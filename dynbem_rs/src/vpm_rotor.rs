// Forward-flight free-wake VPM rotor coupling.
//
// Extends the axial coupling (see examples/vpm_vs_bem_axial.rs) to cyclic
// pitch and crosswind / edgewise flow, following standard unsteady
// lifting-line free-wake practice (Leishman, "Principles of Helicopter
// Aerodynamics"; Bagai & Leishman 1995). See VPM_DESIGN.md.
//
// What this adds over the axial coupling:
//   * Per-blade loading and induced-velocity probing (no azimuthal-symmetry
//     shortcut) -- every blade at its own azimuth.
//   * Cyclic pitch via `cyclic::cyclic_coeffs` (repo swashplate convention).
//   * Full 3D freestream (through-disk + in-plane), so the advancing/
//     retreating asymmetry and the skewed wake fall out naturally.
//   * Shed (temporal) vorticity  -dGamma/dt  in addition to trailed
//     (radial-gradient) vorticity. This is the load-bearing addition: once
//     bound circulation varies with azimuth (any cyclic or crosswind case)
//     Kelvin's theorem requires the shed term, or the wake is wrong.
//   * Hub pitch/roll moments  Mx = sum r dT sin(psi),  My = sum r dT cos(psi)
//     (AGENTS.md hub-frame convention), so cyclic/crosswind are observable.
//
// Frame: NED, CCW-from-above (AGENTS.md). Hub axis +Z, thrust is "up" (-Z).
// Everything here is evaluated in the hub frame (caller supplies the
// freestream already rotated into hub axes).

use crate::aero_io::{AeroResult, RotorInputs};
use crate::aero_model::{AeroModel, IntegrationMethod, RotorStateExt};
use crate::bem_common::{assemble_result, kinematics, Kinematics, PolarTable};
use crate::cyclic::{cyclic_coeffs, ControlGains};
use crate::polar::Polar;
use crate::rotor_definition::RotorDefinition;
use crate::vpm::{
    advect_rk2, advect_rk2_bh, induced_at_points, induced_at_points_bh, ParticleField,
};
use std::f64::consts::PI;

/// Solver resolution / wake settings.
#[derive(Clone, Copy, Debug)]
pub struct VpmRotorConfig {
    pub n_steps_per_rev: usize,
    pub n_wake_rev: usize,
    /// Revolutions to march before averaging the (periodic) loads.
    pub n_settle_rev: usize,
    /// Particle core size (m).
    pub sigma: f32,
    /// Under-relaxation factor on bound circulation.
    pub relax: f64,
    /// Solve the bound circulation with an implicit (Prandtl / Weissinger)
    /// lifting-line, including the within-step trailing near-wake downwash.
    /// When false, falls back to the original pointwise Kutta-Joukowski pass
    /// (no within-step self-induction).
    pub nonlinear_lifting_line: bool,
    /// Cluster the spanwise stations toward the tip (and root) with cosine
    /// spacing, to resolve the bound-circulation roll-off. When false, the
    /// stations are uniformly spaced.
    pub tip_clustering: bool,
    /// Scale each feature's particle / filament core size with the local
    /// radial spacing (smaller cores where stations are clustered), instead
    /// of using the single global `sigma` everywhere.
    pub local_core: bool,
    /// Use the Barnes-Hut O(N log N) wake evaluator instead of the direct
    /// O(N^2) sum (approximate; controlled by `bh_theta`). Only engaged once
    /// the wake reaches `bh_min_particles`.
    pub barnes_hut: bool,
    /// Barnes-Hut opening angle (smaller = more accurate + slower).
    pub bh_theta: f32,
    /// Particle count below which the direct sum is used even with
    /// `barnes_hut` on (the tree overhead is not worth it for small wakes).
    pub bh_min_particles: usize,
}

impl Default for VpmRotorConfig {
    fn default() -> Self {
        Self {
            n_steps_per_rev: 24,
            n_wake_rev: 4,
            n_settle_rev: 6,
            sigma: 0.18,
            relax: 0.35,
            nonlinear_lifting_line: true,
            tip_clustering: true,
            local_core: true,
            barnes_hut: false,
            bh_theta: 0.5,
            bh_min_particles: 2048,
        }
    }
}

impl VpmRotorConfig {
    /// Small, fast preset for unit tests (coarse but sign-correct).
    pub fn fast_test() -> Self {
        Self {
            n_steps_per_rev: 12,
            n_wake_rev: 2,
            n_settle_rev: 3,
            sigma: 0.2,
            relax: 0.4,
            nonlinear_lifting_line: true,
            tip_clustering: true,
            local_core: true,
            barnes_hut: false,
            bh_theta: 0.5,
            bh_min_particles: 2048,
        }
    }
}

/// A single steady flight condition, in hub axes.
#[derive(Clone, Copy, Debug)]
pub struct FlightCondition {
    pub collective_rad: f64,
    pub tilt_lon: f64,
    pub tilt_lat: f64,
    /// Air velocity relative to the hub, hub frame. +Z is through-disk
    /// (axial); X/Y are in-plane (edgewise / crosswind).
    pub v_hub: [f64; 3],
    pub omega_rad_s: f64,
    pub rho: f64,
}

/// Cycle-averaged loads (mean over the final revolution).
#[derive(Clone, Copy, Debug)]
pub struct VpmRotorResult {
    /// Thrust "up" (-Z hub), N.
    pub thrust: f64,
    /// Shaft reaction torque (spin), N*m.
    pub torque: f64,
    /// Hub rolling moment (roll-right positive), N*m.
    pub mx_hub: f64,
    /// Hub pitching moment (pitch-up positive), N*m.
    pub my_hub: f64,
    /// Wake particle count at the end of the run.
    pub n_particles: usize,
    /// Mean wake position (for diagnostics, e.g. skew).
    pub wake_centroid: [f64; 3],
}

/// Forward-flight free-wake VPM rotor.
#[derive(Clone)]
pub struct VpmRotor {
    polar: PolarTable,
    n_blades: usize,
    n_elements: usize,
    /// Radial edges (n+1), possibly cosine-clustered toward the tip.
    r_edge: Vec<f64>,
    /// Station centres (n) -- midpoints of adjacent edges.
    r_mid: Vec<f64>,
    /// Per-station radial width (n) = r_edge[i+1] - r_edge[i].
    dr: Vec<f64>,
    /// Per-station chord (n), sampled at r_mid.
    chord: Vec<f64>,
    /// Per-station twist (rad, n), sampled at r_mid.
    twist: Vec<f64>,
    /// Per-station shed-particle core size (n).
    sigma_mid: Vec<f32>,
    /// Per-edge trailed-particle / near-wake filament core size (n+1).
    sigma_edge: Vec<f32>,
    ctrl: ControlGains,
    config: VpmRotorConfig,
}

#[inline]
fn r_hat(psi: f64) -> [f64; 3] {
    [psi.cos(), -psi.sin(), 0.0]
}
#[inline]
fn t_hat(psi: f64) -> [f64; 3] {
    [-psi.sin(), -psi.cos(), 0.0]
}

/// Regularized Biot-Savart velocity at `p` induced by a straight vortex
/// filament A->B of unit circulation. Van Garrel core smoothing (core radius
/// `core`) keeps it finite near the filament. Right-hand rule about A->B,
/// matching the particle-field sign convention -- multiply the result by the
/// filament circulation.
#[inline]
fn segment_induced(p: [f64; 3], a: [f64; 3], b: [f64; 3], core: f64) -> [f64; 3] {
    let r1 = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let r2 = [p[0] - b[0], p[1] - b[1], p[2] - b[2]];
    let r0 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let cross = [
        r1[1] * r2[2] - r1[2] * r2[1],
        r1[2] * r2[0] - r1[0] * r2[2],
        r1[0] * r2[1] - r1[1] * r2[0],
    ];
    let c2 = cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2];
    let r1n = (r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]).sqrt();
    let r2n = (r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2]).sqrt();
    if c2 < 1e-20 || r1n < 1e-9 || r2n < 1e-9 {
        return [0.0; 3];
    }
    let l2 = r0[0] * r0[0] + r0[1] * r0[1] + r0[2] * r0[2];
    let num = (r0[0] * r1[0] + r0[1] * r1[1] + r0[2] * r1[2]) / r1n
        - (r0[0] * r2[0] + r0[1] * r2[1] + r0[2] * r2[2]) / r2n;
    // Van Garrel regularization: attenuates within ~core of the filament.
    let reg = c2 / (c2 + core * core * l2);
    let s = num * reg / (4.0 * PI * c2);
    [cross[0] * s, cross[1] * s, cross[2] * s]
}

impl VpmRotor {
    pub fn new<P: Polar>(
        defn: &RotorDefinition,
        polar: &P,
        ctrl: ControlGains,
        config: VpmRotorConfig,
    ) -> Self {
        let blade = &defn.blade;
        let n = blade.n_elements;
        let r_root = blade.root_cutout_m;
        let r_tip = blade.radius_m;
        let span = r_tip - r_root;

        // Radial edges: cosine-clustered toward both ends when requested
        // (clusters at the tip and root, resolving the circulation roll-off),
        // otherwise uniform.
        let r_edge: Vec<f64> = (0..=n)
            .map(|k| {
                let f = if config.tip_clustering {
                    0.5 * (1.0 - (PI * k as f64 / n as f64).cos())
                } else {
                    k as f64 / n as f64
                };
                r_root + span * f
            })
            .collect();

        let r_mid: Vec<f64> = (0..n).map(|i| 0.5 * (r_edge[i] + r_edge[i + 1])).collect();
        let dr: Vec<f64> = (0..n).map(|i| r_edge[i + 1] - r_edge[i]).collect();
        let chord: Vec<f64> = r_mid.iter().map(|&r| blade.chord_at(r)).collect();
        let twist: Vec<f64> = r_mid
            .iter()
            .map(|&r| blade.twist_at(r).to_radians())
            .collect();

        // Local core sizing: scale the base sigma by the local radial spacing
        // relative to the uniform spacing, clamped to [0.5, 1.0] x base so the
        // wake stays well-overlapped where stations are coarse and the tip
        // vortex is not over-smoothed where they are fine.
        let dr_uniform = span / n as f64;
        let base = config.sigma;
        let core_for = |width: f64| -> f32 {
            if config.local_core {
                let s = (width / dr_uniform).clamp(0.5, 1.0) as f32;
                base * s
            } else {
                base
            }
        };
        let sigma_mid: Vec<f32> = dr.iter().map(|&w| core_for(w)).collect();
        let sigma_edge: Vec<f32> = (0..=n)
            .map(|k| {
                // Edge width = mean of adjacent station widths.
                let w = if k == 0 {
                    dr[0]
                } else if k == n {
                    dr[n - 1]
                } else {
                    0.5 * (dr[k - 1] + dr[k])
                };
                core_for(w)
            })
            .collect();

        Self {
            polar: PolarTable::from_polar(polar),
            n_blades: blade.n_blades,
            n_elements: n,
            r_edge,
            r_mid,
            dr,
            chord,
            twist,
            sigma_mid,
            sigma_edge,
            ctrl,
            config,
        }
    }

    /// Build the hub-frame flight condition (and the world<->hub kinematics)
    /// from generic `RotorInputs`. Shared by the `AeroModel` load evaluation
    /// (`compute_forces`) and the time-advance (`step`).
    fn flight_condition(&self, inputs: &RotorInputs) -> (FlightCondition, Kinematics) {
        let r_tip = *self.r_edge.last().expect("r_edge always has n+1 >= 2 edges");
        let kin = kinematics(inputs, inputs.omega_rad_s, r_tip);
        let fc = FlightCondition {
            collective_rad: inputs.collective_rad,
            tilt_lon: inputs.tilt_lon,
            tilt_lat: inputs.tilt_lat,
            // Air velocity relative to the hub, hub axes: in-plane (X, Y) from
            // v_inplane_hub, axial (Z, through-disk) from v_climb.
            v_hub: [kin.v_inplane_hub[0], kin.v_inplane_hub[1], kin.v_climb],
            omega_rad_s: inputs.omega_rad_s,
            rho: inputs.rho_kg_m3,
        };
        (fc, kin)
    }

    /// March the free wake from empty to a periodic state and return
    /// cycle-averaged loads (cold start; no persisted wake).
    pub fn simulate(&self, fc: &FlightCondition) -> VpmRotorResult {
        self.march(fc, None).0
    }

    /// Advance the free wake by exactly one azimuthal sub-step
    /// (`dpsi = 2*pi / n_steps_per_rev`).  Call this once per animation
    /// frame; read `out_state.wake` immediately afterwards for the current
    /// particle cloud.
    pub fn step_one(
        &self,
        fc: &FlightCondition,
        state: &VpmRotorState,
    ) -> (VpmRotorResult, VpmRotorState) {
        self.march_window(fc, Some(state), state.psi, 1, 1)
    }

    /// Sub-step duration in seconds (`dpsi / omega`).
    pub fn dt_step(&self, fc: &FlightCondition) -> f64 {
        (2.0 * PI / self.config.n_steps_per_rev as f64) / fc.omega_rad_s.max(1e-9)
    }

    /// Core free-wake march. When `warm` carries a persisted wake and bound
    /// circulation the march *continues* from them instead of starting from an
    /// empty wake; the final wake and circulation are returned in the output
    /// `VpmRotorState`. Runs a full settle (`(n_settle_rev + 1)` revs) and
    /// averages loads over the final revolution. Standard unsteady
    /// lifting-line free-wake time-march.
    pub fn march(
        &self,
        fc: &FlightCondition,
        warm: Option<&VpmRotorState>,
    ) -> (VpmRotorResult, VpmRotorState) {
        let cfg = &self.config;
        let total_steps = (cfg.n_settle_rev + 1) * cfg.n_steps_per_rev;
        let avg_window = cfg.n_steps_per_rev;
        self.march_window(fc, warm, 0.0, total_steps, avg_window)
    }

    /// Free-wake march for an explicit number of sub-steps, averaging loads
    /// over the trailing `avg_window` steps. Each sub-step convects the wake
    /// by `dt = (2*pi / n_steps_per_rev) / omega`. When `warm` carries a
    /// persisted wake / circulation the march continues from it; the final
    /// wake and circulation are returned in the output `VpmRotorState`.
    /// `psi_offset` is the blade-0 azimuth at the start of this window
    /// (use 0.0 for a fresh settle; use `state.psi` when continuing).
    fn march_window(
        &self,
        fc: &FlightCondition,
        warm: Option<&VpmRotorState>,
        psi_offset: f64,
        total_steps: usize,
        avg_window: usize,
    ) -> (VpmRotorResult, VpmRotorState) {
        let n = self.n_elements;
        let nb = self.n_blades;
        let cfg = &self.config;
        let dpsi = 2.0 * PI / cfg.n_steps_per_rev as f64;
        let omega = fc.omega_rad_s;
        let dt = dpsi / omega;
        let (theta_1c, theta_1s) = cyclic_coeffs(fc.tilt_lon, fc.tilt_lat, self.ctrl);

        let mut wake = match warm.and_then(|s| s.wake.clone()) {
            Some(w) => w,
            None => ParticleField::new(),
        };
        // Relaxed bound circulation and its previous value, per (blade, station).
        let mut gamma = vec![vec![0.0f64; n]; nb];
        // Seed the previous-step circulation from the persisted state when its
        // shape matches; otherwise start from zero (cold blade).
        let mut gamma_prev = match warm.and_then(|s| s.gamma.as_ref()) {
            Some(g) if g.len() == nb && g.iter().all(|row| row.len() == n) => g.clone(),
            _ => vec![vec![0.0f64; n]; nb],
        };

        // Trailed edge particles (n+1) + shed station particles (n), per blade.
        let shed_per_step = nb * (2 * n + 1);
        let max_particles = cfg.n_wake_rev * cfg.n_steps_per_rev * shed_per_step;

        // `total_steps` sub-steps are taken; loads are averaged over the final
        // `avg_window` of them (both are caller-supplied).
        let avg_start = total_steps.saturating_sub(avg_window);

        let mut t_acc = 0.0;
        let mut q_acc = 0.0;
        let mut mx_acc = 0.0;
        let mut my_acc = 0.0;
        let mut avg_count = 0usize;

        // Scratch reused each step.
        let mut u_rel = vec![vec![[0.0f64; 3]; n]; nb]; // per-station relative wind

        for step in 0..total_steps {
            let psi0 = psi_offset + step as f64 * dpsi;

            let mut thrust_step = 0.0;
            let mut torque_step = 0.0;
            let mut mx_step = 0.0;
            let mut my_step = 0.0;

            // ---- loads on every blade (nonlinear lifting-line solve) ------
            for b in 0..nb {
                let psi_b = psi0 + b as f64 * 2.0 * PI / nb as f64;
                let rh = r_hat(psi_b);
                let th = t_hat(psi_b);
                let (spsi, cpsi) = psi_b.sin_cos();
                let cyc = theta_1c * cpsi + theta_1s * spsi;

                // Far-field (particle wake) induced velocity at station centers.
                let (tx, ty, tz): (Vec<f32>, Vec<f32>, Vec<f32>) = {
                    let mut xs = Vec::with_capacity(n);
                    let mut ys = Vec::with_capacity(n);
                    let mut zs = Vec::with_capacity(n);
                    for i in 0..n {
                        let r = self.r_mid[i];
                        xs.push((r * rh[0]) as f32);
                        ys.push((r * rh[1]) as f32);
                        zs.push(0.0);
                    }
                    (xs, ys, zs)
                };
                let use_bh = cfg.barnes_hut && wake.len() >= cfg.bh_min_particles;
                let ind = if use_bh {
                    induced_at_points_bh(&wake, &tx, &ty, &tz, cfg.bh_theta)
                } else {
                    induced_at_points(&wake, &tx, &ty, &tz)
                };
                let u_far: Vec<[f64; 3]> = (0..n)
                    .map(|i| [ind[i][0] as f64, ind[i][1] as f64, ind[i][2] as f64])
                    .collect();

                // Background (far-field only) relative wind, used to fix the
                // near-wake trailing-leg geometry for this step.
                let mut urel_bg = vec![[0.0f64; 3]; n];
                for i in 0..n {
                    let vb = omega * self.r_mid[i];
                    urel_bg[i] = [
                        fc.v_hub[0] + u_far[i][0] - vb * th[0],
                        fc.v_hub[1] + u_far[i][1] - vb * th[1],
                        fc.v_hub[2] + u_far[i][2] - vb * th[2],
                    ];
                }

                // Near-wake influence coefficients B[i][k]: 3D velocity at
                // station i's control point induced by a unit-circulation
                // one-step trailing filament at edge k (k = 0..=n). The leg
                // trails along the local relative wind, one convection step
                // long -- exactly the trailed vorticity this step deposits,
                // which becomes a particle next step (so no double count with
                // the far wake). Lift comes from the polar, so the bound
                // vortex carries no self term (classical Prandtl lifting line).
                let mut b_inf: Vec<Vec<[f64; 3]>> = Vec::new();
                if cfg.nonlinear_lifting_line {
                    let mut ep = Vec::with_capacity(n + 1);
                    let mut leg = Vec::with_capacity(n + 1);
                    for k in 0..=n {
                        let r_edge = self.r_edge[k];
                        ep.push([r_edge * rh[0], r_edge * rh[1], 0.0]);
                        let ur = if k == 0 {
                            urel_bg[0]
                        } else if k == n {
                            urel_bg[n - 1]
                        } else {
                            [
                                0.5 * (urel_bg[k - 1][0] + urel_bg[k][0]),
                                0.5 * (urel_bg[k - 1][1] + urel_bg[k][1]),
                                0.5 * (urel_bg[k - 1][2] + urel_bg[k][2]),
                            ]
                        };
                        leg.push([ur[0] * dt, ur[1] * dt, ur[2] * dt]);
                    }
                    b_inf = vec![vec![[0.0f64; 3]; n + 1]; n];
                    for i in 0..n {
                        let r = self.r_mid[i];
                        let cp = [r * rh[0], r * rh[1], 0.0];
                        for k in 0..=n {
                            let a = ep[k];
                            let bb = [a[0] + leg[k][0], a[1] + leg[k][1], a[2] + leg[k][2]];
                            b_inf[i][k] = segment_induced(cp, a, bb, self.sigma_edge[k] as f64);
                        }
                    }
                }

                // Trailing-wake downwash from the current circulation.
                let trail_downwash = |gam: &[f64], i: usize| -> [f64; 3] {
                    let mut u = [0.0f64; 3];
                    if cfg.nonlinear_lifting_line {
                        for k in 0..=n {
                            let g_in = if k == 0 { 0.0 } else { gam[k - 1] };
                            let g_out = if k == n { 0.0 } else { gam[k] };
                            let g_trail = g_in - g_out;
                            u[0] += b_inf[i][k][0] * g_trail;
                            u[1] += b_inf[i][k][1] * g_trail;
                            u[2] += b_inf[i][k][2] * g_trail;
                        }
                    }
                    u
                };

                // Section state at a given circulation: returns (urel, phi,
                // u_mag, cl, cd) using the total induced velocity.
                let section = |gam: &[f64], i: usize| -> ([f64; 3], f64, f64, f64, f64) {
                    let u_ind = trail_downwash(gam, i);
                    let r = self.r_mid[i];
                    let vb = omega * r;
                    let urel = [
                        fc.v_hub[0] + u_far[i][0] + u_ind[0] - vb * th[0],
                        fc.v_hub[1] + u_far[i][1] + u_ind[1] - vb * th[1],
                        fc.v_hub[2] + u_far[i][2] + u_ind[2] - vb * th[2],
                    ];
                    let u_a = urel[2];
                    let u_t = -(urel[0] * th[0] + urel[1] * th[1]);
                    let u_mag = (u_a * u_a + u_t * u_t).sqrt().max(1e-6);
                    let phi = u_a.atan2(u_t);
                    let twist = self.twist[i];
                    let alpha = fc.collective_rad + twist + cyc - phi;
                    let (cl, cd) = self.polar.interp(alpha);
                    (urel, phi, u_mag, cl, cd)
                };

                // Implicit fixed-point solve for bound circulation, seeded from
                // the previous step. With the near wake disabled this collapses
                // to one relaxed Kutta-Joukowski pass (the original scheme).
                let max_iter = if cfg.nonlinear_lifting_line { 30 } else { 1 };
                let mut gam = gamma_prev[b].clone();
                for _ in 0..max_iter {
                    let mut converged = true;
                    for i in 0..n {
                        let (_urel, _phi, u_mag, cl, _cd) = section(&gam, i);
                        let c = self.chord[i];
                        let g_new = 0.5 * u_mag * c * cl;
                        let g_relaxed = gam[i] + cfg.relax * (g_new - gam[i]);
                        if (g_relaxed - gam[i]).abs() > 1e-6 * (g_relaxed.abs() + 1e-6) {
                            converged = false;
                        }
                        gam[i] = g_relaxed;
                    }
                    if converged {
                        break;
                    }
                }

                // Final pass: loads at the converged circulation.
                for i in 0..n {
                    let (urel, phi, u_mag, cl, cd) = section(&gam, i);
                    u_rel[b][i] = urel;
                    let r = self.r_mid[i];
                    let c = self.chord[i];
                    let q_dyn = 0.5 * fc.rho * u_mag * u_mag;
                    let dl = q_dyn * c * cl * self.dr[i];
                    let dd = q_dyn * c * cd * self.dr[i];
                    let d_thrust = dl * phi.cos() - dd * phi.sin(); // up (-Z)
                    thrust_step += d_thrust;
                    torque_step += (dl * phi.sin() + dd * phi.cos()) * r;
                    // Hub moments (AGENTS.md): Mx = r dT sin psi, My = r dT cos psi.
                    mx_step += r * d_thrust * spsi;
                    my_step += r * d_thrust * cpsi;

                    gamma[b][i] = gam[i];
                }
            }

            // ---- shed vorticity from every blade --------------------------
            for b in 0..nb {
                let psi_b = psi0 + b as f64 * 2.0 * PI / nb as f64;
                let rh = r_hat(psi_b);

                // Trailed vorticity: edge j strength = Gamma_{j-1} - Gamma_j,
                // Gamma outside the blade = 0. Segment aligned with the local
                // relative wind (streamwise), length |U_rel| dt.
                for j in 0..=n {
                    let g_in = if j == 0 { 0.0 } else { gamma[b][j - 1] };
                    let g_out = if j == n { 0.0 } else { gamma[b][j] };
                    let g_trail = g_in - g_out;
                    if g_trail == 0.0 {
                        continue;
                    }
                    // Relative wind at the edge = mean of adjacent stations.
                    let ur = if j == 0 {
                        u_rel[b][0]
                    } else if j == n {
                        u_rel[b][n - 1]
                    } else {
                        [
                            0.5 * (u_rel[b][j - 1][0] + u_rel[b][j][0]),
                            0.5 * (u_rel[b][j - 1][1] + u_rel[b][j][1]),
                            0.5 * (u_rel[b][j - 1][2] + u_rel[b][j][2]),
                        ]
                    };
                    let seg = [ur[0] * dt, ur[1] * dt, ur[2] * dt];
                    let r_edge = self.r_edge[j];
                    let pos = [(r_edge * rh[0]) as f32, (r_edge * rh[1]) as f32, 0.0];
                    let a = [
                        (g_trail * seg[0]) as f32,
                        (g_trail * seg[1]) as f32,
                        (g_trail * seg[2]) as f32,
                    ];
                    wake.push(pos, a, self.sigma_edge[j]);
                }

                // Shed (temporal) vorticity: spanwise segment of strength
                // -(Gamma^n - Gamma^{n-1}) along the blade radial, per station.
                for i in 0..n {
                    let d_gamma = gamma[b][i] - gamma_prev[b][i];
                    if d_gamma == 0.0 {
                        continue;
                    }
                    let r = self.r_mid[i];
                    // Spanwise vortex line: strength -dGamma, length dr along r_hat.
                    let mag = -d_gamma * self.dr[i];
                    let pos = [(r * rh[0]) as f32, (r * rh[1]) as f32, 0.0];
                    let a = [
                        (mag * rh[0]) as f32,
                        (mag * rh[1]) as f32,
                        (mag * rh[2]) as f32,
                    ];
                    wake.push(pos, a, self.sigma_mid[i]);
                }
            }

            // Save relaxed circulation for the next step's dGamma/dt.
            for b in 0..nb {
                gamma_prev[b].copy_from_slice(&gamma[b]);
            }

            // ---- convect and truncate the free wake -----------------------
            let freestream = [fc.v_hub[0] as f32, fc.v_hub[1] as f32, fc.v_hub[2] as f32];
            if cfg.barnes_hut && wake.len() >= cfg.bh_min_particles {
                advect_rk2_bh(&mut wake, freestream, dt as f32, cfg.bh_theta);
            } else {
                advect_rk2(&mut wake, freestream, dt as f32);
            }
            if wake.len() > max_particles {
                let excess = wake.len() - max_particles;
                drain_front(&mut wake, excess);
            }

            if step >= avg_start {
                t_acc += thrust_step;
                q_acc += torque_step;
                mx_acc += mx_step;
                my_acc += my_step;
                avg_count += 1;
            }
        }

        let inv = 1.0 / avg_count as f64;
        let centroid = wake_centroid(&wake);
        let result = VpmRotorResult {
            thrust: t_acc * inv,
            torque: q_acc * inv,
            mx_hub: mx_acc * inv,
            my_hub: my_acc * inv,
            n_particles: wake.len(),
            wake_centroid: centroid,
        };
        let out_state = VpmRotorState {
            wake: Some(wake),
            gamma: Some(gamma),
            psi: (psi_offset + total_steps as f64 * dpsi).rem_euclid(2.0 * PI),
        };
        (result, out_state)
    }
}

/// Persistent state for the VPM rotor, carried between `compute_forces` calls.
///
/// The free-wake particle cloud and the relaxed bound circulation live here as
/// plain struct fields so the wake is *continued* across calls instead of being
/// rebuilt from empty each time. They are deliberately NOT exposed through the
/// `RotorStateExt` inflow vector: `get_inflow` returns an empty `Vec`, so the
/// envelope integrator sees zero inflow DOFs and never tries to damp, blend, or
/// serialize the wake (a particle cloud has no meaningful scalar-inflow tau).
#[derive(Clone, Debug, Default)]
pub struct VpmRotorState {
    /// Persisted free-wake particle cloud. `None` on the first call -> the
    /// march starts from an empty wake.
    pub wake: Option<ParticleField>,
    /// Persisted relaxed bound circulation per (blade, station). Seeds the next
    /// call's implicit circulation solve and its dGamma/dt shed vorticity.
    pub gamma: Option<Vec<Vec<f64>>>,
    /// Current blade-0 azimuth (rad). Accumulated by `step_one` so that each
    /// call sheds particles from the correct rotor position rather than always
    /// restarting at psi=0. `march` ignores this field and always starts at 0.
    pub psi: f64,
}

impl RotorStateExt for VpmRotorState {
    // The wake is persisted via the struct fields above, not the inflow vector,
    // so the inflow DOF count is zero.
    fn get_inflow(&self) -> Vec<f64> {
        Vec::new()
    }
    fn set_inflow(&mut self, arr: Vec<f64>) {
        debug_assert!(arr.is_empty());
    }
}

impl AeroModel for VpmRotor {
    type State = VpmRotorState;

    fn initial_state(&self) -> Self::State {
        VpmRotorState::default()
    }

    // inflow_taus: the trait default (all-infinity) is correct -- there are no
    // scalar inflow DOFs (the persistent wake lives in the struct fields).

    fn compute_forces(
        &self,
        _inputs: &RotorInputs,
        _state: &VpmRotorState,
    ) -> (AeroResult, VpmRotorState) {
        // VPM is a time-marching free-wake model: loads cannot be produced from
        // a single instantaneous evaluation the way the BEM-family models can.
        // The wake must be advanced in time, and the settling loop is the
        // caller's responsibility -- drive it with `step(dt)` repeatedly (N
        // sub-steps to settle) and read the loads from each `AeroResult`.
        //
        // This is the "assert false on the instantaneous path" contract: there
        // is no correct single-shot compute_forces for a free wake, so calling
        // it is a programming error. (`panic!` rather than `assert!(false, ..)`
        // only to avoid clippy::assertions_on_constants.)
        panic!(
            "VpmRotor::compute_forces is unsupported: a free-wake VPM rotor has \
             no single-shot instantaneous evaluation. Advance the wake with \
             AeroModel::step(dt) N times (the caller owns the settling loop); \
             for a self-contained steady solve use VpmRotor::simulate / march."
        );
    }

    /// Advance the free wake by `dt` seconds and return the loads averaged over
    /// the advanced window.
    ///
    /// This is the only valid time-advance for a VPM rotor: `compute_forces`
    /// panics because a free wake has no instantaneous single-shot evaluation.
    /// The caller settles the wake by calling `step` repeatedly (N sub-steps).
    ///
    /// VPM carries the whole free wake as its state, not scalar inflow DOFs, so
    /// the generic derivative-integration path in the default `AeroModel::step`
    /// (`lam_{n+1} = lam_n + dt*dlam`) is meaningless here and is deliberately
    /// replaced: the wake itself is marched directly. `method` is ignored for
    /// the same reason. The wake convects on its own fixed sub-step
    /// `dt_wake = (2*pi / n_steps_per_rev) / omega`; this advances the smallest
    /// whole number of sub-steps (>= 1) that covers `dt`.
    fn step(
        &self,
        inputs: &RotorInputs,
        state: &VpmRotorState,
        dt: f64,
        _method: IntegrationMethod,
    ) -> (AeroResult, VpmRotorState) {
        // The derivative-integration path must never run for the VPM model:
        // there are no scalar inflow states to integrate. Enforce the invariant
        // that makes this override the only correct time-advance for the wake.
        assert!(
            state.get_inflow().is_empty(),
            "VPM step: derivative-integration path is invalid -- VpmRotorState \
             carries the free wake, not scalar inflow DOFs"
        );
        let (fc, kin) = self.flight_condition(inputs);
        let dpsi = 2.0 * PI / self.config.n_steps_per_rev as f64;
        let dt_wake = dpsi / fc.omega_rad_s.abs().max(1e-9);
        let n_sub = ((dt / dt_wake).round() as i64).max(1) as usize;
        let (res, out_state) = self.march_window(&fc, Some(state), state.psi, n_sub, n_sub);
        let result = assemble_result(
            res.thrust,
            res.torque,
            res.mx_hub,
            res.my_hub,
            kin.hub_axis,
            &inputs.R_hub,
        );
        (result, out_state)
    }
}

/// Remove the oldest `k` particles (front of the SoA arrays).
fn drain_front(f: &mut ParticleField, k: usize) {
    // `AVec` has no `drain`, so shift the survivors to the front and truncate.
    for v in [
        &mut f.px,
        &mut f.py,
        &mut f.pz,
        &mut f.ax,
        &mut f.ay,
        &mut f.az,
        &mut f.sigma,
    ] {
        let n = v.len();
        v.copy_within(k.., 0);
        v.truncate(n - k);
    }
}

fn wake_centroid(f: &ParticleField) -> [f64; 3] {
    let n = f.len();
    if n == 0 {
        return [0.0; 3];
    }
    let inv = 1.0 / n as f64;
    [
        f.px.iter().map(|&v| v as f64).sum::<f64>() * inv,
        f.py.iter().map(|&v| v as f64).sum::<f64>() * inv,
        f.pz.iter().map(|&v| v as f64).sum::<f64>() * inv,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aero_io::{Mat3, RotorInputs, Vec3};
    use crate::aero_model::AeroModel;
    use crate::polar::LinearPolar;
    use crate::quasi_static_bem::{QuasiStaticBEM, QuasiStaticRotorState};
    use crate::rotor_definition::{
        BladeGeometry, LinearPolarParameters, PitchActuation, RotorDefinition,
    };

    const OMEGA: f64 = 120.0;
    const RHO: f64 = 1.225;

    fn test_rotor(n_elements: usize) -> RotorDefinition {
        RotorDefinition {
            blade: BladeGeometry {
                n_blades: 2,
                radius_m: 1.0,
                root_cutout_m: 0.2,
                chord_m: 0.06,
                twist_deg: 2.0,
                n_elements,
                tip_loss: true,
                r_stations_m: Vec::new(),
                chord_stations_m: Vec::new(),
                twist_stations_deg: Vec::new(),
            },
            airfoil: LinearPolarParameters {
                CL0: 0.0,
                CL_alpha_per_rad: 5.7,
                CD0: 0.01,
                alpha_stall_deg: 15.0,
            },
            control: None,
            pitch_actuation: PitchActuation::DirectMechanical,
            flap: None,
            name: "vpm_rotor_test".to_string(),
            description: String::new(),
        }
    }

    fn polar() -> LinearPolar {
        LinearPolar::new(0.0, 5.7, 0.01, 15.0_f64.to_radians())
    }

    fn rotor(n_elements: usize) -> VpmRotor {
        VpmRotor::new(
            &test_rotor(n_elements),
            &polar(),
            ControlGains::default(),
            VpmRotorConfig::fast_test(),
        )
    }

    fn hover(collective_deg: f64) -> FlightCondition {
        FlightCondition {
            collective_rad: collective_deg.to_radians(),
            tilt_lon: 0.0,
            tilt_lat: 0.0,
            v_hub: [0.0, 0.0, 0.0],
            omega_rad_s: OMEGA,
            rho: RHO,
        }
    }

    /// With no cyclic and no crosswind the forward-flight model reduces to
    /// axial flow and should reproduce a hover thrust in the same ballpark as
    /// the quasi-static BEM (the shed term vanishes at steady state, so this
    /// exercises the per-blade loop and trailed wake). Hub moments should be
    /// ~ 0 by axisymmetry.
    #[test]
    fn axial_reduction_matches_bem_hover() {
        let vpm = rotor(12);
        let res = vpm.simulate(&hover(8.0));

        // BEM anchor.
        let bem = QuasiStaticBEM::build(test_rotor(30), 72, polar());
        let inputs = RotorInputs {
            collective_rad: 8.0_f64.to_radians(),
            tilt_lon: 0.0,
            tilt_lat: 0.0,
            R_hub: Mat3::eye(),
            v_hub_world: Vec3::zero(),
            wind_world: Vec3::zero(),
            rho_kg_m3: RHO,
            omega_rad_s: OMEGA,
        };
        let (bem_res, _) = bem.compute_forces(&inputs, &QuasiStaticRotorState);
        let t_bem = -bem_res.F_world[2];

        assert!(res.thrust > 0.0, "thrust should be positive, got {}", res.thrust);
        let rel = (res.thrust - t_bem).abs() / t_bem;
        assert!(
            rel < 0.30,
            "VPM hover thrust {:.1} N vs BEM {:.1} N ({:.0}% off)",
            res.thrust,
            t_bem,
            rel * 100.0
        );
        // Axisymmetric -> hub moments small relative to thrust*R.
        let scale = res.thrust * 1.0;
        assert!(
            res.mx_hub.abs() < 0.10 * scale && res.my_hub.abs() < 0.10 * scale,
            "hub moments should be ~0 in hover: Mx={:.3} My={:.3}",
            res.mx_hub,
            res.my_hub
        );
    }

    /// Thrust increases with collective pitch.
    #[test]
    fn collective_increases_thrust() {
        let vpm = rotor(10);
        let lo = vpm.simulate(&hover(5.0)).thrust;
        let hi = vpm.simulate(&hover(9.0)).thrust;
        assert!(hi > lo, "thrust should rise with collective: {} -> {}", lo, hi);
    }

    /// Longitudinal cyclic tilt_lon > 0 (forward stick) produces a nose-down
    /// pitching moment: My_hub < 0 (AGENTS.md convention).
    #[test]
    fn longitudinal_cyclic_gives_nose_down_moment() {
        let vpm = rotor(8);
        let base = vpm.simulate(&hover(8.0));
        let mut fc = hover(8.0);
        fc.tilt_lon = 3.0_f64.to_radians();
        let tilted = vpm.simulate(&fc);
        assert!(
            tilted.my_hub < base.my_hub && tilted.my_hub < 0.0,
            "tilt_lon>0 should give nose-down My<0: base {:.3} -> {:.3}",
            base.my_hub,
            tilted.my_hub
        );
    }

    /// Lateral cyclic tilt_lat > 0 produces a roll-right moment: Mx_hub > 0.
    #[test]
    fn lateral_cyclic_gives_roll_right_moment() {
        let vpm = rotor(8);
        let base = vpm.simulate(&hover(8.0));
        let mut fc = hover(8.0);
        fc.tilt_lat = 3.0_f64.to_radians();
        let tilted = vpm.simulate(&fc);
        assert!(
            tilted.mx_hub > base.mx_hub && tilted.mx_hub > 0.0,
            "tilt_lat>0 should give roll-right Mx>0: base {:.3} -> {:.3}",
            base.mx_hub,
            tilted.mx_hub
        );
    }

    /// A crosswind (edgewise in-plane flow) keeps the solution bounded, keeps
    /// thrust positive, develops a hub moment, and skews the wake downstream.
    #[test]
    fn crosswind_stays_bounded_and_skews_wake() {
        let vpm = rotor(10);
        let mut fc = hover(8.0);
        fc.v_hub = [8.0, 0.0, 0.0]; // 8 m/s edgewise along +X
        let res = vpm.simulate(&fc);

        assert!(res.thrust.is_finite() && res.thrust > 0.0, "thrust {}", res.thrust);
        assert!(res.torque.is_finite(), "torque {}", res.torque);
        // Asymmetric loading -> nonzero hub moment.
        let moment = (res.mx_hub.powi(2) + res.my_hub.powi(2)).sqrt();
        assert!(moment > 1e-3, "crosswind should induce a hub moment, got {}", moment);
        // Wake convects downstream (+X) with the freestream.
        assert!(
            res.wake_centroid[0] > 0.05,
            "wake should skew downstream (+X), centroid_x = {}",
            res.wake_centroid[0]
        );
    }
}
