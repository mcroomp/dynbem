// Forward-flight free-wake VPM rotor coupling.
//
// Extends the axial coupling to cyclic pitch and crosswind / edgewise flow,
// following standard unsteady lifting-line free-wake practice (Leishman,
// "Principles of Helicopter Aerodynamics"; Bagai & Leishman 1995). 
// See VPM_DESIGN.md and validation_rs for comprehensive theory validation tests.
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
use crate::bem_common::{assemble_result, kinematics, Kinematics};
use crate::cyclic::{cyclic_coeffs, ControlGains};
use crate::polar::Polar;
use crate::rotor_definition::{
    FlapProperties, PitchActuation, RotorDefinition, ServoFlapActuation,
};
use crate::vpm::{
    advect_rk2, advect_rk2_bh, advect_rk2_seq, advect_rk2_bh_seq,
    induced_at_points, induced_at_points_bh, induced_at_points_seq, induced_at_points_bh_seq,
    ParticleField,
};
use std::f64::consts::PI;

/// Solver resolution / wake settings.
#[derive(Clone, Copy, Debug)]
pub struct VpmRotorConfig {
    /// Maximum number of wake particles retained. Older particles are
    /// discarded (front-truncation) once this cap is reached. Size this to
    /// cover the desired wake age: a rough guide is
    /// `n_blades * (2*n_elements + 1) * steps_per_rev * wake_revs`.
    pub max_particles: usize,
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
    /// Model per-blade rigid flap dynamics (beta, beta_dot) when the rotor
    /// definition supplies `FlapProperties`. When false (or no FlapProperties)
    /// blades are rigid in flap and the wake stays in the disk plane.
    pub flap_dynamics: bool,
    /// Use Rayon parallelism for wake induction and convection. When false,
    /// uses sequential implementations for debugging or single-threaded execution.
    pub use_rayon: bool,
}

impl Default for VpmRotorConfig {
    fn default() -> Self {
        Self {
            max_particles: 4800,
            sigma: 0.18,
            relax: 0.35,
            nonlinear_lifting_line: true,
            tip_clustering: true,
            local_core: true,
            barnes_hut: false,
            bh_theta: 0.5,
            bh_min_particles: 2048,
            flap_dynamics: true,
            use_rayon: true,
        }
    }
}

impl VpmRotorConfig {
    /// Small, fast preset for unit tests (coarse but sign-correct).
    pub fn fast_test() -> Self {
        Self {
            max_particles: 800,
            sigma: 0.2,
            relax: 0.4,
            nonlinear_lifting_line: true,
            tip_clustering: true,
            local_core: true,
            barnes_hut: false,
            bh_theta: 0.5,
            bh_min_particles: 2048,
            flap_dynamics: true,
            use_rayon: true,
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
pub struct VpmRotor<P: Polar> {
    polar: P,
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
    /// Optional rigid-flap properties (inertia + non-rotating flap frequency).
    /// `Some` enables per-blade flap dynamics when `config.flap_dynamics`.
    flap: Option<FlapProperties>,
    /// Optional servo-flap feathering actuation (Kaman path). `Some` enables
    /// the per-blade feathering DOF when the control stiffness is positive.
    feather: Option<ServoFlapActuation>,
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

impl<P: Polar> VpmRotor<P> {
    pub fn new(
        defn: &RotorDefinition,
        polar: P,
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
            polar,
            n_blades: blade.n_blades,
            n_elements: n,
            r_edge,
            r_mid,
            dr,
            chord,
            twist,
            sigma_mid,
            sigma_edge,
            flap: defn.flap.clone(),
            feather: match &defn.pitch_actuation {
                PitchActuation::ServoFlap(act) => Some(act.clone()),
                PitchActuation::DirectMechanical => None,
            },
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
    /// `dt` is the sub-step duration (seconds) -- use the same value as your
    /// controller loop (e.g. 1.0/400.0). `n_steps` is how many sub-steps to
    /// run in total; loads are averaged over the second half.
    pub fn simulate(&self, fc: &FlightCondition, dt: f64, n_steps: usize) -> VpmRotorResult {
        self.march(fc, None, dt, n_steps).0
    }

    /// Advance the free wake by exactly one sub-step of duration `dt`.
    /// Call this once per controller tick; read `out_state.wake` afterwards
    /// for the current particle cloud.
    pub fn step_one(
        &self,
        fc: &FlightCondition,
        state: &VpmRotorState,
        dt: f64,
    ) -> (VpmRotorResult, VpmRotorState) {
        self.march_window(fc, Some(state), state.psi, dt, 1, 1)
    }

    /// Settle the free wake for `n_steps` sub-steps of duration `dt`, then
    /// return loads averaged over the second half. When `warm` carries a
    /// persisted wake the march continues from it. Typical usage: call once
    /// at startup with enough steps to develop a periodic wake (~several
    /// hundred at 400 Hz for a 30 rad/s rotor), then drive with `step()`.
    pub fn march(
        &self,
        fc: &FlightCondition,
        warm: Option<&VpmRotorState>,
        dt: f64,
        n_steps: usize,
    ) -> (VpmRotorResult, VpmRotorState) {
        let avg_window = (n_steps / 2).max(1);
        self.march_window(fc, warm, 0.0, dt, n_steps, avg_window)
    }

    /// Free-wake march for an explicit number of sub-steps, averaging loads
    /// over the trailing `avg_window` steps. Each sub-step convects the wake
    /// by `dt` seconds; `dpsi = omega * dt`. When `warm` carries a persisted
    /// wake / circulation the march continues from it.
    /// `psi_offset` is the blade-0 azimuth at the start of this window
    /// (use 0.0 for a fresh settle; use `state.psi` when continuing).
    fn march_window(
        &self,
        fc: &FlightCondition,
        warm: Option<&VpmRotorState>,
        psi_offset: f64,
        dt: f64,
        total_steps: usize,
        avg_window: usize,
    ) -> (VpmRotorResult, VpmRotorState) {
        let n = self.n_elements;
        let nb = self.n_blades;
        let cfg = &self.config;
        let omega = fc.omega_rad_s;
        // dpsi = omega * dt: angle swept per sub-step at current rotor speed.
        let dpsi = omega * dt;
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

        // Per-blade rigid-flap DOF. Active only when the rotor supplies
        // FlapProperties AND config.flap_dynamics is on. beta > 0 = flap up
        // (tip toward -Z, the thrust side). The aerodynamic flap damping is
        // NOT added as an analytical term here -- it emerges from the loads,
        // because beta_dot feeds the section angle of attack below (the whole
        // reason to model flap in the time-marched VPM rather than as a static
        // hub-moment factor). The ODE integrated per step is therefore the
        // purely structural/inertial one, forced by the aero flap moment:
        //   I_beta * beta'' + I_beta*(Omega^2 + omega_NR^2)*beta = M_flap_aero
        let flap_active = cfg.flap_dynamics && self.flap.is_some();
        let mut beta = match warm.and_then(|s| s.beta.as_ref()) {
            Some(b) if flap_active && b.len() == nb => b.clone(),
            _ => vec![0.0f64; nb],
        };
        let mut beta_dot = match warm.and_then(|s| s.beta_dot.as_ref()) {
            Some(b) if flap_active && b.len() == nb => b.clone(),
            _ => vec![0.0f64; nb],
        };

        // Per-blade feathering DOF (Kaman servo-flap path). Active when the
        // rotor carries a ServoFlapActuation with a positive control stiffness
        // (the pushrod/linkage restoring moment -- feathering has no
        // centrifugal stiffening, so without it the DOF is ill-posed). In
        // servo mode the swashplate collective/cyclic are reinterpreted as
        // flap deflection commands delta_f, which produce the pitching moment
        // M_servo that drives feathering; the feathering angle theta_f then
        // REPLACES the direct swashplate-to-pitch path. The ODE integrated is
        //   I_theta*theta'' + C_theta*theta' + k_ctrl*theta = M_servo
        // with the mechanical damper C_theta the only dissipation (axis at the
        // AC => no aero pitch damping) -- integrated semi-implicitly below.
        let servo_active = match &self.feather {
            Some(act) => act.control_stiffness_Nm_per_rad > 0.0,
            None => false,
        };
        // Swashplate commands reinterpreted as flap deflection harmonics.
        let (delta_f0, delta_f1c, delta_f1s) = (fc.collective_rad, theta_1c, theta_1s);
        let (flap_r_in, flap_r_out, flap_cm_delta) = match &self.feather {
            Some(act) => (
                act.flap.r_inner_m,
                act.flap.r_outer_m,
                act.flap.C_M_delta_per_rad,
            ),
            None => (0.0, 0.0, 0.0),
        };
        let mut theta_f = match warm.and_then(|s| s.theta_f.as_ref()) {
            Some(t) if servo_active && t.len() == nb => t.clone(),
            _ => vec![0.0f64; nb],
        };
        let mut theta_f_dot = match warm.and_then(|s| s.theta_f_dot.as_ref()) {
            Some(t) if servo_active && t.len() == nb => t.clone(),
            _ => vec![0.0f64; nb],
        };

        // Trailed edge particles (n+1) + shed station particles (n), per blade.
        let max_particles = cfg.max_particles;

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

            // Aerodynamic flap moment about the hinge, per blade, accumulated
            // in the loads loop below (M = sum r * dF_z). Drives the flap ODE.
            let mut m_flap = vec![0.0f64; nb];
            // Servo-flap pitching moment about the feathering axis, per blade,
            // accumulated over the flap span below. Drives the feathering ODE.
            let mut m_servo = vec![0.0f64; nb];

            // ---- loads on every blade (nonlinear lifting-line solve) ------
            for b in 0..nb {
                let psi_b = psi0 + b as f64 * 2.0 * PI / nb as f64;
                let rh = r_hat(psi_b);
                let th = t_hat(psi_b);
                let (spsi, cpsi) = psi_b.sin_cos();
                let cyc = theta_1c * cpsi + theta_1s * spsi;
                // Flap DOF for this blade (0 when inactive). beta tilts the
                // blade out of plane (z = -r*beta, up = -Z); beta_dot adds a
                // vertical section velocity that changes the local AoA (the
                // aerodynamic flap damping, captured exactly here).
                let beta_b = beta[b];
                let beta_dot_b = beta_dot[b];
                // Feathering DOF for this blade. In servo mode the blade pitch
                // comes from theta_f (which subsumes collective + cyclic via
                // the servo-flap response); otherwise from the swashplate.
                let theta_f_b = theta_f[b];
                let delta_f_b = if servo_active {
                    delta_f0 + delta_f1c * cpsi + delta_f1s * spsi
                } else {
                    0.0
                };
                let control_pitch = if servo_active {
                    theta_f_b
                } else {
                    fc.collective_rad + cyc
                };

                // Far-field (particle wake) induced velocity at station centers.
                let (tx, ty, tz): (Vec<f32>, Vec<f32>, Vec<f32>) = {
                    let mut xs = Vec::with_capacity(n);
                    let mut ys = Vec::with_capacity(n);
                    let mut zs = Vec::with_capacity(n);
                    for i in 0..n {
                        let r = self.r_mid[i];
                        xs.push((r * rh[0]) as f32);
                        ys.push((r * rh[1]) as f32);
                        zs.push((-r * beta_b) as f32);
                    }
                    (xs, ys, zs)
                };
                let use_bh = cfg.barnes_hut && wake.len() >= cfg.bh_min_particles;
                let ind = if use_bh {
                    if cfg.use_rayon {
                        induced_at_points_bh(&wake, &tx, &ty, &tz, cfg.bh_theta)
                    } else {
                        induced_at_points_bh_seq(&wake, &tx, &ty, &tz, cfg.bh_theta)
                    }
                } else if cfg.use_rayon {
                    induced_at_points(&wake, &tx, &ty, &tz)
                } else {
                    induced_at_points_seq(&wake, &tx, &ty, &tz)
                };
                let u_far: Vec<[f64; 3]> = (0..n)
                    .map(|i| [ind[i][0] as f64, ind[i][1] as f64, ind[i][2] as f64])
                    .collect();

                // Background (far-field only) relative wind, used to fix the
                // near-wake trailing-leg geometry for this step.
                let mut urel_bg = vec![[0.0f64; 3]; n];
                for i in 0..n {
                    let r = self.r_mid[i];
                    let vb = omega * r;
                    // Flap-rate vertical velocity: blade section moves at
                    // -r*beta_dot in Z (up = -Z), so the relative wind gains
                    // +r*beta_dot in Z. Increasing u_a lowers the AoA -> flap
                    // aerodynamic damping.
                    let v_flap = r * beta_dot_b;
                    urel_bg[i] = [
                        fc.v_hub[0] + u_far[i][0] - vb * th[0],
                        fc.v_hub[1] + u_far[i][1] - vb * th[1],
                        fc.v_hub[2] + u_far[i][2] - vb * th[2] + v_flap,
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
                        ep.push([r_edge * rh[0], r_edge * rh[1], -r_edge * beta_b]);
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
                        let cp = [r * rh[0], r * rh[1], -r * beta_b];
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
                    // Flap-rate vertical velocity (see urel_bg above): +r*beta_dot
                    // in Z, so flapping up reduces the AoA (flap damping).
                    let v_flap = r * beta_dot_b;
                    let urel = [
                        fc.v_hub[0] + u_far[i][0] + u_ind[0] - vb * th[0],
                        fc.v_hub[1] + u_far[i][1] + u_ind[1] - vb * th[1],
                        fc.v_hub[2] + u_far[i][2] + u_ind[2] - vb * th[2] + v_flap,
                    ];
                    let u_a = urel[2];
                    let u_t = -(urel[0] * th[0] + urel[1] * th[1]);
                    let u_mag = (u_a * u_a + u_t * u_t).sqrt().max(1e-6);
                    let phi = u_a.atan2(u_t);
                    let twist = self.twist[i];
                    let alpha = twist + control_pitch - phi;
                    let (cl, cd) = self.polar.cl_cd(alpha);
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
                    // Aero flap moment about the hinge: out-of-plane force
                    // (d_thrust, up) at arm r. Positive -> flaps blade up.
                    m_flap[b] += r * d_thrust;
                    // Servo-flap pitching moment about the feathering axis over
                    // the flap span: dM = q_dyn * c * C_M_delta * delta_f * dr,
                    // using the true local dynamic pressure (VPM sees the real
                    // flow, not the Omega*r approximation the BEM solve uses).
                    if servo_active && r >= flap_r_in && r <= flap_r_out {
                        m_servo[b] += q_dyn * c * flap_cm_delta * delta_f_b * self.dr[i];
                    }

                    gamma[b][i] = gam[i];
                }
            }

            // ---- shed vorticity from every blade --------------------------
            for b in 0..nb {
                let psi_b = psi0 + b as f64 * 2.0 * PI / nb as f64;
                let rh = r_hat(psi_b);
                // Out-of-plane flap displacement of this blade (z = -r*beta).
                let beta_b = beta[b];

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
                    let pos = [
                        (r_edge * rh[0]) as f32,
                        (r_edge * rh[1]) as f32,
                        (-r_edge * beta_b) as f32,
                    ];
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
                    let pos = [(r * rh[0]) as f32, (r * rh[1]) as f32, (-r * beta_b) as f32];
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

            // ---- integrate the rigid-flap DOF one sub-step ----------------
            // Structural/inertial ODE forced by the aero flap moment (the aero
            // damping is already inside m_flap via the beta_dot AoA term):
            //   I_beta*beta'' = M_flap - I_beta*(Omega^2 + omega_NR^2)*beta
            // Symplectic (semi-implicit) Euler: advance the rate first, then
            // the angle with the new rate -- stable for the lightly-damped
            // flap oscillator at the resolved sub-step (dpsi ~ 0.1-0.3 rad).
            if flap_active {
                let fp = self.flap.as_ref().expect("flap_active implies Some");
                let i_beta = fp.I_blade_flap_kgm2.max(1e-9);
                let omega_nr = fp.omega_nr_rad_s;
                // Effective rotating flap stiffness / inertia:
                //   K/I = Omega^2 + omega_NR^2  (centrifugal + structural spring)
                //        = Omega^2 * nu_beta^2
                let k_over_i = omega * omega + omega_nr * omega_nr;
                for b in 0..nb {
                    let beta_ddot = m_flap[b] / i_beta - k_over_i * beta[b];
                    beta_dot[b] += dt * beta_ddot;
                    beta[b] += dt * beta_dot[b];
                }
            }

            // ---- integrate the feathering DOF one sub-step ----------------
            // I_theta*theta'' + C_theta*theta' + k_ctrl*theta = M_servo
            // The mechanical damper C_theta is the ONLY dissipation (Kaman axis
            // at the AC => no aero pitch damping), so it is integrated
            // semi-implicitly (implicit on the damping term) for unconditional
            // stability regardless of how stiff the damper is:
            //   theta_dot <- (theta_dot + dt*(M_servo - k*theta)/I) / (1 + dt*C/I)
            //   theta     <- theta + dt*theta_dot
            if servo_active {
                let act = self.feather.as_ref().expect("servo_active implies Some");
                let i_th = act.I_theta_kgm2.max(1e-12);
                let c_th = act.damper_Nms_per_rad;
                let k_ctrl = act.control_stiffness_Nm_per_rad;
                let damp_fac = 1.0 / (1.0 + dt * c_th / i_th);
                for b in 0..nb {
                    let rhs = (m_servo[b] - k_ctrl * theta_f[b]) / i_th;
                    theta_f_dot[b] = (theta_f_dot[b] + dt * rhs) * damp_fac;
                    theta_f[b] += dt * theta_f_dot[b];
                }
            }

            // ---- convect and truncate the free wake -----------------------
            let freestream = [fc.v_hub[0] as f32, fc.v_hub[1] as f32, fc.v_hub[2] as f32];
            if cfg.barnes_hut && wake.len() >= cfg.bh_min_particles {
                if cfg.use_rayon {
                    advect_rk2_bh(&mut wake, freestream, dt as f32, cfg.bh_theta);
                } else {
                    advect_rk2_bh_seq(&mut wake, freestream, dt as f32, cfg.bh_theta);
                }
            } else if cfg.use_rayon {
                advect_rk2(&mut wake, freestream, dt as f32);
            } else {
                advect_rk2_seq(&mut wake, freestream, dt as f32);
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
            beta: if flap_active { Some(beta) } else { None },
            beta_dot: if flap_active { Some(beta_dot) } else { None },
            theta_f: if servo_active { Some(theta_f) } else { None },
            theta_f_dot: if servo_active { Some(theta_f_dot) } else { None },
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
    /// Per-blade flap angle (rad), positive = flap up (tip toward -Z, the
    /// thrust side). `None` when flap dynamics are inactive.
    pub beta: Option<Vec<f64>>,
    /// Per-blade flap rate (rad/s). `None` when flap dynamics are inactive.
    pub beta_dot: Option<Vec<f64>>,
    /// Per-blade feathering angle (rad) about the pitch axis, driven by the
    /// servo-flap moment. `None` when the feathering DOF is inactive
    /// (direct-mechanical pitch or zero control stiffness).
    pub theta_f: Option<Vec<f64>>,
    /// Per-blade feathering rate (rad/s). `None` when inactive.
    pub theta_f_dot: Option<Vec<f64>>,
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

impl<P: Polar> AeroModel for VpmRotor<P> {
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
    /// replaced: the wake itself is marched directly. `method` is ignored.
    /// One sub-step of duration `dt` is taken per call -- align `dt` with your
    /// controller loop (e.g. 1/400 s at 400 Hz).
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
        // One sub-step per call; dt IS the sub-step duration.
        let (res, out_state) = self.march_window(&fc, Some(state), state.psi, dt, 1, 1);
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

// ============================================================================
// Public API for validation / analysis
// ============================================================================

/// Compute induced velocities at arbitrary points in space from a wake
/// particle field. Used for analysis and validation (e.g., disk inflow sampling).
pub fn induced_velocities_at_points(
    wake: &ParticleField,
    tx: &[f32],
    ty: &[f32],
    tz: &[f32],
) -> Vec<[f32; 3]> {
    induced_at_points(wake, tx, ty, tz)
}


