//! The free-wake march loop for [VpmRotor] -- the per-sub-step time advance:
//! blade loads via the nonlinear lifting line, shed/trailed vorticity, the
//! rigid-flap and servo-flap DOF integration, and wake convection. Factored out
//! of mod.rs because it is long and self-contained; the public entry points
//! (step_one, march, simulate, and the AeroModel impl) live in the
//! module root and delegate here via [VpmRotor::march_window].

use super::aging::{core_spread, strength_decay};
use super::common::{
    advect_rk2, advect_rk2_bh, advect_rk2_bh_seq, advect_rk2_nan_check, advect_rk2_seq,
    induced_at_points, induced_at_points_bh, induced_at_points_bh_seq, induced_at_points_nan_check,
    induced_at_points_seq, ParticleField,
};
use super::merge::{merge_particles, MergeOpts};
use super::reformulated::advect_rvpm;
use super::{FlightCondition, VpmRotor, VpmRotorResult, VpmRotorState, WakeEngine};
use crate::cyclic::cyclic_coeffs;
use crate::polar::Polar;
use std::f64::consts::PI;
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
    /// Free-wake march for an explicit number of sub-steps, averaging loads
    /// over the trailing `avg_window` steps. Each sub-step convects the wake
    /// by `dt` seconds; `dpsi = omega * dt`. When `warm` carries a persisted
    /// wake / circulation the march continues from it.
    /// `psi_offset` is the blade-0 azimuth at the start of this window
    /// (use 0.0 for a fresh settle; use `state.psi` when continuing).
    pub(super) fn march_window(
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

        // Per-blade feathering DOF (Kaman servo-flap path). Active whenever the
        // rotor carries a ServoFlapActuation (feathering + damper architecture).
        // In servo mode the swashplate collective/cyclic are reinterpreted as
        // flap deflection commands delta_f, which produce the pitching moment
        // M_servo that drives feathering; the feathering angle theta_f then
        // REPLACES the direct swashplate-to-pitch path. The ODE integrated is
        //   I_theta*theta'' + C_theta*theta' + k_aero*theta = M_servo
        // with the mechanical damper C_theta the only dissipation and k_aero
        // the aerodynamic spring from the AC offset -- integrated
        // semi-implicitly below.
        // Servo-flap mode is active whenever the rotor carries a
        // ServoFlapActuation; there is no stiffness gate.
        let servo_active = self.feather.is_some();
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

        self.dbg_validate_setup(
            &wake, omega, dt, dpsi, n, nb, total_steps, &gamma_prev, theta_1c, theta_1s,
            servo_active, warm.is_some(),
        );

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
                let ind = if cfg.use_scalar_nan_check {
                    induced_at_points_nan_check(&wake, &tx, &ty, &tz)
                } else if use_bh {
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
                self.dbg_check_u_far(step, b, &u_far, wake.len(), omega);

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
                self.dbg_check_kj_inputs(step, b, &gam, n);
                for iter in 0..max_iter {
                    let mut converged = true;
                    for i in 0..n {
                        let (_urel, _phi, u_mag, cl, _cd) = section(&gam, i);
                        let c = self.chord[i];
                        let g_new = 0.5 * u_mag * c * cl;
                        let g_relaxed = gam[i] + cfg.relax * (g_new - gam[i]);
                        if cfg.use_scalar_nan_check {
                            Self::dbg_check_kj_step(
                                step, b, iter, i, u_mag, cl, g_new, g_relaxed, c, &u_far[i], &gam,
                            );
                        }
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
                    if cfg.use_scalar_nan_check {
                        self.dbg_check_final(step, b, i, &urel, u_mag, cl, gam[i], &th, &u_far[i]);
                    }
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
                    // the flap span: dM = q_dyn * c * C_M_delta * delta_f * dr.
                    // Also add the blade zero-lift AC moment over the whole span:
                    // dM = q_dyn * c * blade_Cm_AC * dr (DC only, all stations).
                    // Together these give the physical DC balance: feathering
                    // trims at the collective where flap + blade moments cancel.
                    if servo_active {
                        let act = self.feather.as_ref().unwrap();
                        if r >= flap_r_in && r <= flap_r_out {
                            m_servo[b] += q_dyn * c * flap_cm_delta * delta_f_b * self.dr[i];
                        }
                        if act.blade_Cm_AC.abs() > 1e-12 {
                            m_servo[b] += q_dyn * c * act.blade_Cm_AC * self.dr[i];
                        }
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

                self.dbg_check_shed_inputs(step, b, psi_b, &rh, beta_b, &gamma[b], &u_rel[b]);

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
                    if cfg.use_scalar_nan_check {
                        Self::dbg_check_shed_trail_pre(step, b, j, r_edge, &ur, &seg);
                    }
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
                    if cfg.use_scalar_nan_check {
                        Self::dbg_check_shed_trail_post(
                            step, b, j, &pos, &a, r_edge, &rh, beta_b, g_trail, &seg, &u_rel[b],
                        );
                    }
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
                    if cfg.use_scalar_nan_check {
                        self.dbg_check_shed_span(
                            step, b, i, &pos, &a, r, &rh, beta_b, mag, d_gamma, &gamma[b],
                        );
                    }
                    wake.push(pos, a, self.sigma_mid[i]);
                }
            }

            // Save relaxed circulation for the next step's dGamma/dt.
            for b in 0..nb {
                gamma_prev[b].copy_from_slice(&gamma[b]);
            }

            // Integrate the rigid-flap DOF one sub-step (see `integrate_flap`).
            if flap_active {
                self.integrate_flap(&mut beta, &mut beta_dot, &m_flap, omega, dt);
            }

            // Integrate the servo-flap feathering DOF one sub-step
            // (see `integrate_feather`).
            if servo_active {
                self.integrate_feather(&mut theta_f, &mut theta_f_dot, &m_servo, fc.rho, omega, dt);
            }

            // Convect the wake and apply truncation, optional aging, and
            // optional population control (see `advance_wake`).
            self.advance_wake(&mut wake, fc, dt, dpsi, psi0, step);

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
            theta_f_dot: if servo_active {
                Some(theta_f_dot)
            } else {
                None
            },
        };
        (result, out_state)
    }

    /// Integrate the rigid-flap DOF one sub-step. Structural/inertial ODE
    /// forced by the aero flap moment `m_flap` (the aero damping is already
    /// inside `m_flap` via the beta_dot AoA term):
    ///   `I_beta*beta'' = M_flap - I_beta*(Omega^2 + omega_NR^2)*beta`
    /// Symplectic (semi-implicit) Euler: advance the rate first, then the angle
    /// with the new rate -- stable for the lightly-damped flap oscillator at
    /// the resolved sub-step (dpsi ~ 0.1-0.3 rad). Requires `self.flap`.
    fn integrate_flap(
        &self,
        beta: &mut [f64],
        beta_dot: &mut [f64],
        m_flap: &[f64],
        omega: f64,
        dt: f64,
    ) {
        let fp = self
            .flap
            .as_ref()
            .expect("integrate_flap requires FlapProperties");
        let i_beta = fp.I_blade_flap_kgm2.max(1e-9);
        let omega_nr = fp.omega_nr_rad_s;
        // Effective rotating flap stiffness / inertia:
        //   K/I = Omega^2 + omega_NR^2  (centrifugal + structural spring)
        //        = Omega^2 * nu_beta^2
        let k_over_i = omega * omega + omega_nr * omega_nr;
        for b in 0..beta.len() {
            let beta_ddot = m_flap[b] / i_beta - k_over_i * beta[b];
            beta_dot[b] += dt * beta_ddot;
            beta[b] += dt * beta_dot[b];
        }
    }

    /// Integrate the servo-flap feathering DOF one sub-step:
    ///   `I_theta*theta'' + C_theta*theta' + k_aero*theta = M_servo`
    /// The mechanical damper C_theta is the ONLY dissipation (Kaman axis at the
    /// AC => no aero pitch damping). `k_aero` is the aerodynamic restoring
    /// spring from the AC offset (physical & measurable): a pitch-up makes more
    /// lift at the AC, a distance ac_offset aft of the feathering axis, giving a
    /// nose-down restoring torque. With ac_offset=0 (axis at AC) there is no
    /// spring and the damper alone sets the ~90 deg cyclic phase lag; DC trim
    /// then comes from the blade camber moment (folded into `m_servo`).
    /// Integrated semi-implicitly (implicit damping) for unconditional
    /// stability regardless of damper strength. Requires `self.feather`.
    fn integrate_feather(
        &self,
        theta_f: &mut [f64],
        theta_f_dot: &mut [f64],
        m_servo: &[f64],
        rho: f64,
        omega: f64,
        dt: f64,
    ) {
        let act = self
            .feather
            .as_ref()
            .expect("integrate_feather requires ServoFlapActuation");
        let i_th = act.I_theta_kgm2.max(1e-12);
        let c_th = act.damper_Nms_per_rad;
        // Aerodynamic feathering spring from the AC offset [N*m/rad]:
        //   k_aero = 0.5*rho*omega^2*cl_alpha*ac_offset*Int(c r^2 dr)
        let k_aero = 0.5
            * rho
            * omega
            * omega
            * self.cl_alpha
            * act.ac_offset_m
            * self.feather_span_integral;
        let damp_fac = 1.0 / (1.0 + dt * c_th / i_th);
        for b in 0..theta_f.len() {
            let rhs = (m_servo[b] - k_aero * theta_f[b]) / i_th;
            theta_f_dot[b] = (theta_f_dot[b] + dt * rhs) * damp_fac;
            theta_f[b] += dt * theta_f_dot[b];
        }
    }

    /// Advance the free wake one sub-step: convect it with the selected engine
    /// (classic RK2 direct / Barnes-Hut, or the reformulated VPM), FIFO-truncate
    /// to `max_particles`, then apply the optional wake aging (core spreading /
    /// strength fade) and optional tree-collapse population control. `psi0` and
    /// `dpsi` phase the azimuth-triggered merge pass; `step` is only used in the
    /// debug NaN-check messages.
    fn advance_wake(
        &self,
        wake: &mut ParticleField,
        fc: &FlightCondition,
        dt: f64,
        dpsi: f64,
        psi0: f64,
        step: usize,
    ) {
        let cfg = &self.config;
        let freestream = [fc.v_hub[0] as f32, fc.v_hub[1] as f32, fc.v_hub[2] as f32];
        match cfg.wake_engine {
            WakeEngine::ReformulatedVpm => {
                // rVPM: strengths and cores evolve by vortex stretching.
                // Direct O(N^2); Barnes-Hut / nan-check paths do not apply.
                advect_rvpm(wake, freestream, dt as f32);
            }
            WakeEngine::ClassicVpm => {
                if cfg.use_scalar_nan_check {
                    advect_rk2_nan_check(wake, freestream, dt as f32);
                } else if cfg.barnes_hut && wake.len() >= cfg.bh_min_particles {
                    if cfg.use_rayon {
                        advect_rk2_bh(wake, freestream, dt as f32, cfg.bh_theta);
                    } else {
                        advect_rk2_bh_seq(wake, freestream, dt as f32, cfg.bh_theta);
                    }
                } else if cfg.use_rayon {
                    advect_rk2(wake, freestream, dt as f32);
                } else {
                    advect_rk2_seq(wake, freestream, dt as f32);
                }
            }
        }
        self.dbg_check_advect(step, wake);
        if wake.len() > cfg.max_particles {
            let excess = wake.len() - cfg.max_particles;
            drain_front(wake, excess);
        }

        // Optional wake aging (models the missing viscous decay). Core spreading
        // grows each core (conserves circulation); strength fade decays each
        // strength to 1/e over `strength_decay_tau_rev` revs (non-conservative).
        // Both no-op at their zero defaults.
        if cfg.core_spread_nu > 0.0 {
            core_spread(wake, cfg.core_spread_nu, dt);
        }
        if cfg.strength_decay_tau_rev > 0.0 {
            // Decay fraction over this sub-step: exp(-(dpsi/2pi)/tau_rev).
            let frac_rev = dpsi / (2.0 * PI);
            let factor = (-frac_rev / cfg.strength_decay_tau_rev).exp() as f32;
            strength_decay(wake, factor);
        }

        // Optional population control: collapse small, coherent, far-field wake
        // cells into single equivalent particles. Triggered on accumulated
        // azimuth (not the per-call sub-step index) so it fires identically
        // whether the caller drives one sub-step per `step()` call or a whole
        // window via `march()`. Total vector circulation is conserved exactly.
        if cfg.merge_wake && cfg.merge_every > 0 && wake.len() >= cfg.merge_min_particles {
            let merge_az = cfg.merge_every as f64 * dpsi;
            if merge_az > 0.0 {
                let before = (psi0 / merge_az).floor();
                let after = ((psi0 + dpsi) / merge_az).floor();
                if after > before {
                    let opts = MergeOpts {
                        kappa: cfg.merge_kappa,
                        chi_min: cfg.merge_chi_min,
                        region_dist: cfg.merge_region_dist,
                        min_particles: cfg.merge_min_particles,
                    };
                    *wake = merge_particles(wake, &opts);
                }
            }
        }
    }

    // ---- debug-only NaN / sanity instrumentation --------------------------
    // Hoisted off the hot march loop so it reads as physics. Each is gated on
    // `config.use_scalar_nan_check` (enabled at runtime via NAN_DEBUG=1). The
    // helpers called from the innermost per-element loops keep the `if` guard at
    // the call site so there is no call overhead when disabled; the per-step /
    // per-blade ones self-guard on entry. All are `#[cold]`/`#[inline(never)]`
    // to keep them out of the march loop's codegen.

    /// Debug-only march setup validation: the scalar kinematics and every
    /// geometry / wake-strength entry must be finite before the first sub-step.
    #[cold]
    #[inline(never)]
    #[allow(clippy::too_many_arguments)]
    fn dbg_validate_setup(
        &self,
        wake: &ParticleField,
        omega: f64,
        dt: f64,
        dpsi: f64,
        n: usize,
        nb: usize,
        total_steps: usize,
        gamma_prev: &[Vec<f64>],
        theta_1c: f64,
        theta_1s: f64,
        servo_active: bool,
        warm: bool,
    ) {
        if !self.config.use_scalar_nan_check {
            return;
        }
        eprintln!("march_window: n={} nb={} omega={} dt={} dpsi={} total_steps={} wake_len={} warm={}",
            n, nb, omega, dt, dpsi, total_steps, wake.len(), warm);
        eprintln!("  gamma_prev[0][..3]={:?}", &gamma_prev[0][..3.min(n)]);
        eprintln!(
            "  theta_1c={} theta_1s={} servo_active={}",
            theta_1c, theta_1s, servo_active
        );
        assert!(
            omega.is_finite() && omega > 0.0,
            "march_window: bad omega={}",
            omega
        );
        assert!(dt.is_finite() && dt > 0.0, "march_window: bad dt={}", dt);
        assert!(
            dpsi.is_finite() && dpsi > 0.0,
            "march_window: bad dpsi={}",
            dpsi
        );
        assert!(n > 0, "march_window: n_elements=0");
        for i in 0..n {
            assert!(
                self.r_mid[i].is_finite() && self.r_mid[i] > 0.0,
                "march_window: r_mid[{}]={}",
                i,
                self.r_mid[i]
            );
            assert!(
                self.chord[i].is_finite() && self.chord[i] > 0.0,
                "march_window: chord[{}]={}",
                i,
                self.chord[i]
            );
            assert!(
                self.dr[i].is_finite() && self.dr[i] > 0.0,
                "march_window: dr[{}]={}",
                i,
                self.dr[i]
            );
            assert!(
                self.sigma_mid[i].is_finite() && self.sigma_mid[i] > 0.0,
                "march_window: sigma_mid[{}]={}",
                i,
                self.sigma_mid[i]
            );
        }
        // Check wake particle strengths -- huge strengths cause huge u_far.
        for p in 0..wake.len() {
            let amag = ((wake.ax[p] as f64).powi(2)
                + (wake.ay[p] as f64).powi(2)
                + (wake.az[p] as f64).powi(2))
            .sqrt();
            assert!(
                amag.is_finite() && amag < 1e6,
                "march_window: wake particle {} has |a|={} a=[{},{},{}]",
                p,
                amag,
                wake.ax[p],
                wake.ay[p],
                wake.az[p]
            );
            assert!(
                wake.px[p].is_finite() && wake.py[p].is_finite() && wake.pz[p].is_finite(),
                "march_window: wake particle {} has NaN pos=[{},{},{}]",
                p,
                wake.px[p],
                wake.py[p],
                wake.pz[p]
            );
        }
    }

    /// Debug-only far-field induced-velocity sanity check: flags any non-finite
    /// or absurdly large (`>500x` tip speed) wake-induced velocity at a blade's
    /// stations, and traces the worst-case magnitude.
    #[cold]
    #[inline(never)]
    fn dbg_check_u_far(
        &self,
        step: usize,
        b: usize,
        u_far: &[[f64; 3]],
        wake_len: usize,
        omega: f64,
    ) {
        if !self.config.use_scalar_nan_check {
            return;
        }
        let n = u_far.len();
        let r_tip = *self.r_edge.last().unwrap();
        let u_max = 500.0 * omega * r_tip; // 500x tip speed is absurd
        for i in 0..n {
            for k in 0..3 {
                assert!(u_far[i][k].is_finite() && u_far[i][k].abs() < u_max,
                    "u_far LARGE: step={} b={} i={} k={} u_far={:.3} (limit {:.1}) wake_len={}",
                    step, b, i, k, u_far[i][k], u_max, wake_len);
            }
        }
        // Also eprintln the worst-case magnitude for tracing.
        let u_far_max: f64 = u_far
            .iter()
            .flat_map(|v| v.iter().map(|x| x.abs()))
            .fold(0.0, f64::max);
        if u_far_max > omega * r_tip * 5.0 {
            eprintln!(
                "  u_far warn: step={} b={} u_far_max={:.2} (5x tip_speed={:.1})",
                step,
                b,
                u_far_max,
                omega * r_tip * 5.0
            );
        }
    }

    /// Debug-only pre-solve check: the seed circulation and every geometry entry
    /// feeding the Kutta-Joukowski fixed point must be finite and positive.
    #[cold]
    #[inline(never)]
    fn dbg_check_kj_inputs(&self, step: usize, b: usize, gam: &[f64], n: usize) {
        if !self.config.use_scalar_nan_check {
            return;
        }
        for i in 0..n {
            assert!(
                gam[i].is_finite(),
                "kj nan: step={} b={} gamma_prev[{}]={}",
                step,
                b,
                i,
                gam[i]
            );
            assert!(
                self.r_mid[i].is_finite() && self.r_mid[i] > 0.0,
                "kj nan: step={} b={} r_mid[{}]={}",
                step,
                b,
                i,
                self.r_mid[i]
            );
            assert!(
                self.chord[i].is_finite() && self.chord[i] > 0.0,
                "kj nan: step={} b={} chord[{}]={}",
                step,
                b,
                i,
                self.chord[i]
            );
            assert!(
                self.dr[i].is_finite() && self.dr[i] > 0.0,
                "kj nan: step={} b={} dr[{}]={}",
                step,
                b,
                i,
                self.dr[i]
            );
            assert!(
                self.sigma_mid[i].is_finite() && self.sigma_mid[i] > 0.0,
                "kj nan: step={} b={} sigma_mid[{}]={}",
                step,
                b,
                i,
                self.sigma_mid[i]
            );
            assert!(
                self.sigma_edge[i].is_finite() && self.sigma_edge[i] > 0.0,
                "kj nan: step={} b={} sigma_edge[{}]={}",
                step,
                b,
                i,
                self.sigma_edge[i]
            );
        }
        assert!(
            self.sigma_edge[n].is_finite() && self.sigma_edge[n] > 0.0,
            "kj nan: step={} b={} sigma_edge[n={}]={}",
            step,
            b,
            n,
            self.sigma_edge[n]
        );
    }

    /// Debug-only per-iteration KJ check (guard kept at the call site because
    /// this runs in the innermost fixed-point loop).
    #[cold]
    #[inline(never)]
    #[allow(clippy::too_many_arguments)]
    fn dbg_check_kj_step(
        step: usize,
        b: usize,
        iter: usize,
        i: usize,
        u_mag: f64,
        cl: f64,
        g_new: f64,
        g_relaxed: f64,
        c: f64,
        u_far_i: &[f64; 3],
        gam: &[f64],
    ) {
        assert!(
            u_mag.is_finite() && u_mag < 1e8,
            "kj nan: step={} b={} iter={} i={} u_mag={} u_far={:?} gam={:?}",
            step,
            b,
            iter,
            i,
            u_mag,
            u_far_i,
            gam
        );
        assert!(
            cl.is_finite(),
            "kj nan: step={} b={} iter={} i={} cl={} u_mag={}",
            step,
            b,
            iter,
            i,
            cl,
            u_mag
        );
        assert!(g_relaxed.is_finite() && g_relaxed.abs() < 1e8,
            "kj nan: step={} b={} iter={} i={} g_new={} g_relaxed={} u_mag={} cl={} c={}",
            step, b, iter, i, g_new, g_relaxed, u_mag, cl, c);
    }

    /// Debug-only converged-loads check (guard kept at the call site,
    /// per-element).
    #[cold]
    #[inline(never)]
    #[allow(clippy::too_many_arguments)]
    fn dbg_check_final(
        &self,
        step: usize,
        b: usize,
        i: usize,
        urel: &[f64; 3],
        u_mag: f64,
        cl: f64,
        gam_i: f64,
        th: &[f64; 3],
        u_far_i: &[f64; 3],
    ) {
        for k in 0..3 {
            assert!(urel[k].is_finite() && urel[k].abs() < 1e8,
                "final nan: step={} b={} i={} urel[{}]={} u_mag={} cl={} gam[i]={} r_mid={} th={:?} u_far={:?}",
                step, b, i, k, urel[k], u_mag, cl, gam_i, self.r_mid[i], th, u_far_i);
        }
        assert!(
            gam_i.is_finite() && gam_i.abs() < 1e8,
            "final nan: step={} b={} i={} gam[i]={}",
            step,
            b,
            i,
            gam_i
        );
    }

    /// Debug-only shed-vorticity input check: azimuth, radial unit vector, flap
    /// angle, and every bound-circulation / relative-wind entry must be finite
    /// before shedding this blade's wake.
    #[cold]
    #[inline(never)]
    fn dbg_check_shed_inputs(
        &self,
        step: usize,
        b: usize,
        psi_b: f64,
        rh: &[f64; 3],
        beta_b: f64,
        gamma_b: &[f64],
        u_rel_b: &[[f64; 3]],
    ) {
        if !self.config.use_scalar_nan_check {
            return;
        }
        assert!(
            psi_b.is_finite(),
            "shed nan: step={} b={} psi_b={}",
            step,
            b,
            psi_b
        );
        assert!(
            rh[0].is_finite() && rh[1].is_finite(),
            "shed nan: step={} b={} rh={:?}",
            step,
            b,
            rh
        );
        assert!(
            beta_b.is_finite(),
            "shed nan: step={} b={} beta_b={}",
            step,
            b,
            beta_b
        );
        for i in 0..gamma_b.len() {
            assert!(
                gamma_b[i].is_finite(),
                "shed nan: step={} b={} gamma[{}]={}",
                step,
                b,
                i,
                gamma_b[i]
            );
            for k in 0..3 {
                assert!(
                    u_rel_b[i][k].is_finite(),
                    "shed nan: step={} b={} u_rel[{}][{}]={}",
                    step,
                    b,
                    i,
                    k,
                    u_rel_b[i][k]
                );
            }
        }
    }

    /// Debug-only trailed-particle pre-check (guard kept at call site): edge
    /// radius and streamwise segment must be finite before building the particle.
    #[cold]
    #[inline(never)]
    fn dbg_check_shed_trail_pre(
        step: usize,
        b: usize,
        j: usize,
        r_edge: f64,
        ur: &[f64; 3],
        seg: &[f64; 3],
    ) {
        assert!(
            r_edge.is_finite(),
            "shed trail nan: step={} b={} j={} r_edge={}",
            step,
            b,
            j,
            r_edge
        );
        assert!(
            seg[0].is_finite() && seg[1].is_finite() && seg[2].is_finite(),
            "shed trail nan: step={} b={} j={} ur={:?} seg={:?}",
            step,
            b,
            j,
            ur,
            seg
        );
    }

    /// Debug-only trailed-particle post-check (guard kept at call site): the
    /// built position / strength must be finite and not implausibly large.
    #[cold]
    #[inline(never)]
    #[allow(clippy::too_many_arguments)]
    fn dbg_check_shed_trail_post(
        step: usize,
        b: usize,
        j: usize,
        pos: &[f32; 3],
        a: &[f32; 3],
        r_edge: f64,
        rh: &[f64; 3],
        beta_b: f64,
        g_trail: f64,
        seg: &[f64; 3],
        u_rel_b: &[[f64; 3]],
    ) {
        assert!(pos[0].is_finite() && pos[1].is_finite() && pos[2].is_finite(),
            "shed trail nan: step={} b={} j={} pos={:?} r_edge={} rh={:?} beta_b={}", step, b, j, pos, r_edge, rh, beta_b);
        assert!(
            a[0].is_finite() && a[1].is_finite() && a[2].is_finite(),
            "shed trail nan: step={} b={} j={} a={:?} g_trail={} seg={:?}",
            step,
            b,
            j,
            a,
            g_trail,
            seg
        );
        let amag =
            ((a[0] as f64).powi(2) + (a[1] as f64).powi(2) + (a[2] as f64).powi(2)).sqrt();
        assert!(amag < 1e4,
            "shed trail LARGE: step={} b={} j={} |a|={} a={:?} g_trail={} seg={:?} u_rel_b0={:?}",
            step, b, j, amag, a, g_trail, seg, &u_rel_b[..3.min(u_rel_b.len())]);
    }

    /// Debug-only shed-spanwise particle check (guard kept at call site): the
    /// built position / strength must be finite and bounded.
    #[cold]
    #[inline(never)]
    #[allow(clippy::too_many_arguments)]
    fn dbg_check_shed_span(
        &self,
        step: usize,
        b: usize,
        i: usize,
        pos: &[f32; 3],
        a: &[f32; 3],
        r: f64,
        rh: &[f64; 3],
        beta_b: f64,
        mag: f64,
        d_gamma: f64,
        gamma_b: &[f64],
    ) {
        assert!(
            pos[0].is_finite() && pos[1].is_finite() && pos[2].is_finite(),
            "shed span nan: step={} b={} i={} pos={:?} r={} rh={:?} beta_b={}",
            step,
            b,
            i,
            pos,
            r,
            rh,
            beta_b
        );
        assert!(
            a[0].is_finite() && a[1].is_finite() && a[2].is_finite(),
            "shed span nan: step={} b={} i={} a={:?} mag={} d_gamma={} dr={}",
            step,
            b,
            i,
            a,
            mag,
            d_gamma,
            self.dr[i]
        );
        let amag =
            ((a[0] as f64).powi(2) + (a[1] as f64).powi(2) + (a[2] as f64).powi(2)).sqrt();
        assert!(amag < 1e4,
            "shed span LARGE: step={} b={} i={} |a|={} a={:?} mag={} d_gamma={} gamma={:?}",
            step, b, i, amag, a, mag, d_gamma, gamma_b);
    }

    /// Debug-only post-advection position check.
    #[cold]
    #[inline(never)]
    fn dbg_check_advect(&self, step: usize, wake: &ParticleField) {
        if !self.config.use_scalar_nan_check {
            return;
        }
        for p in 0..wake.len() {
            assert!(
                wake.px[p].is_finite() && wake.py[p].is_finite() && wake.pz[p].is_finite(),
                "advect nan: step={} particle={} pos=[{},{},{}]",
                step,
                p,
                wake.px[p],
                wake.py[p],
                wake.pz[p]
            );
        }
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
