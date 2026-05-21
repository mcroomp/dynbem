// Level 2: Pitt-Peters 3-state dynamic inflow + VRS empirical correction.
// See ../CLAUDE.md "Pitt-Peters inflow ODE" and dynbem/pitt_peters.py.

use std::f64::consts::PI;

use crate::aero_io::{AeroResult, RotorInputs, Vec3};
use crate::aero_model::AeroModel;
use crate::bem_common::RadialGrid;
use crate::common::{
    vrs_lambda1, EPS_DENOM, EPS_OMEGA_R, MU_T_FLOOR, VRS_DESCENT_THRESHOLD, V_T_HOVER_FLOOR_FRAC,
};
use crate::cyclic::cyclic_coeffs;
use crate::polar::{Polar, PolarKind};
use crate::rotor_definition::RotorDefinition;
use crate::rotor_state::PittPetersRotorState;

/// Sum thrust and torque over radial elements in axial flight (mu = 0).
/// Plain scalar loop; LLVM auto-vectorizes the arithmetic between the
/// transcendentals.
#[allow(clippy::too_many_arguments)]
fn axial_forces(
    grid: &RadialGrid,
    col: f64,
    omega: f64,
    omega_r: f64,
    lambda_total: f64,
    rho: f64,
    n_b: usize,
    polar: &PolarKind,
) -> (f64, f64) {
    let mut t_acc = 0.0;
    let mut q_acc = 0.0;
    let v_a = lambda_total * omega_r;
    let n_bf = n_b as f64;
    for i in 0..grid.r_mid.len() {
        let r = grid.r_mid[i];
        let v_t = omega * r;
        if v_t < EPS_DENOM {
            continue;
        }
        let phi = v_a.atan2(v_t);
        let alpha = col + grid.twist_rad[i] - phi;
        let (cl, cd) = polar.cl_cd(alpha);
        let cos_p = phi.cos();
        let sin_p = phi.sin();
        let cn = cl * cos_p - cd * sin_p;
        let ct = cl * sin_p + cd * cos_p;
        let q = 0.5 * rho * (v_a * v_a + v_t * v_t) * grid.chord[i] * grid.dr * n_bf;
        t_acc += q * cn;
        q_acc += q * ct * r;
    }
    (t_acc, q_acc)
}

/// Sum thrust, torque, and in-plane hub moments over (psi, r).
/// Pitt-Peters local inflow:
///     lam_local = lambda_total + x*(lam_c*cos psi + lam_s*sin psi).
/// CCW-from-above; v_t_extra = v_in_hub_x*sin psi + v_in_hub_y*cos psi.
/// Per-element moment: r*dT*[sin psi, cos psi, 0] in hub frame.
#[allow(clippy::too_many_arguments)]
fn fwd_forces(
    grid: &RadialGrid,
    col: f64,
    omega: f64,
    omega_r: f64,
    lambda_total: f64,
    lam_c: f64,
    lam_s: f64,
    rho: f64,
    n_b: usize,
    n_psi: usize,
    v_in_hub_x: f64,
    v_in_hub_y: f64,
    theta_1c: f64,
    theta_1s: f64,
    polar: &PolarKind,
) -> (f64, f64, f64, f64) {
    let mut t_acc = 0.0;
    let mut q_acc = 0.0;
    let mut mx_acc = 0.0;
    let mut my_acc = 0.0;
    let n_bf = n_b as f64;
    let inv_n_psi = 1.0 / (n_psi as f64);

    for i_psi in 0..n_psi {
        let psi = 2.0 * PI * (i_psi as f64) * inv_n_psi;
        let cos_psi = psi.cos();
        let sin_psi = psi.sin();
        let v_t_extra = v_in_hub_x * sin_psi + v_in_hub_y * cos_psi;
        let col_psi = col + theta_1c * cos_psi + theta_1s * sin_psi;
        let mut rdt_sum = 0.0;
        for i in 0..grid.r_mid.len() {
            let r = grid.r_mid[i];
            let v_t = omega * r + v_t_extra;
            if v_t <= 0.0 {
                continue;
            }
            let x = grid.x_mid[i];
            let lam_local = lambda_total + x * (lam_c * cos_psi + lam_s * sin_psi);
            let v_a = lam_local * omega_r;
            let phi = v_a.atan2(v_t);
            let alpha = col_psi + grid.twist_rad[i] - phi;
            let (cl, cd) = polar.cl_cd(alpha);
            let cos_p = phi.cos();
            let sin_p = phi.sin();
            let cn = cl * cos_p - cd * sin_p;
            let ct = cl * sin_p + cd * cos_p;
            let q = 0.5 * rho * (v_a * v_a + v_t * v_t) * grid.chord[i] * grid.dr * n_bf;
            let dt = q * cn;
            t_acc += dt;
            q_acc += q * ct * r;
            rdt_sum += r * dt;
        }
        mx_acc += rdt_sum * sin_psi;
        my_acc += rdt_sum * cos_psi;
    }
    (
        t_acc * inv_n_psi,
        q_acc * inv_n_psi,
        mx_acc * inv_n_psi,
        my_acc * inv_n_psi,
    )
}

#[derive(Clone)]
pub struct PittPetersModel {
    pub defn: RotorDefinition,
    pub n_psi_elements: usize,
    pub polar: PolarKind,
    pub grid: RadialGrid,
}

impl PittPetersModel {
    pub fn build(defn: RotorDefinition, n_psi_elements: usize, polar: PolarKind) -> Self {
        let grid = RadialGrid::from_blade(&defn.blade);
        Self {
            defn,
            n_psi_elements,
            polar,
            grid,
        }
    }
}

impl AeroModel for PittPetersModel {
    type State = PittPetersRotorState;

    fn initial_state(&self) -> Self::State {
        PittPetersRotorState::default()
    }

    fn inflow_taus(&self, inputs: &RotorInputs, state: &Self::State) -> Vec<f64> {
        let omega = state.omega_rad_s;
        let r_tip = self.defn.blade.radius_m;
        let omega_r = omega * r_tip;
        let z_hub = Vec3::new(0.0, 0.0, 1.0);
        let hub_axis = inputs.R_hub * z_hub;
        let v_rel = inputs.wind_world - inputs.v_hub_world;
        let v_climb = v_rel.dot(hub_axis);
        let v_edge = (v_rel - hub_axis * v_climb).norm();
        let v0 = state.lambda_0 * omega_r;
        let v_t_disk = (v_edge * v_edge + (v_climb + v0).powi(2))
            .sqrt()
            .max(V_T_HOVER_FLOOR_FRAC * omega_r.max(1.0));
        let tau_0 = (8.0 * r_tip) / (3.0 * PI * v_t_disk);
        let tau_cs = (16.0 * r_tip) / (45.0 * PI * v_t_disk);
        vec![tau_0, tau_cs, tau_cs, f64::INFINITY, f64::INFINITY]
    }

    fn compute_forces(
        &self,
        inputs: &RotorInputs,
        state: &PittPetersRotorState,
    ) -> (AeroResult, PittPetersRotorState) {
        let blade = &self.defn.blade;
        let omega = state.omega_rad_s;
        let rho = inputs.rho_kg_m3;
        let r_tip = blade.radius_m;
        let area = PI * r_tip * r_tip;

        let omega_r = omega * r_tip;
        let z_hub = Vec3::new(0.0, 0.0, 1.0);
        let hub_axis = inputs.R_hub * z_hub;
        let v_rel = inputs.wind_world - inputs.v_hub_world;
        let v_climb = v_rel.dot(hub_axis);
        let v_inplane = v_rel - hub_axis * v_climb;
        let v_edge = v_inplane.norm();
        let mu = v_edge / omega_r.max(EPS_OMEGA_R);
        let v_inplane_hub = inputs.R_hub.transpose() * v_inplane;

        let lam0 = state.lambda_0;
        let lam_c = state.lambda_c;
        let lam_s = state.lambda_s;
        let lambda_climb = if omega_r > EPS_OMEGA_R {
            v_climb / omega_r
        } else {
            0.0
        };
        let lambda_total = lam0 + lambda_climb;

        let gains = self.defn.control_gains();
        let (theta_1c, theta_1s) = cyclic_coeffs(inputs.tilt_lon, inputs.tilt_lat, gains);
        let has_cyclic = theta_1c.abs() + theta_1s.abs() > EPS_DENOM;

        // ------------------------------------------------------------------
        // Blade element forces
        // ------------------------------------------------------------------
        let (mut t_total, mut q_total, mut mx_hub, mut my_hub) = (0.0, 0.0, 0.0, 0.0);
        if omega_r > EPS_OMEGA_R {
            if (mu > 0.01 || has_cyclic || lam_c.abs() + lam_s.abs() > EPS_DENOM) && omega > 1.0 {
                let (t, q, mx, my) = fwd_forces(
                    &self.grid,
                    inputs.collective_rad,
                    omega,
                    omega_r,
                    lambda_total,
                    lam_c,
                    lam_s,
                    rho,
                    blade.n_blades,
                    self.n_psi_elements,
                    v_inplane_hub[0],
                    v_inplane_hub[1],
                    theta_1c,
                    theta_1s,
                    &self.polar,
                );
                t_total = t;
                q_total = q;
                mx_hub = mx;
                my_hub = my;
            } else {
                let (t, q) = axial_forces(
                    &self.grid,
                    inputs.collective_rad,
                    omega,
                    omega_r,
                    lambda_total,
                    rho,
                    blade.n_blades,
                    &self.polar,
                );
                t_total = t;
                q_total = q;
            }
        }

        // ------------------------------------------------------------------
        // Pitt-Peters L-matrix steady-state targets + ODE
        // ------------------------------------------------------------------
        let t_pos = t_total.max(0.0);
        let v_h = if t_pos > EPS_OMEGA_R {
            (t_pos / (2.0 * rho * area)).sqrt()
        } else {
            0.0
        };

        let v0 = lam0 * omega_r;
        let v_t_disk = (v_edge * v_edge + (v_climb + v0).powi(2))
            .sqrt()
            .max(V_T_HOVER_FLOOR_FRAC * omega_r.max(1.0));

        let v_c = (-v_climb).max(0.0);
        let lam2 = if v_h > VRS_DESCENT_THRESHOLD {
            v_c / v_h
        } else {
            0.0
        };

        let mu_t_eff = (if omega_r > EPS_OMEGA_R {
            v_t_disk / omega_r
        } else {
            0.0
        })
        .max(MU_T_FLOOR);
        let mu_inplane = v_edge / omega_r.max(EPS_OMEGA_R);
        // small_atan2_eps avoids chi flipping between +-pi/2 when lambda_total
        // is exactly zero in pure-edgewise flow.
        let chi = mu_inplane.atan2(lambda_total.abs() + EPS_OMEGA_R);
        let cos_chi = chi.cos();
        let tan_half_chi = (0.5 * chi).tan();
        let l_off = (15.0 * PI / 64.0) * tan_half_chi;
        let l_cc = 4.0 * cos_chi / (1.0 + cos_chi);
        let l_ss = 4.0 / (1.0 + cos_chi);

        let norm = rho * area * omega_r * r_tip * v_t_disk;
        let c_l_hub = if norm > EPS_DENOM { mx_hub / norm } else { 0.0 };
        let c_m_hub = if norm > EPS_DENOM { my_hub / norm } else { 0.0 };

        let lam0_ss = if v_climb < -VRS_DESCENT_THRESHOLD && lam2 > 0.0 && lam2 < 2.0 {
            if omega_r > EPS_OMEGA_R {
                vrs_lambda1(lam2) * v_h / omega_r
            } else {
                0.0
            }
        } else if omega_r > EPS_OMEGA_R {
            t_total / (2.0 * rho * area * v_t_disk * omega_r) + l_off * c_m_hub / mu_t_eff
        } else {
            0.0
        };

        let c_t = if omega_r > EPS_OMEGA_R {
            t_total / (rho * area * omega_r * omega_r)
        } else {
            0.0
        };
        let lam_c_ss = (-l_off * c_t + l_cc * c_m_hub) / mu_t_eff;
        let lam_s_ss = (l_ss * c_l_hub) / mu_t_eff;

        let tau_0 = (8.0 * r_tip) / (3.0 * PI * v_t_disk);
        let tau_cs = (16.0 * r_tip) / (45.0 * PI * v_t_disk);

        let d_lam0 = (lam0_ss - lam0) / tau_0;
        let d_lam_c = (lam_c_ss - lam_c) / tau_cs;
        let d_lam_s = (lam_s_ss - lam_s) / tau_cs;

        // Mechanical
        let i_ode = self.defn.autorotation.I_ode_kgm2.unwrap_or(1.0);
        let d_omega = (-q_total + inputs.motor_torque_Nm) / i_ode;
        let d_spin_angle = omega;

        // Outputs
        let f_world = hub_axis * (-t_total);
        let mxyz_hub = Vec3::new(mx_hub, my_hub, 0.0);
        let m_orbital = inputs.R_hub * mxyz_hub;
        let m_spin = hub_axis * q_total;

        let result = AeroResult {
            F_world: f_world,
            M_orbital: m_orbital,
            Q_spin: q_total,
            M_spin: m_spin,
        };
        let derivative = PittPetersRotorState {
            lambda_0: d_lam0,
            lambda_c: d_lam_c,
            lambda_s: d_lam_s,
            omega_rad_s: d_omega,
            spin_angle_rad: d_spin_angle,
        };
        (result, derivative)
    }
}
