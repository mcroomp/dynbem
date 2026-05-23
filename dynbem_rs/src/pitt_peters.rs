// Level 2: Pitt-Peters 3-state dynamic inflow + VRS empirical correction.
// See ../CLAUDE.md "Pitt-Peters inflow ODE" and dynbem/pitt_peters.py.

use std::f64::consts::PI;

use crate::aero_io::{AeroResult, RotorInputs};
use crate::aero_model::AeroModel;
use crate::bem_common::{
    assemble_result, element_force, kinematics, v_t_disk, vrs_regime, ElementCtx, PsiKernel,
    RadialGrid, SweepCtx,
};
use crate::common::{vrs_lambda1, EPS_DENOM, EPS_OMEGA_R, MAX_BEM_ELEMENTS, MU_T_FLOOR};
use crate::cyclic::cyclic_coeffs;
use crate::polar::PolarKind;
use crate::rotor_definition::RotorDefinition;
use crate::rotor_state::PittPetersRotorState;

/// Sum thrust and torque over radial elements in axial flight (mu = 0).
/// Fast path bypassing the full psi-loop -- builds a single SweepCtx
/// (n_psi=1, zero cyclic, zero in-plane wind) and calls element_force
/// directly with the prescribed v_a = lambda_total * omega_r.
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
    let sweep = SweepCtx {
        grid,
        polar,
        col,
        omega,
        omega_r,
        rho,
        n_b,
        n_psi: 1,
        n_psi_inv: 1.0,
        v_in_hub_x: 0.0,
        v_in_hub_y: 0.0,
        theta_1c: 0.0,
        theta_1s: 0.0,
    };
    let mut t_acc = 0.0;
    let mut q_acc = 0.0;
    let v_a = lambda_total * omega_r;

    let n = grid.n_elements;
    assert!(n < MAX_BEM_ELEMENTS);

    let r_mid = &grid.r_mid;
    let chord = &grid.chord;
    let twist = &grid.twist_rad;
    for i in 0..n {
        let r = r_mid[i];
        let v_t = omega * r;
        if v_t < EPS_DENOM {
            continue;
        }
        let ctx = ElementCtx {
            i,
            cos_psi: 1.0,
            sin_psi: 0.0,
            r,
            chord: chord[i],
            twist: twist[i],
            dr: grid.dr,
            col_psi: col,
            v_t,
        };
        let (dt, dq) = element_force(v_a, &sweep, &ctx);
        t_acc += dt;
        q_acc += dq;
    }
    (t_acc, q_acc)
}

/// Pitt-Peters psi-loop kernel.
///
/// Local inflow expands the three harmonic states (lambda_0 in
/// `lambda_total`, plus lam_c and lam_s) at element i and azimuth psi:
///     lam_local = lambda_total + x*(lam_c*cos psi + lam_s*sin psi).
/// No per-element callback -- PP doesn't need the azimuth-averaged dT/dx.
struct PpKernel<'a> {
    lambda_total: f64,
    lam_c: f64,
    lam_s: f64,
    x_mid: &'a [f64],
}

impl<'a> PsiKernel for PpKernel<'a> {
    #[inline(always)]
    fn lam_local(&self, i: usize, cos_psi: f64, sin_psi: f64) -> f64 {
        self.lambda_total + self.x_mid[i] * (self.lam_c * cos_psi + self.lam_s * sin_psi)
    }
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
        let r_tip = self.defn.blade.radius_m;
        let kin = kinematics(inputs, state.omega_rad_s, r_tip);
        let v0 = state.lambda_0 * kin.omega_r;
        let vt = v_t_disk(kin.v_edge, kin.v_climb, v0, kin.omega_r);
        let tau_0 = (8.0 * r_tip) / (3.0 * PI * vt);
        let tau_cs = (16.0 * r_tip) / (45.0 * PI * vt);
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

        let kin = kinematics(inputs, omega, r_tip);
        let omega_r = kin.omega_r;
        let hub_axis = kin.hub_axis;
        let v_climb = kin.v_climb;
        let v_edge = kin.v_edge;
        let mu = kin.mu;
        let v_inplane_hub = kin.v_inplane_hub;

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
                let mut kernel = PpKernel {
                    lambda_total,
                    lam_c,
                    lam_s,
                    x_mid: &self.grid.x_mid[..self.grid.n_elements],
                };
                let sweep = SweepCtx {
                    grid: &self.grid,
                    polar: &self.polar,
                    col: inputs.collective_rad,
                    omega,
                    omega_r,
                    rho,
                    n_b: blade.n_blades,
                    n_psi: self.n_psi_elements,
                    n_psi_inv: 1.0 / (self.n_psi_elements as f64),
                    v_in_hub_x: v_inplane_hub[0],
                    v_in_hub_y: v_inplane_hub[1],
                    theta_1c,
                    theta_1s,
                };
                let (t, q, mx, my) = sweep.run(&mut kernel);
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
        let vrs = vrs_regime(t_total, v_climb, rho, area);

        let v0 = lam0 * omega_r;
        let vt_disk = v_t_disk(v_edge, v_climb, v0, omega_r);

        let mu_t_eff = (if omega_r > EPS_OMEGA_R {
            vt_disk / omega_r
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

        let norm = rho * area * omega_r * r_tip * vt_disk;
        let c_l_hub = if norm > EPS_DENOM { mx_hub / norm } else { 0.0 };
        let c_m_hub = if norm > EPS_DENOM { my_hub / norm } else { 0.0 };

        let lam0_ss = if vrs.in_vrs {
            if omega_r > EPS_OMEGA_R {
                vrs_lambda1(vrs.lam2) * vrs.v_h / omega_r
            } else {
                0.0
            }
        } else if omega_r > EPS_OMEGA_R {
            t_total / (2.0 * rho * area * vt_disk * omega_r) + l_off * c_m_hub / mu_t_eff
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

        let tau_0 = (8.0 * r_tip) / (3.0 * PI * vt_disk);
        let tau_cs = (16.0 * r_tip) / (45.0 * PI * vt_disk);

        let d_lam0 = (lam0_ss - lam0) / tau_0;
        let d_lam_c = (lam_c_ss - lam_c) / tau_cs;
        let d_lam_s = (lam_s_ss - lam_s) / tau_cs;

        // Mechanical
        let i_ode = self.defn.autorotation.I_ode_kgm2.unwrap_or(1.0);
        let d_omega = (-q_total + inputs.motor_torque_Nm) / i_ode;
        let d_spin_angle = omega;

        // Outputs
        let result = assemble_result(t_total, q_total, mx_hub, my_hub, hub_axis, &inputs.R_hub);
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
