// Level 2: Oye 2-stage annular dynamic inflow.
// See ../CLAUDE.md "Oye 2-stage annular dynamic inflow" and dynbem/oye.py.

use std::f64::consts::PI;

use crate::aero_io::{AeroResult, RotorInputs};
use crate::aero_model::AeroModel;
use crate::bem_common::{
    assemble_result, element_force, kinematics, v_t_disk, vrs_regime, ElementCtx, PsiKernel,
    RadialGrid, SweepCtx,
};
use crate::common::{vrs_lambda1, EPS_OMEGA_R, VRS_DESCENT_THRESHOLD, V_T_HOVER_FLOOR_FRAC};
use crate::cyclic::cyclic_coeffs;
use crate::polar::PolarKind;
use crate::rotor_definition::RotorDefinition;
use crate::rotor_state::OyeRotorState;

/// Oye empirical coupling constant (original Oye 1990, OpenFAST default).
/// Damps the coupling between the quasi-steady target W_qs and the
/// intermediate filter state W_int.
pub const OYE_K: f64 = 0.6;

/// Upper clamp on rotor-mean induction factor `a` in the tau1 formula.
/// Without this, tau1 -> infinity as a -> 1/1.3 ~= 0.77 and the filter
/// stalls. 0.5 is the actuator-disk limit; matches OpenFAST DBEMT_Mod=1.
const A_AVG_CLAMP: f64 = 0.5;

/// Sanity clamp on quasi-steady inflow per annulus. Physical W is small
/// (~0.01-0.1); this catches transient blow-ups from a noisy dT estimate
/// without poisoning later steps.
const W_QS_CLAMP: f64 = 1.0;

/// Oye psi-loop kernel.
///
/// Local inflow is per-annulus and azimuth-independent
/// (`lambda_climb + W[i]`), and we accumulate the azimuth-averaged dT
/// into `dt_avg[i]` for the W_qs solver. Both happen in one `element`
/// override -- splitting them into separate trait methods would force
/// `dt` to be recomputed or carried by side state, with no benefit.
struct OyeKernel<'a> {
    lambda_climb: f64,
    w: &'a [f64],
    dt_avg: &'a mut [f64],
}

impl<'a> PsiKernel for OyeKernel<'a> {
    #[inline(always)]
    fn element(&mut self, sweep: &SweepCtx<'_>, ctx: &ElementCtx) -> (f64, f64) {
        let lam = self.lambda_climb + self.w[ctx.i];
        let v_a = lam * sweep.omega_r;
        let (dt, dq) = element_force(v_a, sweep, ctx);
        self.dt_avg[ctx.i] += dt / (sweep.n_psi as f64);
        (dt, dq)
    }
}

/// Quasi-steady annulus inflow from Glauert mass-flow momentum balance.
///
///     dT = 4*pi*r*dr*rho*V_resultant*v_i
///
/// Linearised in W_qs using a rotor-mean mu_T (computed externally from the
/// converged state). Stays stable in forward flight where the pure
/// axial-momentum form blows up at low lambda_r (autorotation / VRS).
fn solve_w_qs(
    dt_avg: &[f64],
    x_arr: &[f64],
    _dr: f64,
    r_tip: f64,
    omega_r: f64,
    mu_t: f64,
    rho: f64,
) -> Vec<f64> {
    if omega_r < EPS_OMEGA_R {
        return vec![0.0; dt_avg.len()];
    }
    let area = PI * r_tip * r_tip;
    // rho_norm here is (rho * A * Omega_R^2 * dr / R) per Python implementation.
    let rho_norm = (rho * area * omega_r * omega_r * _dr / r_tip).max(1e-30);
    let mu_t_safe = mu_t.max(crate::common::MU_T_FLOOR);
    let mut out = Vec::with_capacity(dt_avg.len());
    for i in 0..dt_avg.len() {
        let d_cdx = dt_avg[i] / rho_norm;
        let x = x_arr[i].max(EPS_OMEGA_R);
        let w = d_cdx / (4.0 * x * mu_t_safe);
        out.push(w.clamp(-W_QS_CLAMP, W_QS_CLAMP));
    }
    out
}

/// Oye time constants per annulus.
///
///     tau1     = 1.1 / (1 - 1.3*min(a, 0.5)) * R / V_inf      (rotor-mean)
///     tau2(r)  = (0.39 - 0.26*(r/R)^2) * tau1                 (radius-dependent)
///
/// The 0.5 clamp on a keeps tau1 finite at the actuator-disk limit.
fn oye_taus(r_tip: f64, x_arr: &[f64], v_inf: f64, a_avg: f64) -> (f64, Vec<f64>) {
    let a_c = a_avg.clamp(0.0, A_AVG_CLAMP);
    let tau1 = 1.1 / (1.0 - 1.3 * a_c) * r_tip / v_inf.max(V_T_HOVER_FLOOR_FRAC);
    let tau2: Vec<f64> = x_arr
        .iter()
        .map(|&x| (0.39 - 0.26 * x * x) * tau1)
        .collect();
    (tau1, tau2)
}

#[derive(Clone)]
pub struct OyeBEMModel {
    pub defn: RotorDefinition,
    pub n_psi_elements: usize,
    pub coupling_k: f64,
    pub polar: PolarKind,
    pub grid: RadialGrid,
}

impl AeroModel for OyeBEMModel {
    type State = OyeRotorState;

    fn initial_state(&self) -> Self::State {
        let n = self.defn.blade.n_elements;
        OyeRotorState {
            W_int: vec![0.0; n],
            W: vec![0.0; n],
            omega_rad_s: 0.0,
            spin_angle_rad: 0.0,
        }
    }

    fn inflow_taus(&self, inputs: &RotorInputs, state: &Self::State) -> Vec<f64> {
        let r_tip = self.defn.blade.radius_m;
        let n = self.defn.blade.n_elements;
        let kin = kinematics(inputs, state.omega_rad_s, r_tip);
        let omega_r = kin.omega_r;
        if omega_r < EPS_OMEGA_R {
            return vec![f64::INFINITY; 2 * n + 2];
        }
        let v_climb = kin.v_climb;
        let v_edge = kin.v_edge;
        let w_mean: f64 = if n > 0 {
            state.W.iter().sum::<f64>() / (n as f64)
        } else {
            0.0
        };
        let v_inf = v_t_disk(v_edge, v_climb, w_mean * omega_r, omega_r);
        let a_avg = if v_inf > VRS_DESCENT_THRESHOLD {
            w_mean * omega_r / v_inf
        } else {
            0.0
        };
        let (tau1, tau2) = oye_taus(r_tip, &self.grid.x_mid, v_inf, a_avg);
        let mut out = Vec::with_capacity(2 * n + 2);
        out.extend(std::iter::repeat(tau1).take(n));
        out.extend(tau2);
        out.push(f64::INFINITY);
        out.push(f64::INFINITY);
        out
    }

    fn compute_forces(
        &self,
        inputs: &RotorInputs,
        state: &OyeRotorState,
    ) -> (AeroResult, OyeRotorState) {
        let blade = &self.defn.blade;
        let omega = state.omega_rad_s;
        let rho = inputs.rho_kg_m3;
        let r_tip = blade.radius_m;
        let n = blade.n_elements;

        let kin = kinematics(inputs, omega, r_tip);
        let omega_r = kin.omega_r;
        let hub_axis = kin.hub_axis;
        let v_climb = kin.v_climb;
        let v_edge = kin.v_edge;
        let v_inplane_hub = kin.v_inplane_hub;

        let gains = self.defn.control_gains();
        let (theta_1c, theta_1s) = cyclic_coeffs(inputs.tilt_lon, inputs.tilt_lat, gains);

        // Blade element forces -- psi-loop with per-annulus W.
        let lambda_climb = if omega_r > EPS_OMEGA_R {
            v_climb / omega_r
        } else {
            0.0
        };
        let mut dt_avg = vec![0.0f64; n];
        let (t_total, q_total, mx_hub, my_hub) = if omega_r > EPS_OMEGA_R && omega > 1.0 {
            let mut kernel = OyeKernel {
                lambda_climb,
                w: &state.W,
                dt_avg: &mut dt_avg,
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
                v_in_hub_x: v_inplane_hub[0],
                v_in_hub_y: v_inplane_hub[1],
                theta_1c,
                theta_1s,
            };
            sweep.run(&mut kernel)
        } else {
            (0.0, 0.0, 0.0, 0.0)
        };

        // W_qs per annulus, with VRS uniform override.
        let w_mean: f64 = if n > 0 {
            state.W.iter().sum::<f64>() / (n as f64)
        } else {
            0.0
        };
        let v0_mean = w_mean * omega_r;
        let vt_disk = v_t_disk(v_edge, v_climb, v0_mean, omega_r);
        let mu_t = vt_disk / omega_r.max(EPS_OMEGA_R);

        let area = PI * r_tip * r_tip;
        let vrs = vrs_regime(t_total, v_climb, rho, area);

        let w_qs: Vec<f64> = if vrs.in_vrs && omega_r > EPS_OMEGA_R {
            let w_uniform = vrs_lambda1(vrs.lam2) * vrs.v_h / omega_r;
            vec![w_uniform; n]
        } else {
            solve_w_qs(
                &dt_avg,
                &self.grid.x_mid,
                self.grid.dr,
                r_tip,
                omega_r,
                mu_t,
                rho,
            )
        };

        // Oye 2-stage filter ODE (Mod=1: dW_qs/dt = 0 across the outer step).
        let a_avg = if vt_disk > VRS_DESCENT_THRESHOLD {
            w_mean * omega_r / vt_disk
        } else {
            0.0
        };
        let (_tau1_scalar, tau2_arr) = oye_taus(r_tip, &self.grid.x_mid, vt_disk, a_avg);
        let tau1_scalar = _tau1_scalar;

        let mut d_w_int = Vec::with_capacity(n);
        let mut d_w = Vec::with_capacity(n);
        for i in 0..n {
            d_w_int.push((w_qs[i] - state.W_int[i]) / tau1_scalar);
            d_w.push((state.W_int[i] - state.W[i]) / tau2_arr[i]);
        }

        let i_ode = self.defn.autorotation.I_ode_kgm2.unwrap_or(1.0);
        let d_omega = (-q_total + inputs.motor_torque_Nm) / i_ode;
        let d_psi = omega;

        let result = assemble_result(t_total, q_total, mx_hub, my_hub, hub_axis, &inputs.R_hub);
        let derivative = OyeRotorState {
            W_int: d_w_int,
            W: d_w,
            omega_rad_s: d_omega,
            spin_angle_rad: d_psi,
        };
        (result, derivative)
    }
}
