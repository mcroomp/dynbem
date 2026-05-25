// Level 2: Oye 2-stage annular dynamic inflow.
// See ../CLAUDE.md "Oye 2-stage annular dynamic inflow" and dynbem/oye.py.

use std::f64::consts::PI;

use crate::aero_io::{AeroResult, RotorInputs};
use crate::aero_model::AeroModel;
use crate::bem_common::{
    assemble_result, element_force, kinematics, v_t_disk, vrs_regime, ElementCtx, PsiKernel,
    RadialGrid, SweepCtx,
};
use crate::common::{
    vrs_lambda1, EPS_OMEGA_R, MIN_LOSS_FACTOR, VRS_DESCENT_THRESHOLD, V_T_HOVER_FLOOR_FRAC,
};
use crate::cyclic::cyclic_coeffs;
use crate::polar::Polar;
use crate::quasi_static_bem::{prandtl_hub_loss, prandtl_tip_loss};
use crate::rotor_definition::RotorDefinition;
use crate::rotor_state::OyeRotorState;

/// Oye empirical coupling constant (original Oye 1990, OpenFAST default).
/// Damps the coupling between the quasi-steady target W_qs and the
/// intermediate filter state W_int.
pub const OYE_K: f64 = 0.6;

/// Oye 2-stage filter step.
#[inline(always)]
pub fn oye_filter_step(
    d_w_int: &mut [f64],
    d_w: &mut [f64],
    w_qs: &[f64],
    w_int: &[f64],
    w: &[f64],
    tau2_arr: &[f64],
    tau1_scalar: f64,
    n_active: usize,
) {
    for i in 0..n_active {
        d_w_int[i] = (w_qs[i] - w_int[i]) / tau1_scalar;
        d_w[i] = (w_int[i] - w[i]) / tau2_arr[i];
    }
}

/// Upper clamp on rotor-mean induction factor `a` in the tau1 formula.
/// Without this, tau1 -> infinity as a -> 1/1.3 ~= 0.77 and the filter
/// stalls. 0.5 is the actuator-disk limit; matches OpenFAST DBEMT_Mod=1.
const A_AVG_CLAMP: f64 = 0.5;

/// Sanity clamp on quasi-steady inflow per annulus. Physical W is small
/// (~0.01-0.1); this catches transient blow-ups from a noisy dT estimate
/// without poisoning later steps.
const W_QS_CLAMP: f64 = 1.0;

/// Oye psi-loop kernel.
struct OyeKernel<'a> {
    lambda_climb: f64,
    w: &'a [f64],
    dt_avg: &'a mut [f64],
}

impl<'a> PsiKernel for OyeKernel<'a> {
    #[inline(always)]
    fn element<P: Polar>(&mut self, sweep: &SweepCtx<'_, P>, ctx: &ElementCtx) -> (f64, f64) {
        let lam = self.lambda_climb + self.w[ctx.i];
        let v_a = lam * sweep.omega_r;
        let (dt, dq) = element_force(v_a, sweep, ctx);
        self.dt_avg[ctx.i] += dt * sweep.n_psi_inv;
        (dt, dq)
    }
}

fn solve_w_qs(
    dt_avg: &[f64],
    x_arr: &[f64],
    f_loss: &[f64],
    n_active: usize,
    _dr: f64,
    r_tip: f64,
    omega_r: f64,
    mu_t: f64,
    rho: f64,
) -> Vec<f64> {
    if omega_r < EPS_OMEGA_R {
        return vec![0.0; n_active];
    }
    let area = PI * r_tip * r_tip;
    let rho_norm = (rho * area * omega_r * omega_r * _dr / r_tip).max(1e-30);
    let mu_t_safe = mu_t.max(crate::common::MU_T_FLOOR);
    let mut out = vec![0.0; n_active];
    for i in 0..n_active {
        let d_cdx = dt_avg[i] / rho_norm;
        let x = x_arr[i].max(EPS_OMEGA_R);
        let f = f_loss[i].max(MIN_LOSS_FACTOR);
        let w = d_cdx / (4.0 * x * f * mu_t_safe);
        out[i] = w.clamp(-W_QS_CLAMP, W_QS_CLAMP);
    }
    out
}

fn oye_taus(r_tip: f64, x_arr: &[f64], v_inf: f64, a_avg: f64, n_active: usize) -> (f64, Vec<f64>) {
    let a_c = a_avg.clamp(0.0, A_AVG_CLAMP);
    let tau1 = 1.1 / (1.0 - 1.3 * a_c) * r_tip / v_inf.max(V_T_HOVER_FLOOR_FRAC);
    let mut tau2 = vec![0.0; n_active];
    for i in 0..n_active {
        let x = x_arr[i];
        tau2[i] = (0.39 - 0.26 * x * x) * tau1;
    }
    (tau1, tau2)
}

#[derive(Clone)]
pub struct OyeBEMModel<P: Polar> {
    pub defn: RotorDefinition,
    pub n_psi_elements: usize,
    pub coupling_k: f64,
    pub polar: P,
    pub grid: RadialGrid,
}

impl<P: Polar + Clone> OyeBEMModel<P> {
    pub fn build(defn: RotorDefinition, n_psi_elements: usize, polar: P) -> Self {
        Self::build_with_k(defn, n_psi_elements, polar, OYE_K)
    }

    pub fn build_with_k(
        defn: RotorDefinition,
        n_psi_elements: usize,
        polar: P,
        coupling_k: f64,
    ) -> Self {
        let grid = RadialGrid::from_blade(&defn.blade);
        Self { defn, n_psi_elements, coupling_k, polar, grid }
    }
}

impl<P: Polar + Clone> AeroModel for OyeBEMModel<P> {
    type State = OyeRotorState;

    fn initial_state(&self) -> Self::State {
        let n = self.defn.blade.n_elements;
        OyeRotorState::zeros(n)
    }

    fn inflow_taus(&self, inputs: &RotorInputs, state: &Self::State) -> Vec<f64> {
        let r_tip = self.defn.blade.radius_m;
        let n = self.defn.blade.n_elements;
        let kin = kinematics(inputs, inputs.omega_rad_s, r_tip);
        let omega_r = kin.omega_r;
        if omega_r < EPS_OMEGA_R {
            return vec![f64::INFINITY; 2 * n];
        }
        let v_climb = kin.v_climb;
        let v_edge = kin.v_edge;
        let w_mean: f64 = if n > 0 {
            state.w_slice().iter().sum::<f64>() / (n as f64)
        } else {
            0.0
        };
        let v_inf = v_t_disk(v_edge, v_climb, w_mean * omega_r, omega_r);
        let a_avg = if v_inf > VRS_DESCENT_THRESHOLD {
            w_mean * omega_r / v_inf
        } else {
            0.0
        };
        let (tau1, tau2) = oye_taus(r_tip, &self.grid.x_mid[..n], v_inf, a_avg, n);
        let mut out = Vec::with_capacity(2 * n);
        out.extend(std::iter::repeat(tau1).take(n));
        out.extend(tau2);
        out
    }

    fn compute_forces(
        &self,
        inputs: &RotorInputs,
        state: &OyeRotorState,
    ) -> (AeroResult, OyeRotorState) {
        let blade = &self.defn.blade;
        let omega = inputs.omega_rad_s;
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

        let lambda_climb = if omega_r > EPS_OMEGA_R {
            v_climb / omega_r
        } else {
            0.0
        };
        let mut dt_avg = vec![0.0; n];
        let (t_total, q_total, mx_hub, my_hub) = if omega_r > EPS_OMEGA_R && omega > 1.0 {
            let mut kernel = OyeKernel {
                lambda_climb,
                w: state.w_slice(),
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
                n_psi_inv: 1.0 / (self.n_psi_elements as f64),
                v_in_hub_x: v_inplane_hub[0],
                v_in_hub_y: v_inplane_hub[1],
                theta_1c,
                theta_1s,
            };
            sweep.run(&mut kernel)
        } else {
            (0.0, 0.0, 0.0, 0.0)
        };

        let w_mean: f64 = if n > 0 {
            state.w_slice().iter().sum::<f64>() / (n as f64)
        } else {
            0.0
        };
        let v0_mean = w_mean * omega_r;
        let vt_disk = v_t_disk(v_edge, v_climb, v0_mean, omega_r);
        let mu_t = vt_disk / omega_r.max(EPS_OMEGA_R);

        let area = PI * r_tip * r_tip;
        let vrs = vrs_regime(t_total, v_climb, rho, area);

        let use_tip_loss = self.defn.blade.tip_loss;
        let mut f_loss = vec![1.0; n];
        if use_tip_loss && omega_r > EPS_OMEGA_R {
            let n_b = self.defn.blade.n_blades;
            let x_hub = self.grid.x_hub;
            for i in 0..n {
                let r = self.grid.r_mid[i];
                let lam_local = lambda_climb + state.W[i];
                let v_a = lam_local * omega_r;
                let v_t = omega * r;
                let phi = v_a.atan2(v_t);
                let x = self.grid.x_mid[i];
                f_loss[i] = (prandtl_tip_loss(n_b, x, phi) * prandtl_hub_loss(n_b, x, x_hub, phi))
                    .max(MIN_LOSS_FACTOR);
            }
        }

        let w_qs: Vec<f64> = if vrs.in_vrs && omega_r > EPS_OMEGA_R {
            let w_uniform = vrs_lambda1(vrs.lam2) * vrs.v_h / omega_r;
            vec![w_uniform; n]
        } else {
            solve_w_qs(
                &dt_avg,
                &self.grid.x_mid[..n],
                &f_loss,
                n,
                self.grid.dr,
                r_tip,
                omega_r,
                mu_t,
                rho,
            )
        };

        let a_avg = if vt_disk > VRS_DESCENT_THRESHOLD {
            w_mean * omega_r / vt_disk
        } else {
            0.0
        };
        let (tau1_scalar, tau2_arr) = oye_taus(r_tip, &self.grid.x_mid[..n], vt_disk, a_avg, n);

        let mut d_w_int = vec![0.0; n];
        let mut d_w = vec![0.0; n];
        oye_filter_step(
            &mut d_w_int,
            &mut d_w,
            &w_qs,
            &state.W_int,
            &state.W,
            &tau2_arr,
            tau1_scalar,
            n,
        );

        let result = assemble_result(t_total, q_total, mx_hub, my_hub, hub_axis, &inputs.R_hub);
        let derivative = OyeRotorState {
            n_elements: n,
            W_int: d_w_int,
            W: d_w,
        };
        (result, derivative)
    }
}
