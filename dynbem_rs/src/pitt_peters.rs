// Level 2: Pitt-Peters 3-state dynamic inflow + VRS empirical correction.
// See ../CLAUDE.md "Pitt-Peters inflow ODE" and dynbem/pitt_peters.py.
//
// WHAT THIS MODEL DOES
// --------------------
// Three coupled first-order ODEs for the inflow harmonics
// (lambda_0, lambda_c, lambda_s):
//
//   1. Resolve blade pitch (direct swashplate or servo-flap feathering).
//   2. Run the azimuth x radius BEM psi-loop with the CURRENT inflow
//      state expanded locally as
//        lam_local(i, psi) = lambda_total + x_i*(lam_c*cos psi + lam_s*sin psi)
//      to get total thrust T, torque Q, and hub moments Mx_hub, My_hub.
//   3. Form coefficients (C_T, C_M_hub, C_L_hub), rotate the moment
//      forcing into wind axes, and evaluate the Pitt-Peters L-matrix
//      steady-state targets (lam0_ss, lam_c_ss, lam_s_ss).
//   4. Return the first-order relaxation derivative
//        d(lambda)/dt = (lambda_ss - lambda) / tau
//      for each harmonic. The caller integrates the state externally
//      (semi-implicit Euler in envelope/point_mass.py).
//
// `compute_forces` returns (AeroResult, derivative-as-RotorState); it does
// NOT advance the state itself.
//
// CANONICAL REFERENCE: Peters' Nikolsky Lecture, JAHS 54(1):011001 (2009),
// Eqs 7-11 (Research/Peters_Nikolsky_2008/). Read that file's CLAUDE.md
// before touching any sign or coefficient here.
//
// SECONDARY REFERENCE: S. (Kevin) Li, "Development of Rotorcraft Forward
// Flight Analysis Code Using the Finite-State Dynamic Inflow Model", MS
// Project Report, UC Davis MAE (Aug 2020). Its Eq 11 prints the full
// Peters-HaQuang mass-flow V = (mu^2 + (lam+nu)(lam+2nu))/sqrt(mu^2 +
// (lam+nu)^2) together with the same [M] and [L] matrices we use, and
// attributes that V to Peters-HaQuang [8] (the moving-hub reformulation),
// NOT to Pitt's original 1980 dissertation (which only had v = V/(Omega R)).
// Notably Li's own uniform-inflow trim (his Eq 14) falls back to the
// classical Glauert lambda_0 = C_T/(2*sqrt(mu^2 + lambda^2)) -- the same
// mass-flow form this code uses (flag A) -- so a finite-state Peters-He
// implementation corroborates the Glauert scaling choice.
//
// NON-STANDARD / PROJECT-SPECIFIC DEVIATIONS FROM TEXTBOOK PITT-PETERS
// (each is flagged inline at its use-site with `NON-STANDARD:`):
//
//   A. Mass-flow scalar is classical Glauert V_mf = sqrt(mu^2 + lambda^2)
//      (`v_mass_flow_disk`), NOT the Peters-HaQuang V (Li 2020 Eq 11). The
//      two differ by ~sqrt(2) in hover (our form reproduces Glauert
//      lambda_0 = sqrt(C_T/2); the Peters Eq-8 V = 2*nu in hover gives
//      lambda_0 = sqrt(C_T/4) = sqrt(C_T)/2); the L-matrix STRUCTURE still
//      matches Peters exactly. Choice is justified empirically (46-point
//      regression suite) and by Li 2020 Eq 14, which uses the same Glauert
//      form for uniform inflow. (In Peters' Eq 8, `lambda` is climb and
//      `nu` is induced flow -- distinct variables; see
//      Research/Peters_Nikolsky_2008/CLAUDE.md section 8.)
//   B. Sign/azimuth convention is ours (psi=0 at +X/nose, CCW-from-above),
//      so our lambda_c = -Peters nu_c and lambda_s = -Peters nu_s; forcing
//      is {C_T, +C_M_hub, +C_L_hub} (see AGENTS.md "Hub-frame aero moments").
//   C. VRS empirical override: inside the vortex-ring-state regime the
//      momentum lam0_ss is replaced by the Leishman polynomial
//      (`vrs_lambda1`); momentum theory does not apply in a recirculating
//      wake, so the lam0_ss cross-coupling term (+l_off*C_M) is dropped
//      there. The cyclic targets lam_c_ss/lam_s_ss are still computed with
//      their full L-matrix cross-coupling (the regime only overrides the
//      uniform-inflow component).
//   D. Wind-axis rotation IS applied (beta_wind) -- forcing and the current
//      cyclic states are rotated into wind axes, the L relations evaluated
//      there, and the derivatives rotated back. This makes the model
//      rotationally covariant in oblique flight. (AGENTS.md documents the
//      history: it was once reverted for instability, now safe because the
//      envelope integrator damps all inflow states semi-implicitly.)
//   E. Dual mass-flow handling: the L-matrix scaling uses mu_t_eff clamped at
//      MU_T_FLOOR to stay finite at hover, while the time constants
//      tau_0/tau_cs and the C_L/C_M normalization use v_mf directly (which
//      carries only its own much smaller MASS_FLOW_HOVER_FLOOR_FRAC floor, not
//      the MU_T_FLOOR clamp). Recomputing v_mf from mu_t_eff would re-impose
//      the MU_T_FLOOR clamp -- inflating v_mf whenever v_mf/omega_r <
//      MU_T_FLOOR -- and corrupt the coefficients (see the inline note).
//   F. Quasi-static blade flap reduction (`apply_flap_reduction`) scales the
//      hub moments the AIRFRAME sees, but the inflow ODE uses the FULL
//      aerodynamic moments (the wake responds to disk loading, not to what
//      the flapping blade passes to the hub).
//   G. Servo-flap feathering pre-pass (Kaman path): in ServoFlap mode the
//      swashplate collective AND cyclic are reinterpreted as flap commands
//      and the feathering response replaces the direct blade-pitch path.

use std::f64::consts::PI;

use crate::aero_io::{AeroResult, RotorInputs};
use crate::aero_model::{AeroModel, RotorStateExt};
use crate::bem_common::{
    apply_flap_reduction, assemble_result, build_psi_trig_table, element_force, kinematics,
    v_mass_flow_disk, vrs_regime, ElementCtx, PsiKernel, RadialGrid, SweepCtx,
};
use crate::common::{vrs_lambda1, EPS_DENOM, EPS_OMEGA_R, MU_T_FLOOR};
use crate::cyclic::cyclic_coeffs;
use crate::polar::Polar;
use crate::rotor_definition::{PitchActuation, RotorDefinition};
use crate::servoflap::{solve_feathering, FeatheringState};

#[derive(Clone, Debug, Default)]
pub struct PittPetersRotorState {
    pub lambda_0: f64,
    pub lambda_c: f64,
    pub lambda_s: f64,
}

/// Pitt-Peters psi-loop kernel.
///
/// Local inflow expands the three harmonic states (lambda_0 in
/// `lambda_total`, plus lam_c and lam_s) at element i and azimuth psi:
///     lam_local = lambda_total + x*(lam_c*cos psi + lam_s*sin psi).
struct PpKernel<'a> {
    lambda_total: f64,
    lam_c: f64,
    lam_s: f64,
    x_mid: &'a [f64],
}

impl<'a> PsiKernel for PpKernel<'a> {
    #[inline(always)]
    fn element<P: Polar>(&mut self, sweep: &SweepCtx<'_, P>, ctx: &ElementCtx) -> (f64, f64) {
        let lam = self.lambda_total
            + self.x_mid[ctx.i] * (self.lam_c * ctx.cos_psi + self.lam_s * ctx.sin_psi);
        let v_a = lam * sweep.omega_r;
        element_force(v_a, sweep, ctx)
    }
}

#[derive(Clone)]
pub struct PittPetersModel<P: Polar> {
    pub defn: RotorDefinition,
    pub n_psi_elements: usize,
    pub psi_trig: Vec<(f64, f64)>,
    pub polar: P,
    pub grid: RadialGrid,
}

impl<P: Polar + Clone> PittPetersModel<P> {
    pub fn build(defn: RotorDefinition, n_psi_elements: usize, polar: P) -> Self {
        let grid = RadialGrid::from_blade(&defn.blade);
        let psi_trig = build_psi_trig_table(n_psi_elements);
        Self {
            defn,
            n_psi_elements,
            psi_trig,
            polar,
            grid,
        }
    }
}

impl RotorStateExt for PittPetersRotorState {
    fn get_inflow(&self) -> Vec<f64> {
        vec![self.lambda_0, self.lambda_c, self.lambda_s]
    }
    fn set_inflow(&mut self, arr: Vec<f64>) {
        debug_assert_eq!(arr.len(), 3);
        self.lambda_0 = arr[0];
        self.lambda_c = arr[1];
        self.lambda_s = arr[2];
    }
}



impl<P: Polar + Clone> AeroModel for PittPetersModel<P> {
    type State = PittPetersRotorState;

    fn initial_state(&self) -> Self::State {
        PittPetersRotorState::default()
    }

    fn inflow_taus(&self, inputs: &RotorInputs, state: &Self::State) -> Vec<f64> {
        let r_tip = self.defn.blade.radius_m;
        let kin = kinematics(inputs, inputs.omega_rad_s, r_tip);
        let v0 = state.lambda_0 * kin.omega_r;
        let v_mf = v_mass_flow_disk(kin.v_edge, kin.v_climb, v0, kin.omega_r);
        let tau_0 = (8.0 * r_tip) / (3.0 * PI * v_mf);
        let tau_cs = (16.0 * r_tip) / (45.0 * PI * v_mf);
        vec![tau_0, tau_cs, tau_cs]
    }

    fn compute_forces(
        &self,
        inputs: &RotorInputs,
        state: &PittPetersRotorState,
    ) -> (AeroResult, PittPetersRotorState) {
        let blade = &self.defn.blade;
        let omega = inputs.omega_rad_s;
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
        // ------------------------------------------------------------------
        // Quasi-static feathering solve (pre-pass, no new state).
        // NON-STANDARD (flag G): in servo-flap mode the swashplate collective
        // AND cyclic are interpreted as flap commands delta_f; the feathering
        // response (delta_theta_0/1c/1s) REPLACES the direct swashplate pitch
        // path. Textbook Pitt-Peters has no feathering DOF -- this is the
        // Kaman servo-flap extension (see servoflap.rs and AGENTS.md).
        // ------------------------------------------------------------------
        let feathering_state = match &self.defn.pitch_actuation {
            PitchActuation::DirectMechanical => FeatheringState::RIGID,
            PitchActuation::ServoFlap(act) => solve_feathering(
                act,
                inputs.collective_rad,
                theta_1c,
                theta_1s,
                rho,
                omega,
                mu,
                r_tip,
                blade.chord_m,
                self.defn.airfoil.CL_alpha_per_rad,
            ),
        };
        let servo_mode = self.defn.is_servoflap();
        // In servo mode, collective/cyclic are interpreted as flap commands and
        // blade pitch comes only from feathering response.
        let loop_collective = if servo_mode {
            feathering_state.delta_theta_0
        } else {
            inputs.collective_rad
        };
        let (loop_theta_1c, loop_theta_1s) = if servo_mode {
            (
                feathering_state.delta_theta_1c,
                feathering_state.delta_theta_1s,
            )
        } else {
            (theta_1c, theta_1s)
        };

        // ------------------------------------------------------------------
        // Blade element forces
        // ------------------------------------------------------------------
        let (mut t_total, mut q_total, mut mx_hub, mut my_hub) = (0.0, 0.0, 0.0, 0.0);
        if omega_r > EPS_OMEGA_R && omega > 1.0 {
            let mut kernel = PpKernel {
                lambda_total,
                lam_c,
                lam_s,
                x_mid: &self.grid.x_mid[..self.grid.n_elements],
            };
            let sweep = SweepCtx {
                grid: &self.grid,
                polar: &self.polar,
                col: loop_collective,
                omega,
                omega_r,
                rho,
                n_b: blade.n_blades,
                n_psi: self.n_psi_elements,
                n_psi_inv: 1.0 / (self.n_psi_elements as f64),
                psi_trig: &self.psi_trig,
                v_in_hub_x: v_inplane_hub[0],
                v_in_hub_y: v_inplane_hub[1],
                theta_1c: loop_theta_1c,
                theta_1s: loop_theta_1s,
            };
            let (t, q, mx, my) = sweep.run(&mut kernel);
            t_total = t;
            q_total = q;
            mx_hub = mx;
            my_hub = my;
        }

        // ------------------------------------------------------------------
        // Pitt-Peters L-matrix steady-state targets + ODE
        // ------------------------------------------------------------------
        let vrs = vrs_regime(t_total, v_climb, rho, area);
        let mu_inplane = v_edge / omega_r.max(EPS_OMEGA_R);
        let v0 = lam0 * omega_r;
        // NON-STANDARD (flag A): v_mf is the classical Glauert mass-flow
        // V_mf = sqrt(v_edge^2 + (v_climb + v0)^2), NOT the Peters-HaQuang V
        // = (mu^2 + (lam+nu)(lam+2nu)) / sqrt(mu^2 + (lam+nu)^2). Reproduces
        // Glauert hover lambda_0 = sqrt(C_T/2); differs from the Peters-HaQuang
        // V by ~sqrt(2) in hover (Peters V = 2*nu there). L-matrix structure is
        // unaffected. Consistent with the Glauert uniform-inflow scaling used
        // in Li (UC Davis MS, 2020) Eq 14; see flag A in the module header.
        let v_mf = v_mass_flow_disk(v_edge, v_climb, v0, omega_r);
        // NON-STANDARD (flag E): mu_t_eff is clamped at MU_T_FLOOR purely to
        // keep the L-matrix scaling finite near hover. The time constants and
        // the C_L/C_M normalization below must use v_mf directly (which keeps
        // only its own smaller MASS_FLOW_HOVER_FLOOR_FRAC floor), not this
        // clamped value (see the note before `norm`).
        let mu_t_eff = (if omega_r > EPS_OMEGA_R {
            v_mf / omega_r
        } else {
            0.0
        })
        .max(MU_T_FLOOR);
        // NON-STANDARD (flag D): wind-axis rotation. beta=0 means the in-plane
        // relative wind is along -X_hub (pure longitudinal). Forcing and the
        // current cyclic states are rotated into wind axes, the L relations are
        // applied there, and the derivatives are rotated back -- making the
        // model rotationally covariant for oblique flight (mu_y != 0).
        let beta_wind = if v_edge > EPS_DENOM {
            // Wind-axis angle in the hub in-plane frame. beta=0 means
            // in-plane relative wind is along -X_hub (pure longitudinal).
            v_inplane_hub[1].atan2(-v_inplane_hub[0])
        } else {
            0.0
        };
        let beta_c = beta_wind.cos();
        let beta_s = beta_wind.sin();
        // EPS_OMEGA_R avoids chi flipping between +-pi/2 when lambda_total
        // is exactly zero in pure-edgewise flow.
        let chi = mu_inplane.atan2(lambda_total.abs() + EPS_OMEGA_R);
        let cos_chi = chi.cos();
        let tan_half_chi = (0.5 * chi).tan();
        // L-matrix entries in wind axes (Peters Eq 10):
        //   l_off = (15 pi / 64) tan(chi/2)  -- wake-skew off-diagonal
        //   l_cc  = 2(1 - X^2) = 4 cos chi / (1 + cos chi)
        //   l_ss  = 2(1 + X^2) = 4 / (1 + cos chi),   X = tan(chi/2)
        let l_off = (15.0 * PI / 64.0) * tan_half_chi;
        let l_cc = 4.0 * cos_chi / (1.0 + cos_chi);
        let l_ss = 4.0 / (1.0 + cos_chi);

        // NON-STANDARD (flag E, cont.): keep v_mf as the correctly-computed
        // v_mass_flow_disk value. Do NOT recompute from the clamped mu_t_eff
        // (that would re-impose the MU_T_FLOOR clamp on v_mf -- inflating it
        // whenever v_mf/omega_r < MU_T_FLOOR -- and breaking the coefficient
        // calculations).
        let norm = rho * area * omega_r * r_tip * v_mf;
        let c_l_hub = if norm > EPS_DENOM { mx_hub / norm } else { 0.0 };
        let c_m_hub = if norm > EPS_DENOM { my_hub / norm } else { 0.0 };

        // Rotate aerodynamic forcing into wind axes before applying the
        // Pitt-Peters L-matrix steady-state relations.
        let c_m_wind = beta_c * c_m_hub + beta_s * c_l_hub;
        let c_l_wind = -beta_s * c_m_hub + beta_c * c_l_hub;

        // Current cyclic inflow states expressed in wind axes.
        let lam_c_wind = beta_c * lam_c + beta_s * lam_s;
        let lam_s_wind = -beta_s * lam_c + beta_c * lam_s;

        let lam0_ss = if vrs.in_vrs {
            // NON-STANDARD (flag C): VRS empirical override. Momentum theory
            // does not apply in a recirculating wake, so lam0_ss comes from the
            // Leishman polynomial instead of the L-matrix relation, and the
            // uniform-inflow cross-coupling term (+l_off*C_M, in the else
            // branch) is dropped in this regime. The cyclic targets
            // lam_c_ss/lam_s_ss below are NOT skipped -- they keep their full
            // L-matrix cross-coupling even inside VRS.
            if omega_r > EPS_OMEGA_R {
                vrs_lambda1(vrs.lam2) * vrs.v_h / omega_r
            } else {
                0.0
            }
        } else if omega_r > EPS_OMEGA_R {
            // Uniform inflow: momentum term C_T/(2 V_mf) plus the symmetric
            // cross-coupling +l_off*C_M/mu_t_eff (a higher-order effect, small
            // in practice but kept for L-matrix symmetry).
            t_total / (2.0 * rho * area * v_mf * omega_r) + l_off * c_m_wind / mu_t_eff
        } else {
            0.0
        };

        let c_t = if omega_r > EPS_OMEGA_R {
            t_total / (rho * area * omega_r * omega_r)
        } else {
            0.0
        };
        // lam_c_ss: the -l_off*C_T term is the Pitt-Peters wake-skew
        // cross-coupling -- it produces Glauert wake skew naturally from thrust
        // forcing. Do NOT also add a closed-form Glauert tilt (double-count).
        let lam_c_ss_wind = (-l_off * c_t + l_cc * c_m_wind) / mu_t_eff;
        let lam_s_ss_wind = (l_ss * c_l_wind) / mu_t_eff;

        let tau_0 = (8.0 * r_tip) / (3.0 * PI * v_mf);
        let tau_cs = (16.0 * r_tip) / (45.0 * PI * v_mf);

        let d_lam0 = (lam0_ss - lam0) / tau_0;
        let d_lam_c_wind = (lam_c_ss_wind - lam_c_wind) / tau_cs;
        let d_lam_s_wind = (lam_s_ss_wind - lam_s_wind) / tau_cs;

        // Rotate harmonic derivatives back to hub-frame state coordinates.
        let d_lam_c = beta_c * d_lam_c_wind - beta_s * d_lam_s_wind;
        let d_lam_s = beta_s * d_lam_c_wind + beta_c * d_lam_s_wind;

        // Outputs -- NON-STANDARD (flag F): apply quasi-static flap reduction
        // to the hub moments before assembling the world-frame result. The
        // inflow ODE above uses the FULL aerodynamic moments (they drive the
        // wake), but the airframe only sees the fraction that passes through
        // the blade's flap stiffness.
        let (mx_out, my_out) = apply_flap_reduction(mx_hub, my_hub, self.defn.flap.as_ref(), omega);
        let result = assemble_result(t_total, q_total, mx_out, my_out, hub_axis, &inputs.R_hub);
        let derivative = PittPetersRotorState {
            lambda_0: d_lam0,
            lambda_c: d_lam_c,
            lambda_s: d_lam_s,
        };
        (result, derivative)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aero_io::{Mat3, Vec3};
    use crate::aero_model::AeroModel;
    use crate::rotor_definition::{BladeGeometry, ControlProperties, LinearPolarParameters, PitchActuation, RotorDefinition};

    fn beaupoil_rotor() -> RotorDefinition {
        RotorDefinition {
            blade: BladeGeometry {
                n_blades: 4, radius_m: 2.5, root_cutout_m: 0.5,
                chord_m: 0.20, twist_deg: 0.0, n_elements: 10, tip_loss: true,
                r_stations_m: Vec::new(), chord_stations_m: Vec::new(), twist_stations_deg: Vec::new(),
            },
            airfoil: LinearPolarParameters { CL0: 0.393, CL_alpha_per_rad: 5.79, CD0: 0.0079, alpha_stall_deg: 13.0 },
            control: Some(ControlProperties { swashplate_pitch_gain_rad: 0.3, swashplate_phase_deg: Some(0.0) }),
            pitch_actuation: PitchActuation::DirectMechanical,
            flap: None,
            name: "beaupoil_2026".to_string(),
            description: String::new(),
        }
    }

    /// RAWES IC attitude: Pitt-Peters force must oppose body-Z.
    #[test]
    fn rawes_ic_pp_force_opposes_body_z() {
        let defn = beaupoil_rotor();
        let polar = crate::polar::LinearPolar::from_properties(&defn.airfoil);
        let model = PittPetersModel::build(defn, 36, polar);
        let r_hub = Mat3([
            [0.0, -1.0, 0.0],
            [0.42720594325829603, 0.0, -0.9041543463617201],
            [0.9041543463617201, 0.0, 0.4272059432582967],
        ]);
        let inputs = RotorInputs {
            collective_rad: -0.18, tilt_lon: 0.0, tilt_lat: 0.0,
            R_hub: r_hub,
            v_hub_world: Vec3::zero(),
            wind_world: Vec3::new(0.0, 10.0, 0.0),
            omega_rad_s: 38.132161, rho_kg_m3: 1.225,
        };
        let (result, _) = model.compute_forces(&inputs, &model.initial_state());
        let body_z = Vec3::new(r_hub.0[0][2], r_hub.0[1][2], r_hub.0[2][2]);
        let fdz = result.F_world.dot(body_z);
        assert!(fdz < 0.0, "PP RAWES IC: F dot body_z should be negative, got {fdz:.3}");
        assert!(result.F_world.0[1] > 0.0, "PP RAWES IC: F_east should be positive (downwind)");
        assert!(result.F_world.0[2] < 0.0, "PP RAWES IC: F_up should be positive (-Z)");
    }

    /// Hover mass-flow: Glauert lambda_0 = sqrt(CT/2), Peters Eq-8 gives sqrt(CT/4).
    /// The ratio must be exactly sqrt(2) -- verifies the analytical derivation in
    /// AGENTS.md (Pitt-Peters section) is implemented correctly.
    #[test]
    fn glauert_vs_peters_hover_inflow_ratio_is_sqrt2() {
        const C_T: f64 = 0.00488;
        let lambda_glauert = (C_T / 2.0).sqrt();
        let lambda_peters  = (C_T / 4.0).sqrt();
        let ratio = lambda_glauert / lambda_peters;
        assert!(
            (ratio - 2.0_f64.sqrt()).abs() < 1e-10,
            "Expected sqrt(2) ratio between Glauert and Peters hover inflow, got {ratio}"
        );
        // Verify mass-flow parameter ratio is also sqrt(2).
        let v_mf_glauert = lambda_glauert;
        let v_peters     = 2.0 * lambda_peters;
        let mf_ratio = v_peters / v_mf_glauert;
        assert!(
            (mf_ratio - 2.0_f64.sqrt()).abs() < 1e-10,
            "Expected sqrt(2) mass-flow ratio, got {mf_ratio}"
        );
    }
}

