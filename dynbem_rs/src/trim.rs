// Steady-state cyclic trim solver. Generic over any AeroModel: the
// integrator only needs get_inflow/set_inflow on State, and compute_forces +
// inflow_taus on the model.
//
// The Python facade in the glue crate adds an AeroAny enum to dispatch
// across the three concrete model types at the Python boundary -- pyo3
// can't dispatch over Rust generics from Python, so the enum lives
// there and resolves to one of these generic functions internally.

use crate::aero_io::{Mat3, RotorInputs};
use crate::aero_model::{AeroModel, RotorStateExt};

/// One semi-implicit Euler step: damp dynamic-inflow states by
/// 1/(1 + dt/tau); explicit Euler on inflow states (tau = inf).
fn semi_implicit_step<M: AeroModel>(
    aero: &M,
    state: &M::State,
    derivative: &M::State,
    inputs: &RotorInputs,
    dt: f64,
) -> M::State {
    let taus = aero.inflow_taus(inputs, state);
    let arr = state.get_inflow();
    let darr = derivative.get_inflow();
    let n = arr.len();
    debug_assert_eq!(arr.len(), n);
    debug_assert_eq!(darr.len(), n);
    debug_assert_eq!(taus.len(), n);
    let mut new_arr = Vec::with_capacity(n);
    for i in 0..n {
        let damp = if taus[i].is_finite() {
            1.0 / (1.0 + dt / taus[i])
        } else {
            1.0
        };
        new_arr.push(arr[i] + dt * darr[i] * damp);
    }
    let mut out = state.clone();
    out.set_inflow(new_arr);
    out
}

/// Advance the state to quasi-steady inflow at fixed inputs.
pub fn relax_inflow<M: AeroModel>(
    aero: &M,
    mut state: M::State,
    inputs: &RotorInputs,
    n_steps: usize,
    dt: f64,
) -> M::State {
    for _ in 0..n_steps {
        let (_, deriv) = aero.compute_forces(inputs, &state);
        state = semi_implicit_step(aero, &state, &deriv, inputs, dt);
    }
    state
}

/// Hub-frame moment (Mx, My) at the given cyclic inputs, plus the state
/// derivative for the semi-implicit step.
fn eval_moment<M: AeroModel>(
    aero: &M,
    state: &M::State,
    inputs: &RotorInputs,
    r_hub: &Mat3,
    target_x: f64,
    target_y: f64,
) -> (f64, f64, M::State) {
    let (result, deriv) = aero.compute_forces(inputs, state);
    // R_hub.T @ M_orbital -- hub-frame moment.
    let m_hub = r_hub.transpose() * result.M_hub_world;
    (m_hub[0] - target_x, m_hub[1] - target_y, deriv)
}

#[derive(Clone, Debug)]
pub struct TrimOutcome<S: RotorStateExt> {
    pub tilt_lon: f64,
    pub tilt_lat: f64,
    pub mx_residual: f64,
    pub my_residual: f64,
    pub iterations: usize,
    pub converged: bool,
    pub final_state: S,
}

#[allow(clippy::too_many_arguments)]
pub fn solve_trim_cyclic<M: AeroModel>(
    aero: &M,
    mut state: M::State,
    base_inputs: &RotorInputs,
    target_x: f64,
    target_y: f64,
    tilt_lon_init: f64,
    tilt_lat_init: f64,
    tilt_min: f64,
    tilt_max: f64,
    tolerance_n_m: f64,
    max_iterations: usize,
    probe_rad: f64,
    dt_relax: f64,
    n_inflow_relax: usize,
    n_settle: usize,
) -> TrimOutcome<M::State> {
    let mut tilt_lon = tilt_lon_init.clamp(tilt_min, tilt_max);
    let mut tilt_lat = tilt_lat_init.clamp(tilt_min, tilt_max);

    let make_inputs = |tlon: f64, tlat: f64| {
        let mut i = base_inputs.clone();
        i.tilt_lon = tlon;
        i.tilt_lat = tlat;
        i
    };

    let relax = |s: M::State, tlon: f64, tlat: f64, n: usize| -> M::State {
        let inp = make_inputs(tlon, tlat);
        relax_inflow(aero, s, &inp, n, dt_relax)
    };

    if n_settle > 0 {
        state = relax(state, tilt_lon, tilt_lat, n_settle);
    }
    state = relax(state, tilt_lon, tilt_lat, n_inflow_relax);

    let mut inp = make_inputs(tilt_lon, tilt_lat);
    let (mut mx, mut my, _) =
        eval_moment(aero, &state, &inp, &base_inputs.R_hub, target_x, target_y);
    let mut converged = mx.abs() < tolerance_n_m && my.abs() < tolerance_n_m;
    let mut iter = 0usize;

    for k in 1..=max_iterations {
        iter = k;
        if converged {
            break;
        }
        // probe d(My)/d(tilt_lon)
        let inp_p = make_inputs(tilt_lon + probe_rad, tilt_lat);
        let (_, my_p, _) =
            eval_moment(aero, &state, &inp_p, &base_inputs.R_hub, target_x, target_y);
        let d_my_dlon = (my_p - my) / probe_rad;
        // probe d(Mx)/d(tilt_lat)
        let inp_p = make_inputs(tilt_lon, tilt_lat + probe_rad);
        let (mx_p, _, _) =
            eval_moment(aero, &state, &inp_p, &base_inputs.R_hub, target_x, target_y);
        let d_mx_dlat = (mx_p - mx) / probe_rad;

        if d_my_dlon.abs() > 1e-6 {
            tilt_lon = (tilt_lon - 0.5 * my / d_my_dlon).clamp(tilt_min, tilt_max);
        }
        if d_mx_dlat.abs() > 1e-6 {
            tilt_lat = (tilt_lat - 0.5 * mx / d_mx_dlat).clamp(tilt_min, tilt_max);
        }
        state = relax(state, tilt_lon, tilt_lat, n_inflow_relax);
        inp = make_inputs(tilt_lon, tilt_lat);
        let (mxi, myi, _) = eval_moment(aero, &state, &inp, &base_inputs.R_hub, target_x, target_y);
        mx = mxi;
        my = myi;
        converged = mx.abs() < tolerance_n_m && my.abs() < tolerance_n_m;
    }
    TrimOutcome {
        tilt_lon,
        tilt_lat,
        mx_residual: mx,
        my_residual: my,
        iterations: iter,
        converged,
        final_state: state,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aero_model::AeroModel;
    use crate::oye::OyeBEMModel;
    use crate::pitt_peters::PittPetersModel;
    use crate::polar::LinearPolar;
    use crate::rotor_definition::{
        BladeGeometry, ControlProperties, LinearPolarParameters, PitchActuation, RotorDefinition,
    };

    const OMEGA: f64 = 28.0;
    const COLLECTIVE: f64 = -9.0_f64.to_radians();
    const TOL_OYE: f64 = 0.05;
    const TOL_PITT: f64 = 1.6;

    fn beaupoil_rotor() -> RotorDefinition {
        RotorDefinition {
            blade: BladeGeometry {
                n_blades: 4,
                radius_m: 2.5,
                root_cutout_m: 0.5,
                chord_m: 0.20,
                twist_deg: 0.0,
                n_elements: 10,
                tip_loss: true,
                r_stations_m: Vec::new(),
                chord_stations_m: Vec::new(),
                twist_stations_deg: Vec::new(),
            },
            airfoil: LinearPolarParameters {
                CL0: 0.393,
                CL_alpha_per_rad: 5.79,
                CD0: 0.0079,
                alpha_stall_deg: 13.0,
            },
            control: Some(ControlProperties {
                swashplate_pitch_gain_rad: 0.3,
                swashplate_phase_deg: Some(0.0),
            }),
            pitch_actuation: PitchActuation::DirectMechanical,
            flap: None,
            name: "beaupoil_2026".to_string(),
            description: String::new(),
        }
    }

    fn base_inputs(wind_world: [f64; 3]) -> RotorInputs {
        use crate::aero_io::{Mat3, Vec3};
        RotorInputs {
            collective_rad: COLLECTIVE,
            tilt_lon: 0.0,
            tilt_lat: 0.0,
            R_hub: Mat3::eye(),
            v_hub_world: Vec3::zero(),
            wind_world: Vec3::new(wind_world[0], wind_world[1], wind_world[2]),
            rho_kg_m3: 1.225,
            omega_rad_s: OMEGA,
        }
    }

    fn moments_hub<M: AeroModel>(
        aero: &M,
        state: &M::State,
        tilt_lon: f64,
        tilt_lat: f64,
        wind_world: [f64; 3],
    ) -> (f64, f64) {
        let mut inputs = base_inputs(wind_world);
        inputs.collective_rad = COLLECTIVE;
        inputs.tilt_lon = tilt_lon;
        inputs.tilt_lat = tilt_lat;
        let (res, _) = aero.compute_forces(&inputs, state);
        let m_hub = inputs.R_hub.transpose() * res.M_hub_world;
        (m_hub[0], m_hub[1])
    }

    fn run_trim<M: AeroModel>(
        aero: &M,
        state: M::State,
        wind_world: [f64; 3],
        target_x: f64,
        target_y: f64,
        tol: f64,
    ) -> TrimOutcome<M::State> {
        let inputs = base_inputs(wind_world);
        solve_trim_cyclic(
            aero,
            state,
            &inputs,
            target_x,
            target_y,
            0.0,
            0.0,
            -0.261_799_387_799_149_4,
            0.261_799_387_799_149_4,
            tol,
            50,
            0.008_726_646_259_971_648,
            0.005,
            100,
            0,
        )
    }

    #[test]
    fn oye_hover_trim_is_near_zero() {
        let defn = beaupoil_rotor();
        let polar = LinearPolar::from_properties(&defn.airfoil);
        let aero = OyeBEMModel::build(defn, 36, polar);
        let out = run_trim(
            &aero,
            aero.initial_state(),
            [0.0, 0.0, 0.0],
            0.0,
            0.0,
            TOL_OYE,
        );
        assert!(
            out.converged,
            "Oye hover trim did not converge: iters={} mx={:.4} my={:.4}",
            out.iterations, out.mx_residual, out.my_residual
        );
        assert!(out.tilt_lon.abs() < 0.5_f64.to_radians());
        assert!(out.tilt_lat.abs() < 0.5_f64.to_radians());
    }

    #[test]
    fn pitt_forward_trim_residual_below_tolerance() {
        let defn = beaupoil_rotor();
        let polar = LinearPolar::from_properties(&defn.airfoil);
        let aero = PittPetersModel::build(defn, 36, polar);
        let out = run_trim(
            &aero,
            aero.initial_state(),
            [0.0, 10.0, 0.0],
            0.0,
            0.0,
            TOL_PITT,
        );
        assert!(
            out.converged,
            "Pitt forward trim did not converge: iters={} mx={:.4} my={:.4}",
            out.iterations, out.mx_residual, out.my_residual
        );
        assert!(out.mx_residual.abs() < TOL_PITT);
        assert!(out.my_residual.abs() < TOL_PITT);
    }

    #[test]
    fn trim_residual_matches_direct_evaluation_pitt() {
        let defn = beaupoil_rotor();
        let polar = LinearPolar::from_properties(&defn.airfoil);
        let aero = PittPetersModel::build(defn, 36, polar);
        let out = run_trim(
            &aero,
            aero.initial_state(),
            [0.0, 10.0, 0.0],
            0.0,
            0.0,
            TOL_PITT,
        );
        let (mx, my) = moments_hub(
            &aero,
            &out.final_state,
            out.tilt_lon,
            out.tilt_lat,
            [0.0, 10.0, 0.0],
        );
        assert!((mx - out.mx_residual).abs() < 1e-6);
        assert!((my - out.my_residual).abs() < 1e-6);
    }

    #[test]
    fn trim_to_nonzero_target_moment_pitt() {
        let defn = beaupoil_rotor();
        let polar = LinearPolar::from_properties(&defn.airfoil);
        let aero = PittPetersModel::build(defn, 36, polar);
        let m_target = 5.0;
        let out = run_trim(
            &aero,
            aero.initial_state(),
            [0.0, 10.0, 0.0],
            0.0,
            m_target,
            TOL_PITT,
        );
        assert!(
            out.converged,
            "Pitt target-moment trim did not converge: iters={} mx={:.4} my={:.4}",
            out.iterations, out.mx_residual, out.my_residual
        );
        let (mx, my) = moments_hub(
            &aero,
            &out.final_state,
            out.tilt_lon,
            out.tilt_lat,
            [0.0, 10.0, 0.0],
        );
        assert!(mx.abs() < TOL_PITT, "Mx={mx:.4} should be near 0");
        assert!(
            (my - m_target).abs() < TOL_PITT,
            "My={my:.4} should be near {m_target:.4}"
        );
    }

    #[test]
    fn relax_inflow_settles_to_fixed_point_for_both_models() {
        use crate::aero_model::RotorStateExt;
        let defn = beaupoil_rotor();
        let pp = PittPetersModel::build(
            defn.clone(),
            36,
            LinearPolar::from_properties(&defn.airfoil),
        );
        let oye = OyeBEMModel::build(
            defn.clone(),
            36,
            LinearPolar::from_properties(&defn.airfoil),
        );
        let inputs = base_inputs([0.0, 10.0, 0.0]);

        let s1_pp = relax_inflow(&pp, pp.initial_state(), &inputs, 500, 0.005);
        let s2_pp = relax_inflow(&pp, s1_pp.clone(), &inputs, 500, 0.005);
        let d_pp: f64 = s1_pp
            .get_inflow()
            .iter()
            .zip(s2_pp.get_inflow().iter())
            .map(|(a, b)| (b - a) * (b - a))
            .sum::<f64>()
            .sqrt();
        assert!(d_pp < 1e-4, "Pitt inflow not settled: delta={d_pp:.4e}");

        let s1_oye = relax_inflow(&oye, oye.initial_state(), &inputs, 500, 0.005);
        let s2_oye = relax_inflow(&oye, s1_oye.clone(), &inputs, 500, 0.005);
        let d_oye: f64 = s1_oye
            .get_inflow()
            .iter()
            .zip(s2_oye.get_inflow().iter())
            .map(|(a, b)| (b - a) * (b - a))
            .sum::<f64>()
            .sqrt();
        assert!(d_oye < 1e-4, "Oye inflow not settled: delta={d_oye:.4e}");
    }

    #[test]
    fn solver_reduces_baseline_disturbance_pitt() {
        use crate::aero_model::RotorStateExt;
        let defn = beaupoil_rotor();
        let polar = LinearPolar::from_properties(&defn.airfoil);
        let aero = PittPetersModel::build(defn, 36, polar);
        let wind = [0.0, 10.0, 0.0];

        let mut state = aero.initial_state();
        let inputs = base_inputs(wind);
        for _ in 0..200 {
            let (_, deriv) = aero.compute_forces(&inputs, &state);
            let arr: Vec<f64> = state
                .get_inflow()
                .iter()
                .zip(deriv.get_inflow().iter())
                .map(|(x, dx)| x + 0.005 * dx)
                .collect();
            state.set_inflow(arr);
        }

        let (mx0, my0) = moments_hub(&aero, &state, 0.0, 0.0, wind);
        let baseline = (mx0 * mx0 + my0 * my0).sqrt();
        assert!(baseline > 10.0, "baseline too small: {baseline:.2}");

        let out = run_trim(&aero, state, wind, 0.0, 0.0, TOL_PITT);
        let trim_mag =
            (out.mx_residual * out.mx_residual + out.my_residual * out.my_residual).sqrt();
        assert!(
            trim_mag < baseline / 70.0,
            "solver did not cancel disturbance: baseline={baseline:.2} trim={trim_mag:.4}"
        );
    }

    #[test]
    fn trim_clips_to_bounds_oye() {
        let defn = beaupoil_rotor();
        let polar = LinearPolar::from_properties(&defn.airfoil);
        let aero = OyeBEMModel::build(defn, 36, polar);
        let tight = 1.0_f64.to_radians();
        let out = solve_trim_cyclic(
            &aero,
            aero.initial_state(),
            &base_inputs([0.0, 10.0, 0.0]),
            0.0,
            0.0,
            0.0,
            0.0,
            -tight,
            tight,
            0.01,
            20,
            0.008_726_646_259_971_648,
            0.005,
            100,
            0,
        );
        assert!(out.tilt_lon >= -tight && out.tilt_lon <= tight);
        assert!(out.tilt_lat >= -tight && out.tilt_lat <= tight);
    }
}
