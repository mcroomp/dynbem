// Steady-state cyclic trim solver. Generic over any AeroModel: the
// integrator only needs to_vec/from_vec on State, and compute_forces +
// inflow_taus on the model.
//
// The Python facade in the glue crate adds an AeroAny enum to dispatch
// across the three concrete model types at the Python boundary -- pyo3
// can't dispatch over Rust generics from Python, so the enum lives
// there and resolves to one of these generic functions internally.

use crate::aero_io::{Mat3, RotorInputs};
use crate::aero_model::{AeroModel, RotorStateExt};

/// One semi-implicit Euler step: damp dynamic-inflow states by
/// 1/(1 + dt/tau); explicit Euler on mechanical (tau = inf). Optionally
/// clamps omega to a fixed value (state convention: omega is the
/// second-to-last entry).
fn semi_implicit_step<M: AeroModel>(
    aero: &M,
    state: &M::State,
    derivative: &M::State,
    inputs: &RotorInputs,
    dt: f64,
    fix_omega_to: Option<f64>,
) -> M::State {
    let taus = aero.inflow_taus(inputs, state);
    let arr = state.to_vec();
    let darr = derivative.to_vec();
    let n = arr.len();
    let mut new_arr = Vec::with_capacity(n);
    for i in 0..n {
        let damp = if taus[i].is_finite() {
            1.0 / (1.0 + dt / taus[i])
        } else {
            1.0
        };
        new_arr.push(arr[i] + dt * darr[i] * damp);
    }
    if let Some(om) = fix_omega_to {
        new_arr[n - 2] = om;
    }
    M::State::from_vec(&new_arr)
}

/// Advance the state to quasi-steady inflow at fixed inputs.
pub fn relax_inflow<M: AeroModel>(
    aero: &M,
    mut state: M::State,
    inputs: &RotorInputs,
    n_steps: usize,
    dt: f64,
    fix_omega: bool,
) -> M::State {
    let fix_to = if fix_omega { Some(state.omega()) } else { None };
    for _ in 0..n_steps {
        let (_, deriv) = aero.compute_forces(inputs, &state);
        state = semi_implicit_step(aero, &state, &deriv, inputs, dt, fix_to);
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
    let m_hub = r_hub.transpose() * result.M_orbital;
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
    fix_omega: bool,
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
        relax_inflow(aero, s, &inp, n, dt_relax, fix_omega)
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
