// Internal traits unifying the three aero models.
//
// Lets the trim solver (and future generic code) be written once over
// any model type. Not exposed to Python -- pyo3 can't dispatch over
// generic Rust traits from Python. Instead `AeroAny` (an enum wrapping
// the three concrete models) is the Python-facing entry point that
// resolves to a trait-using generic function internally.

use crate::aero_io::{AeroResult, RotorInputs};

/// Common aero-model interface. Each implementor caches a polar table
/// and a radial grid in its struct; compute_forces is the hot path.
pub trait AeroModel {
    type State: RotorStateExt;

    /// Forces + state derivative for one timestep.
    fn compute_forces(
        &self,
        inputs: &RotorInputs,
        state: &Self::State,
    ) -> (AeroResult, Self::State);

    /// Time constants per inflow state DOF. Used by the semi-implicit damping
    /// in the trim integrator. Default = all-infinity (no dynamic-inflow lags);
    fn inflow_taus(&self, _inputs: &RotorInputs, state: &Self::State) -> Vec<f64> {
        vec![f64::INFINITY; state.inflow_dof()]
    }

    /// Zero state at the right shape for this model.
    fn initial_state(&self) -> Self::State;
}

// Rotor state types: quasi-static, Pitt-Peters, Oye.
// States carry inflow DOFs only.  omega_rad_s lives in RotorInputs and is
// supplied by the caller on every compute_forces call; the mechanical ODE
// is NOT part of this state -- the caller owns and integrates omega.
//
// Canonical integration pattern (explicit Euler, same dt as inflow loop):
//
//   let (result, dstate) = aero.compute_forces(&inputs, &state);
//   state = step_state(state, dstate, dt);                     // inflow
//   omega += dt * (motor_torque - result.Q_spin) / I_kgm2;     // spin
//   inputs.omega_rad_s = omega;                                 // feed back
//
// The Python helper `dynbem.mechanical.omega_derivative` and
// `euler_step_omega` implement the scalar spin ODE.
// The inflow states may be stiff (tau << dt); use the semi-implicit stepper
// in `dynbem.mechanical` / `envelope.point_mass._step_state_semi_implicit`
// when dt is large relative to the inflow time constants.


/// Inflow-state serialization for generic integrators.
/// Omega is part of RotorInputs; state carries only inflow DOFs.
/// Each concrete state's impl lives in its aero module (quasi_static_bem,
/// pitt_peters, oye).
pub trait RotorStateExt: Clone {
    fn get_inflow(&self) -> Vec<f64>;
    fn set_inflow(&mut self, arr: Vec<f64>);
    fn inflow_dof(&self) -> usize {
        self.get_inflow().len()
    }
}
