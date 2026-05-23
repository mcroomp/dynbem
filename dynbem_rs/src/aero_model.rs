// Internal traits unifying the three aero models.
//
// Lets the trim solver (and future generic code) be written once over
// any model type. Not exposed to Python -- pyo3 can't dispatch over
// generic Rust traits from Python. Instead `AeroAny` (an enum wrapping
// the three concrete models) is the Python-facing entry point that
// resolves to a trait-using generic function internally.

use crate::aero_io::{AeroResult, RotorInputs};
use crate::rotor_state::{OyeRotorState, PittPetersRotorState, QuasiStaticRotorState};

/// Inflow-state serialization for generic integrators.
/// Mechanical state (omega/spin) stays explicit on the typed state.
pub trait RotorStateExt: Clone {
    fn get_inflow(&self) -> Vec<f64>;
    fn set_inflow(&mut self, arr: Vec<f64>);
    fn omega(&self) -> f64;
    fn set_omega(&mut self, omega: f64);
    fn spin(&self) -> f64;
    fn set_spin(&mut self, spin: f64);
    fn inflow_dof(&self) -> usize {
        self.get_inflow().len()
    }
}

impl RotorStateExt for QuasiStaticRotorState {
    fn get_inflow(&self) -> Vec<f64> {
        Vec::new()
    }
    fn set_inflow(&mut self, arr: Vec<f64>) {
        debug_assert!(arr.is_empty());
    }
    fn omega(&self) -> f64 {
        self.omega_rad_s
    }
    fn set_omega(&mut self, omega: f64) {
        self.omega_rad_s = omega;
    }
    fn spin(&self) -> f64 {
        self.spin_angle_rad
    }
    fn set_spin(&mut self, spin: f64) {
        self.spin_angle_rad = spin;
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
    fn omega(&self) -> f64 {
        self.omega_rad_s
    }
    fn set_omega(&mut self, omega: f64) {
        self.omega_rad_s = omega;
    }
    fn spin(&self) -> f64 {
        self.spin_angle_rad
    }
    fn set_spin(&mut self, spin: f64) {
        self.spin_angle_rad = spin;
    }
}

impl RotorStateExt for OyeRotorState {
    fn get_inflow(&self) -> Vec<f64> {
        let n = self.n_elements;
        let mut v = Vec::with_capacity(2 * n);
        v.extend_from_slice(self.w_int_slice());
        v.extend_from_slice(self.w_slice());
        v
    }
    fn set_inflow(&mut self, arr: Vec<f64>) {
        let n = self.n_elements;
        debug_assert_eq!(arr.len(), 2 * n);
        self.W_int[..n].copy_from_slice(&arr[..n]);
        self.W[..n].copy_from_slice(&arr[n..2 * n]);
    }
    fn omega(&self) -> f64 {
        self.omega_rad_s
    }
    fn set_omega(&mut self, omega: f64) {
        self.omega_rad_s = omega;
    }
    fn spin(&self) -> f64 {
        self.spin_angle_rad
    }
    fn set_spin(&mut self, spin: f64) {
        self.spin_angle_rad = spin;
    }
}

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

    /// Time constants per state DOF (infinite for mechanical / quasi-static
    /// states). Used by the semi-implicit damping in the trim integrator.
    /// Default = all-infinity (no dynamic-inflow lags); dynamic-inflow
    /// models override this with their per-state lag formulas.
    fn inflow_taus(&self, _inputs: &RotorInputs, state: &Self::State) -> Vec<f64> {
        vec![f64::INFINITY; state.inflow_dof()]
    }

    /// Zero state at the right shape for this model.
    fn initial_state(&self) -> Self::State;
}
