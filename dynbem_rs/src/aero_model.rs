// Internal traits unifying the three aero models.
//
// Lets the trim solver (and future generic code) be written once over
// any model type. Not exposed to Python -- pyo3 can't dispatch over
// generic Rust traits from Python. Instead `AeroAny` (an enum wrapping
// the three concrete models) is the Python-facing entry point that
// resolves to a trait-using generic function internally.

use crate::aero_io::{AeroResult, RotorInputs};
use crate::rotor_state::{OyeRotorState, PittPetersRotorState, QuasiStaticRotorState};

/// State vector serialization for the ODE integrator.
/// Convention: omega_rad_s is the **second-to-last** entry, spin_angle_rad
/// is the last (so the trim solver can clamp omega without knowing the
/// state layout).
pub trait RotorStateExt: Clone {
    fn to_vec(&self) -> Vec<f64>;
    fn from_vec(arr: &[f64]) -> Self;
    fn omega(&self) -> f64;
    fn n_dof(&self) -> usize;
}

impl RotorStateExt for QuasiStaticRotorState {
    fn to_vec(&self) -> Vec<f64> {
        vec![self.omega_rad_s, self.spin_angle_rad]
    }
    fn from_vec(arr: &[f64]) -> Self {
        QuasiStaticRotorState {
            omega_rad_s: arr[0],
            spin_angle_rad: arr[1],
        }
    }
    fn omega(&self) -> f64 {
        self.omega_rad_s
    }
    fn n_dof(&self) -> usize {
        2
    }
}

impl RotorStateExt for PittPetersRotorState {
    fn to_vec(&self) -> Vec<f64> {
        vec![
            self.lambda_0,
            self.lambda_c,
            self.lambda_s,
            self.omega_rad_s,
            self.spin_angle_rad,
        ]
    }
    fn from_vec(arr: &[f64]) -> Self {
        PittPetersRotorState {
            lambda_0: arr[0],
            lambda_c: arr[1],
            lambda_s: arr[2],
            omega_rad_s: arr[3],
            spin_angle_rad: arr[4],
        }
    }
    fn omega(&self) -> f64 {
        self.omega_rad_s
    }
    fn n_dof(&self) -> usize {
        5
    }
}

impl RotorStateExt for OyeRotorState {
    fn to_vec(&self) -> Vec<f64> {
        let n = self.W_int.len();
        let mut v = Vec::with_capacity(2 * n + 2);
        v.extend_from_slice(&self.W_int);
        v.extend_from_slice(&self.W);
        v.push(self.omega_rad_s);
        v.push(self.spin_angle_rad);
        v
    }
    fn from_vec(arr: &[f64]) -> Self {
        let n_total = arr.len();
        let n = (n_total - 2) / 2;
        OyeRotorState {
            W_int: arr[..n].to_vec(),
            W: arr[n..2 * n].to_vec(),
            omega_rad_s: arr[n_total - 2],
            spin_angle_rad: arr[n_total - 1],
        }
    }
    fn omega(&self) -> f64 {
        self.omega_rad_s
    }
    fn n_dof(&self) -> usize {
        2 * self.W_int.len() + 2
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
    fn inflow_taus(&self, inputs: &RotorInputs, state: &Self::State) -> Vec<f64>;

    /// Zero state at the right shape for this model.
    fn initial_state(&self) -> Self::State;
}
