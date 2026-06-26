// Internal traits unifying the three aero models.
//
// Lets the trim solver (and future generic code) be written once over
// any model type. Not exposed to Python -- pyo3 can't dispatch over
// generic Rust traits from Python. Instead `AeroAny` (an enum wrapping
// the three concrete models) is the Python-facing entry point that
// resolves to a trait-using generic function internally.

use crate::aero_io::{AeroResult, RotorInputs};

/// Integration method for inflow state updates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IntegrationMethod {
    /// Explicit Euler: lam_{n+1} = lam_n + dt * dlam_n
    ExplicitEuler,
    /// Semi-implicit Euler with tau damping: (lam + dt*dlam)/(1 + dt/tau)
    SemiImplicitEuler,
    /// Exact first-order relaxation update (frozen target over dt):
    /// lam_{n+1} = lam_n + tau * (1 - exp(-dt/tau)) * dlam_n
    ExponentialRelaxation,
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

    /// Time constants per inflow state DOF. Used by the semi-implicit damping
    /// in the trim integrator. Default = all-infinity (no dynamic-inflow lags);
    fn inflow_taus(&self, _inputs: &RotorInputs, state: &Self::State) -> Vec<f64> {
        vec![f64::INFINITY; state.inflow_dof()]
    }

    /// Zero state at the right shape for this model.
    fn initial_state(&self) -> Self::State;

    /// Compute forces and advance the inflow state by `dt` seconds using the
    /// selected integration method.
    ///
    /// For each inflow DOF i:
    /// - ExplicitEuler:     new_lam[i] = lam[i] + dt * dlam[i]
    /// - SemiImplicitEuler: new_lam[i] = (lam[i] + dt * dlam[i]) / (1 + dt / tau[i])
    /// - ExponentialRelaxation: new_lam[i] = lam[i] + tau[i] * (1 - exp(-dt/tau[i])) * dlam[i]
    ///
    /// For SemiImplicitEuler, when tau[i] == f64::INFINITY (quasi-static DOF)
    /// this reduces to plain explicit Euler. When dt >> tau the scheme is
    /// unconditionally stable and damps to the steady-state target without
    /// oscillation.
    ///
    /// ExponentialRelaxation is exact for first-order lag states when the
    /// steady-state target is treated as frozen over the step.
    ///
    /// Returns (AeroResult at the start-of-step state, new integrated state).
    fn step(
        &self,
        inputs: &RotorInputs,
        state: &Self::State,
        dt: f64,
        method: IntegrationMethod,
    ) -> (AeroResult, Self::State) {
        let (result, d_state) = self.compute_forces(inputs, state);
        let taus = self.inflow_taus(inputs, state);
        let inflow = state.get_inflow();
        let d_inflow = d_state.get_inflow();
        let new_inflow: Vec<f64> = inflow
            .iter()
            .zip(d_inflow.iter())
            .zip(taus.iter())
            .map(|((lam, dlam), tau)| {
                let explicit = lam + dt * dlam;
                match method {
                    IntegrationMethod::ExplicitEuler => explicit,
                    IntegrationMethod::SemiImplicitEuler => {
                        if tau.is_finite() {
                            explicit / (1.0 + dt / tau)
                        } else {
                            explicit
                        }
                    }
                    IntegrationMethod::ExponentialRelaxation => {
                        if tau.is_finite() && tau.abs() > f64::EPSILON {
                            let gain = -(-dt / tau).exp_m1();
                            lam + tau * gain * dlam
                        } else {
                            explicit
                        }
                    }
                }
            })
            .collect();
        let mut new_state = state.clone();
        new_state.set_inflow(new_inflow);
        (result, new_state)
    }
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

#[cfg(test)]
mod tests {
    use super::{AeroModel, IntegrationMethod, RotorStateExt};
    use crate::aero_io::{AeroResult, Mat3, RotorInputs, Vec3};

    #[derive(Clone, Debug)]
    struct DummyState(Vec<f64>);

    impl RotorStateExt for DummyState {
        fn get_inflow(&self) -> Vec<f64> {
            self.0.clone()
        }

        fn set_inflow(&mut self, arr: Vec<f64>) {
            self.0 = arr;
        }
    }

    struct DummyModel {
        d_inflow: Vec<f64>,
        taus: Option<Vec<f64>>,
    }

    impl AeroModel for DummyModel {
        type State = DummyState;

        fn compute_forces(
            &self,
            _inputs: &RotorInputs,
            _state: &Self::State,
        ) -> (AeroResult, Self::State) {
            (
                AeroResult {
                    F_world: Vec3::new(1.0, 2.0, 3.0),
                    M_hub_world: Vec3::new(4.0, 5.0, 6.0),
                    Q_spin: 7.0,
                    M_spin: Vec3::new(8.0, 9.0, 10.0),
                },
                DummyState(self.d_inflow.clone()),
            )
        }

        fn inflow_taus(&self, _inputs: &RotorInputs, state: &Self::State) -> Vec<f64> {
            self.taus
                .clone()
                .unwrap_or_else(|| vec![f64::INFINITY; state.inflow_dof()])
        }

        fn initial_state(&self) -> Self::State {
            DummyState(vec![])
        }
    }

    fn dummy_inputs() -> RotorInputs {
        RotorInputs {
            collective_rad: 0.0,
            tilt_lon: 0.0,
            tilt_lat: 0.0,
            R_hub: Mat3::eye(),
            v_hub_world: Vec3::zero(),
            wind_world: Vec3::zero(),
            rho_kg_m3: 1.225,
            omega_rad_s: 100.0,
        }
    }

    #[test]
    fn step_uses_explicit_update_for_infinite_tau() {
        let model = DummyModel {
            d_inflow: vec![2.0],
            taus: None,
        };
        let state = DummyState(vec![1.0]);
        let (result, next_state) = model.step(
            &dummy_inputs(),
            &state,
            0.5,
            IntegrationMethod::SemiImplicitEuler,
        );

        assert_eq!(next_state.get_inflow(), vec![2.0]);
        assert_eq!(result.Q_spin, 7.0);
    }

    #[test]
    fn step_uses_explicit_method_update() {
        let model = DummyModel {
            d_inflow: vec![2.0],
            taus: Some(vec![0.1]),
        };
        let state = DummyState(vec![1.0]);
        let (_, next_state) = model.step(
            &dummy_inputs(),
            &state,
            0.5,
            IntegrationMethod::ExplicitEuler,
        );

        // explicit = 1.0 + 0.5 * 2.0 = 2.0
        assert!((next_state.get_inflow()[0] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn step_uses_semi_implicit_update_for_finite_tau() {
        let model = DummyModel {
            d_inflow: vec![1.0],
            taus: Some(vec![0.1]),
        };
        let state = DummyState(vec![1.0]);
        let (_, next_state) = model.step(
            &dummy_inputs(),
            &state,
            0.1,
            IntegrationMethod::SemiImplicitEuler,
        );

        // explicit = 1.0 + 0.1 * 1.0 = 1.1
        // denominator = 1 + dt/tau = 1 + 0.1/0.1 = 2
        // new = 1.1 / 2 = 0.55
        assert!((next_state.get_inflow()[0] - 0.55).abs() < 1e-12);
    }

    #[test]
    fn step_uses_exponential_relaxation_for_finite_tau() {
        let model = DummyModel {
            d_inflow: vec![1.0],
            taus: Some(vec![0.1]),
        };
        let state = DummyState(vec![1.0]);
        let (_, next_state) = model.step(
            &dummy_inputs(),
            &state,
            0.1,
            IntegrationMethod::ExponentialRelaxation,
        );

        // lam_new = lam + tau * (1 - exp(-dt/tau)) * dlam
        //         = 1 + 0.1 * (1 - exp(-1)) * 1 = 1.0632120558828558
        assert!((next_state.get_inflow()[0] - 1.0632120558828558).abs() < 1e-12);
    }
}
