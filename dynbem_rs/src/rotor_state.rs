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

#[derive(Clone, Debug, Default)]
pub struct QuasiStaticRotorState;

#[derive(Clone, Debug, Default)]
pub struct PittPetersRotorState {
    pub lambda_0: f64,
    pub lambda_c: f64,
    pub lambda_s: f64,
}

#[derive(Clone, Debug)]
#[allow(non_snake_case)]
pub struct OyeRotorState {
    pub n_elements: usize,
    pub W_int: Vec<f64>,
    pub W: Vec<f64>,
}

impl OyeRotorState {
    pub fn zeros(n_elements: usize) -> Self {
        Self {
            n_elements,
            W_int: vec![0.0; n_elements],
            W: vec![0.0; n_elements],
        }
    }

    #[inline(always)]
    pub fn w_int_slice(&self) -> &[f64] {
        &self.W_int[..self.n_elements]
    }

    #[inline(always)]
    pub fn w_slice(&self) -> &[f64] {
        &self.W[..self.n_elements]
    }
}
