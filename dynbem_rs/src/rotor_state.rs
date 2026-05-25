// Rotor state types: quasi-static, Pitt-Peters, Oye.
// States carry inflow DOFs only.  omega_rad_s lives in RotorInputs and is
// supplied by the caller on every compute_forces call; the mechanical ODE
// (d_omega/dt) is integrated externally by the caller.

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
