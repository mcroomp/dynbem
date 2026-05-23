// Rotor state types: quasi-static, Pitt-Peters, Oye.
// State-vector layout (including omega/spin indices) is reported by
// each RotorStateExt implementation in aero_model.rs.

#[derive(Clone, Debug, Default)]
pub struct QuasiStaticRotorState {
    pub omega_rad_s: f64,
    pub spin_angle_rad: f64,
}

#[derive(Clone, Debug, Default)]
pub struct PittPetersRotorState {
    pub lambda_0: f64,
    pub lambda_c: f64,
    pub lambda_s: f64,
    pub omega_rad_s: f64,
    pub spin_angle_rad: f64,
}

#[derive(Clone, Debug)]
#[allow(non_snake_case)]
pub struct OyeRotorState {
    pub n_elements: usize,
    pub W_int: Vec<f64>,
    pub W: Vec<f64>,
    pub omega_rad_s: f64,
    pub spin_angle_rad: f64,
}

impl OyeRotorState {
    pub fn zeros(n_elements: usize, omega_rad_s: f64) -> Self {
        Self {
            n_elements,
            W_int: vec![0.0; n_elements],
            W: vec![0.0; n_elements],
            omega_rad_s,
            spin_angle_rad: 0.0,
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
