// Rotor state types: quasi-static, Pitt-Peters, Oye. Mechanical states
// (omega, spin_angle) are ALWAYS the last two entries in the state vector
// serialization (see aero_model::RotorStateExt).

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
    pub W_int: Vec<f64>,
    pub W: Vec<f64>,
    pub omega_rad_s: f64,
    pub spin_angle_rad: f64,
}

impl OyeRotorState {
    pub fn zeros(n_elements: usize, omega_rad_s: f64) -> Self {
        Self {
            W_int: vec![0.0; n_elements],
            W: vec![0.0; n_elements],
            omega_rad_s,
            spin_angle_rad: 0.0,
        }
    }
}
