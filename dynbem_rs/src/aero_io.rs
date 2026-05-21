// RotorInputs and AeroResult: the common compute_forces() boundary types.
// See ../../../CLAUDE.md for sign conventions and the NED frame.

// ---------------------------------------------------------------------------
// Vec3 / Mat3: tiny fixed-size newtypes wrapping plain f64 arrays.
//
// Same codegen as raw [f64; 3] / [[f64; 3]; 3] (Copy, no indirection), but
// the operator impls let call sites read like math: `R * v`, `a - b`,
// `s * v`. No external dependency; LLVM still autovectorizes the bodies.
// ---------------------------------------------------------------------------

#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Vec3(pub [f64; 3]);

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Mat3(pub [[f64; 3]; 3]);

impl Vec3 {
    #[inline]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3([x, y, z])
    }
    #[inline]
    pub const fn zero() -> Self {
        Vec3([0.0; 3])
    }
    #[inline]
    pub fn dot(self, other: Vec3) -> f64 {
        self.0[0] * other.0[0] + self.0[1] * other.0[1] + self.0[2] * other.0[2]
    }
    #[inline]
    pub fn norm(self) -> f64 {
        self.dot(self).sqrt()
    }
}

impl std::ops::Index<usize> for Vec3 {
    type Output = f64;
    #[inline]
    fn index(&self, i: usize) -> &f64 {
        &self.0[i]
    }
}

impl std::ops::Add for Vec3 {
    type Output = Vec3;
    #[inline]
    fn add(self, r: Vec3) -> Vec3 {
        Vec3([self.0[0] + r.0[0], self.0[1] + r.0[1], self.0[2] + r.0[2]])
    }
}
impl std::ops::Sub for Vec3 {
    type Output = Vec3;
    #[inline]
    fn sub(self, r: Vec3) -> Vec3 {
        Vec3([self.0[0] - r.0[0], self.0[1] - r.0[1], self.0[2] - r.0[2]])
    }
}
impl std::ops::Neg for Vec3 {
    type Output = Vec3;
    #[inline]
    fn neg(self) -> Vec3 {
        Vec3([-self.0[0], -self.0[1], -self.0[2]])
    }
}
impl std::ops::Mul<f64> for Vec3 {
    type Output = Vec3;
    #[inline]
    fn mul(self, s: f64) -> Vec3 {
        Vec3([self.0[0] * s, self.0[1] * s, self.0[2] * s])
    }
}
impl std::ops::Mul<Vec3> for f64 {
    type Output = Vec3;
    #[inline]
    fn mul(self, v: Vec3) -> Vec3 {
        v * self
    }
}

impl Mat3 {
    #[inline]
    pub const fn eye() -> Self {
        Mat3([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    }
    #[inline]
    pub fn transpose(&self) -> Mat3 {
        let m = &self.0;
        Mat3([
            [m[0][0], m[1][0], m[2][0]],
            [m[0][1], m[1][1], m[2][1]],
            [m[0][2], m[1][2], m[2][2]],
        ])
    }
}

impl std::ops::Mul<Vec3> for Mat3 {
    type Output = Vec3;
    #[inline]
    fn mul(self, v: Vec3) -> Vec3 {
        let m = &self.0;
        let x = &v.0;
        Vec3([
            m[0][0] * x[0] + m[0][1] * x[1] + m[0][2] * x[2],
            m[1][0] * x[0] + m[1][1] * x[1] + m[1][2] * x[2],
            m[2][0] * x[0] + m[2][1] * x[1] + m[2][2] * x[2],
        ])
    }
}

// ---------------------------------------------------------------------------
// Pure-Rust I/O structs. The Python-facing wrappers in the glue crate hold
// these by value and forward.
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
#[allow(non_snake_case)]
pub struct RotorInputs {
    pub collective_rad: f64,
    pub tilt_lon: f64,
    pub tilt_lat: f64,
    pub R_hub: Mat3,
    pub v_hub_world: Vec3,
    pub wind_world: Vec3,
    pub t: f64,
    pub rho_kg_m3: f64,
    pub motor_torque_Nm: f64,
}

#[derive(Clone, Debug)]
#[allow(non_snake_case)]
pub struct AeroResult {
    pub F_world: Vec3,
    pub M_orbital: Vec3,
    pub Q_spin: f64,
    pub M_spin: Vec3,
}
