// Cyclic pitch / swashplate mapping. See ../CLAUDE.md "Cyclic pitch convention".
//
//   theta(psi) = collective + theta_1c*cos(psi) + theta_1s*sin(psi)
//
// With helicopter-standard signs (tilt_lon>0 -> nose-down, tilt_lat>0 -> roll right):
//   theta_1c = gain * (-tilt_lon*cos(phi) - tilt_lat*sin(phi))
//   theta_1s = gain * (-tilt_lon*sin(phi) + tilt_lat*cos(phi))
// phi defaults to 0; gain defaults to 1.

#[derive(Clone, Copy, Debug)]
pub struct ControlGains {
    pub gain: f64,
    pub phase_rad: f64,
}

impl Default for ControlGains {
    fn default() -> Self {
        Self {
            gain: 1.0,
            phase_rad: 0.0,
        }
    }
}

#[inline]
pub fn cyclic_coeffs(tilt_lon: f64, tilt_lat: f64, ctrl: ControlGains) -> (f64, f64) {
    let cos_phi = ctrl.phase_rad.cos();
    let sin_phi = ctrl.phase_rad.sin();
    let theta_1c = ctrl.gain * (-tilt_lon * cos_phi - tilt_lat * sin_phi);
    let theta_1s = ctrl.gain * (-tilt_lon * sin_phi + tilt_lat * cos_phi);
    (theta_1c, theta_1s)
}
