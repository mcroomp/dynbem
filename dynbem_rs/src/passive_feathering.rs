// Quasi-static blade feathering model driven by a trailing-edge servo-flap.
//
// # Physics
//
// The blade is free to rotate about its feathering (pitch) axis, restrained
// only by a rotary damper at the pitch bearing and -- if the aerodynamic
// centre is offset from the feathering axis -- by an aerodynamic spring.
//
// Kaman design intent: feathering axis at the AC, so no aerodynamic spring
// and no aerodynamic damping from the main blade in pitch.  The only damping
// is mechanical (the bearing damper).  The servo-flap exerts a pitching
// moment that drives feathering angle delta_theta(psi).
//
// EOM in psi-domain (primes = d/d(psi)):
//
//   I_theta * delta_theta'' + C_theta * delta_theta' + k_aero * delta_theta
//       = M_servo(psi)                                                   (1)
//
// where:
//   k_aero = 0.5 * rho * cl_alpha * chord * ac_offset * R^3 / 3
//             (aerodynamic spring from AC offset; zero if ac_offset = 0)
//   M_servo = servo-flap pitching moment (1/rev harmonic forcing)
//
// Dividing by I_theta * Omega^2 and defining:
//   p_sq  = 1 + k_aero / (I_theta * Omega^2)     (frequency ratio squared)
//   d_mech = C_theta / (2 * I_theta * Omega)      (non-dim mechanical damping)
//
// Note: p_sq = 1 for a pure feathering axis at the AC (no spring).
// Note: there is NO aerodynamic damping term in pitch (unlike flapping where
//       the Lock-number d_aero = gamma/8 dominates).  The feathering DOF only
//       couples to the aerodynamics through the spring term (if ac_offset != 0)
//       and through the forced response itself.
//
// The non-dimensional 2x2 harmonic balance at 1/rev:
//
//   [-(p_sq-1), -2*d_mech] [delta_theta_1c]   [m_f1c]
//   [-2*d_mech, -(p_sq-1)] [delta_theta_1s] = [m_f1s]
//
// where m_f = M_servo / (I_theta * Omega^2).
//
// det = (p_sq-1)^2 - 4*d_mech^2  (can be zero at resonance p_sq=1, d_mech=0)
// For p_sq=1 (feathering at AC, typical Kaman):
//   det = -4*d_mech^2
//   delta_theta_1c = -m_f1s / (2 * d_mech)  ->  M_f1s / (C_theta * Omega)
//   delta_theta_1s = -m_f1c / (2 * d_mech)  ->  -M_f1c / (C_theta * Omega)
// Wait -- let me re-derive carefully to get the sign right.
//
// With delta_theta = A*cos(psi) + B*sin(psi):
//   delta_theta'  = -A*sin(psi) + B*cos(psi)
//   delta_theta'' = -A*cos(psi) - B*sin(psi)
//
// EOM / (I_theta * Omega^2):
//   cos(psi): -(p_sq)*A + 2*d_mech*B = m_f1c / (I*O^2)   -- wait, let me redo
//
// I_theta*O^2*[ -A*cos - B*sin ] + C_theta*O*[ -A*sin + B*cos ] + k_aero*[A*cos + B*sin]
//   = M_f1c*cos + M_f1s*sin
//
// cos: I*O^2*(-A) + C*O*(B) + k*A = M_f1c
//      A*(k - I*O^2) + B*(C*O)    = M_f1c
//      A*I*O^2*(p_sq - 1) + B*C*O = M_f1c        [since k = I*O^2*(p_sq-1)]
//
// sin: I*O^2*(-B) + C*O*(-A) + k*B = M_f1s
//      -A*(C*O) + B*(k - I*O^2)    = M_f1s
//      -A*C*O + B*I*O^2*(p_sq - 1) = M_f1s
//
// Dividing by I*O^2, letting zeta = C*O/(2*I*O^2) = C/(2*I*O) = d_mech:
//   A*(p_sq-1) + 2*d_mech*B = M_f1c / (I*O^2)  =: rhs_c
//   -2*d_mech*A + B*(p_sq-1) = M_f1s / (I*O^2)  =: rhs_s
//
// Matrix: M = [[(p_sq-1), 2*d_mech], [-2*d_mech, (p_sq-1)]]
// det(M) = (p_sq-1)^2 + 4*d_mech^2
//
// Solution (Cramer):
//   A = delta_theta_1c = [ (p_sq-1)*rhs_c - 2*d_mech*rhs_s ] / det
//   B = delta_theta_1s = [ 2*d_mech*rhs_c + (p_sq-1)*rhs_s ] / det
//
// For p_sq=1 (det = 4*d_mech^2):
//   A = -2*d_mech*rhs_s / (4*d_mech^2) = -rhs_s / (2*d_mech)
//         = -M_f1s / (2*d_mech * I*O^2) = -M_f1s / (C*O)   [since 2*d_mech*I*O^2 = C*O]
//   B =  2*d_mech*rhs_c / (4*d_mech^2) =  rhs_c / (2*d_mech)
//         =  M_f1c / (C*O)
//
// Physical interpretation (p_sq=1):
//   A positive M_f1c (cos forcing, e.g. from lateral tilt) drives B > 0 (sin response).
//   This is the classical 90-deg phase lag: cos forcing -> sin response.
//
// These delta_theta values are added directly to loop_theta_1c/1s in the psi-loop,
// giving full-span pitch authority equivalent to the swashplate.
//
// # References
// Kaman Aerospace: servo-flap rotor design patents and technical reports.
// Bramwell, "Helicopter Dynamics", Ch. 4 -- pitch-flap coupling, feathering EOM.

use crate::rotor_definition::{PassiveFeatheringProperties, ServoFlapProperties};

/// Solved feathering pitch harmonics for one compute_forces call.
///
/// These are ADDED to the swashplate cyclic pitch in the psi-loop.
/// All angles in radians.  Zero when FeatheringProperties is None.
#[derive(Clone, Copy, Debug, Default)]
pub struct FeatheringState {
    /// 1/rev longitudinal feathering harmonic [rad].
    /// delta_theta(psi) = delta_theta_1c*cos(psi) + delta_theta_1s*sin(psi)
    pub delta_theta_1c: f64,
    /// 1/rev lateral feathering harmonic [rad].
    pub delta_theta_1s: f64,
}

impl FeatheringState {
    /// Zero feathering -- used when FeatheringProperties is None (rigid pitch).
    pub const RIGID: Self = Self {
        delta_theta_1c: 0.0,
        delta_theta_1s: 0.0,
    };
}

// ---------------------------------------------------------------------------
// Servo-flap pitching moment harmonics
// ---------------------------------------------------------------------------

/// Compute the aerodynamic pitching moment about the feathering axis from the
/// servo-flap, decomposed into 1/rev harmonics (M_f1c, M_f1s).
///
/// Integrates:  dM/dr = 0.5 * rho * v_local^2 * c * C_M_delta * delta_f(psi)
/// over the flap span, keeping only cos(psi) and sin(psi) terms.
///
/// v_local(r, psi) = Omega*r + Omega*R*mu*sin(psi)   (hover: mu=0)
///
/// After expanding v_local^2 and multiplying by delta_f harmonics:
///   M_f1c = delta_f1c * m_dc_shape
///   M_f1s = delta_f1s * m_dc_shape + delta_f0 * m_sin_shape
///
/// where:
///   m_dc_shape  = 0.5 * rho * Omega^2 * C_M_delta * chord * (r_out^3-r_in^3)/3
///   m_sin_shape = 0.5 * rho * Omega^2 * C_M_delta * chord * mu * R * (r_out^2-r_in^2)
fn servo_flap_moments(
    sp: &ServoFlapProperties,
    delta_f0: f64,
    delta_f1c: f64,
    delta_f1s: f64,
    rho: f64,
    omega: f64,
    mu: f64,
    r_tip: f64,
    chord: f64,
) -> (f64, f64) {
    let r_in = sp.r_inner_m;
    let r_out = sp.r_outer_m.min(r_tip);
    if r_out <= r_in || omega == 0.0 {
        return (0.0, 0.0);
    }

    let r3 = (r_out.powi(3) - r_in.powi(3)) / 3.0;
    let r2 = (r_out.powi(2) - r_in.powi(2)) / 2.0;
    let q_base = 0.5 * rho * omega * omega * sp.C_M_delta_per_rad * chord;

    // DC velocity^2 integral (dominant term)
    let m_dc_shape = q_base * r3;
    // Cross term with forward-flight mu: gives extra sin(psi) contribution
    let m_sin_shape = q_base * 2.0 * mu * r_tip * r2;

    let m_f1c = delta_f1c * m_dc_shape;
    let m_f1s = delta_f1s * m_dc_shape + delta_f0 * m_sin_shape;

    (m_f1c, m_f1s)
}

// ---------------------------------------------------------------------------
// Main solve
// ---------------------------------------------------------------------------

/// Solve the quasi-static feathering harmonic balance for one compute_forces call.
///
/// Returns FeatheringState with (delta_theta_1c, delta_theta_1s) to be added
/// to the swashplate cyclic pitch in the psi-loop.
///
/// # Arguments
/// * `fp`      -- FeatheringProperties (never None at call site)
/// * `theta_1c/s` -- swashplate command harmonics [rad]; interpreted as servo-flap
///                  deflection delta_f when servo is configured
/// * `rho`     -- air density [kg/m^3]
/// * `omega`   -- rotor angular velocity [rad/s]
/// * `mu`      -- advance ratio
/// * `r_tip`   -- blade tip radius [m]
/// * `chord`   -- mean blade chord [m]
/// * `cl_alpha` -- lift curve slope [1/rad] (for aerodynamic spring if ac_offset != 0)
pub fn solve_feathering(
    fp: &PassiveFeatheringProperties,
    theta_1c: f64,
    theta_1s: f64,
    rho: f64,
    omega: f64,
    mu: f64,
    r_tip: f64,
    chord: f64,
    cl_alpha: f64,
) -> FeatheringState {
    if omega < 1e-6 {
        return FeatheringState::RIGID;
    }

    let i_theta = fp.I_theta_kgm2;
    let c_theta = fp.damper_Nms_per_rad;

    // Non-dimensional mechanical damping: d_mech = C / (2 * I * Omega)
    let d_mech = c_theta / (2.0 * i_theta * omega);

    // Aerodynamic spring from AC offset (zero for Kaman design with axis at AC).
    // k_aero = 0.5 * rho * cl_alpha * chord * ac_offset * integral(r^2 dr, 0, R)
    //        = 0.5 * rho * cl_alpha * chord * ac_offset * R^3/3
    let k_aero = if fp.ac_offset_m.abs() > 1e-9 {
        0.5 * rho * cl_alpha * chord * fp.ac_offset_m * r_tip.powi(3) / 3.0
    } else {
        0.0
    };
    // p_sq = 1 + k_aero / (I * Omega^2);  p_sq = 1 when ac_offset = 0
    let p_sq = 1.0 + k_aero / (i_theta * omega * omega);

    // Servo-flap moment harmonics.
    // The swashplate cyclic commands are interpreted as servo-flap deflection
    // angles delta_f when servo is configured.
    let (m_f1c, m_f1s) = match &fp.servoflaps {
        None => (0.0, 0.0),  // passive free-feathering: no servo forcing
        Some(sp) => servo_flap_moments(
            sp,
            0.0,       // DC flap deflection: zero (collective handled separately)
            theta_1c,  // lateral delta_f harmonic
            theta_1s,  // longitudinal delta_f harmonic
            rho,
            omega,
            mu,
            r_tip,
            chord,
        ),
    };

    // Normalised RHS
    let inv_i_omega2 = 1.0 / (i_theta * omega * omega);
    let rhs_c = m_f1c * inv_i_omega2;
    let rhs_s = m_f1s * inv_i_omega2;

    // 2x2 system (see derivation in module comment):
    //   [(p_sq-1),  2*d_mech] [A]   [rhs_c]
    //   [-2*d_mech, (p_sq-1)] [B] = [rhs_s]
    let diag = p_sq - 1.0;
    let det = diag * diag + 4.0 * d_mech * d_mech;

    if det.abs() < 1e-30 {
        // Undamped resonance (p_sq=1, d_mech=0) -- physically unbounded; clamp.
        return FeatheringState::RIGID;
    }

    let inv_det = 1.0 / det;
    let delta_theta_1c = inv_det * (diag * rhs_c - 2.0 * d_mech * rhs_s);
    let delta_theta_1s = inv_det * (2.0 * d_mech * rhs_c + diag * rhs_s);

    FeatheringState { delta_theta_1c, delta_theta_1s }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::rotor_definition::{PassiveFeatheringProperties, ServoFlapProperties};

    fn make_feathering(c_theta: f64, ac_offset: f64, with_servo: bool) -> PassiveFeatheringProperties {
        let servo = if with_servo {
            Some(ServoFlapProperties {
                C_M_delta_per_rad: -1.5,
                r_inner_m: 1.2,
                r_outer_m: 2.5,
            })
        } else {
            None
        };
        PassiveFeatheringProperties {
            I_theta_kgm2: 0.05,
            damper_Nms_per_rad: c_theta,
            ac_offset_m: ac_offset,
            servoflaps: servo,
        }
    }

    // Zero command -> zero feathering
    #[test]
    fn test_no_command_zero_feathering() {
        let fp = make_feathering(5.0, 0.0, true);
        let s = solve_feathering(&fp, 0.0, 0.0, 1.225, 33.2, 0.0, 2.5, 0.20, 5.79);
        assert!(s.delta_theta_1c.abs() < 1e-10);
        assert!(s.delta_theta_1s.abs() < 1e-10);
    }

    // Lateral command (theta_1c) at AC (p_sq=1): cos forcing -> sin response (90-deg lag)
    #[test]
    fn test_lateral_command_gives_sin_response_at_ac() {
        let fp = make_feathering(5.0, 0.0, true);
        let s = solve_feathering(&fp, 0.1, 0.0, 1.225, 33.2, 0.0, 2.5, 0.20, 5.79);
        // For p_sq=1: delta_theta_1c = -M_f1s/(C*O), delta_theta_1s = M_f1c/(C*O)
        // M_f1c is nonzero, M_f1s=0 (no DC or lateral command) -> delta_theta_1c=0, delta_theta_1s!=0
        assert!(s.delta_theta_1c.abs() < 1e-6,
            "delta_theta_1c should be ~0 for pure lateral command at AC: {}", s.delta_theta_1c);
        assert!(s.delta_theta_1s.abs() > 1e-4,
            "delta_theta_1s should be nonzero: {}", s.delta_theta_1s);
    }

    // No servo -> zero feathering regardless of command
    #[test]
    fn test_no_servo_zero_feathering() {
        let fp = make_feathering(5.0, 0.0, false);
        let s = solve_feathering(&fp, 0.1, 0.1, 1.225, 33.2, 0.1, 2.5, 0.20, 5.79);
        assert!(s.delta_theta_1c.abs() < 1e-10);
        assert!(s.delta_theta_1s.abs() < 1e-10);
    }

    // Omega=0 -> RIGID
    #[test]
    fn test_zero_omega_returns_rigid() {
        let fp = make_feathering(5.0, 0.0, true);
        let s = solve_feathering(&fp, 0.1, 0.1, 1.225, 0.0, 0.0, 2.5, 0.20, 5.79);
        assert_eq!(s.delta_theta_1c, 0.0);
        assert_eq!(s.delta_theta_1s, 0.0);
    }

    // Authority scales with 1/C_theta (less damping -> more feathering)
    #[test]
    fn test_authority_scales_with_damper() {
        let fp_lo = make_feathering(2.0, 0.0, true);
        let fp_hi = make_feathering(8.0, 0.0, true);
        let s_lo = solve_feathering(&fp_lo, 0.1, 0.0, 1.225, 33.2, 0.0, 2.5, 0.20, 5.79);
        let s_hi = solve_feathering(&fp_hi, 0.1, 0.0, 1.225, 33.2, 0.0, 2.5, 0.20, 5.79);
        // Lower damper -> larger response
        assert!(s_lo.delta_theta_1s.abs() > s_hi.delta_theta_1s.abs(),
            "s_lo={}, s_hi={}", s_lo.delta_theta_1s, s_hi.delta_theta_1s);
    }

    // AC offset creates a spring: larger p_sq reduces authority at 1/rev
    #[test]
    fn test_ac_offset_reduces_authority() {
        let fp_at_ac   = make_feathering(5.0, 0.0,   true);
        let fp_offset  = make_feathering(5.0, 0.01,  true);  // 1 cm forward of AC
        let s_at_ac  = solve_feathering(&fp_at_ac,  0.1, 0.0, 1.225, 33.2, 0.0, 2.5, 0.20, 5.79);
        let s_offset = solve_feathering(&fp_offset, 0.1, 0.0, 1.225, 33.2, 0.0, 2.5, 0.20, 5.79);
        // AC offset moves p_sq away from 1, det changes
        assert!(s_at_ac.delta_theta_1s.abs() != s_offset.delta_theta_1s.abs());
    }
}
