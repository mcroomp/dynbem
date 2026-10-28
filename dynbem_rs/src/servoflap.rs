// Quasi-static 1/rev blade feathering response to a trailing-edge servo-flap.
//
// The time-domain model used by the VPM is
//
//   I_theta * theta_ddot + C_theta * theta_dot
//       + (I_theta * Omega^2 + k_aero) * theta = M_servo + M_camber.
//
// For a thin blade, I_theta * Omega^2 is the centrifugal restoring stiffness
// about its feathering axis. k_aero is the additional lift-arm stiffness when
// the aerodynamic centre is offset from that axis. Substituting
// theta = theta_0 + A*cos(psi) + B*sin(psi), psi_dot = Omega, gives the 1/rev
// matrix below. The centrifugal term cancels the 1/rev inertial term, but it
// remains in the DC trim equation.
//
// This BEM approximation keeps only the mean dynamic pressure and its linear
// first-harmonic perturbation from tangential in-plane wind. The VPM instead
// computes local dynamic pressure and advances the same physical feathering
// DOF in time. Neither model includes torsional blade flexibility.

use crate::aero_io::Vec3;
use crate::bem_common::RadialGrid;
use crate::rotor_definition::{ServoFlapActuation, ServoFlapGeometry};

#[derive(Clone, Copy, Debug, Default)]
struct FeatheringIntegrals {
    lift_arm_r2: f64,
    blade_moment_r2: f64,
    blade_moment_r1: f64,
    flap_moment_r2: f64,
    flap_moment_r1: f64,
}

impl FeatheringIntegrals {
    fn from_grid(grid: &RadialGrid, flap: &ServoFlapGeometry) -> Self {
        let mut integrals = Self::default();
        let half_dr = 0.5 * grid.dr;

        for i in 0..grid.n_elements {
            let r = grid.r_mid[i];
            let chord = grid.chord[i];
            let lo = r - half_dr;
            let hi = r + half_dr;
            let int_r2 = (hi.powi(3) - lo.powi(3)) / 3.0;
            let int_r1 = (hi.powi(2) - lo.powi(2)) / 2.0;
            integrals.lift_arm_r2 += chord * int_r2;
            integrals.blade_moment_r2 += chord * chord * int_r2;
            integrals.blade_moment_r1 += chord * chord * int_r1;

            let flap_lo = lo.max(flap.r_inner_m);
            let flap_hi = hi.min(flap.r_outer_m);
            if flap_hi > flap_lo {
                let flap_r2 = (flap_hi.powi(3) - flap_lo.powi(3)) / 3.0;
                let flap_r1 = (flap_hi.powi(2) - flap_lo.powi(2)) / 2.0;
                integrals.flap_moment_r2 += chord * chord * flap_r2;
                integrals.flap_moment_r1 += chord * chord * flap_r1;
            }
        }

        integrals
    }
}

#[inline]
pub(crate) fn centrifugal_pitch_stiffness(i_theta: f64, omega: f64) -> f64 {
    i_theta * omega * omega
}

#[inline]
pub(crate) fn aerodynamic_pitch_stiffness(
    rho: f64,
    omega: f64,
    cl_alpha: f64,
    ac_offset: f64,
    lift_arm_r2: f64,
) -> f64 {
    0.5 * rho * omega * omega * cl_alpha * ac_offset * lift_arm_r2
}

/// Solved feathering pitch harmonics for one compute_forces call.
///
/// In servo-flap mode these REPLACE the swashplate collective/cyclic pitch in
/// the psi-loop (the direct swashplate-to-pitch path is disabled there).
/// All angles in radians.  Zero when FeatheringProperties is None.
#[derive(Clone, Copy, Debug, Default)]
pub struct FeatheringState {
    /// Mean (DC) feathering offset [rad], driven by the collective command.
    pub delta_theta_0: f64,
    /// 1/rev longitudinal feathering harmonic [rad].
    /// delta_theta(psi) = delta_theta_1c*cos(psi) + delta_theta_1s*sin(psi)
    pub delta_theta_1c: f64,
    /// 1/rev lateral feathering harmonic [rad].
    pub delta_theta_1s: f64,
}

impl FeatheringState {
    /// Zero feathering -- used when FeatheringProperties is None (rigid pitch).
    pub const RIGID: Self = Self {
        delta_theta_0: 0.0,
        delta_theta_1c: 0.0,
        delta_theta_1s: 0.0,
    };
}

// ---------------------------------------------------------------------------
// Servo-flap pitching moment harmonics
// ---------------------------------------------------------------------------

/// Compute the aerodynamic pitching moment about the feathering axis from the
/// servo-flap, decomposed into mean + 1/rev harmonics (M_f0, M_f1c, M_f1s).
///
/// Integrates `dM/dr = q * chord^2 * C_M_delta * delta_f` over the flap span.
/// The first-harmonic dynamic-pressure terms use the project hub-frame
/// tangential wind components and are linearized about `Omega*r`.
fn servo_flap_moments(
    sp: &ServoFlapGeometry,
    delta_f0: f64,
    delta_f1c: f64,
    delta_f1s: f64,
    rho: f64,
    omega: f64,
    v_inplane_hub: Vec3,
    integrals: FeatheringIntegrals,
) -> (f64, f64, f64) {
    if omega == 0.0 || integrals.flap_moment_r2 == 0.0 {
        return (0.0, 0.0, 0.0);
    }

    let m_dc_shape = 0.5 * rho * omega * omega * sp.C_M_delta_per_rad * integrals.flap_moment_r2;
    let m_c_shape =
        rho * omega * sp.C_M_delta_per_rad * v_inplane_hub.0[1] * integrals.flap_moment_r1;
    let m_s_shape =
        rho * omega * sp.C_M_delta_per_rad * v_inplane_hub.0[0] * integrals.flap_moment_r1;

    let m_f0 = delta_f0 * m_dc_shape + 0.5 * (delta_f1c * m_c_shape + delta_f1s * m_s_shape);
    let m_f1c = delta_f1c * m_dc_shape + delta_f0 * m_c_shape;
    let m_f1s = delta_f1s * m_dc_shape + delta_f0 * m_s_shape;

    (m_f0, m_f1c, m_f1s)
}

// ---------------------------------------------------------------------------
// Main solve
// ---------------------------------------------------------------------------

/// Solve the quasi-static feathering harmonic balance for one compute_forces call.
///
/// Returns FeatheringState (delta_theta_0, delta_theta_1c, delta_theta_1s). In
/// servo-flap mode these REPLACE the swashplate collective/cyclic pitch in the
/// psi-loop; the direct swashplate-to-pitch path is disabled there.
///
/// # Arguments
/// * `act`     -- ServoFlapActuation (the rotor is in servo-flap mode)
/// * `theta_0` -- collective command [rad], mapped to DC flap command in
///                servo-flap mode
/// * `theta_1c/s` -- swashplate command harmonics [rad], passed through to
///                  flap harmonics without internal phase compensation.
/// * `rho`     -- air density [kg/m^3]
/// * `omega`   -- rotor angular velocity [rad/s]
/// * `v_inplane_hub` -- hub-frame in-plane relative-air velocity [m/s]
/// * `grid`    -- blade radial stations and per-station chord
/// * `cl_alpha` -- lift curve slope [1/rad] (for aerodynamic spring if ac_offset != 0)
pub fn solve_feathering(
    act: &ServoFlapActuation,
    theta_0: f64,
    theta_1c: f64,
    theta_1s: f64,
    rho: f64,
    omega: f64,
    v_inplane_hub: Vec3,
    grid: &RadialGrid,
    cl_alpha: f64,
) -> FeatheringState {
    if omega < 1e-6 {
        return FeatheringState::RIGID;
    }

    let i_theta = act.I_theta_kgm2;
    let c_theta = act.damper_Nms_per_rad;
    let integrals = FeatheringIntegrals::from_grid(grid, &act.flap);
    assert!(
        i_theta.is_finite() && i_theta > 0.0,
        "servo-flap pitch inertia must be finite and positive"
    );
    assert!(
        c_theta.is_finite() && c_theta >= 0.0,
        "servo-flap bearing damping must be finite and nonnegative"
    );

    // Non-dimensional mechanical damping: d_mech = C / (2 * I * Omega)
    let d_mech = c_theta / (2.0 * i_theta * omega);

    // Lift-arm spring about the feathering axis. The thin-blade centrifugal
    // spring is I_theta*Omega^2, so DC stiffness remains at zero AC offset.
    let omega2 = omega * omega;
    let k_aero =
        aerodynamic_pitch_stiffness(rho, omega, cl_alpha, act.ac_offset_m, integrals.lift_arm_r2);
    let k_centrifugal = centrifugal_pitch_stiffness(i_theta, omega);
    let k_total = k_centrifugal + k_aero;
    assert!(
        k_total.is_finite() && k_total != 0.0,
        "servo-flap has no finite static feathering trim at this rotor speed"
    );

    // Servo-flap moment harmonics (raw physical mapping, no pre-rotation).
    let (mut m_f0, mut m_f1c, mut m_f1s) = servo_flap_moments(
        &act.flap,
        theta_0,
        theta_1c,
        theta_1s,
        rho,
        omega,
        v_inplane_hub,
        integrals,
    );

    // Blade zero-lift pitching moment over the full blade, using the same
    // chord-squared moment-coefficient normalization as the flap moment.
    if act.blade_Cm_AC.abs() > 1e-12 {
        m_f0 += 0.5 * rho * omega2 * act.blade_Cm_AC * integrals.blade_moment_r2;
        m_f1c += rho * omega * act.blade_Cm_AC * v_inplane_hub.0[1] * integrals.blade_moment_r1;
        m_f1s += rho * omega * act.blade_Cm_AC * v_inplane_hub.0[0] * integrals.blade_moment_r1;
    }

    // Static trim balances the servo and camber moments against both
    // centrifugal and aerodynamic pitch stiffness.
    let delta_theta_0 = m_f0 / k_total;

    // Normalised RHS
    let inv_i_omega2 = 1.0 / (i_theta * omega2);
    let rhs_c = m_f1c * inv_i_omega2;
    let rhs_s = m_f1s * inv_i_omega2;

    // At 1/rev, the centrifugal stiffness I*Omega^2 cancels the inertial
    // term -I*Omega^2. The remaining diagonal is k_aero/(I*Omega^2).
    //   [diag,  2*d] [A] = [rhs_c]
    //   [-2*d, diag] [B]   [rhs_s]
    let diag = k_aero * inv_i_omega2;
    let det = diag * diag + 4.0 * d_mech * d_mech;

    assert!(
        det.is_finite() && det > 0.0,
        "servo-flap has an undamped 1/rev feathering resonance"
    );

    let inv_det = 1.0 / det;
    let delta_theta_1c = inv_det * (diag * rhs_c - 2.0 * d_mech * rhs_s);
    let delta_theta_1s = inv_det * (2.0 * d_mech * rhs_c + diag * rhs_s);

    FeatheringState {
        delta_theta_0,
        delta_theta_1c,
        delta_theta_1s,
    }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::rotor_definition::{BladeGeometry, ServoFlapActuation, ServoFlapGeometry};

    fn make_grid(chord: f64) -> RadialGrid {
        RadialGrid::from_blade(&BladeGeometry::uniform(2, 2.5, 0.0, chord, 0.0, 10))
    }

    fn make_actuation(c_theta: f64, ac_offset: f64) -> ServoFlapActuation {
        ServoFlapActuation {
            I_theta_kgm2: 0.05,
            damper_Nms_per_rad: c_theta,
            ac_offset_m: ac_offset,
            blade_Cm_AC: 0.0,
            flap: ServoFlapGeometry {
                C_M_delta_per_rad: -1.5,
                r_inner_m: 1.2,
                r_outer_m: 2.5,
            },
        }
    }

    // Zero command -> zero feathering
    #[test]
    fn test_no_command_zero_feathering() {
        let fp = make_actuation(5.0, 0.0);
        let s = solve_feathering(
            &fp,
            0.0,
            0.0,
            0.0,
            1.225,
            33.2,
            Vec3::zero(),
            &make_grid(0.20),
            5.79,
        );
        assert!(s.delta_theta_0.abs() < 1e-10);
        assert!(s.delta_theta_1c.abs() < 1e-10);
        assert!(s.delta_theta_1s.abs() < 1e-10);
    }

    // Lateral command (theta_1c) at AC (p_sq=1): cos forcing -> sin response.
    #[test]
    fn test_theta_1c_command_cross_couples_to_theta_1s_at_ac() {
        let fp = make_actuation(5.0, 0.0);
        let grid = make_grid(0.20);
        let s = solve_feathering(&fp, 0.0, 0.1, 0.0, 1.225, 33.2, Vec3::zero(), &grid, 5.79);
        let integrals = FeatheringIntegrals::from_grid(&grid, &fp.flap);
        let moment = 0.5
            * 1.225
            * 33.2_f64.powi(2)
            * fp.flap.C_M_delta_per_rad
            * integrals.flap_moment_r2
            * 0.1;
        assert!(s.delta_theta_1s.abs() > 1e-4);
        assert!((s.delta_theta_1s - moment / (fp.damper_Nms_per_rad * 33.2)).abs() < 1e-12);
        assert!(
            s.delta_theta_1c.abs() < 1e-6,
            "delta_theta_1c should be ~0 for pure theta_1c command: {}",
            s.delta_theta_1c
        );
    }

    // Omega=0 -> RIGID
    #[test]
    fn test_zero_omega_returns_rigid() {
        let fp = make_actuation(5.0, 0.0);
        let s = solve_feathering(
            &fp,
            0.0,
            0.1,
            0.1,
            1.225,
            0.0,
            Vec3::zero(),
            &make_grid(0.20),
            5.79,
        );
        assert_eq!(s.delta_theta_1c, 0.0);
        assert_eq!(s.delta_theta_1s, 0.0);
    }

    // Authority scales with 1/C_theta (less damping -> more feathering)
    #[test]
    fn test_authority_scales_with_damper() {
        let fp_lo = make_actuation(2.0, 0.0);
        let fp_hi = make_actuation(8.0, 0.0);
        let grid = make_grid(0.20);
        let s_lo = solve_feathering(
            &fp_lo,
            0.0,
            0.1,
            0.0,
            1.225,
            33.2,
            Vec3::zero(),
            &grid,
            5.79,
        );
        let s_hi = solve_feathering(
            &fp_hi,
            0.0,
            0.1,
            0.0,
            1.225,
            33.2,
            Vec3::zero(),
            &grid,
            5.79,
        );
        // Lower damper -> larger response
        assert!(
            s_lo.delta_theta_1s.abs() > s_hi.delta_theta_1s.abs(),
            "s_lo={}, s_hi={}",
            s_lo.delta_theta_1s,
            s_hi.delta_theta_1s
        );
    }

    // AC offset creates a spring: larger p_sq reduces authority at 1/rev
    #[test]
    fn test_ac_offset_reduces_authority() {
        let fp_at_ac = make_actuation(5.0, 0.0);
        let fp_offset = make_actuation(5.0, 0.01); // 1 cm forward of AC
        let grid = make_grid(0.20);
        let s_at_ac = solve_feathering(
            &fp_at_ac,
            0.0,
            0.1,
            0.0,
            1.225,
            33.2,
            Vec3::zero(),
            &grid,
            5.79,
        );
        let s_offset = solve_feathering(
            &fp_offset,
            0.0,
            0.1,
            0.0,
            1.225,
            33.2,
            Vec3::zero(),
            &grid,
            5.79,
        );
        // AC offset moves p_sq away from 1, det changes
        assert!(s_at_ac.delta_theta_1s.abs() != s_offset.delta_theta_1s.abs());
    }

    #[test]
    fn test_nonzero_offset_dc_balance_is_independent_of_omega() {
        let theta_0 = 0.001;
        let cl_alpha = 5.75;
        let grid = make_grid(0.10);

        for ac_offset in [-0.01, 0.01] {
            let fp = make_actuation(0.5, ac_offset);
            let integrals = FeatheringIntegrals::from_grid(&grid, &fp.flap);
            let expected =
                0.5 * 1.225 * fp.flap.C_M_delta_per_rad * theta_0 * integrals.flap_moment_r2
                    / (fp.I_theta_kgm2
                        + 0.5 * 1.225 * cl_alpha * ac_offset * integrals.lift_arm_r2);

            for omega in [14.135, 28.27, 56.54] {
                let s = solve_feathering(
                    &fp,
                    theta_0,
                    0.0,
                    0.0,
                    1.225,
                    omega,
                    Vec3::zero(),
                    &grid,
                    cl_alpha,
                );
                assert!(
                    (s.delta_theta_0 - expected).abs() < 1e-12,
                    "offset={ac_offset}, omega={omega}: theta_0={}, expected={expected}",
                    s.delta_theta_0,
                );
            }
        }
    }

    #[test]
    fn test_dc_collective_uses_centrifugal_stiffness_not_pass_through() {
        let grid = make_grid(0.10);
        let fp = make_actuation(5.0, 0.0);
        let command = 0.07;
        let integrals = FeatheringIntegrals::from_grid(&grid, &fp.flap);
        let expected = 0.5 * 1.225 * fp.flap.C_M_delta_per_rad * command * integrals.flap_moment_r2
            / fp.I_theta_kgm2;
        let response = solve_feathering(
            &fp,
            command,
            0.0,
            0.0,
            1.225,
            33.2,
            Vec3::zero(),
            &grid,
            5.79,
        );
        assert!((response.delta_theta_0 - expected).abs() < 1e-12);

        let mut no_flap_moment = fp;
        no_flap_moment.flap.C_M_delta_per_rad = 0.0;
        let response = solve_feathering(
            &no_flap_moment,
            command,
            0.0,
            0.0,
            1.225,
            33.2,
            Vec3::zero(),
            &grid,
            5.79,
        );
        assert_eq!(response.delta_theta_0, 0.0);
    }

    #[test]
    fn test_flap_moment_has_chord_squared_scaling() {
        let fp = make_actuation(5.0, 0.0);
        let grid = make_grid(0.10);
        let integrals = FeatheringIntegrals::from_grid(&grid, &fp.flap);
        let expected_integral =
            0.10_f64.powi(2) * (fp.flap.r_outer_m.powi(3) - fp.flap.r_inner_m.powi(3)) / 3.0;
        assert!((integrals.flap_moment_r2 - expected_integral).abs() < 1e-12);
        let command = 0.01;
        let moment = servo_flap_moments(
            &fp.flap,
            command,
            0.0,
            0.0,
            1.225,
            33.2,
            Vec3::zero(),
            integrals,
        )
        .0;
        let expected = 0.5
            * 1.225
            * 33.2_f64.powi(2)
            * command
            * fp.flap.C_M_delta_per_rad
            * integrals.flap_moment_r2;
        assert!((moment - expected).abs() < 1e-12);
        assert!(moment < 0.0);
    }

    #[test]
    fn test_forward_flight_harmonics_follow_hub_wind_direction() {
        let fp = make_actuation(5.0, 0.0);
        let grid = make_grid(0.10);
        let x_wind = solve_feathering(
            &fp,
            0.001,
            0.0,
            0.0,
            1.225,
            33.2,
            Vec3::new(5.0, 0.0, 0.0),
            &grid,
            5.79,
        );
        let y_wind = solve_feathering(
            &fp,
            0.001,
            0.0,
            0.0,
            1.225,
            33.2,
            Vec3::new(0.0, 5.0, 0.0),
            &grid,
            5.79,
        );
        assert!(x_wind.delta_theta_1c.abs() > 0.0);
        assert_eq!(x_wind.delta_theta_1s, 0.0);
        assert!(y_wind.delta_theta_1s.abs() > 0.0);
        assert_eq!(y_wind.delta_theta_1c, 0.0);
    }
}
