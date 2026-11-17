// Directional checks for servo-flap (Kaman-style) feathering dynamics.
//
// Moved from dynbem_rs/src/vpm/mod.rs::tests. Three checks:
//   1. Zero swashplate command -> zero feathering on all blades.
//   2. Collective command drives feathering to a bounded, settled angle.
//   3. Cyclic command changes the hub-moment vector from the no-command
//      baseline. Its expected axis and sign are checked separately in
//      cyclic_phase_servo.
//
// The bounded/settled DC feathering in check 2 requires a restoring spring.
// This is the aerodynamic spring from a nonzero ac_offset (AC aft of the
// feathering axis) -- a physical, measurable restoring moment.

use crate::helpers::*;
use crate::report::Report;
use dynbem_rs::rotor_definition::{PitchActuation, ServoFlapActuation, ServoFlapGeometry};

fn servo_actuation(ac_offset_m: f64) -> ServoFlapActuation {
    ServoFlapActuation {
        I_theta_kgm2: 0.02,
        damper_Nms_per_rad: 0.5,
        ac_offset_m,
        blade_Cm_AC: 0.0,
        flap: ServoFlapGeometry {
            C_M_delta_per_rad: -1.0,
            r_inner_m: R_ROOT,
            r_outer_m: 0.9 * R_TIP,
        },
    }
}

pub fn check_servo_flap(r: &mut Report) {
    r.begin_module(
        "servo_flap",
        "Directional: Kaman servo-flap feathering dynamics in VPM",
    );

    let mut defn = theory_rotor(10, 2.0);
    defn.pitch_actuation = PitchActuation::ServoFlap(servo_actuation(0.10));
    let rotor = make_fast_rotor(&defn);
    let dt = 1.0 / (OMEGA * 2.0 / std::f64::consts::PI);

    // --- 1. Zero command -> zero feathering ---
    use dynbem_rs::vpm::FlightCondition;
    let fc_zero = FlightCondition {
        collective_rad: 0.0,
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        v_hub: [0.0, 0.0, 0.0],
        omega_rad_s: OMEGA,
        rho: RHO,
    };
    let (_res_z, state_z) = rotor.march(&fc_zero, None, dt, 13 * STEPS_PER_REV);
    let theta_f_z = state_z.theta_f.expect("servo active -> theta_f present");
    let max_theta = theta_f_z
        .iter()
        .cloned()
        .fold(0.0_f64, |a, v| a.max(v.abs()));
    r.info("zero_cmd", "max_theta_f_rad", max_theta, 0.0);
    r.assert_bool(
        "zero_cmd",
        "feathering_stays_zero",
        max_theta,
        0.0,
        max_theta < 1e-6,
        &format!("zero command -> zero feathering, max={:.2e} rad", max_theta),
    );

    // --- 2. Collective command -> feathering settles to bounded nonzero angle ---
    let fc_col = hover_fc(8.0);
    let (_res_c, state_c) = rotor.march(&fc_col, None, dt, 26 * STEPS_PER_REV);
    let theta_f = state_c.theta_f.expect("servo active -> theta_f present");
    let theta_f_dot = state_c
        .theta_f_dot
        .expect("servo active -> theta_f_dot present");
    for (i, (&t, &td)) in theta_f.iter().zip(&theta_f_dot).enumerate() {
        let blade = format!("collective_blade{i}");
        r.info(&blade, "theta_f_rad", t, f64::NAN);
        r.info(&blade, "theta_f_dot_rad_per_s", td, f64::NAN);
        r.assert_bool(
            &blade,
            "feathering_nonzero",
            t.abs(),
            1e-3,
            t.abs() > 1e-3,
            &format!("collective should drive feathering, got {:.4} rad", t),
        );
        r.assert_bool(
            &blade,
            "feathering_bounded",
            t.abs(),
            0.7,
            t.abs() < 0.7,
            &format!("feathering implausibly large: {:.3} rad", t),
        );
        r.assert_bool(
            &blade,
            "feathering_settled",
            td.abs(),
            2.0,
            td.abs() < 2.0,
            &format!(
                "feathering should settle (theta_dot~0), got {:.3} rad/s",
                td
            ),
        );
    }

    // At zero AC offset there is no aerodynamic spring, but the rotating
    // blade's centrifugal stiffness still gives a finite DC equilibrium.
    let mut defn_centered = theory_rotor(10, 2.0);
    defn_centered.pitch_actuation = PitchActuation::ServoFlap(servo_actuation(0.0));
    let rotor_centered = make_fast_rotor(&defn_centered);
    let (_res_centered, state_centered) =
        rotor_centered.march(&fc_col, None, dt, 13 * STEPS_PER_REV);
    let theta_centered = state_centered
        .theta_f
        .expect("servo active -> theta_f present");
    let theta_dot_centered = state_centered
        .theta_f_dot
        .expect("servo active -> theta_f_dot present");
    for (i, (&theta, &theta_dot)) in theta_centered.iter().zip(&theta_dot_centered).enumerate() {
        r.assert_bool(
            &format!("centered_blade{i}"),
            "finite_dc_trim",
            theta.abs(),
            0.7,
            theta.is_finite() && theta.abs() > 1e-3 && theta.abs() < 0.7,
            &format!(
                "centrifugal stiffness should give bounded DC trim at zero offset, got {theta:.4} rad"
            ),
        );
        r.assert_bool(
            &format!("centered_blade{i}"),
            "dc_trim_settled",
            theta_dot.abs(),
            2.0,
            theta_dot.is_finite() && theta_dot.abs() < 2.0,
            &format!("zero-offset DC trim should settle, got {theta_dot:.3} rad/s"),
        );
    }

    // --- 3. Cyclic command -> changes the hub-moment vector ---
    let n_steps = 21 * STEPS_PER_REV;
    let (no_cyc, _) = rotor.march(&hover_fc(8.0), None, dt, n_steps);
    let mut fc_cyc = hover_fc(8.0);
    fc_cyc.tilt_lon = 4.0_f64.to_radians();
    let (with_cyc, _) = rotor.march(&fc_cyc, None, dt, n_steps);

    let m_none = no_cyc.my_hub.hypot(no_cyc.mx_hub);
    let m_cyc = with_cyc.my_hub.hypot(with_cyc.mx_hub);
    let m_delta = (with_cyc.my_hub - no_cyc.my_hub).hypot(with_cyc.mx_hub - no_cyc.mx_hub);
    r.info("cyclic_cmd", "hub_moment_no_cyc_Nm", m_none, f64::NAN);
    r.info("cyclic_cmd", "hub_moment_with_cyc_Nm", m_cyc, f64::NAN);
    r.info("cyclic_cmd", "hub_moment_change_Nm", m_delta, f64::NAN);
    r.assert_bool(
        "cyclic_cmd",
        "hub_moment_changed",
        m_delta,
        1e-3,
        m_delta > 1e-3,
        &format!(
            "cyclic feathering should change the hub-moment vector by more than 1e-3 N*m, got {:.3}",
            m_delta
        ),
    );
}
