// Directional checks for blade flap dynamics (coning and hub-moment relief).
//
// Moved from dynbem_rs/src/vpm/mod.rs::tests. Two checks:
//   1. With flap dynamics on, blades cone up (beta > 0) in hover and the
//      coning is equal across blades (axisymmetry).
//   2. A freely-hinged blade transmits less hub moment under cyclic than a
//      rigid blade (flap relief).

use crate::helpers::*;
use crate::report::Report;

pub fn check_flap_directional(r: &mut Report) {
    r.begin_module(
        "flap_directional",
        "Directional: VPM blade flap coning in hover and hub-moment relief under cyclic",
    );

    const I_BETA: f64 = 0.1; // freely hinged, omega_NR = 0
    let dt = 1.0 / (OMEGA * 2.0 / std::f64::consts::PI);
    let n_hover = 21 * STEPS_PER_REV;

    // --- 1. Hover coning: beta > 0, spread < 0.02 rad ---
    let defn_flap = theory_rotor_flap(10, I_BETA);
    let rotor_flap = make_fast_rotor(&defn_flap);
    let fc = hover_fc(8.0);
    let (_res, state) = rotor_flap.march(&fc, None, dt, n_hover);
    let beta = state.beta.expect("flap active -> beta present");

    let hi = beta.iter().cloned().fold(f64::MIN, f64::max);
    let lo = beta.iter().cloned().fold(f64::MAX, f64::min);
    r.info("hover_coning", "beta_min_rad", lo, f64::NAN);
    r.info("hover_coning", "beta_max_rad", hi, f64::NAN);
    r.info("hover_coning", "spread_rad", hi - lo, f64::NAN);
    r.assert_bool(
        "hover_coning",
        "coning_positive",
        lo,
        0.0,
        lo > 0.0,
        &format!("all blades should cone up (beta>0), min={:.4} rad", lo),
    );
    r.assert_bool(
        "hover_coning",
        "coning_not_implausible",
        hi,
        0.35,
        hi < 0.35,
        &format!("coning implausibly large: {:.3} rad", hi),
    );
    r.assert_bool(
        "hover_coning",
        "coning_axisymmetric",
        hi - lo,
        0.02,
        hi - lo < 0.02,
        &format!(
            "hover coning should be ~equal per blade, spread {:.4} rad",
            hi - lo
        ),
    );

    // --- 2. Flap relief: hinged rotor has smaller |M_hub| than rigid under cyclic ---
    let n_cyclic = 21 * STEPS_PER_REV;
    let mut fc_cyc = hover_fc(8.0);
    fc_cyc.tilt_lon = 4.0_f64.to_radians();

    // Rigid rotor (no flap properties).
    let defn_rigid = theory_rotor(10, 2.0);
    let rotor_rigid = make_fast_rotor(&defn_rigid);
    let (rigid, _) = rotor_rigid.march(&fc_cyc, None, dt, n_cyclic);

    // Hinged rotor.
    let (flex, _) = rotor_flap.march(&fc_cyc, None, dt, n_cyclic);

    let m_rigid = rigid.my_hub.hypot(rigid.mx_hub);
    let m_flex = flex.my_hub.hypot(flex.mx_hub);
    r.info("flap_relief", "m_rigid_Nm", m_rigid, f64::NAN);
    r.info("flap_relief", "m_flex_Nm", m_flex, f64::NAN);
    r.assert_bool(
        "flap_relief",
        "flex_lt_rigid",
        m_flex,
        m_rigid,
        m_flex < m_rigid,
        &format!(
            "flapping should reduce hub moment: rigid {:.3} vs flex {:.3} N*m",
            m_rigid, m_flex
        ),
    );
}
