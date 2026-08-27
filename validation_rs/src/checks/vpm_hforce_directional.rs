// Directional check: the VPM free-wake rotor accumulates an in-plane hub
// force ("H-force") from the blade-element loads, matching the BEM-family
// convention (see checks/h_force.rs for the quasi-static-BEM analogue and
// dynbem_rs/src/vpm/march.rs for the projection).
//
// Physical claim under test: with edgewise flow across a level disk the
// rotor's in-plane hub force points along the incoming freestream
// (downwind), the same directional convention as the BEM level-disk
// crosswind check. FlightCondition::v_hub is the freestream air velocity in
// the hub frame (relative wind = v_hub - blade velocity), so a +X freestream
// yields a +X in-plane force. In hover the loading is axisymmetric and the
// in-plane force vanishes.
//
// The VPM marches the flap DOF in the time domain, so the flapping-tilt part
// of the H-force emerges from the instantaneous disk geometry -- no separate
// harmonic flap solve is needed (unlike the BEM-family models). This check
// runs with a rigid rotor to isolate the always-present profile/induced-drag
// term; the flapping-tilt contribution is exercised in flap_directional.

use crate::helpers::*;
use crate::report::Report;
use std::f64::consts::PI;

pub fn check_vpm_hforce_directional(r: &mut Report) {
    r.begin_module(
        "vpm_hforce_directional",
        "Directional: VPM in-plane hub force (H-force) points downwind in edgewise flow; ~0 in hover",
    );

    let dt = (2.0 * PI / OMEGA) / STEPS_PER_REV as f64;
    let n_revs = 21 * STEPS_PER_REV;

    let defn = theory_rotor(10, 0.0);
    let rotor = make_fast_rotor(&defn);

    // --- Hover baseline: axisymmetric loading -> no in-plane force ---
    let (hover, _) = rotor.march(&hover_fc(8.0), None, dt, n_revs);
    let h_hover = hover.fx_hub.hypot(hover.fy_hub);
    r.info("hover", "fx_hub_N", hover.fx_hub, f64::NAN);
    r.info("hover", "fy_hub_N", hover.fy_hub, f64::NAN);
    r.info("hover", "h_force_N", h_hover, f64::NAN);
    r.info("hover", "thrust_N", hover.thrust, f64::NAN);
    r.assert_bool(
        "hover",
        "no_in_plane_force",
        h_hover,
        0.0,
        h_hover < 0.02 * hover.thrust.abs().max(1.0),
        &format!(
            "hover loading is axisymmetric -> in-plane force ~0, got |H|={:.4} N (T={:.2} N)",
            h_hover, hover.thrust
        ),
    );

    // --- Edgewise flow along +X: H-force points downwind (+X) ---
    let (fwd, _) = rotor.march(&forward_fc(8.0, 0.15), None, dt, n_revs);
    r.info("edgewise_flow", "fx_hub_N", fwd.fx_hub, f64::NAN);
    r.info("edgewise_flow", "fy_hub_N", fwd.fy_hub, f64::NAN);
    r.info("edgewise_flow", "thrust_N", fwd.thrust, f64::NAN);
    r.assert_bool(
        "edgewise_flow",
        "h_force_downwind",
        fwd.fx_hub,
        0.0,
        fwd.fx_hub > 0.05,
        &format!(
            "+X freestream should give a downwind (+X) H-force, got Fx_hub={:.4} N",
            fwd.fx_hub
        ),
    );
    r.assert_bool(
        "edgewise_flow",
        "mostly_aligned_with_flow",
        fwd.fx_hub.abs(),
        fwd.fy_hub.abs(),
        fwd.fx_hub.abs() > fwd.fy_hub.abs(),
        &format!(
            "H-force should be mostly along the flow axis: |Fx|={:.4} N vs |Fy|={:.4} N",
            fwd.fx_hub.abs(),
            fwd.fy_hub.abs()
        ),
    );
    r.assert_bool(
        "edgewise_flow",
        "larger_than_hover",
        fwd.fx_hub.hypot(fwd.fy_hub),
        h_hover,
        fwd.fx_hub.hypot(fwd.fy_hub) > h_hover,
        &format!(
            "edgewise flow should raise the H-force above the hover baseline: {:.4} N vs {:.4} N",
            fwd.fx_hub.hypot(fwd.fy_hub),
            h_hover
        ),
    );
}
