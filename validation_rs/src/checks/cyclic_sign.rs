// Directional checks for cyclic sign conventions and collective monotonicity.
//
// These replace the long-running unit tests that were formerly in
// dynbem_rs/src/vpm/mod.rs::tests. They exercise the same physics at
// the same tolerance level, but live here so the cargo test suite for
// dynbem_rs stays fast.
//
// All checks are purely directional (boolean pass/fail). No empirical data.

use crate::helpers::*;
use crate::report::Report;

pub fn check_cyclic_sign(r: &mut Report) {
    r.begin_module(
        "cyclic_sign",
        "Directional: collective monotone and cyclic sign conventions (AGENTS.md)",
    );

    let defn = theory_rotor(8, 2.0);
    let rotor = make_fast_rotor(&defn);
    let dt = 1.0 / (OMEGA * 2.0 / std::f64::consts::PI);
    // ~400 steps at STEPS_PER_REV
    let n_steps = 17 * STEPS_PER_REV;

    // --- collective monotone: 5 deg < 9 deg ---
    let fc_lo = hover_fc(5.0);
    let fc_hi = hover_fc(9.0);
    let (res_lo, _) = rotor.march(&fc_lo, None, dt, n_steps);
    let (res_hi, _) = rotor.march(&fc_hi, None, dt, n_steps);
    r.info("collective_5deg", "thrust_N", res_lo.thrust, f64::NAN);
    r.info("collective_9deg", "thrust_N", res_hi.thrust, f64::NAN);
    r.assert_bool(
        "collective_monotone",
        "hi_gt_lo",
        res_hi.thrust,
        res_lo.thrust,
        res_hi.thrust > res_lo.thrust,
        &format!(
            "thrust should rise with collective: {:.1} -> {:.1} N",
            res_lo.thrust, res_hi.thrust
        ),
    );

    // --- longitudinal cyclic: tilt_lon > 0 -> My < 0 (nose-down) ---
    let fc_base = hover_fc(8.0);
    let mut fc_lon = hover_fc(8.0);
    fc_lon.tilt_lon = 3.0_f64.to_radians();

    let (base, _) = rotor.march(&fc_base, None, dt, n_steps);
    let (tilted_lon, _) = rotor.march(&fc_lon, None, dt, n_steps);
    r.info("lon_cyclic_base", "my_hub_Nm", base.my_hub, f64::NAN);
    r.info(
        "lon_cyclic_tilted",
        "my_hub_Nm",
        tilted_lon.my_hub,
        f64::NAN,
    );
    r.assert_bool(
        "lon_cyclic",
        "my_lt_base_and_negative",
        tilted_lon.my_hub,
        base.my_hub,
        tilted_lon.my_hub < base.my_hub && tilted_lon.my_hub < 0.0,
        &format!(
            "tilt_lon>0 should give nose-down My<0: base {:.3} -> {:.3} N*m",
            base.my_hub, tilted_lon.my_hub
        ),
    );

    // --- lateral cyclic: tilt_lat > 0 -> Mx > 0 (roll-right) ---
    let mut fc_lat = hover_fc(8.0);
    fc_lat.tilt_lat = 3.0_f64.to_radians();

    let (tilted_lat, _) = rotor.march(&fc_lat, None, dt, n_steps);
    r.info("lat_cyclic_base", "mx_hub_Nm", base.mx_hub, f64::NAN);
    r.info(
        "lat_cyclic_tilted",
        "mx_hub_Nm",
        tilted_lat.mx_hub,
        f64::NAN,
    );
    r.assert_bool(
        "lat_cyclic",
        "mx_gt_base_and_positive",
        tilted_lat.mx_hub,
        base.mx_hub,
        tilted_lat.mx_hub > base.mx_hub && tilted_lat.mx_hub > 0.0,
        &format!(
            "tilt_lat>0 should give roll-right Mx>0: base {:.3} -> {:.3} N*m",
            base.mx_hub, tilted_lat.mx_hub
        ),
    );

    // --- crosswind: thrust bounded, hub moment nonzero, wake skews downstream ---
    let mut fc_cross = hover_fc(8.0);
    fc_cross.v_hub = [8.0, 0.0, 0.0]; // 8 m/s edgewise along +X
    let (cross, _) = rotor.march(&fc_cross, None, dt, n_steps);
    let moment = cross.mx_hub.hypot(cross.my_hub);
    r.info("crosswind", "thrust_N", cross.thrust, f64::NAN);
    r.info("crosswind", "hub_moment_Nm", moment, f64::NAN);
    r.info(
        "crosswind",
        "wake_centroid_x_m",
        cross.wake_centroid[0],
        f64::NAN,
    );
    r.assert_bool(
        "crosswind",
        "thrust_positive",
        cross.thrust,
        0.0,
        cross.thrust.is_finite() && cross.thrust > 0.0,
        &format!(
            "crosswind: thrust should be positive, got {:.1}",
            cross.thrust
        ),
    );
    r.assert_bool(
        "crosswind",
        "hub_moment_nonzero",
        moment,
        0.0,
        moment > 1e-3,
        &format!("crosswind should induce hub moment, got {:.4}", moment),
    );
    r.assert_bool(
        "crosswind",
        "wake_skews_downstream",
        cross.wake_centroid[0],
        0.0,
        cross.wake_centroid[0] > 0.05,
        &format!(
            "wake should skew downstream (+X), centroid_x = {:.3}",
            cross.wake_centroid[0]
        ),
    );
}
