// Directional check: the rotor's net in-plane hub force ("H-force") under
// a pure crosswind, with the disk held level and collective at ~0 (no
// thrust). See AGENTS.md "Rotor rotation direction" for the hub-frame
// convention and dynbem_rs/src/bem_common.rs::assemble_result for how
// Fx_hub/Fy_hub are assembled into F_world.
//
// Physical claim under test: a horizontal rotor disk is not "invisible" to
// a crosswind just because collective is zero. Blade profile drag varies
// with azimuth (advancing/retreating asymmetry from the edgewise flow), and
// that asymmetric tangential loading sums to a net in-plane force that
// points along the wind -- the same mechanism that would make a freely
// flying rotor drift downwind. This is separate from (and in addition to)
// the hub *moment* (Mx_hub/My_hub) already checked in cyclic_sign.rs.

use crate::helpers::*;
use crate::report::Report;
use dynbem_rs::aero_io::{Mat3, RotorInputs, Vec3};
use dynbem_rs::aero_model::AeroModel;
use dynbem_rs::quasi_static_bem::QuasiStaticBEM;

pub fn check_h_force(r: &mut Report) {
    r.begin_module(
        "h_force",
        "Directional: level-disk crosswind produces an in-plane hub force along the wind",
    );

    let defn = theory_rotor(8, 0.0);
    let polar = theory_polar();
    let rotor = QuasiStaticBEM::build(defn, 36, polar);

    // Reasonable operating RPM (same as the rest of the theory suite), disk
    // perfectly level (R_hub = identity), collective ~0 so there is no
    // vertical thrust -- isolates the in-plane term.
    let omega = OMEGA;
    let wind_speed = 8.0; // m/s crosswind
    let wind = Vec3::new(0.0, wind_speed, 0.0); // pure +Y (East) crosswind

    let inputs = RotorInputs {
        collective_rad: 0.0,
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        R_hub: Mat3::eye(),
        v_hub_world: Vec3::zero(),
        wind_world: wind,
        rho_kg_m3: RHO,
        omega_rad_s: omega,
    };

    let (result, _) = rotor.compute_forces(&inputs, &rotor.initial_state());
    let f = result.F_world.0;

    r.info("level_disk_crosswind", "thrust_z_N", f[2], f64::NAN);
    r.info("level_disk_crosswind", "f_north_N", f[0], f64::NAN);
    r.info("level_disk_crosswind", "f_east_N", f[1], f64::NAN);

    // Collective ~0 -> no vertical thrust.
    r.assert_bool(
        "level_disk_crosswind",
        "no_vertical_thrust",
        f[2],
        0.0,
        f[2].abs() < 1.0,
        &format!(
            "collective=0 should give ~zero vertical force, got F_z={:.3} N",
            f[2]
        ),
    );

    // The wind is purely +Y (East): the horizontal force should point the
    // same way (downwind), not across or against it.
    r.assert_bool(
        "level_disk_crosswind",
        "force_along_wind",
        f[1],
        0.0,
        f[1] > 0.05,
        &format!(
            "crosswind +Y should push the disk +Y (downwind), got F_east={:.3} N",
            f[1]
        ),
    );
    r.assert_bool(
        "level_disk_crosswind",
        "no_cross_axis_force",
        f[0],
        0.0,
        f[0].abs() < 0.1 * f[1].abs().max(1.0),
        &format!(
            "pure +Y crosswind should not push North/South, got F_north={:.3} N (F_east={:.3} N)",
            f[0], f[1]
        ),
    );

    // Reversing the wind should reverse the force.
    let inputs_rev = RotorInputs {
        wind_world: Vec3::new(0.0, -wind_speed, 0.0),
        ..inputs.clone()
    };
    let (result_rev, _) = rotor.compute_forces(&inputs_rev, &rotor.initial_state());
    let f_rev = result_rev.F_world.0;
    r.info("reversed_crosswind", "f_east_N", f_rev[1], f64::NAN);
    r.assert_bool(
        "reversed_crosswind",
        "force_flips_with_wind",
        f_rev[1],
        0.0,
        f_rev[1] < -0.05,
        &format!(
            "reversing wind to -Y should reverse the force to F_east<0, got {:.3} N",
            f_rev[1]
        ),
    );

    // No wind at all (hover, zero collective) -> no in-plane force either.
    let inputs_calm = RotorInputs {
        wind_world: Vec3::zero(),
        ..inputs.clone()
    };
    let (result_calm, _) = rotor.compute_forces(&inputs_calm, &rotor.initial_state());
    let f_calm = result_calm.F_world.0;
    r.info("no_wind", "f_east_N", f_calm[1], f64::NAN);
    r.assert_bool(
        "no_wind",
        "no_force_without_wind",
        f_calm[1],
        0.0,
        f_calm[1].abs() < 0.5,
        &format!(
            "zero wind + zero collective should give ~zero in-plane force, got F_east={:.3} N",
            f_calm[1]
        ),
    );

    // --- disk tilting toward the wind should shrink the H-force to zero ---
    //
    // Roll the disk about the North (X) axis by phi: hub_axis rotates from
    // [0,0,1] (level) toward [0,1,0] (pointing straight into the +Y wind).
    // As phi -> 90 deg the relative wind becomes purely axial (v_edge -> 0),
    // and the in-plane hub force -- which only exists because of the
    // azimuthal (edgewise-flow) asymmetry -- should vanish with it.
    fn hub_force_mag(rotor: &QuasiStaticBEM<dynbem_rs::polar::LinearPolar>, phi: f64) -> f64 {
        let (s, c) = phi.sin_cos();
        let r_hub = Mat3([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]]);
        let inputs = RotorInputs {
            collective_rad: 0.0,
            tilt_lon: 0.0,
            tilt_lat: 0.0,
            R_hub: r_hub,
            v_hub_world: Vec3::zero(),
            wind_world: Vec3::new(0.0, 8.0, 0.0),
            rho_kg_m3: RHO,
            omega_rad_s: OMEGA,
        };
        let (result, _) = rotor.compute_forces(&inputs, &rotor.initial_state());
        // Undo the rotation to read the in-plane hub-frame force back out.
        let f_hub = r_hub.transpose() * result.F_world;
        f_hub.0[0].hypot(f_hub.0[1])
    }

    let mag_level = hub_force_mag(&rotor, 0.0);
    let mag_mid = hub_force_mag(&rotor, 45.0_f64.to_radians());
    let mag_aligned = hub_force_mag(&rotor, 89.0_f64.to_radians());

    r.info("tilt_toward_wind", "h_force_level_N", mag_level, f64::NAN);
    r.info("tilt_toward_wind", "h_force_45deg_N", mag_mid, f64::NAN);
    r.info("tilt_toward_wind", "h_force_89deg_N", mag_aligned, f64::NAN);
    r.assert_bool(
        "tilt_toward_wind",
        "shrinks_monotonically",
        mag_mid,
        mag_level,
        mag_mid < mag_level && mag_aligned < mag_mid,
        &format!(
            "H-force should shrink as the disk tilts into the wind: level={:.4} N, 45deg={:.4} N, 89deg={:.4} N",
            mag_level, mag_mid, mag_aligned
        ),
    );
    r.assert_bool(
        "tilt_toward_wind",
        "vanishes_when_aligned",
        mag_aligned,
        0.0,
        mag_aligned < 0.02,
        &format!(
            "H-force should be ~0 once the disk axis aligns with the wind, got {:.4} N",
            mag_aligned
        ),
    );

    // --- textbook H-force: drag opposing the direction of flight ---
    //
    // Classical rotor theory defines H-force as a drag-like force that
    // opposes the aircraft's forward flight, not "whatever direction the
    // wind happens to blow a stationary rotor". Model that directly: the
    // hub itself moves at +X through still air (wind_world = 0), so the
    // relative wind seen by the disk is v_rel = wind - v_hub = -X (apparent
    // headwind from the nose). The in-plane force should point rearward
    // (-X, opposing the direction of flight), matching the standard
    // definition of rotor profile-drag H-force.
    let inputs_forward_flight = RotorInputs {
        collective_rad: 0.0,
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        R_hub: Mat3::eye(),
        v_hub_world: Vec3::new(8.0, 0.0, 0.0),
        wind_world: Vec3::zero(),
        rho_kg_m3: RHO,
        omega_rad_s: OMEGA,
    };
    let (result_ff, _) = rotor.compute_forces(&inputs_forward_flight, &rotor.initial_state());
    let f_ff = result_ff.F_world.0;
    r.info("forward_flight_still_air", "f_north_N", f_ff[0], f64::NAN);
    r.assert_bool(
        "forward_flight_still_air",
        "drag_opposes_flight_direction",
        f_ff[0],
        0.0,
        f_ff[0] < -0.05,
        &format!(
            "flying forward (+X) through still air should give a rearward (-X) H-force, got F_north={:.3} N",
            f_ff[0]
        ),
    );
}
