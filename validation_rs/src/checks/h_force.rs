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

    // --- flapping-tilt H-force: flapback drag ---
    //
    // Classical rotor theory's second (usually larger) H-force term comes
    // from blade flapping: in forward flight the advancing/retreating lift
    // asymmetry drives a 1/rev flap response that, with the ~90 deg
    // aerodynamic-damping phase lag, tilts the disk AFT (flapback). The
    // thrust vector tilting rearward is an extra rearward in-plane force.
    //
    // Model it directly with a freely-hinged flap blade (omega_NR=0, so no
    // hub moment is transmitted -- this isolates the flapping-tilt H-force
    // from the transmitted-moment path) at a realistic Lock number, in
    // forward flight with real thrust (nonzero collective). Compare against
    // the identical rigid-blade rotor: enabling flap must make the H-force
    // MORE rearward, and the flapping contribution should dominate the
    // small profile-drag term.
    let i_beta = 0.03; // -> Lock number ~8 (see info below)
    let defn_flap = theory_rotor_flap(8, i_beta);
    let rotor_flap = QuasiStaticBEM::build(defn_flap, 36, theory_polar());
    let rotor_rigid = QuasiStaticBEM::build(theory_rotor(8, 0.0), 36, theory_polar());

    let collective = 8.0_f64.to_radians();
    let ff = |c: f64| RotorInputs {
        collective_rad: c,
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        R_hub: Mat3::eye(),
        v_hub_world: Vec3::new(8.0, 0.0, 0.0), // forward flight +X, still air
        wind_world: Vec3::zero(),
        rho_kg_m3: RHO,
        omega_rad_s: OMEGA,
    };

    let (res_flap, _) = rotor_flap.compute_forces(&ff(collective), &rotor_flap.initial_state());
    let (res_rigid, _) = rotor_rigid.compute_forces(&ff(collective), &rotor_rigid.initial_state());
    let f_flap = res_flap.F_world.0;
    let f_rigid = res_rigid.F_world.0;

    r.info("flapback", "lock_number", lock_number(i_beta), f64::NAN);
    r.info("flapback", "thrust_z_N", -f_flap[2], f64::NAN);
    r.info("flapback", "f_north_flap_N", f_flap[0], f64::NAN);
    r.info("flapback", "f_north_rigid_N", f_rigid[0], f64::NAN);

    // Flapback must be rearward (-X).
    r.assert_bool(
        "flapback",
        "flapping_h_force_is_rearward",
        f_flap[0],
        0.0,
        f_flap[0] < -0.05,
        &format!(
            "flapping in forward flight (+X) should tilt the disk aft -> rearward (-X) H-force, got F_north={:.3} N",
            f_flap[0]
        ),
    );
    // Flapping adds to (dominates) the rigid-blade profile-drag H-force.
    r.assert_bool(
        "flapback",
        "flapping_adds_rearward_drag",
        f_flap[0],
        f_rigid[0],
        f_flap[0] < f_rigid[0] - 0.05,
        &format!(
            "enabling flap should make the H-force more rearward: flap F_north={:.3} N vs rigid {:.3} N",
            f_flap[0], f_rigid[0]
        ),
    );
    // A freely-hinged blade transmits no hub moment: the flapping shows up
    // in the H-force, not as an airframe roll/pitch moment.
    let m_flap = res_flap.M_hub_world.0;
    r.info("flapback", "hub_moment_x_N_m", m_flap[0], f64::NAN);
    r.info("flapback", "hub_moment_y_N_m", m_flap[1], f64::NAN);
    r.assert_bool(
        "flapback",
        "free_hinge_transmits_no_moment",
        m_flap[0].hypot(m_flap[1]),
        0.0,
        m_flap[0].hypot(m_flap[1]) < 1e-6,
        &format!(
            "freely-hinged flap should transmit ~zero hub moment, got |M|={:.3e} N*m",
            m_flap[0].hypot(m_flap[1])
        ),
    );
}
