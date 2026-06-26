// Port of tests/test_rawes_ic_aero_sign.py
//
// Regression tests for the RAWES IC aero sign: the rotor force should point
// roughly -body_z (downwind and upward) at this kite tethered-rotor attitude.
//
// RAWES clean IC attitude from windpower test_generate_ic.py:
//   pos      = [0, 90.415, -42.721]  NED
//   wind     = [0, 10, 0]            NED (+East)
//   body_z   = R_hub[:, 2] = [0, -0.9042, 0.4272]
//
// For a thrust-like force: F dot body_z < 0, F[1] > 0 (downwind), F[2] < 0 (up).

mod common;

use dynbem_rs::aero_io::{Mat3, RotorInputs, Vec3};
use dynbem_rs::aero_model::AeroModel;
use dynbem_rs::pitt_peters::PittPetersModel;
use dynbem_rs::polar::LinearPolar;

// Clean wind-aligned RAWES IC R_hub (columns are NED body axes).
// body_z = R_hub[:, 2] = [0.0, -0.9041543463617201, 0.4272059432582967]
fn r_rawes_ic() -> Mat3 {
    Mat3([
        [0.0, -1.0, 0.0],
        [0.42720594325829603, 0.0, -0.9041543463617201],
        [0.9041543463617201, 0.0, 0.4272059432582967],
    ])
}

fn rawes_ic_inputs(omega: f64, tilt_lon: f64, tilt_lat: f64) -> RotorInputs {
    RotorInputs {
        collective_rad: -0.18,
        tilt_lon,
        tilt_lat,
        R_hub: r_rawes_ic(),
        v_hub_world: Vec3::new(0.0, 0.0, 0.0),
        wind_world: Vec3::new(0.0, 10.0, 0.0),
        omega_rad_s: omega,
        rho_kg_m3: 1.225,
    }
}

/// Returns (F_dot_body_z, F_world_east, F_world_up).
fn force_projections(f_world: Vec3, r_hub: &Mat3) -> (f64, f64, f64) {
    let body_z = Vec3::new(r_hub.0[0][2], r_hub.0[1][2], r_hub.0[2][2]);
    (f_world.dot(body_z), f_world.0[1], -f_world.0[2])
}

// ---------------------------------------------------------------------------
// Quasi-static BEM: zero cyclic
// ---------------------------------------------------------------------------

#[test]
fn quasi_static_rawes_ic_force_points_against_body_z() {
    let model = common::qs_model();
    let inputs = rawes_ic_inputs(53.161687, 0.0, 0.0);
    let (result, _) = model.compute_forces(&inputs, &model.initial_state());
    let (f_dot_bz, downwind, up) = force_projections(result.F_world, &inputs.R_hub);
    assert!(
        f_dot_bz < 0.0,
        "QS RAWES IC: F dot body_z should be negative, got {f_dot_bz:.3}  \
         F_world={:?}",
        result.F_world.0,
    );
    assert!(
        downwind > 0.0,
        "QS RAWES IC: expected F downwind > 0, got {downwind:.3}"
    );
    assert!(up > 0.0, "QS RAWES IC: expected F up > 0, got {up:.3}");
}

// ---------------------------------------------------------------------------
// Pitt-Peters: reference model, zero cyclic
// ---------------------------------------------------------------------------

#[test]
fn pitt_peters_rawes_ic_force_has_expected_sign() {
    let defn = common::beaupoil_rotor();
    let polar = LinearPolar::from_properties(&defn.airfoil);
    let model = PittPetersModel::build(defn, 36, polar);
    let inputs = rawes_ic_inputs(38.132161, 0.0, 0.0);
    let (result, _) = model.compute_forces(&inputs, &model.initial_state());
    let (f_dot_bz, downwind, up) = force_projections(result.F_world, &inputs.R_hub);
    assert!(
        f_dot_bz < 0.0,
        "PP RAWES IC: F dot body_z should be negative, got {f_dot_bz:.3}  \
         F_world={:?}",
        result.F_world.0,
    );
    assert!(
        downwind > 0.0,
        "PP RAWES IC: expected F downwind > 0, got {downwind:.3}"
    );
    assert!(up > 0.0, "PP RAWES IC: expected F up > 0, got {up:.3}");
}

// ---------------------------------------------------------------------------
// Quasi-static BEM: captured trim cyclic
// ---------------------------------------------------------------------------

#[test]
fn quasi_static_rawes_ic_sign_with_trim_cyclic() {
    // tilt_lon=0.0, tilt_lat=0.022616 -- captured trim cyclic from windpower
    let model = common::qs_model();
    let inputs = rawes_ic_inputs(53.161687, 0.0, 0.022616);
    let (result, _) = model.compute_forces(&inputs, &model.initial_state());
    let (f_dot_bz, _, _) = force_projections(result.F_world, &inputs.R_hub);
    assert!(
        f_dot_bz < 0.0,
        "QS RAWES IC with trim cyclic: F dot body_z should be negative, got {f_dot_bz:.3}  \
         F_world={:?}",
        result.F_world.0,
    );
}
