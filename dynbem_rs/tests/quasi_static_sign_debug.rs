// Port of tests/test_quasi_static_sign_debug.py
//
// Four cases with identity R_hub (hub_axis = NED +Z = downward) that isolate
// the windmill psi-loop sign bug.  With R_hub = I:
//   v_climb         = wind[2]   (positive = helicopter, negative = windmill)
//   v_inplane_hub   = wind[0:2] (in-plane wind triggers psi-loop when large)
//
// Case A: hover              (mu=0, v_climb=0)  -- axial path
// Case B: in-plane only      (mu>0, v_climb=0)  -- psi-loop
// Case C: in-plane + axial   (mu>0, v_climb<0)  -- psi-loop, RAWES IC analog
// Case D: axial windmill     (mu=0, v_climb<0)  -- axial windmill path

mod common;

use dynbem_rs::aero_io::{Mat3, RotorInputs, Vec3};
use dynbem_rs::aero_model::AeroModel;

fn identity_inputs(wind: [f64; 3], omega: f64) -> RotorInputs {
    RotorInputs {
        collective_rad: -0.18,
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        R_hub: Mat3([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        v_hub_world: Vec3::new(0.0, 0.0, 0.0),
        wind_world: Vec3::new(wind[0], wind[1], wind[2]),
        omega_rad_s: omega,
        t: 0.0,
        rho_kg_m3: 1.225,
    }
}

// QS trim omega from the original RAWES IC scenario
const OMEGA_QS: f64 = 53.161687;
// PP trim omega from the same scenario
const OMEGA_PP: f64 = 38.132161;
// Decomposed RAWES IC in-plane and axial components (hub_axis = NED +Z)
const V_CLIMB: f64 = -9.04;
const V_INPLANE: f64 = 4.27;

// ---------------------------------------------------------------------------
// Case A: hover -- force magnitude must be non-trivially nonzero
// ---------------------------------------------------------------------------

#[test]
fn case_a_hover_thrust_nontrivial_at_qs_omega() {
    let model = common::qs_model();
    let (result, _) = model.compute_forces(&identity_inputs([0.0, 0.0, 0.0], OMEGA_QS), &model.initial_state());
    let f2 = result.F_world.0[2];
    assert!(f2.abs() > 10.0, "Case A hover OMEGA_QS: near-zero thrust: {f2:.3}");
}

#[test]
fn case_a_hover_thrust_nontrivial_at_pp_omega() {
    let model = common::qs_model();
    let (result, _) = model.compute_forces(&identity_inputs([0.0, 0.0, 0.0], OMEGA_PP), &model.initial_state());
    let f2 = result.F_world.0[2];
    assert!(f2.abs() > 10.0, "Case A hover OMEGA_PP: near-zero thrust: {f2:.3}");
}

// ---------------------------------------------------------------------------
// Case B: in-plane wind only -- sign must agree at both omegas
// ---------------------------------------------------------------------------

#[test]
fn case_b_inplane_sign_consistent() {
    let model = common::qs_model();
    let (r_high, _) = model.compute_forces(&identity_inputs([V_INPLANE, 0.0, 0.0], OMEGA_QS), &model.initial_state());
    let (r_low, _)  = model.compute_forces(&identity_inputs([V_INPLANE, 0.0, 0.0], OMEGA_PP), &model.initial_state());
    let f2_high = r_high.F_world.0[2];
    let f2_low  = r_low.F_world.0[2];
    assert_eq!(
        f2_high > 0.0, f2_low > 0.0,
        "Case B: sign flips between omegas -- OMEGA_QS={f2_high:.1} OMEGA_PP={f2_low:.1}",
    );
}

// ---------------------------------------------------------------------------
// Case C: psi-loop + axial (RAWES IC analog) -- must match Case D sign
// ---------------------------------------------------------------------------

#[test]
fn case_c_psiloop_matches_axial_windmill_sign_at_high_omega() {
    let model = common::qs_model();
    let (r_c, _) = model.compute_forces(&identity_inputs([V_INPLANE, 0.0, V_CLIMB], OMEGA_QS), &model.initial_state());
    let (r_d, _) = model.compute_forces(&identity_inputs([0.0, 0.0, V_CLIMB], OMEGA_QS), &model.initial_state());
    let f2_c = r_c.F_world.0[2];
    let f2_d = r_d.F_world.0[2];
    assert_eq!(
        f2_c > 0.0, f2_d > 0.0,
        "Case C vs D at OMEGA_QS: sign mismatch -- C={f2_c:.1} D={f2_d:.1}",
    );
}

#[test]
fn case_c_sign_consistent_across_omegas() {
    let model = common::qs_model();
    let (r_high, _) = model.compute_forces(&identity_inputs([V_INPLANE, 0.0, V_CLIMB], OMEGA_QS), &model.initial_state());
    let (r_low, _)  = model.compute_forces(&identity_inputs([V_INPLANE, 0.0, V_CLIMB], OMEGA_PP), &model.initial_state());
    let f2_high = r_high.F_world.0[2];
    let f2_low  = r_low.F_world.0[2];
    assert_eq!(
        f2_high > 0.0, f2_low > 0.0,
        "Case C: sign flips between omegas (windmill psi-loop regression) -- \
         OMEGA_QS={f2_high:.1} OMEGA_PP={f2_low:.1}",
    );
}

// ---------------------------------------------------------------------------
// Case D: axial windmill -- F_world[2] < 0 (upward force) at both omegas
// ---------------------------------------------------------------------------

#[test]
fn case_d_axial_windmill_force_is_upward_at_qs_omega() {
    let model = common::qs_model();
    let (result, _) = model.compute_forces(&identity_inputs([0.0, 0.0, V_CLIMB], OMEGA_QS), &model.initial_state());
    let f2 = result.F_world.0[2];
    assert!(f2 < 0.0, "Case D OMEGA_QS: expected F_world[2] < 0, got {f2:.3}");
}

#[test]
fn case_d_axial_windmill_force_is_upward_at_pp_omega() {
    let model = common::qs_model();
    let (result, _) = model.compute_forces(&identity_inputs([0.0, 0.0, V_CLIMB], OMEGA_PP), &model.initial_state());
    let f2 = result.F_world.0[2];
    assert!(f2 < 0.0, "Case D OMEGA_PP: expected F_world[2] < 0, got {f2:.3}");
}
