// Port of tests/test_trim.py core trim-solver behavior into Rust integration tests.
//
// Python now keeps lightweight API smoke checks; detailed trim behavior and
// solver convergence lives here for faster local iteration.

mod common;

use dynbem_rs::aero_io::{Mat3, RotorInputs, Vec3};
use dynbem_rs::aero_model::{AeroModel, RotorStateExt};
use dynbem_rs::oye::OyeBEMModel;
use dynbem_rs::pitt_peters::PittPetersModel;
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::trim::{relax_inflow, solve_trim_cyclic};

const OMEGA: f64 = 28.0;
const COLLECTIVE: f64 = -9.0_f64.to_radians();
const TOL_OYE: f64 = 0.05;
const TOL_PITT: f64 = 1.6;

fn level_r() -> Mat3 {
    Mat3::eye()
}

fn base_inputs(wind_world: Vec3) -> RotorInputs {
    RotorInputs {
        collective_rad: COLLECTIVE,
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        R_hub: level_r(),
        v_hub_world: Vec3::zero(),
        wind_world,
        rho_kg_m3: 1.225,
        omega_rad_s: OMEGA,
    }
}

fn moments_hub<M: AeroModel>(
    aero: &M,
    state: &M::State,
    collective: f64,
    tilt_lon: f64,
    tilt_lat: f64,
    wind_world: Vec3,
) -> (f64, f64) {
    let mut inputs = base_inputs(wind_world);
    inputs.collective_rad = collective;
    inputs.tilt_lon = tilt_lon;
    inputs.tilt_lat = tilt_lat;
    let (res, _) = aero.compute_forces(&inputs, state);
    let m_hub = inputs.R_hub.transpose() * res.M_hub_world;
    (m_hub[0], m_hub[1])
}

fn run_trim<M: AeroModel>(
    aero: &M,
    state: M::State,
    wind_world: Vec3,
    target_x: f64,
    target_y: f64,
    tol: f64,
) -> dynbem_rs::trim::TrimOutcome<M::State> {
    let inputs = base_inputs(wind_world);
    solve_trim_cyclic(
        aero,
        state,
        &inputs,
        target_x,
        target_y,
        0.0,
        0.0,
        -0.261_799_387_799_149_4,
        0.261_799_387_799_149_4,
        tol,
        50,
        0.008_726_646_259_971_648,
        0.005,
        100,
        0,
    )
}

#[test]
fn oye_hover_trim_is_near_zero() {
    let defn = common::beaupoil_rotor();
    let polar = LinearPolar::from_properties(&defn.airfoil);
    let aero = OyeBEMModel::build(defn, 36, polar);
    let out = run_trim(&aero, aero.initial_state(), Vec3::zero(), 0.0, 0.0, TOL_OYE);
    assert!(
        out.converged,
        "Oye hover trim did not converge: iters={} mx={:.4} my={:.4}",
        out.iterations,
        out.mx_residual,
        out.my_residual
    );
    assert!(out.tilt_lon.abs() < 0.5_f64.to_radians());
    assert!(out.tilt_lat.abs() < 0.5_f64.to_radians());
}

#[test]
fn pitt_forward_trim_residual_below_tolerance() {
    let defn = common::beaupoil_rotor();
    let polar = LinearPolar::from_properties(&defn.airfoil);
    let aero = PittPetersModel::build(defn, 36, polar);
    let out = run_trim(
        &aero,
        aero.initial_state(),
        Vec3::new(0.0, 10.0, 0.0),
        0.0,
        0.0,
        TOL_PITT,
    );
    assert!(
        out.converged,
        "Pitt forward trim did not converge: iters={} mx={:.4} my={:.4}",
        out.iterations,
        out.mx_residual,
        out.my_residual
    );
    assert!(out.mx_residual.abs() < TOL_PITT);
    assert!(out.my_residual.abs() < TOL_PITT);
}

#[test]
fn trim_residual_matches_direct_evaluation_pitt() {
    let defn = common::beaupoil_rotor();
    let polar = LinearPolar::from_properties(&defn.airfoil);
    let aero = PittPetersModel::build(defn, 36, polar);
    let out = run_trim(
        &aero,
        aero.initial_state(),
        Vec3::new(0.0, 10.0, 0.0),
        0.0,
        0.0,
        TOL_PITT,
    );

    let (mx, my) = moments_hub(
        &aero,
        &out.final_state,
        COLLECTIVE,
        out.tilt_lon,
        out.tilt_lat,
        Vec3::new(0.0, 10.0, 0.0),
    );

    assert!((mx - out.mx_residual).abs() < 1e-6);
    assert!((my - out.my_residual).abs() < 1e-6);
}

#[test]
fn trim_to_nonzero_target_moment_pitt() {
    let defn = common::beaupoil_rotor();
    let polar = LinearPolar::from_properties(&defn.airfoil);
    let aero = PittPetersModel::build(defn, 36, polar);
    let m_target = 5.0;

    let out = run_trim(
        &aero,
        aero.initial_state(),
        Vec3::new(0.0, 10.0, 0.0),
        0.0,
        m_target,
        TOL_PITT,
    );

    assert!(
        out.converged,
        "Pitt target-moment trim did not converge: iters={} mx={:.4} my={:.4}",
        out.iterations,
        out.mx_residual,
        out.my_residual
    );

    let (mx, my) = moments_hub(
        &aero,
        &out.final_state,
        COLLECTIVE,
        out.tilt_lon,
        out.tilt_lat,
        Vec3::new(0.0, 10.0, 0.0),
    );
    assert!(mx.abs() < TOL_PITT, "Mx={mx:.4} should be near 0");
    assert!(
        (my - m_target).abs() < TOL_PITT,
        "My={my:.4} should be near {m_target:.4}"
    );
}

#[test]
fn relax_inflow_settles_to_fixed_point_for_both_models() {
    let defn = common::beaupoil_rotor();
    let polar_pp = LinearPolar::from_properties(&defn.airfoil);
    let polar_oye = LinearPolar::from_properties(&defn.airfoil);
    let pp = PittPetersModel::build(defn.clone(), 36, polar_pp);
    let oye = OyeBEMModel::build(defn, 36, polar_oye);
    let inputs = base_inputs(Vec3::new(0.0, 10.0, 0.0));

    let s1_pp = relax_inflow(&pp, pp.initial_state(), &inputs, 500, 0.005);
    let s2_pp = relax_inflow(&pp, s1_pp.clone(), &inputs, 500, 0.005);
    let d_pp: f64 = s1_pp
        .get_inflow()
        .iter()
        .zip(s2_pp.get_inflow().iter())
        .map(|(a, b)| (b - a) * (b - a))
        .sum::<f64>()
        .sqrt();
    assert!(d_pp < 1e-4, "Pitt inflow not settled: delta={d_pp:.4e}");

    let s1_oye = relax_inflow(&oye, oye.initial_state(), &inputs, 500, 0.005);
    let s2_oye = relax_inflow(&oye, s1_oye.clone(), &inputs, 500, 0.005);
    let d_oye: f64 = s1_oye
        .get_inflow()
        .iter()
        .zip(s2_oye.get_inflow().iter())
        .map(|(a, b)| (b - a) * (b - a))
        .sum::<f64>()
        .sqrt();
    assert!(d_oye < 1e-4, "Oye inflow not settled: delta={d_oye:.4e}");
}

#[test]
fn solver_reduces_baseline_disturbance_pitt() {
    let defn = common::beaupoil_rotor();
    let polar = LinearPolar::from_properties(&defn.airfoil);
    let aero = PittPetersModel::build(defn, 36, polar);
    let wind = Vec3::new(0.0, 10.0, 0.0);

    let mut state = aero.initial_state();
    let inputs = base_inputs(wind);
    for _ in 0..200 {
        let (_, deriv) = aero.compute_forces(&inputs, &state);
        let arr = state
            .get_inflow()
            .iter()
            .zip(deriv.get_inflow().iter())
            .map(|(x, dx)| x + 0.005 * dx)
            .collect::<Vec<f64>>();
        state.set_inflow(arr);
    }

    let (mx0, my0) = moments_hub(&aero, &state, COLLECTIVE, 0.0, 0.0, wind);
    let baseline = (mx0 * mx0 + my0 * my0).sqrt();
    assert!(baseline > 10.0, "baseline too small: {baseline:.2}");

    let out = run_trim(&aero, state, wind, 0.0, 0.0, TOL_PITT);
    let trim_mag = (out.mx_residual * out.mx_residual + out.my_residual * out.my_residual).sqrt();

    // Keep this robust to moderate solver/mode variation while still enforcing
    // strong disturbance cancellation.
    assert!(
        trim_mag < baseline / 70.0,
        "solver did not sufficiently cancel disturbance: baseline={baseline:.2} trim={trim_mag:.4}"
    );
}

#[test]
fn trim_clips_to_bounds_oye() {
    let defn = common::beaupoil_rotor();
    let polar = LinearPolar::from_properties(&defn.airfoil);
    let aero = OyeBEMModel::build(defn, 36, polar);
    let tight = 1.0_f64.to_radians();
    let out = solve_trim_cyclic(
        &aero,
        aero.initial_state(),
        &base_inputs(Vec3::new(0.0, 10.0, 0.0)),
        0.0,
        0.0,
        0.0,
        0.0,
        -tight,
        tight,
        0.01,
        20,
        0.008_726_646_259_971_648,
        0.005,
        100,
        0,
    );
    assert!(out.tilt_lon >= -tight && out.tilt_lon <= tight);
    assert!(out.tilt_lat >= -tight && out.tilt_lat <= tight);
}
