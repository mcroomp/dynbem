// VPM empirical test: forward-flight autorotation vs Wheatley & Hood
// (NACA TR 515, 1935), Tables III and IV.
//
// Context: this is the RAWES operating regime -- a rotor autorotating in
// forward/crosswind flight, extracting power from the wind.  The PCA-2
// autogiro was held at fixed (mu, alpha, Omega, pitch) in the NACA full-scale
// wind tunnel; at equilibrium Q_aero = 0 by construction.  We compare the VPM
// predicted CL against the measured value at each operating point.
//
// CL (airplane axes) = lift / (q * pi * R^2) where q = 0.5 * rho * V^2.
// The VPM returns axial thrust T (along hub axis -Z_hub); with shaft angle
// alpha the lift contribution is T * cos(alpha).  In-plane H-force is not
// modelled, so only CL (not CD) is asserted.
//
// To refresh the tolerance CSV after a model change (release mode, ~60 s):
//   $env:REWRITE_VPM_EMPIRICAL=1
//   cargo test --test vpm_forward_flight_empirical --release -- --nocapture

mod common;

use dynbem_rs::cyclic::ControlGains;
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::rotor_definition::{
    BladeGeometry, LinearPolarParameters, PitchActuation, RotorDefinition,
};
use dynbem_rs::vpm_rotor::{FlightCondition, VpmRotor, VpmRotorConfig};
use serde::{Deserialize, Serialize};
use std::env;
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// PCA-2 geometry constants (from rotors/wheatley_pca2/rotor.yaml)
// ---------------------------------------------------------------------------
const R: f64 = 6.85;   // tip radius (m)
const RHO: f64 = 1.225;

// Steps per revolution and settle / averaging lengths.
// 48 steps/rev = 7.5 deg/step.  4 settle + 2 avg = 6 revolutions per point.
const STEPS_PER_REV: usize = 48;
const N_SETTLE_REV: usize = 4;
const N_AVG_REV: usize = 2;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn pca2_rotor() -> VpmRotor<LinearPolar> {
    let defn = RotorDefinition {
        blade: BladeGeometry {
            n_blades: 4,
            radius_m: R,
            root_cutout_m: 1.20,
            chord_m: 0.55,
            twist_deg: 0.0,
            n_elements: 12,  // 12 stations -- fast; use 30 for production sweeps
            tip_loss: true,
            r_stations_m: Vec::new(),
            chord_stations_m: Vec::new(),
            twist_stations_deg: Vec::new(),
        },
        airfoil: LinearPolarParameters {
            CL0: 0.0,
            CL_alpha_per_rad: 5.73,   // 2*pi for symmetric thin airfoil
            CD0: 0.0098,
            alpha_stall_deg: 15.0,
        },
        control: None,
        pitch_actuation: PitchActuation::DirectMechanical,
        flap: None,
        name: "pca2_vpm_fwd".to_string(),
        description: String::new(),
    };
    let polar = LinearPolar::from_properties(&defn.airfoil);
    let config = VpmRotorConfig {
        max_particles: 2000,
        sigma: 0.50,          // ~chord scale for the PCA-2 (chord = 0.55 m)
        relax: 0.35,
        nonlinear_lifting_line: true,
        tip_clustering: true,
        local_core: true,
        barnes_hut: true,
        bh_theta: 0.5,
        bh_min_particles: 200,  // engage BH early -- wake fills fast at 4 blades
    };
    VpmRotor::new(&defn, polar, ControlGains::default(), config)
}

/// Build a FlightCondition for a Wheatley/Hood operating point.
///
/// The hub is tilted back by `alpha_deg` from vertical so that the rotor disk
/// faces into the oncoming wind.  The wind is along +X in NED.  Rotating the
/// world wind vector into hub frame gives:
///   v_hub = [V * cos(a), 0, -V * sin(a)]
/// where `a` = alpha_deg in radians and `V = omega * R * mu / cos(a)`.
fn wheatley_fc(pitch_deg: f64, mu: f64, alpha_deg: f64, n_rpm: f64) -> FlightCondition {
    let omega = n_rpm * PI / 30.0;
    let a = alpha_deg.to_radians();
    let v = omega * R * mu / a.cos();
    FlightCondition {
        collective_rad: pitch_deg.to_radians(),
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        v_hub: [v * a.cos(), 0.0, -v * a.sin()],
        omega_rad_s: omega,
        rho: RHO,
    }
}

/// CL (airplane axes) from VPM thrust.
/// VPM returns T along -Z_hub; projected onto the lift direction gives T*cos(a).
fn cl_model(thrust: f64, mu: f64, alpha_deg: f64, n_rpm: f64) -> f64 {
    let omega = n_rpm * PI / 30.0;
    let a = alpha_deg.to_radians();
    let v = omega * R * mu / a.cos();
    let q = 0.5 * RHO * v * v;
    let area = PI * R * R;
    thrust * a.cos() / (q * area)
}

/// CQ = Q / (rho * A * (omega*R)^2 * R) -- tip-speed normalized.
fn cq_model(torque: f64, n_rpm: f64) -> f64 {
    let omega = n_rpm * PI / 30.0;
    let area = PI * R * R;
    torque / (RHO * area * (omega * R).powi(2) * R)
}

// ---------------------------------------------------------------------------
// CSV record
// ---------------------------------------------------------------------------
#[derive(Serialize, Deserialize, Clone)]
struct VpmFwdRecord {
    label: String,
    table: String,
    pitch_deg: f64,
    mu: f64,
    alpha_deg: f64,
    n_rpm: f64,
    cl_meas: f64,
    cl_max_err: f64,
}

// ---------------------------------------------------------------------------
// Test
// ---------------------------------------------------------------------------
#[test]
fn vpm_wheatley_forward_flight_cl() {
    let csv_data = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/vpm_forward_flight_empirical.csv"
    ));

    let rewrite_mode = env::var("REWRITE_VPM_EMPIRICAL").is_ok();

    let mut rdr = common::csv_reader_with_comments(csv_data.as_bytes());
    let mut records: Vec<VpmFwdRecord> = Vec::new();
    let mut failures: Vec<String> = Vec::new();

    let rotor = pca2_rotor();

    for result in rdr.deserialize() {
        let record: VpmFwdRecord = result.expect("failed to deserialize CSV row");

        let omega = record.n_rpm * PI / 30.0;
        let t_rev = 2.0 * PI / omega;
        let dt = t_rev / STEPS_PER_REV as f64;
        let n_total = (N_SETTLE_REV + N_AVG_REV) * STEPS_PER_REV;

        let fc = wheatley_fc(record.pitch_deg, record.mu, record.alpha_deg, record.n_rpm);

        let t0 = std::time::Instant::now();
        let res = rotor.simulate(&fc, dt, n_total);
        let elapsed = t0.elapsed().as_secs_f64();

        let cl = cl_model(res.thrust, record.mu, record.alpha_deg, record.n_rpm);
        let cq = cq_model(res.torque, record.n_rpm);
        let err = (cl - record.cl_meas).abs() / record.cl_meas;

        let max_err = record.cl_max_err;

        eprintln!(
            "{} (table {}, pitch={:.1} deg, mu={:.3}, alpha={:.1} deg, N={:.0} rpm): \
             CL_vpm={:.4}  CL_meas={:.4}  err={:.1}%  \
             CQ={:+.5}  T={:.1} N  {:.1} s",
            record.label, record.table, record.pitch_deg,
            record.mu, record.alpha_deg, record.n_rpm,
            cl, record.cl_meas, err * 100.0,
            cq, res.thrust, elapsed,
        );

        let mut r = record.clone();
        if rewrite_mode {
            r.cl_max_err = (err + 0.05).max(err * 3.0);  // observed error + 5 pp margin
        } else {
            if err > max_err {
                failures.push(format!(
                    "{}: CL err {:.1}% exceeds ceiling {:.1}%  \
                     (CL_vpm={:.4} vs CL_meas={:.4})",
                    record.label, err * 100.0, max_err * 100.0, cl, record.cl_meas
                ));
            }
        }
        records.push(r);
    }

    if rewrite_mode {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/vpm_forward_flight_empirical.csv"
        );
        let mut wtr = csv::Writer::from_path(path).unwrap();
        for r in &records {
            wtr.serialize(r).unwrap();
        }
        wtr.flush().unwrap();
        eprintln!("Rewrote vpm_forward_flight_empirical.csv with updated error ceilings.");
        return;
    }

    if !failures.is_empty() {
        panic!(
            "VPM forward-flight CL validation failures:\n{}",
            failures.join("\n")
        );
    }
}
