// Shared fixtures for integration tests.
#![allow(dead_code)]

use dynbem_rs::polar::LinearPolar;
use dynbem_rs::aero_io::{Mat3, RotorInputs, Vec3};
use dynbem_rs::quasi_static_bem::QuasiStaticBEM;
use dynbem_rs::rotor_definition::{
    BladeGeometry, ControlProperties, LinearPolarParameters, PitchActuation, RotorDefinition,
};
use csv;

/// Beaupoil 2026 / RAWES rotor definition (matches rotors/beaupoil_2026/rotor.yaml).
pub fn beaupoil_rotor() -> RotorDefinition {
    RotorDefinition {
        blade: BladeGeometry {
            n_blades: 4,
            radius_m: 2.5,
            root_cutout_m: 0.5,
            chord_m: 0.20,
            twist_deg: 0.0,
            n_elements: 10,
            tip_loss: true,
            r_stations_m: Vec::new(),
            chord_stations_m: Vec::new(),
            twist_stations_deg: Vec::new(),
        },
        airfoil: LinearPolarParameters {
            CL0: 0.393,
            CL_alpha_per_rad: 5.79,
            CD0: 0.0079,
            alpha_stall_deg: 13.0,
        },
        control: Some(ControlProperties {
            swashplate_pitch_gain_rad: 0.3,
            swashplate_phase_deg: Some(0.0),
        }),
        pitch_actuation: PitchActuation::DirectMechanical,
        flap: None,
        name: "beaupoil_2026".to_string(),
        description: String::new(),
    }
}

pub fn qs_model() -> QuasiStaticBEM<LinearPolar> {
    let defn = beaupoil_rotor();
    let polar = LinearPolar::from_properties(&defn.airfoil);
    QuasiStaticBEM::build(defn, 36, polar)
}

/// Castles and Gray (1951) NACA TN-2474 rotor definition.
pub fn castles_gray_rotor() -> RotorDefinition {
    RotorDefinition {
        blade: BladeGeometry::uniform(3, 0.914, 0.155, 0.0479, 0.0, 30),
        airfoil: LinearPolarParameters {
            CL0: 0.0,
            CL_alpha_per_rad: 5.90,
            CD0: 0.01046,
            alpha_stall_deg: 15.5,
        },
        control: None,
        pitch_actuation: PitchActuation::DirectMechanical,
        flap: None,
        name: "castles_gray_6ft".to_string(),
        description: String::new(),
    }
}

pub fn hover_inputs(theta_deg: f64, rpm: f64) -> RotorInputs {
    let omega = rpm * std::f64::consts::PI / 30.0;
    RotorInputs {
        collective_rad: theta_deg.to_radians(),
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        R_hub: Mat3::eye(),
        v_hub_world: Vec3::zero(),
        wind_world: Vec3::zero(),
        rho_kg_m3: 1.225,
        omega_rad_s: omega,
    }
}

pub fn ct_from_result(defn: &RotorDefinition, rpm: f64, thrust_n: f64) -> f64 {
    let omega = rpm * std::f64::consts::PI / 30.0;
    let r = defn.blade.radius_m;
    let a = std::f64::consts::PI * r * r;
    thrust_n / (1.225 * a * (omega * r).powi(2))
}

/// CQ = Q_spin / (rho * A * (omega*R)^2 * R).
/// Q_spin > 0 opposes rotation (helicopter mode); Q_spin < 0 drives it (autorotation).
pub fn cq_from_result(defn: &RotorDefinition, rpm: f64, q_spin: f64) -> f64 {
    let omega = rpm * std::f64::consts::PI / 30.0;
    let r = defn.blade.radius_m;
    let a = std::f64::consts::PI * r * r;
    q_spin / (1.225 * a * (omega * r).powi(2) * r)
}

/// RotorInputs for axial descent.  v_descent_m_s > 0 = aircraft moving downward
/// in NED, which gives v_climb = -v_descent_m_s (upward relative air through disk).
pub fn descent_inputs(theta_deg: f64, rpm: f64, v_descent_m_s: f64) -> RotorInputs {
    let omega = rpm * std::f64::consts::PI / 30.0;
    RotorInputs {
        collective_rad: theta_deg.to_radians(),
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        R_hub: Mat3::eye(),
        v_hub_world: Vec3::new(0.0, 0.0, v_descent_m_s),
        wind_world: Vec3::zero(),
        rho_kg_m3: 1.225,
        omega_rad_s: omega,
    }
}

/// Compute error band bounds from a max threshold.
/// Returns (min_err, max_err) where min_err has epsilon margin for floating-point precision.
pub fn error_band(max_err: f64) -> (f64, f64) {
    let min_err = (max_err - 0.011).max(0.0);
    (min_err, max_err)
}

/// Create a CSV reader with comment line support (lines starting with #).
/// Caller provides the CSV data as a byte slice.
pub fn csv_reader_with_comments(data: &[u8]) -> csv::Reader<&[u8]> {
    csv::ReaderBuilder::new()
        .comment(Some(b'#'))
        .from_reader(data)
}

