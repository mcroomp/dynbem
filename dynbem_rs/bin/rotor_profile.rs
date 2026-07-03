//! Rotor aerodynamic model profiler: timing comparisons across VPM (direct/BH),
//! Pitt-Peters, Øye, and quasi-static BEM models at varying particle counts.
//!
//! Usage:
//!   cargo build --release -p dynbem_rs --features parallel
//!   ./target/release/rotor_profile [MODELS] [OPTIONS]
//!
//! Models (comma-separated, case-insensitive):
//!   vpm-direct     - VPM direct O(N^2) velocity evaluation
//!   vpm-bh         - VPM Barnes-Hut O(N log N) tree evaluation (theta=0.5)
//!   pitt-peters    - Global 3-state inflow model (L-matrix)
//!   oye            - Per-annulus 2-stage filter model
//!   bem            - Quasi-static BEM with Brent root finder
//!   all            - All models (default)
//!
//! Particle counts (comma-separated integers):
//!   Default: 500,2000,5000,10000,20000 for VPM; 1 for other models (unused)
//!
//! Options:
//!   --max-particles <N>     Cap VPM max_particles per config (default: 20000)
//!   --duration <secs>       Wall-clock target per model/count (default: 10.0 s)
//!   --output <fmt>          csv, md, or text (default: text)
//!   --seq                   Force single-threaded evaluation (no Rayon)
//!   --help                  Show this message
//!
//! Example:
//!   rotor_profile vpm-direct,vpm-bh --duration 30
//!   rotor_profile all --output csv > profile.csv
//!   rotor_profile pitt-peters,oye --seq



use dynbem_rs::aero_io::{RotorInputs, Vec3, Mat3};
use dynbem_rs::aero_model::AeroModel;
use dynbem_rs::cyclic::ControlGains;
use dynbem_rs::oye::OyeBEMModel;
use dynbem_rs::pitt_peters::PittPetersModel;
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::quasi_static_bem::QuasiStaticBEM;
use dynbem_rs::rotor_definition::{
    BladeGeometry, ControlProperties, LinearPolarParameters, PitchActuation,
    RotorDefinition,
};
use dynbem_rs::vpm_rotor::{FlightCondition, VpmRotor, VpmRotorConfig};
use std::env;
use std::f64::consts::PI;
use std::time::{Duration, Instant};

// -----------

const ROTOR_RADIUS_M: f64 = 3.0;
const ROTOR_RPM: f64 = 120.0;

fn test_rotor_definition() -> RotorDefinition {
    RotorDefinition {
        blade: BladeGeometry {
            n_blades: 2,
            radius_m: ROTOR_RADIUS_M,
            root_cutout_m: 0.3,
            chord_m: 0.3,
            twist_deg: 0.0,
            n_elements: 20,
            tip_loss: true,
            r_stations_m: Vec::new(),
            chord_stations_m: Vec::new(),
            twist_stations_deg: Vec::new(),
        },
        airfoil: LinearPolarParameters {
            CL0: 0.0,
            CL_alpha_per_rad: 5.73,
            CD0: 0.01,
            alpha_stall_deg: 15.0,
        },
        control: Some(ControlProperties {
            swashplate_pitch_gain_rad: 1.0,
            swashplate_phase_deg: Some(0.0),
        }),
        pitch_actuation: PitchActuation::DirectMechanical,
        flap: None,
        name: "profile_test".to_string(),
        description: String::new(),
    }
}

fn test_flight_condition() -> FlightCondition {
    let omega = ROTOR_RPM * PI / 30.0;
    FlightCondition {
        omega_rad_s: omega,
        v_hub: [12.0, 5.0, -2.0], // forward flight + edgewise + descent
        collective_rad: 0.1,
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        rho: 1.225,
    }
}

fn flight_condition_to_rotor_inputs(fc: &FlightCondition) -> RotorInputs {
    // Identity hub orientation (no tilt); v_hub == v_hub_world
    RotorInputs {
        collective_rad: fc.collective_rad,
        tilt_lon: fc.tilt_lon,
        tilt_lat: fc.tilt_lat,
        R_hub: Mat3::eye(),
        v_hub_world: Vec3::new(fc.v_hub[0], fc.v_hub[1], fc.v_hub[2]),
        wind_world: Vec3::zero(),
        rho_kg_m3: fc.rho,
        omega_rad_s: fc.omega_rad_s,
    }
}

#[derive(Debug, Clone)]
#[allow(dead_code)]
struct ProfileResult {
    model: String,
    particle_count: usize,
    iterations: u64,
    elapsed_s: f64,
    ms_per_step: f64,
}

// -----

enum ModelKind {
    VpmDirect,
    VpmBh,
    PittPeters,
    Oye,
    Bem,
}

impl ModelKind {
    #[allow(dead_code)]
    fn name(&self) -> &'static str {
        match self {
            ModelKind::VpmDirect => "vpm-direct",
            ModelKind::VpmBh => "vpm-bh",
            ModelKind::PittPeters => "pitt-peters",
            ModelKind::Oye => "oye",
            ModelKind::Bem => "bem",
        }
    }

    fn from_str(s: &str) -> Option<Self> {
        match s.trim().to_lowercase().as_str() {
            "vpm-direct" | "vpm_direct" => Some(ModelKind::VpmDirect),
            "vpm-bh" | "vpm_bh" => Some(ModelKind::VpmBh),
            "pitt-peters" | "pitt_peters" | "pp" => Some(ModelKind::PittPeters),
            "oye" => Some(ModelKind::Oye),
            "bem" => Some(ModelKind::Bem),
            "all" => None, // special case: caller should expand
            _ => None,
        }
    }
}

// -----

fn profile_vpm(
    max_particles: usize,
    barnes_hut: bool,
    _force_seq: bool,
    duration: Duration,
) -> ProfileResult {
    let defn = test_rotor_definition();
    let polar = LinearPolar::from_properties(&defn.airfoil);
    let config = VpmRotorConfig {
        max_particles,
        sigma: 0.15,
        relax: 0.35,
        nonlinear_lifting_line: true,
        tip_clustering: true,
        local_core: true,
        barnes_hut,
        bh_theta: 0.5,
        bh_min_particles: 200,
        flap_dynamics: false,
        use_rayon: true,
    };
    let rotor = VpmRotor::new(&defn, polar, ControlGains::default(), config);
    let fc = test_flight_condition();

    // One full revolution; with 18 steps/rev we get 18 substeps.
    let omega = fc.omega_rad_s;
    let dpsi_per_step = 2.0 * PI / 18.0;
    let dt = dpsi_per_step / omega;

    let model_label = if barnes_hut {
        "vpm-bh".to_string()
    } else {
        "vpm-direct".to_string()
    };

    eprintln!("Profiling {}...", model_label);

    // Warm up
    let (_, mut state) = rotor.march(&fc, None, dt, 18);
    let mut checksum = 0.0;

    // Measure
    let start = Instant::now();
    let mut iters = 0u64;
    while start.elapsed() < duration {
        let (_, new_state) = rotor.step_one(&fc, &state, dt);
        state = new_state;
        checksum += state.psi;
        iters += 1;
    }
    let elapsed = start.elapsed().as_secs_f64();
    let ms_per_step = 1e3 * elapsed / iters as f64;

    eprintln!(
        "{}: {} steps in {:.2}s = {:.2} ms/step (checksum {:.6})",
        model_label, iters, elapsed, ms_per_step, checksum
    );

    ProfileResult {
        model: model_label,
        particle_count: state.wake.as_ref().map(|w| w.len()).unwrap_or(0),
        iterations: iters,
        elapsed_s: elapsed,
        ms_per_step,
    }
}

fn profile_pitt_peters(_force_seq: bool, duration: Duration) -> ProfileResult {
    let defn = test_rotor_definition();
    let polar = LinearPolar::from_properties(&defn.airfoil);
    let model = PittPetersModel::build(defn, 72, polar);
    let fc = flight_condition_to_rotor_inputs(&test_flight_condition());
    let mut state = model.initial_state();

    eprintln!("Profiling pitt-peters...");

    // Warm up
    let _ = model.compute_forces(&fc, &state);

    let mut checksum = 0.0;
    let start = Instant::now();
    let mut iters = 0u64;
    while start.elapsed() < duration {
        let (aero_result, new_state) = model.compute_forces(&fc, &state);
        state = new_state;
        checksum += aero_result.F_world[0];
        iters += 1;
    }
    let elapsed = start.elapsed().as_secs_f64();
    let ms_per_step = 1e3 * elapsed / iters as f64;

    eprintln!(
        "pitt-peters: {} steps in {:.2}s = {:.2} ms/step (checksum {:.6})",
        iters, elapsed, ms_per_step, checksum
    );

    ProfileResult {
        model: "pitt-peters".to_string(),
        particle_count: 0,
        iterations: iters,
        elapsed_s: elapsed,
        ms_per_step,
    }
}

fn profile_oye(_force_seq: bool, duration: Duration) -> ProfileResult {
    let defn = test_rotor_definition();
    let polar = LinearPolar::from_properties(&defn.airfoil);
    let model = OyeBEMModel::build(defn, 72, polar);
    let fc = flight_condition_to_rotor_inputs(&test_flight_condition());
    let mut state = model.initial_state();

    eprintln!("Profiling oye...");

    // Warm up
    let _ = model.compute_forces(&fc, &state);

    let mut checksum = 0.0;
    let start = Instant::now();
    let mut iters = 0u64;
    while start.elapsed() < duration {
        let (aero_result, new_state) = model.compute_forces(&fc, &state);
        state = new_state;
        checksum += aero_result.F_world[0];
        iters += 1;
    }
    let elapsed = start.elapsed().as_secs_f64();
    let ms_per_step = 1e3 * elapsed / iters as f64;

    eprintln!(
        "oye: {} steps in {:.2}s = {:.2} ms/step (checksum {:.6})",
        iters, elapsed, ms_per_step, checksum
    );

    ProfileResult {
        model: "oye".to_string(),
        particle_count: 0,
        iterations: iters,
        elapsed_s: elapsed,
        ms_per_step,
    }
}

fn profile_bem(_force_seq: bool, duration: Duration) -> ProfileResult {
    let defn = test_rotor_definition();
    let polar = LinearPolar::from_properties(&defn.airfoil);
    let model = QuasiStaticBEM::build(defn, 72, polar);
    let fc = flight_condition_to_rotor_inputs(&test_flight_condition());
    let mut state = model.initial_state();

    eprintln!("Profiling bem...");

    // Warm up
    let _ = model.compute_forces(&fc, &state);

    let mut checksum = 0.0;
    let start = Instant::now();
    let mut iters = 0u64;
    while start.elapsed() < duration {
        let (aero_result, new_state) = model.compute_forces(&fc, &state);
        state = new_state;
        checksum += aero_result.F_world[0];
        iters += 1;
    }
    let elapsed = start.elapsed().as_secs_f64();
    let ms_per_step = 1e3 * elapsed / iters as f64;

    eprintln!(
        "bem: {} steps in {:.2}s = {:.2} ms/step (checksum {:.6})",
        iters, elapsed, ms_per_step, checksum
    );

    ProfileResult {
        model: "bem".to_string(),
        particle_count: 0,
        iterations: iters,
        elapsed_s: elapsed,
        ms_per_step,
    }
}

// -----

fn print_text(results: &[ProfileResult]) {
    println!();
    println!("Rotor profiling results");
    println!("Rotor: R={:.1}m, rpm={:.0}, blades={}", 
             ROTOR_RADIUS_M, ROTOR_RPM, 2);
    println!();
    println!(
        "{:<15} {:<15} {:>12} {:>12}",
        "Model", "Particles", "Iterations", "ms/step"
    );
    println!("{}", "-".repeat(60));
    for r in results {
        let part_str = if r.particle_count == 0 {
            "-".to_string()
        } else {
            format!("{}", r.particle_count)
        };
        println!(
            "{:<15} {:<15} {:>12} {:>12.3}",
            r.model, part_str, r.iterations, r.ms_per_step
        );
    }
    println!();
}

fn print_csv(results: &[ProfileResult]) {
    println!("Model,Particles,Iterations,ms_per_step");
    for r in results {
        println!(
            "{},{},{},{}",
            r.model,
            if r.particle_count == 0 { 0 } else { r.particle_count },
            r.iterations,
            r.ms_per_step
        );
    }
}

fn print_markdown(results: &[ProfileResult]) {
    println!();
    println!("| Model | Particles | Iterations | ms/step |");
    println!("|-------|-----------|------------|---------|");
    for r in results {
        let part_str = if r.particle_count == 0 {
            "-".to_string()
        } else {
            format!("{}", r.particle_count)
        };
        println!(
            "| {} | {} | {} | {:.3} |",
            r.model, part_str, r.iterations, r.ms_per_step
        );
    }
    println!();
}

fn print_help() {
    eprintln!(
        r#"rotor_profile: Aerodynamic model profiler for VPM (direct/BH), Pitt-Peters, Øye, and BEM

USAGE:
  rotor_profile [MODELS] [OPTIONS]

MODELS (comma-separated, case-insensitive):
  vpm-direct     - VPM direct O(N²) velocity evaluation
  vpm-bh         - VPM Barnes-Hut O(N log N) tree (theta=0.5)
  pitt-peters    - Global 3-state Pitt-Peters inflow model
  oye            - Per-annulus 2-stage Øye filter model
  bem            - Quasi-static BEM with Brent solver
  all            - All models (default if none specified)

OPTIONS:
  --max-particles <N>    Cap VPM max_particles (default: 20000)
  --duration <secs>      Wall-clock time per test (default: 10.0 s)
  --output <fmt>         csv, md, or text (default: text)
  --seq                  Force single-threaded evaluation
  --help, -h             Show this message

EXAMPLES:
  rotor_profile all --duration 10
  rotor_profile vpm-direct,vpm-bh --duration 30 --output csv
  rotor_profile pitt-peters,oye,bem --output csv > results.csv
"#
    );
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let mut models = vec![];
    let particle_counts = vec![500, 2000, 5000, 10000, 20000];
    let mut max_particles = 20000;
    let mut duration = Duration::from_secs_f64(10.0);
    let mut output_fmt = "text";
    let mut force_seq = false;

    let mut i = 1;
    while i < args.len() {
        let arg = &args[i];

        if arg == "--help" || arg == "-h" {
            print_help();
            return;
        } else if arg == "--seq" {
            force_seq = true;
            i += 1;
        } else if arg == "--max-particles" {
            i += 1;
            if i < args.len() {
                max_particles = args[i].parse().unwrap_or(20000);
            }
            i += 1;
        } else if arg == "--duration" {
            i += 1;
            if i < args.len() {
                if let Ok(secs) = args[i].parse::<f64>() {
                    duration = Duration::from_secs_f64(secs);
                }
            }
            i += 1;
        } else if arg == "--output" {
            i += 1;
            if i < args.len() {
                output_fmt = match args[i].as_str() {
                    "csv" => "csv",
                    "md" | "markdown" => "markdown",
                    _ => "text",
                };
            }
            i += 1;
        } else if !arg.starts_with("--") {
            // Parse model list
            if let Some(m) = ModelKind::from_str(arg) {
                models.push(m);
            } else if arg.to_lowercase() == "all" {
                models = vec![
                    ModelKind::VpmDirect,
                    ModelKind::VpmBh,
                    ModelKind::PittPeters,
                    ModelKind::Oye,
                    ModelKind::Bem,
                ];
            } else {
                // Try parsing as comma-separated list
                for part in arg.split(',') {
                    if let Some(m) = ModelKind::from_str(part) {
                        models.push(m);
                    }
                }
            }
            i += 1;
        } else {
            i += 1;
        }
    }

    // Default to all models
    if models.is_empty() {
        models = vec![
            ModelKind::VpmDirect,
            ModelKind::VpmBh,
            ModelKind::PittPeters,
            ModelKind::Oye,
            ModelKind::Bem,
        ];
    }

    let mut results = vec![];

    for model in models {
        match model {
            ModelKind::VpmDirect => {
                for &n in &particle_counts {
                    let r = profile_vpm(n.min(max_particles), false, force_seq, duration);
                    results.push(r);
                }
            }
            ModelKind::VpmBh => {
                for &n in &particle_counts {
                    let r = profile_vpm(n.min(max_particles), true, force_seq, duration);
                    results.push(r);
                }
            }
            ModelKind::PittPeters => {
                let r = profile_pitt_peters(force_seq, duration);
                results.push(r);
            }
            ModelKind::Oye => {
                let r = profile_oye(force_seq, duration);
                results.push(r);
            }
            ModelKind::Bem => {
                let r = profile_bem(force_seq, duration);
                results.push(r);
            }
        }
    }

    match output_fmt {
        "csv" => print_csv(&results),
        "markdown" => print_markdown(&results),
        _ => print_text(&results),
    }
}
