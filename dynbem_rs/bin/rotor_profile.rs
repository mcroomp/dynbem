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
//!   --seq                   Time only the sequential path (default: time both
//!                           Rayon-parallel `-par` and sequential `-seq`, so the
//!                           multi-core speedup shows up as paired VPM rows)
//!   --par                   Time only the Rayon-parallel path (fast; avoids the
//!                           slow sequential O(N^2) direct sum at large N)
//!   --long                  Extend the particle sweep to large N
//!                           (2k,5k,10k,16k,32k) to show direct vs BH growth;
//!                           raises --max-particles to 32000
//!   --help                  Show this message
//!
//! Example:
//!   rotor_profile vpm-direct,vpm-bh --duration 30
//!   rotor_profile all --output csv > profile.csv
//!   rotor_profile pitt-peters,oye --seq

use dynbem_rs::aero_io::{Mat3, RotorInputs, Vec3};
use dynbem_rs::aero_model::AeroModel;
use dynbem_rs::cyclic::ControlGains;
use dynbem_rs::oye::OyeBEMModel;
use dynbem_rs::pitt_peters::PittPetersModel;
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::quasi_static_bem::QuasiStaticBEM;
use dynbem_rs::rotor_definition::{
    BladeGeometry, ControlProperties, LinearPolarParameters, PitchActuation, RotorDefinition,
};
use dynbem_rs::vpm_rotor::{FlightCondition, VpmRotor, VpmRotorConfig, VpmRotorState};
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

/// Build a VPM rotor at the given particle cap. `barnes_hut` selects the
/// evaluator and `use_rayon` selects the parallel (multi-core) vs sequential
/// velocity evaluation; every other parameter is identical, so all runs are a
/// like-for-like comparison.
fn build_vpm_rotor(
    max_particles: usize,
    barnes_hut: bool,
    use_rayon: bool,
) -> VpmRotor<LinearPolar> {
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
        use_rayon,
        use_scalar_nan_check: false,
    };
    VpmRotor::new(&defn, polar, ControlGains::default(), config)
}

/// Convection step size: one full revolution in 18 steps.
fn vpm_dt(fc: &FlightCondition) -> f64 {
    (2.0 * PI / 18.0) / fc.omega_rad_s
}

/// Settle a wake up to the `max_particles` cap so timing runs at a fixed N.
/// Marches until the wake reaches the cap or stops growing (whichever first),
/// then returns the steady FIFO state. Using the direct (exact) evaluator here
/// means the settled wake is identical regardless of which evaluator is timed.
fn settle_vpm_to_cap(
    rotor: &VpmRotor<LinearPolar>,
    fc: &FlightCondition,
    dt: f64,
    max_particles: usize,
) -> VpmRotorState {
    let wake_len = |s: &VpmRotorState| s.wake.as_ref().map(|w| w.len()).unwrap_or(0);
    let (_, mut state) = rotor.march(fc, None, dt, 4);
    let mut last_len = wake_len(&state);
    let mut stall = 0u32;
    // Hard cap on settle steps so an unreachable target can't loop forever.
    for _ in 0..5000 {
        if wake_len(&state) >= max_particles {
            break;
        }
        let (_, s) = rotor.step_one(fc, &state, dt);
        state = s;
        let len = wake_len(&state);
        if len == last_len {
            stall += 1;
            if stall > 20 {
                break; // wake has plateaued below the cap
            }
        } else {
            stall = 0;
        }
        last_len = len;
    }
    state
}

/// Time one VPM evaluator starting from a pre-settled wake `state` (held at the
/// cap, so N is constant through the run). Returns the per-step cost.
fn time_vpm_from(
    rotor: &VpmRotor<LinearPolar>,
    fc: &FlightCondition,
    dt: f64,
    mut state: VpmRotorState,
    duration: Duration,
    model_label: &str,
) -> ProfileResult {
    eprintln!("Profiling {}...", model_label);
    let mut checksum = 0.0;
    let start = Instant::now();
    let mut iters = 0u64;
    while start.elapsed() < duration {
        let (_, new_state) = rotor.step_one(fc, &state, dt);
        state = new_state;
        checksum += state.psi;
        iters += 1;
    }
    let elapsed = start.elapsed().as_secs_f64();
    let ms_per_step = 1e3 * elapsed / iters as f64;
    let particle_count = state.wake.as_ref().map(|w| w.len()).unwrap_or(0);

    eprintln!(
        "{}: {} steps in {:.2}s = {:.2} ms/step (N={}, checksum {:.6})",
        model_label, iters, elapsed, ms_per_step, particle_count, checksum
    );

    ProfileResult {
        model: model_label.to_string(),
        particle_count,
        iterations: iters,
        elapsed_s: elapsed,
        ms_per_step,
    }
}

/// Which parallelism variants to time for each VPM evaluator.
#[derive(Clone, Copy, PartialEq)]
enum TimingMode {
    /// Sequential only (`-seq` rows).
    SeqOnly,
    /// Rayon-parallel only (`-par` rows).
    ParOnly,
    /// Both, parallel first then sequential (default), so the multi-core
    /// speedup is a `-par`/`-seq` row pair.
    Both,
}

impl TimingMode {
    /// (use_rayon, label suffix) pairs to time for this mode.
    fn variants(self) -> &'static [(bool, &'static str)] {
        match self {
            TimingMode::SeqOnly => &[(false, "-seq")],
            TimingMode::ParOnly => &[(true, "-par")],
            TimingMode::Both => &[(true, "-par"), (false, "-seq")],
        }
    }
}

/// Profile the requested VPM evaluators at one particle cap. The wake is
/// settled ONCE (with the exact direct evaluator) to the cap, then each
/// requested evaluator is timed from that same state -- so all rows are
/// measured on an identical wake at an identical N. `mode` selects sequential,
/// parallel, or both, so the multi-core speedup is visible in the table.
fn profile_vpm_at(
    max_particles: usize,
    want_direct: bool,
    want_bh: bool,
    mode: TimingMode,
    duration: Duration,
) -> Vec<ProfileResult> {
    let fc = test_flight_condition();
    let dt = vpm_dt(&fc);

    // Settle with the fast parallel direct evaluator; the wake is identical
    // regardless of how the subsequent timing runs are evaluated.
    let settler = build_vpm_rotor(max_particles, false, true);
    let settled = settle_vpm_to_cap(&settler, &fc, dt, max_particles);

    // (barnes_hut, want) pairs to time.
    let evaluators = [
        (false, want_direct, "vpm-direct"),
        (true, want_bh, "vpm-bh"),
    ];
    let modes = mode.variants();

    let mut out = Vec::new();
    for &(barnes_hut, want, base) in &evaluators {
        if !want {
            continue;
        }
        for &(use_rayon, suffix) in modes {
            let rotor = build_vpm_rotor(max_particles, barnes_hut, use_rayon);
            let label = format!("{base}{suffix}");
            out.push(time_vpm_from(
                &rotor,
                &fc,
                dt,
                settled.clone(),
                duration,
                &label,
            ));
        }
    }
    out
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
    println!(
        "Rotor: R={:.1}m, rpm={:.0}, blades={}",
        ROTOR_RADIUS_M, ROTOR_RPM, 2
    );
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
            if r.particle_count == 0 {
                0
            } else {
                r.particle_count
            },
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
  --seq                  Time only the sequential path (default: both -par and -seq)
  --par                  Time only the Rayon-parallel path (fast at large N)
  --long                 Extend sweep to large N (2k,5k,10k,16k,32k); cap -> 32000
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
    let mut particle_counts = vec![500, 2000, 5000, 10000, 20000];
    let mut max_particles = 20000;
    let mut duration = Duration::from_secs_f64(10.0);
    let mut output_fmt = "text";
    let mut force_seq = false;
    let mut force_par = false;
    let mut long_sweep = false;

    let mut i = 1;
    while i < args.len() {
        let arg = &args[i];

        if arg == "--help" || arg == "-h" {
            print_help();
            return;
        } else if arg == "--seq" {
            force_seq = true;
            i += 1;
        } else if arg == "--par" {
            force_par = true;
            i += 1;
        } else if arg == "--long" {
            long_sweep = true;
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

    // `--long`: extend the sweep to large N (16k, 32k) so the O(N^2) direct sum
    // vs O(N log N) Barnes-Hut growth is visible. Raises the cap to match unless
    // the user set a larger one explicitly.
    if long_sweep {
        particle_counts = vec![2000, 5000, 10000, 16000, 32000];
        max_particles = max_particles.max(32000);
    }

    let mut results = vec![];

    // Parallelism variants to time: --seq or --par restrict to one; default both.
    let timing_mode = match (force_seq, force_par) {
        (true, false) => TimingMode::SeqOnly,
        (false, true) => TimingMode::ParOnly,
        _ => TimingMode::Both,
    };

    // VPM: settle each wake once (to the cap) and time the requested evaluators
    // on that same state, so vpm-direct and vpm-bh compare at an identical N.
    let want_direct = models.iter().any(|m| matches!(m, ModelKind::VpmDirect));
    let want_bh = models.iter().any(|m| matches!(m, ModelKind::VpmBh));
    if want_direct || want_bh {
        for &n in &particle_counts {
            let mut rows = profile_vpm_at(
                n.min(max_particles),
                want_direct,
                want_bh,
                timing_mode,
                duration,
            );
            results.append(&mut rows);
        }
    }

    for model in models {
        match model {
            // VPM handled above (paired settle/compare).
            ModelKind::VpmDirect | ModelKind::VpmBh => {}
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
