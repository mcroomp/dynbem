// Standalone kernel profiling harness for external profilers.
// Build: cargo build --release -p dynbem_rs
// Run: ./target/release/profile_kernels.exe [benchmark_name]
//
// Examples:
//   profile_kernels.exe solve_bem_element
//   profile_kernels.exe pitt_peters
//   profile_kernels.exe oye
//   profile_kernels.exe sweep

use dynbem_rs::aero_io::{Mat3, RotorInputs, Vec3};
use dynbem_rs::aero_model::AeroModel;
use dynbem_rs::bem_common::RadialGrid;
use dynbem_rs::oye::{OYE_K, OyeBEMModel};
use dynbem_rs::pitt_peters::PittPetersModel;
use dynbem_rs::polar::{LinearPolar, PolarKind};
use dynbem_rs::quasi_static_bem::{QuasiStaticBEM, solve_bem_element};
use dynbem_rs::rotor_definition::{
    AirfoilProperties, AutorotationProperties, BladeGeometry, InertiaProperties, RotorDefinition,
};
use dynbem_rs::rotor_state::{OyeRotorState, PittPetersRotorState, QuasiStaticRotorState};
use std::env;
use std::time::Instant;

fn make_rotor_definition(n_elements: usize) -> RotorDefinition {
    RotorDefinition {
        blade: BladeGeometry {
            n_blades: 2,
            radius_m: 1.0,
            root_cutout_m: 0.2,
            chord_m: 0.06,
            twist_deg: 2.0,
            n_elements,
            r_stations_m: Vec::new(),
            chord_stations_m: Vec::new(),
            twist_stations_deg: Vec::new(),
        },
        airfoil: AirfoilProperties {
            Re_design: 1_000_000,
            CL0: 0.0,
            CL_alpha_per_rad: 5.7,
            CD0: 0.01,
            alpha_stall_deg: 15.0,
            tip_loss: true,
            name: "bench_linear".to_string(),
            source: "profile_kernels".to_string(),
            polar_csv: None,
            CD_structural: 0.0,
            Re_operating: None,
        },
        control: None,
        inertia: InertiaProperties::default(),
        autorotation: AutorotationProperties::default(),
        name: "bench_rotor".to_string(),
        description: "standalone harness rotor".to_string(),
    }
}

fn make_inputs() -> RotorInputs {
    RotorInputs {
        collective_rad: 8.0_f64.to_radians(),
        tilt_lon: 1.0_f64.to_radians(),
        tilt_lat: -0.8_f64.to_radians(),
        R_hub: Mat3::eye(),
        v_hub_world: Vec3::new(0.0, 0.0, 0.0),
        wind_world: Vec3::new(8.0, 1.5, -1.0),
        t: 0.0,
        rho_kg_m3: 1.225,
        motor_torque_Nm: 0.0,
    }
}

fn bench_solve_bem_element(iterations: usize) {
    let polar = PolarKind::Linear(LinearPolar::new(0.0, 5.7, 0.01, 15.0_f64.to_radians()));
    let start = Instant::now();
    for _ in 0..iterations {
        let _ = solve_bem_element(
            0.85,
            0.02,
            0.06,
            2.0_f64.to_radians(),
            8.0_f64.to_radians(),
            120.0,
            -1.0,
            1.225,
            2,
            1.0,
            &polar,
            true,
            0.0,
            0.2,
        );
    }
    let elapsed = start.elapsed();
    println!(
        "solve_bem_element: {} iterations in {:.2}ms ({:.3}us per)",
        iterations,
        elapsed.as_secs_f64() * 1000.0,
        elapsed.as_secs_f64() * 1_000_000.0 / iterations as f64
    );
}

fn bench_pitt_peters(iterations: usize) {
    let defn = make_rotor_definition(30);
    let polar = PolarKind::Linear(LinearPolar::new(0.0, 5.7, 0.01, 15.0_f64.to_radians()));
    let inputs = make_inputs();

    let pp = PittPetersModel::build(defn, 72, polar);
    let pp_state = PittPetersRotorState {
        lambda_0: 0.06,
        lambda_c: 0.01,
        lambda_s: -0.008,
        omega_rad_s: 120.0,
        spin_angle_rad: 0.0,
    };

    let start = Instant::now();
    for _ in 0..iterations {
        let _ = pp.compute_forces(&inputs, &pp_state);
    }
    let elapsed = start.elapsed();
    println!(
        "pitt_peters: {} iterations in {:.2}ms ({:.3}us per)",
        iterations,
        elapsed.as_secs_f64() * 1000.0,
        elapsed.as_secs_f64() * 1_000_000.0 / iterations as f64
    );
}

fn bench_oye(iterations: usize) {
    let defn = make_rotor_definition(30);
    let polar = PolarKind::Linear(LinearPolar::new(0.0, 5.7, 0.01, 15.0_f64.to_radians()));
    let inputs = make_inputs();

    let oye = OyeBEMModel {
        defn: defn.clone(),
        n_psi_elements: 72,
        coupling_k: OYE_K,
        polar,
        grid: RadialGrid::from_blade(&defn.blade),
    };

    let oye_state = OyeRotorState::zeros(defn.blade.n_elements, 120.0);

    let start = Instant::now();
    for _ in 0..iterations {
        let _ = oye.compute_forces(&inputs, &oye_state);
    }
    let elapsed = start.elapsed();
    println!(
        "oye: {} iterations in {:.2}ms ({:.3}us per)",
        iterations,
        elapsed.as_secs_f64() * 1000.0,
        elapsed.as_secs_f64() * 1_000_000.0 / iterations as f64
    );
}

fn bench_sweep(iterations: usize) {
    let defn = make_rotor_definition(30);
    let polar = PolarKind::Linear(LinearPolar::new(0.0, 5.7, 0.01, 15.0_f64.to_radians()));
    let inputs = make_inputs();

    let bem = QuasiStaticBEM::build(defn, 72, polar);
    let bem_state = QuasiStaticRotorState {
        omega_rad_s: 120.0,
        spin_angle_rad: 0.0,
    };

    let start = Instant::now();
    for _ in 0..iterations {
        let _ = bem.compute_forces(&inputs, &bem_state);
    }
    let elapsed = start.elapsed();
    println!(
        "sweep (bem): {} iterations in {:.2}ms ({:.3}us per)",
        iterations,
        elapsed.as_secs_f64() * 1000.0,
        elapsed.as_secs_f64() * 1_000_000.0 / iterations as f64
    );
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let bench_name = if args.len() > 1 {
        args[1].as_str()
    } else {
        "all"
    };

    // Iteration counts chosen to give ~5-10 seconds of runtime per benchmark.
    // Adjust these based on your profiler's sampling window needs.
    match bench_name {
        "solve_bem_element" => bench_solve_bem_element(100_000),
        "pitt_peters" => bench_pitt_peters(100_000),
        "oye" => bench_oye(50_000),
        "sweep" => bench_sweep(5_000),
        "all" => {
            bench_solve_bem_element(100_000);
            bench_pitt_peters(10_000);
            bench_oye(5_000);
            bench_sweep(5_000);
        }
        _ => {
            eprintln!("Unknown benchmark: {}", bench_name);
            eprintln!("Available: solve_bem_element, pitt_peters, oye, sweep, all");
            std::process::exit(1);
        }
    }
}
