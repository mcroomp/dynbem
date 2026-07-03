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
use dynbem_rs::oye::OyeRotorState;
use dynbem_rs::oye::{OyeBEMModel, OYE_K};
use dynbem_rs::pitt_peters::PittPetersModel;
use dynbem_rs::pitt_peters::PittPetersRotorState;
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::quasi_static_bem::QuasiStaticRotorState;
use dynbem_rs::quasi_static_bem::{solve_bem_element, BEMElementGeometry, QuasiStaticBEM};
use dynbem_rs::rotor_definition::{
    BladeGeometry, LinearPolarParameters, PitchActuation, RotorDefinition,
};
use dynbem_rs::vpm_rotor::{advect_rk2, induced_velocities, ParticleField};
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
            tip_loss: true,
            r_stations_m: Vec::new(),
            chord_stations_m: Vec::new(),
            twist_stations_deg: Vec::new(),
        },
        airfoil: LinearPolarParameters {
            CL0: 0.0,
            CL_alpha_per_rad: 5.7,
            CD0: 0.01,
            alpha_stall_deg: 15.0,
        },
        control: None,
        pitch_actuation: PitchActuation::DirectMechanical,
        flap: None,
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
        rho_kg_m3: 1.225,
        omega_rad_s: 120.0,
    }
}

fn bench_solve_bem_element(iterations: usize) {
    let polar = LinearPolar::new(0.0, 5.7, 0.01, 15.0_f64.to_radians());
    let geom = BEMElementGeometry::new(
        0.85,
        0.02,
        0.06,
        2.0_f64.to_radians(),
        120.0,
        1.225,
        2,
        1.0,
        &polar,
        true,
        0.2,
    );
    let start = Instant::now();
    for _ in 0..iterations {
        let _ = solve_bem_element(&geom, 8.0_f64.to_radians(), -1.0, 0.0);
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
    let polar = LinearPolar::new(0.0, 5.7, 0.01, 15.0_f64.to_radians());
    let inputs = make_inputs();

    let pp = PittPetersModel::build(defn, 72, polar);
    let pp_state = PittPetersRotorState {
        lambda_0: 0.06,
        lambda_c: 0.01,
        lambda_s: -0.008,
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
    let polar = LinearPolar::new(0.0, 5.7, 0.01, 15.0_f64.to_radians());
    let inputs = make_inputs();

    let oye = OyeBEMModel::build_with_k(defn.clone(), 72, polar, OYE_K);

    let oye_state = OyeRotorState::zeros(defn.blade.n_elements);

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
    let polar = LinearPolar::new(0.0, 5.7, 0.01, 15.0_f64.to_radians());
    let inputs = make_inputs();

    let bem = QuasiStaticBEM::build(defn, 72, polar);
    let bem_state = QuasiStaticRotorState;

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

/// Deterministic pseudo-random particle cloud (LCG, no rand dependency).
fn seed_cloud(n: usize) -> ParticleField {
    let mut state = 0x9E37_79B9u32;
    let mut rng = || {
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        (state >> 8) as f32 / (1u32 << 24) as f32 - 0.5
    };
    let mut f = ParticleField::with_capacity(n);
    for _ in 0..n {
        f.push(
            [rng() * 6.0, rng() * 6.0, rng() * 6.0],
            [rng() * 0.1, rng() * 0.1, rng() * 0.1],
            0.15,
        );
    }
    f
}

fn bench_vpm() {
    // One "solver step" = one advect_rk2 (two O(N^2) velocity evaluations).
    // Report per-step wall time across a range of particle counts, plus the
    // per-velocity-evaluation cost. Iteration count shrinks as N grows to
    // keep each measurement ~1s.
    println!("VPM direct O(N^2), SIMD (wide::f32x8), single-threaded:");
    println!(
        "{:>8}  {:>12}  {:>14}  {:>16}",
        "N", "steps", "us/eval", "ms/advect_step"
    );
    for &n in &[500usize, 2_000, 5_000, 10_000, 20_000] {
        let mut field = seed_cloud(n);
        // Scale iterations so total work ~ constant.
        let iters = (2_000_000_000 / (n * n).max(1)).clamp(3, 2000);

        // Time a bare velocity evaluation.
        let start = Instant::now();
        for _ in 0..iters {
            let v = induced_velocities(&field);
            std::hint::black_box(&v);
        }
        let per_eval_us = start.elapsed().as_secs_f64() * 1e6 / iters as f64;

        // Time a full RK2 advect step (two evaluations + position updates).
        let start = Instant::now();
        for _ in 0..iters {
            advect_rk2(&mut field, [8.0, 0.0, -1.0], 1e-4);
        }
        let per_step_ms = start.elapsed().as_secs_f64() * 1e3 / iters as f64;

        println!(
            "{:>8}  {:>12}  {:>14.2}  {:>16.4}",
            n, iters, per_eval_us, per_step_ms
        );
    }
}

fn bench_compare() {
    // One call of each BEM-family model (fixed cost, independent of wake
    // history) against one VPM velocity evaluation at several N. These are
    // different units of work -- a BEM "call" is a full converged rotor
    // solution, a VPM "eval" is one wake-convection velocity field -- so read
    // the table as "cost of one solver step", not as equal fidelity.
    let n_iter = 20_000;
    println!("Per-call cost of one solver step (single-threaded):\n");

    let defn = make_rotor_definition(30);
    let polar = LinearPolar::new(0.0, 5.7, 0.01, 15.0_f64.to_radians());
    let inputs = make_inputs();

    let bem = QuasiStaticBEM::build(defn.clone(), 72, polar.clone());
    let bem_state = QuasiStaticRotorState;
    let t = Instant::now();
    for _ in 0..n_iter {
        std::hint::black_box(bem.compute_forces(&inputs, &bem_state));
    }
    println!(
        "  BEM (quasi-static, 30 elem x 72 psi)   {:>10.3} us/call",
        t.elapsed().as_secs_f64() * 1e6 / n_iter as f64
    );

    let pp = PittPetersModel::build(defn.clone(), 72, polar.clone());
    let pp_state = PittPetersRotorState {
        lambda_0: 0.06,
        lambda_c: 0.01,
        lambda_s: -0.008,
    };
    let t = Instant::now();
    for _ in 0..n_iter {
        std::hint::black_box(pp.compute_forces(&inputs, &pp_state));
    }
    println!(
        "  Pitt-Peters (30 elem x 72 psi)         {:>10.3} us/call",
        t.elapsed().as_secs_f64() * 1e6 / n_iter as f64
    );

    let oye = OyeBEMModel::build_with_k(defn.clone(), 72, polar, OYE_K);
    let oye_state = OyeRotorState::zeros(defn.blade.n_elements);
    let t = Instant::now();
    for _ in 0..n_iter {
        std::hint::black_box(oye.compute_forces(&inputs, &oye_state));
    }
    println!(
        "  Oye (30 elem x 72 psi)                 {:>10.3} us/call",
        t.elapsed().as_secs_f64() * 1e6 / n_iter as f64
    );

    println!();
    for &n in &[500usize, 2_000, 5_000, 10_000] {
        let field = seed_cloud(n);
        let iters = (2_000_000_000 / (n * n).max(1)).clamp(3, 2000);
        let t = Instant::now();
        for _ in 0..iters {
            std::hint::black_box(induced_velocities(&field));
        }
        println!(
            "  VPM velocity eval, N = {:>6}          {:>10.3} us/eval",
            n,
            t.elapsed().as_secs_f64() * 1e6 / iters as f64
        );
    }
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
        "vpm" => bench_vpm(),
        "compare" => bench_compare(),
        "all" => {
            bench_solve_bem_element(100_000);
            bench_pitt_peters(10_000);
            bench_oye(5_000);
            bench_sweep(5_000);
            bench_vpm();
        }
        _ => {
            eprintln!("Unknown benchmark: {}", bench_name);
            eprintln!(
                "Available: solve_bem_element, pitt_peters, oye, sweep, vpm, compare, all"
            );
            std::process::exit(1);
        }
    }
}
