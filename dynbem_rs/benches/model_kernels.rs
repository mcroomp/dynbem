use criterion::{black_box, criterion_group, criterion_main, Criterion};
use dynbem_rs::aero_io::{Mat3, RotorInputs, Vec3};
use dynbem_rs::aero_model::AeroModel;
use dynbem_rs::bem_common::{element_force, ElementCtx, RadialGrid, SweepCtx};
use dynbem_rs::oye::{OyeBEMModel, OYE_K};
use dynbem_rs::pitt_peters::PittPetersModel;
use dynbem_rs::polar::{LinearPolar, TabulatedPolar};
use dynbem_rs::quasi_static_bem::{solve_bem_element, QuasiStaticBEM};
use dynbem_rs::rotor_definition::{BladeGeometry, LinearPolarParameters, RotorDefinition};
use dynbem_rs::rotor_state::{OyeRotorState, PittPetersRotorState, QuasiStaticRotorState};

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
        name: "bench_rotor".to_string(),
        description: "criterion baseline rotor".to_string(),
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
        omega_rad_s: 120.0,
    }
}

fn bench_solve_bem_element(c: &mut Criterion) {
    let polar = LinearPolar::new(0.0, 5.7, 0.01, 15.0_f64.to_radians());
    c.bench_function("solve_bem_element", |b| {
        b.iter(|| {
            black_box(solve_bem_element(
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
            ))
        })
    });
}

fn bench_element_force_only(c: &mut Criterion) {
    let defn = make_rotor_definition(30);
    let grid = RadialGrid::from_blade(&defn.blade);
    let linear = LinearPolar::new(0.0, 5.7, 0.01, 15.0_f64.to_radians());

    let n_tab = 121usize;
    let amin = -0.5;
    let amax = 0.5;
    let step = (amax - amin) / ((n_tab - 1) as f64);
    let mut alpha = Vec::with_capacity(n_tab);
    let mut cl = Vec::with_capacity(n_tab);
    let mut cd = Vec::with_capacity(n_tab);
    for i in 0..n_tab {
        let a = amin + step * (i as f64);
        alpha.push(a);
        if a.abs() < linear.alpha_stall_rad {
            cl.push(linear.CL0 + linear.CL_alpha_per_rad * a);
            cd.push(linear.CD0);
        } else {
            let cl_mag = linear.CL0 + linear.CL_alpha_per_rad * linear.alpha_stall_rad;
            cl.push(cl_mag.copysign(a));
            cd.push(linear.CD0 + (a.abs() - linear.alpha_stall_rad));
        }
    }
    let tabulated =
        TabulatedPolar::new(alpha, cl, cd).expect("tabulated benchmark polar must be valid");

    let omega = 120.0;
    let r = grid.r_mid[10];
    let chord = grid.chord[10];
    let twist = grid.twist_rad[10];
    let dr = grid.dr;
    let cos_psi = 0.766044443118978;
    let sin_psi = 0.642787609686539;
    let theta_1c = 1.0_f64.to_radians();
    let theta_1s = -0.8_f64.to_radians();
    let col = 8.0_f64.to_radians();
    let col_psi = col + theta_1c * cos_psi + theta_1s * sin_psi;
    let v_t_extra = 8.0 * sin_psi + 1.5 * cos_psi;
    let ctx = ElementCtx {
        i: 10,
        cos_psi,
        sin_psi,
        r,
        chord,
        twist,
        dr,
        col_psi,
        v_t: omega * r + v_t_extra,
    };
    let v_a = 0.065 * (omega * defn.blade.radius_m);

    let sweep_linear = SweepCtx {
        grid: &grid,
        polar: &linear,
        col,
        omega,
        omega_r: omega * defn.blade.radius_m,
        rho: 1.225,
        n_b: defn.blade.n_blades,
        n_psi: 1,
        n_psi_inv: 1.0,
        psi_trig: &[(cos_psi, sin_psi)],
        v_in_hub_x: 8.0,
        v_in_hub_y: 1.5,
        theta_1c,
        theta_1s,
    };
    let sweep_tab = SweepCtx {
        grid: &grid,
        polar: &tabulated,
        col,
        omega,
        omega_r: omega * defn.blade.radius_m,
        rho: 1.225,
        n_b: defn.blade.n_blades,
        n_psi: 1,
        n_psi_inv: 1.0,
        psi_trig: &[(cos_psi, sin_psi)],
        v_in_hub_x: 8.0,
        v_in_hub_y: 1.5,
        theta_1c,
        theta_1s,
    };

    let mut group = c.benchmark_group("element_force_only");
    group.bench_function("linear", |b| {
        b.iter(|| {
            black_box(element_force(
                black_box(v_a),
                black_box(&sweep_linear),
                black_box(&ctx),
            ))
        })
    });
    group.bench_function("tabulated", |b| {
        b.iter(|| {
            black_box(element_force(
                black_box(v_a),
                black_box(&sweep_tab),
                black_box(&ctx),
            ))
        })
    });
    group.finish();
}

fn bench_model_compute_forces(c: &mut Criterion) {
    let defn = make_rotor_definition(30);
    let polar = LinearPolar::new(0.0, 5.7, 0.01, 15.0_f64.to_radians());
    let inputs = make_inputs();

    let bem = QuasiStaticBEM::build(defn.clone(), 72, polar.clone());
    let pp = PittPetersModel::build(defn.clone(), 72, polar.clone());
    let oye = OyeBEMModel::build_with_k(defn.clone(), 72, polar.clone(), OYE_K);

    let bem_state = QuasiStaticRotorState;
    let pp_state = PittPetersRotorState {
        lambda_0: 0.06,
        lambda_c: 0.01,
        lambda_s: -0.008,
    };
    let oye_state = OyeRotorState::zeros(defn.blade.n_elements);

    let mut group = c.benchmark_group("models_compute_forces");
    group.bench_function("bem", |b| {
        b.iter(|| black_box(bem.compute_forces(&inputs, &bem_state)))
    });
    group.bench_function("pitt_peters", |b| {
        b.iter(|| black_box(pp.compute_forces(&inputs, &pp_state)))
    });
    group.bench_function("oye", |b| {
        b.iter(|| black_box(oye.compute_forces(&inputs, &oye_state)))
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_solve_bem_element,
    bench_element_force_only,
    bench_model_compute_forces
);
criterion_main!(benches);
