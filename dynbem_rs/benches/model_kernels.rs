use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use dynbem_rs::aero_io::{Mat3, RotorInputs, Vec3};
use dynbem_rs::aero_model::AeroModel;
use dynbem_rs::bem_common::{element_force, ElementCtx, PsiKernel, RadialGrid, SweepCtx};
use dynbem_rs::oye::{OyeBEMModel, OYE_K};
use dynbem_rs::pitt_peters::PittPetersModel;
use dynbem_rs::polar::{LinearPolar, PolarKind};
use dynbem_rs::quasi_static_bem::{solve_bem_element, QuasiStaticBEM};
use dynbem_rs::rotor_definition::{
    AirfoilProperties, AutorotationProperties, BladeGeometry, RotorDefinition,
};
use dynbem_rs::rotor_state::{OyeRotorState, PittPetersRotorState, QuasiStaticRotorState};

struct PrescribedKernel<'a> {
    lambda_total: f64,
    lam_c: f64,
    lam_s: f64,
    x_mid: &'a [f64],
}

impl<'a> PsiKernel for PrescribedKernel<'a> {
    #[inline(always)]
    fn element(&mut self, sweep: &SweepCtx<'_>, ctx: &ElementCtx) -> (f64, f64) {
        let lam = self.lambda_total
            + self.x_mid[ctx.i] * (self.lam_c * ctx.cos_psi + self.lam_s * ctx.sin_psi);
        element_force(lam * sweep.omega_r, sweep, ctx)
    }
}

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
            CL0: 0.0,
            CL_alpha_per_rad: 5.7,
            CD0: 0.01,
            alpha_stall_deg: 15.0,
            tip_loss: true,
        },
        control: None,
        autorotation: AutorotationProperties::default(),
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
        motor_torque_Nm: 0.0,
    }
}

fn bench_solve_bem_element(c: &mut Criterion) {
    let polar = PolarKind::Linear(LinearPolar::new(0.0, 5.7, 0.01, 15.0_f64.to_radians()));
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

fn bench_scalar_sweep(c: &mut Criterion) {
    let polar = PolarKind::Linear(LinearPolar::new(0.0, 5.7, 0.01, 15.0_f64.to_radians()));
    let mut group = c.benchmark_group("sweep_scalar");

    for &(n_r, n_psi) in &[(30usize, 36usize), (30, 72), (80, 72)] {
        let defn = make_rotor_definition(n_r);
        let grid = RadialGrid::from_blade(&defn.blade);
        let omega = 120.0;
        let sweep = SweepCtx {
            grid: &grid,
            polar: &polar,
            col: 8.0_f64.to_radians(),
            omega,
            omega_r: omega * defn.blade.radius_m,
            rho: 1.225,
            n_b: defn.blade.n_blades,
            n_psi,
            n_psi_inv: 1.0 / (n_psi as f64),
            v_in_hub_x: 8.0,
            v_in_hub_y: 1.5,
            theta_1c: 1.0_f64.to_radians(),
            theta_1s: -0.8_f64.to_radians(),
        };
        let id = format!("npsi{}_nr{}", n_psi, n_r);
        group.bench_function(BenchmarkId::new("prescribed", id), |b| {
            b.iter(|| {
                let mut kernel = PrescribedKernel {
                    lambda_total: 0.065,
                    lam_c: 0.012,
                    lam_s: -0.009,
                    x_mid: &grid.x_mid[..grid.n_elements],
                };
                black_box(sweep.run(&mut kernel))
            })
        });
    }

    group.finish();
}

fn bench_model_compute_forces(c: &mut Criterion) {
    let defn = make_rotor_definition(30);
    let polar = PolarKind::Linear(LinearPolar::new(0.0, 5.7, 0.01, 15.0_f64.to_radians()));
    let inputs = make_inputs();

    let bem = QuasiStaticBEM::build(defn.clone(), 72, polar.clone());
    let pp = PittPetersModel::build(defn.clone(), 72, polar.clone());
    let oye = OyeBEMModel {
        defn: defn.clone(),
        n_psi_elements: 72,
        coupling_k: OYE_K,
        polar: polar.clone(),
        grid: RadialGrid::from_blade(&defn.blade),
    };

    let bem_state = QuasiStaticRotorState {
        omega_rad_s: 120.0,
        spin_angle_rad: 0.0,
    };
    let pp_state = PittPetersRotorState {
        lambda_0: 0.06,
        lambda_c: 0.01,
        lambda_s: -0.008,
        omega_rad_s: 120.0,
        spin_angle_rad: 0.0,
    };
    let oye_state = OyeRotorState::zeros(defn.blade.n_elements, 120.0);

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
    bench_scalar_sweep,
    bench_model_compute_forces
);
criterion_main!(benches);
