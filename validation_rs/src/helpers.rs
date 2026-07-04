//! Rotor geometry, flight-condition builders, and theory helpers for validation checks.

use dynbem_rs::cyclic::ControlGains;
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::rotor_definition::{
    BladeGeometry, FlapProperties, LinearPolarParameters, PitchActuation, RotorDefinition,
};
pub use dynbem_rs::rotor_definition::{ServoFlapActuation, ServoFlapGeometry};
use dynbem_rs::vpm_rotor::{
    induced_velocities_at_points, FlightCondition, VpmRotor, VpmRotorConfig, VpmRotorResult,
    VpmRotorState,
};
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

pub const R_TIP: f64 = 0.914;
pub const R_ROOT: f64 = 0.155;
pub const CHORD: f64 = 0.0479;
pub const N_BLADES: usize = 3;
pub const CL_ALPHA: f64 = 5.90;
pub const CD0: f64 = 0.01046;
pub const ALPHA_STALL_DEG: f64 = 15.5;
pub const RHO: f64 = 1.225;
pub const OMEGA: f64 = 125.6637; // 1200 rpm
pub const STEPS_PER_REV: usize = 24;
pub const PCA2_R: f64 = 6.85;

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

pub fn solidity() -> f64 {
    N_BLADES as f64 * CHORD / (PI * R_TIP)
}

pub fn tip_speed() -> f64 {
    OMEGA * R_TIP
}

pub fn omega_from_rpm(rpm: f64) -> f64 {
    rpm * PI / 30.0
}

pub fn ct_at(thrust: f64, omega: f64, r: f64) -> f64 {
    thrust / (RHO * PI * r * r * (omega * r).powi(2))
}

pub fn cq_at(torque: f64, omega: f64, r: f64) -> f64 {
    torque / (RHO * PI * r * r * (omega * r).powi(2) * r)
}

pub fn c_t(thrust: f64) -> f64 {
    ct_at(thrust, OMEGA, R_TIP)
}

pub fn c_q(torque: f64) -> f64 {
    cq_at(torque, OMEGA, R_TIP)
}

pub fn lock_number(i_beta: f64) -> f64 {
    RHO * CL_ALPHA * CHORD * R_TIP.powi(4) / i_beta
}

// ---------------------------------------------------------------------------
// Polar and blade geometry
// ---------------------------------------------------------------------------

pub fn base_polar_params() -> LinearPolarParameters {
    LinearPolarParameters {
        CL0: 0.0,
        CL_alpha_per_rad: CL_ALPHA,
        CD0,
        alpha_stall_deg: ALPHA_STALL_DEG,
    }
}

pub fn polar_for(params: &LinearPolarParameters) -> LinearPolar {
    LinearPolar::from_properties(params)
}

pub fn theory_polar() -> LinearPolar {
    LinearPolar::new(0.0, CL_ALPHA, CD0, ALPHA_STALL_DEG.to_radians())
}

pub fn theory_blade(n_elements: usize, twist_deg: f64, tip_loss: bool) -> BladeGeometry {
    BladeGeometry {
        n_blades: N_BLADES,
        radius_m: R_TIP,
        root_cutout_m: R_ROOT,
        chord_m: CHORD,
        twist_deg,
        n_elements,
        tip_loss,
        r_stations_m: Vec::new(),
        chord_stations_m: Vec::new(),
        twist_stations_deg: Vec::new(),
    }
}

pub fn theory_rotor(n_elements: usize, twist_deg: f64) -> RotorDefinition {
    RotorDefinition {
        blade: theory_blade(n_elements, twist_deg, true),
        airfoil: base_polar_params(),
        control: None,
        pitch_actuation: PitchActuation::DirectMechanical,
        flap: None,
        name: "theory_rotor".to_string(),
        description: String::new(),
    }
}

pub fn theory_rotor_flap(n_elements: usize, i_beta: f64) -> RotorDefinition {
    let mut defn = theory_rotor(n_elements, 0.0);
    defn.flap = Some(FlapProperties {
        I_blade_flap_kgm2: i_beta,
        omega_nr_rad_s: 0.0,
    });
    defn
}

pub fn castles_gray_rotor(n_elements: usize) -> RotorDefinition {
    theory_rotor(n_elements, 0.0)
}

pub fn pca2_rotor(n_elements: usize) -> RotorDefinition {
    RotorDefinition {
        blade: BladeGeometry {
            n_blades: 4,
            radius_m: PCA2_R,
            root_cutout_m: 1.20,
            chord_m: 0.55,
            twist_deg: 0.0,
            n_elements,
            tip_loss: true,
            r_stations_m: Vec::new(),
            chord_stations_m: Vec::new(),
            twist_stations_deg: Vec::new(),
        },
        airfoil: LinearPolarParameters {
            CL0: 0.0,
            CL_alpha_per_rad: 5.73,
            CD0: 0.0098,
            alpha_stall_deg: 15.0,
        },
        control: None,
        pitch_actuation: PitchActuation::DirectMechanical,
        flap: None,
        name: "pca2_wheatley".to_string(),
        description: String::new(),
    }
}

// ---------------------------------------------------------------------------
// VPM rotor creation
// ---------------------------------------------------------------------------

pub fn make_rotor(defn: &RotorDefinition) -> VpmRotor<LinearPolar> {
    let cfg = VpmRotorConfig {
        max_particles: 14000,
        sigma: (1.5 * defn.blade.chord_m) as f32,
        barnes_hut: true,
        bh_theta: 0.5,
        bh_min_particles: 2000,
        ..VpmRotorConfig::default()
    };
    VpmRotor::new(defn, polar_for(&defn.airfoil), ControlGains::default(), cfg)
}

/// Lightweight preset for directional (sign/monotone) checks that do not need
/// quantitative accuracy. Uses VpmRotorConfig::fast_test() -- no Barnes-Hut,
/// 800 particles. Much faster than make_rotor; suitable for checks that just
/// need the right sign, not converged magnitudes.
pub fn make_fast_rotor(defn: &RotorDefinition) -> VpmRotor<LinearPolar> {
    VpmRotor::new(
        defn,
        polar_for(&defn.airfoil),
        ControlGains::default(),
        VpmRotorConfig::fast_test(),
    )
}

pub fn make_pca2_rotor(defn: &RotorDefinition) -> VpmRotor<LinearPolar> {
    let cfg = VpmRotorConfig {
        max_particles: 2000,
        sigma: 0.50,
        relax: 0.35,
        nonlinear_lifting_line: true,
        tip_clustering: true,
        local_core: true,
        barnes_hut: true,
        bh_theta: 0.5,
        bh_min_particles: 200,
        flap_dynamics: false,
        ..VpmRotorConfig::default()
    };
    VpmRotor::new(defn, polar_for(&defn.airfoil), ControlGains::default(), cfg)
}

// ---------------------------------------------------------------------------
// Flight conditions
// ---------------------------------------------------------------------------

pub fn hover_fc(collective_deg: f64) -> FlightCondition {
    hover_fc_omega(collective_deg, OMEGA)
}

pub fn hover_fc_omega(collective_deg: f64, omega: f64) -> FlightCondition {
    FlightCondition {
        collective_rad: collective_deg.to_radians(),
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        v_hub: [0.0, 0.0, 0.0],
        omega_rad_s: omega,
        rho: RHO,
    }
}

pub fn climb_fc(collective_deg: f64, v_climb: f64) -> FlightCondition {
    FlightCondition {
        collective_rad: collective_deg.to_radians(),
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        v_hub: [0.0, 0.0, v_climb],
        omega_rad_s: OMEGA,
        rho: RHO,
    }
}

pub fn forward_fc(collective_deg: f64, mu: f64) -> FlightCondition {
    FlightCondition {
        collective_rad: collective_deg.to_radians(),
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        v_hub: [mu * tip_speed(), 0.0, 0.0],
        omega_rad_s: OMEGA,
        rho: RHO,
    }
}

pub fn wheatley_fc(pitch_deg: f64, mu: f64, alpha_deg: f64, n_rpm: f64) -> FlightCondition {
    let omega = omega_from_rpm(n_rpm);
    let a = alpha_deg.to_radians();
    let v = omega * PCA2_R * mu / a.cos();
    FlightCondition {
        collective_rad: pitch_deg.to_radians(),
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        v_hub: [v * a.cos(), 0.0, -v * a.sin()],
        omega_rad_s: omega,
        rho: RHO,
    }
}

// ---------------------------------------------------------------------------
// Inflow and wake sampling
// ---------------------------------------------------------------------------

pub fn settle(
    rotor: &VpmRotor<LinearPolar>,
    fc: &FlightCondition,
    revs: usize,
) -> (VpmRotorResult, VpmRotorState) {
    let dt = (2.0 * PI / fc.omega_rad_s) / STEPS_PER_REV as f64;
    rotor.march(fc, None, dt, revs * STEPS_PER_REV)
}

pub fn mean_disk_inflow(state: &VpmRotorState) -> f64 {
    let wake = match &state.wake {
        Some(w) if w.len() > 0 => w,
        _ => return 0.0,
    };
    let mut tx = Vec::new();
    let mut ty = Vec::new();
    let mut tz = Vec::new();
    for &rf in &[0.4_f64, 0.5, 0.6, 0.7, 0.8] {
        let r = rf * R_TIP;
        for k in 0..12 {
            let psi = 2.0 * PI * k as f64 / 12.0;
            tx.push((r * psi.cos()) as f32);
            ty.push((r * -psi.sin()) as f32);
            tz.push(0.0_f32);
        }
    }
    let ind = induced_velocities_at_points(wake, &tx, &ty, &tz);
    let mean_vz: f64 = ind.iter().map(|v| v[2] as f64).sum::<f64>() / ind.len() as f64;
    mean_vz / tip_speed()
}

pub fn wake_skew_angle(state: &VpmRotorState, result: &VpmRotorResult) -> f64 {
    let _ = state;
    let c = result.wake_centroid;
    let horiz = (c[0] * c[0] + c[1] * c[1]).sqrt();
    horiz.atan2(c[2].abs().max(1e-9))
}

pub fn sample_flap(
    rotor: &VpmRotor<LinearPolar>,
    fc: &FlightCondition,
    state: &VpmRotorState,
) -> Vec<(f64, f64)> {
    let dt = (2.0 * PI / fc.omega_rad_s) / STEPS_PER_REV as f64;
    let mut s = state.clone();
    let mut out = Vec::with_capacity(STEPS_PER_REV);
    for _ in 0..STEPS_PER_REV {
        let (_r, ns) = rotor.step_one(fc, &s, dt);
        s = ns;
        let beta0 = s.beta.as_ref().map(|b| b[0]).unwrap_or(0.0);
        out.push((s.psi, beta0));
    }
    out
}

pub fn fourier_flap(samples: &[(f64, f64)]) -> (f64, f64, f64) {
    let n = samples.len() as f64;
    let a0 = samples.iter().map(|&(_, b)| b).sum::<f64>() / n;
    let sc: f64 = samples.iter().map(|&(psi, b)| b * psi.cos()).sum::<f64>() / n;
    let ss: f64 = samples.iter().map(|&(psi, b)| b * psi.sin()).sum::<f64>() / n;
    (a0, -2.0 * sc, -2.0 * ss)
}

// ---------------------------------------------------------------------------
// Closed-form theory helpers
// ---------------------------------------------------------------------------

pub fn bemt_hover_lambda(theta_rad: f64) -> f64 {
    let sa = solidity() * CL_ALPHA;
    (sa / 16.0) * ((1.0 + 32.0 * theta_rad / sa).sqrt() - 1.0)
}

pub fn bemt_hover_ct(theta_rad: f64, lambda: f64) -> f64 {
    let sa = solidity() * CL_ALPHA;
    (sa / 2.0) * (theta_rad / 3.0 - lambda / 2.0)
}

pub fn glauert_lambda(ct: f64, mu: f64, alpha_rad: f64) -> f64 {
    let mut lam = mu * alpha_rad.tan() + (ct / 2.0).max(1e-8).sqrt();
    for _ in 0..500 {
        let ln = mu * alpha_rad.tan() + ct / (2.0 * (mu * mu + lam * lam).sqrt());
        if (ln - lam).abs() < 1e-11 {
            lam = ln;
            break;
        }
        lam = ln;
    }
    lam
}

pub fn wheatley_cl_from_thrust(thrust: f64, mu: f64, alpha_deg: f64, n_rpm: f64) -> f64 {
    let omega = omega_from_rpm(n_rpm);
    let a = alpha_deg.to_radians();
    let v = omega * PCA2_R * mu / a.cos();
    let q = 0.5 * RHO * v * v;
    thrust * a.cos() / (q * PI * PCA2_R * PCA2_R)
}

// ---------------------------------------------------------------------------
// BEM model helpers (QS / Pitt-Peters / Oye against Castles-Gray data)
// ---------------------------------------------------------------------------

use dynbem_rs::aero_io::{Mat3, RotorInputs, Vec3};
use dynbem_rs::aero_model::{AeroModel, IntegrationMethod};
use dynbem_rs::oye::OyeBEMModel;
use dynbem_rs::pitt_peters::PittPetersModel;
use dynbem_rs::quasi_static_bem::QuasiStaticBEM;

pub fn bem_hover_inputs(theta_deg: f64, rpm: f64) -> RotorInputs {
    let omega = omega_from_rpm(rpm);
    RotorInputs {
        collective_rad: theta_deg.to_radians(),
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        R_hub: Mat3::eye(),
        v_hub_world: Vec3::zero(),
        wind_world: Vec3::zero(),
        rho_kg_m3: RHO,
        omega_rad_s: omega,
    }
}

pub fn bem_descent_inputs(theta_deg: f64, rpm: f64, v_descent_m_s: f64) -> RotorInputs {
    let omega = omega_from_rpm(rpm);
    RotorInputs {
        collective_rad: theta_deg.to_radians(),
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        R_hub: Mat3::eye(),
        v_hub_world: Vec3::new(0.0, 0.0, v_descent_m_s),
        wind_world: Vec3::zero(),
        rho_kg_m3: RHO,
        omega_rad_s: omega,
    }
}

pub fn bem_ct(thrust_n: f64, rpm: f64) -> f64 {
    let omega = omega_from_rpm(rpm);
    thrust_n / (RHO * PI * R_TIP * R_TIP * (omega * R_TIP).powi(2))
}

pub fn bem_cq(q_spin: f64, rpm: f64) -> f64 {
    let omega = omega_from_rpm(rpm);
    q_spin / (RHO * PI * R_TIP * R_TIP * (omega * R_TIP).powi(2) * R_TIP)
}

/// Run a BEM model for n_steps with ExplicitEuler and return the final CT.
pub fn run_ct<M: AeroModel>(model: &M, inp: &RotorInputs, rpm: f64, n_steps: usize) -> f64 {
    let mut state = model.initial_state();
    let mut ct = 0.0;
    for _ in 0..n_steps {
        let (res, next) = model.step(inp, &state, 0.001, IntegrationMethod::ExplicitEuler);
        state = next;
        ct = bem_ct(-res.F_world[2], rpm);
    }
    ct
}

/// Run a BEM model for n_steps with ExplicitEuler and return the final CQ.
pub fn run_cq<M: AeroModel>(model: &M, inp: &RotorInputs, rpm: f64, n_steps: usize) -> f64 {
    let mut state = model.initial_state();
    let mut cq = 0.0;
    for _ in 0..n_steps {
        let (res, next) = model.step(inp, &state, 0.001, IntegrationMethod::ExplicitEuler);
        state = next;
        cq = bem_cq(res.Q_spin, rpm);
    }
    cq
}

pub fn castles_gray_qs(n_elem: usize) -> QuasiStaticBEM<LinearPolar> {
    let defn = castles_gray_rotor(n_elem);
    let polar = polar_for(&defn.airfoil);
    QuasiStaticBEM::build(defn, n_elem * 3, polar)
}

pub fn castles_gray_pp(n_elem: usize) -> PittPetersModel<LinearPolar> {
    let defn = castles_gray_rotor(n_elem);
    let polar = polar_for(&defn.airfoil);
    PittPetersModel::build(defn, n_elem * 3, polar)
}

pub fn castles_gray_oye(n_elem: usize) -> OyeBEMModel<LinearPolar> {
    let defn = castles_gray_rotor(n_elem);
    let polar = polar_for(&defn.airfoil);
    OyeBEMModel::build(defn, n_elem * 3, polar)
}

/// Parse a CSV (with optional # comment lines) and return owned string records.
pub fn csv_rows(data: &str) -> Vec<csv::StringRecord> {
    let mut rdr = csv::ReaderBuilder::new()
        .comment(Some(b'#'))
        .from_reader(data.as_bytes());
    rdr.records().map(|r| r.expect("CSV parse error")).collect()
}
