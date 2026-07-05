// Exploratory comparison: VPM rigid-flap DOF vs classical closed-form
// flapping theory for a centrally-hinged (nu_beta = 1) rotor.
//
// This is NOT a validation test -- it is a "how close are we?" print. The
// classical theory (Bramwell / Seddon / Prouty, no cyclic pitch, untwisted,
// constant chord, LINEAR aero, UNIFORM inflow, no tip loss / reverse flow)
// gives closed forms for the flap harmonics:
//
//   beta(psi) = a_0 - a_1 cos(psi) - b_1 sin(psi)
//   a_0 = gamma [ theta_0/8 (1 + mu^2) - lambda/6 ]        (coning)
//   a_1 = 2 mu (4/3 theta_0 - lambda) / (1 - mu^2/2)       (longitudinal)
//   b_1 = 4/3 mu a_0 / (1 + mu^2/2)                        (lateral)
//
// with Lock number gamma = rho * a * c * R^4 / I_beta.
//
// The VPM has a free-wake (non-uniform) inflow, a real polar, tip effects and
// (at high mu) reverse flow -- so exact agreement is NOT expected. To keep the
// inflow apples-to-apples, the theory's lambda is taken from the VPM's OWN
// thrust via Glauert momentum. a_0 (coning) and the 1/rev amplitude
// sqrt(a_1^2 + b_1^2) are convention-robust; the individual a_1/b_1 split and
// phase depend on the azimuth origin (ours: psi=0 at +X; classical: psi=0
// downwind), so those are printed for inspection, not asserted.
//
// Run: cargo run --release -p dynbem_rs --example vpm_flapping_vs_theory

use dynbem_rs::cyclic::ControlGains;
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::rotor_definition::{
    BladeGeometry, FlapProperties, LinearPolarParameters, PitchActuation, RotorDefinition,
};
use dynbem_rs::vpm::{FlightCondition, VpmRotor, VpmRotorConfig};
use std::f64::consts::PI;

const N_BLADES: usize = 2;
const R_TIP: f64 = 1.0;
const R_ROOT: f64 = 0.15;
const CHORD: f64 = 0.06;
const CL_ALPHA: f64 = 5.7;
const CD0: f64 = 0.008;
const ALPHA_STALL_DEG: f64 = 15.0;
const RHO: f64 = 1.225;
const OMEGA: f64 = 120.0;
const N_STATIONS: usize = 14;

// Blade flap inertia chosen for a realistic Lock number (~6-8).
const I_BETA: f64 = 0.05;

fn rotor_definition() -> RotorDefinition {
    RotorDefinition {
        blade: BladeGeometry {
            n_blades: N_BLADES,
            radius_m: R_TIP,
            root_cutout_m: R_ROOT,
            chord_m: CHORD,
            twist_deg: 0.0,
            n_elements: N_STATIONS,
            tip_loss: true,
            r_stations_m: Vec::new(),
            chord_stations_m: Vec::new(),
            twist_stations_deg: Vec::new(),
        },
        airfoil: LinearPolarParameters {
            CL0: 0.0,
            CL_alpha_per_rad: CL_ALPHA,
            CD0,
            alpha_stall_deg: ALPHA_STALL_DEG,
        },
        control: None,
        pitch_actuation: PitchActuation::DirectMechanical,
        flap: Some(FlapProperties {
            I_blade_flap_kgm2: I_BETA,
            omega_nr_rad_s: 0.0, // freely hinged -> nu_beta = 1 (matches theory)
        }),
        name: "flap_theory_rotor".to_string(),
        description: "VPM flap DOF vs classical theory".to_string(),
    }
}

fn polar() -> LinearPolar {
    LinearPolar::new(0.0, CL_ALPHA, CD0, ALPHA_STALL_DEG.to_radians())
}

/// Glauert forward-flight induced inflow ratio for a horizontal disk
/// (freestream in-plane, no through-disk component): solve
///   lambda = C_T / (2 sqrt(mu^2 + lambda^2))
fn glauert_lambda(c_t: f64, mu: f64) -> f64 {
    let mut lam = (c_t / 2.0).max(1e-6).sqrt();
    for _ in 0..100 {
        let lam_new = c_t / (2.0 * (mu * mu + lam * lam).sqrt());
        if (lam_new - lam).abs() < 1e-9 {
            lam = lam_new;
            break;
        }
        lam = lam_new;
    }
    lam
}

fn main() {
    let defn = rotor_definition();
    let cfg = VpmRotorConfig {
        max_particles: 3000,
        sigma: 0.12,
        ..VpmRotorConfig::default()
    };
    let rotor = VpmRotor::new(&defn, polar(), ControlGains::default(), cfg);

    let area = PI * R_TIP * R_TIP;
    let tip_speed = OMEGA * R_TIP;
    let gamma = RHO * CL_ALPHA * CHORD * R_TIP.powi(4) / I_BETA;

    let collective_deg = 8.0;
    let theta_0 = (collective_deg as f64).to_radians();

    println!(
        "PCA-style flap comparison: R={R_TIP} m, {N_BLADES} blades, chord={CHORD} m, \
         Omega={OMEGA} rad/s, collective={collective_deg} deg"
    );
    println!(
        "Lock number gamma = {:.2}, I_beta = {} kg*m^2, nu_beta = 1",
        gamma, I_BETA
    );
    println!();
    println!(
        "{:>5} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}",
        "mu", "a0_vpm", "a0_th", "amp_vpm", "amp_th", "a1_vpm", "a1_th", "b1_vpm", "b1_th"
    );

    let steps_per_rev = 48usize;
    let dt = (2.0 * PI / OMEGA) / steps_per_rev as f64;

    for &mu in &[0.10, 0.15, 0.20, 0.25, 0.30] {
        let v_x = mu * tip_speed;
        let fc = FlightCondition {
            collective_rad: theta_0,
            tilt_lon: 0.0,
            tilt_lat: 0.0,
            v_hub: [v_x, 0.0, 0.0], // pure in-plane forward flight (horizontal disk)
            omega_rad_s: OMEGA,
            rho: RHO,
        };

        // Settle the wake + flap DOF (~10 revs).
        let n_settle = 10 * steps_per_rev;
        let (res, mut state) = rotor.march(&fc, None, dt, n_settle);
        let c_t = res.thrust / (RHO * area * tip_speed * tip_speed);

        // Sample beta_0(psi) over one revolution.
        let mut sum = 0.0;
        let mut sum_c = 0.0;
        let mut sum_s = 0.0;
        let m = steps_per_rev;
        for _ in 0..m {
            let (_r, s) = rotor.step_one(&fc, &state, dt);
            state = s;
            let psi = state.psi;
            let beta0 = state.beta.as_ref().map(|b| b[0]).unwrap_or(0.0);
            sum += beta0;
            sum_c += beta0 * psi.cos();
            sum_s += beta0 * psi.sin();
        }
        let a0_vpm = sum / m as f64;
        // beta ~ a0 + A1c cos + A1s sin ; classical a1 = -A1c, b1 = -A1s.
        let a1c = 2.0 * sum_c / m as f64;
        let a1s = 2.0 * sum_s / m as f64;
        let a1_vpm = -a1c;
        let b1_vpm = -a1s;
        let amp_vpm = (a1_vpm * a1_vpm + b1_vpm * b1_vpm).sqrt();

        // Theory, with lambda from the VPM's own thrust (Glauert).
        let lambda = glauert_lambda(c_t, mu);
        let a0_th = gamma * (theta_0 / 8.0 * (1.0 + mu * mu) - lambda / 6.0);
        let a1_th = 2.0 * mu * (4.0 / 3.0 * theta_0 - lambda) / (1.0 - 0.5 * mu * mu);
        let b1_th = (4.0 / 3.0) * mu * a0_th / (1.0 + 0.5 * mu * mu);
        let amp_th = (a1_th * a1_th + b1_th * b1_th).sqrt();

        println!(
            "{:>5.2} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4}",
            mu,
            a0_vpm.to_degrees(),
            a0_th.to_degrees(),
            amp_vpm.to_degrees(),
            amp_th.to_degrees(),
            a1_vpm.to_degrees(),
            a1_th.to_degrees(),
            b1_vpm.to_degrees(),
            b1_th.to_degrees(),
        );
    }
    println!();
    println!("(angles in degrees; a0=coning, amp=sqrt(a1^2+b1^2) 1/rev amplitude)");
    println!("a0 and amp are convention-robust; a1/b1 split depends on azimuth origin.");
}
