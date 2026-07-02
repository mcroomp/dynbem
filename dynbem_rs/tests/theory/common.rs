// Shared harness for the VPM-vs-standard-rotor-theory tests.
//
// One clean "theory rotor" (constant chord, linear polar, optional twist) so
// the classical closed forms actually apply, plus reusable extractors that
// turn a settled VPM run into the observables the theory is stated in:
// C_T, C_Q/C_P, mean disk inflow, figure of merit, flap harmonics, wake skew.
#![allow(dead_code)]

use dynbem_rs::cyclic::ControlGains;
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::rotor_definition::{
    BladeGeometry, FlapProperties, LinearPolarParameters, PitchActuation, RotorDefinition,
};
use dynbem_rs::vpm::induced_at_points;
use dynbem_rs::vpm_rotor::{
    FlightCondition, VpmRotor, VpmRotorConfig, VpmRotorResult, VpmRotorState,
};
use std::f64::consts::PI;

// ---- theory rotor geometry --------------------------------------------------
//
// All theory tests use the realistic Castles & Gray NACA TN-2474 6-ft rotor
// (constant chord, so BEMT closed forms still apply) rather than a synthetic
// low-drag rotor -- a real CD0 / solidity keeps every derived quantity
// (figure of merit, power, flapping) physical. Hover is validated against
// this rotor's measured C_T/C_Q directly.

pub const R_TIP: f64 = 0.914;
pub const R_ROOT: f64 = 0.155;
pub const CHORD: f64 = 0.0479;
pub const N_BLADES: usize = 3;
pub const CL_ALPHA: f64 = 5.90; // lift-curve slope a [1/rad]
pub const CD0: f64 = 0.01046;
pub const ALPHA_STALL_DEG: f64 = 15.5;
pub const RHO: f64 = 1.225;
pub const OMEGA: f64 = 125.6637; // 1200 rpm (matches the CG hover test speed)
pub const STEPS_PER_REV: usize = 24;

/// Rotor solidity sigma = N_b c / (pi R).
pub fn solidity() -> f64 {
    N_BLADES as f64 * CHORD / (PI * R_TIP)
}

/// Tip speed Omega*R.
pub fn tip_speed() -> f64 {
    OMEGA * R_TIP
}

/// Lock number gamma = rho a c R^4 / I_beta (needs a flap inertia).
pub fn lock_number(i_beta: f64) -> f64 {
    RHO * CL_ALPHA * CHORD * R_TIP.powi(4) / i_beta
}

pub fn theory_rotor(n_elements: usize, twist_deg: f64) -> RotorDefinition {
    RotorDefinition {
        blade: BladeGeometry {
            n_blades: N_BLADES,
            radius_m: R_TIP,
            root_cutout_m: R_ROOT,
            chord_m: CHORD,
            twist_deg,
            n_elements,
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
        flap: None,
        name: "theory_rotor".to_string(),
        description: String::new(),
    }
}

/// Theory rotor with a freely-hinged (nu_beta = 1) flap DOF for the flapping
/// module.
pub fn theory_rotor_flap(n_elements: usize, i_beta: f64) -> RotorDefinition {
    let mut defn = theory_rotor(n_elements, 0.0);
    defn.flap = Some(FlapProperties {
        I_blade_flap_kgm2: i_beta,
        omega_nr_rad_s: 0.0,
    });
    defn
}

pub fn polar() -> LinearPolar {
    LinearPolar::new(0.0, CL_ALPHA, CD0, ALPHA_STALL_DEG.to_radians())
}

/// Castles & Gray (1951) NACA TN-2474 6-ft rotor -- the realistic rotor all
/// theory tests share (identical to `theory_rotor` with zero twist). Named
/// separately for the hover-vs-measured module that references the dataset.
pub fn castles_gray_rotor(n_elements: usize) -> RotorDefinition {
    theory_rotor(n_elements, 0.0)
}

/// Linear polar matching a given rotor definition's airfoil parameters.
pub fn polar_for(defn: &RotorDefinition) -> LinearPolar {
    LinearPolar::from_properties(&defn.airfoil)
}

pub fn omega_from_rpm(rpm: f64) -> f64 {
    rpm * PI / 30.0
}

/// C_T at an arbitrary tip speed: T / (rho A (Omega R)^2).
pub fn ct_at(thrust: f64, omega: f64, r: f64) -> f64 {
    thrust / (RHO * PI * r * r * (omega * r).powi(2))
}

/// C_Q = C_P at an arbitrary tip speed: Q / (rho A (Omega R)^2 R).
pub fn cq_at(torque: f64, omega: f64, r: f64) -> f64 {
    torque / (RHO * PI * r * r * (omega * r).powi(2) * r)
}

pub fn make_rotor(defn: &RotorDefinition) -> VpmRotor<LinearPolar> {
    // Large particle budget + Barnes-Hut so the hover/forward wake develops
    // enough age for the induced inflow to converge. Core scaled to the chord.
    let cfg = VpmRotorConfig {
        max_particles: 14000,
        sigma: (1.5 * defn.blade.chord_m) as f32,
        barnes_hut: true,
        bh_theta: 0.5,
        bh_min_particles: 2000,
        ..VpmRotorConfig::default()
    };
    VpmRotor::new(defn, polar_for(defn), ControlGains::default(), cfg)
}

// ---- flight conditions -----------------------------------------------------

pub fn hover_fc(collective_deg: f64) -> FlightCondition {
    hover_fc_omega(collective_deg, OMEGA)
}

/// Hover at an arbitrary rotor speed (for the measured-rpm datasets).
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

/// Axial climb: `v_climb > 0` = climbing, i.e. air passes down through the
/// disk faster (+Z through-disk in NED, on top of the induced flow).
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

/// Edgewise forward flight, horizontal disk (no through-disk freestream):
/// `v_hub = [mu * Omega R, 0, 0]`. Advance ratio `mu`.
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

// ---- run + observables -----------------------------------------------------

pub fn dt() -> f64 {
    (2.0 * PI / OMEGA) / STEPS_PER_REV as f64
}

/// Settle the wake + DOFs for `revs` revolutions and return the cycle-averaged
/// loads plus the final state. Uses the flight condition's own rotor speed.
pub fn settle(
    rotor: &VpmRotor<LinearPolar>,
    fc: &FlightCondition,
    revs: usize,
) -> (VpmRotorResult, VpmRotorState) {
    let dt = (2.0 * PI / fc.omega_rad_s) / STEPS_PER_REV as f64;
    rotor.march(fc, None, dt, revs * STEPS_PER_REV)
}

/// Thrust coefficient C_T = T / (rho A (Omega R)^2).
pub fn c_t(thrust: f64) -> f64 {
    thrust / (RHO * PI * R_TIP * R_TIP * tip_speed().powi(2))
}

/// Torque / power coefficient C_Q = C_P = Q / (rho A (Omega R)^2 R).
pub fn c_q(torque: f64) -> f64 {
    torque / (RHO * PI * R_TIP * R_TIP * tip_speed().powi(2) * R_TIP)
}

/// Mean through-disk induced inflow ratio lambda = mean(v_z) / (Omega R),
/// sampled from the particle wake on a ring of disk points (z = 0). This is
/// the direct momentum-theory observable (thrust <-> inflow).
pub fn mean_disk_inflow(state: &VpmRotorState) -> f64 {
    let wake = match &state.wake {
        Some(w) if w.len() > 0 => w,
        _ => return 0.0,
    };
    let mut tx = Vec::new();
    let mut ty = Vec::new();
    let mut tz = Vec::new();
    // r_hat(psi) = [cos psi, -sin psi, 0] (project convention).
    for &rf in &[0.4_f64, 0.5, 0.6, 0.7, 0.8] {
        let r = rf * R_TIP;
        for k in 0..12 {
            let psi = 2.0 * PI * k as f64 / 12.0;
            tx.push((r * psi.cos()) as f32);
            ty.push((r * -psi.sin()) as f32);
            tz.push(0.0);
        }
    }
    let ind = induced_at_points(wake, &tx, &ty, &tz);
    let mean_vz: f64 = ind.iter().map(|v| v[2] as f64).sum::<f64>() / ind.len() as f64;
    mean_vz / tip_speed()
}

/// Wake-skew angle chi = atan2(horizontal centroid offset, |vertical offset|)
/// -- the angle the wake leans from the disk normal. Compare to the Glauert
/// chi = atan2(mu, |lambda|).
pub fn wake_skew_angle(state: &VpmRotorState, result: &VpmRotorResult) -> f64 {
    let c = result.wake_centroid;
    let horiz = (c[0] * c[0] + c[1] * c[1]).sqrt();
    horiz.atan2(c[2].abs().max(1e-9))
}

// ---- closed-form helpers ---------------------------------------------------

/// Glauert forward-flight inflow: solve lambda = mu tan(alpha) +
/// C_T / (2 sqrt(mu^2 + lambda^2)) by fixed-point iteration.
pub fn glauert_lambda(c_t: f64, mu: f64, alpha_rad: f64) -> f64 {
    let mut lam = mu * alpha_rad.tan() + (c_t / 2.0).max(1e-8).sqrt();
    for _ in 0..500 {
        let ln = mu * alpha_rad.tan() + c_t / (2.0 * (mu * mu + lam * lam).sqrt());
        if (ln - lam).abs() < 1e-11 {
            lam = ln;
            break;
        }
        lam = ln;
    }
    lam
}

/// Combined-BEMT hover inflow (uniform, no tip loss):
/// lambda = (sigma a / 16) (sqrt(1 + 32 theta / (sigma a)) - 1).
pub fn bemt_hover_lambda(theta_rad: f64) -> f64 {
    let sa = solidity() * CL_ALPHA;
    (sa / 16.0) * ((1.0 + 32.0 * theta_rad / sa).sqrt() - 1.0)
}

/// Blade-element hover thrust with uniform inflow:
/// C_T = (sigma a / 2)(theta_0/3 - lambda/2).
pub fn bemt_hover_ct(theta_rad: f64, lambda: f64) -> f64 {
    let sa = solidity() * CL_ALPHA;
    (sa / 2.0) * (theta_rad / 3.0 - lambda / 2.0)
}

/// Fourier-fit beta(psi) samples to (a0, a1, b1) for
/// beta = a0 - a1 cos(psi) - b1 sin(psi).
pub fn fourier_flap(samples: &[(f64, f64)]) -> (f64, f64, f64) {
    let n = samples.len() as f64;
    let a0 = samples.iter().map(|(_, b)| b).sum::<f64>() / n;
    let sc = samples.iter().map(|(p, b)| b * p.cos()).sum::<f64>() / n;
    let ss = samples.iter().map(|(p, b)| b * p.sin()).sum::<f64>() / n;
    (a0, -2.0 * sc, -2.0 * ss)
}

/// Sample blade-0 flap angle over one revolution at a settled state.
pub fn sample_flap(
    rotor: &VpmRotor<LinearPolar>,
    fc: &FlightCondition,
    state: &VpmRotorState,
) -> Vec<(f64, f64)> {
    let mut s = state.clone();
    let mut out = Vec::with_capacity(STEPS_PER_REV);
    for _ in 0..STEPS_PER_REV {
        let (_r, ns) = rotor.step_one(fc, &s, dt());
        s = ns;
        let beta0 = s.beta.as_ref().map(|b| b[0]).unwrap_or(0.0);
        out.push((s.psi, beta0));
    }
    out
}
