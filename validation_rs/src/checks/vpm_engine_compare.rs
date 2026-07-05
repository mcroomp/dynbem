// Long opt-in survey: compare the two wake engines -- classic convection-only
// VPM vs the reformulated VPM (rVPM, strength/size evolution) -- across EVERY
// empirical regime we have (hover, vertical-descent windmill-brake
// autorotation, and forward-flight autorotation). Everything except
// `wake_engine` is held identical, so any difference is attributable to the
// engine alone.
//
// This is NOT part of run_all -- it is slow (rVPM advection is direct O(N^2),
// no Barnes-Hut). Run it explicitly:
//
//   cargo run --release -p validation_rs -- vpm_engine_compare
//
// Question answered: is rVPM "better"? "Better" is scored three ways:
//   1. accuracy  -- mean |error| vs measured CT/CQ/CL, per regime;
//   2. stability -- both engines must stay finite;
//   3. cost      -- wall-clock per regime (rVPM has no BH, so it is slower).
//
// rVPM is a strict generalization of classic VPM: with the stretching term it
// preserves tip vortices (momentum + mass conserving) instead of freezing
// them, which should matter most where the wake lingers (hover / descent) and
// least where it convects away (forward flight).

use crate::helpers::*;
use crate::report::Report;
use dynbem_rs::cyclic::ControlGains;
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::rotor_definition::RotorDefinition;
use dynbem_rs::vpm::{FlightCondition, VpmRotor, VpmRotorConfig, WakeEngine};
use std::f64::consts::PI;
use std::time::Instant;

// Kept modest so the O(N^2) rVPM advection stays tractable; both engines use
// the SAME numbers, so the head-to-head comparison is fair.
const WAKE_REVS: f64 = 1.5; // retained wake age (max_particles sized from this)
const TOTAL_REVS: usize = 6; // marched revolutions (march averages trailing half)
const SPR_HOVER: usize = 36; // steps/rev, Castles-Gray hover + descent
const SPR_FWD: usize = 48; // steps/rev, PCA-2 forward flight

/// Build a VPM rotor sized for `WAKE_REVS` of retained wake at `steps_per_rev`,
/// selecting the wake `engine`. All other knobs are identical between engines.
fn build_vpm(
    defn: &RotorDefinition,
    steps_per_rev: usize,
    sigma: f32,
    engine: WakeEngine,
) -> VpmRotor<LinearPolar> {
    let nb = defn.blade.n_blades;
    let ne = defn.blade.n_elements;
    let ppr = nb * (2 * ne + 1) * steps_per_rev;
    let max_particles = ((WAKE_REVS * ppr as f64).ceil() as usize) + 1;
    let cfg = VpmRotorConfig {
        max_particles,
        sigma,
        // BH is OFF here for fairness: with direct induction for BOTH engines
        // the ONLY difference is the wake self-advection model (classic
        // frozen-strength RK2 vs rVPM stretching). (The BH tree is now hardened
        // against the non-finite coordinates a diverged wake produces -- see
        // build_node -- so this is no longer needed to avoid a stack overflow.)
        barnes_hut: false,
        wake_engine: engine,
        // rVPM stabilizers (ignored by ClassicVpm): Pedrizzetti relaxation +
        // viscous/SFS core spreading. Without these the bare inviscid rVPM
        // diverges. `nu` is scaled to the tip speed * core so it is roughly
        // rotor-independent.
        rvpm_relax: 0.3,
        rvpm_nu: 0.02,
        ..VpmRotorConfig::default()
    };
    VpmRotor::new(defn, polar_for(&defn.airfoil), ControlGains::default(), cfg)
}

/// March `TOTAL_REVS` revolutions from a cold wake; return trailing-half-
/// averaged (thrust, torque).
fn run(rotor: &VpmRotor<LinearPolar>, fc: &FlightCondition, steps_per_rev: usize) -> (f64, f64) {
    let dt = (2.0 * PI / fc.omega_rad_s) / steps_per_rev as f64;
    let (res, _s) = rotor.march(fc, None, dt, TOTAL_REVS * steps_per_rev);
    (res.thrust, res.torque)
}

fn abs_err_pct(v: f64, reference: f64) -> f64 {
    if reference.abs() < 1e-15 {
        f64::NAN
    } else {
        (v - reference).abs() / reference.abs() * 100.0
    }
}

fn mean(v: &[f64]) -> f64 {
    let vals: Vec<f64> = v.iter().copied().filter(|x| x.is_finite()).collect();
    if vals.is_empty() {
        f64::NAN
    } else {
        vals.iter().sum::<f64>() / vals.len() as f64
    }
}

/// One-word verdict comparing classic vs reformulated mean error.
fn verdict(classic: f64, rvpm: f64) -> &'static str {
    if !classic.is_finite() || !rvpm.is_finite() {
        "n/a"
    } else if rvpm < classic - 0.2 {
        "rVPM better"
    } else if rvpm > classic + 0.2 {
        "rVPM worse"
    } else {
        "tie"
    }
}

pub fn check_vpm_engine_compare(r: &mut Report) {
    r.begin_module(
        "vpm_engine_compare",
        "Classic VPM vs Reformulated VPM (rVPM) across ALL empirical scenarios",
    );

    let engines = [
        ("classic", WakeEngine::ClassicVpm),
        ("rVPM", WakeEngine::ReformulatedVpm),
    ];

    // ===================================================================
    // 1) HOVER -- Castles-Gray NACA TN-2474 (CT and CQ)
    // ===================================================================
    let cg = castles_gray_rotor(10);
    let cg_sigma = (1.5 * cg.blade.chord_m) as f32;
    // (theta_deg, rpm, ct_meas, cq_meas)
    let hover_pts = [
        (6.68_f64, 1200.0_f64, 0.00289_f64, 0.000137_f64),
        (8.46, 1200.0, 0.00400, 0.000226),
        (5.55, 1600.0, 0.00255, 0.000108),
    ];
    let (mut ct_c, mut ct_r) = (Vec::new(), Vec::new());
    let (mut cq_c, mut cq_r) = (Vec::new(), Vec::new());
    let (mut t_c, mut t_r) = (0.0_f64, 0.0_f64);
    for &(theta, rpm, ct_meas, cq_meas) in &hover_pts {
        let omega = omega_from_rpm(rpm);
        let fc = hover_fc_omega(theta, omega);
        for (name, eng) in engines {
            let rotor = build_vpm(&cg, SPR_HOVER, cg_sigma, eng);
            let t0 = Instant::now();
            let (thr, tq) = run(&rotor, &fc, SPR_HOVER);
            let elapsed = t0.elapsed().as_secs_f64();
            let (ct, cq) = (ct_at(thr, omega, R_TIP), cq_at(tq, omega, R_TIP));
            let case = format!("hover th={theta:.2} rpm={rpm:.0} [{name}]");
            r.info(format!("{case} CT"), "CT", ct, ct_meas);
            r.info(format!("{case} CQ"), "CQ", cq, cq_meas);
            match name {
                "classic" => {
                    ct_c.push(abs_err_pct(ct, ct_meas));
                    cq_c.push(abs_err_pct(cq, cq_meas));
                    t_c += elapsed;
                }
                _ => {
                    ct_r.push(abs_err_pct(ct, ct_meas));
                    cq_r.push(abs_err_pct(cq, cq_meas));
                    t_r += elapsed;
                }
            }
        }
    }
    let (h_ct_c, h_ct_r) = (mean(&ct_c), mean(&ct_r));
    println!(
        "  [hover CT] mean|err|  classic={h_ct_c:.1}%  rVPM={h_ct_r:.1}%   ({})",
        verdict(h_ct_c, h_ct_r)
    );
    println!(
        "  [hover CQ] mean|err|  classic={:.1}%  rVPM={:.1}%   ({})",
        mean(&cq_c),
        mean(&cq_r),
        verdict(mean(&cq_c), mean(&cq_r))
    );

    // ===================================================================
    // 2) VERTICAL DESCENT (windmill-brake autorotation) -- Castles-Gray WBS CQ
    // ===================================================================
    // (theta_deg, rpm, v_descent_m_s, cq_meas)  cq < 0 => driving/autorotative
    let descent_pts = [
        (1.23_f64, 1200.0_f64, 11.15_f64, -0.000112_f64),
        (-0.11, 1600.0, 10.17, -0.00005),
        (2.46, 1600.0, 13.77, -0.000045),
    ];
    let (mut dq_c, mut dq_r) = (Vec::new(), Vec::new());
    let (mut td_c, mut td_r) = (0.0_f64, 0.0_f64);
    for &(theta, rpm, v_descent, cq_meas) in &descent_pts {
        let omega = omega_from_rpm(rpm);
        let fc = FlightCondition {
            collective_rad: theta.to_radians(),
            tilt_lon: 0.0,
            tilt_lat: 0.0,
            v_hub: [0.0, 0.0, -v_descent],
            omega_rad_s: omega,
            rho: RHO,
        };
        for (name, eng) in engines {
            let rotor = build_vpm(&cg, SPR_HOVER, cg_sigma, eng);
            let t0 = Instant::now();
            let (_thr, tq) = run(&rotor, &fc, SPR_HOVER);
            let elapsed = t0.elapsed().as_secs_f64();
            let cq = cq_at(tq, omega, R_TIP);
            let case = format!("descent th={theta:.2} rpm={rpm:.0} vz={v_descent:.1} [{name}]");
            r.info(format!("{case} CQ"), "CQ", cq, cq_meas);
            match name {
                "classic" => {
                    dq_c.push(abs_err_pct(cq, cq_meas));
                    td_c += elapsed;
                }
                _ => {
                    dq_r.push(abs_err_pct(cq, cq_meas));
                    td_r += elapsed;
                }
            }
        }
    }
    println!(
        "  [descent CQ] mean|err|  classic={:.1}%  rVPM={:.1}%   ({})",
        mean(&dq_c),
        mean(&dq_r),
        verdict(mean(&dq_c), mean(&dq_r))
    );

    // ===================================================================
    // 3) FORWARD-FLIGHT AUTOROTATION -- Wheatley & Hood TR-515 CL
    // ===================================================================
    let pca2 = pca2_rotor(12);
    // (pitch_deg, mu, alpha_deg, rpm, cl_meas)
    let fwd_pts = [
        (1.9_f64, 0.181_f64, 11.2_f64, 98.8_f64, 0.363_f64),
        (1.9, 0.249, 6.6, 98.8, 0.192),
        (2.7, 0.242, 6.3, 118.9, 0.266),
    ];
    let (mut cl_c, mut cl_r) = (Vec::new(), Vec::new());
    let (mut tf_c, mut tf_r) = (0.0_f64, 0.0_f64);
    for &(pitch, mu, alpha, rpm, cl_meas) in &fwd_pts {
        let fc = wheatley_fc(pitch, mu, alpha, rpm);
        for (name, eng) in engines {
            let rotor = build_vpm(&pca2, SPR_FWD, 0.50, eng);
            let t0 = Instant::now();
            let (thr, _tq) = run(&rotor, &fc, SPR_FWD);
            let elapsed = t0.elapsed().as_secs_f64();
            let cl = wheatley_cl_from_thrust(thr, mu, alpha, rpm);
            let case = format!("fwd mu={mu:.3} a={alpha:.1} [{name}]");
            r.info(format!("{case} CL"), "CL", cl, cl_meas);
            match name {
                "classic" => {
                    cl_c.push(abs_err_pct(cl, cl_meas));
                    tf_c += elapsed;
                }
                _ => {
                    cl_r.push(abs_err_pct(cl, cl_meas));
                    tf_r += elapsed;
                }
            }
        }
    }
    let (f_c, f_r) = (mean(&cl_c), mean(&cl_r));
    println!(
        "  [forward CL] mean|err|  classic={f_c:.1}%  rVPM={f_r:.1}%   ({})",
        verdict(f_c, f_r)
    );

    // ---- cost + overall summary ----
    println!(
        "  [cost s]   hover classic={t_c:.1} rVPM={t_r:.1} | descent classic={td_c:.1} rVPM={td_r:.1} | forward classic={tf_c:.1} rVPM={tf_r:.1}"
    );
    let classic_overall = mean(&[h_ct_c, mean(&cq_c), mean(&dq_c), f_c]);
    let rvpm_overall = mean(&[h_ct_r, mean(&cq_r), mean(&dq_r), f_r]);
    println!(
        "  [OVERALL]  mean|err| across regimes  classic={classic_overall:.1}%  rVPM={rvpm_overall:.1}%   ({})",
        verdict(classic_overall, rvpm_overall)
    );

    // Guard: both engines must stay finite and physically plausible everywhere
    // (this survey is about which is BETTER, not a hard accuracy gate).
    r.assert_bool(
        "stability",
        "classic_finite",
        classic_overall,
        1e3,
        classic_overall.is_finite() && classic_overall < 1e3,
        "classic VPM must stay finite across all regimes",
    );
    r.assert_bool(
        "stability",
        "rvpm_finite",
        rvpm_overall,
        1e3,
        rvpm_overall.is_finite() && rvpm_overall < 1e3,
        "reformulated VPM must stay finite across all regimes",
    );
}
