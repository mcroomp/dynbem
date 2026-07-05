// Long opt-in survey: run the VPM with the "ideal" wake-aging config across
// EVERY empirical dataset we have (hover, vertical-descent windmill-brake
// autorotation, and forward-flight autorotation) and compare aging ON vs OFF.
//
// This is NOT part of run_all -- it is slow (many VPM operating points, each
// marched several revolutions). Run it explicitly:
//
//   cargo run --release -p validation_rs -- vpm_aging_survey
//
// Purpose: earlier work showed that adding wake aging (strength fade / core
// spreading) removes the "more wake = worse" degradation and lets a longer
// wake improve accuracy on the Wheatley forward-flight case. This survey asks
// whether that same single aging config also helps (or at least does not hurt)
// on the hover and descent datasets -- i.e. whether it generalizes.
//
// The "ideal" config uses the dimensionless strength-fade knob
// (`strength_decay_tau_rev`), which is RPM- and rotor-independent, so one value
// applies to every scenario.

use crate::helpers::*;
use crate::report::Report;
use dynbem_rs::cyclic::ControlGains;
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::rotor_definition::RotorDefinition;
use dynbem_rs::vpm::{FlightCondition, VpmRotor, VpmRotorConfig};
use std::f64::consts::PI;

// ---- ideal aging config ----
const TAU_REV: f64 = 1.0; // strength fade: decay to 1/e over 1 revolution
const WAKE_REVS: f64 = 2.5; // retained wake age (max_particles sized from this)
const TOTAL_REVS: usize = 8; // marched revolutions (march averages trailing half)
const SPR_HOVER: usize = 36; // steps/rev, Castles-Gray hover + descent
const SPR_FWD: usize = 48; // steps/rev, PCA-2 forward flight

/// Build a VPM rotor sized for `wake_revs` of retained wake at `steps_per_rev`,
/// with the strength-fade aging knob set to `tau_rev` (0 = aging off).
fn build_vpm(
    defn: &RotorDefinition,
    steps_per_rev: usize,
    wake_revs: f64,
    sigma: f32,
    tau_rev: f64,
) -> VpmRotor<LinearPolar> {
    let nb = defn.blade.n_blades;
    let ne = defn.blade.n_elements;
    let ppr = nb * (2 * ne + 1) * steps_per_rev;
    let max_particles = ((wake_revs * ppr as f64).ceil() as usize) + 1;
    let cfg = VpmRotorConfig {
        max_particles,
        sigma,
        barnes_hut: true,
        bh_theta: 0.5,
        bh_min_particles: 512,
        strength_decay_tau_rev: tau_rev,
        ..VpmRotorConfig::default()
    };
    VpmRotor::new(defn, polar_for(&defn.airfoil), ControlGains::default(), cfg)
}

/// March `total_revs` revolutions from a cold wake; return trailing-half-
/// averaged (thrust, torque).
fn run(
    rotor: &VpmRotor<LinearPolar>,
    fc: &FlightCondition,
    steps_per_rev: usize,
    total_revs: usize,
) -> (f64, f64) {
    let dt = (2.0 * PI / fc.omega_rad_s) / steps_per_rev as f64;
    let (res, _s) = rotor.march(fc, None, dt, total_revs * steps_per_rev);
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

pub fn check_vpm_aging_survey(r: &mut Report) {
    r.begin_module(
        "vpm_aging_survey",
        "VPM ideal aging config (strength-fade tau=1 rev) across ALL empirical scenarios",
    );

    // ===================================================================
    // 1) HOVER -- Castles-Gray NACA TN-2474 (CT and CQ), 1200 & 1600 rpm
    // ===================================================================
    let cg = castles_gray_rotor(10);
    let cg_sigma = (1.5 * cg.blade.chord_m) as f32;
    // (theta_deg, rpm, ct_meas, cq_meas)
    let hover_pts = [
        (4.91_f64, 1200.0_f64, 0.00168_f64, 0.00007_f64),
        (6.68, 1200.0, 0.00289, 0.000137),
        (8.46, 1200.0, 0.00400, 0.000226),
        (10.29, 1200.0, 0.00488, 0.000342),
        (3.96, 1600.0, 0.00160, 0.000053),
        (5.55, 1600.0, 0.00255, 0.000108),
        (7.18, 1600.0, 0.00346, 0.000194),
    ];
    let (mut ct_off, mut ct_on) = (Vec::new(), Vec::new());
    let (mut cq_off, mut cq_on) = (Vec::new(), Vec::new());
    for &(theta, rpm, ct_meas, cq_meas) in &hover_pts {
        let omega = omega_from_rpm(rpm);
        let fc = hover_fc_omega(theta, omega);

        let rotor_off = build_vpm(&cg, SPR_HOVER, WAKE_REVS, cg_sigma, 0.0);
        let (t0, q0) = run(&rotor_off, &fc, SPR_HOVER, TOTAL_REVS);
        let (ct0, cq0) = (ct_at(t0, omega, R_TIP), cq_at(q0, omega, R_TIP));

        let rotor_on = build_vpm(&cg, SPR_HOVER, WAKE_REVS, cg_sigma, TAU_REV);
        let (t1, q1) = run(&rotor_on, &fc, SPR_HOVER, TOTAL_REVS);
        let (ct1, cq1) = (ct_at(t1, omega, R_TIP), cq_at(q1, omega, R_TIP));

        let case = format!("hover th={theta:.2} rpm={rpm:.0}");
        r.info(format!("{case} CT_off"), "CT", ct0, ct_meas);
        r.info(format!("{case} CT_on"), "CT", ct1, ct_meas);
        r.info(format!("{case} CQ_off"), "CQ", cq0, cq_meas);
        r.info(format!("{case} CQ_on"), "CQ", cq1, cq_meas);

        ct_off.push(abs_err_pct(ct0, ct_meas));
        ct_on.push(abs_err_pct(ct1, ct_meas));
        cq_off.push(abs_err_pct(cq0, cq_meas));
        cq_on.push(abs_err_pct(cq1, cq_meas));
    }
    let (h_ct_off, h_ct_on) = (mean(&ct_off), mean(&ct_on));
    println!(
        "  [hover CT] mean|err|  OFF={h_ct_off:.1}%  ON={h_ct_on:.1}%   (aging {})",
        verdict(h_ct_off, h_ct_on)
    );
    println!(
        "  [hover CQ] mean|err|  OFF={:.1}%  ON={:.1}%   (aging {})",
        mean(&cq_off),
        mean(&cq_on),
        verdict(mean(&cq_off), mean(&cq_on))
    );
    // Regime finding: aging (a convecting-wake model) should NOT help hover,
    // where the wake lingers under the disk and its persistence IS the inflow.
    // Guard the conclusion: aging-off hover CT is good AND aging does not
    // improve it.
    r.assert_bool(
        "hover_regime",
        "aging_off_ct_accurate",
        h_ct_off,
        15.0,
        h_ct_off < 15.0,
        "aging-off hover CT mean|err| should stay accurate (<15%)",
    );
    r.assert_bool(
        "hover_regime",
        "aging_does_not_help_hover",
        h_ct_on,
        h_ct_off,
        h_ct_on >= h_ct_off - 0.2,
        "strength-fade should not improve hover CT (wake must persist in hover)",
    );

    // ===================================================================
    // 2) VERTICAL DESCENT (windmill-brake autorotation) -- Castles-Gray WBS CQ
    // ===================================================================
    // (theta_deg, rpm, v_descent_m_s, cq_meas)  cq < 0 => driving/autorotative
    let descent_pts = [
        (1.23_f64, 1200.0_f64, 11.15_f64, -0.000112_f64),
        (1.18, 1200.0, 11.15, -0.000116),
        (-0.11, 1600.0, 10.17, -0.00005),
        (-1.32, 1600.0, 11.23, -0.000072),
        (-1.66, 1600.0, 11.91, -0.000084),
        (2.46, 1600.0, 13.77, -0.000045),
    ];
    let (mut dq_off, mut dq_on) = (Vec::new(), Vec::new());
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
        let rotor_off = build_vpm(&cg, SPR_HOVER, WAKE_REVS, cg_sigma, 0.0);
        let (_t0, q0) = run(&rotor_off, &fc, SPR_HOVER, TOTAL_REVS);
        let cq0 = cq_at(q0, omega, R_TIP);

        let rotor_on = build_vpm(&cg, SPR_HOVER, WAKE_REVS, cg_sigma, TAU_REV);
        let (_t1, q1) = run(&rotor_on, &fc, SPR_HOVER, TOTAL_REVS);
        let cq1 = cq_at(q1, omega, R_TIP);

        let case = format!("descent th={theta:.2} rpm={rpm:.0} vz={v_descent:.1}");
        r.info(format!("{case} CQ_off"), "CQ", cq0, cq_meas);
        r.info(format!("{case} CQ_on"), "CQ", cq1, cq_meas);
        dq_off.push(abs_err_pct(cq0, cq_meas));
        dq_on.push(abs_err_pct(cq1, cq_meas));
    }
    println!(
        "  [descent CQ] mean|err|  OFF={:.1}%  ON={:.1}%   (aging {})",
        mean(&dq_off),
        mean(&dq_on),
        verdict(mean(&dq_off), mean(&dq_on))
    );

    // ===================================================================
    // 3) FORWARD-FLIGHT AUTOROTATION -- Wheatley & Hood TR-515 CL
    // ===================================================================
    let pca2 = pca2_rotor(12);
    // (pitch_deg, mu, alpha_deg, rpm, cl_meas)
    let fwd_pts = [
        (1.9_f64, 0.181_f64, 11.2_f64, 98.8_f64, 0.363_f64),
        (1.9, 0.249, 6.6, 98.8, 0.192),
        (1.9, 0.315, 4.3, 98.8, 0.116),
        (2.7, 0.242, 6.3, 118.9, 0.266),
    ];
    let (mut cl_off, mut cl_on) = (Vec::new(), Vec::new());
    for &(pitch, mu, alpha, rpm, cl_meas) in &fwd_pts {
        let fc = wheatley_fc(pitch, mu, alpha, rpm);

        let rotor_off = build_vpm(&pca2, SPR_FWD, WAKE_REVS, 0.50, 0.0);
        let (t0, _q0) = run(&rotor_off, &fc, SPR_FWD, TOTAL_REVS);
        let cl0 = wheatley_cl_from_thrust(t0, mu, alpha, rpm);

        let rotor_on = build_vpm(&pca2, SPR_FWD, WAKE_REVS, 0.50, TAU_REV);
        let (t1, _q1) = run(&rotor_on, &fc, SPR_FWD, TOTAL_REVS);
        let cl1 = wheatley_cl_from_thrust(t1, mu, alpha, rpm);

        let case = format!("fwd mu={mu:.3} a={alpha:.1}");
        r.info(format!("{case} CL_off"), "CL", cl0, cl_meas);
        r.info(format!("{case} CL_on"), "CL", cl1, cl_meas);
        cl_off.push(abs_err_pct(cl0, cl_meas));
        cl_on.push(abs_err_pct(cl1, cl_meas));
    }
    let (f_off, f_on) = (mean(&cl_off), mean(&cl_on));
    println!(
        "  [forward CL] mean|err|  OFF={f_off:.1}%  ON={f_on:.1}%   (aging {})",
        verdict(f_off, f_on)
    );
    // Regime finding: in forward flight the wake convects downstream, so the
    // aged far wake is spurious and fading it HELPS. Guard both the direction
    // and the achieved accuracy.
    r.assert_bool(
        "forward_regime",
        "aging_helps_forward",
        f_on,
        f_off,
        f_on < f_off - 0.2,
        "strength-fade should improve forward-flight CL (far wake convects away)",
    );
    r.assert_bool(
        "forward_regime",
        "aging_on_cl_accurate",
        f_on,
        15.0,
        f_on < 15.0,
        "aging-on forward CL mean|err| should be accurate (<15%)",
    );
}

/// One-word verdict comparing OFF vs ON mean error.
fn verdict(off: f64, on: f64) -> &'static str {
    if !off.is_finite() || !on.is_finite() {
        "n/a"
    } else if on < off - 0.2 {
        "HELPS"
    } else if on > off + 0.2 {
        "HURTS"
    } else {
        "neutral"
    }
}
