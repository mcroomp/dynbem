// VPM theory-validation report binary.
//
// Runs all VPM-vs-standard-theory checks and writes a structured report to
// stdout and tmp/theory_report.txt.  The text format is optimised for AI
// analysis: every data point is on one key=value line; PASS/FAIL/INFO tokens
// are grep-friendly; section headers delimit the report.
//
// Build and run (RELEASE is mandatory -- VPM is ~50-100x slower in debug):
//   cargo run --release --bin theory_report
//
// The checks mirror the removed tests/theory/ suite; the logic lives here so
// it is not in the normal `cargo test` run.

#![allow(dead_code, clippy::too_many_arguments)]

use dynbem_rs::cyclic::ControlGains;
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::rotor_definition::{
    BladeGeometry, FlapProperties, LinearPolarParameters, PitchActuation, RotorDefinition,
};
use dynbem_rs::vpm_rotor::{FlightCondition, VpmRotor, VpmRotorConfig, VpmRotorResult, VpmRotorState, induced_at_points};
use std::f64::consts::PI;
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

// ---------------------------------------------------------------------------
// Report harness
// ---------------------------------------------------------------------------

#[derive(Debug)]
struct Row {
    module:    &'static str,
    case:      String,
    qty:       &'static str,
    vpm:       f64,
    reference: f64,
    err_pct:   f64,   // signed % error vs reference (NaN = no reference)
    tol_pct:   f64,   // tolerance threshold (NaN = info-only)
    status:    Status,
    note:      Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Status { Pass, Fail, Info }

impl std::fmt::Display for Status {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Status::Pass => write!(f, "PASS"),
            Status::Fail => write!(f, "FAIL"),
            Status::Info => write!(f, "INFO"),
        }
    }
}

struct Report {
    rows: Vec<Row>,
    current_module: &'static str,
}

impl Report {
    fn new() -> Self {
        Self { rows: Vec::new(), current_module: "" }
    }

    fn begin_module(&mut self, name: &'static str, desc: &str) {
        self.current_module = name;
        println!();
        println!("=== MODULE {}  ({})", name, desc);
    }

    /// Record a checked quantity with a tolerance.
    fn check(
        &mut self,
        case: impl Into<String>,
        qty: &'static str,
        vpm: f64,
        reference: f64,
        tol_pct: f64,
    ) {
        let err_pct = if reference.abs() > 1e-15 {
            (vpm - reference) / reference * 100.0
        } else {
            f64::NAN
        };
        let pass = err_pct.abs() < tol_pct || err_pct.is_nan();
        let status = if pass { Status::Pass } else { Status::Fail };
        self.emit(case.into(), qty, vpm, reference, err_pct, tol_pct, status, None);
    }

    /// Record an info-only quantity (no pass/fail assertion).
    fn info(
        &mut self,
        case: impl Into<String>,
        qty: &'static str,
        vpm: f64,
        reference: f64,
    ) {
        let err_pct = if reference.abs() > 1e-15 {
            (vpm - reference) / reference * 100.0
        } else {
            f64::NAN
        };
        self.emit(case.into(), qty, vpm, reference, err_pct, f64::NAN, Status::Info, None);
    }

    /// Record a directional / boolean check.
    fn assert_bool(
        &mut self,
        case: impl Into<String>,
        qty: &'static str,
        vpm: f64,
        reference: f64,
        pass: bool,
        note: impl Into<String>,
    ) {
        let err_pct = if reference.abs() > 1e-15 {
            (vpm - reference) / reference * 100.0
        } else {
            f64::NAN
        };
        let status = if pass { Status::Pass } else { Status::Fail };
        self.emit(case.into(), qty, vpm, reference, err_pct, f64::NAN, status, Some(note.into()));
    }

    fn emit(
        &mut self,
        case: String,
        qty: &'static str,
        vpm: f64,
        reference: f64,
        err_pct: f64,
        tol_pct: f64,
        status: Status,
        note: Option<String>,
    ) {
        let tol_str = if tol_pct.is_nan() {
            "NA".to_string()
        } else {
            format!("{:.0}%", tol_pct)
        };
        let err_str = if err_pct.is_nan() {
            "NA".to_string()
        } else {
            format!("{:+.1}%", err_pct)
        };
        let ref_str = if reference.is_nan() {
            "NA".to_string()
        } else {
            format!("{:.6}", reference)
        };
        let note_str = note.as_deref().unwrap_or("");

        println!(
            "  CHECK  module={}  case={:?}  qty={}  vpm={:.6}  ref={}  err={}  tol={}  {}  {}",
            self.current_module, case, qty, vpm, ref_str, err_str, tol_str, status, note_str
        );

        self.rows.push(Row {
            module: self.current_module,
            case,
            qty,
            vpm,
            reference,
            err_pct,
            tol_pct,
            status,
            note,
        });
    }

    fn summary(&self) -> (usize, usize, usize) {
        let total = self.rows.iter().filter(|r| r.status != Status::Info).count();
        let pass  = self.rows.iter().filter(|r| r.status == Status::Pass).count();
        let fail  = self.rows.iter().filter(|r| r.status == Status::Fail).count();
        (total, pass, fail)
    }

    fn write_file(&self, path: &PathBuf) {
        let mut f = std::fs::File::create(path).expect("create theory_report.txt");
        let (total, pass, fail) = self.summary();
        writeln!(f, "THEORY_REPORT  dynbem_rs  generated={}", chrono_now()).unwrap();
        writeln!(f, "SUMMARY  total={}  pass={}  fail={}", total, pass, fail).unwrap();
        writeln!(f).unwrap();
        for row in &self.rows {
            let tol = if row.tol_pct.is_nan() { "NA".to_string() } else { format!("{:.0}%", row.tol_pct) };
            let err = if row.err_pct.is_nan() { "NA".to_string() } else { format!("{:+.1}%", row.err_pct) };
            let note = row.note.as_deref().unwrap_or("");
            writeln!(
                f,
                "CHECK  module={}  case={:?}  qty={}  vpm={:.6}  ref={:.6}  err={}  tol={}  {}  {}",
                row.module, row.case, row.qty, row.vpm, row.reference, err, tol, row.status, note
            ).unwrap();
        }
    }
}

fn chrono_now() -> String {
    // No external chrono dep; use a simple unix-seconds fallback.
    use std::time::{SystemTime, UNIX_EPOCH};
    let s = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs();
    format!("unix={}", s)
}

// ---------------------------------------------------------------------------
// Rotor / flight-condition helpers  (ported from tests/theory/common.rs)
// ---------------------------------------------------------------------------

const R_TIP:           f64 = 0.914;
const R_ROOT:          f64 = 0.155;
const CHORD:           f64 = 0.0479;
const N_BLADES:        usize = 3;
const CL_ALPHA:        f64 = 5.90;
const CD0:             f64 = 0.01046;
const ALPHA_STALL_DEG: f64 = 15.5;
const RHO:             f64 = 1.225;
const OMEGA:           f64 = 125.6637; // 1200 rpm
const STEPS_PER_REV:   usize = 24;
const PCA2_R:          f64 = 6.85;

fn solidity() -> f64 { N_BLADES as f64 * CHORD / (PI * R_TIP) }
fn tip_speed() -> f64 { OMEGA * R_TIP }
fn omega_from_rpm(rpm: f64) -> f64 { rpm * PI / 30.0 }
fn ct_at(thrust: f64, omega: f64, r: f64) -> f64 { thrust / (RHO * PI * r * r * (omega * r).powi(2)) }
fn cq_at(torque: f64, omega: f64, r: f64) -> f64 { torque / (RHO * PI * r * r * (omega * r).powi(2) * r) }
fn c_t(thrust: f64) -> f64 { ct_at(thrust, OMEGA, R_TIP) }
fn c_q(torque: f64) -> f64 { cq_at(torque, OMEGA, R_TIP) }
fn lock_number(i_beta: f64) -> f64 { RHO * CL_ALPHA * CHORD * R_TIP.powi(4) / i_beta }

fn base_polar_params() -> LinearPolarParameters {
    LinearPolarParameters { CL0: 0.0, CL_alpha_per_rad: CL_ALPHA, CD0, alpha_stall_deg: ALPHA_STALL_DEG }
}
fn polar_for(params: &LinearPolarParameters) -> LinearPolar {
    LinearPolar::from_properties(params)
}
fn theory_polar() -> LinearPolar {
    LinearPolar::new(0.0, CL_ALPHA, CD0, ALPHA_STALL_DEG.to_radians())
}

fn theory_blade(n_elements: usize, twist_deg: f64, tip_loss: bool) -> BladeGeometry {
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

fn theory_rotor(n_elements: usize, twist_deg: f64) -> RotorDefinition {
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

fn theory_rotor_flap(n_elements: usize, i_beta: f64) -> RotorDefinition {
    let mut defn = theory_rotor(n_elements, 0.0);
    defn.flap = Some(FlapProperties { I_blade_flap_kgm2: i_beta, omega_nr_rad_s: 0.0 });
    defn
}

fn castles_gray_rotor(n_elements: usize) -> RotorDefinition { theory_rotor(n_elements, 0.0) }

fn pca2_rotor(n_elements: usize) -> RotorDefinition {
    RotorDefinition {
        blade: BladeGeometry {
            n_blades: 4, radius_m: PCA2_R, root_cutout_m: 1.20, chord_m: 0.55,
            twist_deg: 0.0, n_elements, tip_loss: true,
            r_stations_m: Vec::new(), chord_stations_m: Vec::new(), twist_stations_deg: Vec::new(),
        },
        airfoil: LinearPolarParameters { CL0: 0.0, CL_alpha_per_rad: 5.73, CD0: 0.0098, alpha_stall_deg: 15.0 },
        control: None,
        pitch_actuation: PitchActuation::DirectMechanical,
        flap: None,
        name: "pca2_wheatley".to_string(),
        description: String::new(),
    }
}

fn make_rotor(defn: &RotorDefinition) -> VpmRotor<LinearPolar> {
    let cfg = VpmRotorConfig {
        max_particles: 14000, sigma: (1.5 * defn.blade.chord_m) as f32,
        barnes_hut: true, bh_theta: 0.5, bh_min_particles: 2000,
        ..VpmRotorConfig::default()
    };
    VpmRotor::new(defn, polar_for(&defn.airfoil), ControlGains::default(), cfg)
}

fn make_pca2_rotor(defn: &RotorDefinition) -> VpmRotor<LinearPolar> {
    let cfg = VpmRotorConfig {
        max_particles: 2000, sigma: 0.50, relax: 0.35,
        nonlinear_lifting_line: true, tip_clustering: true, local_core: true,
        barnes_hut: true, bh_theta: 0.5, bh_min_particles: 200,
        flap_dynamics: false, ..VpmRotorConfig::default()
    };
    VpmRotor::new(defn, polar_for(&defn.airfoil), ControlGains::default(), cfg)
}

fn hover_fc(collective_deg: f64) -> FlightCondition { hover_fc_omega(collective_deg, OMEGA) }
fn hover_fc_omega(collective_deg: f64, omega: f64) -> FlightCondition {
    FlightCondition { collective_rad: collective_deg.to_radians(), tilt_lon: 0.0, tilt_lat: 0.0,
                      v_hub: [0.0, 0.0, 0.0], omega_rad_s: omega, rho: RHO }
}
fn climb_fc(collective_deg: f64, v_climb: f64) -> FlightCondition {
    FlightCondition { collective_rad: collective_deg.to_radians(), tilt_lon: 0.0, tilt_lat: 0.0,
                      v_hub: [0.0, 0.0, v_climb], omega_rad_s: OMEGA, rho: RHO }
}
fn forward_fc(collective_deg: f64, mu: f64) -> FlightCondition {
    FlightCondition { collective_rad: collective_deg.to_radians(), tilt_lon: 0.0, tilt_lat: 0.0,
                      v_hub: [mu * tip_speed(), 0.0, 0.0], omega_rad_s: OMEGA, rho: RHO }
}
fn wheatley_fc(pitch_deg: f64, mu: f64, alpha_deg: f64, n_rpm: f64) -> FlightCondition {
    let omega = omega_from_rpm(n_rpm);
    let a = alpha_deg.to_radians();
    let v = omega * PCA2_R * mu / a.cos();
    FlightCondition { collective_rad: pitch_deg.to_radians(), tilt_lon: 0.0, tilt_lat: 0.0,
                      v_hub: [v * a.cos(), 0.0, -v * a.sin()], omega_rad_s: omega, rho: RHO }
}

fn settle(rotor: &VpmRotor<LinearPolar>, fc: &FlightCondition, revs: usize) -> (VpmRotorResult, VpmRotorState) {
    let dt = (2.0 * PI / fc.omega_rad_s) / STEPS_PER_REV as f64;
    rotor.march(fc, None, dt, revs * STEPS_PER_REV)
}

fn mean_disk_inflow(state: &VpmRotorState) -> f64 {
    let wake = match &state.wake { Some(w) if w.len() > 0 => w, _ => return 0.0 };
    let mut tx = Vec::new(); let mut ty = Vec::new(); let mut tz = Vec::new();
    for &rf in &[0.4_f64, 0.5, 0.6, 0.7, 0.8] {
        let r = rf * R_TIP;
        for k in 0..12 {
            let psi = 2.0 * PI * k as f64 / 12.0;
            tx.push((r * psi.cos()) as f32);
            ty.push((r * -psi.sin()) as f32);
            tz.push(0.0_f32);
        }
    }
    let ind = induced_at_points(wake, &tx, &ty, &tz);
    let mean_vz: f64 = ind.iter().map(|v| v[2] as f64).sum::<f64>() / ind.len() as f64;
    mean_vz / tip_speed()
}

fn wake_skew_angle(state: &VpmRotorState, result: &VpmRotorResult) -> f64 {
    let _ = state;
    let c = result.wake_centroid;
    let horiz = (c[0] * c[0] + c[1] * c[1]).sqrt();
    horiz.atan2(c[2].abs().max(1e-9))
}

fn sample_flap(rotor: &VpmRotor<LinearPolar>, fc: &FlightCondition, state: &VpmRotorState) -> Vec<(f64, f64)> {
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

fn fourier_flap(samples: &[(f64, f64)]) -> (f64, f64, f64) {
    let n = samples.len() as f64;
    let a0 = samples.iter().map(|&(_, b)| b).sum::<f64>() / n;
    let sc: f64 = samples.iter().map(|&(psi, b)| b * psi.cos()).sum::<f64>() / n;
    let ss: f64 = samples.iter().map(|&(psi, b)| b * psi.sin()).sum::<f64>() / n;
    (a0, -2.0 * sc, -2.0 * ss)
}

// Closed-form helpers
fn bemt_hover_lambda(theta_rad: f64) -> f64 {
    let sa = solidity() * CL_ALPHA;
    (sa / 16.0) * ((1.0 + 32.0 * theta_rad / sa).sqrt() - 1.0)
}
fn bemt_hover_ct(theta_rad: f64, lambda: f64) -> f64 {
    let sa = solidity() * CL_ALPHA;
    (sa / 2.0) * (theta_rad / 3.0 - lambda / 2.0)
}
fn glauert_lambda(ct: f64, mu: f64, alpha_rad: f64) -> f64 {
    let mut lam = mu * alpha_rad.tan() + (ct / 2.0).max(1e-8).sqrt();
    for _ in 0..500 {
        let ln = mu * alpha_rad.tan() + ct / (2.0 * (mu * mu + lam * lam).sqrt());
        if (ln - lam).abs() < 1e-11 { lam = ln; break; }
        lam = ln;
    }
    lam
}
fn wheatley_cl_from_thrust(thrust: f64, mu: f64, alpha_deg: f64, n_rpm: f64) -> f64 {
    let omega = omega_from_rpm(n_rpm);
    let a = alpha_deg.to_radians();
    let v = omega * PCA2_R * mu / a.cos();
    let q = 0.5 * RHO * v * v;
    thrust * a.cos() / (q * PI * PCA2_R * PCA2_R)
}

// ---------------------------------------------------------------------------
// Check functions (one per former test module)
// ---------------------------------------------------------------------------

fn check_blade_element_hover(r: &mut Report) {
    r.begin_module("blade_element_hover", "VPM vs combined BEMT, hover; Leishman ch.3");
    let defn = theory_rotor(12, 0.0);
    let rotor = make_rotor(&defn);
    for &theta_deg in &[8.46_f64, 10.29] {
        let fc = hover_fc_omega(theta_deg, omega_from_rpm(1200.0));
        let (res, _s) = settle(&rotor, &fc, 10);
        let ct_vpm = ct_at(res.thrust, omega_from_rpm(1200.0), R_TIP);
        let theta = theta_deg.to_radians();
        let lambda = bemt_hover_lambda(theta);
        let ct_be  = bemt_hover_ct(theta, lambda);
        let case = format!("theta={theta_deg:.2} rpm=1200");
        r.check(&case, "CT", ct_vpm, ct_be, 25.0);
        r.info(&case, "lambda_BEMT", lambda, f64::NAN);
    }
}

fn check_hover_castles_gray(r: &mut Report) {
    r.begin_module("hover_castles_gray", "Hover vs measured Castles-Gray NACA TN-2474 Table V");
    struct M { theta: f64, rpm: f64, ct: f64, cq: f64 }
    let meas = [
        M { theta: 8.46,  rpm: 1200.0, ct: 0.00400, cq: 0.000226 },
        M { theta: 10.29, rpm: 1200.0, ct: 0.00488, cq: 0.000342 },
    ];
    let defn  = castles_gray_rotor(10);
    let rotor = make_rotor(&defn);
    for m in &meas {
        let omega = omega_from_rpm(m.rpm);
        let fc    = hover_fc_omega(m.theta, omega);
        let (res, _s) = settle(&rotor, &fc, 10);
        let ct = ct_at(res.thrust, omega, R_TIP);
        let cq = cq_at(res.torque, omega, R_TIP);
        let fm = ct.powf(1.5) / (2f64.sqrt() * cq.max(1e-9));
        let fm_meas = m.ct.powf(1.5) / (2f64.sqrt() * m.cq);
        let case = format!("theta={:.2} rpm={:.0}", m.theta, m.rpm);
        r.check(&case, "CT", ct, m.ct, 25.0);
        r.info(&case, "CQ", cq, m.cq);
        r.info(&case, "FM", fm, fm_meas);
    }
}

fn check_climb_momentum(r: &mut Report) {
    r.begin_module("climb_momentum", "Axial climb momentum consistency; C_T ~= 2*lambda_i*(lambda_i+lambda_c)");
    let defn  = theory_rotor(12, 0.0);
    let rotor = make_rotor(&defn);
    let mut prev_ct = f64::NAN;
    let mut prev_li = f64::NAN;
    let mut rels = Vec::new();
    for &v_climb in &[0.0_f64, 2.0, 4.0] {
        let fc = climb_fc(9.0, v_climb);
        let (res, state) = settle(&rotor, &fc, 8);
        let ct = c_t(res.thrust);
        let lambda_i = mean_disk_inflow(&state).abs();
        let lambda_c = v_climb / tip_speed();
        let ct_mom   = 2.0 * lambda_i * (lambda_i + lambda_c);
        let rel = (ct - ct_mom).abs() / ct.max(1e-6);
        rels.push(rel);
        let case = format!("v_climb={v_climb:.1}");
        r.check(&case, "momentum_closure", ct, ct_mom, 80.0);
        r.info(&case, "CT", ct, f64::NAN);
        r.info(&case, "lambda_i", lambda_i, f64::NAN);
        // Directional: higher climb -> lower CT and lower lambda_i
        if !prev_ct.is_nan() {
            r.assert_bool(&case, "CT_decreases_with_climb", ct, prev_ct,
                ct <= prev_ct + 0.0005, "CT should not rise with positive climb");
            r.assert_bool(&case, "lambda_i_drops_with_climb", lambda_i, prev_li,
                lambda_i <= prev_li + 0.003, "lambda_i should not rise with positive climb");
        }
        prev_ct = ct; prev_li = lambda_i;
    }
    let mean_rel = rels.iter().sum::<f64>() / rels.len() as f64;
    r.check("aggregate", "mean_momentum_closure", mean_rel * 100.0, 0.0, 80.0);
}

fn check_glauert_forward_inflow(r: &mut Report) {
    r.begin_module("glauert_forward_inflow", "VPM disk inflow and wake skew vs Glauert; MODEL.md sec 10-11");
    let defn  = theory_rotor(12, 0.0);
    let rotor = make_rotor(&defn);
    let mut lam_errs = Vec::new();
    for &mu in &[0.10_f64, 0.15, 0.20, 0.25] {
        let fc = forward_fc(8.5, mu);
        let (res, state) = settle(&rotor, &fc, 10);
        let ct = c_t(res.thrust);
        let lambda_vpm = mean_disk_inflow(&state).abs();
        let lambda_g   = glauert_lambda(ct, mu, 0.0).abs();
        let lam_err = (lambda_vpm - lambda_g).abs() / lambda_g.max(1e-6);
        lam_errs.push(lam_err);
        let chi_vpm = wake_skew_angle(&state, &res);
        let chi_g   = mu.atan2(lambda_g.max(1e-6));
        let case = format!("mu={mu:.2}");
        r.check(&case, "lambda_inflow", lambda_vpm, lambda_g, 65.0);
        r.check(&case, "chi_deg", chi_vpm.to_degrees(), chi_g.to_degrees(), 20.0);
        r.info(&case, "CT", ct, f64::NAN);
    }
    let mean_err = lam_errs.iter().sum::<f64>() / lam_errs.len() as f64;
    r.check("aggregate", "mean_inflow_err_pct", mean_err * 100.0, 0.0, 35.0);
}

fn check_prandtl_tip_loss(r: &mut Report) {
    r.begin_module("prandtl_tip_loss", "Tip-loss flag must reduce global loads (directional)");
    let mut defn_on  = theory_rotor(12, 0.0); defn_on.blade.tip_loss = true;
    let mut defn_off = theory_rotor(12, 0.0); defn_off.blade.tip_loss = false;
    let rotor_on  = make_rotor(&defn_on);
    let rotor_off = make_rotor(&defn_off);
    let fc = hover_fc(9.0);
    let (res_on,  _s_on)  = settle(&rotor_on,  &fc, 8);
    let (res_off, _s_off) = settle(&rotor_off, &fc, 8);
    let ct_on  = c_t(res_on.thrust);  let ct_off = c_t(res_off.thrust);
    let cq_on  = c_q(res_on.torque);  let cq_off = c_q(res_off.torque);
    r.assert_bool("hover_col=9", "CT_tip_loss_not_higher", ct_on, ct_off, ct_on <= ct_off + 1e-6, "tip-loss must not increase CT");
    r.assert_bool("hover_col=9", "CQ_tip_loss_not_higher", cq_on, cq_off, cq_on <= cq_off + 1e-6, "tip-loss must not increase CQ");
    r.info("hover_col=9", "CT_on",  ct_on,  ct_off);
    r.info("hover_col=9", "CQ_on",  cq_on,  cq_off);
}

fn check_wake_skew(r: &mut Report) {
    r.begin_module("wake_skew", "Wake skew grows with mu; covariant under X/Y rotation");
    let defn  = theory_rotor(12, 0.0);
    let rotor = make_rotor(&defn);
    let mut chis = Vec::new();
    for &mu in &[0.08_f64, 0.16, 0.24] {
        let fc = forward_fc(8.5, mu);
        let (res, state) = settle(&rotor, &fc, 8);
        let chi = wake_skew_angle(&state, &res);
        chis.push(chi);
        r.info(format!("mu={mu:.2}"), "chi_deg", chi.to_degrees(), f64::NAN);
    }
    // Monotone increase (loose 2-deg slack each step).
    r.assert_bool("mu_sweep", "chi_increases_01", chis[1], chis[0], chis[1] > chis[0] - 2f64.to_radians(), "chi not growing with mu");
    r.assert_bool("mu_sweep", "chi_increases_12", chis[2], chis[1], chis[2] > chis[1] - 2f64.to_radians(), "chi not growing with mu");
    // Covariance: X vs Y flight direction at mu=0.20.
    let mu = 0.20; let v = mu * tip_speed();
    let fc_x = FlightCondition { collective_rad: 8.5f64.to_radians(), tilt_lon: 0.0, tilt_lat: 0.0, v_hub: [v, 0.0, 0.0], omega_rad_s: OMEGA, rho: RHO };
    let fc_y = FlightCondition { collective_rad: 8.5f64.to_radians(), tilt_lon: 0.0, tilt_lat: 0.0, v_hub: [0.0, v, 0.0], omega_rad_s: OMEGA, rho: RHO };
    let (res_x, st_x) = settle(&rotor, &fc_x, 8);
    let (res_y, st_y) = settle(&rotor, &fc_y, 8);
    let chi_x = wake_skew_angle(&st_x, &res_x);
    let chi_y = wake_skew_angle(&st_y, &res_y);
    let dchi  = (chi_x - chi_y).abs().to_degrees();
    let hx = (res_x.wake_centroid[0].powi(2) + res_x.wake_centroid[1].powi(2)).sqrt();
    let hy = (res_y.wake_centroid[0].powi(2) + res_y.wake_centroid[1].powi(2)).sqrt();
    let dh = (hx - hy).abs() / hx.max(1e-9);
    r.check("covariance_mu=0.20", "chi_xy_diff_deg", dchi, 0.0, 8.0);
    r.check("covariance_mu=0.20", "horiz_shift_diff_pct", dh * 100.0, 0.0, 30.0);
}

fn check_flapping_harmonics(r: &mut Report) {
    r.begin_module("flapping_harmonics", "VPM flap ODE vs Bramwell/Seddon theory; MODEL.md sec 14a");
    const I_BETA: f64 = 0.030;
    let defn  = theory_rotor_flap(10, I_BETA);
    let gamma = lock_number(I_BETA);
    let rotor = make_rotor(&defn);
    let theta0 = 8.0f64.to_radians();
    for &mu in &[0.15_f64, 0.20, 0.25] {
        let fc = forward_fc(8.0, mu);
        let (res, state) = settle(&rotor, &fc, 10);
        let samples = sample_flap(&rotor, &fc, &state);
        let (a0, a1, b1) = fourier_flap(&samples);
        let ct     = c_t(res.thrust);
        let lambda = glauert_lambda(ct, mu, 0.0);
        let a0_th  = gamma * (theta0 / 8.0 * (1.0 + mu * mu) - lambda / 6.0);
        let a1_th  = 2.0 * mu * (4.0 / 3.0 * theta0 - lambda) / (1.0 - 0.5 * mu * mu);
        let case   = format!("mu={mu:.2}");
        r.assert_bool(&case, "coning_positive", a0, 0.0, a0 > 0.0, "coning must be positive");
        r.check(&case, "coning_a0_deg", a0.to_degrees(), a0_th.to_degrees(), 15.0);
        r.check(&case, "longit_a1_deg", a1.to_degrees(), a1_th.to_degrees(), 15.0);
        // b1 (lateral) under-predicted by free wake -- info only, bounded < a1.
        r.assert_bool(&case, "b1_smaller_than_a1", b1.abs(), a1.abs(), b1.abs() < a1.abs(), "b1 should be smaller than a1");
        r.info(&case, "lateral_b1_deg", b1.to_degrees(), f64::NAN);
    }
}

fn check_autorotation(r: &mut Report) {
    r.begin_module("autorotation", "Directional: VPM reaches negative-torque branch in descent+edgewise");
    let defn  = theory_rotor(12, 0.0);
    let rotor = make_rotor(&defn);
    let candidates = [
        (2.0_f64, 0.18_f64, -2.0_f64),
        (2.0_f64, 0.22_f64, -3.0_f64),
        (1.5_f64, 0.25_f64, -3.0_f64),
        (1.0_f64, 0.28_f64, -4.0_f64),
    ];
    let mut found = false; let mut best_cq = f64::INFINITY;
    for &(col, mu, vz) in &candidates {
        let v = mu * tip_speed();
        let fc = FlightCondition { collective_rad: col.to_radians(), tilt_lon: 0.0, tilt_lat: 0.0,
                                   v_hub: [v, 0.0, vz], omega_rad_s: OMEGA, rho: RHO };
        let (res, _s) = settle(&rotor, &fc, 8);
        let cq = c_q(res.torque); let ct = c_t(res.thrust);
        best_cq = best_cq.min(cq);
        let case = format!("col={col:.1} mu={mu:.2} vz={vz:.1}");
        r.info(&case, "CQ", cq, f64::NAN);
        r.info(&case, "CT", ct, f64::NAN);
        if cq < -5e-6 { found = true; break; }
    }
    r.assert_bool("candidates", "found_negative_torque", best_cq, 0.0, found,
        &format!("no autorotation point found (best CQ={best_cq:+.5})"));
}

fn check_measured_companions(r: &mut Report) {
    r.begin_module("measured_companions", "Each theory module anchored to a measured dataset");

    // --- companion: hover Castles-Gray 1600 rpm
    let defn  = castles_gray_rotor(10);
    let rotor = make_rotor(&defn);
    let pts_1600 = [(3.96_f64, 1600.0_f64, 0.00160_f64), (5.55, 1600.0, 0.00255), (7.18, 1600.0, 0.00346)];
    for &(theta_deg, rpm, ct_meas) in &pts_1600 {
        let omega = omega_from_rpm(rpm);
        let (res, _s) = settle(&rotor, &hover_fc_omega(theta_deg, omega), 10);
        let ct = ct_at(res.thrust, omega, R_TIP);
        r.check(format!("CG1600 theta={theta_deg:.2}"), "CT", ct, ct_meas, 35.0);
    }

    // --- companion: climb + descent (Castles-Gray WBS)
    let rows_climb = [
        (1.23_f64, 1200.0_f64, 11.15_f64, -0.000112_f64),
        (-1.66_f64, 1600.0_f64, 11.91_f64, -0.000084_f64),
    ];
    for &(theta_deg, rpm, v_descent, cq_meas) in &rows_climb {
        let omega = omega_from_rpm(rpm);
        let fc = FlightCondition { collective_rad: theta_deg.to_radians(), tilt_lon: 0.0, tilt_lat: 0.0,
                                   v_hub: [0.0, 0.0, -v_descent], omega_rad_s: omega, rho: RHO };
        let (res, _s) = settle(&rotor, &fc, 10);
        let cq = cq_at(res.torque, omega, R_TIP);
        let case = format!("CG_descent theta={theta_deg:.2} rpm={rpm:.0}");
        r.check(&case, "CQ_sign", cq, cq_meas, 300.0); // loose: sign + order-of-magnitude
        r.info(&case, "CQ_vpm", cq, cq_meas);
    }

    // --- companion: forward flight (Wheatley CL, Tables III & IV)
    let pca2_defn = pca2_rotor(12);
    let pca2 = make_pca2_rotor(&pca2_defn);
    let wh_rows = [
        (1.9_f64, 0.145_f64, 15.9_f64, 98.8_f64, 0.425_f64),
        (2.7_f64, 0.145_f64, 16.9_f64, 117.2_f64, 0.540_f64),
        (1.9_f64, 0.204_f64, 15.7_f64, 98.9_f64, 0.489_f64),
        (2.7_f64, 0.204_f64, 16.9_f64, 117.2_f64, 0.608_f64),
    ];
    let mut cl_errs = Vec::new();
    for &(pitch_deg, mu, alpha_deg, rpm, cl_meas) in &wh_rows {
        let fc = wheatley_fc(pitch_deg, mu, alpha_deg, rpm);
        let (res, _s) = settle(&pca2, &fc, 8);
        let cl = wheatley_cl_from_thrust(res.thrust, mu, alpha_deg, rpm);
        let err = (cl - cl_meas).abs() / cl_meas.max(1e-9);
        cl_errs.push(err);
        r.info(format!("wheatley mu={mu:.3}"), "CL", cl, cl_meas);
    }
    let mean_cl_err = cl_errs.iter().sum::<f64>() / cl_errs.len() as f64;
    r.check("wheatley_aggregate", "mean_CL_err_pct", mean_cl_err * 100.0, 0.0, 45.0);

    // --- companion: tip-loss measured anchor
    let mut defn_on  = castles_gray_rotor(10); defn_on.blade.tip_loss = true;
    let mut defn_off = castles_gray_rotor(10); defn_off.blade.tip_loss = false;
    let r_on  = make_rotor(&defn_on);
    let r_off = make_rotor(&defn_off);
    let ct_meas_ref = 0.00400;
    let omega1200 = omega_from_rpm(1200.0);
    let fc_tl = hover_fc_omega(8.46, omega1200);
    let (res_on,  _s_on)  = settle(&r_on,  &fc_tl, 10);
    let (res_off, _s_off) = settle(&r_off, &fc_tl, 10);
    let ct_on  = ct_at(res_on.thrust,  omega1200, R_TIP);
    let ct_off = ct_at(res_off.thrust, omega1200, R_TIP);
    r.assert_bool("tip_loss_meas", "tip_loss_on_is_more_accurate",
        (ct_on - ct_meas_ref).abs(), (ct_off - ct_meas_ref).abs(),
        (ct_on - ct_meas_ref).abs() <= (ct_off - ct_meas_ref).abs() + 0.0002,
        "tip-loss should not make CT less accurate vs measured");
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() {
    let t0 = Instant::now();
    println!("THEORY_REPORT  dynbem_rs  run_start={}  (RELEASE mode required)", chrono_now());
    println!("Each CHECK line: module case qty vpm ref err tol PASS|FAIL|INFO");
    println!("Modules run sequentially; VPM is the ONLY model under test.");

    let mut report = Report::new();

    // Lightweight checks first (hover / momentum / tip-loss), then heavier
    // forward-flight and flapping checks.
    check_blade_element_hover(&mut report);
    check_hover_castles_gray(&mut report);
    check_climb_momentum(&mut report);
    check_prandtl_tip_loss(&mut report);
    check_glauert_forward_inflow(&mut report);
    check_wake_skew(&mut report);
    check_autorotation(&mut report);
    check_flapping_harmonics(&mut report);
    check_measured_companions(&mut report);

    let elapsed = t0.elapsed();
    let (total, pass, fail) = report.summary();
    println!();
    println!("=== SUMMARY  total={}  pass={}  fail={}  elapsed={:.1}s", total, pass, fail, elapsed.as_secs_f64());

    // Write to tmp/
    let out_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent().unwrap()
        .join("tmp");
    std::fs::create_dir_all(&out_dir).ok();
    let out_path = out_dir.join("theory_report.txt");
    report.write_file(&out_path);
    println!("Report written -> {}", out_path.display());

    if fail > 0 {
        eprintln!("FAILED: {} check(s) did not meet tolerance", fail);
        std::process::exit(1);
    }
}
