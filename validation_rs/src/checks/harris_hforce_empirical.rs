// Harris NASA CR-2008-215370 PCA-2 rotor in-plane hub force (H-force) vs
// advance ratio. Primary tabulated data (report pages 504-505, extracted from
// the text layer) lives in
// Research/Harris_CR-2008-215370/pages_504_505_appendix_pca2_coeffs.md.
//
// Method: the PCA-2 rotor autorotates, so blade pitch is not tabulated. For
// each operating point (mu, shaft alpha, RPM) the collective is trimmed so the
// quasi-static BEM matches the tabulated Rotor CT, then the resulting Rotor CH
// is compared against the tabulated value. Trimming to CT isolates the H-force
// model from the unknown pitch (and, since flapping does not change thrust in
// our model, the trim is independent of the flap DOF).
//
// The blade is modelled freely hinged (omega_NR = 0) so the flapping-tilt
// H-force is captured; for a free hinge that force is independent of blade
// inertia (see pca2_rotor_flap docs). Coefficients use the Harris/American
// normalization C = force / (rho * pi * R^2 * (Omega*R)^2).
//
// Per-case error ceilings live in harris_hforce_empirical.csv. To re-baseline
// after a model change:
//   REWRITE_EMPIRICAL_CSV=1 cargo run --release -p validation_rs -- harris_hforce

use crate::helpers::*;
use crate::report::Report;
use dynbem_rs::aero_model::AeroModel;
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::quasi_static_bem::QuasiStaticBEM;
use serde::{Deserialize, Serialize};

#[derive(Deserialize, Serialize)]
struct Row {
    label: String,
    mu: f64,
    alpha_deg: f64,
    n_rpm: f64,
    ct_meas: f64,
    ch_meas: f64,
    ch_max_err: f64,
}

/// Compute (CT, CH) for the PCA-2 rotor at a given collective pitch and
/// operating point. CH is the in-plane hub-force coefficient (magnitude of the
/// disk-plane force), CT the thrust coefficient, both Harris-normalized.
fn ct_ch_at(rotor: &QuasiStaticBEM<LinearPolar>, pitch_deg: f64, row: &Row) -> (f64, f64) {
    let inp = harris_pca2_inputs(pitch_deg, row.mu, row.alpha_deg, row.n_rpm);
    let (res, _) = rotor.compute_forces(&inp, &rotor.initial_state());
    let f = res.F_world.0;
    let norm = pca2_coeff_norm(row.n_rpm);
    let ct = -f[2] / norm;
    let ch = f[0].hypot(f[1]) / norm;
    (ct, ch)
}

/// Bisection trim: find the collective (deg) whose BEM CT matches ct_target.
/// CT is monotone-increasing in pitch over the bracket. Returns the trimmed
/// pitch and the (CT, CH) achieved there.
fn trim_to_ct(rotor: &QuasiStaticBEM<LinearPolar>, row: &Row) -> (f64, f64, f64) {
    let mut lo = -10.0_f64;
    let mut hi = 20.0_f64;
    let (ct_lo, ch_lo) = ct_ch_at(rotor, lo, row);
    let (ct_hi, ch_hi) = ct_ch_at(rotor, hi, row);
    // Guard: if the target is outside the bracket, clamp to the nearer edge.
    if row.ct_meas <= ct_lo {
        return (lo, ct_lo, ch_lo);
    }
    if row.ct_meas >= ct_hi {
        return (hi, ct_hi, ch_hi);
    }
    let mut mid = 0.5 * (lo + hi);
    let (mut ct_mid, mut ch_mid) = ct_ch_at(rotor, mid, row);
    for _ in 0..50 {
        if (ct_mid - row.ct_meas).abs() < 1e-6 {
            break;
        }
        if ct_mid < row.ct_meas {
            lo = mid;
        } else {
            hi = mid;
        }
        mid = 0.5 * (lo + hi);
        let pair = ct_ch_at(rotor, mid, row);
        ct_mid = pair.0;
        ch_mid = pair.1;
    }
    (mid, ct_mid, ch_mid)
}

pub fn check_harris_hforce(r: &mut Report) {
    r.begin_module(
        "harris_hforce_empirical",
        "QS BEM Rotor CH vs Harris CR-2008-215370 PCA-2 (trim-to-CT, per-case sweep)",
    );

    let csv_data = include_str!("harris_hforce_empirical.csv");
    let rewrite = std::env::var("REWRITE_EMPIRICAL_CSV").is_ok();
    let mut new_rows: Vec<Row> = Vec::new();

    // Freely-hinged PCA-2 blade; i_beta value is irrelevant for a free hinge
    // (see pca2_rotor_flap docs), Lock ~5 chosen only as a plausible number.
    let i_beta = 1700.0;
    let defn = pca2_rotor_flap(16, i_beta);
    let polar = polar_for(&defn.airfoil);
    let rotor = QuasiStaticBEM::build(defn, 48, polar);

    for rec in csv_rows(csv_data) {
        let row: Row = rec
            .deserialize(None)
            .expect("harris_hforce_empirical.csv parse");
        let (pitch, ct, ch) = trim_to_ct(&rotor, &row);
        let err_pct = (ch - row.ch_meas).abs() / row.ch_meas * 100.0;

        let case = format!("{} (mu={:.3})", row.label, row.mu);
        r.info(&case, "trim_pitch_deg", pitch, f64::NAN);
        r.info(&case, "CT_trimmed", ct, row.ct_meas);
        r.info(&case, "CH_bem", ch, row.ch_meas);
        r.check(&case, "CH_err_pct", err_pct, 0.0, row.ch_max_err * 100.0);

        if rewrite {
            let new_tol = (err_pct / 100.0 + 0.03).max(0.03);
            new_rows.push(Row {
                ch_max_err: new_tol,
                ..row
            });
        }
    }

    if rewrite {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("src/checks/harris_hforce_empirical.csv");
        let mut wtr = csv::Writer::from_path(&path).expect("open CSV for rewrite");
        for row in &new_rows {
            wtr.serialize(row).expect("write row");
        }
        wtr.flush().expect("flush");
        eprintln!("Rewrote {}", path.display());
    }
}
