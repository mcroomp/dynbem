// VPM forward-flight autorotation vs Wheatley & Hood NACA TR-515 Tables III & IV.
// Full per-case sweep with error ceilings from vpm_forward_flight_empirical.csv.
// Release mode required -- each case runs ~6 rotor revolutions (~60 s total).
//
// To regenerate ceilings after a model change:
//   REWRITE_VPM_EMPIRICAL=1 cargo run --release -p validation_rs

use crate::helpers::*;
use crate::report::Report;
use serde::{Deserialize, Serialize};

#[derive(Deserialize, Serialize)]
struct Row {
    label: String,
    table: String,
    pitch_deg: f64,
    mu: f64,
    alpha_deg: f64,
    n_rpm: f64,
    cl_meas: f64,
    cl_max_err: f64,
}

pub fn check_vpm_empirical(r: &mut Report) {
    r.begin_module(
        "vpm_forward_flight_empirical",
        "VPM CL vs Wheatley & Hood NACA TR-515 Tables III & IV (full per-case sweep)",
    );

    let csv_data = include_str!("vpm_forward_flight_empirical.csv");
    let rewrite = std::env::var("REWRITE_VPM_EMPIRICAL").is_ok();
    let mut new_rows: Vec<Row> = Vec::new();

    let defn = pca2_rotor(12);
    let rotor = make_pca2_rotor(&defn);

    const STEPS_PER_REV: usize = 48;
    const N_SETTLE: usize = 4;
    const N_AVG: usize = 2;

    for rec in csv_rows(csv_data) {
        let row: Row = rec.deserialize(None).expect("vpm_forward_flight_empirical.csv parse");
        let omega = omega_from_rpm(row.n_rpm);
        let dt = (2.0 * std::f64::consts::PI / omega) / STEPS_PER_REV as f64;
        let n_total = (N_SETTLE + N_AVG) * STEPS_PER_REV;

        let fc = wheatley_fc(row.pitch_deg, row.mu, row.alpha_deg, row.n_rpm);
        let (res, _) = rotor.march(&fc, None, dt, n_total);
        let cl = wheatley_cl_from_thrust(res.thrust, row.mu, row.alpha_deg, row.n_rpm);
        let err_pct = (cl - row.cl_meas).abs() / row.cl_meas * 100.0;

        let case = format!("{} (table {}, mu={:.3})", row.label, row.table, row.mu);
        r.info(&case, "CL_vpm", cl, row.cl_meas);
        r.check(&case, "CL_err_pct", err_pct, 0.0, row.cl_max_err * 100.0);

        if rewrite {
            let new_tol = (err_pct / 100.0 + 0.05).max(err_pct / 100.0 * 3.0);
            new_rows.push(Row { cl_max_err: new_tol, ..row });
        }
    }

    if rewrite {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("src/checks/vpm_forward_flight_empirical.csv");
        let mut wtr = csv::Writer::from_path(&path).expect("open CSV for rewrite");
        for row in &new_rows {
            wtr.serialize(row).expect("write row");
        }
        wtr.flush().expect("flush");
        eprintln!("Rewrote {}", path.display());
    }
}
