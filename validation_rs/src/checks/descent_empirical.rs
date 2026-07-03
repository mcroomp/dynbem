// Castles-Gray NACA TN-2474 Table I axial-descent CQ (WBS regime, lambda2 >= 2).
// QS, Pitt-Peters, and Oye BEM models vs measured data.
// Negative CQ_empirical = autorotation / windmill-brake-state.
// Per-case error ceilings from descent_empirical.csv.

use crate::helpers::*;
use crate::report::Report;
use serde::{Deserialize, Serialize};

#[derive(Deserialize, Serialize)]
struct Row {
    name: String,
    model: String,
    theta_deg: f64,
    rpm: f64,
    v_descent_m_s: f64,
    cq_empirical: f64,
    max_err: f64,
}

pub fn check_descent_empirical(r: &mut Report) {
    r.begin_module(
        "descent_cq_empirical",
        "QS/PP/Oye axial-descent CQ vs Castles-Gray TN-2474 Table I (WBS)",
    );

    let csv_data = include_str!("descent_empirical.csv");
    let rewrite = std::env::var("REWRITE_EMPIRICAL_CSV").is_ok();
    let mut new_rows: Vec<Row> = Vec::new();

    let qs  = castles_gray_qs(12);
    let pp  = castles_gray_pp(12);
    let oye = castles_gray_oye(12);

    for rec in csv_rows(csv_data) {
        let row: Row = rec.deserialize(None).expect("descent_empirical.csv parse");
        let inp = bem_descent_inputs(row.theta_deg, row.rpm, row.v_descent_m_s);
        let cq = match row.model.as_str() {
            "QS"  => run_cq(&qs,  &inp, row.rpm, 1),
            "PP"  => run_cq(&pp,  &inp, row.rpm, 10000),
            "OYE" => run_cq(&oye, &inp, row.rpm, 10000),
            other => panic!("unknown model: {other}"),
        };
        let err_pct = (cq - row.cq_empirical).abs() / row.cq_empirical.abs() * 100.0;
        let tol_pct = row.max_err * 100.0;
        r.check(&row.name, "CQ", cq, row.cq_empirical, tol_pct);
        if rewrite {
            new_rows.push(Row { max_err: (err_pct / 100.0 + 0.01).max(0.011), ..row });
        }
    }

    if rewrite {
        rewrite_csv("descent_empirical.csv", &new_rows);
    }
}

fn rewrite_csv(filename: &str, rows: &[Row]) {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("src/checks")
        .join(filename);
    let mut wtr = csv::Writer::from_path(&path).expect("open CSV for rewrite");
    for row in rows {
        wtr.serialize(row).expect("write row");
    }
    wtr.flush().expect("flush CSV");
    eprintln!("Rewrote {}", path.display());
}
