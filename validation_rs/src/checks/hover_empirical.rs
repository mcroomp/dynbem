// Castles-Gray NACA TN-2474 Table V hover CT: QS, Pitt-Peters, and Oye BEM
// models vs measured data.  Per-case error ceilings from hover_empirical.csv.
//
// To regenerate bounds after a model change, run:
//   REWRITE_EMPIRICAL_CSV=1 cargo run --release -p validation_rs

use crate::helpers::*;
use crate::report::Report;
use serde::{Deserialize, Serialize};

#[derive(Deserialize, Serialize)]
struct Row {
    name: String,
    model: String,
    theta_deg: f64,
    rpm: f64,
    ct_empirical: f64,
    max_err: f64,
}

pub fn check_hover_empirical(r: &mut Report) {
    r.begin_module(
        "hover_ct_empirical",
        "QS/PP/Oye hover CT vs Castles-Gray TN-2474 Table V",
    );

    let csv_data = include_str!("hover_empirical.csv");
    let rewrite = std::env::var("REWRITE_EMPIRICAL_CSV").is_ok();
    let mut new_rows: Vec<Row> = Vec::new();

    let qs  = castles_gray_qs(12);
    let pp  = castles_gray_pp(12);
    let oye = castles_gray_oye(12);

    for rec in csv_rows(csv_data) {
        let row: Row = rec.deserialize(None).expect("hover_empirical.csv parse");
        let inp = bem_hover_inputs(row.theta_deg, row.rpm);
        let ct = match row.model.as_str() {
            "QS"  => run_ct(&qs,  &inp, row.rpm, 1),
            "PP"  => run_ct(&pp,  &inp, row.rpm, 10000),
            "OYE" => run_ct(&oye, &inp, row.rpm, 10000),
            other => panic!("unknown model: {other}"),
        };
        let err_pct = (ct - row.ct_empirical).abs() / row.ct_empirical * 100.0;
        let tol_pct = row.max_err * 100.0;
        r.check(&row.name, "CT", ct, row.ct_empirical, tol_pct);
        if rewrite {
            new_rows.push(Row { max_err: (err_pct / 100.0 + 0.01).max(0.011), ..row });
        }
    }

    if rewrite {
        rewrite_csv("hover_empirical.csv", &new_rows);
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
