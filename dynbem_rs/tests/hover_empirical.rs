mod common;

use dynbem_rs::aero_model::{AeroModel, IntegrationMethod};
use dynbem_rs::oye::OyeBEMModel;
use dynbem_rs::pitt_peters::PittPetersModel;
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::quasi_static_bem::QuasiStaticBEM;
use serde::{Deserialize, Serialize};
use std::env;

#[derive(Serialize, Deserialize, Clone)]
struct HoverRecord {
    name: String,
    model: String,
    theta_deg: f64,
    rpm: f64,
    ct_empirical: f64,
    max_err: f64,
}

fn run_model_loop<M: AeroModel>(
    model: &M,
    inp: &dynbem_rs::aero_io::RotorInputs,
    defn: &dynbem_rs::rotor_definition::RotorDefinition,
    rpm: f64,
    n_steps: usize,
    dt: f64,
    method: IntegrationMethod,
) -> f64 {
    let mut state = model.initial_state();
    let mut ct = 0.0;
    for _ in 0..n_steps {
        let (res, next_state) = model.step(inp, &state, dt, method);
        state = next_state;
        ct = common::ct_from_result(defn, rpm, -res.F_world[2]);
    }
    ct
}

fn ct_hover(model_name: &str, theta_deg: f64, rpm: f64) -> f64 {
    let defn = common::castles_gray_rotor();
    let inp = common::hover_inputs(theta_deg, rpm);
    let polar = LinearPolar::from_properties(&defn.airfoil);
    match model_name {
        "QS" => {
            let model = QuasiStaticBEM::build(defn.clone(), 36, polar);
            run_model_loop(&model, &inp, &defn, rpm, 1, 0.001, IntegrationMethod::ExplicitEuler)
        }
        "PP" => {
            let model = PittPetersModel::build(defn.clone(), 36, polar);
            run_model_loop(&model, &inp, &defn, rpm, 10000, 0.001, IntegrationMethod::ExplicitEuler)
        }
        "OYE" => {
            let model = OyeBEMModel::build(defn.clone(), 36, polar);
            run_model_loop(&model, &inp, &defn, rpm, 10000, 0.001, IntegrationMethod::ExplicitEuler)
        }
        _ => panic!("unknown model: {}", model_name),
    }
}

// Castles-Gray empirical CT values from Table V (page 36).
// Cases and per-case error ceilings are loaded from hover_empirical.csv.
// Each ceiling is the observed error + ~1 percentage point; a lower bound
// at ceiling - 1.1% catches unexpected improvements.
//
// To rewrite CSV with actual errors, set REWRITE_EMPIRICAL_CSV=1:
//   $env:REWRITE_EMPIRICAL_CSV=1; cargo test --test hover_empirical -- --nocapture
#[test]
fn hover_ct_models_vs_empirical() {
    let csv_data = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/hover_empirical.csv"
    ));
    
    let rewrite_mode = env::var("REWRITE_EMPIRICAL_CSV").is_ok();
    
    let mut rdr = common::csv_reader_with_comments(csv_data.as_bytes());
    let mut records: Vec<HoverRecord> = Vec::new();
    let mut failures: Vec<String> = Vec::new();

    for result in rdr.deserialize() {
        let record: HoverRecord = result.expect("failed to deserialize CSV record");
        let ct_model = ct_hover(&record.model, record.theta_deg, record.rpm);
        let err = (ct_model - record.ct_empirical).abs() / record.ct_empirical;

        let (min_err, max_err) = common::error_band(record.max_err);

        eprintln!(
            "{}: model={} theta={} rpm={}: \
             CT_model={:.5}, CT_empirical={:.5}, \
             err={:.1}% (band=[{:.1}%, {:.1}%])",
            record.name,
            record.model,
            record.theta_deg,
            record.rpm,
            ct_model,
            record.ct_empirical,
            err * 100.0,
            min_err * 100.0,
            max_err * 100.0
        );

        if rewrite_mode {
            let new_max_err = (err + 0.01).max(0.011);
            let mut updated = record.clone();
            updated.max_err = new_max_err;
            records.push(updated);
        } else {
            if err >= max_err {
                failures.push(format!(
                    "{}: err={:.1}% exceeds max {:.1}%",
                    record.name,
                    err * 100.0,
                    max_err * 100.0
                ));
            } else if err < min_err {
                failures.push(format!(
                    "{}: err={:.1}% is below min {:.1}% (improvement detected - please update threshold)",
                    record.name,
                    err * 100.0,
                    min_err * 100.0
                ));
            }
        }
    }

    if rewrite_mode {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let csv_path = format!("{}/tests/hover_empirical.csv", manifest_dir);
        let mut wtr = csv::Writer::from_path(&csv_path).expect("failed to open CSV for writing");

        for record in records {
            wtr.serialize(record).expect("failed to write record");
        }

        wtr.flush().expect("failed to flush CSV");
        eprintln!("\n=== CSV rewritten ===");
    } else if !failures.is_empty() {
        panic!("hover_empirical failures:\n{}", failures.join("\n"));
    }
}
