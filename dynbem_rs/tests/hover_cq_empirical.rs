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
    cq_empirical: f64,
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
    let mut cq = 0.0;
    for _ in 0..n_steps {
        let (res, next_state) = model.step(inp, &state, dt, method);
        state = next_state;
        cq = common::cq_from_result(defn, rpm, res.Q_spin);
    }
    cq
}

fn cq_hover(model_name: &str, theta_deg: f64, rpm: f64) -> f64 {
    let defn = common::castles_gray_rotor();
    let inp = common::hover_inputs(theta_deg, rpm);
    match model_name {
        "QS" => {
            let polar = LinearPolar::from_properties(&defn.airfoil);
            let model = QuasiStaticBEM::build(defn.clone(), 36, polar);
            run_model_loop(&model, &inp, &defn, rpm, 1, 0.001, IntegrationMethod::ExplicitEuler)
        }
        "PP" => {
            let polar = LinearPolar::from_properties(&defn.airfoil);
            let model = PittPetersModel::build(defn.clone(), 36, polar);
            run_model_loop(&model, &inp, &defn, rpm, 10000, 0.001, IntegrationMethod::ExplicitEuler)
        }
        "OYE" => {
            let polar = LinearPolar::from_properties(&defn.airfoil);
            let model = OyeBEMModel::build(defn.clone(), 36, polar);
            run_model_loop(&model, &inp, &defn, rpm, 10000, 0.001, IntegrationMethod::ExplicitEuler)
        }
        _ => panic!("unknown model: {}", model_name),
    }
}

// Castles-Gray TN-2474 Table V hover CQ vs model predictions.
// Same operating points as hover_empirical (CT); both quantities from the same experiment.
// Per-case ceilings loaded from hover_cq_empirical.csv.
//
// To rewrite CSV with actual errors, set REWRITE_EMPIRICAL_CSV=1:
//   $env:REWRITE_EMPIRICAL_CSV=1; cargo test --test hover_cq_empirical -- --nocapture
#[test]
fn hover_cq_models_vs_empirical() {
    let csv_data = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/hover_cq_empirical.csv"
    ));
    
    let rewrite_mode = env::var("REWRITE_EMPIRICAL_CSV").is_ok();
    
    let mut rdr = common::csv_reader_with_comments(csv_data.as_bytes());
    let mut records: Vec<HoverRecord> = Vec::new();
    let mut failures: Vec<String> = Vec::new();

    for result in rdr.deserialize() {
        let record: HoverRecord = result.expect("failed to deserialize CSV record");
        let cq_model = cq_hover(&record.model, record.theta_deg, record.rpm);
        let err = (cq_model - record.cq_empirical).abs() / record.cq_empirical.abs();

        let (min_err, max_err) = common::error_band(record.max_err);

        eprintln!(
            "{}: model={} theta={} rpm={}: \
             CQ_model={:.7}, CQ_empirical={:.7}, \
             err={:.1}% (band=[{:.1}%, {:.1}%])",
            record.name,
            record.model,
            record.theta_deg,
            record.rpm,
            cq_model,
            record.cq_empirical,
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
        let csv_path = format!("{}/tests/hover_cq_empirical.csv", manifest_dir);
        let mut wtr = csv::Writer::from_path(&csv_path).expect("failed to open CSV for writing");

        for record in records {
            wtr.serialize(record).expect("failed to write record");
        }

        wtr.flush().expect("failed to flush CSV");
        eprintln!("\n=== CSV rewritten ===");
    } else if !failures.is_empty() {
        panic!("hover_cq_empirical failures:\n{}", failures.join("\n"));
    }
}
