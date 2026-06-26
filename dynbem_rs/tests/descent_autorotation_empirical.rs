mod common;

use dynbem_rs::aero_model::{AeroModel, IntegrationMethod};
use dynbem_rs::oye::OyeBEMModel;
use dynbem_rs::pitt_peters::PittPetersModel;
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::quasi_static_bem::QuasiStaticBEM;
use serde::{Deserialize, Serialize};
use std::env;

#[derive(Serialize, Deserialize, Clone)]
struct DescentRecord {
    name: String,
    model: String,
    theta_deg: f64,
    rpm: f64,
    v_descent_m_s: f64,
    cq_empirical: f64,
    max_err: f64,
}

fn run_descent_step<M: AeroModel>(
    model: &M,
    inp: &dynbem_rs::aero_io::RotorInputs,
    rpm: f64,
    n_steps: usize,
) -> f64 {
    let mut state = model.initial_state();
    let mut cq = 0.0;
    let defn = common::castles_gray_rotor();
    for _ in 0..n_steps {
        let (res, next_state) = model.step(inp, &state, 0.001, IntegrationMethod::ExplicitEuler);
        state = next_state;
        cq = common::cq_from_result(&defn, rpm, res.Q_spin);
    }
    cq
}

fn cq_descent(model_name: &str, theta_deg: f64, rpm: f64, v_descent_m_s: f64) -> f64 {
    let defn = common::castles_gray_rotor();
    let inp = common::descent_inputs(theta_deg, rpm, v_descent_m_s);
    match model_name {
        "QS" => {
            let polar = LinearPolar::from_properties(&defn.airfoil);
            let model = QuasiStaticBEM::build(defn, 36, polar);
            run_descent_step(&model, &inp, rpm, 1)
        }
        "PP" => {
            let polar = LinearPolar::from_properties(&defn.airfoil);
            let model = PittPetersModel::build(defn, 36, polar);
            run_descent_step(&model, &inp, rpm, 10000)
        }
        "OYE" => {
            let polar = LinearPolar::from_properties(&defn.airfoil);
            let model = OyeBEMModel::build(defn, 36, polar);
            run_descent_step(&model, &inp, rpm, 10000)
        }
        _ => panic!("unknown model: {}", model_name),
    }
}

// Castles-Gray TN-2474 Table I axial-descent data, WBS regime (lambda2 >= 2.0).
// All three aero models (QS, PP, OYE); single step for QS, 10000 steps for PP/OYE.
// Negative CQ_empirical means the air drives the rotor (autorotation / WBS).
// Operating points and per-case ceilings are loaded from descent_autorotation_empirical.csv.
//
// To rewrite CSV with actual errors, set REWRITE_EMPIRICAL_CSV=1:
//   $env:REWRITE_EMPIRICAL_CSV=1; cargo test --test descent_autorotation_empirical -- --nocapture
#[test]
fn descent_autorotation_cq_vs_empirical() {
    let csv_data = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/descent_autorotation_empirical.csv"
    ));
    
    let rewrite_mode = env::var("REWRITE_EMPIRICAL_CSV").is_ok();
    
    let mut rdr = common::csv_reader_with_comments(csv_data.as_bytes());
    let mut records: Vec<DescentRecord> = Vec::new();
    let mut failures: Vec<String> = Vec::new();

    for result in rdr.deserialize() {
        let record: DescentRecord = result.expect("failed to deserialize CSV record");
        let err = (cq_descent(
            &record.model,
            record.theta_deg,
            record.rpm,
            record.v_descent_m_s,
        ) - record.cq_empirical)
            .abs()
            / record.cq_empirical.abs();

        let (min_err, max_err) = common::error_band(record.max_err);

        eprintln!(
            "{}: model={} theta={} rpm={} v_descent={:.2}: \
             err={:.1}% (band=[{:.1}%, {:.1}%])",
            record.name,
            record.model,
            record.theta_deg,
            record.rpm,
            record.v_descent_m_s,
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
        let csv_path = format!("{}/tests/descent_autorotation_empirical.csv", manifest_dir);
        let mut wtr = csv::Writer::from_path(&csv_path).expect("failed to open CSV for writing");

        for record in records {
            wtr.serialize(record).expect("failed to write record");
        }

        wtr.flush().expect("failed to flush CSV");
        eprintln!("\n=== CSV rewritten ===");
    } else if !failures.is_empty() {
        panic!("descent_autorotation_empirical failures:\n{}", failures.join("\n"));
    }
}
