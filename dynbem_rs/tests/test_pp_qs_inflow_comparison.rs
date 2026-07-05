mod common;

use dynbem_rs::aero_model::{AeroModel, IntegrationMethod, RotorStateExt};
use dynbem_rs::pitt_peters::PittPetersModel;
use dynbem_rs::quasi_static_bem::QuasiStaticBEM;
use dynbem_rs::polar::LinearPolar;

#[test]
fn pp_vs_qs_steady_state_inflow_4_91deg() {
    // The question: why does PP settle on a different inflow level than QS?
    // In pure hover, both should produce similar thrust.
    
    let defn = common::castles_gray_rotor();
    let inp = common::hover_inputs(4.91, 1200.0);
    let ct_empirical = 0.00168;
    
    let polar_qs = LinearPolar::from_properties(&defn.airfoil);
    let model_qs = QuasiStaticBEM::build(defn.clone(), 36, polar_qs);
    let state_qs = model_qs.initial_state();
    let (res_qs, _) = model_qs.step(&inp, &state_qs, 0.001, IntegrationMethod::ExplicitEuler);
    let ct_qs = common::ct_from_result(&defn, 1200.0, -res_qs.F_world[2]);
    
    eprintln!("QuasiStatic:");
    eprintln!("  CT={:.5}", ct_qs);
    eprintln!("  (no inflow state - QS is algebraic)");
    eprintln!("  Thrust: {:.5} N", -res_qs.F_world[2]);
    eprintln!("  Torque: {:.5} N*m", res_qs.Q_spin);
    
    // Now run PP to steady state
    let polar_pp = LinearPolar::from_properties(&defn.airfoil);
    let model_pp = PittPetersModel::build(defn.clone(), 36, polar_pp);
    let mut state_pp = model_pp.initial_state();
    
    for step in 0..10000 {
        let (res, next_state) = model_pp.step(&inp, &state_pp, 0.001, IntegrationMethod::ExplicitEuler);
        state_pp = next_state;
    }
    
    let (res_pp, _) = model_pp.compute_forces(&inp, &state_pp);
    let ct_pp = common::ct_from_result(&defn, 1200.0, -res_pp.F_world[2]);
    let inflow_pp = state_pp.get_inflow();
    
    eprintln!("\nPittPeters (after 10000 steps):");
    eprintln!("  CT={:.5}", ct_pp);
    eprintln!("  Inflow: lambda_0={:.6}, lambda_c={:.6}, lambda_s={:.6}", inflow_pp[0], inflow_pp[1], inflow_pp[2]);
    eprintln!("  Thrust: {:.5} N", -res_pp.F_world[2]);
    eprintln!("  Torque: {:.5} N*m", res_pp.Q_spin);
    
    eprintln!("\nComparison:");
    eprintln!("  CT difference: {:.5} ({:.1}%)", (ct_pp - ct_qs).abs(), ((ct_pp - ct_qs).abs() / ct_qs) * 100.0);
    eprintln!("  PP inflow lambda_0: {:.6}", inflow_pp[0]);
    let expected_lambda = ((ct_empirical) / 2.0_f64).sqrt();
    eprintln!("  Expected (hover momentum): {:.6}", expected_lambda);
    eprintln!("  Empirical CT: {:.5}", ct_empirical);
}
