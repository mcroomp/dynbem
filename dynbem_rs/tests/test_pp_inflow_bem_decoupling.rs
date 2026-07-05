mod common;

use dynbem_rs::aero_model::{AeroModel, IntegrationMethod, RotorStateExt};
use dynbem_rs::pitt_peters::PittPetersModel;
use dynbem_rs::polar::LinearPolar;

#[test]
fn pp_inflow_bem_decoupling_4_91deg() {
    // Check if PP's inflow ODE and BEM are decoupling.
    // In hover, both should converge to: lambda_0 = sqrt(CT/2) and CT from momentum.
    
    let defn = common::castles_gray_rotor();
    let inp = common::hover_inputs(4.91, 1200.0);
    let ct_empirical = 0.00168;
    let expected_lambda = ((ct_empirical) / 2.0_f64).sqrt();
    
    let polar_pp = LinearPolar::from_properties(&defn.airfoil);
    let model_pp = PittPetersModel::build(defn.clone(), 36, polar_pp);
    let mut state_pp = model_pp.initial_state();
    
    eprintln!("Step  | lambda_0  | CT_model  | CT_theory | Difference");
    eprintln!("------|-----------|-----------|-----------|----------");
    
    for step in 0..10000 {
        let (res, next_state) = model_pp.step(&inp, &state_pp, 0.001, IntegrationMethod::ExplicitEuler);
        state_pp = next_state;
        
        let steps_to_check = [1, 10, 100, 500, 1000, 5000, 10000];
        if steps_to_check.contains(&(step + 1)) {
            let inflow = state_pp.get_inflow();
            let ct_model = common::ct_from_result(&defn, 1200.0, -res.F_world[2]);
            
            // From momentum theory: if lambda_0 is the inflow, then CT = 2*lambda_0^2 (in hover, mu_T ≈ lambda_0)
            // But that's assuming the BEM equilibrium satisfies momentum theory.
            // Actually: CT = 2 * lambda_0 * mu_T, and in hover mu_T ≈ lambda_0 for the solution
            // But let's compute what the "theory" CT would be from the inflow
            let lambda_0 = inflow[0];
            let ct_from_momentum = 2.0 * lambda_0 * lambda_0;
            
            eprintln!(
                "{:5} | {:.6} | {:.6} | {:.6} | {:.6}",
                step + 1,
                lambda_0,
                ct_model,
                ct_from_momentum,
                (ct_model - ct_from_momentum).abs()
            );
        }
    }
    
    eprintln!("\nExpected CT from momentum theory: {:.6}", ct_empirical);
    eprintln!("Expected lambda_0 from momentum theory: {:.6}", expected_lambda);
    
    let (res_final, _) = model_pp.compute_forces(&inp, &state_pp);
    let ct_final = common::ct_from_result(&defn, 1200.0, -res_final.F_world[2]);
    let inflow_final = state_pp.get_inflow();
    
    eprintln!("\nFinal PP state:");
    eprintln!("  CT_model: {:.6}", ct_final);
    eprintln!("  lambda_0: {:.6}", inflow_final[0]);
    eprintln!("  Error (model vs empirical): {:.1}%", ((ct_final - ct_empirical).abs() / ct_empirical) * 100.0);
}
