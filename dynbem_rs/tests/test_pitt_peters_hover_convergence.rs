mod common;

use dynbem_rs::aero_model::{AeroModel, IntegrationMethod, RotorStateExt};
use dynbem_rs::pitt_peters::PittPetersModel;
use dynbem_rs::quasi_static_bem::QuasiStaticBEM;
use dynbem_rs::polar::LinearPolar;

#[test]
fn pp_vs_qs_convergence_4_91deg_1200rpm() {
    // Case: 4.91 deg, 1200 RPM
    // QS gets 16% error, PP gets 45.4% error - what's different?
    
    let defn = common::castles_gray_rotor();
    let inp = common::hover_inputs(4.91, 1200.0);
    let ct_empirical = 0.00168;
    
    // Run QuasiStatic
    let polar_qs = LinearPolar::from_properties(&defn.airfoil);
    let model_qs = QuasiStaticBEM::build(defn.clone(), 36, polar_qs);
    let mut state_qs = model_qs.initial_state();
    let (res_qs, _) = model_qs.step(&inp, &state_qs, 0.001, IntegrationMethod::ExplicitEuler);
    let ct_qs = common::ct_from_result(&defn, 1200.0, -res_qs.F_world[2]);
    let err_qs = (ct_qs - ct_empirical).abs() / ct_empirical;
    
    eprintln!("QS: CT={:.5}, err={:.1}%", ct_qs, err_qs * 100.0);
    eprintln!("QS thrust: {:.5} N", -res_qs.F_world[2]);
    eprintln!("QS torque: {:.5} N*m", res_qs.Q_spin);
    
    // Run PittPeters with convergence tracking
    let polar_pp = LinearPolar::from_properties(&defn.airfoil);
    let model_pp = PittPetersModel::build(defn.clone(), 36, polar_pp);
    let mut state_pp = model_pp.initial_state();
    
    let n_steps = 10000;
    let dt = 0.001;
    let check_steps = [1, 10, 100, 1000, 10000];
    
    for step in 0..n_steps {
        let (res, next_state) = model_pp.step(&inp, &state_pp, dt, IntegrationMethod::ExplicitEuler);
        state_pp = next_state;
        
        if check_steps.contains(&(step + 1)) {
            let inflow = state_pp.get_inflow();
            let ct_pp = common::ct_from_result(&defn, 1200.0, -res.F_world[2]);
            let err_pp = (ct_pp - ct_empirical).abs() / ct_empirical;
            
            eprintln!(
                "PP step {}: CT={:.5}, err={:.1}%, lambda_0={:.6}, lambda_c={:.6}, lambda_s={:.6}",
                step + 1,
                ct_pp,
                err_pp * 100.0,
                inflow[0],
                inflow[1],
                inflow[2]
            );
        }
    }
    
    let (res_pp, _) = model_pp.compute_forces(&inp, &state_pp);
    let ct_pp = common::ct_from_result(&defn, 1200.0, -res_pp.F_world[2]);
    let err_pp = (ct_pp - ct_empirical).abs() / ct_empirical;
    
    eprintln!("\nFinal comparison:");
    eprintln!("QS: CT={:.5} ({:.1}% error)", ct_qs, err_qs * 100.0);
    eprintln!("PP: CT={:.5} ({:.1}% error)", ct_pp, err_pp * 100.0);
    eprintln!("PP thrust: {:.5} N", -res_pp.F_world[2]);
    eprintln!("PP torque: {:.5} N*m", res_pp.Q_spin);
    eprintln!("Empirical CT: {:.5}", ct_empirical);
}
