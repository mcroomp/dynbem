use crate::helpers::*;
use crate::report::Report;

pub fn check_blade_element_hover(r: &mut Report) {
    r.begin_module(
        "blade_element_hover",
        "VPM vs combined BEMT, hover; Leishman ch.3",
    );
    let defn = theory_rotor(12, 0.0);
    let rotor = make_rotor(&defn);
    for &theta_deg in &[8.46_f64, 10.29] {
        let fc = hover_fc_omega(theta_deg, omega_from_rpm(1200.0));
        let (res, _s) = settle(&rotor, &fc, 10);
        let ct_vpm = ct_at(res.thrust, omega_from_rpm(1200.0), R_TIP);
        let theta = theta_deg.to_radians();
        let lambda = bemt_hover_lambda(theta);
        let ct_be = bemt_hover_ct(theta, lambda);
        let case = format!("theta={theta_deg:.2} rpm=1200");
        r.check(&case, "CT", ct_vpm, ct_be, 25.0);
        r.info(&case, "lambda_BEMT", lambda, f64::NAN);
    }
}
