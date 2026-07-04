use crate::helpers::*;
use crate::report::Report;

pub fn check_glauert_forward_inflow(r: &mut Report) {
    r.begin_module(
        "glauert_forward_inflow",
        "VPM disk inflow and wake skew vs Glauert; BEM_COMMON.md sec 10-11",
    );
    let defn = theory_rotor(12, 0.0);
    let rotor = make_rotor(&defn);
    let mut lam_errs = Vec::new();
    for &mu in &[0.10_f64, 0.15, 0.20, 0.25] {
        let fc = forward_fc(8.5, mu);
        let (res, state) = settle(&rotor, &fc, 10);
        let ct = c_t(res.thrust);
        let lambda_vpm = mean_disk_inflow(&state).abs();
        let lambda_g = glauert_lambda(ct, mu, 0.0).abs();
        let lam_err = (lambda_vpm - lambda_g).abs() / lambda_g.max(1e-6);
        lam_errs.push(lam_err);
        let chi_vpm = wake_skew_angle(&state, &res);
        let chi_g = mu.atan2(lambda_g.max(1e-6));
        let case = format!("mu={mu:.2}");
        r.check(&case, "lambda_inflow", lambda_vpm, lambda_g, 65.0);
        r.check(
            &case,
            "chi_deg",
            chi_vpm.to_degrees(),
            chi_g.to_degrees(),
            20.0,
        );
        r.info(&case, "CT", ct, f64::NAN);
    }
    let mean_err = lam_errs.iter().sum::<f64>() / lam_errs.len() as f64;
    r.check(
        "aggregate",
        "mean_inflow_err_pct",
        mean_err * 100.0,
        0.0,
        35.0,
    );
}
