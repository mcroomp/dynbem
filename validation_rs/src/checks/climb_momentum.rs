use crate::helpers::*;
use crate::report::Report;

pub fn check_climb_momentum(r: &mut Report) {
    r.begin_module(
        "climb_momentum",
        "Axial climb momentum consistency; C_T ~= 2*lambda_i*(lambda_i+lambda_c)",
    );
    let defn = theory_rotor(12, 0.0);
    let rotor = make_rotor(&defn);
    let mut prev_ct = f64::NAN;
    let mut prev_li = f64::NAN;
    let mut rels = Vec::new();
    for &v_climb in &[0.0_f64, 2.0, 4.0] {
        let fc = climb_fc(9.0, v_climb);
        let (res, state) = settle(&rotor, &fc, 8);
        let ct = c_t(res.thrust);
        let lambda_i = mean_disk_inflow(&state).abs();
        let lambda_c = v_climb / tip_speed();
        let ct_mom = 2.0 * lambda_i * (lambda_i + lambda_c);
        let rel = (ct - ct_mom).abs() / ct.max(1e-6);
        rels.push(rel);
        let case = format!("v_climb={v_climb:.1}");
        r.check(&case, "momentum_closure", ct, ct_mom, 80.0);
        r.info(&case, "CT", ct, f64::NAN);
        r.info(&case, "lambda_i", lambda_i, f64::NAN);
        // Directional: higher climb -> lower CT and lower lambda_i
        if !prev_ct.is_nan() {
            r.assert_bool(
                &case,
                "CT_decreases_with_climb",
                ct,
                prev_ct,
                ct <= prev_ct + 0.0005,
                "CT should not rise with positive climb",
            );
            r.assert_bool(
                &case,
                "lambda_i_drops_with_climb",
                lambda_i,
                prev_li,
                lambda_i <= prev_li + 0.003,
                "lambda_i should not rise with positive climb",
            );
        }
        prev_ct = ct;
        prev_li = lambda_i;
    }
    let mean_rel = rels.iter().sum::<f64>() / rels.len() as f64;
    r.check(
        "aggregate",
        "mean_momentum_closure",
        mean_rel * 100.0,
        0.0,
        80.0,
    );
}
