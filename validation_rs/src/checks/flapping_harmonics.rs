use crate::helpers::*;
use crate::report::Report;

pub fn check_flapping_harmonics(r: &mut Report) {
    r.begin_module(
        "flapping_harmonics",
        "VPM flap ODE vs Bramwell/Seddon theory; MODEL.md sec 14a",
    );
    const I_BETA: f64 = 0.030;
    let defn = theory_rotor_flap(10, I_BETA);
    let gamma = lock_number(I_BETA);
    let rotor = make_rotor(&defn);
    let theta0 = 8.0f64.to_radians();
    for &mu in &[0.15_f64, 0.20, 0.25] {
        let fc = forward_fc(8.0, mu);
        let (res, state) = settle(&rotor, &fc, 10);
        let samples = sample_flap(&rotor, &fc, &state);
        let (a0, a1, b1) = fourier_flap(&samples);
        let ct = c_t(res.thrust);
        let lambda = glauert_lambda(ct, mu, 0.0);
        let a0_th = gamma * (theta0 / 8.0 * (1.0 + mu * mu) - lambda / 6.0);
        let a1_th = 2.0 * mu * (4.0 / 3.0 * theta0 - lambda) / (1.0 - 0.5 * mu * mu);
        let case = format!("mu={mu:.2}");
        r.assert_bool(
            &case,
            "coning_positive",
            a0,
            0.0,
            a0 > 0.0,
            "coning must be positive",
        );
        r.check(&case, "coning_a0_deg", a0.to_degrees(), a0_th.to_degrees(), 15.0);
        r.check(&case, "longit_a1_deg", a1.to_degrees(), a1_th.to_degrees(), 15.0);
        // b1 (lateral) under-predicted by free wake -- info only, bounded < a1.
        r.assert_bool(
            &case,
            "b1_smaller_than_a1",
            b1.abs(),
            a1.abs(),
            b1.abs() < a1.abs(),
            "b1 should be smaller than a1",
        );
        r.info(&case, "lateral_b1_deg", b1.to_degrees(), f64::NAN);
    }
}
