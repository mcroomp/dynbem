// Flap harmonics -- VPM flap DOF vs classical closed-form flapping theory.
//
// Standard theory (Bramwell / Seddon / Prouty), centrally-hinged (nu_beta = 1),
// untwisted, no cyclic, with the inflow lambda taken from the VPM's own thrust
// (Glauert) so the comparison is apples-to-apples on inflow:
//
//   a_0 = gamma [ theta_0/8 (1+mu^2) - lambda/6 ]        (coning)
//   a_1 = 2 mu (4/3 theta_0 - lambda) / (1 - mu^2/2)     (longitudinal)
//   b_1 = 4/3 mu a_0 / (1 + mu^2/2)                      (lateral)
//
// Validates the flap ODE, the aero flap-moment forcing, and the ~90 deg flap
// phase lag (MODEL.md sec 14a is the quasi-static analogue). Coning a_0 and
// longitudinal a_1 are asserted; lateral b_1 is a documented under-prediction
// (VPM_DESIGN sec 8.4 / 8.5) and only sanity-bounded here.

use crate::common::*;

const COLLECTIVE_DEG: f64 = 8.0;
const I_BETA: f64 = 0.030; // Lock number gamma ~ 8 for the Castles-Gray rotor

#[test]
fn flap_harmonics_vs_theory() {
    let defn = theory_rotor_flap(10, I_BETA);
    let rotor = make_rotor(&defn);
    let gamma = lock_number(I_BETA);
    let theta0 = COLLECTIVE_DEG.to_radians();

    for &mu in &[0.15_f64, 0.20, 0.25] {
        let fc = forward_fc(COLLECTIVE_DEG, mu);
        let (res, state) = settle(&rotor, &fc, 10);

        let samples = sample_flap(&rotor, &fc, &state);
        let (a0, a1, b1) = fourier_flap(&samples);

        let ct = c_t(res.thrust);
        let lambda = glauert_lambda(ct, mu, 0.0);
        let a0_th = gamma * (theta0 / 8.0 * (1.0 + mu * mu) - lambda / 6.0);
        let a1_th = 2.0 * mu * (4.0 / 3.0 * theta0 - lambda) / (1.0 - 0.5 * mu * mu);

        let a0_err = (a0 - a0_th).abs() / a0_th.abs();
        let a1_err = (a1 - a1_th).abs() / a1_th.abs();

        eprintln!(
            "flap mu={mu:.2}: a0={:.2} (th {:.2}, {:.0}%)  a1={:.2} (th {:.2}, {:.0}%)  \
             b1={:.2} deg  [C_T={ct:.4}, lambda={lambda:.4}]",
            a0.to_degrees(),
            a0_th.to_degrees(),
            a0_err * 100.0,
            a1.to_degrees(),
            a1_th.to_degrees(),
            a1_err * 100.0,
            b1.to_degrees(),
        );

        assert!(
            a0 > 0.0,
            "coning must be positive (blade cones up), got {a0} rad"
        );
        assert!(
            a0_err < 0.15,
            "coning a0 {:.2} deg vs theory {:.2} deg ({:.0}% off, > 15%)",
            a0.to_degrees(),
            a0_th.to_degrees(),
            a0_err * 100.0
        );
        assert!(
            a1_err < 0.15,
            "longitudinal a1 {:.2} deg vs theory {:.2} deg ({:.0}% off, > 15%)",
            a1.to_degrees(),
            a1_th.to_degrees(),
            a1_err * 100.0
        );
        // Lateral b1 is under-predicted by the free wake vs uniform-inflow
        // theory (documented open item); only bound it as small vs a1.
        assert!(
            b1.abs() < a1.abs(),
            "lateral b1 {:.2} deg should be smaller than longitudinal a1 {:.2} deg",
            b1.to_degrees(),
            a1.to_degrees()
        );
    }
}
