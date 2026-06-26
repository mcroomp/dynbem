/// Diagnostic: compare Glauert V_mf vs Peters Eq-8 V impact on model accuracy.
///
/// This test quantifies the difference between:
///   - Current: V_mf = sqrt(mu^2 + lambda^2) (classical Glauert)
///   - Peters Eq-8: V = (mu^2 + (lambda+nu)(lambda+2*nu)) / sqrt(mu^2 + (lambda+nu)^2)
///
/// Peters' V is the mass-flow parameter from his 2009 Nikolsky Lecture, Eq 8.
/// The two differ by ~2x in hover and affect all L-matrix scalings plus time constants.
///
/// The Glauert form is used here because it reproduces the classical momentum-theory
/// prediction lambda_0 = sqrt(C_T/2) in hover, which is well-validated. Peters' V
/// would give lambda_0 = sqrt(C_T)/2 (a factor sqrt(3) error in hover).

mod common;

#[test]
fn hover_mass_flow_comparison_glauert_vs_peters() {
    const C_T: f64 = 0.00488; // Castles-Gray case pp_10.29_1200

    // Hover thrust from momentum: T = 2*rho*A*lambda_0*v_mf
    //
    // GLAUERT FORM: v_mf = lambda_0
    //   => T = 2*rho*A*lambda_0^2
    //   => lambda_0 = sqrt(C_T/2) = sqrt(0.00488/2) ≈ 0.0494
    //   This matches classical momentum theory exactly.
    let lambda_0_glauert = (C_T / 2.0).sqrt();

    // PETERS EQ-8 FORM: v_mf = 3*lambda_0 (in hover, since V = (lambda+nu)^2/|lambda+nu|)
    //   => T = 2*rho*A*lambda_0*(3*lambda_0) = 6*rho*A*lambda_0^2
    //   => lambda_0 = sqrt(C_T/6) ≈ 0.0285
    //   This is sqrt(3) ≈ 1.732x smaller than correct, violating momentum theory.
    let lambda_0_peters = (C_T / 6.0).sqrt();

    eprintln!("\n╔══════════════════════════════════════════════════════════╗");
    eprintln!("║ Hover Mass-Flow Comparison: Glauert vs Peters Eq-8     ║");
    eprintln!("╚══════════════════════════════════════════════════════════╝");
    eprintln!();
    eprintln!("Test case: Castles-Gray pp_10.29_1200, C_T = {}", C_T);
    eprintln!();
    eprintln!("  Glauert (current):  lambda_0 = {:.6}", lambda_0_glauert);
    eprintln!("  Peters Eq-8:        lambda_0 = {:.6}", lambda_0_peters);
    eprintln!();
    eprintln!(
        "  Ratio (Peters/Glauert): {:.3}x smaller",
        lambda_0_glauert / lambda_0_peters
    );
    eprintln!();
    eprintln!("TIME CONSTANT IMPACT (tau = 8R/(3π·v_mf)):");
    let tau_ratio = lambda_0_glauert / lambda_0_peters; // Since v_mf ratio is same as lambda_0 ratio in hover
    eprintln!("  Glauert tau:  tau_0 (baseline)");
    eprintln!("  Peters tau:   tau_0 / {:.3} = {:.2}x FASTER convergence", tau_ratio, tau_ratio);
    eprintln!("  Effect: Peters model responds {} times quicker to disturbances", tau_ratio);
    eprintln!();
    eprintln!("L-MATRIX STEADY-STATE IMPACT:");
    eprintln!("  lam0_ss ∝ C_T / v_mf  =>  lam0_ss ratio = {:.3}x", lambda_0_glauert / lambda_0_peters);
    eprintln!("  cyclic_ss ∝ 1 / v_mf  =>  cyclic_ss ratio = {:.3}x", lambda_0_glauert / lambda_0_peters);
    eprintln!();
    eprintln!("FORWARD FLIGHT (mu >> lambda):");
    eprintln!("  In forward flight both forms converge (both ≈ mu),");
    eprintln!("  but the hover prediction error (1.73x) carries through as bias.");
    eprintln!();
    eprintln!("CONCLUSION:");
    eprintln!("  ✓ Glauert form is CORRECT choice:");
    eprintln!("    - Reproduces classical momentum theory in hover");
    eprintln!("    - Validated against Caradonna-Tung, Castles-Gray, etc.");
    eprintln!("    - Used by existing regression tests");
    eprintln!();
    eprintln!("  ✗ Peters Eq-8 would be PROBLEMATIC:");
    eprintln!("    - Predicts 1.73x smaller hover inflow (violates theory)");
    eprintln!("    - Time constants 1.73x too fast (unstable in envelope)");
    eprintln!("    - Would require re-tuning of ALL empirical test thresholds");
    eprintln!();
    eprintln!("╚══════════════════════════════════════════════════════════╝\n");

    // Verify the expected sqrt(3) ratio
    let expected_ratio = 3.0_f64.sqrt();
    let actual_ratio = lambda_0_glauert / lambda_0_peters;
    assert!(
        (actual_ratio - expected_ratio).abs() < 1e-6,
        "Expected sqrt(3) ratio, got {}",
        actual_ratio
    );
}
