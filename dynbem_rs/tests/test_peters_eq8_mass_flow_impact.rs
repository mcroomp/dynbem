/// Diagnostic: compare Glauert V_mf vs Peters Eq-8 V impact on model accuracy.
///
/// This test quantifies the difference between:
///   - Current: V_mf = sqrt(mu^2 + lambda^2) (classical Glauert)
///   - Peters Eq-8: V = (mu^2 + (lambda+nu)(lambda+2*nu)) / sqrt(mu^2 + (lambda+nu)^2)
///
/// Peters' V is the mass-flow parameter from his 2009 Nikolsky Lecture, Eq 8
/// (see Research/Peters_Nikolsky_2008/CLAUDE.md). In Eq 8, `lambda` is the
/// CLIMB ratio (= 0 in hover) and `nu` is the induced flow -- they are DISTINCT
/// variables. In hover (mu=0, lambda=0) the formula reduces to
///   V = (nu)(2*nu) / sqrt(nu^2) = 2*nu.
///
/// The Glauert form is used in the model because it reproduces the classical
/// momentum-theory hover inflow lambda_0 = sqrt(C_T/2). Peters' V gives
/// lambda_0 = sqrt(C_T/4) = sqrt(C_T)/2 in hover -- a factor sqrt(2) smaller.

mod common;

#[test]
fn hover_mass_flow_comparison_glauert_vs_peters() {
    const C_T: f64 = 0.00488; // Castles-Gray case pp_10.29_1200

    // Hover thrust from momentum: T = 2*rho*A*lambda_0*v_mf, i.e.
    // C_T = 2*lambda_0*(v_mf / (Omega*R)). In hover the non-dim mass-flow
    // equals lambda_0 (Glauert) or 2*lambda_0 (Peters Eq-8).
    //
    // GLAUERT FORM: v_mf = lambda_0
    //   => C_T = 2*lambda_0^2
    //   => lambda_0 = sqrt(C_T/2) = sqrt(0.00488/2) ~ 0.0494
    //   This matches classical momentum theory exactly.
    let lambda_0_glauert = (C_T / 2.0).sqrt();

    // PETERS EQ-8 FORM: in hover V = 2*nu = 2*lambda_0
    //   => C_T = 2*lambda_0*(2*lambda_0) = 4*lambda_0^2
    //   => lambda_0 = sqrt(C_T/4) = sqrt(C_T)/2 ~ 0.0349
    //   This is sqrt(2) ~ 1.414x smaller than the Glauert/momentum value.
    let lambda_0_peters = (C_T / 4.0).sqrt();

    eprintln!();
    eprintln!("=== Hover Mass-Flow Comparison: Glauert vs Peters Eq-8 ===");
    eprintln!();
    eprintln!("Test case: Castles-Gray pp_10.29_1200, C_T = {}", C_T);
    eprintln!();
    eprintln!("  Glauert (current):  lambda_0 = {:.6}", lambda_0_glauert);
    eprintln!("  Peters Eq-8:        lambda_0 = {:.6}", lambda_0_peters);
    eprintln!();
    eprintln!(
        "  Ratio (Glauert/Peters): {:.3}x  (expected sqrt(2) ~ 1.414)",
        lambda_0_glauert / lambda_0_peters
    );
    eprintln!();
    // Mass-flow parameter ratio in hover: v_mf_glauert = lambda_0_glauert,
    // V_peters = 2*lambda_0_peters = sqrt(C_T). So V_peters / v_mf_glauert
    // = sqrt(C_T) / sqrt(C_T/2) = sqrt(2): Peters' mass flow is LARGER.
    let v_mf_glauert = lambda_0_glauert;
    let v_peters = 2.0 * lambda_0_peters;
    let massflow_ratio = v_peters / v_mf_glauert;
    eprintln!("TIME CONSTANT IMPACT (tau = 8R/(3*pi*v_mf)):");
    eprintln!("  v_mf ratio (Peters/Glauert): {:.3}x", massflow_ratio);
    eprintln!(
        "  Peters tau = tau_0 / {:.3} => {:.2}x FASTER convergence",
        massflow_ratio, massflow_ratio
    );
    eprintln!();
    eprintln!("L-MATRIX STEADY-STATE IMPACT:");
    eprintln!(
        "  lam0_ss propto C_T / v_mf  => Peters value {:.3}x smaller",
        massflow_ratio
    );
    eprintln!(
        "  cyclic_ss propto 1 / v_mf  => Peters value {:.3}x smaller",
        massflow_ratio
    );
    eprintln!();
    eprintln!("FORWARD FLIGHT (mu >> lambda):");
    eprintln!("  Both forms converge to V ~ mu, so the difference is a");
    eprintln!("  low-speed/hover effect that biases the hover prediction.");
    eprintln!();
    eprintln!("CONCLUSION:");
    eprintln!("  Glauert form reproduces classical momentum theory in hover");
    eprintln!("  (lambda_0 = sqrt(C_T/2)) and is calibrated against the");
    eprintln!("  46-point Castles-Gray regression suite. Peters' Eq-8 V would");
    eprintln!("  give sqrt(2)x smaller hover inflow and sqrt(2)x faster time");
    eprintln!("  constants, requiring a full re-validation.");
    eprintln!();

    // Verify the expected sqrt(2) ratio between the two hover inflow values.
    let expected_ratio = 2.0_f64.sqrt();
    let actual_ratio = lambda_0_glauert / lambda_0_peters;
    assert!(
        (actual_ratio - expected_ratio).abs() < 1e-6,
        "Expected sqrt(2) ratio, got {}",
        actual_ratio
    );

    // And verify the mass-flow parameter ratio is also sqrt(2).
    assert!(
        (massflow_ratio - expected_ratio).abs() < 1e-6,
        "Expected sqrt(2) mass-flow ratio, got {}",
        massflow_ratio
    );
}
