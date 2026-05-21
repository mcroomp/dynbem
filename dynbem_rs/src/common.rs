// Shared low-level numerics. Pure Rust, no PyO3 here.

// ---------------------------------------------------------------------------
// Numerical thresholds used across the BEM models.
// Names express WHAT the value guards, not just its magnitude.
// ---------------------------------------------------------------------------

/// Denominator-safety floor for an O(1) f64 quantity. Used wherever we'd
/// otherwise divide by, take a ratio of, or compare to "tiny" -- velocities
/// (m/s on rotors of size ~1 m operating at ~100 rad/s), coefficients (~1),
/// flow angles (rad), etc.
pub const EPS_DENOM: f64 = 1e-9;

/// Below this, the rotor is considered not spinning -- bypass all BEM math.
pub const EPS_OMEGA_R: f64 = 1e-6;

/// Floor for the combined Prandtl tip+hub loss factor. Prevents the
/// 1/F division in the momentum-BEM quadratic from blowing up at the
/// tip / hub-cutout.
pub const MIN_LOSS_FACTOR: f64 = 1e-4;

/// V_T (mass-flow speed at the disk) floor expressed as a fraction of
/// max(Omega*R, 1). Stops the Pitt-Peters / Oye time constants going
/// infinite as the disk approaches VRS / hover at zero thrust.
pub const V_T_HOVER_FLOOR_FRAC: f64 = 1e-2;

/// Threshold below which |v_climb| is considered "definitively in
/// descent" -- guards the VRS-region detection from chattering at hover.
pub const VRS_DESCENT_THRESHOLD: f64 = 1e-3;

/// Floor for the non-dim mass-flow parameter mu_T in the Pitt-Peters L
/// matrix denominator. Without this the cyclic targets blow up at
/// hover/climb where mu_T -> 0. 0.05 corresponds to V_T ~= 5% Omega_R,
/// well below any practical operating point. Empirical, matches the
/// Python reference.
pub const MU_T_FLOOR: f64 = 0.05;

/// Leishman VRS empirical polynomial coefficients (NACA TN-2474 fit).
const VRS_C: [f64; 4] = [1.125, -1.372, 1.718, -0.655];

#[inline]
pub fn vrs_lambda1(k: f64) -> f64 {
    // Horner form: lets LLVM emit a tight FMA chain.
    let k2 = k * k;
    let k3 = k2 * k;
    let k4 = k2 * k2;
    1.0 + VRS_C[0] * k + VRS_C[1] * k2 + VRS_C[2] * k3 + VRS_C[3] * k4
}
