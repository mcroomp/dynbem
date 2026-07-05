//! Wake aging: core spreading and strength fade.
//!
//! The classic engine is inviscid and non-decaying, so a retained wake stays
//! artificially coherent. Two optional, cheap aging mechanisms model the
//! missing physical decay; both are O(N) in-place passes applied once per
//! sub-step by [`crate::vpm`].
//!
//! * [`core_spread`] grows each particle's core by the viscous law
//!   `d(sigma^2)/dt = 2*nu`, i.e. `sigma <- sqrt(sigma^2 + 2*nu*dt)`. This
//!   CONSERVES circulation (alpha is unchanged) -- it only spreads the
//!   vorticity. Note the regularized kernel's far field (r >> sigma) is
//!   sigma-independent, so core spreading mainly softens the NEAR wake.
//! * [`strength_decay`] multiplies every particle strength by a per-step factor
//!   < 1, modelling turbulent breakdown of the aged wake. This is
//!   NON-conservative and directly weakens the far-field influence (which
//!   scales with alpha), so it targets the spurious far-wake bias that
//!   conservative merging and core spreading cannot.

use super::common::ParticleField;

/// Grow every particle core by one viscous step: `sigma^2 += 2*nu*dt`.
/// No-op when `nu <= 0`. Conserves circulation.
pub fn core_spread(field: &mut ParticleField, nu: f64, dt: f64) {
    if nu <= 0.0 || dt <= 0.0 {
        return;
    }
    let d = (2.0 * nu * dt) as f32;
    for s in field.sigma.iter_mut() {
        *s = (*s * *s + d).sqrt();
    }
}

/// Multiply every particle strength by `factor` (a per-step decay < 1).
/// No-op when `factor` is not in (0, 1). Non-conservative (models wake
/// breakdown).
pub fn strength_decay(field: &mut ParticleField, factor: f32) {
    if !(factor > 0.0 && factor < 1.0) {
        return;
    }
    for a in field.ax.iter_mut() {
        *a *= factor;
    }
    for a in field.ay.iter_mut() {
        *a *= factor;
    }
    for a in field.az.iter_mut() {
        *a *= factor;
    }
}
