//! Reformulated VPM (rVPM) wake evolution -- the alternative wake engine to the
//! classic convection-only path in [`crate::vpm`].
//!
//! # Relationship to the classic engine
//!
//! rVPM (Alvarez & Ning 2020/2022) is a strict generalization of the classic
//! VPM: same regularized particle field ([`ParticleField`]) and same algebraic
//! Biot-Savart kernel, but the wake particles evolve their vector strength
//! `Gamma` and core size `sigma` according to vortex stretching, instead of
//! being frozen. Setting the stretching term to zero recovers the classic
//! convection step exactly.
//!
//! This module therefore reuses [`ParticleField`] from [`crate::vpm`] and adds
//! only the two things classic VPM lacks:
//!   1. the analytic velocity gradient (vortex stretching `S = (Gamma . grad) u`),
//!   2. the rVPM ODE integration of `x`, `Gamma`, `sigma`.
//!
//! The rotor coupling (shedding, lifting line) in [`crate::vpm`] is
//! unchanged; it simply calls [`advect_rvpm`] instead of `vpm::advect_rk2` when
//! the caller selects `WakeEngine::ReformulatedVpm`.
//!
//! # Governing equations (inviscid, momentum + mass conserving: g = 1/5, f = 0)
//!
//! With `S_p = (Gamma_p . grad) u(x_p)` the vortex stretching at particle `p`:
//!
//! ```text
//!   dx_p/dt     = u(x_p)
//!   dGamma_p/dt = S_p - (3/5) (S_p . Ghat_p) Ghat_p
//!   dsigma_p/dt = -(1/5) sigma_p (S_p . Gamma_p) / |Gamma_p|^2
//! ```
//!
//! The SFS/LES turbulence term and the Pedrizzetti relaxation of the FLOWVPM
//! formulation are intentionally omitted from this first implementation; this
//! is the "inviscid rVPM" that provides the strength/size evolution (and the
//! conservation properties) without the turbulence closure.
//!
//! # Analytic gradient of the algebraic kernel
//!
//! The classic velocity is `u = (1/4pi) sum_p K (alpha_p x d)` with
//! `d = x - x_p`, `s = |d|^2 / sigma_p^2`, and
//! `K(s) = (s + 5/2) / (sigma_p^3 (s + 1)^{5/2})`. Differentiating along the
//! target strength `Gamma` gives the per-source stretching contribution
//!
//! ```text
//!   S_contrib = 2 (dK/ds) (Gamma . d) / sigma_p^2 (alpha_p x d)
//!             + K (alpha_p x Gamma)
//!   dK/ds     = (-3/2 s - 21/4) / (sigma_p^3 (s + 1)^{7/2})
//! ```
//!
//! The self term (`p` == target) vanishes: `d = 0` kills the first term and
//! `alpha_p x Gamma = Gamma x Gamma = 0` kills the second.

use super::common::ParticleField;
use rayon::prelude::*;

/// 1 / (4 pi), the Biot-Savart prefactor (f64).
const INV_4PI_F64: f64 = 0.079_577_471_545_947_67;

/// Lower bound on `sigma` after a step (cores must stay positive; keeps
/// `sigma^3` finite in the kernel).
const SIGMA_FLOOR: f64 = 1.0e-4;

/// Below this `|Gamma|^2` a particle is treated as strengthless: the stretching
/// reorientation and the `sigma` evolution (both divide by `|Gamma|^2`) are
/// skipped and the raw stretching is used for `dGamma`.
const GAMMA2_FLOOR: f64 = 1.0e-20;

/// Induced velocity `u` and vortex stretching `S = (Gamma . grad) u` at every
/// particle, evaluated pair-by-pair in double precision. `Gamma` at the target
/// is that particle's own strength. Direct O(N^2); the outer (target) loop runs
/// on the Rayon pool. Returns one `(u, S)` pair per particle, in field order.
pub fn vel_and_stretch(field: &ParticleField) -> Vec<([f64; 3], [f64; 3])> {
    let n = field.len();
    if n == 0 {
        return Vec::new();
    }
    (0..n)
        .into_par_iter()
        .map(|j| {
            let xj = field.px[j] as f64;
            let yj = field.py[j] as f64;
            let zj = field.pz[j] as f64;
            // Target strength Gamma (for the directional derivative).
            let gx = field.ax[j] as f64;
            let gy = field.ay[j] as f64;
            let gz = field.az[j] as f64;

            let mut ux = 0.0;
            let mut uy = 0.0;
            let mut uz = 0.0;
            let mut sx = 0.0;
            let mut sy = 0.0;
            let mut sz = 0.0;

            for p in 0..n {
                let dx = xj - field.px[p] as f64;
                let dy = yj - field.py[p] as f64;
                let dz = zj - field.pz[p] as f64;
                let r2 = dx * dx + dy * dy + dz * dz;

                let sig = field.sigma[p] as f64;
                let sig2 = sig * sig;
                let sig3 = sig2 * sig;
                let s = r2 / sig2;
                let base = s + 1.0;
                let base_sqrt = base.sqrt();
                let base25 = base * base * base_sqrt; // (s+1)^{5/2}
                let base35 = base25 * base; // (s+1)^{7/2}
                let k = (s + 2.5) / (sig3 * base25);
                let dkds = (-1.5 * s - 5.25) / (sig3 * base35);

                let apx = field.ax[p] as f64;
                let apy = field.ay[p] as f64;
                let apz = field.az[p] as f64;

                // alpha_src x d  (the velocity cross product).
                let cadx = apy * dz - apz * dy;
                let cady = apz * dx - apx * dz;
                let cadz = apx * dy - apy * dx;
                ux += k * cadx;
                uy += k * cady;
                uz += k * cadz;

                // Stretching: coef*(alpha_src x d) + K*(alpha_src x Gamma_tgt).
                let gdotd = gx * dx + gy * dy + gz * dz;
                let coef = 2.0 * dkds * gdotd / sig2;
                let cagx = apy * gz - apz * gy;
                let cagy = apz * gx - apx * gz;
                let cagz = apx * gy - apy * gx;
                sx += coef * cadx + k * cagx;
                sy += coef * cady + k * cagy;
                sz += coef * cadz + k * cagz;
            }

            (
                [ux * INV_4PI_F64, uy * INV_4PI_F64, uz * INV_4PI_F64],
                [sx * INV_4PI_F64, sy * INV_4PI_F64, sz * INV_4PI_F64],
            )
        })
        .collect()
}

/// rVPM strength/size derivatives for one particle given its strength `gamma`,
/// core `sigma`, and the stretching `s = (gamma . grad) u`.
/// Returns `(dGamma[3], dSigma)`.
#[inline]
fn rvpm_deriv(gamma: [f64; 3], sigma: f64, s: [f64; 3]) -> ([f64; 3], f64) {
    let g2 = gamma[0] * gamma[0] + gamma[1] * gamma[1] + gamma[2] * gamma[2];
    if g2 < GAMMA2_FLOOR {
        // Strengthless particle: no well-defined orientation. Take the raw
        // stretching for dGamma and leave sigma fixed.
        return (s, 0.0);
    }
    let sdotg = s[0] * gamma[0] + s[1] * gamma[1] + s[2] * gamma[2];
    // dGamma = S - (3/5)(S.Ghat)Ghat = S - (3/5)(S.G/|G|^2) G
    let c = 0.6 * sdotg / g2;
    let dgamma = [s[0] - c * gamma[0], s[1] - c * gamma[1], s[2] - c * gamma[2]];
    // dSigma = -(1/5) sigma (S.G) / |G|^2
    let dsigma = -0.2 * sigma * sdotg / g2;
    (dgamma, dsigma)
}

/// Advance the free wake one step of size `dt` with the reformulated VPM:
/// midpoint (RK2) integration of position (convection + freestream), vector
/// strength (stretching + reorientation), and core size (mass-conserving
/// evolution). Drop-in replacement for `vpm::advect_rk2` selected by
/// `WakeEngine::ReformulatedVpm`.
pub fn advect_rvpm(field: &mut ParticleField, freestream: [f32; 3], dt: f32) {
    let n = field.len();
    if n == 0 {
        return;
    }
    let dt64 = dt as f64;
    let half = 0.5 * dt64;
    let fs = [freestream[0] as f64, freestream[1] as f64, freestream[2] as f64];

    // ---- Stage 1: derivatives at the current state ----
    let k1 = vel_and_stretch(field);
    let mut mid = field.clone();
    for j in 0..n {
        let (u, s) = k1[j];
        let gamma = [field.ax[j] as f64, field.ay[j] as f64, field.az[j] as f64];
        let sigma = field.sigma[j] as f64;
        let (dg, dsig) = rvpm_deriv(gamma, sigma, s);
        mid.px[j] = (xf(field.px[j]) + half * (u[0] + fs[0])) as f32;
        mid.py[j] = (xf(field.py[j]) + half * (u[1] + fs[1])) as f32;
        mid.pz[j] = (xf(field.pz[j]) + half * (u[2] + fs[2])) as f32;
        mid.ax[j] = (gamma[0] + half * dg[0]) as f32;
        mid.ay[j] = (gamma[1] + half * dg[1]) as f32;
        mid.az[j] = (gamma[2] + half * dg[2]) as f32;
        mid.sigma[j] = (sigma + half * dsig).max(SIGMA_FLOOR) as f32;
    }

    // ---- Stage 2: derivatives at the midpoint, applied over the full step ----
    let k2 = vel_and_stretch(&mid);
    for j in 0..n {
        let (u, s) = k2[j];
        let gamma_mid = [mid.ax[j] as f64, mid.ay[j] as f64, mid.az[j] as f64];
        let sigma_mid = mid.sigma[j] as f64;
        let (dg, dsig) = rvpm_deriv(gamma_mid, sigma_mid, s);
        field.px[j] = (xf(field.px[j]) + dt64 * (u[0] + fs[0])) as f32;
        field.py[j] = (xf(field.py[j]) + dt64 * (u[1] + fs[1])) as f32;
        field.pz[j] = (xf(field.pz[j]) + dt64 * (u[2] + fs[2])) as f32;
        field.ax[j] = (field.ax[j] as f64 + dt64 * dg[0]) as f32;
        field.ay[j] = (field.ay[j] as f64 + dt64 * dg[1]) as f32;
        field.az[j] = (field.az[j] as f64 + dt64 * dg[2]) as f32;
        field.sigma[j] = (field.sigma[j] as f64 + dt64 * dsig).max(SIGMA_FLOOR) as f32;
    }
}

#[inline(always)]
fn xf(v: f32) -> f64 {
    v as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    fn one_particle() -> ParticleField {
        let mut f = ParticleField::new();
        f.push([1.0, 2.0, 3.0], [0.5, -0.2, 0.1], 0.3);
        f
    }

    /// A single particle has no other vorticity to stretch it: rVPM must reduce
    /// to pure convection by the freestream, leaving strength and core fixed.
    #[test]
    fn rvpm_single_particle_is_pure_convection() {
        let mut f = one_particle();
        let (a0, s0) = ([f.ax[0], f.ay[0], f.az[0]], f.sigma[0]);
        let p0 = [f.px[0], f.py[0], f.pz[0]];
        let fs = [2.0f32, -1.0, 0.5];
        let dt = 0.01f32;
        advect_rvpm(&mut f, fs, dt);
        // Position advanced by exactly freestream*dt (self-induction is zero).
        assert!((f.px[0] - (p0[0] + fs[0] * dt)).abs() < 1e-5);
        assert!((f.py[0] - (p0[1] + fs[1] * dt)).abs() < 1e-5);
        assert!((f.pz[0] - (p0[2] + fs[2] * dt)).abs() < 1e-5);
        // Strength and core unchanged (no stretching).
        assert!((f.ax[0] - a0[0]).abs() < 1e-6);
        assert!((f.ay[0] - a0[1]).abs() < 1e-6);
        assert!((f.az[0] - a0[2]).abs() < 1e-6);
        assert!((f.sigma[0] - s0).abs() < 1e-6);
    }

    /// A random cloud advances one rVPM step with all-finite state and positive
    /// cores (basic stability / no-NaN guard).
    #[test]
    fn rvpm_step_stays_finite() {
        let mut state = 0x1234_5678u32;
        let mut rng = || {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            (state >> 8) as f32 / (1u32 << 24) as f32 - 0.5
        };
        let mut f = ParticleField::new();
        for _ in 0..200 {
            f.push(
                [rng() * 2.0, rng() * 2.0, rng() * 2.0],
                [rng(), rng(), rng()],
                0.15 + 0.1 * (rng() + 0.5).abs(),
            );
        }
        advect_rvpm(&mut f, [1.0, 0.0, 0.0], 0.02);
        for i in 0..f.len() {
            assert!(f.px[i].is_finite() && f.py[i].is_finite() && f.pz[i].is_finite());
            assert!(f.ax[i].is_finite() && f.ay[i].is_finite() && f.az[i].is_finite());
            assert!(f.sigma[i].is_finite() && f.sigma[i] > 0.0);
        }
    }

    /// The stretching evaluator must produce the same induced velocity as the
    /// classic engine (the `u` half of `(u, S)` is the identical Biot-Savart
    /// sum). Guards against a kernel-transcription error in the rVPM path.
    #[test]
    fn rvpm_velocity_matches_classic() {
        let mut state = 0x9E37_79B9u32;
        let mut rng = || {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            (state >> 8) as f32 / (1u32 << 24) as f32 - 0.5
        };
        let mut f = ParticleField::new();
        for _ in 0..300 {
            f.push(
                [rng() * 4.0, rng() * 4.0, rng() * 4.0],
                [rng(), rng(), rng()],
                0.2 + 0.1 * (rng() + 0.5).abs(),
            );
        }
        let classic = crate::vpm::common::induced_velocities_ref(&f);
        let rvpm = vel_and_stretch(&f);
        let peak = classic
            .iter()
            .map(|u| (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt())
            .fold(0.0f64, f64::max);
        for j in 0..f.len() {
            for k in 0..3 {
                assert!(
                    (rvpm[j].0[k] - classic[j][k]).abs() <= 1e-9 * (1.0 + peak),
                    "velocity mismatch at {j},{k}: rvpm {} vs classic {}",
                    rvpm[j].0[k],
                    classic[j][k]
                );
            }
        }
    }
}
