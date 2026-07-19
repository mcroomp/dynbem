// Level 1 BEM: helicopter momentum quadratic + wind-turbine windmill (Brent).
// See ../CLAUDE.md, dynbem/bem.py for the full physics.

use std::f64::consts::PI;

use crate::aero_io::{AeroResult, RotorInputs};
use crate::aero_model::{AeroModel, RotorStateExt};
use crate::bem_common::{
    apply_flap_reduction, assemble_result, build_psi_trig_table, kinematics, ElementCtx, PsiKernel,
    RadialGrid, SweepCtx,
};
use crate::common::{EPS_DENOM, EPS_OMEGA_R, MIN_LOSS_FACTOR};
use crate::cyclic::cyclic_coeffs;
use crate::polar::Polar;
use crate::rotor_definition::{PitchActuation, RotorDefinition};
use crate::servoflap::{solve_feathering, FeatheringState};

const MAX_BEM_ITER: usize = 60;
const BEM_TOL: f64 = 1e-7;

/// Minimum |lambda_climb| (tip-referenced axial inflow ratio, v_climb /
/// (omega*radius_m)) required before the dedicated windmill Brent solver is
/// even attempted. Below this, the local axial flow is too small relative
/// to blade speed for the windmill solver's assumptions to hold reliably
/// (its `lam_local = omega*r/u_up` term grows large and the solved root
/// becomes numerically noisy/inconsistent from one element or azimuth to
/// the next -- see tests/test_bem_windmill_boundary.py in the
/// windpower-repo history). Typical hover/climb induced inflow ratios are
/// O(0.02-0.08), so this threshold sits at the low end of that range: well
/// above numerical noise, comfortably below where genuine windmill-brake
/// descent physics take over. Below this threshold, all
/// elements/azimuths uniformly fall back to solve_bem_element (helicopter
/// momentum quadratic), which is continuous through lambda_climb == 0.
const MIN_LAMBDA_CLIMB_WINDMILL: f64 = 0.02;

#[derive(Clone, Debug, Default)]
pub struct QuasiStaticRotorState;

// ---------------------------------------------------------------------------
// Prandtl tip / hub losses
// ---------------------------------------------------------------------------

#[inline]
pub fn prandtl_tip_loss(n_blades: usize, x: f64, phi_rad: f64) -> f64 {
    prandtl_tip_loss_from_sin_abs(n_blades, x, phi_rad.sin().abs())
}

#[inline]
pub fn prandtl_tip_loss_from_sin_abs(n_blades: usize, x: f64, sin_phi_abs: f64) -> f64 {
    if sin_phi_abs < EPS_DENOM || x >= 1.0 {
        return 1.0;
    }
    let f = (n_blades as f64) / 2.0 * (1.0 - x) / (x * sin_phi_abs);
    (2.0 / PI) * (1.0_f64.min((-f).exp())).acos()
}

#[inline]
pub fn prandtl_hub_loss(n_blades: usize, x: f64, x_hub: f64, phi_rad: f64) -> f64 {
    prandtl_hub_loss_from_sin_abs(n_blades, x, x_hub, phi_rad.sin().abs())
}

#[inline]
pub fn prandtl_hub_loss_from_sin_abs(n_blades: usize, x: f64, x_hub: f64, sin_phi_abs: f64) -> f64 {
    if sin_phi_abs < EPS_DENOM || x <= x_hub || x_hub <= 0.0 {
        return 1.0;
    }
    let f = (n_blades as f64) / 2.0 * (x - x_hub) / (x_hub * sin_phi_abs);
    (2.0 / PI) * (1.0_f64.min((-f).exp())).acos()
}

// ---------------------------------------------------------------------------
// BEM element result (internal)
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, Default)]
pub struct BEMElementResult {
    pub d_t: f64,               // thrust contribution [N]
    pub d_q: f64,               // torque contribution [N.m]
    pub lambda_r: f64,          // converged axial inflow ratio v_a / (Omega * R)
    pub a_prime: f64,           // converged tangential induction factor
    pub momentum_residual: f64, // |4*F*lambda_r*(lambda_r - lambda_climb) - sigma_r*cn*(lambda_r^2 + x^2)|
}

// ---------------------------------------------------------------------------
// BEM element geometry (constant-per-element cached values)
// ---------------------------------------------------------------------------

/// Precomputed geometry and constants for one blade element.
/// Holds the radial/geometric quantities that remain constant across multiple
/// solver calls at the same station. Only the time-varying aerodynamic
/// quantities (v_climb, collective_rad per azimuth, v_t_extra) are passed
/// separately to the solver functions.
pub struct BEMElementGeometry<'a, P: Polar> {
    // Blade/rotor geometry
    pub r: f64,
    pub dr: f64,
    pub chord: f64,
    pub twist_rad: f64,
    pub omega: f64,
    pub rho: f64,
    pub n_blades: usize,
    pub radius_m: f64,
    pub polar: &'a P,
    pub use_tip_loss: bool,
    pub root_cutout_m: f64,
    // Precomputed constants
    pub omega_r: f64,
    pub inv_omega_r: f64,
    pub inv_r: f64,
    pub x: f64,
    pub x_hub: f64,
    pub sigma_r: f64,
}

impl<'a, P: Polar> BEMElementGeometry<'a, P> {
    pub fn new(
        r: f64,
        dr: f64,
        chord: f64,
        twist_rad: f64,
        omega: f64,
        rho: f64,
        n_blades: usize,
        radius_m: f64,
        polar: &'a P,
        use_tip_loss: bool,
        root_cutout_m: f64,
    ) -> Self {
        let omega_r = omega * radius_m;
        let inv_omega_r = if omega_r > EPS_OMEGA_R {
            1.0 / omega_r
        } else {
            1.0
        };
        let inv_r = if r > 0.0 { 1.0 / r } else { 0.0 };
        let inv_radius_m = if radius_m > 0.0 { 1.0 / radius_m } else { 0.0 };
        let x = r * inv_radius_m;
        let x_hub = root_cutout_m * inv_radius_m;
        let sigma_r = (n_blades as f64) * chord * inv_r / (2.0 * PI);

        Self {
            r,
            dr,
            chord,
            twist_rad,
            omega,
            rho,
            n_blades,
            radius_m,
            polar,
            use_tip_loss,
            root_cutout_m,
            omega_r,
            inv_omega_r,
            inv_r,
            x,
            x_hub,
            sigma_r,
        }
    }
}

// ---------------------------------------------------------------------------
// Helicopter momentum quadratic
// ---------------------------------------------------------------------------

/// Helicopter momentum-BEM solver for one annulus.
///
/// Fixed-point iteration on (lambda_r, a_prime) with 50% under-relaxation;
/// the converged root of the quadratic is always the climb/hover branch
/// (see the seeding comment in the function body for why). Reverse-flow
/// region (v_t < 0) breaks out early and returns zero forces -- caller
/// is responsible for the surrounding psi-loop's reverse-flow skip.
pub fn solve_bem_element<P: Polar>(
    geom: &BEMElementGeometry<P>,
    collective_rad: f64,
    v_climb: f64,
    v_t_extra: f64,
) -> BEMElementResult {
    if geom.omega_r < EPS_OMEGA_R {
        return BEMElementResult::default();
    }
    let theta = collective_rad + geom.twist_rad;
    let lambda_climb = v_climb * geom.inv_omega_r;

    // Seed from climb branch for hover/climb and for the near-hover descent
    // band (|v_climb| < MIN_LAMBDA_CLIMB_WINDMILL * tip_speed). In that band
    // the windmill solver has already declined (its Brent bracket doesn't
    // exist), and hover is the v_climb -> 0 limit of the *climb* branch --
    // seeding from the descent branch there previously caused a large
    // spurious discontinuity right at lambda_climb == 0.
    //
    // For genuine deep descent (outside the windmill threshold) this function
    // may be called as a fallback when the windmill solver found no bracket
    // for a particular azimuth/element; in that regime the descent-branch
    // root is physically correct and is used instead.
    let near_hover = v_climb >= -MIN_LAMBDA_CLIMB_WINDMILL * geom.omega * geom.radius_m;
    let mut lambda_r = if near_hover || lambda_climb >= 0.0 {
        (lambda_climb + 0.03).max(0.02)
    } else {
        (lambda_climb * 0.85).min(-0.02)
    };
    let mut a_prime: f64 = 0.0;

    for _ in 0..MAX_BEM_ITER {
        let v_a = lambda_r * geom.omega_r;
        let v_t = geom.omega * geom.r * (1.0 + a_prime) + v_t_extra;
        if v_t < EPS_DENOM {
            break;
        }

        let phi = v_a.atan2(v_t);
        let alpha = theta - phi;
        let (cl, cd) = geom.polar.cl_cd(alpha);

        let cos_p = phi.cos();
        let sin_p = phi.sin();
        let sin_phi_abs = sin_p.abs();

        let f_loss = if geom.use_tip_loss {
            (prandtl_tip_loss_from_sin_abs(geom.n_blades, geom.x, sin_phi_abs)
                * prandtl_hub_loss_from_sin_abs(geom.n_blades, geom.x, geom.x_hub, sin_phi_abs))
            .max(MIN_LOSS_FACTOR)
        } else {
            1.0
        };
        let cn = cl * cos_p - cd * sin_p;
        let ct = cl * sin_p + cd * cos_p;

        let k = geom.sigma_r * cn / (4.0 * f_loss);

        let lambda_r_new = if (k - 1.0).abs() > 1e-6 {
            let disc =
                (lambda_climb * lambda_climb - 4.0 * (k - 1.0) * k * geom.x * geom.x).max(0.0);
            let sq = disc.sqrt();
            let denom = 2.0 * (k - 1.0);
            let r1 = (-lambda_climb + sq) / denom;
            let r2 = (-lambda_climb - sq) / denom;
            // Mirror the seeding choice: climb branch near hover, descent
            // branch for genuine deep descent.
            if near_hover || lambda_climb >= 0.0 {
                if r2 > 0.0 {
                    r2
                } else {
                    r1
                }
            } else {
                if r1 < 0.0 {
                    r1
                } else {
                    r2
                }
            }
        } else if lambda_climb.abs() > 1e-8 {
            -k * geom.x * geom.x / lambda_climb
        } else {
            geom.x * k.max(0.0).sqrt()
        };
        let lambda_r_new = lambda_r_new.clamp(-2.0, 2.0);

        let sc = sin_p * cos_p;
        let a_prime_new = if sc.abs() > 1e-8 && ct.abs() > 1e-10 {
            let ap_denom = 4.0 * f_loss * sc / (geom.sigma_r * ct) - 1.0;
            let v = if ap_denom.abs() > 1e-8 {
                1.0 / ap_denom
            } else {
                0.0
            };
            v.clamp(-0.5, 0.5)
        } else {
            0.0
        };

        let converged =
            (lambda_r_new - lambda_r).abs() < BEM_TOL && (a_prime_new - a_prime).abs() < BEM_TOL;
        lambda_r = 0.5 * lambda_r + 0.5 * lambda_r_new;
        a_prime = 0.5 * a_prime + 0.5 * a_prime_new;
        if converged {
            break;
        }
    }

    let v_a = lambda_r * geom.omega_r;
    let v_t = geom.omega * geom.r * (1.0 + a_prime) + v_t_extra;
    let v_rel_sq = v_a * v_a + v_t * v_t;
    let phi = v_a.atan2(v_t);
    let alpha = theta - phi;
    let (cl, cd) = geom.polar.cl_cd(alpha);
    let cos_p = phi.cos();
    let sin_p = phi.sin();
    let sin_phi_abs = sin_p.abs();
    let cn = cl * cos_p - cd * sin_p;
    let ct = cl * sin_p + cd * cos_p;
    let q = 0.5 * geom.rho * v_rel_sq * geom.chord * geom.dr * (geom.n_blades as f64);

    // Momentum-balance residual at the converged state:
    //   |4*F*lambda_r*(lambda_r - lambda_climb) - sigma_r*cn*(lambda_r^2 + x^2)|
    let f_loss = if geom.use_tip_loss {
        (prandtl_tip_loss_from_sin_abs(geom.n_blades, geom.x, sin_phi_abs)
            * prandtl_hub_loss_from_sin_abs(geom.n_blades, geom.x, geom.x_hub, sin_phi_abs))
        .max(MIN_LOSS_FACTOR)
    } else {
        1.0
    };
    let momentum_residual = (4.0 * f_loss * lambda_r * (lambda_r - lambda_climb)
        - geom.sigma_r * cn * (lambda_r * lambda_r + geom.x * geom.x))
        .abs();

    BEMElementResult {
        d_t: q * cn,
        d_q: q * ct * geom.r,
        lambda_r,
        a_prime,
        momentum_residual,
    }
}

// ---------------------------------------------------------------------------
// Brent's method (van Wijngaarden-Dekker-Brent) for the windmill solver
// ---------------------------------------------------------------------------

/// Brent's method root finder on a sign-changing bracket.
///
/// Returns None when f(a)*f(b) > 0 (no sign change), or when the residual
/// goes non-finite mid-iteration. Otherwise converges in at most max_iter
/// steps; returns the current best estimate after that.
fn brentq<F: FnMut(f64) -> f64>(
    mut f: F,
    mut a: f64,
    mut b: f64,
    xtol: f64,
    max_iter: usize,
) -> Option<f64> {
    let mut fa = f(a);
    let mut fb = f(b);
    if !fa.is_finite() || !fb.is_finite() {
        return None;
    }
    if fa * fb > 0.0 {
        return None;
    }
    if fa.abs() < fb.abs() {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut fa, &mut fb);
    }
    let mut c = a;
    let mut fc = fa;
    let mut mflag = true;
    let mut d = 0.0f64;

    for _ in 0..max_iter {
        if fb == 0.0 || (b - a).abs() < xtol {
            return Some(b);
        }

        let s = if fa != fc && fb != fc {
            // inverse quadratic interpolation
            a * fb * fc / ((fa - fb) * (fa - fc))
                + b * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb))
        } else {
            // secant
            b - fb * (b - a) / (fb - fa)
        };

        let bound_lo = (3.0 * a + b) / 4.0;
        let cond1 = (s - bound_lo) * (s - b) > 0.0;
        let cond2 = mflag && (s - b).abs() >= (b - c).abs() / 2.0;
        let cond3 = !mflag && (s - b).abs() >= (c - d).abs() / 2.0;
        let cond4 = mflag && (b - c).abs() < xtol;
        let cond5 = !mflag && (c - d).abs() < xtol;
        let use_bisect = cond1 || cond2 || cond3 || cond4 || cond5;
        let s = if use_bisect { (a + b) / 2.0 } else { s };
        mflag = use_bisect;

        let fs = f(s);
        if !fs.is_finite() {
            return None;
        }
        d = c;
        c = b;
        fc = fb;
        if fa * fs < 0.0 {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        if fa.abs() < fb.abs() {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut fa, &mut fb);
        }
    }
    Some(b)
}

// ---------------------------------------------------------------------------
// Wind-turbine windmill solver (Ning 2014 Brent-on-phi)
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct WindmillInductions {
    a: f64,
    ap: f64,
    cn: f64,
    ct: f64,
}

fn induction_at_phi<P: Polar>(
    phi: f64,
    theta: f64,
    lam_local: f64,
    geom: &BEMElementGeometry<P>,
) -> Option<WindmillInductions> {
    let sin_phi = phi.sin();
    let cos_phi = phi.cos();
    if sin_phi.abs() < 1e-12 {
        return None;
    }
    let alpha = theta - phi;
    let (cl, cd) = geom.polar.cl_cd(alpha);
    let cn = cl * cos_phi - cd * sin_phi;
    let ct = cl * sin_phi + cd * cos_phi;
    let sin_phi_abs = sin_phi.abs();
    let f_loss = if geom.use_tip_loss {
        prandtl_tip_loss_from_sin_abs(geom.n_blades, geom.x, sin_phi_abs)
            * prandtl_hub_loss_from_sin_abs(geom.n_blades, geom.x, geom.x_hub, sin_phi_abs)
    } else {
        1.0
    };
    let f_loss = f_loss.max(MIN_LOSS_FACTOR);
    let sin2 = sin_phi * sin_phi;
    if cn <= EPS_DENOM {
        return None;
    }
    let k_axial = geom.sigma_r * cn / (4.0 * f_loss * sin2);
    let mut a = k_axial / (1.0 + k_axial);
    if a > 0.4 {
        // Buhl quadratic in the turbulent-wake state
        let k = geom.sigma_r * cn / sin2;
        let aa = 50.0 / 9.0 - 4.0 * f_loss - k;
        let bb = 4.0 * f_loss - 40.0 / 9.0 + 2.0 * k;
        let cc = 8.0 / 9.0 - k;
        let disc = bb * bb - 4.0 * aa * cc;
        if disc < 0.0 || aa.abs() < 1e-12 {
            return None;
        }
        let sq = disc.sqrt();
        let r1 = (-bb - sq) / (2.0 * aa);
        let r2 = (-bb + sq) / (2.0 * aa);
        let cands = [r1, r2];
        let mut found: Option<f64> = None;
        for &cand in &cands {
            if (0.4..=1.0).contains(&cand) {
                found = Some(match found {
                    Some(p) => p.min(cand),
                    None => cand,
                });
            }
        }
        let Some(a_buhl) = found else { return None };
        a = a_buhl;
    }
    let sc = sin_phi * cos_phi;
    let ap = if ct.abs() > EPS_DENOM && sc.abs() > EPS_DENOM {
        let k_tan = geom.sigma_r * ct / (4.0 * f_loss * sc);
        let v = if (1.0 - k_tan).abs() > EPS_DENOM {
            k_tan / (1.0 - k_tan)
        } else {
            0.0
        };
        v.clamp(-0.5, 0.5)
    } else {
        0.0
    };
    let _ = lam_local;
    Some(WindmillInductions { a, ap, cn, ct })
}

/// Wind-turbine BEM solver for one annulus.
///
/// Ning 2014 reformulation: residual-on-phi with Brent over (-pi/2, 0),
/// plus the Buhl quadratic for the turbulent-wake state (a > 0.4) where
/// classical momentum theory breaks down. Returns None when no sign-change
/// bracket exists or the iteration leaves the valid windmill regime
/// (0 < a < 1 and Cn > 0) -- caller falls back to the helicopter
/// quadratic in those cases.
///
/// `v_t_extra` accounts for in-plane wind at this azimuth (same meaning as
/// in `solve_bem_element`). The geometric residual uses
/// `(1+ap)*lam_local + v_t_extra/u_up` for the tangential term so that
/// the converged phi reflects the full tangential velocity. Forces are
/// computed with `v_t = (1+ap)*omega*r + v_t_extra`.
fn solve_bem_element_windmill<P: Polar>(
    geom: &BEMElementGeometry<P>,
    collective_rad: f64,
    v_climb: f64,
    v_t_extra: f64,
) -> Option<BEMElementResult> {
    if geom.omega_r < EPS_OMEGA_R {
        return None;
    }
    // Rotor-level normalized threshold (see MIN_LAMBDA_CLIMB_WINDMILL doc
    // comment): reject when the *tip*-referenced inflow ratio is too small
    // for the windmill solver's assumptions to hold. Deliberately uses the
    // tip speed (geom.omega * geom.radius_m), not this element's own local
    // omega*r, so every element/azimuth in a given call switches branch
    // together -- gating per-element would let inboard and outboard
    // stations disagree on which solver to use at the same rotor-level
    // v_climb, producing element-to-element inconsistency across the
    // integrated blade.
    if v_climb >= -MIN_LAMBDA_CLIMB_WINDMILL * geom.omega * geom.radius_m {
        return None;
    }
    let u_up = -v_climb;
    let inv_u_up = 1.0 / u_up;
    let theta = collective_rad + geom.twist_rad;
    let lam_local = geom.omega * geom.r * inv_u_up;
    // Tangential wind contribution normalised by u_up (for the geometric residual).
    let lam_extra = v_t_extra * inv_u_up;

    let residual = |phi: f64| -> f64 {
        match induction_at_phi(phi, theta, lam_local, geom) {
            None => 1e3,
            // Geometric consistency: sin(phi)*v_t/u_up + cos(phi)*(1-a) = 0
            // v_t/u_up = (1+ap)*lam_local + lam_extra
            Some(ind) => {
                phi.sin() * ((1.0 + ind.ap) * lam_local + lam_extra) + phi.cos() * (1.0 - ind.a)
            }
        }
    };

    let phi_lo = -0.5 * PI + 1e-4;
    let phi_hi = -1e-4;
    let r_lo = residual(phi_lo);
    let r_hi = residual(phi_hi);
    if !(r_lo.is_finite() && r_hi.is_finite()) || r_lo * r_hi >= 0.0 {
        return None;
    }
    let phi_star = brentq(residual, phi_lo, phi_hi, 1e-8, 80)?;
    let ind = induction_at_phi(phi_star, theta, lam_local, geom)?;
    if ind.cn <= 0.0 || !(0.0..=1.0).contains(&ind.a) {
        return None;
    }
    let v_a = -(1.0 - ind.a) * u_up;
    let v_t = (1.0 + ind.ap) * geom.omega * geom.r + v_t_extra;
    let v_rel_sq = v_a * v_a + v_t * v_t;
    let q = 0.5 * geom.rho * v_rel_sq * geom.chord * geom.dr * (geom.n_blades as f64);
    // Windmill solver works in (a, a') space, not (lambda_r, a_prime); reconstruct
    // lambda_r consistent with the rest of the pipeline (axial inflow ratio).
    let lambda_r = v_a * geom.inv_omega_r;
    Some(BEMElementResult {
        d_t: q * ind.cn,
        d_q: q * ind.ct * geom.r,
        lambda_r,
        a_prime: ind.ap,
        momentum_residual: 0.0, // converged via Brent-on-phi, residual is on phi not lambda_r
    })
}

// ---------------------------------------------------------------------------
// BEM PsiKernel: overrides element() entirely so the shared psi-loop runs
// solve_bem_element (iterative quadratic) per (psi, r) instead of the
// prescribed-inflow path PP and Oye use.
// ---------------------------------------------------------------------------

/// BEM-specific state for the shared psi-loop.
struct BemKernel {
    v_climb: f64,
    r_tip: f64,
    root_cutout_m: f64,
    use_tip_loss: bool,
}

impl PsiKernel for BemKernel {
    #[inline(always)]
    fn element<P: Polar>(&mut self, sweep: &SweepCtx<'_, P>, ctx: &ElementCtx) -> (f64, f64) {
        let v_t_extra = ctx.v_t - sweep.omega * ctx.r;

        // Construct the geometry once, reuse for both windmill and helicopter paths
        let geom = BEMElementGeometry::new(
            ctx.r,
            ctx.dr,
            ctx.chord,
            ctx.twist,
            sweep.omega,
            sweep.rho,
            sweep.n_b,
            self.r_tip,
            sweep.polar,
            self.use_tip_loss,
            self.root_cutout_m,
        );

        // When v_climb < 0 (upflow through disk), always try the Ning 2014
        // windmill Brent solver first.  The helicopter momentum quadratic
        // cannot find the windmill-brake root in this regime: its two roots
        // for lambda_climb < 0 land outside (lambda_climb, 0) and the chosen
        // root is increasingly wrong as |lambda_climb| shrinks.
        //
        // The bracket-existence check inside solve_bem_element_windmill is the
        // natural physics-based filter:
        //
        // * Windmill / steep-descent (negative collective, theta < 0):
        //   at phi_hi = -1e-4, alpha = theta - phi_hi ≈ theta < 0 → CL < 0 →
        //   cn < 0 → induction_at_phi returns None → residual = 1e3 (positive).
        //   phi_lo gives large negative.  Sign change → bracket found → solver
        //   returns the correct windmill root.
        //
        // * Autorotation / forward-flight descent (positive collective):
        //   at phi_hi, alpha ≈ theta > 0 → cn > 0 → a → 1 → (1-a) → 0 →
        //   residual ≈ small negative.  phi_lo also negative.  No sign change
        //   → windmill solver returns None → falls back to helicopter quadratic.
        if self.v_climb < -EPS_DENOM {
            if let Some(wm) =
                solve_bem_element_windmill(&geom, ctx.col_psi, self.v_climb, v_t_extra)
            {
                return (wm.d_t, wm.d_q);
            }
        }
        let elem = solve_bem_element(&geom, ctx.col_psi, self.v_climb, v_t_extra);
        (elem.d_t, elem.d_q)
    }
}

// ---------------------------------------------------------------------------
// QuasiStaticBEM: pyclass holding cached radial grid + polar
// ---------------------------------------------------------------------------

#[derive(Clone)]
pub struct QuasiStaticBEM<P: Polar> {
    pub defn: RotorDefinition,
    pub n_psi_elements: usize,
    pub psi_trig: Vec<(f64, f64)>,
    pub polar: P,
    pub grid: RadialGrid,
}

impl<P: Polar + Clone> QuasiStaticBEM<P> {
    pub fn build(defn: RotorDefinition, n_psi_elements: usize, polar: P) -> Self {
        let grid = RadialGrid::from_blade(&defn.blade);
        let psi_trig = build_psi_trig_table(n_psi_elements);
        Self {
            defn,
            n_psi_elements,
            psi_trig,
            polar,
            grid,
        }
    }
}

impl RotorStateExt for QuasiStaticRotorState {
    fn get_inflow(&self) -> Vec<f64> {
        Vec::new()
    }
    fn set_inflow(&mut self, arr: Vec<f64>) {
        debug_assert!(arr.is_empty());
    }
}

impl<P: Polar + Clone> AeroModel for QuasiStaticBEM<P> {
    type State = QuasiStaticRotorState;

    fn initial_state(&self) -> Self::State {
        QuasiStaticRotorState::default()
    }

    // inflow_taus: trait default (all-infinity) is correct for the
    // quasi-static BEM model; no override needed.

    fn compute_forces(
        &self,
        inputs: &RotorInputs,
        _state: &QuasiStaticRotorState,
    ) -> (AeroResult, QuasiStaticRotorState) {
        let blade = &self.defn.blade;
        let omega = inputs.omega_rad_s;
        let rho = inputs.rho_kg_m3;
        let r_tip = blade.radius_m;
        let n_blades = blade.n_blades;
        let use_tip_loss = self.defn.blade.tip_loss;
        let grid = &self.grid;

        let kin = kinematics(inputs, omega, r_tip);
        let omega_r = kin.omega_r;
        let hub_axis = kin.hub_axis;
        let v_climb = kin.v_climb;
        let v_inplane_hub = kin.v_inplane_hub;
        let mu = kin.mu;

        let n = blade.n_elements;
        let r_arr = &grid.r_mid;
        let chord = &grid.chord;
        let twist = &grid.twist_rad;
        let dr = grid.dr;

        // Cyclic pitch -> theta_1c, theta_1s
        let gains = self.defn.control_gains();
        let (theta_1c, theta_1s) = cyclic_coeffs(inputs.tilt_lon, inputs.tilt_lat, gains);
        let has_cyclic = theta_1c.abs() + theta_1s.abs() > 1e-12;

        // Quasi-static feathering solve (pre-pass)
        let feathering_state = match &self.defn.pitch_actuation {
            PitchActuation::DirectMechanical => FeatheringState::RIGID,
            PitchActuation::ServoFlap(act) => solve_feathering(
                act,
                inputs.collective_rad,
                theta_1c,
                theta_1s,
                rho,
                omega,
                mu,
                r_tip,
                blade.chord_m,
                self.defn.airfoil.CL_alpha_per_rad,
            ),
        };
        let has_feathering_cyclic =
            feathering_state.delta_theta_1c.abs() + feathering_state.delta_theta_1s.abs() > 1e-12;
        let servo_mode = self.defn.is_servoflap();
        let loop_collective = if servo_mode {
            feathering_state.delta_theta_0
        } else {
            inputs.collective_rad
        };
        let (loop_theta_1c, loop_theta_1s) = if servo_mode {
            (
                feathering_state.delta_theta_1c,
                feathering_state.delta_theta_1s,
            )
        } else {
            (theta_1c, theta_1s)
        };

        let mut t_total: f64 = 0.0;
        let mut q_total: f64 = 0.0;
        let mut mx_hub: f64 = 0.0;
        let mut my_hub: f64 = 0.0;
        let mut fx_hub: f64 = 0.0;
        let mut fy_hub: f64 = 0.0;

        if (mu > 0.01 || has_cyclic || has_feathering_cyclic) && omega > 1.0 {
            let mut kernel = BemKernel {
                v_climb,
                r_tip,
                root_cutout_m: blade.root_cutout_m,
                use_tip_loss,
            };
            let sweep = SweepCtx {
                grid,
                polar: &self.polar,
                col: loop_collective,
                omega,
                omega_r,
                rho,
                n_b: n_blades,
                n_psi: self.n_psi_elements,
                n_psi_inv: 1.0 / (self.n_psi_elements as f64),
                psi_trig: &self.psi_trig,
                v_in_hub_x: v_inplane_hub[0],
                v_in_hub_y: v_inplane_hub[1],
                theta_1c: loop_theta_1c,
                theta_1s: loop_theta_1s,
            };
            let (t, q, mx, my, fx, fy) = sweep.run(&mut kernel);
            t_total = t;
            q_total = q;
            mx_hub = mx;
            my_hub = my;
            fx_hub = fx;
            fy_hub = fy;
        } else {
            // Axial: try wind-turbine windmill solver first when v_climb < 0,
            // fall back to helicopter quadratic per element.
            // (Fx_hub/Fy_hub stay 0: azimuth-symmetric loading has no net
            // in-plane force.)
            for i_r in 0..n {
                let r = r_arr[i_r];
                let geom = BEMElementGeometry::new(
                    r,
                    dr,
                    chord[i_r],
                    twist[i_r],
                    omega,
                    rho,
                    n_blades,
                    r_tip,
                    &self.polar,
                    use_tip_loss,
                    blade.root_cutout_m,
                );
                let mut elem: Option<BEMElementResult> = None;
                if v_climb < -EPS_DENOM {
                    elem = solve_bem_element_windmill(&geom, loop_collective, v_climb, 0.0);
                }
                let elem =
                    elem.unwrap_or_else(|| solve_bem_element(&geom, loop_collective, v_climb, 0.0));
                t_total += elem.d_t;
                q_total += elem.d_q;
            }
        }

        let (mx_out, my_out) =
            apply_flap_reduction(mx_hub, my_hub, self.defn.flap.as_ref(), inputs.omega_rad_s);
        let result = assemble_result(
            t_total,
            q_total,
            mx_out,
            my_out,
            fx_hub,
            fy_hub,
            hub_axis,
            &inputs.R_hub,
        );

        let derivative = QuasiStaticRotorState;
        // suppress unused warning (use_tip_loss read inside windmill helper).
        let _ = use_tip_loss;
        (result, derivative)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aero_io::{Mat3, Vec3};
    use crate::rotor_definition::{ControlProperties, LinearPolarParameters, PitchActuation};

    // RAWES clean IC attitude (NED, +East wind = 10 m/s).
    // body_z = R_hub[:,2] = [0, -0.9042, 0.4272]
    fn r_rawes_ic() -> Mat3 {
        Mat3([
            [0.0, -1.0, 0.0],
            [0.42720594325829603, 0.0, -0.9041543463617201],
            [0.9041543463617201, 0.0, 0.4272059432582967],
        ])
    }

    fn beaupoil_rotor() -> RotorDefinition {
        use crate::rotor_definition::BladeGeometry;
        RotorDefinition {
            blade: BladeGeometry {
                n_blades: 4,
                radius_m: 2.5,
                root_cutout_m: 0.5,
                chord_m: 0.20,
                twist_deg: 0.0,
                n_elements: 10,
                tip_loss: true,
                r_stations_m: Vec::new(),
                chord_stations_m: Vec::new(),
                twist_stations_deg: Vec::new(),
            },
            airfoil: LinearPolarParameters {
                CL0: 0.393,
                CL_alpha_per_rad: 5.79,
                CD0: 0.0079,
                alpha_stall_deg: 13.0,
            },
            control: Some(ControlProperties {
                swashplate_pitch_gain_rad: 0.3,
                swashplate_phase_deg: Some(0.0),
            }),
            pitch_actuation: PitchActuation::DirectMechanical,
            flap: None,
            name: "beaupoil_2026".to_string(),
            description: String::new(),
        }
    }

    fn rawes_ic_inputs(omega: f64, tilt_lon: f64, tilt_lat: f64) -> RotorInputs {
        RotorInputs {
            collective_rad: -0.18,
            tilt_lon,
            tilt_lat,
            R_hub: r_rawes_ic(),
            v_hub_world: Vec3::zero(),
            wind_world: Vec3::new(0.0, 10.0, 0.0),
            omega_rad_s: omega,
            rho_kg_m3: 1.225,
        }
    }

    fn f_dot_body_z(result: &AeroResult, r_hub: &Mat3) -> f64 {
        let body_z = Vec3::new(r_hub.0[0][2], r_hub.0[1][2], r_hub.0[2][2]);
        result.F_world.dot(body_z)
    }

    fn qs_beaupoil() -> QuasiStaticBEM<crate::polar::LinearPolar> {
        let defn = beaupoil_rotor();
        let polar = crate::polar::LinearPolar::from_properties(&defn.airfoil);
        QuasiStaticBEM::build(defn, 36, polar)
    }

    /// QS BEM at the RAWES IC: force must oppose body-Z (thrust direction).
    #[test]
    fn rawes_ic_qs_force_opposes_body_z() {
        let model = qs_beaupoil();
        let inputs = rawes_ic_inputs(53.161687, 0.0, 0.0);
        let (result, _) = model.compute_forces(&inputs, &model.initial_state());
        let fdz = f_dot_body_z(&result, &inputs.R_hub);
        assert!(
            fdz < 0.0,
            "QS RAWES IC: F dot body_z should be negative, got {fdz:.3}"
        );
        assert!(
            result.F_world.0[1] > 0.0,
            "QS RAWES IC: F_east should be positive (downwind)"
        );
        assert!(
            result.F_world.0[2] < 0.0,
            "QS RAWES IC: F_up should be positive (-Z)"
        );
    }

    /// Same test with trim cyclic applied.
    #[test]
    fn rawes_ic_qs_with_trim_cyclic_opposes_body_z() {
        let model = qs_beaupoil();
        let inputs = rawes_ic_inputs(53.161687, 0.0, 0.022616);
        let (result, _) = model.compute_forces(&inputs, &model.initial_state());
        let fdz = f_dot_body_z(&result, &inputs.R_hub);
        assert!(
            fdz < 0.0,
            "QS RAWES IC trim cyclic: F dot body_z should be negative, got {fdz:.3}"
        );
    }

    /// QS BEM on logged RAWES flight row 122 (oblique descent + crosswind).
    /// Force must still oppose body-Z at this highly oblique operating point.
    #[test]
    fn rawes_row122_force_opposes_body_z() {
        use crate::aero_model::AeroModel;
        let defn = beaupoil_rotor();
        let polar = crate::polar::LinearPolar::from_properties(&defn.airfoil);
        let model = QuasiStaticBEM::build(defn, 36, polar);
        let inputs = RotorInputs {
            collective_rad: -0.18146020320832462,
            tilt_lon: 0.012614825392536453,
            tilt_lat: 0.035447368174067954,
            R_hub: Mat3([
                [
                    -0.007232662001129507,
                    -0.9995813722829041,
                    0.02801372495419454,
                ],
                [0.684832857996051, -0.025365018889292906, -0.728258589019448],
                [0.7286642884136498, 0.013917471093409295, 0.684729624547885],
            ]),
            v_hub_world: Vec3::new(
                -0.5280840184646872,
                -0.17333482520461627,
                -0.4766214883089747,
            ),
            wind_world: Vec3::new(0.0, 10.0, 0.0),
            rho_kg_m3: 1.225,
            omega_rad_s: 37.02311435435481,
        };
        let (result, _) = model.compute_forces(&inputs, &model.initial_state());
        let body_z = Vec3::new(
            inputs.R_hub.0[0][2],
            inputs.R_hub.0[1][2],
            inputs.R_hub.0[2][2],
        );
        let minus_f_dot_bz = -result.F_world.dot(body_z);
        assert!(
            minus_f_dot_bz > 0.0,
            "RAWES row122: force along +body_z: -F.body_z = {minus_f_dot_bz:+.6}"
        );
    }

    // -----------------------------------------------------------------------
    // Prandtl tip/hub loss formula verification
    // -----------------------------------------------------------------------

    fn prandtl_expected_tip(n: usize, x: f64, phi: f64) -> f64 {
        let f = (n as f64 / 2.0) * (1.0 - x) / (x * phi.sin().abs());
        (2.0 / PI) * f64::acos(f64::exp(-f).min(1.0))
    }

    fn prandtl_expected_hub(n: usize, x: f64, x_hub: f64, phi: f64) -> f64 {
        if x_hub <= 0.0 || (x - x_hub).abs() < 1e-12 {
            return 1.0;
        }
        let f = (n as f64 / 2.0) * (x - x_hub) / (x_hub * phi.sin().abs());
        (2.0 / PI) * f64::acos(f64::exp(-f).min(1.0))
    }

    #[test]
    fn test_prandtl_tip_loss_matches_formula() {
        let cases = [
            (2usize, 0.90_f64, 5.0_f64),
            (2, 0.95, 3.0),
            (4, 0.90, 5.0),
            (2, 0.80, 8.0),
            (3, 0.95, 4.0),
        ];
        for (n, x, phi_deg) in cases {
            let phi = phi_deg.to_radians();
            let expected = prandtl_expected_tip(n, x, phi);
            let got = prandtl_tip_loss(n, x, phi);
            assert!(
                (got - expected).abs() < 1e-12,
                "tip n={n} x={x} phi={phi_deg}deg: got {got:.10} expected {expected:.10}"
            );
        }
    }

    #[test]
    fn test_prandtl_tip_loss_boundary_cases() {
        // Far from tip: F -> 1
        assert!((prandtl_tip_loss(2, 0.3, 5_f64.to_radians()) - 1.0).abs() < 1e-4);
        // phi = 0: F = 1
        assert_eq!(prandtl_tip_loss(2, 0.9, 0.0), 1.0);
        // x = 1 (at tip): F = 1
        assert_eq!(prandtl_tip_loss(2, 1.0, 5_f64.to_radians()), 1.0);
        // More blades -> less tip loss
        let phi = 5_f64.to_radians();
        assert!(prandtl_tip_loss(4, 0.95, phi) > prandtl_tip_loss(2, 0.95, phi));
        // Larger phi -> more loss
        assert!(
            prandtl_tip_loss(2, 0.95, 8_f64.to_radians())
                < prandtl_tip_loss(2, 0.95, 2_f64.to_radians())
        );
    }

    #[test]
    fn test_prandtl_hub_loss_matches_formula() {
        let cases = [
            (2usize, 0.25_f64, 0.15_f64, 5.0_f64),
            (3, 0.20, 0.17, 4.0),
            (2, 0.30, 0.10, 8.0),
            (4, 0.18, 0.15, 5.0),
        ];
        for (n, x, x_hub, phi_deg) in cases {
            let phi = phi_deg.to_radians();
            let expected = prandtl_expected_hub(n, x, x_hub, phi);
            let got = prandtl_hub_loss(n, x, x_hub, phi);
            assert!(
                (got - expected).abs() < 1e-12,
                "hub n={n} x={x} x_hub={x_hub} phi={phi_deg}deg: got {got:.10} expected {expected:.10}"
            );
        }
    }

    #[test]
    fn test_prandtl_hub_loss_boundary_cases() {
        // Far from hub: F -> 1
        assert!((prandtl_hub_loss(2, 0.8, 0.1, 5_f64.to_radians()) - 1.0).abs() < 1e-4);
        // At hub: F = 1 (degenerate guard)
        assert_eq!(prandtl_hub_loss(2, 0.15, 0.15, 5_f64.to_radians()), 1.0);
        // phi = 0: F = 1
        assert_eq!(prandtl_hub_loss(2, 0.3, 0.15, 0.0), 1.0);
        // x_hub = 0: F = 1 (no hub cutout)
        assert_eq!(prandtl_hub_loss(2, 0.3, 0.0, 5_f64.to_radians()), 1.0);
        // More blades -> less hub loss
        let phi = 5_f64.to_radians();
        assert!(prandtl_hub_loss(4, 0.18, 0.15, phi) > prandtl_hub_loss(2, 0.18, 0.15, phi));
        // Closer to hub -> more loss
        assert!(prandtl_hub_loss(2, 0.16, 0.15, phi) < prandtl_hub_loss(2, 0.25, 0.15, phi));
    }

    // -----------------------------------------------------------------------
    // Near-hover boundary: solve_bem_element must not jump across v_climb = 0
    // -----------------------------------------------------------------------

    fn hover_polar() -> crate::polar::LinearPolar {
        crate::polar::LinearPolar::from_properties(&LinearPolarParameters {
            CL0: 0.0,
            CL_alpha_per_rad: 5.7,
            CD0: 0.01,
            alpha_stall_deg: 15.0,
        })
    }

    /// lambda_r must be continuous across v_climb = 0 within the near-hover
    /// band (|v_climb| < MIN_LAMBDA_CLIMB_WINDMILL * tip_speed). Previously,
    /// the sign-based seeding caused a large jump right at v_climb = 0.
    #[test]
    fn test_hover_boundary_continuity() {
        let polar = hover_polar();
        let omega = 1250.0 * std::f64::consts::PI / 30.0;
        let radius_m = 1.143;
        let collective = f64::to_radians(8.0);
        let eps = 1e-4; // tiny v_climb -- well inside near-hover band

        let geom = BEMElementGeometry::new(
            0.8 * radius_m,
            0.05,
            0.1905,
            0.0,
            omega,
            1.225,
            2,
            radius_m,
            &polar,
            false,
            0.0,
        );
        let above = solve_bem_element(&geom, collective, eps, 0.0);
        let below = solve_bem_element(&geom, collective, -eps, 0.0);

        let jump = (above.lambda_r - below.lambda_r).abs();
        assert!(
            jump < 0.01,
            "lambda_r jumps {jump:.4} across v_climb = 0; near-hover continuity broken"
        );
    }

    /// Deep descent (well outside windmill threshold) must return a negative
    /// lambda_r (net upward inflow) -- the helicopter quadratic descent root.
    #[test]
    fn test_deep_descent_negative_lambda_r() {
        let polar = hover_polar();
        let omega = 50.0; // slow spin so v_climb dominates
        let radius_m = 1.143;
        let geom = BEMElementGeometry::new(
            0.8 * radius_m,
            0.05,
            0.1905,
            0.0,
            omega,
            1.225,
            2,
            radius_m,
            &polar,
            false,
            0.0,
        );
        let elem = solve_bem_element(&geom, f64::to_radians(5.0), -15.0, 0.0);
        assert!(
            elem.lambda_r < 0.0,
            "deep descent should give lambda_r < 0, got {:.4}",
            elem.lambda_r
        );
    }

    // -----------------------------------------------------------------------
    // Element-level physics (ported from Python TestBEMElementConvergence)
    // -----------------------------------------------------------------------

    fn ct_polar() -> crate::polar::LinearPolar {
        crate::polar::LinearPolar::from_properties(&LinearPolarParameters {
            CL0: 0.0,
            CL_alpha_per_rad: 2.0 * std::f64::consts::PI,
            CD0: 0.008,
            alpha_stall_deg: 15.0,
        })
    }

    fn ct_geom<'a>(
        omega: f64,
        v_climb: f64,
        polar: &'a crate::polar::LinearPolar,
        use_tip_loss: bool,
    ) -> (BEMElementGeometry<'a, crate::polar::LinearPolar>, f64) {
        let radius_m = 1.143_f64;
        let r = 0.8 * radius_m;
        let dr = 0.05;
        let geom = BEMElementGeometry::new(
            r,
            dr,
            0.1905,
            0.0,
            omega,
            1.225,
            2,
            radius_m,
            polar,
            use_tip_loss,
            0.0,
        );
        let _ = v_climb; // consumed by caller
        (geom, dr)
    }

    /// Momentum balance residual must be near zero at convergence for all
    /// conditions: hover, autorotation (upward wind), and climb.
    #[test]
    fn test_momentum_balance_residual_hover_and_climb() {
        let polar = ct_polar();
        let cases: &[(f64, f64, f64)] = &[
            (8.0, 1250.0, 0.0),   // hover
            (5.0, 1250.0, 0.0),   // hover low pitch
            (12.0, 1250.0, 0.0),  // hover high pitch
            (8.0, 1000.0, 0.0),   // hover lower RPM
            (5.0, 1000.0, -10.0), // autorotation (upward wind)
            (8.0, 1250.0, 5.0),   // climbing
        ];
        for &(coll_deg, omega_rpm, v_climb) in cases {
            let omega = omega_rpm * std::f64::consts::PI / 30.0;
            let radius_m = 1.143_f64;
            let r = 0.8 * radius_m;
            let dr = 0.05_f64;
            let geom = BEMElementGeometry::new(
                r, dr, 0.1905, 0.0, omega, 1.225, 2, radius_m, &polar, true, 0.0,
            );
            let elem = solve_bem_element(&geom, f64::to_radians(coll_deg), v_climb, 0.0);
            let scale = (elem.d_t / dr).abs().max(1.0);
            assert!(
                elem.momentum_residual / scale < 1e-3,
                "coll={coll_deg} rpm={omega_rpm} v_climb={v_climb}: \
                 momentum residual {:.2e} / scale {:.2e} = {:.2e}",
                elem.momentum_residual,
                scale,
                elem.momentum_residual / scale
            );
        }
    }

    /// Stopped rotor (omega=0) must produce zero forces.
    #[test]
    fn test_zero_omega_gives_zero_forces() {
        let polar = ct_polar();
        let geom =
            BEMElementGeometry::new(0.8, 0.05, 0.2, 0.0, 0.0, 1.225, 2, 1.0, &polar, true, 0.0);
        let elem = solve_bem_element(&geom, f64::to_radians(8.0), 0.0, 0.0);
        assert_eq!(elem.d_t, 0.0, "stopped rotor: dT must be zero");
        assert_eq!(elem.d_q, 0.0, "stopped rotor: dQ must be zero");
    }

    /// Hover must produce downward induction (lambda_r > 0).
    #[test]
    fn test_hover_lambda_r_positive() {
        let polar = ct_polar();
        let omega = 1250.0 * std::f64::consts::PI / 30.0;
        let radius_m = 1.143_f64;
        let geom = BEMElementGeometry::new(
            0.8 * radius_m,
            0.05,
            0.1905,
            0.0,
            omega,
            1.225,
            2,
            radius_m,
            &polar,
            false,
            0.0,
        );
        let elem = solve_bem_element(&geom, f64::to_radians(8.0), 0.0, 0.0);
        assert!(
            elem.lambda_r > 0.0,
            "hover induced flow must be downward (lambda_r > 0), got {:.4}",
            elem.lambda_r
        );
    }

    /// At hover, dT must match momentum theory: dT = 4*pi*r*dr*rho*vi^2 (F=1, no tip loss).
    #[test]
    fn test_hover_thrust_matches_momentum_theory() {
        let polar = ct_polar();
        let omega = 1250.0 * std::f64::consts::PI / 30.0;
        let radius_m = 1.143_f64;
        let r = 0.8 * radius_m;
        let dr = 0.02_f64;
        let rho = 1.225_f64;
        let geom = BEMElementGeometry::new(
            r, dr, 0.1905, 0.0, omega, rho, 2, radius_m, &polar, false, 0.0,
        );
        let elem = solve_bem_element(&geom, f64::to_radians(8.0), 0.0, 0.0);
        let v_i = elem.lambda_r * omega * radius_m;
        let dt_momentum = 4.0 * std::f64::consts::PI * r * dr * rho * v_i * v_i;
        let rel_err = (elem.d_t - dt_momentum).abs() / dt_momentum;
        assert!(
            rel_err < 0.02,
            "hover dT {:.4e} vs momentum {:.4e}: rel_err {:.2}% > 2%",
            elem.d_t,
            dt_momentum,
            rel_err * 100.0
        );
    }

    /// Omega decay continuity: rotor spinning down to zero with no wind/climb
    /// must be smooth (no force/torque jumps), and torque must always resist
    /// rotation (Q_spin > 0) since there's no wind to drive it -- second law
    /// of thermodynamics: with no energy input from the air, the rotor can
    /// only ever lose rotational energy, never gain it. This must hold for
    /// every collective setting (positive AND negative pitch), not just one.
    ///
    /// Q_spin convention (see AeroResult docs / aero_io.rs): positive opposes
    /// rotation (drag, motor must supply +torque to sustain); negative means
    /// the airflow is driving the rotor (autorotation). With no wind, only
    /// Q_spin > 0 is physically possible.
    #[test]
    fn test_omega_decay_continuity_to_zero() {
        use crate::rotor_definition::BladeGeometry;

        // Caradonna-Tung rotor: 2 blades, R=1.143 m, NACA 0012, no twist
        let defn = RotorDefinition {
            blade: BladeGeometry {
                n_blades: 2,
                radius_m: 1.143,
                root_cutout_m: 0.1,
                chord_m: 0.1905,
                twist_deg: 0.0,
                n_elements: 20,
                tip_loss: true,
                r_stations_m: Vec::new(),
                chord_stations_m: Vec::new(),
                twist_stations_deg: Vec::new(),
            },
            airfoil: LinearPolarParameters {
                CL0: 0.0,
                CL_alpha_per_rad: 2.0 * PI,
                CD0: 0.008,
                alpha_stall_deg: 15.0,
            },
            control: None,
            pitch_actuation: PitchActuation::DirectMechanical,
            flap: None,
            name: "Caradonna-Tung".to_string(),
            description: String::new(),
        };

        let polar = crate::polar::LinearPolar::from_properties(&defn.airfoil);
        let bem = QuasiStaticBEM::build(defn, 36, polar);

        let rho = 1.225;
        let base_omega = 1250.0 * PI / 30.0;

        // Sweep across positive and negative collectives -- no wind means
        // the rotor must decelerate (Q_spin > 0) regardless of blade pitch
        // sign.
        for &collective_deg in &[-10.0_f64, -5.0, -2.0, 0.0, 2.0, 5.0, 8.0, 12.0] {
            // Geometric decay (constant ratio per step) so that thrust/torque
            // coefficients -- which are what should stay continuous, not raw
            // magnitudes (those naturally scale like omega^2 by momentum
            // theory) -- can be compared step-to-step on an equal footing.
            // Stop at a small but nonzero omega; exact omega=0 is covered by
            // test_zero_omega_gives_zero_forces.
            let mut omega = base_omega;
            let mut prev_ct = None;
            let mut prev_cq = None;

            for _ in 0..30 {
                let inputs = RotorInputs {
                    collective_rad: collective_deg.to_radians(),
                    tilt_lon: 0.0,
                    tilt_lat: 0.0,
                    R_hub: Mat3::eye(),
                    v_hub_world: Vec3::zero(),
                    wind_world: Vec3::zero(),
                    omega_rad_s: omega,
                    rho_kg_m3: rho,
                };

                let (res, _) = bem.compute_forces(&inputs, &QuasiStaticRotorState::default());

                let thrust = res.F_world[2].abs();
                let torque_raw = res.Q_spin;

                // No wind driving the rotor: torque must always resist
                // rotation (never accelerate it) for every collective.
                assert!(
                    torque_raw >= 0.0,
                    "collective={collective_deg} deg, omega={omega:.4}: Q_spin must be \
                     >= 0 with no wind (rotor can only decelerate), got {torque_raw:.3e} N.m"
                );

                // Non-dimensionalize by omega^2 (momentum theory: T, Q ~
                // omega^2 at fixed collective/geometry) so we're checking
                // for genuine regime-switch discontinuities, not the
                // expected quadratic falloff of raw magnitude.
                let ct = thrust / (omega * omega);
                let cq = torque_raw / (omega * omega);

                if let Some(prev) = prev_ct {
                    if prev > 1e-6 {
                        let ratio = ct / prev;
                        assert!(
                            0.5 < ratio && ratio < 2.0,
                            "collective={collective_deg} deg: thrust coefficient \
                             (T/omega^2) jumped at omega={omega:.4}: {prev:.3e} to \
                             {ct:.3e} (ratio {ratio:.2})"
                        );
                    }
                }
                if let Some(prev) = prev_cq {
                    if prev > 1e-8 {
                        let ratio = cq / prev;
                        assert!(
                            0.5 < ratio && ratio < 2.0,
                            "collective={collective_deg} deg: torque coefficient \
                             (Q/omega^2) jumped at omega={omega:.4}: {prev:.3e} to \
                             {cq:.3e} (ratio {ratio:.2})"
                        );
                    }
                }

                prev_ct = Some(ct);
                prev_cq = Some(cq);

                omega *= 0.7; // ~30 steps from 1250 RPM down to ~1e-4 RPM
            }
        }
    }

    /// True time-integration spindown: uses the built-in `AeroModel::step()`
    /// trait method (aero_model.rs) for the inflow-state update, plus the
    /// caller-owned mechanical spin ODE via `crate::mechanical::omega_derivative`
    /// / `crate::mechanical::step_omega` -- the single canonical place this
    /// ODE is evaluated and stepped (see mechanical.rs module docs).
    /// Integrates omega from a nonzero start, with no wind and no motor
    /// torque, and checks it decays toward zero monotonically and smoothly --
    /// exactly the caller pattern documented in AeroModel::step's doc
    /// comment.
    ///
    /// Note: with pure aerodynamic (quadratic) drag and zero other losses,
    /// the continuous-time ODE only reaches exactly omega=0 asymptotically
    /// (1/t tail) -- so this checks decay to a small fraction of the start
    /// speed, not literal zero (test_zero_omega_gives_zero_forces already
    /// covers the omega==0 endpoint directly).
    #[test]
    fn test_omega_spindown_time_integration() {
        use crate::aero_model::IntegrationMethod;
        use crate::rotor_definition::BladeGeometry;

        let defn = RotorDefinition {
            blade: BladeGeometry {
                n_blades: 2,
                radius_m: 1.143,
                root_cutout_m: 0.1,
                chord_m: 0.1905,
                twist_deg: 0.0,
                n_elements: 20,
                tip_loss: true,
                r_stations_m: Vec::new(),
                chord_stations_m: Vec::new(),
                twist_stations_deg: Vec::new(),
            },
            airfoil: LinearPolarParameters {
                CL0: 0.0,
                CL_alpha_per_rad: 2.0 * PI,
                CD0: 0.008,
                alpha_stall_deg: 15.0,
            },
            control: None,
            pitch_actuation: PitchActuation::DirectMechanical,
            flap: None,
            name: "Caradonna-Tung".to_string(),
            description: String::new(),
        };

        let polar = crate::polar::LinearPolar::from_properties(&defn.airfoil);
        let bem = QuasiStaticBEM::build(defn, 36, polar);

        let rho = 1.225;
        let i_ode_kgm2 = 1.0_f64;
        let omega_start = 1250.0 * PI / 30.0;
        let max_steps = 2000;

        for &collective_deg in &[-10.0_f64, -5.0, 5.0, 8.0] {
            let mut omega = omega_start;
            let mut state = bem.initial_state();
            let mut dt = 1e-3_f64;
            let mut reached_target = false;

            for step in 0..max_steps {
                let inputs = RotorInputs {
                    collective_rad: collective_deg.to_radians(),
                    tilt_lon: 0.0,
                    tilt_lat: 0.0,
                    R_hub: Mat3::eye(),
                    v_hub_world: Vec3::zero(),
                    wind_world: Vec3::zero(),
                    omega_rad_s: omega,
                    rho_kg_m3: rho,
                };

                // Built-in trait method: advances the inflow state (a
                // no-op here -- QuasiStaticBEM carries zero inflow DOF).
                let (res, new_state) =
                    bem.step(&inputs, &state, dt, IntegrationMethod::ExplicitEuler);
                state = new_state;

                // Mechanical spin ODE is caller-owned by design (see
                // AeroModel::step doc comment) -- call into the single
                // canonical mechanical ODE functions rather than
                // re-deriving the formula locally.
                let motor_torque_nm = 0.0;
                let d_omega = crate::mechanical::omega_derivative(
                    omega,
                    res.Q_spin,
                    motor_torque_nm,
                    i_ode_kgm2,
                    0.0,
                );

                // No wind/motor driving the rotor: it must never
                // accelerate (Q_spin >= 0 always, per the physics
                // established in test_omega_decay_continuity_to_zero).
                assert!(
                    d_omega <= 1e-9,
                    "collective={collective_deg} deg, step {step}: omega must \
                     not accelerate with no wind/motor torque, d_omega={d_omega:.3e}"
                );

                // Adapt dt to keep the per-step relative change bounded
                // (~5%) -- necessary because the aerodynamic torque
                // shrinks as omega^2, so a fixed dt would either take
                // forever near omega_start or blow past zero near the
                // end.
                if d_omega.abs() > 1e-12 {
                    dt = (0.05 * omega / d_omega.abs()).clamp(1e-6, 5.0);
                }

                let new_omega = crate::mechanical::step_omega(
                    omega,
                    res.Q_spin,
                    motor_torque_nm,
                    i_ode_kgm2,
                    dt,
                    0.0,
                );
                assert!(
                    new_omega <= omega + 1e-9,
                    "collective={collective_deg} deg, step {step}: omega \
                     increased from {omega:.6} to {new_omega:.6} during spindown"
                );

                omega = new_omega;

                if omega <= 1e-3 * omega_start {
                    reached_target = true;
                    break;
                }
            }

            assert!(
                reached_target,
                "collective={collective_deg} deg: omega failed to decay to \
                 0.1% of its starting value within {max_steps} steps (stuck \
                 at {omega:.4} rad/s -- possible spindown regression)"
            );
        }
    }
}
