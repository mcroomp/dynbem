// Level 1 BEM: helicopter momentum quadratic + wind-turbine windmill (Brent).
// See ../CLAUDE.md, dynbem/bem.py for the full physics.

use std::f64::consts::PI;

use crate::aero_io::{AeroResult, RotorInputs};
use crate::aero_model::AeroModel;
use crate::bem_common::{assemble_result, kinematics, ElementCtx, PsiKernel, RadialGrid, SweepCtx};
use crate::common::{EPS_DENOM, EPS_OMEGA_R, MIN_LOSS_FACTOR};
use crate::cyclic::cyclic_coeffs;
use crate::polar::{Polar, PolarKind};
use crate::rotor_definition::RotorDefinition;
use crate::rotor_state::QuasiStaticRotorState;

const MAX_BEM_ITER: usize = 60;
const BEM_TOL: f64 = 1e-7;

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
    pub momentum_residual: f64, // |4*F*lambda_r*(lambda_r - lambda_c) - sigma_r*cn*(lambda_r^2 + x^2)|
}

// ---------------------------------------------------------------------------
// Helicopter momentum quadratic
// ---------------------------------------------------------------------------

/// Helicopter momentum-BEM solver for one annulus.
///
/// Fixed-point iteration on (lambda_r, a_prime) with 50% under-relaxation;
/// the converged root of the quadratic is selected explicitly by sign of
/// lambda_c (climb -> positive root, descent -> negative). Reverse-flow
/// region (v_t < 0) breaks out early and returns zero forces -- caller
/// is responsible for the surrounding ψ-loop's reverse-flow skip.
#[allow(clippy::too_many_arguments)]
pub fn solve_bem_element(
    r: f64,
    dr: f64,
    chord: f64,
    twist_rad: f64,
    collective_rad: f64,
    omega: f64,
    v_climb: f64,
    rho: f64,
    n_blades: usize,
    radius_m: f64,
    polar: &PolarKind,
    use_tip_loss: bool,
    v_t_extra: f64,
    root_cutout_m: f64,
) -> BEMElementResult {
    let omega_r = omega * radius_m;
    if omega_r < EPS_OMEGA_R {
        return BEMElementResult::default();
    }
    let inv_radius_m = if radius_m > 0.0 { 1.0 / radius_m } else { 0.0 };
    let inv_omega_r = 1.0 / omega_r;
    let inv_r = if r > 0.0 { 1.0 / r } else { 0.0 };
    let x = r * inv_radius_m;
    let x_hub = if radius_m > 0.0 {
        root_cutout_m * inv_radius_m
    } else {
        0.0
    };
    let sigma_r = (n_blades as f64) * chord * inv_r / (2.0 * PI);
    let theta = collective_rad + twist_rad;
    let lambda_c = v_climb * inv_omega_r;

    let mut lambda_r = if lambda_c >= 0.0 {
        (lambda_c + 0.03).max(0.02)
    } else {
        (lambda_c * 0.85).min(-0.02)
    };
    let mut a_prime: f64 = 0.0;

    for _ in 0..MAX_BEM_ITER {
        let v_a = lambda_r * omega_r;
        let v_t = omega * r * (1.0 + a_prime) + v_t_extra;
        if v_t < EPS_DENOM {
            break;
        }

        let phi = v_a.atan2(v_t);
        let alpha = theta - phi;
        let (cl, cd) = polar.cl_cd(alpha);

        let cos_p = phi.cos();
        let sin_p = phi.sin();
        let sin_phi_abs = sin_p.abs();

        let f_loss = if use_tip_loss {
            (prandtl_tip_loss_from_sin_abs(n_blades, x, sin_phi_abs)
                * prandtl_hub_loss_from_sin_abs(n_blades, x, x_hub, sin_phi_abs))
            .max(MIN_LOSS_FACTOR)
        } else {
            1.0
        };
        let cn = cl * cos_p - cd * sin_p;
        let ct = cl * sin_p + cd * cos_p;

        let k = sigma_r * cn / (4.0 * f_loss);

        let lambda_r_new = if (k - 1.0).abs() > 1e-6 {
            let disc = (lambda_c * lambda_c - 4.0 * (k - 1.0) * k * x * x).max(0.0);
            let sq = disc.sqrt();
            let denom = 2.0 * (k - 1.0);
            let r1 = (-lambda_c + sq) / denom;
            let r2 = (-lambda_c - sq) / denom;
            if lambda_c >= 0.0 {
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
        } else if lambda_c.abs() > 1e-8 {
            -k * x * x / lambda_c
        } else {
            x * k.max(0.0).sqrt()
        };
        let lambda_r_new = lambda_r_new.clamp(-2.0, 2.0);

        let sc = sin_p * cos_p;
        let a_prime_new = if sc.abs() > 1e-8 && ct.abs() > 1e-10 {
            let ap_denom = 4.0 * f_loss * sc / (sigma_r * ct) - 1.0;
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

    let v_a = lambda_r * omega_r;
    let v_t = omega * r * (1.0 + a_prime) + v_t_extra;
    let v_rel_sq = v_a * v_a + v_t * v_t;
    let phi = v_a.atan2(v_t);
    let alpha = theta - phi;
    let (cl, cd) = polar.cl_cd(alpha);
    let cos_p = phi.cos();
    let sin_p = phi.sin();
    let sin_phi_abs = sin_p.abs();
    let cn = cl * cos_p - cd * sin_p;
    let ct = cl * sin_p + cd * cos_p;
    let q = 0.5 * rho * v_rel_sq * chord * dr * (n_blades as f64);

    // Momentum-balance residual at the converged state:
    //   |4*F*lambda_r*(lambda_r - lambda_c) - sigma_r*cn*(lambda_r^2 + x^2)|
    // Same definition as the legacy Python BEMElementResult.momentum_residual.
    let f_loss = if use_tip_loss {
        (prandtl_tip_loss_from_sin_abs(n_blades, x, sin_phi_abs)
            * prandtl_hub_loss_from_sin_abs(n_blades, x, x_hub, sin_phi_abs))
        .max(MIN_LOSS_FACTOR)
    } else {
        1.0
    };
    let momentum_residual = (4.0 * f_loss * lambda_r * (lambda_r - lambda_c)
        - sigma_r * cn * (lambda_r * lambda_r + x * x))
        .abs();

    BEMElementResult {
        d_t: q * cn,
        d_q: q * ct * r,
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

#[allow(clippy::too_many_arguments)]
fn induction_at_phi(
    phi: f64,
    theta: f64,
    lam_local: f64,
    sigma_r: f64,
    n_blades: usize,
    x: f64,
    x_hub: f64,
    polar: &PolarKind,
    use_tip_loss: bool,
) -> Option<WindmillInductions> {
    let sin_phi = phi.sin();
    let cos_phi = phi.cos();
    if sin_phi.abs() < 1e-12 {
        return None;
    }
    let alpha = theta - phi;
    let (cl, cd) = polar.cl_cd(alpha);
    let cn = cl * cos_phi - cd * sin_phi;
    let ct = cl * sin_phi + cd * cos_phi;
    let f_loss = if use_tip_loss {
        prandtl_tip_loss(n_blades, x, phi) * prandtl_hub_loss(n_blades, x, x_hub, phi)
    } else {
        1.0
    };
    let f_loss = f_loss.max(MIN_LOSS_FACTOR);
    let sin2 = sin_phi * sin_phi;
    if cn <= EPS_DENOM {
        return None;
    }
    let k_axial = sigma_r * cn / (4.0 * f_loss * sin2);
    let mut a = k_axial / (1.0 + k_axial);
    if a > 0.4 {
        // Buhl quadratic in the turbulent-wake state
        let k = sigma_r * cn / sin2;
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
        let k_tan = sigma_r * ct / (4.0 * f_loss * sc);
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

/// Wind-turbine BEM solver for one annulus (axial upflow only).
///
/// Ning 2014 reformulation: residual-on-phi with Brent over (-pi/2, 0),
/// plus the Buhl quadratic for the turbulent-wake state (a > 0.4) where
/// classical momentum theory breaks down. Returns None when no sign-change
/// bracket exists or the iteration leaves the valid windmill regime
/// (0 < a < 1 and Cn > 0) -- caller falls back to the helicopter
/// quadratic in those cases.
#[allow(clippy::too_many_arguments)]
fn solve_bem_element_windmill(
    r: f64,
    dr: f64,
    chord: f64,
    twist_rad: f64,
    collective_rad: f64,
    omega: f64,
    v_climb: f64,
    rho: f64,
    n_blades: usize,
    radius_m: f64,
    polar: &PolarKind,
    use_tip_loss: bool,
    root_cutout_m: f64,
) -> Option<BEMElementResult> {
    if v_climb >= -EPS_DENOM {
        return None;
    }
    let u_up = -v_climb;
    let omega_r = omega * radius_m;
    if omega_r < EPS_OMEGA_R {
        return None;
    }
    let inv_radius_m = if radius_m > 0.0 { 1.0 / radius_m } else { 0.0 };
    let inv_u_up = 1.0 / u_up;
    let inv_omega_r = 1.0 / omega_r;
    let inv_r = if r > 0.0 { 1.0 / r } else { 0.0 };
    let x = r * inv_radius_m;
    let x_hub = if radius_m > 0.0 {
        root_cutout_m * inv_radius_m
    } else {
        0.0
    };
    let sigma_r = (n_blades as f64) * chord * inv_r / (2.0 * PI);
    let theta = collective_rad + twist_rad;
    let lam_local = omega * r * inv_u_up;

    let residual = |phi: f64| -> f64 {
        match induction_at_phi(
            phi,
            theta,
            lam_local,
            sigma_r,
            n_blades,
            x,
            x_hub,
            polar,
            use_tip_loss,
        ) {
            None => 1e3,
            Some(ind) => phi.sin() * (1.0 + ind.ap) * lam_local + phi.cos() * (1.0 - ind.a),
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
    let ind = induction_at_phi(
        phi_star,
        theta,
        lam_local,
        sigma_r,
        n_blades,
        x,
        x_hub,
        polar,
        use_tip_loss,
    )?;
    if ind.cn <= 0.0 || !(0.0..=1.0).contains(&ind.a) {
        return None;
    }
    let v_a = -(1.0 - ind.a) * u_up;
    let v_t = (1.0 + ind.ap) * omega * r;
    let v_rel_sq = v_a * v_a + v_t * v_t;
    let q = 0.5 * rho * v_rel_sq * chord * dr * (n_blades as f64);
    // Windmill solver works in (a, a') space, not (lambda_r, a_prime); reconstruct
    // lambda_r consistent with the rest of the pipeline (axial inflow ratio).
    let lambda_r = v_a * inv_omega_r;
    Some(BEMElementResult {
        d_t: q * ind.cn,
        d_q: q * ind.ct * r,
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

/// BEM-specific state for the shared psi-loop. Holds only what's *not*
/// already in SweepCtx or ElementCtx -- so no per-call constant is
/// duplicated. Everything here is either a model parameter (r_tip,
/// root_cutout_m, use_tip_loss) or a flight-state value that no other
/// model reads through the kernel interface (v_climb).
struct BemKernel {
    v_climb: f64,
    r_tip: f64,
    root_cutout_m: f64,
    use_tip_loss: bool,
}

impl PsiKernel for BemKernel {
    #[inline(always)]
    fn element(&mut self, sweep: &SweepCtx<'_>, ctx: &ElementCtx) -> (f64, f64) {
        // run_psi_loop already filtered v_t > 0; reconstruct v_t_extra from
        // the convention v_t = omega*r + v_t_extra so we don't add another
        // parameter to ElementCtx just for the BEM solver.
        let v_t_extra = ctx.v_t - sweep.omega * ctx.r;
        let elem = solve_bem_element(
            ctx.r,
            ctx.dr,
            ctx.chord,
            ctx.twist,
            ctx.col_psi,
            sweep.omega,
            self.v_climb,
            sweep.rho,
            sweep.n_b,
            self.r_tip,
            sweep.polar,
            self.use_tip_loss,
            v_t_extra,
            self.root_cutout_m,
        );
        (elem.d_t, elem.d_q)
    }
}

// ---------------------------------------------------------------------------
// QuasiStaticBEM: pyclass holding cached radial grid + polar
// ---------------------------------------------------------------------------

#[derive(Clone)]
pub struct QuasiStaticBEM {
    pub defn: RotorDefinition,
    pub n_psi_elements: usize,
    pub polar: PolarKind,
    pub grid: RadialGrid,
}

impl QuasiStaticBEM {
    pub fn build(defn: RotorDefinition, n_psi_elements: usize, polar: PolarKind) -> Self {
        let grid = RadialGrid::from_blade(&defn.blade);
        Self {
            defn,
            n_psi_elements,
            polar,
            grid,
        }
    }
}

impl AeroModel for QuasiStaticBEM {
    type State = QuasiStaticRotorState;

    fn initial_state(&self) -> Self::State {
        QuasiStaticRotorState::default()
    }

    // inflow_taus: trait default (all-infinity) is correct for the
    // quasi-static BEM model; no override needed.

    fn compute_forces(
        &self,
        inputs: &RotorInputs,
        state: &QuasiStaticRotorState,
    ) -> (AeroResult, QuasiStaticRotorState) {
        let blade = &self.defn.blade;
        let omega = inputs.omega_rad_s;
        let rho = inputs.rho_kg_m3;
        let r_tip = blade.radius_m;
        let n_blades = blade.n_blades;
        let use_tip_loss = self.defn.airfoil.tip_loss;
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

        let mut t_total: f64 = 0.0;
        let mut q_total: f64 = 0.0;
        let mut mx_hub: f64 = 0.0;
        let mut my_hub: f64 = 0.0;

        if (mu > 0.01 || has_cyclic) && omega > 1.0 {
            // Forward / cyclic: shared psi-loop with the iterative solver
            // injected via BemKernel::element override.
            let mut kernel = BemKernel {
                v_climb,
                r_tip,
                root_cutout_m: blade.root_cutout_m,
                use_tip_loss,
            };
            let sweep = SweepCtx {
                grid,
                polar: &self.polar,
                col: inputs.collective_rad,
                omega,
                omega_r,
                rho,
                n_b: n_blades,
                n_psi: self.n_psi_elements,
                n_psi_inv: 1.0 / (self.n_psi_elements as f64),
                v_in_hub_x: v_inplane_hub[0],
                v_in_hub_y: v_inplane_hub[1],
                theta_1c,
                theta_1s,
            };
            let (t, q, mx, my) = sweep.run(&mut kernel);
            t_total = t;
            q_total = q;
            mx_hub = mx;
            my_hub = my;
        } else {
            // Axial: try wind-turbine windmill solver first when v_climb < 0,
            // fall back to helicopter quadratic per element.
            for i_r in 0..n {
                let r = r_arr[i_r];
                let mut elem: Option<BEMElementResult> = None;
                if v_climb < -EPS_DENOM {
                    elem = solve_bem_element_windmill(
                        r,
                        dr,
                        chord[i_r],
                        twist[i_r],
                        inputs.collective_rad,
                        omega,
                        v_climb,
                        rho,
                        n_blades,
                        r_tip,
                        &self.polar,
                        use_tip_loss,
                        blade.root_cutout_m,
                    );
                }
                let elem = elem.unwrap_or_else(|| {
                    solve_bem_element(
                        r,
                        dr,
                        chord[i_r],
                        twist[i_r],
                        inputs.collective_rad,
                        omega,
                        v_climb,
                        rho,
                        n_blades,
                        r_tip,
                        &self.polar,
                        use_tip_loss,
                        0.0,
                        blade.root_cutout_m,
                    )
                });
                t_total += elem.d_t;
                q_total += elem.d_q;
            }
        }

        let result = assemble_result(t_total, q_total, mx_hub, my_hub, hub_axis, &inputs.R_hub);

        let derivative = QuasiStaticRotorState;
        // suppress unused warning (use_tip_loss read inside windmill helper).
        let _ = use_tip_loss;
        (result, derivative)
    }
}
