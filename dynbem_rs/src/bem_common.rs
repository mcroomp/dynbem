// Shared infrastructure for the BEM models.
// See ../CLAUDE.md "Shared BEM infrastructure".

use crate::aero_io::{AeroResult, Mat3, RotorInputs, Vec3};
use crate::common::{EPS_OMEGA_R, VRS_DESCENT_THRESHOLD, V_T_HOVER_FLOOR_FRAC};
use crate::polar::{Polar, PolarKind};
use crate::rotor_definition::BladeGeometry;
use std::f64::consts::PI;

/// Cached fixed radial geometry for a BEM kernel.
///
/// Per-station chord/twist arrays mean we can support BladeGeometry with
/// radial-station arrays (wind-turbine blades) AND scalar chord_m/twist_deg
/// (helicopters) through a single uniform interface in the inner loop --
/// no branches per element.
#[derive(Clone, Debug)]
pub struct RadialGrid {
    pub dr: f64,
    pub r_mid: Vec<f64>,     // n
    pub x_mid: Vec<f64>,     // n = r_mid / R
    pub x_hub: f64,          // root_cutout / R
    pub chord: Vec<f64>,     // n  -- per-station chord (m)
    pub twist_rad: Vec<f64>, // n  -- per-station twist (rad)
}

impl RadialGrid {
    pub fn from_blade(blade: &BladeGeometry) -> Self {
        let r_root = blade.root_cutout_m;
        let r_tip = blade.radius_m;
        let n = blade.n_elements;
        let dr = (r_tip - r_root) / (n as f64);
        let mut r_mid = Vec::with_capacity(n);
        let mut x_mid = Vec::with_capacity(n);
        let mut chord = Vec::with_capacity(n);
        let mut twist_rad = Vec::with_capacity(n);
        for i in 0..n {
            let r = r_root + (i as f64 + 0.5) * dr;
            r_mid.push(r);
            x_mid.push(if r_tip > 0.0 { r / r_tip } else { 0.0 });
            chord.push(blade.chord_at(r));
            twist_rad.push(blade.twist_at(r).to_radians());
        }
        let x_hub = if r_tip > 0.0 { r_root / r_tip } else { 0.0 };
        Self {
            dr,
            r_mid,
            x_mid,
            x_hub,
            chord,
            twist_rad,
        }
    }
}

/// Tabulate any polar onto contiguous arrays for the JIT-equivalent inner
/// loop. TabulatedPolar passes its arrays through; analytical polars get
/// sampled to 4001 points over [-pi/2, pi/2] (matching the Python version).
#[derive(Clone, Debug)]
pub struct PolarTable {
    pub alpha: Vec<f64>,
    pub cl: Vec<f64>,
    pub cd: Vec<f64>,
}

impl PolarTable {
    pub fn from_polar(polar: &PolarKind) -> Self {
        match polar {
            PolarKind::Tabulated(p) => Self {
                alpha: p.alpha.clone(),
                cl: p.cl.clone(),
                cd: p.cd.clone(),
            },
            PolarKind::Linear(_) => {
                let n = 4001usize;
                let mut alpha = Vec::with_capacity(n);
                let mut cl = vec![0.0; n];
                let mut cd = vec![0.0; n];
                let amin = -0.5 * PI;
                let amax = 0.5 * PI;
                let step = (amax - amin) / ((n - 1) as f64);
                for i in 0..n {
                    alpha.push(amin + (i as f64) * step);
                }
                polar.cl_cd_into(&alpha, &mut cl, &mut cd);
                Self { alpha, cl, cd }
            }
        }
    }

    /// Scalar interp at one alpha; same semantics as numpy.interp + the
    /// binary search in the Python _interp_polar.
    #[inline]
    pub fn interp(&self, alpha: f64) -> (f64, f64) {
        let a = &self.alpha[..];
        let n = a.len();
        if alpha <= a[0] {
            return (self.cl[0], self.cd[0]);
        }
        if alpha >= a[n - 1] {
            return (self.cl[n - 1], self.cd[n - 1]);
        }
        let mut lo = 0usize;
        let mut hi = n - 1;
        while hi - lo > 1 {
            let mid = (lo + hi) >> 1;
            if a[mid] <= alpha {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let a_lo = a[lo];
        let a_hi = a[hi];
        let t = (alpha - a_lo) / (a_hi - a_lo);
        let cl = self.cl[lo] + t * (self.cl[hi] - self.cl[lo]);
        let cd = self.cd[lo] + t * (self.cd[hi] - self.cd[lo]);
        (cl, cd)
    }
}

// ---------------------------------------------------------------------------
// Once-per-call kinematics. Identical across BEM / Pitt-Peters / Oye; runs
// outside any inner loop, so abstracting it has no autovectorization cost.
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
pub struct Kinematics {
    pub omega_r: f64,
    pub hub_axis: Vec3,
    pub v_climb: f64,
    pub v_inplane: Vec3,
    pub v_edge: f64,
    pub v_inplane_hub: Vec3,
    /// Advance ratio mu = v_edge / max(omega_r, EPS_OMEGA_R).
    pub mu: f64,
}

#[inline]
pub fn kinematics(inputs: &RotorInputs, omega: f64, r_tip: f64) -> Kinematics {
    let omega_r = omega * r_tip;
    let hub_axis = inputs.R_hub * Vec3::new(0.0, 0.0, 1.0);
    let v_rel = inputs.wind_world - inputs.v_hub_world;
    let v_climb = v_rel.dot(hub_axis);
    let v_inplane = v_rel - hub_axis * v_climb;
    let v_edge = v_inplane.norm();
    let v_inplane_hub = inputs.R_hub.transpose() * v_inplane;
    let mu = v_edge / omega_r.max(EPS_OMEGA_R);
    Kinematics {
        omega_r,
        hub_axis,
        v_climb,
        v_inplane,
        v_edge,
        v_inplane_hub,
        mu,
    }
}

// ---------------------------------------------------------------------------
// VRS regime detection. v_h is the hover induced velocity (positive sqrt
// of T/(2 rho A)); lam2 = V_descent / V_h is the descent-positive ratio used
// in Leishman's polynomial. in_vrs picks out 0 < lam2 < 2 while in descent.
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
pub struct VrsRegime {
    pub v_h: f64,
    pub lam2: f64,
    pub in_vrs: bool,
}

/// Mass-flow speed at the disk (the Glauert V_T scalar).
///
///     V_T = sqrt(v_edge^2 + (v_climb + v0_axial)^2)
///
/// `v0_axial` is the axial (induced) component in m/s. Floored at
/// `V_T_HOVER_FLOOR_FRAC * max(omega_r, 1)` to keep the Pitt-Peters / Oye
/// time constants finite at hover / zero thrust.
#[inline]
pub fn v_t_disk(v_edge: f64, v_climb: f64, v0_axial: f64, omega_r: f64) -> f64 {
    (v_edge * v_edge + (v_climb + v0_axial).powi(2))
        .sqrt()
        .max(V_T_HOVER_FLOOR_FRAC * omega_r.max(1.0))
}

#[inline]
pub fn vrs_regime(t_total: f64, v_climb: f64, rho: f64, area: f64) -> VrsRegime {
    let t_pos = t_total.max(0.0);
    let v_h = if t_pos > EPS_OMEGA_R {
        (t_pos / (2.0 * rho * area)).sqrt()
    } else {
        0.0
    };
    let v_c = (-v_climb).max(0.0);
    let lam2 = if v_h > VRS_DESCENT_THRESHOLD {
        v_c / v_h
    } else {
        0.0
    };
    let in_vrs = v_climb < -VRS_DESCENT_THRESHOLD && lam2 > 0.0 && lam2 < 2.0;
    VrsRegime { v_h, lam2, in_vrs }
}

// ---------------------------------------------------------------------------
// AeroResult assembly. Same translation from (T, Q, Mx_hub, My_hub) to the
// world-frame outputs for every model.
// ---------------------------------------------------------------------------

#[inline]
pub fn assemble_result(
    t_total: f64,
    q_total: f64,
    mx_hub: f64,
    my_hub: f64,
    hub_axis: Vec3,
    r_hub: &Mat3,
) -> AeroResult {
    let f_world = hub_axis * (-t_total);
    let mxyz_hub = Vec3::new(mx_hub, my_hub, 0.0);
    let m_orbital = *r_hub * mxyz_hub;
    let m_spin = hub_axis * q_total;
    AeroResult {
        F_world: f_world,
        M_orbital: m_orbital,
        Q_spin: q_total,
        M_spin: m_spin,
    }
}

// ---------------------------------------------------------------------------
// Per-element BEM integrand: given the prescribed axial velocity `v_a` plus
// the (sweep, element) contexts, return the element's (dT, dQ).
//
// `#[inline(always)]` preserves the autovectorization the per-model loops
// had before extraction -- LLVM sees the same arithmetic, just routed
// through one function. The opaque polar.cl_cd call is the same
// vectorization barrier it was before.
// ---------------------------------------------------------------------------

#[inline(always)]
pub fn element_force(v_a: f64, sweep: &SweepCtx<'_>, ctx: &ElementCtx) -> (f64, f64) {
    let v_t = ctx.v_t;
    let phi = v_a.atan2(v_t);
    let alpha = ctx.col_psi + ctx.twist - phi;
    let (cl, cd) = sweep.polar.cl_cd(alpha);
    let cos_p = phi.cos();
    let sin_p = phi.sin();
    let cn = cl * cos_p - cd * sin_p;
    let ct = cl * sin_p + cd * cos_p;
    let q = 0.5 * sweep.rho * (v_a * v_a + v_t * v_t) * ctx.chord * ctx.dr * (sweep.n_b as f64);
    (q * cn, q * ct * ctx.r)
}

// ---------------------------------------------------------------------------
// Shared psi-loop kernel: monomorphized via the PsiKernel trait so each model
// gets its own specialised copy, with #[inline(always)] callbacks. Codegen
// is bit-identical to the prior per-model hand-rolled loops -- the trait is
// used as a *static interface* (generic K), not runtime dispatch.
//
// Two override points, with sensible defaults:
//   - element(ctx): compute (dT, dQ) for one element. Default = the
//     prescribed-inflow path used by Pitt-Peters and Oye -- evaluate
//     lam_local, then run element_force. BEM overrides this to call its
//     iterative solve_bem_element instead.
//   - lam_local(i, cos_psi, sin_psi): the local inflow formula
//     (Pitt-Peters: lambda_total + x*(lam_c*cos psi + lam_s*sin psi);
//      Oye: lambda_climb + W[i]). Unused if element() is overridden.
//   - on_element(i, dt, inv_n_psi): per-element side effect (no-op for
//     PP/BEM; Oye uses it to accumulate the azimuth-averaged dT/dx the
//     W_qs solver needs downstream).
// ---------------------------------------------------------------------------

/// Per-element transients passed to a PsiKernel. Call-invariants
/// (omega_r, rho, n_b, polar) live in `SweepCtx` instead -- the kernel
/// receives both contexts and reads each field from its natural home.
pub struct ElementCtx {
    pub i: usize,
    pub cos_psi: f64,
    pub sin_psi: f64,
    pub r: f64,
    pub chord: f64,
    pub twist: f64,
    pub dr: f64,
    pub col_psi: f64,
    /// Tangential velocity at this (r, psi): `omega * r + v_t_extra`.
    pub v_t: f64,
}

pub trait PsiKernel {
    /// Element-level force computation. Default = prescribed-inflow path
    /// for Pitt-Peters (and any future model that just needs to define a
    /// `lam_local` formula). Override for models that need their own
    /// solver per element (BEM) or that need a per-element side effect
    /// (Oye accumulating azimuth-averaged dT/dx alongside the force
    /// computation -- a separate callback would just force the kernel
    /// to recompute `dt` or split state setup, so we collapse both into
    /// one override point).
    #[inline(always)]
    fn element(&mut self, sweep: &SweepCtx<'_>, ctx: &ElementCtx) -> (f64, f64) {
        let lam = self.lam_local(ctx.i, ctx.cos_psi, ctx.sin_psi);
        let v_a = lam * sweep.omega_r;
        element_force(v_a, sweep, ctx)
    }

    /// Local axial inflow ratio at element i, azimuth (cos_psi, sin_psi).
    /// Unused when the kernel overrides `element` directly; provide a
    /// stub default so override-everything kernels (BEM, Oye) don't have
    /// to implement it.
    #[inline(always)]
    #[allow(unused_variables)]
    fn lam_local(&self, i: usize, cos_psi: f64, sin_psi: f64) -> f64 {
        0.0
    }
}

/// Whole-sweep configuration: the call-invariant inputs that describe one
/// psi x radial pass over the rotor disk. Built once by each model's
/// compute_forces, then handed to run_psi_loop, which derives per-element
/// `ElementCtx` values from these fields each iteration. Symmetric with
/// `ElementCtx` (which describes one element rather than one sweep).
pub struct SweepCtx<'a> {
    pub grid: &'a RadialGrid,
    pub polar: &'a PolarKind,
    /// Base collective pitch (rad). Per-azimuth pitch is `col + theta_1c*cos psi + theta_1s*sin psi`.
    pub col: f64,
    pub omega: f64,
    pub omega_r: f64,
    pub rho: f64,
    pub n_b: usize,
    pub n_psi: usize,
    /// In-plane wind in hub frame; `v_t_extra = v_in_hub_x*sin psi + v_in_hub_y*cos psi`.
    pub v_in_hub_x: f64,
    pub v_in_hub_y: f64,
    pub theta_1c: f64,
    pub theta_1s: f64,
}

impl<'a> SweepCtx<'a> {
    /// Run one full psi x radial sweep with the given kernel. Returns the
    /// azimuth-averaged (T, Q, Mx_hub, My_hub) over the rotor disk.
    ///
    /// `self.omega > 0` is assumed (caller filters out the not-spinning case
    /// before invoking). Reverse-flow region (`v_t <= 0`) is skipped
    /// per-element.
    #[inline(always)]
    pub fn run<K: PsiKernel>(&self, kernel: &mut K) -> (f64, f64, f64, f64) {
        let mut t_acc = 0.0;
        let mut q_acc = 0.0;
        let mut mx_acc = 0.0;
        let mut my_acc = 0.0;
        let inv_n_psi = 1.0 / (self.n_psi as f64);
        let grid = self.grid;
        let n_r = grid.r_mid.len();
        for i_psi in 0..self.n_psi {
            let psi = 2.0 * PI * (i_psi as f64) * inv_n_psi;
            let cos_psi = psi.cos();
            let sin_psi = psi.sin();
            let v_t_extra = self.v_in_hub_x * sin_psi + self.v_in_hub_y * cos_psi;
            let col_psi = self.col + self.theta_1c * cos_psi + self.theta_1s * sin_psi;
            let mut rdt_sum = 0.0;
            for i in 0..n_r {
                let r = grid.r_mid[i];
                let v_t = self.omega * r + v_t_extra;
                if v_t <= 0.0 {
                    continue;
                }
                let ctx = ElementCtx {
                    i,
                    cos_psi,
                    sin_psi,
                    r,
                    chord: grid.chord[i],
                    twist: grid.twist_rad[i],
                    dr: grid.dr,
                    col_psi,
                    v_t,
                };
                let (dt, dq) = kernel.element(self, &ctx);
                t_acc += dt;
                q_acc += dq;
                rdt_sum += r * dt;
            }
            mx_acc += rdt_sum * sin_psi;
            my_acc += rdt_sum * cos_psi;
        }
        (
            t_acc * inv_n_psi,
            q_acc * inv_n_psi,
            mx_acc * inv_n_psi,
            my_acc * inv_n_psi,
        )
    }
}
