// Shared infrastructure for the BEM models.
// See ../CLAUDE.md "Shared BEM infrastructure".

use crate::aero_io::{AeroResult, Mat3, RotorInputs, Vec3};
use crate::common::{
    AlignedBEMF64, EPS_OMEGA_R, MASS_FLOW_HOVER_FLOOR_FRAC, MAX_BEM_ELEMENTS, VRS_DESCENT_THRESHOLD,
};
use crate::polar::Polar;
use crate::rotor_definition::BladeGeometry;
use std::f64::consts::PI;

/// Build a fixed trig table for a given azimuth resolution.
///
/// Entry i stores (cos(psi_i), sin(psi_i)) where
/// psi_i = 2*pi*i/n_psi.
pub fn build_psi_trig_table(n_psi: usize) -> Vec<(f64, f64)> {
    assert!(n_psi > 0);

    let inv_n_psi = 1.0 / (n_psi as f64);
    let mut out = Vec::with_capacity(n_psi);
    for i in 0..n_psi {
        let psi = 2.0 * PI * (i as f64) * inv_n_psi;
        out.push((psi.cos(), psi.sin()));
    }
    out
}

/// Cached fixed radial geometry for a BEM kernel.
///
/// Per-station chord/twist arrays mean we can support BladeGeometry with
/// radial-station arrays (wind-turbine blades) AND scalar chord_m/twist_deg
/// (helicopters) through a single uniform interface in the inner loop --
/// no branches per element.
#[derive(Clone, Debug)]
pub struct RadialGrid {
    pub n_elements: usize,
    pub dr: f64,
    pub r_mid: AlignedBEMF64,     // n
    pub x_mid: AlignedBEMF64,     // n = r_mid / R
    pub x_hub: f64,               // root_cutout / R
    pub chord: AlignedBEMF64,     // n  -- per-station chord (m)
    pub twist_rad: AlignedBEMF64, // n  -- per-station twist (rad)
}

impl RadialGrid {
    pub fn from_blade(blade: &BladeGeometry) -> Self {
        let r_root = blade.root_cutout_m;
        let r_tip = blade.radius_m;
        let n = blade.n_elements;
        assert!(n <= MAX_BEM_ELEMENTS);
        let dr = (r_tip - r_root) / (n as f64);
        let mut r_mid = aligned::Aligned([0.0; MAX_BEM_ELEMENTS]);
        let mut x_mid = aligned::Aligned([0.0; MAX_BEM_ELEMENTS]);
        let mut chord = aligned::Aligned([0.0; MAX_BEM_ELEMENTS]);
        let mut twist_rad = aligned::Aligned([0.0; MAX_BEM_ELEMENTS]);
        for i in 0..n {
            let r = r_root + (i as f64 + 0.5) * dr;
            r_mid[i] = r;
            x_mid[i] = if r_tip > 0.0 { r / r_tip } else { 0.0 };
            chord[i] = blade.chord_at(r);
            twist_rad[i] = blade.twist_at(r).to_radians();
        }
        let x_hub = if r_tip > 0.0 { r_root / r_tip } else { 0.0 };
        Self {
            n_elements: n,
            dr,
            r_mid,
            x_mid,
            x_hub,
            chord,
            twist_rad,
        }
    }
}

/// Tabulate any polar onto contiguous arrays for the fast inner
/// loop. TabulatedPolar passes its arrays through; analytical polars get
/// sampled to 4001 points over [-pi/2, pi/2] (matching the Python version).
#[derive(Clone, Debug)]
pub struct PolarTable {
    pub alpha: Vec<f64>,
    pub cl: Vec<f64>,
    pub cd: Vec<f64>,
}

impl PolarTable {
    pub fn from_polar<P: Polar>(polar: &P) -> Self {
        if let Some((alpha, cl, cd)) = polar.table_data() {
            return Self {
                alpha: alpha.to_vec(),
                cl: cl.to_vec(),
                cd: cd.to_vec(),
            };
        }
        // Analytical polar: sample 4001 points over [-pi/2, pi/2].
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

/// Mass-flow speed at the disk (the Glauert mass-flow scalar V_mf).
///
/// This is the resultant flow speed *through* the rotor disk used to scale
/// the dynamic-inflow time constants -- NOT the blade tip speed `omega_r`
/// and NOT the blade-element tangential velocity `v_t`.
///
/// ```text
/// V_mf = sqrt(v_edge^2 + (v_climb + v0_axial)^2)
/// ```
///
/// `v0_axial` is the axial (induced) component in m/s. Floored at
/// `MASS_FLOW_HOVER_FLOOR_FRAC * max(omega_r, 1)` to keep the Pitt-Peters / Oye
/// time constants finite at hover / zero thrust.
#[inline]
pub fn v_mass_flow_disk(v_edge: f64, v_climb: f64, v0_axial: f64, omega_r: f64) -> f64 {
    (v_edge * v_edge + (v_climb + v0_axial).powi(2))
        .sqrt()
        .max(MASS_FLOW_HOVER_FLOOR_FRAC * omega_r.max(1.0))
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

/// Result of the 1/rev blade-flapping solve.
///
/// `mx_out`/`my_out` are the hub moments that reach the airframe (the flap
/// deflection redistributes the aerodynamic reaction between blade
/// centrifugal restoring and hub structure). `dfx_hub`/`dfy_hub` are the
/// additional in-plane hub force (the flapping-tilt contribution to
/// H-force): as the blade flaps by beta(psi) the thrust vector tilts into
/// the disk plane. Both are in the hub frame and are added to the
/// profile-drag H-force before `assemble_result`.
pub struct FlapOutput {
    pub mx_out: f64,
    pub my_out: f64,
    pub dfx_hub: f64,
    pub dfy_hub: f64,
}

/// Phase-correct quasi-static blade-flapping solve (1/rev harmonic balance
/// with aerodynamic flap damping).
///
/// The flap DOF responds to the RAW aerodynamic hub moment (the caller must
/// pass the un-reduced `mx_hub`/`my_hub`; the inflow ODE still uses those
/// full moments because the wake responds to disk loading, not to what the
/// airframe sees). The solve returns both the reduced hub moments and the
/// flapping-tilt H-force.
///
/// ## Convention (see AGENTS.md "blade flapping" section)
///
/// - Azimuth psi=0 at +X, CCW from above; blade flap `beta > 0` = tip up
///   (toward -z, the thrust direction). `beta(psi) = beta_0 + beta_1c cos
///   psi + beta_1s sin psi`.
/// - The sweep's azimuth-averaged hub moments are the 1/rev aerodynamic
///   flap-moment harmonics: `mx_hub = N_b M1s / 2`, `my_hub = N_b M1c / 2`.
/// - Flap equation (per rev, ' = d/dpsi):
///   `I_b Omega^2 (beta'' + d beta' + nu^2 beta) = M_beta(psi)` with
///   aerodynamic damping `d = gamma/8`, Lock number `gamma = rho a c R^4 /
///   I_b`. Computed here from the actual radial geometry to allow taper:
///   `d = 0.5 rho a S3 / I_b`, `S3 = sum c(r) r^3 dr` (Omega cancels).
/// - 1/rev balance gives the damped transfer
///   `[[a, d],[-d, a]] [beta_1c; beta_1s] = (1/(I_b Omega^2)) [M1c; M1s]`,
///   `a = nu^2 - 1`. At `nu ~ 1` the damping dominates and produces the
///   classical ~90 deg flap lag (flapback) rather than an in-phase response.
/// - Transmitted hub moment (Johnson): `M_hub = (N_b/2) I_b Omega^2 (nu^2-1)
///   beta_1`. With no damping this reduces to the full aero moment.
/// - Flapping-tilt H-force: `Fx = -<T beta cos psi> = -T beta_1c/2`,
///   `Fy = +<T beta sin psi> = +T beta_1s/2` (T = mean disk thrust).
///
/// The H-force sign is pinned by `validation_rs/src/checks/h_force.rs`
/// (forward flight with collective must give a rearward, flow-opposing
/// force that adds to the profile-drag term).
#[inline]
pub fn apply_flap_dynamics(
    t_total: f64,
    mx_hub: f64,
    my_hub: f64,
    flap: Option<&crate::rotor_definition::FlapProperties>,
    grid: &RadialGrid,
    cl_alpha: f64,
    rho: f64,
    n_b: usize,
    omega: f64,
) -> FlapOutput {
    let fp = match flap {
        Some(fp) => fp,
        None => {
            return FlapOutput {
                mx_out: mx_hub,
                my_out: my_hub,
                dfx_hub: 0.0,
                dfy_hub: 0.0,
            }
        }
    };

    let i_b = fp.I_blade_flap_kgm2;
    if !(i_b > 0.0) || omega.abs() < 1e-6 {
        // Degenerate: no inertia or not spinning -> no meaningful flap solve.
        return FlapOutput {
            mx_out: mx_hub,
            my_out: my_hub,
            dfx_hub: 0.0,
            dfy_hub: 0.0,
        };
    }

    let nu2 = fp.nu_beta_sq(omega);
    let a_stiff = nu2 - 1.0;

    // Aerodynamic flap damping (nondimensional, = gamma/8). Integrated over
    // the radial grid so taper is handled: C_beta = 0.5 rho a Omega * S3,
    // S3 = sum c(r) r^3 dr, and gamma/8 = C_beta/(I_b Omega) = 0.5 rho a
    // S3 / I_b (Omega cancels).
    let mut s3 = 0.0;
    for i in 0..grid.n_elements {
        let r = grid.r_mid[i];
        s3 += grid.chord[i] * r * r * r * grid.dr;
    }
    let d_damp = 0.5 * rho * cl_alpha * s3 / i_b;

    // 1/rev aerodynamic flap-moment harmonics (per blade).
    let n_bf = n_b as f64;
    let m1c = 2.0 * my_hub / n_bf;
    let m1s = 2.0 * mx_hub / n_bf;

    // Damped 1/rev flap response:
    //   [[a, d],[-d, a]] [b1c; b1s] = (1/(I_b Omega^2)) [M1c; M1s]
    // => [b1c; b1s] = (1/(I_b Omega^2 D)) [[a, -d],[d, a]] [M1c; M1s]
    let om2 = omega * omega;
    let det = a_stiff * a_stiff + d_damp * d_damp;
    let inv = 1.0 / (i_b * om2 * det);
    let b1c = inv * (a_stiff * m1c - d_damp * m1s);
    let b1s = inv * (d_damp * m1c + a_stiff * m1s);

    // Transmitted hub moment (Johnson): M_hub = (N_b/2) I_b Omega^2 (nu^2-1) beta1.
    let k = 0.5 * n_bf * i_b * om2 * a_stiff;
    let my_out = k * b1c;
    let mx_out = k * b1s;

    // Flapping-tilt H-force from the thrust vector tilting with beta(psi).
    let dfx_hub = -0.5 * t_total * b1c;
    let dfy_hub = 0.5 * t_total * b1s;

    FlapOutput {
        mx_out,
        my_out,
        dfx_hub,
        dfy_hub,
    }
}

///
/// `fx_hub, fy_hub` are the net in-plane (H-force) components in the hub
/// frame -- the world-frame reaction of the blades' tangential aerodynamic
/// force, summed over the disk (see `SweepCtx::run`'s `fx_acc`/`fy_acc`).
/// This is zero whenever the tangential loading is azimuth-independent
/// (pure axial flow / hover) and grows with edgewise flow (`v_edge`),
/// vanishing again as the hub axis re-aligns with the relative wind.
/// Pass `0.0, 0.0` for callers that don't yet compute it (e.g. VPM).
///
/// The caller is expected to have already folded in the flapping-tilt
/// contribution to H-force (`FlapOutput::dfx_hub`/`dfy_hub` from
/// `apply_flap_dynamics`) when a `FlapProperties` is configured; the
/// profile/induced-drag part comes from `SweepCtx::run`.
#[inline]
pub fn assemble_result(
    t_total: f64,
    q_total: f64,
    mx_hub: f64,
    my_hub: f64,
    fx_hub: f64,
    fy_hub: f64,
    hub_axis: Vec3,
    r_hub: &Mat3,
) -> AeroResult {
    let f_hub = Vec3::new(fx_hub, fy_hub, -t_total);
    let f_world = *r_hub * f_hub;
    let mxyz_hub = Vec3::new(mx_hub, my_hub, 0.0);
    let m_orbital = *r_hub * mxyz_hub;
    let m_spin = hub_axis * q_total;
    AeroResult {
        F_world: f_world,
        M_hub_world: m_orbital,
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
pub fn element_force<P: Polar>(v_a: f64, sweep: &SweepCtx<'_, P>, ctx: &ElementCtx) -> (f64, f64) {
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
// One override point:
//   - element(ctx): compute (dT, dQ) for one element.
//     Pitt-Peters and Oye both use the prescribed-inflow path
//     (compute lam_local, then call element_force).
//     BEM overrides this to call its iterative solve_bem_element instead.
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
    /// Element-level force computation.
    ///
    /// Models with prescribed inflow (Pitt-Peters, Oye) compute `lam`
    /// directly and then call `element_force`; BEM-style models can run
    /// their own per-element solver here.
    fn element<P: Polar>(&mut self, sweep: &SweepCtx<'_, P>, ctx: &ElementCtx) -> (f64, f64);
}

/// Whole-sweep configuration: the call-invariant inputs that describe one
/// psi x radial pass over the rotor disk. Built once by each model's
/// compute_forces, then handed to run_psi_loop, which derives per-element
/// `ElementCtx` values from these fields each iteration. Symmetric with
/// `ElementCtx` (which describes one element rather than one sweep).
pub struct SweepCtx<'a, P: Polar> {
    pub grid: &'a RadialGrid,
    pub polar: &'a P,
    /// Base collective pitch (rad). Per-azimuth pitch is `col + theta_1c*cos psi + theta_1s*sin psi`.
    pub col: f64,
    pub omega: f64,
    pub omega_r: f64,
    pub rho: f64,
    pub n_b: usize,
    pub n_psi: usize,
    pub n_psi_inv: f64,
    pub psi_trig: &'a [(f64, f64)],
    /// In-plane wind in hub frame; `v_t_extra = v_in_hub_x*sin psi + v_in_hub_y*cos psi`.
    pub v_in_hub_x: f64,
    pub v_in_hub_y: f64,
    pub theta_1c: f64,
    pub theta_1s: f64,
}

impl<'a, P: Polar> SweepCtx<'a, P> {
    /// Run one full psi x radial sweep with the given kernel. Returns the
    /// azimuth-averaged (T, Q, Mx_hub, My_hub, Fx_hub, Fy_hub) over the rotor
    /// disk.
    ///
    /// `Fx_hub`/`Fy_hub` are the net in-plane hub force (H-force): each
    /// element's tangential aerodynamic force is `dFt = dQ / r` (the same
    /// force whose moment produces `dQ`; recovered by undoing the `* r`
    /// weighting rather than recomputing it), projected onto the fixed hub
    /// x/y axes with the same `(sin psi, cos psi)` pairing used everywhere
    /// else in this sweep (`v_t_extra = v_in_hub_x*sin psi + v_in_hub_y*cos
    /// psi`, `Mx_hub`/`My_hub` below). This is zero for azimuth-independent
    /// (pure axial) loading and grows with edgewise flow -- see
    /// `assemble_result`.
    ///
    /// `self.omega > 0` is assumed (caller filters out the not-spinning case
    /// before invoking). Reverse-flow region (`v_t <= 0`) is skipped
    /// per-element.
    #[inline(always)]
    pub fn run<K: PsiKernel>(&self, kernel: &mut K) -> (f64, f64, f64, f64, f64, f64) {
        let mut t_acc = 0.0;
        let mut q_acc = 0.0;
        let mut mx_acc = 0.0;
        let mut my_acc = 0.0;
        let mut fx_acc = 0.0;
        let mut fy_acc = 0.0;
        let inv_n_psi = self.n_psi_inv;
        let grid = self.grid;
        let n_r = grid.n_elements;
        let n_psi = self.n_psi;

        assert!(n_r < MAX_BEM_ELEMENTS);
        assert_eq!(self.psi_trig.len(), n_psi);

        let r_mid = &grid.r_mid;
        let chord = &grid.chord;
        let twist = &grid.twist_rad;

        for i_psi in 0..n_psi {
            let (cos_psi, sin_psi) = self.psi_trig[i_psi];
            let v_t_extra = self.v_in_hub_x * sin_psi + self.v_in_hub_y * cos_psi;
            let col_psi = self.col + self.theta_1c * cos_psi + self.theta_1s * sin_psi;
            let mut rdt_sum = 0.0;
            for i in 0..n_r {
                let r = r_mid[i];
                let v_t = self.omega * r + v_t_extra;
                if v_t <= 0.0 {
                    continue;
                }
                let ctx = ElementCtx {
                    i,
                    cos_psi,
                    sin_psi,
                    r,
                    chord: chord[i],
                    twist: twist[i],
                    dr: grid.dr,
                    col_psi,
                    v_t,
                };
                let (dt, dq) = kernel.element(self, &ctx);
                t_acc += dt;
                q_acc += dq;
                rdt_sum += r * dt;
                let d_ft = dq / r;
                fx_acc += d_ft * sin_psi;
                fy_acc += d_ft * cos_psi;
            }
            mx_acc += rdt_sum * sin_psi;
            my_acc += rdt_sum * cos_psi;
        }
        (
            t_acc * inv_n_psi,
            q_acc * inv_n_psi,
            mx_acc * inv_n_psi,
            my_acc * inv_n_psi,
            fx_acc * inv_n_psi,
            fy_acc * inv_n_psi,
        )
    }
}
