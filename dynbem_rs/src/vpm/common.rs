// Vortex Particle Method (VPM) core -- free-wake particle engine.
//
// This is the first, deliberately un-abstracted implementation: a
// self-contained particle field, a regularized Biot-Savart velocity
// evaluator (direct O(N^2), SIMD via the `wide` crate), and a free-wake
// convection step. No rotor coupling, no vortex stretching, no viscous
// diffusion yet -- those are the documented next steps below. Getting the
// convection kinematics correct and fast is the load-bearing part, so it
// comes first and is validated on its own (see the vortex-ring test).
//
// ---------------------------------------------------------------------------
// Standard formulation (following current VPM literature)
// ---------------------------------------------------------------------------
//
// The wake vorticity field is discretized as a sum of regularized particles
//
//     omega(x) = sum_p  alpha_p  zeta_sigma(x - x_p)
//
// where `alpha_p` is the particle's vector-valued strength (the integral of
// vorticity it carries, units m^3/s = circulation * length) and
// `zeta_sigma` is a radially-symmetric regularization (smoothing) kernel of
// core size `sigma_p`.
//
// The induced velocity is the regularized Biot-Savart law
//
//     u(x) = (1 / 4 pi) sum_p  g(rho_p)  [ alpha_p x (x - x_p) ] / r_p^3
//
// with r_p = |x - x_p|, rho_p = r_p / sigma_p, and `g` the low-order
// algebraic regularization function of Winckelmans & Leonard (1993):
//
//     g(rho) = rho^3 (rho^2 + 5/2) / (rho^2 + 1)^(5/2)
//
// (this is the smoothing paired with zeta(rho) = (15 / 8 pi) (rho^2 + 1)^(-7/2);
// the same kernel used by modern rotorcraft VPM codes, e.g. FLOWVPM /
// Alvarez & Ning 2020). As rho -> infinity, g -> 1 and the singular
// Biot-Savart law is recovered; as rho -> 0, g ~ (5/2) rho^3 so the velocity
// stays finite.
//
// Singularity-free rearrangement (what the code actually evaluates)
// -----------------------------------------------------------------
// We never divide by r. Writing rho^2 = r^2 / sigma^2 and folding
// rho^3 = r^3 / sigma^3 into g / r^3 gives a kernel that depends only on
// rho^2 and sigma:
//
//     K(rho, sigma) = g(rho) / r^3
//                   = (rho^2 + 5/2) / [ sigma^3 (rho^2 + 1)^(5/2) ]
//
//     u(x) = (1 / 4 pi) sum_p  K(rho_p, sigma_p)  [ alpha_p x (x - x_p) ]
//
// K is finite for all r (K(0, sigma) = (5/2) / sigma^3), and the self term
// (x == x_p) contributes exactly zero because the cross product
// alpha_p x (x - x_p) vanishes there. So there are no branches, no masks,
// and no division by r in the hot loop -- ideal for SIMD.
//
// Sign check: a particle strength along +z, alpha = (0, 0, A) with A > 0,
// induces velocity in +y at a point on the +x axis, i.e. counter-clockwise
// circulation about +z -- the right-hand rule for vorticity along +z.
//
// ---------------------------------------------------------------------------
// NOT yet implemented (next steps, in rough priority order)
// ---------------------------------------------------------------------------
//   1. Vortex stretching  d alpha / dt = (alpha . grad) u  -- needed for 3D
//      circulation conservation once the wake distorts. The algebraic kernel
//      has an analytic velocity gradient; add it as a second accumulator.
//   2. Viscous diffusion (core spreading or PSE) so long-lived hover/VRS
//      wakes don't stay artificially coherent.
//   3. Rotor coupling: shed particles from blade-element bound circulation
//      each azimuth step, and feed the induced velocity at the disk back into
//      the blade-element loads.
//
// Done since first draft: the rotor coupling (this module's parent, `vpm`) and a
// Barnes-Hut O(N log N) evaluator (`induced_at_points_bh` / `advect_rk2_bh`)
// for when N outgrows the direct path.

use aligned_vec::AVec;
use bytemuck::cast_slice;
use rayon::prelude::*;
use wide::f32x8;

/// 1 / (4 pi), the Biot-Savart prefactor.
const INV_4PI: f32 = 0.079_577_47;
const INV_4PI_F64: f64 = 0.079_577_471_545_947_67;

/// Memory alignment (bytes) for the SoA particle arrays. A cache line (64)
/// also satisfies the 32-byte alignment `wide::f32x8` needs for aligned
/// vector loads, so the SIMD loops can read each eight-wide chunk with an
/// aligned move.
const SIMD_ALIGN: usize = 64;

/// A cloud of regularized vortex particles in structure-of-arrays layout.
///
/// SoA (not array-of-structs) is deliberate: it lets each SIMD lane load
/// eight contiguous particle components at once. All quantities are `f32`
/// -- the wake spans meters and the velocity evaluation is compute-bound,
/// so single precision is the standard performance choice for VPM. The
/// backing storage is [`AVec`] aligned to [`SIMD_ALIGN`] so the eight-wide
/// reads are aligned vector loads.
#[derive(Clone, Debug)]
pub struct ParticleField {
    /// Positions.
    pub px: AVec<f32>,
    pub py: AVec<f32>,
    pub pz: AVec<f32>,
    /// Vector strengths alpha (integral of vorticity), units m^3/s.
    pub ax: AVec<f32>,
    pub ay: AVec<f32>,
    pub az: AVec<f32>,
    /// Per-particle regularization core size sigma (m).
    pub sigma: AVec<f32>,
}

impl Default for ParticleField {
    fn default() -> Self {
        Self {
            px: AVec::new(SIMD_ALIGN),
            py: AVec::new(SIMD_ALIGN),
            pz: AVec::new(SIMD_ALIGN),
            ax: AVec::new(SIMD_ALIGN),
            ay: AVec::new(SIMD_ALIGN),
            az: AVec::new(SIMD_ALIGN),
            sigma: AVec::new(SIMD_ALIGN),
        }
    }
}

impl ParticleField {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_capacity(n: usize) -> Self {
        Self {
            px: AVec::with_capacity(SIMD_ALIGN, n),
            py: AVec::with_capacity(SIMD_ALIGN, n),
            pz: AVec::with_capacity(SIMD_ALIGN, n),
            ax: AVec::with_capacity(SIMD_ALIGN, n),
            ay: AVec::with_capacity(SIMD_ALIGN, n),
            az: AVec::with_capacity(SIMD_ALIGN, n),
            sigma: AVec::with_capacity(SIMD_ALIGN, n),
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.px.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.px.is_empty()
    }

    /// Append one particle. `sigma` must be > 0.
    #[inline]
    pub fn push(&mut self, pos: [f32; 3], strength: [f32; 3], sigma: f32) {
        debug_assert!(sigma > 0.0, "particle core size sigma must be positive");
        self.px.push(pos[0]);
        self.py.push(pos[1]);
        self.pz.push(pos[2]);
        self.ax.push(strength[0]);
        self.ay.push(strength[1]);
        self.az.push(strength[2]);
        self.sigma.push(sigma);
    }

    pub fn clear(&mut self) {
        self.px.clear();
        self.py.clear();
        self.pz.clear();
        self.ax.clear();
        self.ay.clear();
        self.az.clear();
        self.sigma.clear();
    }

    /// Iterate over all particles as `(pos, strength, sigma)` tuples.
    /// `pos` and `strength` are `[f32; 3]`; `sigma` is the core radius in metres.
    pub fn particles(&self) -> impl Iterator<Item = ([f32; 3], [f32; 3], f32)> + '_ {
        (0..self.len()).map(move |i| {
            (
                [self.px[i], self.py[i], self.pz[i]],
                [self.ax[i], self.ay[i], self.az[i]],
                self.sigma[i],
            )
        })
    }

    /// Remove the oldest `k` particles from the front of the SoA arrays.
    ///
    /// `AVec` does not implement `drain`, so we rotate the survivors to the
    /// front (an in-place O(n) slice operation) and then truncate.
    pub fn drain_front(&mut self, k: usize) {
        let n = self.px.len();
        if k == 0 || k > n {
            return;
        }
        let keep = n - k;
        macro_rules! rotate_trunc {
            ($field:expr) => {
                $field.rotate_left(k);
                $field.truncate(keep);
            };
        }
        rotate_trunc!(self.px);
        rotate_trunc!(self.py);
        rotate_trunc!(self.pz);
        rotate_trunc!(self.ax);
        rotate_trunc!(self.ay);
        rotate_trunc!(self.az);
        rotate_trunc!(self.sigma);
    }
}

/// Regularized-Biot-Savart velocity induced by the whole field at every
/// particle location. Direct O(N^2) evaluation, SIMD-accelerated over the
/// source particles (eight at a time via `wide::f32x8`).
///
/// The self term contributes exactly zero (the cross product vanishes), so
/// `out[j]` is the velocity induced on particle `j` by all *other*
/// particles -- i.e. the free-wake convection velocity of particle `j`.
///
/// Returns one `[vx, vy, vz]` per particle, in field order.
pub fn induced_velocities(field: &ParticleField) -> Vec<[f32; 3]> {
    induced_at_points(field, &field.px, &field.py, &field.pz)
}

/// Evaluate the Biot-Savart kernel at a single target point given the
/// pre-built `Chunks` view of the source field. Shared by the sequential
/// and parallel paths in [`induced_at_points`].
#[inline]
fn eval_target(xj: f32, yj: f32, zj: f32, src: &Chunks<'_>, n_chunks: usize) -> [f32; 3] {
    let xjv = f32x8::splat(xj);
    let yjv = f32x8::splat(yj);
    let zjv = f32x8::splat(zj);
    let mut ux = f32x8::splat(0.0);
    let mut uy = f32x8::splat(0.0);
    let mut uz = f32x8::splat(0.0);
    accumulate_chunks(xjv, yjv, zjv, src, 0, n_chunks, &mut ux, &mut uy, &mut uz);
    reduce_velocity(ux, uy, uz)
}

/// Regularized-Biot-Savart velocity induced by `field` at a set of arbitrary
/// probe points (e.g. blade-element control points). Same kernel and
/// singularity-free formulation as [`induced_velocities`]; a probe point that
/// coincides with a particle still gets a finite, self-free contribution
/// because the cross product vanishes there.
///
/// When the `parallel` crate feature is enabled the outer target loop runs on
/// the Rayon thread pool (one target per task, sources vectorized eight-wide
/// per lane as in the sequential path). The number of threads is controlled by
/// Rayon's global pool (defaults to the logical core count).
pub fn induced_at_points(
    field: &ParticleField,
    tx: &[f32],
    ty: &[f32],
    tz: &[f32],
) -> Vec<[f32; 3]> {
    let n = field.len();
    let m = tx.len();
    let mut out = vec![[0.0f32; 3]; m];
    if n == 0 || m == 0 {
        return out;
    }

    // Pad the source arrays up to a multiple of 8 so the SIMD loop has no
    // remainder. Padding particles carry zero strength (zero contribution)
    // and sigma = 1 (keeps sigma^3 finite; never actually used). The buffers
    // are aligned so they can be reinterpreted as `&[f32x8]` and read with
    // aligned vector loads.
    let n_pad = n.div_ceil(8) * 8;
    let pad = |src: &[f32], fill: f32| -> AVec<f32> {
        let mut v = AVec::<f32>::with_capacity(SIMD_ALIGN, n_pad);
        v.extend_from_slice(src);
        v.resize(n_pad, fill);
        v
    };
    let px = pad(&field.px, 0.0);
    let py = pad(&field.py, 0.0);
    let pz = pad(&field.pz, 0.0);
    let ax = pad(&field.ax, 0.0);
    let ay = pad(&field.ay, 0.0);
    let az = pad(&field.az, 0.0);
    let sg = pad(&field.sigma, 1.0);

    // Reinterpret each aligned, multiple-of-8 array as a slice of eight-wide
    // vectors; indexing these is an aligned load.
    let src = Chunks::from_avecs(&px, &py, &pz, &ax, &ay, &az, &sg);
    let n_chunks = n_pad / 8;

    out.par_iter_mut().enumerate().for_each(|(j, o)| {
        *o = eval_target(tx[j], ty[j], tz[j], &src, n_chunks);
    });

    out
}

/// Sequential version of [`induced_at_points`] -- always single-threaded.
/// Useful for benchmarking or when the caller manages its own parallelism
/// at a higher level.
pub fn induced_at_points_seq(
    field: &ParticleField,
    tx: &[f32],
    ty: &[f32],
    tz: &[f32],
) -> Vec<[f32; 3]> {
    let n = field.len();
    let m = tx.len();
    let mut out = vec![[0.0f32; 3]; m];
    if n == 0 || m == 0 {
        return out;
    }
    let n_pad = n.div_ceil(8) * 8;
    let pad = |src: &[f32], fill: f32| -> AVec<f32> {
        let mut v = AVec::<f32>::with_capacity(SIMD_ALIGN, n_pad);
        v.extend_from_slice(src);
        v.resize(n_pad, fill);
        v
    };
    let px = pad(&field.px, 0.0);
    let py = pad(&field.py, 0.0);
    let pz = pad(&field.pz, 0.0);
    let ax = pad(&field.ax, 0.0);
    let ay = pad(&field.ay, 0.0);
    let az = pad(&field.az, 0.0);
    let sg = pad(&field.sigma, 1.0);
    let src = Chunks::from_avecs(&px, &py, &pz, &ax, &ay, &az, &sg);
    let n_chunks = n_pad / 8;
    for j in 0..m {
        out[j] = eval_target(tx[j], ty[j], tz[j], &src, n_chunks);
    }
    out
}

/// Sequential self-evaluation variant (wake on itself). Always single-threaded.
pub fn induced_velocities_seq(field: &ParticleField) -> Vec<[f32; 3]> {
    induced_at_points_seq(field, &field.px, &field.py, &field.pz)
}

/// Scalar `f64` reference implementation of [`induced_velocities`].
///
/// Same formulation, evaluated pair-by-pair with double-precision
/// accumulation. This is the readable specification and the ground truth the
/// SIMD path is validated against; it is not the production hot path.
#[allow(dead_code)] // reference/spec impl; used by the SIMD-vs-scalar tests
pub fn induced_velocities_ref(field: &ParticleField) -> Vec<[f64; 3]> {
    let n = field.len();
    let mut out = vec![[0.0f64; 3]; n];
    for j in 0..n {
        let xj = field.px[j] as f64;
        let yj = field.py[j] as f64;
        let zj = field.pz[j] as f64;
        let mut ux = 0.0f64;
        let mut uy = 0.0f64;
        let mut uz = 0.0f64;
        for s in 0..n {
            let dx = xj - field.px[s] as f64;
            let dy = yj - field.py[s] as f64;
            let dz = zj - field.pz[s] as f64;
            let sax = field.ax[s] as f64;
            let say = field.ay[s] as f64;
            let saz = field.az[s] as f64;
            let sigma = field.sigma[s] as f64;

            let r2 = dx * dx + dy * dy + dz * dz;
            let sigma2 = sigma * sigma;
            let sigma3 = sigma2 * sigma;
            let rho2 = r2 / sigma2;
            let base = rho2 + 1.0;
            let denom = base * base * base.sqrt();
            let k = (rho2 + 2.5) / (sigma3 * denom);

            ux += k * (say * dz - saz * dy);
            uy += k * (saz * dx - sax * dz);
            uz += k * (sax * dy - say * dx);
        }
        out[j] = [ux * INV_4PI_F64, uy * INV_4PI_F64, uz * INV_4PI_F64];
    }
    out
}

/// Scalar f64 reference implementation of [`induced_at_points`] with NaN
/// assertions. Evaluates each target against all sources pair-by-pair in f64
/// and panics if any output is non-finite, printing which target index and
/// which source index (if detectable) caused the issue.
///
/// Intended for debugging only -- O(N^2) and non-vectorized. Enable via
/// `VpmRotorConfig::use_scalar_nan_check`.
pub fn induced_at_points_nan_check(
    field: &ParticleField,
    tx: &[f32],
    ty: &[f32],
    tz: &[f32],
) -> Vec<[f32; 3]> {
    let n = field.len();
    let m = tx.len();
    let mut out = vec![[0.0f32; 3]; m];
    if n == 0 || m == 0 {
        return out;
    }

    // Sanity-check the source field first.
    for s in 0..n {
        let sigma = field.sigma[s] as f64;
        assert!(
            sigma > 0.0 && sigma.is_finite(),
            "nan_check: source particle {} has bad sigma={}",
            s,
            sigma
        );
        let pos = [field.px[s] as f64, field.py[s] as f64, field.pz[s] as f64];
        let str_ = [field.ax[s] as f64, field.ay[s] as f64, field.az[s] as f64];
        for k in 0..3 {
            assert!(
                pos[k].is_finite(),
                "nan_check: source particle {} pos[{}]={} is non-finite",
                s,
                k,
                pos[k]
            );
            assert!(
                str_[k].is_finite(),
                "nan_check: source particle {} strength[{}]={} is non-finite",
                s,
                k,
                str_[k]
            );
        }
    }

    for j in 0..m {
        let xj = tx[j] as f64;
        let yj = ty[j] as f64;
        let zj = tz[j] as f64;
        let mut ux = 0.0f64;
        let mut uy = 0.0f64;
        let mut uz = 0.0f64;
        for s in 0..n {
            let dx = xj - field.px[s] as f64;
            let dy = yj - field.py[s] as f64;
            let dz = zj - field.pz[s] as f64;
            let sax = field.ax[s] as f64;
            let say = field.ay[s] as f64;
            let saz = field.az[s] as f64;
            let sigma = field.sigma[s] as f64;

            let r2 = dx * dx + dy * dy + dz * dz;
            let sigma2 = sigma * sigma;
            let sigma3 = sigma2 * sigma;
            let rho2 = r2 / sigma2;
            let base = rho2 + 1.0;
            let denom = base * base * base.sqrt();
            let k = (rho2 + 2.5) / (sigma3 * denom);

            let dux = k * (say * dz - saz * dy);
            let duy = k * (saz * dx - sax * dz);
            let duz = k * (sax * dy - say * dx);

            assert!(
                dux.is_finite() && duy.is_finite() && duz.is_finite(),
                "nan_check: NaN at target j={} from source s={}: \
                 pos_j=[{},{},{}] pos_s=[{},{},{}] alpha_s=[{},{},{}] sigma={} \
                 r2={} rho2={} denom={} k={} du=[{},{},{}]",
                j,
                s,
                xj,
                yj,
                zj,
                field.px[s],
                field.py[s],
                field.pz[s],
                sax,
                say,
                saz,
                sigma,
                r2,
                rho2,
                denom,
                k,
                dux,
                duy,
                duz
            );

            ux += dux;
            uy += duy;
            uz += duz;
        }
        let vx = (ux * INV_4PI_F64) as f32;
        let vy = (uy * INV_4PI_F64) as f32;
        let vz = (uz * INV_4PI_F64) as f32;
        assert!(
            vx.is_finite() && vy.is_finite() && vz.is_finite(),
            "nan_check: accumulated NaN at target j={}: v=[{},{},{}]",
            j,
            vx,
            vy,
            vz
        );
        out[j] = [vx, vy, vz];
    }
    out
}

/// Shared RK2 midpoint free-wake step. `eval` supplies the induced velocity at
/// every particle (direct or Barnes-Hut); the integration is otherwise
/// identical. Monomorphized per call site, so the closure inlines away.
#[inline]
fn advect_rk2_with(
    field: &mut ParticleField,
    freestream: [f32; 3],
    dt: f32,
    mut eval: impl FnMut(&ParticleField) -> Vec<[f32; 3]>,
) {
    let n = field.len();
    if n == 0 {
        return;
    }

    // Stage 1: velocity at the current positions.
    let u1 = eval(field);

    // Midpoint positions: x + 0.5 dt (u1 + freestream).
    let half = 0.5 * dt;
    let mut mid = field.clone();
    for i in 0..n {
        mid.px[i] += half * (u1[i][0] + freestream[0]);
        mid.py[i] += half * (u1[i][1] + freestream[1]);
        mid.pz[i] += half * (u1[i][2] + freestream[2]);
    }

    // Stage 2: velocity at the midpoint, used for the full step.
    let u2 = eval(&mid);
    for i in 0..n {
        field.px[i] += dt * (u2[i][0] + freestream[0]);
        field.py[i] += dt * (u2[i][1] + freestream[1]);
        field.pz[i] += dt * (u2[i][2] + freestream[2]);
    }
}

/// Advance the free wake by one step of size `dt` using midpoint (RK2)
/// integration. Each particle convects at its own induced velocity plus the
/// uniform `freestream`. Strengths are held constant (convection only --
/// stretching and diffusion are not modelled yet).
pub fn advect_rk2(field: &mut ParticleField, freestream: [f32; 3], dt: f32) {
    advect_rk2_with(field, freestream, dt, induced_velocities);
}

/// Sequential (non-parallelized) variant of `advect_rk2`. Uses the sequential
/// induced-velocity evaluator for single-threaded execution.
pub fn advect_rk2_seq(field: &mut ParticleField, freestream: [f32; 3], dt: f32) {
    advect_rk2_with(field, freestream, dt, induced_velocities_seq);
}

/// Debug variant of `advect_rk2` that routes through the scalar NaN-asserting
/// induced-velocity path. Panics at the first non-finite contribution with
/// source/target indices and values. Very slow -- O(N^2) non-vectorized.
pub fn advect_rk2_nan_check(field: &mut ParticleField, freestream: [f32; 3], dt: f32) {
    advect_rk2_with(field, freestream, dt, |f| {
        induced_at_points_nan_check(f, &f.px, &f.py, &f.pz)
    });
}

// ---------------------------------------------------------------------------
// Barnes-Hut O(N log N) evaluator
// ---------------------------------------------------------------------------
//
// The direct evaluator above is O(N^2): every target sums every source. The
// Barnes-Hut tree groups distant sources into a single "super-particle" -- the
// lumped vector strength Sum(alpha) placed at the strength-weighted centre of
// the cell -- whenever the cell is far enough that the monopole approximation
// is accurate. Accuracy is set by the opening angle `theta`: a cell of width
// `s` seen at distance `d` is accepted as far when s / d < theta. Smaller
// theta -> more direct sums -> more accurate and slower.
//
// The far field of the regularized kernel IS the singular Biot-Savart law
// (K -> 1/r^3 as r >> sigma), so a far cell is evaluated with the *same*
// kernel as a real particle: a super-particle carrying Sum(alpha) and a
// representative sigma, so the same eight-wide kernel serves both. The tree is
// flattened into a compact pre-order node array with escape pointers (no child
// arrays, no recursion) and the particles are reordered into leaf-contiguous,
// 8-padded blocks at build time. Per target the stackless walk runs the kernel
// directly over each near leaf's packed block and batches the far-cell
// monopoles into one small padded pass -- there is no per-target near-field
// gather. Traversal is scalar and cheap (O(N log N)); the arithmetic stays
// vectorized.
//
// The expansion centre is weighted by |alpha| so the strongest particles are
// represented at their own location (this cancels the |vorticity| dipole and
// keeps the monopole error small). This is a monopole (p = 0) Barnes-Hut; a
// dipole term could be added later for more accuracy at a given theta.

// Leaves hold up to this many particles. Larger => shallower tree, fewer nodes
// for the walk to visit, and bigger SIMD leaf blocks (more work handed to the
// cheap vector kernel); the build pays a little more to save the hot walk.
const BH_LEAF_MAX: usize = 32;
pub(crate) const BH_MIN_HALF: f32 = 1e-5;

/// Transient node used only while building the tree recursively. It carries
/// the |alpha|-weighting accumulators and a per-leaf index list; after the
/// flatten pass these are discarded and only [`FlatNode`] + the packed source
/// arrays survive.
pub(crate) struct TmpNode {
    pub(crate) half: f32,
    /// Strength-weighted expansion centre (falls back to the cube centre
    /// when the cell carries no net strength).
    pub(crate) center: [f32; 3],
    /// Sum(alpha) over the cell -- the monopole used for far interactions.
    pub(crate) sum_a: [f32; 3],
    /// Representative core size (|alpha|-weighted mean sigma).
    pub(crate) sigma_rep: f32,
    /// |alpha|-weighting accumulators, kept so parents can combine children.
    pub(crate) wsum: f32,
    wpos: [f32; 3],
    wsig: f32,
    pub(crate) children: [i32; 8],
    /// Particle indices, populated only for leaves.
    pub(crate) particles: Vec<u32>,
    pub(crate) is_leaf: bool,
}

impl TmpNode {
    fn placeholder(cube_center: [f32; 3], half: f32) -> Self {
        Self {
            half,
            center: cube_center,
            sum_a: [0.0; 3],
            sigma_rep: 1.0,
            wsum: 0.0,
            wpos: [0.0; 3],
            wsig: 0.0,
            children: [-1; 8],
            particles: Vec::new(),
            is_leaf: false,
        }
    }
}

/// Compact octree node laid out in DFS pre-order.
///
/// The walk reads `center` + `four_half2` on every visit (the opening-angle
/// test) and `sum_a` + `sigma_rep` only when a cell is accepted as far. There
/// is no child array: in pre-order the first child is always the next node, so
/// a near internal node descends to `i + 1`, while a far or leaf node jumps to
/// `escape` (the index just past this whole subtree). Leaves reference a
/// contiguous, 8-padded block of the packed source arrays via
/// `[leaf_start, leaf_start + 8 * leaf_chunks)`; `leaf_chunks == 0` marks an
/// internal node.
#[derive(Clone)]
struct FlatNode {
    center: [f32; 3],
    /// (2 * half)^2, precomputed so the far test is `four_half2 < theta2 * d2`.
    four_half2: f32,
    sum_a: [f32; 3],
    sigma_rep: f32,
    escape: u32,
    leaf_start: u32,
    leaf_chunks: u32,
}

/// Flattened octree plus the reordered, leaf-contiguous source arrays.
///
/// Built once per evaluation. The particles are permuted into leaf order and
/// each leaf is padded to a multiple of eight (zero strength, sigma = 1), so a
/// leaf's sources are a whole number of aligned `f32x8` chunks the kernel can
/// consume directly -- no per-target gather for the near field.
struct FlatTree {
    nodes: Vec<FlatNode>,
    sx: AVec<f32>,
    sy: AVec<f32>,
    sz: AVec<f32>,
    sax: AVec<f32>,
    say: AVec<f32>,
    saz: AVec<f32>,
    sg: AVec<f32>,
}

impl FlatTree {
    fn build(field: &ParticleField) -> Self {
        let mut tree = Self {
            nodes: Vec::new(),
            sx: AVec::new(SIMD_ALIGN),
            sy: AVec::new(SIMD_ALIGN),
            sz: AVec::new(SIMD_ALIGN),
            sax: AVec::new(SIMD_ALIGN),
            say: AVec::new(SIMD_ALIGN),
            saz: AVec::new(SIMD_ALIGN),
            sg: AVec::new(SIMD_ALIGN),
        };
        let n = field.len();
        if n == 0 {
            return tree;
        }
        // Per-axis bounds via eight-wide min/max reductions (each axis is its
        // own contiguous array, so the whole coordinate stream vectorizes).
        let (lo_x, hi_x) = simd_min_max(&field.px);
        let (lo_y, hi_y) = simd_min_max(&field.py);
        let (lo_z, hi_z) = simd_min_max(&field.pz);
        let lo = [lo_x, lo_y, lo_z];
        let hi = [hi_x, hi_y, hi_z];
        let center = [
            0.5 * (lo[0] + hi[0]),
            0.5 * (lo[1] + hi[1]),
            0.5 * (lo[2] + hi[2]),
        ];
        let mut half = 0.0f32;
        for k in 0..3 {
            half = half.max(0.5 * (hi[k] - lo[k]));
        }
        half = (half * 1.0001).max(BH_MIN_HALF);

        // Recursive temporary tree (aggregates + leaf index lists), then a
        // single pre-order pass that flattens it and emits the packed leaves.
        let idx: Vec<u32> = (0..n as u32).collect();
        let mut tmp = Vec::new();
        let root = build_node(&mut tmp, field, idx, center, half);
        tree.nodes.reserve(tmp.len());
        tree.flatten(&tmp, root, field);
        tree
    }

    /// Emit `tmp[ti]` and its subtree in DFS pre-order, filling escape links
    /// and (for leaves) the packed, 8-padded source block. Escape of each node
    /// is the node index just past its subtree.
    fn flatten(&mut self, tmp: &[TmpNode], ti: usize, field: &ParticleField) {
        let my = self.nodes.len();
        let src = &tmp[ti];
        let w = 2.0 * src.half;
        self.nodes.push(FlatNode {
            center: src.center,
            four_half2: w * w,
            sum_a: src.sum_a,
            sigma_rep: src.sigma_rep,
            escape: 0,
            leaf_start: 0,
            leaf_chunks: 0,
        });

        if src.is_leaf {
            let start = self.sx.len() as u32;
            for &p in &src.particles {
                let pi = p as usize;
                self.sx.push(field.px[pi]);
                self.sy.push(field.py[pi]);
                self.sz.push(field.pz[pi]);
                self.sax.push(field.ax[pi]);
                self.say.push(field.ay[pi]);
                self.saz.push(field.az[pi]);
                self.sg.push(field.sigma[pi]);
            }
            // Pad the block to a multiple of eight (zero strength, sigma = 1).
            while self.sx.len() % 8 != 0 {
                self.sx.push(0.0);
                self.sy.push(0.0);
                self.sz.push(0.0);
                self.sax.push(0.0);
                self.say.push(0.0);
                self.saz.push(0.0);
                self.sg.push(1.0);
            }
            let chunks = (self.sx.len() as u32 - start) / 8;
            self.nodes[my].leaf_start = start;
            self.nodes[my].leaf_chunks = chunks;
        } else {
            for o in 0..8 {
                let c = src.children[o];
                if c >= 0 {
                    self.flatten(tmp, c as usize, field);
                }
            }
        }
        self.nodes[my].escape = self.nodes.len() as u32;
    }

    /// Evaluate the induced velocity at one target with a stackless pre-order
    /// walk: near leaves run the SIMD kernel directly over their packed block,
    /// far cells are batched into `far` and kernelled in one pass at the end.
    #[inline]
    fn evaluate(&self, tx: f32, ty: f32, tz: f32, theta2: f32, far: &mut Scratch) -> [f32; 3] {
        let xj = f32x8::splat(tx);
        let yj = f32x8::splat(ty);
        let zj = f32x8::splat(tz);
        let mut ux = f32x8::splat(0.0);
        let mut uy = f32x8::splat(0.0);
        let mut uz = f32x8::splat(0.0);

        // Aligned eight-wide views of the packed leaf sources (the whole
        // buffer is a multiple of eight; each leaf owns a chunk sub-range).
        let src = Chunks::from_avecs(
            &self.sx, &self.sy, &self.sz, &self.sax, &self.say, &self.saz, &self.sg,
        );

        far.clear();
        let end = self.nodes.len();
        let mut i = 0usize;
        while i < end {
            let node = &self.nodes[i];
            let dx = tx - node.center[0];
            let dy = ty - node.center[1];
            let dz = tz - node.center[2];
            let d2 = dx * dx + dy * dy + dz * dz;
            if node.four_half2 < theta2 * d2 {
                // Far cell: batch its monopole for the tail kernel pass.
                far.push(
                    node.center[0],
                    node.center[1],
                    node.center[2],
                    node.sum_a[0],
                    node.sum_a[1],
                    node.sum_a[2],
                    node.sigma_rep,
                );
                i = node.escape as usize;
            } else if node.leaf_chunks != 0 {
                // Near leaf: kernel straight over its packed block.
                let c0 = (node.leaf_start / 8) as usize;
                let c1 = c0 + node.leaf_chunks as usize;
                accumulate_chunks(xj, yj, zj, &src, c0, c1, &mut ux, &mut uy, &mut uz);
                i = node.escape as usize;
            } else {
                // Near internal node: descend to the first child (next node).
                i += 1;
            }
        }

        // Far monopoles gathered above, in one padded kernel pass.
        far.pad_to_8();
        let far_chunks = far.len / 8;
        if far_chunks != 0 {
            let fsrc = Chunks::from_avecs(
                &far.sx, &far.sy, &far.sz, &far.sax, &far.say, &far.saz, &far.sg,
            );
            accumulate_chunks(xj, yj, zj, &fsrc, 0, far_chunks, &mut ux, &mut uy, &mut uz);
        }

        reduce_velocity(ux, uy, uz)
    }
}

/// Build a subtree over `idx` inside the cube (`cc`, `half`); returns the node
/// index. Recursive; leaves hold up to `BH_LEAF_MAX` particles.
pub(crate) fn build_node(
    nodes: &mut Vec<TmpNode>,
    field: &ParticleField,
    idx: Vec<u32>,
    cc: [f32; 3],
    half: f32,
) -> usize {
    let my = nodes.len();
    nodes.push(TmpNode::placeholder(cc, half));

    let finalize = |wsum: f32, wpos: [f32; 3], wsig: f32| -> ([f32; 3], f32) {
        if wsum > 1e-20 {
            (
                [wpos[0] / wsum, wpos[1] / wsum, wpos[2] / wsum],
                wsig / wsum,
            )
        } else {
            (cc, 1.0)
        }
    };

    if idx.len() <= BH_LEAF_MAX || half < BH_MIN_HALF {
        let mut wsum = 0.0f32;
        let mut wpos = [0.0f32; 3];
        let mut wsig = 0.0f32;
        let mut sum_a = [0.0f32; 3];
        for &p in &idx {
            let pi = p as usize;
            let a = [field.ax[pi], field.ay[pi], field.az[pi]];
            let w = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
            wsum += w;
            wsig += w * field.sigma[pi];
            for k in 0..3 {
                sum_a[k] += a[k];
            }
            wpos[0] += w * field.px[pi];
            wpos[1] += w * field.py[pi];
            wpos[2] += w * field.pz[pi];
        }
        let (center, sigma_rep) = finalize(wsum, wpos, wsig);
        let node = &mut nodes[my];
        node.is_leaf = true;
        node.particles = idx;
        node.wsum = wsum;
        node.wpos = wpos;
        node.wsig = wsig;
        node.sum_a = sum_a;
        node.center = center;
        node.sigma_rep = sigma_rep;
        return my;
    }

    // Partition into the eight octants of this cube.
    let mut buckets: [Vec<u32>; 8] = std::array::from_fn(|_| Vec::new());
    for &p in &idx {
        let pi = p as usize;
        let mut o = 0usize;
        if field.px[pi] >= cc[0] {
            o |= 1;
        }
        if field.py[pi] >= cc[1] {
            o |= 2;
        }
        if field.pz[pi] >= cc[2] {
            o |= 4;
        }
        buckets[o].push(p);
    }

    let ch = half * 0.5;
    let mut wsum = 0.0f32;
    let mut wpos = [0.0f32; 3];
    let mut wsig = 0.0f32;
    let mut sum_a = [0.0f32; 3];
    let mut children = [-1i32; 8];
    for o in 0..8 {
        if buckets[o].is_empty() {
            continue;
        }
        let ccc = [
            cc[0] + if o & 1 != 0 { ch } else { -ch },
            cc[1] + if o & 2 != 0 { ch } else { -ch },
            cc[2] + if o & 4 != 0 { ch } else { -ch },
        ];
        let bucket = std::mem::take(&mut buckets[o]);
        let ci = build_node(nodes, field, bucket, ccc, ch);
        children[o] = ci as i32;
        let child = &nodes[ci];
        wsum += child.wsum;
        wsig += child.wsig;
        for k in 0..3 {
            wpos[k] += child.wpos[k];
            sum_a[k] += child.sum_a[k];
        }
    }
    let (center, sigma_rep) = finalize(wsum, wpos, wsig);
    let node = &mut nodes[my];
    node.children = children;
    node.wsum = wsum;
    node.wpos = wpos;
    node.wsig = wsig;
    node.sum_a = sum_a;
    node.center = center;
    node.sigma_rep = sigma_rep;
    my
}

/// Reused per-target gather buffer (near particles + far cells as sources).
///
/// The backing arrays are sized once to the worst-case gather (every source
/// direct) and never reallocated or freed between targets: `clear` just resets
/// the write cursor `len`, and `push` writes through that cursor with a plain
/// aligned store instead of a bounds-growing `Vec::push`. `sg` is initialised
/// to 1.0 so any never-written padding slot has a finite sigma.
struct Scratch {
    sx: AVec<f32>,
    sy: AVec<f32>,
    sz: AVec<f32>,
    sax: AVec<f32>,
    say: AVec<f32>,
    saz: AVec<f32>,
    sg: AVec<f32>,
    /// Number of source slots written so far (the SIMD kernel reads `len / 8`
    /// eight-wide chunks). Reset by `clear`, never shrinks the allocation.
    len: usize,
}

impl Scratch {
    /// Allocate buffers big enough for `n` gathered sources (rounded up to a
    /// multiple of eight for the padded SIMD tail). Done once per evaluation;
    /// the storage is reused for every target.
    fn with_capacity(n: usize) -> Self {
        let cap = n.div_ceil(8) * 8;
        let mk = |fill: f32| -> AVec<f32> {
            let mut v = AVec::<f32>::with_capacity(SIMD_ALIGN, cap);
            v.resize(cap, fill);
            v
        };
        Self {
            sx: mk(0.0),
            sy: mk(0.0),
            sz: mk(0.0),
            sax: mk(0.0),
            say: mk(0.0),
            saz: mk(0.0),
            sg: mk(1.0),
            len: 0,
        }
    }

    #[inline]
    fn clear(&mut self) {
        self.len = 0;
    }

    #[inline]
    #[allow(clippy::too_many_arguments)]
    fn push(&mut self, x: f32, y: f32, z: f32, ax: f32, ay: f32, az: f32, sig: f32) {
        let k = self.len;
        self.sx[k] = x;
        self.sy[k] = y;
        self.sz[k] = z;
        self.sax[k] = ax;
        self.say[k] = ay;
        self.saz[k] = az;
        self.sg[k] = sig;
        self.len += 1;
    }

    /// Pad up to a multiple of eight with zero-strength (sigma = 1) sources so
    /// the SIMD loop has no remainder.
    #[inline]
    fn pad_to_8(&mut self) {
        while self.len % 8 != 0 {
            let k = self.len;
            self.sx[k] = 0.0;
            self.sy[k] = 0.0;
            self.sz[k] = 0.0;
            self.sax[k] = 0.0;
            self.say[k] = 0.0;
            self.saz[k] = 0.0;
            self.sg[k] = 1.0;
            self.len += 1;
        }
    }
}

/// Per-axis min/max of a coordinate array via eight-wide reductions, with a
/// scalar tail for the final < 8 elements. Used for the octree bounding box.
pub(crate) fn simd_min_max(v: &[f32]) -> (f32, f32) {
    if v.is_empty() {
        return (0.0, 0.0);
    }
    let n8 = v.len() / 8 * 8;
    let (head, tail) = v.split_at(n8);
    let mut vmin = f32x8::splat(f32::INFINITY);
    let mut vmax = f32x8::splat(f32::NEG_INFINITY);
    if !head.is_empty() {
        let head8: &[f32x8] = cast_slice(head);
        for &x in head8 {
            vmin = vmin.min(x);
            vmax = vmax.max(x);
        }
    }
    let mut lo = vmin.to_array().into_iter().fold(f32::INFINITY, f32::min);
    let mut hi = vmax
        .to_array()
        .into_iter()
        .fold(f32::NEG_INFINITY, f32::max);
    for &r in tail {
        if r < lo {
            lo = r;
        }
        if r > hi {
            hi = r;
        }
    }
    (lo, hi)
}

/// Aligned eight-wide views of a set of SoA source arrays (positions, vector
/// strengths, core sizes). Built by reinterpreting multiple-of-eight, aligned
/// `AVec<f32>` buffers as `f32x8` slices -- shared by the direct, near-leaf,
/// and far-cell kernel paths.
struct Chunks<'a> {
    x: &'a [f32x8],
    y: &'a [f32x8],
    z: &'a [f32x8],
    ax: &'a [f32x8],
    ay: &'a [f32x8],
    az: &'a [f32x8],
    sg: &'a [f32x8],
}

impl<'a> Chunks<'a> {
    #[inline]
    #[allow(clippy::too_many_arguments)]
    fn from_avecs(
        x: &'a AVec<f32>,
        y: &'a AVec<f32>,
        z: &'a AVec<f32>,
        ax: &'a AVec<f32>,
        ay: &'a AVec<f32>,
        az: &'a AVec<f32>,
        sg: &'a AVec<f32>,
    ) -> Self {
        Self {
            x: cast_slice(x),
            y: cast_slice(y),
            z: cast_slice(z),
            ax: cast_slice(ax),
            ay: cast_slice(ay),
            az: cast_slice(az),
            sg: cast_slice(sg),
        }
    }
}

/// Apply the 1/4pi prefactor and horizontally reduce the eight-wide velocity
/// accumulators to a single `[vx, vy, vz]`.
#[inline]
fn reduce_velocity(ux: f32x8, uy: f32x8, uz: f32x8) -> [f32; 3] {
    let inv = f32x8::splat(INV_4PI);
    let ux = (ux * inv).to_array();
    let uy = (uy * inv).to_array();
    let uz = (uz * inv).to_array();
    [
        ux.iter().sum::<f32>(),
        uy.iter().sum::<f32>(),
        uz.iter().sum::<f32>(),
    ]
}

/// Accumulate the regularized Biot-Savart velocity from source chunks
/// `[c0, c1)` of `src` into the running `f32x8` accumulators. This is the one
/// SIMD kernel used by every path: the direct O(N^2) sum, the near-leaf packed
/// tree blocks (read in place), and the batched far-cell monopoles.
#[inline]
#[allow(clippy::too_many_arguments)]
fn accumulate_chunks(
    xj: f32x8,
    yj: f32x8,
    zj: f32x8,
    src: &Chunks,
    c0: usize,
    c1: usize,
    ux: &mut f32x8,
    uy: &mut f32x8,
    uz: &mut f32x8,
) {
    // Bounds-check elision: assert the range fits every source array once, up
    // front, so LLVM can prove `c < c1 <= len` for each `[c]` in the hot loop
    // and drop the per-element checks. (The arrays are all the same padded
    // length, so these never trip in practice.)
    assert!(c0 <= c1);
    assert!(c1 <= src.x.len());
    assert!(c1 <= src.y.len());
    assert!(c1 <= src.z.len());
    assert!(c1 <= src.ax.len());
    assert!(c1 <= src.ay.len());
    assert!(c1 <= src.az.len());
    assert!(c1 <= src.sg.len());

    for c in c0..c1 {
        let sx = src.x[c];
        let sy = src.y[c];
        let sz = src.z[c];
        let sax = src.ax[c];
        let say = src.ay[c];
        let saz = src.az[c];
        let sigma = src.sg[c];

        // d = x_target - x_source
        let dx = xj - sx;
        let dy = yj - sy;
        let dz = zj - sz;

        let r2 = dx * dx + dy * dy + dz * dz;
        let sigma2 = sigma * sigma;
        let sigma3 = sigma2 * sigma;
        let rho2 = r2 / sigma2;
        let base = rho2 + f32x8::splat(1.0);
        // (rho2 + 1)^(5/2) = base^2 * sqrt(base)
        let denom = base * base * base.sqrt();
        // K = (rho2 + 5/2) / (sigma^3 * (rho2 + 1)^(5/2))
        let k = (rho2 + f32x8::splat(2.5)) / (sigma3 * denom);

        // cross = alpha x d
        let cx = say * dz - saz * dy;
        let cy = saz * dx - sax * dz;
        let cz = sax * dy - say * dx;

        *ux += k * cx;
        *uy += k * cy;
        *uz += k * cz;
    }
}

/// Barnes-Hut version of [`induced_at_points`]. Same kernel and sign
/// convention; approximate to the opening angle `theta` (0 recovers the
/// direct sum, larger is faster and coarser; `theta ~ 0.5` is the usual
/// accuracy/speed sweet spot).
pub fn induced_at_points_bh(
    field: &ParticleField,
    tx: &[f32],
    ty: &[f32],
    tz: &[f32],
    theta: f32,
) -> Vec<[f32; 3]> {
    let n = field.len();
    let m = tx.len();
    if n == 0 || m == 0 {
        return vec![[0.0f32; 3]; m];
    }
    let tree = FlatTree::build(field);
    let theta2 = theta * theta;
    let cap = tree.nodes.len();

    let mut out = vec![[0.0f32; 3]; m];

    // Parallel: each Rayon task handles one target; per-thread Scratch
    // eliminates cross-thread mutable state.  FlatTree is read-only (Sync).
    // Using par_chunks_mut(CHUNK) gives each thread CHUNK targets per
    // task -- enough work to amortize Rayon dispatch overhead, while
    // staying cache-friendly (each thread's Scratch stays warm between
    // consecutive targets in the same chunk).
    const CHUNK: usize = 64;
    out.par_chunks_mut(CHUNK)
        .zip(tx.par_chunks(CHUNK))
        .zip(ty.par_chunks(CHUNK))
        .zip(tz.par_chunks(CHUNK))
        .for_each(|(((oc, tc_x), tc_y), tc_z)| {
            let mut far = Scratch::with_capacity(cap);
            for k in 0..tc_x.len() {
                oc[k] = tree.evaluate(tc_x[k], tc_y[k], tc_z[k], theta2, &mut far);
            }
        });

    out
}

/// Barnes-Hut version of [`induced_velocities`] (targets are the particles).
pub fn induced_velocities_bh(field: &ParticleField, theta: f32) -> Vec<[f32; 3]> {
    induced_at_points_bh(field, &field.px, &field.py, &field.pz, theta)
}

/// Sequential (single-threaded) BH velocity evaluation, for benchmarking
/// the parallel speedup of the BH path.  Bypasses the Rayon dispatch in
/// `induced_at_points_bh` and runs the per-target tree walk on one thread.
pub fn induced_velocities_bh_seq(field: &ParticleField, theta: f32) -> Vec<[f32; 3]> {
    let n = field.len();
    if n == 0 {
        return vec![];
    }
    let tree = FlatTree::build(field);
    let theta2 = theta * theta;
    let cap = tree.nodes.len();
    let mut far = Scratch::with_capacity(cap);
    let mut out = vec![[0.0f32; 3]; n];
    for j in 0..n {
        out[j] = tree.evaluate(field.px[j], field.py[j], field.pz[j], theta2, &mut far);
    }
    out
}

/// Compute Barnes-Hut induced velocities at arbitrary points (sequential).
/// This is the single-threaded variant of `induced_at_points_bh`.
pub fn induced_at_points_bh_seq(
    field: &ParticleField,
    tx: &[f32],
    ty: &[f32],
    tz: &[f32],
    theta: f32,
) -> Vec<[f32; 3]> {
    let n = field.len();
    let m = tx.len();
    if n == 0 || m == 0 {
        return vec![[0.0f32; 3]; m];
    }
    let tree = FlatTree::build(field);
    let theta2 = theta * theta;
    let cap = tree.nodes.len();

    let mut out = vec![[0.0f32; 3]; m];
    let mut far = Scratch::with_capacity(cap);
    for j in 0..m {
        out[j] = tree.evaluate(tx[j], ty[j], tz[j], theta2, &mut far);
    }

    out
}

/// Barnes-Hut version of [`advect_rk2`]. Identical midpoint integration, but
/// both velocity evaluations use the tree (rebuilt each stage, since the
/// particles move).
pub fn advect_rk2_bh(field: &mut ParticleField, freestream: [f32; 3], dt: f32, theta: f32) {
    advect_rk2_with(field, freestream, dt, |f| induced_velocities_bh(f, theta));
}

/// Sequential (non-parallelized) variant of `advect_rk2_bh`. Uses the
/// sequential Barnes-Hut evaluator for single-threaded execution.
pub fn advect_rk2_bh_seq(field: &mut ParticleField, freestream: [f32; 3], dt: f32, theta: f32) {
    advect_rk2_with(field, freestream, dt, |f| {
        induced_velocities_bh_seq(f, theta)
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    // The merge/aging tests below exercise the sibling submodules directly.
    use crate::vpm::aging::{core_spread, strength_decay};
    use crate::vpm::merge::{merge_particles, MergeOpts};
    use std::f64::consts::PI;

    /// Single particle with strength along +z induces velocity in +y at a
    /// point on the +x axis (right-hand rule), and in the far field
    /// (rho >> 1, g -> 1) the magnitude matches the singular Biot-Savart
    /// value |u| = A / (4 pi d^2).
    #[test]
    fn far_field_matches_singular_biot_savart() {
        let mut f = ParticleField::new();
        let a = 1.0f32;
        f.push([0.0, 0.0, 0.0], [0.0, 0.0, a], 0.01); // thin core
        f.push([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.01); // massless probe on +x

        let u = induced_velocities(&f);
        let up = u[1];
        let expected = (a / (4.0 * std::f32::consts::PI * 1.0 * 1.0)) as f32;

        assert!(up[0].abs() < 1e-4, "u_x should vanish, got {}", up[0]);
        assert!(up[2].abs() < 1e-4, "u_z should vanish, got {}", up[2]);
        assert!(
            up[1] > 0.0,
            "u_y should be +y (CCW about +z), got {}",
            up[1]
        );
        assert!(
            (up[1] - expected).abs() / expected < 1e-2,
            "u_y = {}, expected ~ {}",
            up[1],
            expected
        );
    }

    /// A lone particle induces zero velocity on itself (the regularized self
    /// term is exactly zero).
    #[test]
    fn self_induced_velocity_is_zero() {
        let mut f = ParticleField::new();
        f.push([0.3, -0.2, 1.1], [0.5, -0.7, 0.9], 0.05);
        let u = induced_velocities(&f);
        assert!(
            u[0].iter().all(|c| c.abs() < 1e-6),
            "self velocity {:?}",
            u[0]
        );
    }

    /// The SIMD path agrees with the scalar f64 reference on a random cloud.
    #[test]
    fn simd_matches_reference() {
        // Deterministic pseudo-random cloud (simple LCG, no rand dependency).
        let mut state = 0x1234_5678u32;
        let mut rng = || {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            (state >> 8) as f32 / (1u32 << 24) as f32 - 0.5
        };
        let mut f = ParticleField::new();
        for _ in 0..137 {
            f.push(
                [rng() * 4.0, rng() * 4.0, rng() * 4.0],
                [rng(), rng(), rng()],
                0.1 + 0.05 * (rng() + 0.5).abs(),
            );
        }

        let simd = induced_velocities(&f);
        let refv = induced_velocities_ref(&f);
        for j in 0..f.len() {
            for k in 0..3 {
                let a = simd[j][k] as f64;
                let b = refv[j][k];
                assert!(
                    (a - b).abs() <= 1e-3 * (1.0 + b.abs()),
                    "particle {} comp {}: simd {} vs ref {}",
                    j,
                    k,
                    a,
                    b
                );
            }
        }
    }

    /// A thin vortex ring self-propagates along its axis at the Kelvin
    /// thin-core speed U = Gamma / (4 pi R) [ln(8 R / a) - 1/4]. Discretize a
    /// ring in the xy-plane with vorticity tangent (-sin, cos, 0); it should
    /// translate in +z (its particles move at ~U), with negligible transverse
    /// velocity by symmetry.
    #[test]
    fn vortex_ring_self_propagation() {
        let radius = 1.0f64;
        let circulation = 1.0f64; // Gamma
        let n = 200usize;
        let sigma = 0.1f64; // core ~ a
        let seg = 2.0 * PI * radius / n as f64; // segment length

        let mut f = ParticleField::new();
        for i in 0..n {
            let theta = 2.0 * PI * (i as f64) / (n as f64);
            let (s, c) = theta.sin_cos();
            let pos = [(radius * c) as f32, (radius * s) as f32, 0.0];
            // alpha = Gamma * seg * tangent, tangent = (-sin, cos, 0)
            let strength = [
                (circulation * seg * -s) as f32,
                (circulation * seg * c) as f32,
                0.0,
            ];
            f.push(pos, strength, sigma as f32);
        }

        let u = induced_velocities(&f);
        // Mean axial velocity of the ring particles = ring translation speed.
        let mut vz = 0.0f64;
        let mut vtrans = 0.0f64;
        for uu in &u {
            vz += uu[2] as f64;
            vtrans += ((uu[0] as f64).powi(2) + (uu[1] as f64).powi(2)).sqrt();
        }
        vz /= n as f64;
        vtrans /= n as f64;

        let u_kelvin = circulation / (4.0 * PI * radius) * ((8.0 * radius / sigma).ln() - 0.25);

        assert!(vz > 0.0, "ring should propagate in +z, got vz = {}", vz);
        assert!(
            (vz - u_kelvin).abs() / u_kelvin < 0.30,
            "ring speed {} vs Kelvin {} (>30% off)",
            vz,
            u_kelvin
        );
        assert!(
            vtrans < 0.1 * vz,
            "transverse velocity {} should be small vs axial {}",
            vtrans,
            vz
        );
    }

    /// One RK2 step moves a thin ring downstream in +z by ~ U * dt.
    #[test]
    fn advect_moves_ring_downstream() {
        let radius = 1.0f32;
        let n = 120usize;
        let sigma = 0.1f32;
        let seg = 2.0 * std::f32::consts::PI * radius / n as f32;
        let mut f = ParticleField::new();
        for i in 0..n {
            let theta = 2.0 * std::f32::consts::PI * (i as f32) / (n as f32);
            let (s, c) = theta.sin_cos();
            f.push(
                [radius * c, radius * s, 0.0],
                [seg * -s, seg * c, 0.0],
                sigma,
            );
        }
        let z0: f32 = f.pz.iter().sum::<f32>() / n as f32;
        advect_rk2(&mut f, [0.0, 0.0, 0.0], 0.1);
        let z1: f32 = f.pz.iter().sum::<f32>() / n as f32;
        assert!(z1 > z0, "ring mean z should increase: {} -> {}", z0, z1);
    }

    /// Deterministic pseudo-random particle cloud (LCG, no rand dependency).
    fn random_cloud(n: usize) -> ParticleField {
        let mut state = 0x9E37_79B9u32;
        let mut rng = || {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            (state >> 8) as f32 / (1u32 << 24) as f32 - 0.5
        };
        let mut f = ParticleField::new();
        for _ in 0..n {
            f.push(
                [rng() * 6.0, rng() * 6.0, rng() * 6.0],
                [rng(), rng(), rng()],
                0.15 + 0.1 * (rng() + 0.5).abs(),
            );
        }
        f
    }

    /// Merging conserves total vector circulation exactly and reduces the
    /// particle count. With a large kappa the far coherent cells collapse, so
    /// the count must drop; `sum(alpha)` is invariant to f32 rounding.
    #[test]
    fn merge_conserves_circulation_and_shrinks() {
        let f = random_cloud(2000);
        let sum_before = [
            f.ax.iter().map(|&v| v as f64).sum::<f64>(),
            f.ay.iter().map(|&v| v as f64).sum::<f64>(),
            f.az.iter().map(|&v| v as f64).sum::<f64>(),
        ];
        // Aggressive, coherence gate off, merge everywhere.
        let opts = MergeOpts {
            kappa: 100.0,
            chi_min: 0.0,
            region_dist: 0.0,
            min_particles: 0,
        };
        let m = merge_particles(&f, &opts);
        let sum_after = [
            m.ax.iter().map(|&v| v as f64).sum::<f64>(),
            m.ay.iter().map(|&v| v as f64).sum::<f64>(),
            m.az.iter().map(|&v| v as f64).sum::<f64>(),
        ];
        assert!(
            m.len() < f.len(),
            "merge should reduce particle count: {} -> {}",
            f.len(),
            m.len()
        );
        for k in 0..3 {
            let tol = 1e-3 * (1.0 + sum_before[k].abs());
            assert!(
                (sum_before[k] - sum_after[k]).abs() <= tol,
                "circulation comp {} not conserved: {} -> {}",
                k,
                sum_before[k],
                sum_after[k]
            );
        }
    }

    /// The coherence gate protects sheared cells: a cloud whose strengths
    /// cancel (low |sum_a|/wsum) is left untouched when chi_min is high.
    #[test]
    fn merge_coherence_gate_preserves_sheared_wake() {
        let f = random_cloud(1500);
        // High coherence requirement + huge kappa: only near-aligned cells
        // would merge, but random strengths cancel, so few/none collapse.
        let opts = MergeOpts {
            kappa: 100.0,
            chi_min: 0.99,
            region_dist: 0.0,
            min_particles: 0,
        };
        let m = merge_particles(&f, &opts);
        // Random-direction strengths rarely satisfy |sum_a| >= 0.99*wsum, so
        // the field is essentially preserved.
        assert!(
            m.len() >= f.len() * 9 / 10,
            "coherence gate should preserve most particles: {} -> {}",
            f.len(),
            m.len()
        );
    }

    /// Core spreading grows every core (sigma^2 += 2*nu*dt) and leaves the
    /// strengths (circulation) untouched.
    #[test]
    fn core_spread_grows_sigma_conserves_circulation() {
        let mut f = random_cloud(500);
        let sig_before: Vec<f32> = f.sigma.iter().copied().collect();
        let sum_before: f64 = f.ax.iter().map(|&v| v as f64).sum();
        let nu = 5.0;
        let dt = 0.01;
        core_spread(&mut f, nu, dt);
        let expect = (2.0 * nu * dt) as f32;
        for i in 0..f.len() {
            let grown = (sig_before[i] * sig_before[i] + expect).sqrt();
            assert!((f.sigma[i] - grown).abs() <= 1e-4 * (1.0 + grown));
        }
        let sum_after: f64 = f.ax.iter().map(|&v| v as f64).sum();
        assert!((sum_before - sum_after).abs() <= 1e-6 * (1.0 + sum_before.abs()));
    }

    /// Strength fade scales every strength by the factor and leaves positions
    /// and cores untouched.
    #[test]
    fn strength_decay_scales_alpha_only() {
        let mut f = random_cloud(500);
        let ax0: Vec<f32> = f.ax.iter().copied().collect();
        let px0: Vec<f32> = f.px.iter().copied().collect();
        let sig0: Vec<f32> = f.sigma.iter().copied().collect();
        let factor = 0.9f32;
        strength_decay(&mut f, factor);
        for i in 0..f.len() {
            assert!((f.ax[i] - ax0[i] * factor).abs() <= 1e-6 * (1.0 + ax0[i].abs()));
            assert_eq!(f.px[i], px0[i]);
            assert_eq!(f.sigma[i], sig0[i]);
        }
    }

    /// The Barnes-Hut evaluator agrees with the direct O(N^2) sum to within
    /// the opening-angle tolerance on a random cloud (monopole, theta = 0.5).
    #[test]
    fn barnes_hut_matches_direct() {
        let f = random_cloud(800);
        let direct = induced_velocities(&f);
        let bh = induced_velocities_bh(&f, 0.5);

        // Peak velocity magnitude sets the error scale.
        let peak = direct
            .iter()
            .map(|u| (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt())
            .fold(0.0f32, f32::max);
        let max_err = direct
            .iter()
            .zip(&bh)
            .map(|(d, b)| {
                let e = [b[0] - d[0], b[1] - d[1], b[2] - d[2]];
                (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt()
            })
            .fold(0.0f32, f32::max);
        assert!(
            max_err / peak < 0.05,
            "Barnes-Hut theta=0.5 max error {:.4} vs peak {:.4} ({:.1}%)",
            max_err,
            peak,
            100.0 * max_err / peak
        );
    }

    /// As theta -> 0 the tree forces direct sums everywhere, so Barnes-Hut
    /// converges to the exact result (bit-level differences only from a
    /// different summation order).
    #[test]
    fn barnes_hut_small_theta_is_near_exact() {
        let f = random_cloud(300);
        let direct = induced_velocities(&f);
        let bh = induced_velocities_bh(&f, 0.05);
        let peak = direct
            .iter()
            .map(|u| (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt())
            .fold(0.0f32, f32::max);
        let max_err = direct
            .iter()
            .zip(&bh)
            .map(|(d, b)| {
                let e = [b[0] - d[0], b[1] - d[1], b[2] - d[2]];
                (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt()
            })
            .fold(0.0f32, f32::max);
        assert!(
            max_err / peak < 1e-3,
            "small-theta Barnes-Hut should match direct: err {:.2e} vs peak {:.2e}",
            max_err,
            peak
        );
    }

    /// The Barnes-Hut advect step moves a thin ring downstream just like the
    /// direct integrator (same physics, approximate evaluator).
    #[test]
    fn barnes_hut_advect_matches_direct_ring() {
        let radius = 1.0f32;
        let n = 240usize;
        let sigma = 0.1f32;
        let seg = 2.0 * std::f32::consts::PI * radius / n as f32;
        let build = || {
            let mut f = ParticleField::new();
            for i in 0..n {
                let theta = 2.0 * std::f32::consts::PI * (i as f32) / (n as f32);
                let (s, c) = theta.sin_cos();
                f.push(
                    [radius * c, radius * s, 0.0],
                    [seg * -s, seg * c, 0.0],
                    sigma,
                );
            }
            f
        };
        let mut fd = build();
        let mut fb = build();
        advect_rk2(&mut fd, [0.0, 0.0, 0.0], 0.1);
        advect_rk2_bh(&mut fb, [0.0, 0.0, 0.0], 0.1, 0.4);
        let zd: f32 = fd.pz.iter().sum::<f32>() / n as f32;
        let zb: f32 = fb.pz.iter().sum::<f32>() / n as f32;
        assert!(zb > 0.0, "BH ring should move in +z, got {}", zb);
        assert!(
            (zd - zb).abs() / zd < 0.05,
            "BH ring displacement {} vs direct {} (>5% off)",
            zb,
            zd
        );
    }
}
