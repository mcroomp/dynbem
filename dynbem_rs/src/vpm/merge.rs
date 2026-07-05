//! Population control: tree-collapse particle merging.
//!
//! The Barnes-Hut octree (see [`super::barnes_hut`]) already computes, per
//! cell, exactly what a merge needs: `sum_a` (the exact vector-strength sum),
//! `center` (the |alpha|-weighted expansion centre -- the same one that cancels
//! the dipole and minimises the monopole error), `sigma_rep` (a representative
//! core), and `wsum` (the unsigned strength sum). Collapsing a cell into a
//! single particle `(center, sum_a, sigma_rep)` is therefore identical to
//! "accepting the cell as a far monopole" -- the same approximation the
//! evaluator makes, made permanent.
//!
//! A cell is collapsed when all three hold:
//!   * size vs core:   `2*half <= kappa * sigma_rep` (particles within ~a core
//!     radius, so the lumped near field barely changes);
//!   * coherence:      `|sum_a| >= chi_min * wsum` (strengths nearly aligned; a
//!     cell straddling a thin tip vortex has cancelling vorticity, low
//!     `|sum_a|/wsum`, and is NOT merged);
//!   * region gate:    `|center| >= region_dist` (only far/old wake; the
//!     near-blade wake stays at full resolution).
//!
//! Total vector circulation is conserved exactly; the linear-impulse error per
//! merge is bounded by the cell's dipole moment, i.e. the same order as the
//! Barnes-Hut opening-angle error -- keep `kappa` consistent with `bh_theta`.

use super::common::{build_node, simd_min_max, ParticleField, TmpNode, BH_MIN_HALF};

/// Options for [`merge_particles`].
#[derive(Clone, Copy, Debug)]
pub struct MergeOpts {
    /// Collapse a cell only when `2*half <= kappa * sigma_rep`. Larger merges
    /// more aggressively (coarser). ~0.5-1.0 keeps error near the BH floor.
    pub kappa: f32,
    /// Collapse a cell only when `|sum_a| >= chi_min * wsum` (coherence). Near
    /// 1.0 protects thin/sheared vortices; 0 disables the coherence gate.
    pub chi_min: f32,
    /// Only collapse cells whose centre is at least this far from the hub
    /// origin (m). 0 merges everywhere; a positive value preserves the
    /// near-blade wake.
    pub region_dist: f32,
    /// Never merge when the field has this many particles or fewer (leave
    /// small wakes untouched).
    pub min_particles: usize,
}

impl Default for MergeOpts {
    fn default() -> Self {
        Self {
            kappa: 0.7,
            chi_min: 0.9,
            region_dist: 0.0,
            min_particles: 0,
        }
    }
}

/// Whether a tree node may be collapsed into a single equivalent particle.
#[inline]
fn merge_ok(node: &TmpNode, opts: &MergeOpts) -> bool {
    // Size vs core: cell must be small relative to its representative core.
    if 2.0 * node.half > opts.kappa * node.sigma_rep {
        return false;
    }
    // Coherence: strengths must be nearly aligned (protects tip vortices).
    let mag = (node.sum_a[0] * node.sum_a[0]
        + node.sum_a[1] * node.sum_a[1]
        + node.sum_a[2] * node.sum_a[2])
        .sqrt();
    if mag < opts.chi_min * node.wsum {
        return false;
    }
    // Region gate: only far/old wake (distance from hub origin).
    if opts.region_dist > 0.0 {
        let d2 = node.center[0] * node.center[0]
            + node.center[1] * node.center[1]
            + node.center[2] * node.center[2];
        if d2 < opts.region_dist * opts.region_dist {
            return false;
        }
    }
    true
}

/// Emit `tmp[ti]` and its subtree into `out`, collapsing any cell that
/// satisfies [`merge_ok`] into a single equivalent particle. Cells that fail
/// the test are descended; leaves that fail emit their original particles
/// unchanged (exact).
fn collapse(tmp: &[TmpNode], ti: usize, field: &ParticleField, opts: &MergeOpts, out: &mut ParticleField) {
    let node = &tmp[ti];
    if merge_ok(node, opts) {
        // One equivalent particle carrying the exact strength sum. sigma_rep
        // is always > 0 (build_node's finalize falls back to 1.0).
        out.push(node.center, node.sum_a, node.sigma_rep.max(f32::MIN_POSITIVE));
        return;
    }
    if node.is_leaf {
        for &p in &node.particles {
            let pi = p as usize;
            out.push(
                [field.px[pi], field.py[pi], field.pz[pi]],
                [field.ax[pi], field.ay[pi], field.az[pi]],
                field.sigma[pi],
            );
        }
        return;
    }
    for o in 0..8 {
        let c = node.children[o];
        if c >= 0 {
            collapse(tmp, c as usize, field, opts, out);
        }
    }
}

/// Merge the wake for population control by collapsing small, coherent,
/// far-field octree cells into single equivalent particles. Reuses the same
/// octree build (`build_node`) as the Barnes-Hut evaluator, so a merge pass is
/// an `O(N log N)` build plus an `O(N)` collapse walk. Total vector
/// circulation is conserved exactly (see module notes above). Returns a new,
/// smaller `ParticleField`; the input is left untouched.
pub fn merge_particles(field: &ParticleField, opts: &MergeOpts) -> ParticleField {
    let n = field.len();
    if n <= opts.min_particles || n == 0 {
        return field.clone();
    }
    // Bounding cube (same construction as FlatTree::build).
    let (lo_x, hi_x) = simd_min_max(&field.px);
    let (lo_y, hi_y) = simd_min_max(&field.py);
    let (lo_z, hi_z) = simd_min_max(&field.pz);
    let center = [
        0.5 * (lo_x + hi_x),
        0.5 * (lo_y + hi_y),
        0.5 * (lo_z + hi_z),
    ];
    let mut half = 0.0f32;
    half = half.max(0.5 * (hi_x - lo_x));
    half = half.max(0.5 * (hi_y - lo_y));
    half = half.max(0.5 * (hi_z - lo_z));
    half = (half * 1.0001).max(BH_MIN_HALF);

    let idx: Vec<u32> = (0..n as u32).collect();
    let mut tmp = Vec::new();
    let root = build_node(&mut tmp, field, idx, center, half);
    let mut out = ParticleField::with_capacity(n);
    collapse(&tmp, root, field, opts, &mut out);
    out
}
