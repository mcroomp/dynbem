// Shared infrastructure for the BEM models.
// See ../CLAUDE.md "Shared BEM infrastructure".

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
