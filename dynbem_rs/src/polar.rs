// Airfoil polars: maps angle of attack -> (CL, CD).
//
// Two implementations:
//   LinearPolar     - analytical CL = CL0 + CL_alpha*alpha, clipped at stall.
//   TabulatedPolar  - alpha/cl/cd arrays with binary-search interp (matches
//                     np.interp). Batched path loops elementwise; the table
//                     is small (~50 entries) so log2 search is cheap and
//                     LLVM auto-vectorizes the arithmetic between hits.
//
// Both expose:
//   cl_cd(alpha) -> (cl, cd)                              scalar
//   cl_cd_into(&[alpha], &mut[cl], &mut[cd])              batched, no alloc

pub trait Polar: Send + Sync {
    fn cl_cd(&self, alpha: f64) -> (f64, f64);

    /// Batched scalar fallback. Override for vectorized implementations.
    fn cl_cd_into(&self, alpha: &[f64], cl: &mut [f64], cd: &mut [f64]) {
        debug_assert_eq!(alpha.len(), cl.len());
        debug_assert_eq!(alpha.len(), cd.len());
        for i in 0..alpha.len() {
            let (l, d) = self.cl_cd(alpha[i]);
            cl[i] = l;
            cd[i] = d;
        }
    }
}

// ---------------------------------------------------------------------------
// LinearPolar
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
#[allow(non_snake_case)]
pub struct LinearPolar {
    pub CL0: f64,
    pub CL_alpha_per_rad: f64,
    pub CD0: f64,
    pub alpha_stall_rad: f64,
}

impl LinearPolar {
    #[allow(non_snake_case)]
    pub fn new(CL0: f64, CL_alpha_per_rad: f64, CD0: f64, alpha_stall_rad: f64) -> Self {
        Self {
            CL0,
            CL_alpha_per_rad,
            CD0,
            alpha_stall_rad,
        }
    }
}

impl Polar for LinearPolar {
    #[inline]
    fn cl_cd(&self, alpha: f64) -> (f64, f64) {
        if alpha.abs() < self.alpha_stall_rad {
            (self.CL0 + self.CL_alpha_per_rad * alpha, self.CD0)
        } else {
            let cl_mag = self.CL0 + self.CL_alpha_per_rad * self.alpha_stall_rad;
            let cl = cl_mag.copysign(alpha);
            let cd = self.CD0 + (alpha.abs() - self.alpha_stall_rad);
            (cl, cd)
        }
    }

    fn cl_cd_into(&self, alpha: &[f64], cl: &mut [f64], cd: &mut [f64]) {
        debug_assert_eq!(alpha.len(), cl.len());
        debug_assert_eq!(alpha.len(), cd.len());
        let cl0 = self.CL0;
        let cla = self.CL_alpha_per_rad;
        let cd0 = self.CD0;
        let astall = self.alpha_stall_rad;
        let cl_stall_mag = cl0 + cla * astall;
        for i in 0..alpha.len() {
            let a = alpha[i];
            let abs_a = a.abs();
            let stalled = abs_a >= astall;
            cl[i] = if stalled {
                cl_stall_mag.copysign(a)
            } else {
                cl0 + cla * a
            };
            cd[i] = if stalled { cd0 + (abs_a - astall) } else { cd0 };
        }
    }
}

// ---------------------------------------------------------------------------
// TabulatedPolar
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct TabulatedPolar {
    pub alpha: Vec<f64>,
    pub cl: Vec<f64>,
    pub cd: Vec<f64>,
}

impl TabulatedPolar {
    pub fn new(alpha: Vec<f64>, cl: Vec<f64>, cd: Vec<f64>) -> Result<Self, String> {
        if alpha.len() != cl.len() || alpha.len() != cd.len() {
            return Err("alpha, cl, cd must have the same length".into());
        }
        if alpha.is_empty() {
            return Err("polar table is empty".into());
        }
        Ok(Self { alpha, cl, cd })
    }

    /// Binary-search interpolation, matches numpy.interp (clamp at endpoints).
    #[inline]
    fn interp_at(&self, alpha: f64) -> (f64, f64) {
        let a = &self.alpha[..];
        let n = a.len();
        if n == 0 {
            return (0.0, 0.0);
        }
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
        let cl_lo = self.cl[lo];
        let cd_lo = self.cd[lo];
        let cl = cl_lo + t * (self.cl[hi] - cl_lo);
        let cd = cd_lo + t * (self.cd[hi] - cd_lo);
        (cl, cd)
    }
}

impl Polar for TabulatedPolar {
    #[inline]
    fn cl_cd(&self, alpha: f64) -> (f64, f64) {
        self.interp_at(alpha)
    }
}

// ---------------------------------------------------------------------------
// PolarKind: enum used internally by models so they hold the polar by value
// without a Box<dyn>. Keeps the inner loop call site monomorphic.
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
pub enum PolarKind {
    Linear(LinearPolar),
    Tabulated(TabulatedPolar),
}

impl Polar for PolarKind {
    #[inline]
    fn cl_cd(&self, alpha: f64) -> (f64, f64) {
        match self {
            PolarKind::Linear(p) => p.cl_cd(alpha),
            PolarKind::Tabulated(p) => p.cl_cd(alpha),
        }
    }
    fn cl_cd_into(&self, alpha: &[f64], cl: &mut [f64], cd: &mut [f64]) {
        match self {
            PolarKind::Linear(p) => p.cl_cd_into(alpha, cl, cd),
            PolarKind::Tabulated(p) => p.cl_cd_into(alpha, cl, cd),
        }
    }
}
