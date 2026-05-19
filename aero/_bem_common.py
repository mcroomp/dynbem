"""Shared infrastructure for the Level-2 dynamic-inflow BEM models.

Holds the pieces that ``PittPetersModel(JIT)`` and ``OyeBEMModel`` would
otherwise duplicate or cross-import:

  * ``vrs_lambda1``  — Leishman VRS empirical polynomial
  * ``_interp_polar`` — JIT polar (cl, cd) lookup
  * ``build_polar_arrays`` — one-time sampling of any AirfoilPolar onto
    contiguous numba arrays
  * ``radial_grid`` — one-time radial geometry caching

Hot-path kinematics and result-assembly are deliberately *not* here —
they're a few cheap numpy ops and live inline in each model's
``compute_forces`` to avoid the Python function-call overhead.

Lets the two models stay as peers (Øye no longer imports from
Pitt-Peters) and the ψ-loop kernels stay model-specific (each has its
own ``lam_local(r, ψ)`` formula and can't share a numba-compiled body
cleanly — closures don't compose under ``@njit``).
"""
from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
from numba import njit

from .polar import AirfoilPolar, TabulatedPolar

if TYPE_CHECKING:
    from .rotor_definition import BladeGeometry


# ---------------------------------------------------------------------------
# VRS empirical polynomial (Leishman 2000, §12.7)
# ---------------------------------------------------------------------------
# λ_1/V_h = 1 + C[0]·λ₂ + C[1]·λ₂² + C[2]·λ₂³ + C[3]·λ₂⁴
# where λ₂ = V_descent / V_h > 0.  Valid for 0 ≤ λ₂ ≤ 2.
# Fit to Castles-Gray (NACA TN-2474) and Coleman (1945) measured data.
_VRS_C = (1.125, -1.372, 1.718, -0.655)


def vrs_lambda1(lambda2: float) -> float:
    """Normalised induced velocity λ₁ = v_i/V_h from Leishman VRS polynomial.

    lambda2  V_descent / V_h, must be in [0, 2].
    Returns  v_i / V_h  (= 1.0 at λ₂=0 hover; ≈1.0 at λ₂=2 WBS boundary).
    """
    k = lambda2
    return 1.0 + _VRS_C[0]*k + _VRS_C[1]*k**2 + _VRS_C[2]*k**3 + _VRS_C[3]*k**4


# ---------------------------------------------------------------------------
# JIT polar interpolator (cl, cd at angle α)
# ---------------------------------------------------------------------------

@njit(cache=True, fastmath=True)
def _interp_polar(alpha, alpha_tab, cl_tab, cd_tab):
    """Linear-interp (cl, cd) at angle alpha (rad).  Binary-search lookup;
    clamps to endpoints outside the tabulated range.  Matches np.interp.
    """
    n = alpha_tab.shape[0]
    if alpha <= alpha_tab[0]:
        return cl_tab[0], cd_tab[0]
    if alpha >= alpha_tab[n - 1]:
        return cl_tab[n - 1], cd_tab[n - 1]
    lo = 0
    hi = n - 1
    while hi - lo > 1:
        mid = (lo + hi) >> 1
        if alpha_tab[mid] <= alpha:
            lo = mid
        else:
            hi = mid
    a_lo = alpha_tab[lo]
    a_hi = alpha_tab[hi]
    t = (alpha - a_lo) / (a_hi - a_lo)
    return (
        cl_tab[lo] + t * (cl_tab[hi] - cl_tab[lo]),
        cd_tab[lo] + t * (cd_tab[hi] - cd_tab[lo]),
    )


# ---------------------------------------------------------------------------
# Polar tabulation for the JIT kernels
# ---------------------------------------------------------------------------

def build_polar_arrays(polar: AirfoilPolar
                       ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sample any AirfoilPolar onto contiguous numba arrays for _interp_polar.

    TabulatedPolar passes through its existing arrays; analytical polars
    (LinearPolar etc.) are sampled to 4001 points over [−π/2, π/2].
    """
    if isinstance(polar, TabulatedPolar):
        return (
            np.ascontiguousarray(polar._alpha),
            np.ascontiguousarray(polar._cl),
            np.ascontiguousarray(polar._cd),
        )
    n = 4001
    a  = np.linspace(-math.pi / 2, math.pi / 2, n)
    cl = np.empty(n)
    cd = np.empty(n)
    for i in range(n):
        cl_i, cd_i = polar.cl_cd(float(a[i]))
        cl[i] = cl_i
        cd[i] = cd_i
    return (
        np.ascontiguousarray(a),
        np.ascontiguousarray(cl),
        np.ascontiguousarray(cd),
    )


# ---------------------------------------------------------------------------
# Radial grid caching
# ---------------------------------------------------------------------------

def radial_grid(blade: "BladeGeometry"
                ) -> tuple[float, np.ndarray, np.ndarray, float, float]:
    """Cache the fixed radial geometry for a JIT BEM kernel.

    Returns
    -------
    dr          width of each radial element (m)
    r_mid       (n,) midpoint radius per element (m), contiguous
    x_mid       (n,) midpoint r/R, contiguous
    x_hub       root-cutout/R (dimensionless)
    twist_rad   uniform twist (rad).  Per-section twist not yet supported.
    """
    R  = blade.radius_m
    n  = blade.n_elements
    r0 = blade.root_cutout_m
    dr = (R - r0) / n
    r_mid = np.ascontiguousarray(
        np.linspace(r0 + 0.5 * dr, R - 0.5 * dr, n)
    )
    x_mid     = np.ascontiguousarray(r_mid / R)
    x_hub     = float(r0 / R) if R > 0.0 else 0.0
    twist_rad = math.radians(blade.twist_deg)
    return dr, r_mid, x_mid, x_hub, twist_rad
