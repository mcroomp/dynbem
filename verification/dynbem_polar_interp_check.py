"""Validate dynbem.TabulatedPolar against numpy.interp on the same data.

This is the *spec-compliance* check: dynbem's TabulatedPolar
interpolates with a hand-rolled binary search + linear interp; the
docstring claims it matches np.interp's behaviour (clamp-at-endpoints,
linear-between-knots, no periodic wrap-around). This script confirms
that on a real polar -- the S809 table used for NREL Phase VI --
across:

  1. Each tabulated alpha (should be exact, machine-precision match).
  2. 4001 evenly-spaced alphas spanning [alpha_min - 5deg, alpha_max + 5deg]
     (a mix of within-table interp and outside-table clamping).
  3. 10,000 random alphas drawn uniformly from the same range.

A separate script (dynbem_polar_vs_aerodyn_nrel_phase_vi.py) does the
cross-implementation check: does dynbem's polar lookup at a given
alpha agree with what AeroDyn computed at the same alpha during a
real simulation?

Reports max and mean absolute error vs np.interp for both CL and CD.
Also smoke-checks dynbem.LinearPolar against its closed-form formula.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dynbem import LinearPolar, TabulatedPolar                          # noqa: E402


S809_YAML = ROOT / "verification" / "ccblade_docker" / "inputs" / "nrel_phase_vi.yaml"


def load_s809_polar() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (alpha_rad, cl, cd) for the S809 polar."""
    with S809_YAML.open() as fh:
        cfg = yaml.safe_load(fh)
    pts = cfg["airfoil_polar"]["points"]
    a = np.array([math.radians(p["alpha_deg"]) for p in pts], dtype=np.float64)
    cl = np.array([p["cl"] for p in pts], dtype=np.float64)
    cd = np.array([p["cd"] for p in pts], dtype=np.float64)
    return a, cl, cd


def check_tabulated_vs_np_interp() -> int:
    a_tab, cl_tab, cd_tab = load_s809_polar()
    p = TabulatedPolar(alpha_rad=a_tab, cl=cl_tab, cd=cd_tab)
    print(f"S809 polar: {len(a_tab)} rows, alpha = "
          f"[{math.degrees(a_tab[0]):+.2f}, {math.degrees(a_tab[-1]):+.2f}] deg")

    # 1) Exact-knot accuracy
    cl_db = np.array([p.cl_cd(a)[0] for a in a_tab])
    cd_db = np.array([p.cl_cd(a)[1] for a in a_tab])
    knot_cl_err = np.abs(cl_db - cl_tab).max()
    knot_cd_err = np.abs(cd_db - cd_tab).max()
    print(f"\nAt {len(a_tab)} tabulated knots:")
    print(f"  max |Cl_dynbem - Cl_table| = {knot_cl_err:.2e}")
    print(f"  max |Cd_dynbem - Cd_table| = {knot_cd_err:.2e}")

    # 2) Dense alpha grid extending 5 deg beyond either end of the table
    alpha_min = a_tab[0] - math.radians(5)
    alpha_max = a_tab[-1] + math.radians(5)
    grid = np.linspace(alpha_min, alpha_max, 4001)

    cl_db_grid = np.array([p.cl_cd(a)[0] for a in grid])
    cd_db_grid = np.array([p.cl_cd(a)[1] for a in grid])
    cl_np_grid = np.interp(grid, a_tab, cl_tab)
    cd_np_grid = np.interp(grid, a_tab, cd_tab)

    cl_err = np.abs(cl_db_grid - cl_np_grid)
    cd_err = np.abs(cd_db_grid - cd_np_grid)
    print(f"\nOn {len(grid)} alphas in [{math.degrees(alpha_min):+.1f}, "
          f"{math.degrees(alpha_max):+.1f}] deg vs np.interp:")
    print(f"  Cl:  max abs err = {cl_err.max():.2e}   mean = {cl_err.mean():.2e}")
    print(f"  Cd:  max abs err = {cd_err.max():.2e}   mean = {cd_err.mean():.2e}")

    # 3) Random samples (catches order-of-evaluation bugs in the binary search)
    rng = np.random.default_rng(seed=42)
    rand_alpha = rng.uniform(alpha_min, alpha_max, size=10_000)
    cl_db_rand = np.array([p.cl_cd(a)[0] for a in rand_alpha])
    cd_db_rand = np.array([p.cl_cd(a)[1] for a in rand_alpha])
    cl_np_rand = np.interp(rand_alpha, a_tab, cl_tab)
    cd_np_rand = np.interp(rand_alpha, a_tab, cd_tab)
    rcl_err = np.abs(cl_db_rand - cl_np_rand).max()
    rcd_err = np.abs(cd_db_rand - cd_np_rand).max()
    print(f"\nOn 10,000 random alphas in the same range vs np.interp:")
    print(f"  Cl:  max abs err = {rcl_err:.2e}")
    print(f"  Cd:  max abs err = {rcd_err:.2e}")

    # 4) Verify clamping behaviour outside the table
    cl_below, cd_below = p.cl_cd(alpha_min - math.radians(20))
    cl_above, cd_above = p.cl_cd(alpha_max + math.radians(20))
    print(f"\nClamping at alpha << first knot ({math.degrees(alpha_min) - 20:.0f} deg):")
    print(f"  (Cl, Cd) dynbem = ({cl_below:+.4f}, {cd_below:.4f})   "
          f"table[0] = ({cl_tab[0]:+.4f}, {cd_tab[0]:.4f})")
    print(f"Clamping at alpha >> last knot ({math.degrees(alpha_max) + 20:.0f} deg):")
    print(f"  (Cl, Cd) dynbem = ({cl_above:+.4f}, {cd_above:.4f})   "
          f"table[-1] = ({cl_tab[-1]:+.4f}, {cd_tab[-1]:.4f})")

    # Tolerance: dynbem and np.interp both use double; remaining drift is
    # from the (alpha-a_lo)/(a_hi-a_lo) division vs np.interp's slope-based
    # form. We expect <1e-12.
    TOL = 1e-12
    ok = (cl_err.max() < TOL and cd_err.max() < TOL and
          rcl_err < TOL and rcd_err < TOL and
          knot_cl_err < TOL and knot_cd_err < TOL)
    if ok:
        print(f"\nPASS: all errors < {TOL:.0e} (machine precision).")
        return 0
    print(f"\nFAIL: at least one error >= {TOL:.0e}.")
    return 1


def check_linear_polar_formula() -> int:
    """Check LinearPolar matches the documented closed-form expression."""
    p = LinearPolar(CL0=0.05, CL_alpha_per_rad=5.7,
                    CD0=0.012, alpha_stall_rad=math.radians(14.0))
    print()
    print("LinearPolar(CL0=0.05, CL_alpha=5.7/rad, CD0=0.012, stall=14 deg)")

    # Within-stall: linear in alpha
    alpha = np.linspace(-math.radians(13), math.radians(13), 401)
    expected_cl = p.CL0 + p.CL_alpha_per_rad * alpha
    expected_cd = np.full_like(alpha, p.CD0)
    got_cl = np.array([p.cl_cd(a)[0] for a in alpha])
    got_cd = np.array([p.cl_cd(a)[1] for a in alpha])
    cl_err = np.abs(got_cl - expected_cl).max()
    cd_err = np.abs(got_cd - expected_cd).max()
    print(f"  Within |alpha| < 13 deg:   max |dCl| = {cl_err:.2e}, "
          f"max |dCd| = {cd_err:.2e}")

    # Beyond stall: Cl held constant (with sign), Cd grows linearly
    a_pos = math.radians(25)
    cl_stall_mag = p.CL0 + p.CL_alpha_per_rad * p.alpha_stall_rad
    cl_pos, cd_pos = p.cl_cd(a_pos)
    cl_neg, cd_neg = p.cl_cd(-a_pos)
    expected_cl_pos = cl_stall_mag
    expected_cl_neg = -cl_stall_mag
    expected_cd_excess = a_pos - p.alpha_stall_rad
    ok = (abs(cl_pos - expected_cl_pos) < 1e-12 and
          abs(cl_neg - expected_cl_neg) < 1e-12 and
          abs(cd_pos - (p.CD0 + expected_cd_excess)) < 1e-12 and
          abs(cd_neg - (p.CD0 + expected_cd_excess)) < 1e-12)
    print(f"  At alpha = +/-{math.degrees(a_pos):.0f} deg (stalled):")
    print(f"    (Cl, Cd)+ = ({cl_pos:+.4f}, {cd_pos:.4f})   "
          f"expected = ({expected_cl_pos:+.4f}, {p.CD0 + expected_cd_excess:.4f})")
    print(f"    (Cl, Cd)- = ({cl_neg:+.4f}, {cd_neg:.4f})   "
          f"expected = ({expected_cl_neg:+.4f}, {p.CD0 + expected_cd_excess:.4f})")
    return 0 if ok and cl_err < 1e-12 and cd_err < 1e-12 else 1


def main() -> int:
    rc1 = check_tabulated_vs_np_interp()
    rc2 = check_linear_polar_formula()
    return rc1 | rc2


if __name__ == "__main__":
    sys.exit(main())
