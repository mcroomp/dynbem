"""Whole-dataset validation of dynbem.bem spanwise sectional CL against
Caradonna & Tung (1981) NASA TM-81232.

Companion to tests/test_caradonna_spanwise.py, which spot-checks three
operating points with hard assertions. This script sweeps all 32 CT
tables and prints a per-station error breakdown plus aggregate stats --
the numbers quoted in that test's docstring come from running this.
Re-run after any change to the BEM or to the Research/csv/CaradonnaTung
data; if the aggregates shift, update the test docstring.

For each (theta_c, Omega) operating point in C-T Tables 1-32, the
paper measured airfoil section CL at five radial stations
(r/R = 0.50, 0.68, 0.80, 0.89, 0.96) by integrating chordwise Cp.
That data lives in Research/csv/CaradonnaTung/page_NN_table_*__cl.csv .

We call solve_bem_element directly at each station (much faster than
running the full BEM and post-processing), read the converged
lambda_r/a_prime, and recover the local airfoil CL = polar.cl(alpha).

Caveats
-------
- C-T tested at tip Mach up to 0.89; our BEM polar is incompressible.
  Compare at the lowest-Mach operating points to minimize bias
  (Tables 27/28 at 1250 rpm M_tip ~ 0.43; Tables 26/27 at 650 rpm
  M_tip ~ 0.22).
- C-T Appendix A scans for tables 1-32 are known to have OCR damage
  on some pages (see CaradonnaTung/CLAUDE.md).  The extracted CL rows
  may carry spurious values for tables flagged in that index.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


def _sample_evenly(items: list, n: int | None) -> list:
    """Return n evenly-spaced items from ``items``. ``n=None`` returns all."""
    if n is None or n <= 0 or n >= len(items):
        return items
    idx = np.linspace(0, len(items) - 1, n).round().astype(int)
    seen: set[int] = set()
    out = []
    for i in idx:
        i = int(i)
        if i not in seen:
            seen.add(i)
            out.append(items[i])
    return out

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dynbem.bem import solve_bem_element
from dynbem.rotor_definition import (
    AirfoilProperties, AutorotationProperties, BladeGeometry, RotorDefinition)
from dynbem.factory import build_polar

CSV_DIR = ROOT / "Research" / "csv" / "CaradonnaTung"

# Caradonna-Tung rotor (NASA TM-81232 page_01.png)
ROTOR = RotorDefinition(
    blade=BladeGeometry(
        n_blades=2, radius_m=1.143, root_cutout_m=0.1,
        chord_m=0.1905, twist_deg=0.0, n_elements=30),
    airfoil=AirfoilProperties(
        Re_design=1_000_000, CL0=0.0,
        CL_alpha_per_rad=2 * math.pi, CD0=0.008,
        alpha_stall_deg=15.0, tip_loss=True),
    autorotation=AutorotationProperties(I_ode_kgm2=1.0),
    name="Caradonna-Tung",
)
POLAR = build_polar(ROTOR.airfoil)

# Operating points: one entry per Cp table (table 1-32) in C-T.
# (table_num, collective_deg, omega_rpm, M_tip) - the M_tip column is for
# annotation only; we don't use it in the BEM (incompressible polar).
OPERATING_POINTS = [
    ( 1,  0.0, 1500, 0.520),
    ( 2,  2.0, 1250, 0.436),
    ( 3,  2.0, 1500, 0.520),
    ( 4,  2.0, 1750, 0.607),
    ( 5,  2.0, 2062, 0.723),
    ( 6,  2.0, 2365, 0.796),
    ( 7,  4.0, 1334, 0.415),
    ( 8,  4.0, 1600, 0.830),
    ( 9,  5.0,  450, 0.328),
    (10,  5.0, 1250, 0.433),
    (11,  5.0, 1500, 0.520),
    (12,  5.0, 1750, 0.607),
    (13,  5.0, 2062, 0.723),
    (14,  5.0, 2365, 0.794),
    (15,  7.0, 2500, 0.877),
    (16,  7.0, 2500, 0.877),
    (17,  8.0, 1250, 0.439),
    (18,  8.0, 1250, 0.326),
    (19,  8.0, 1750, 0.612),
    (20,  8.0, 2050, 0.717),
    (21,  8.0, 2250, 0.794),
    (22,  8.0, 2500, 0.813),
    (23,  8.0, 2500, 0.827),
    (24,  8.0, 2500, 0.865),
    (25,  8.0, 2574, 0.890),
    (26, 10.0,  650, 0.226),
    (27, 12.0,  650, 0.226),
    (28, 12.0, 1250, 0.433),
    (29, 12.0, 1750, 0.330),
    (30, 12.0, 1750, 0.432),
    (31, 12.0, 2074, 0.723),
    (32, 12.0, 2280, 0.784),
]


@dataclass
class Comparison:
    table_num: int
    theta_deg: float
    rpm: float
    m_tip: float
    r_over_R: float
    cl_meas: float
    cl_bem: float

    @property
    def err(self) -> float:
        if abs(self.cl_meas) <= 0.01:
            return float("inf")
        return abs(self.cl_bem - self.cl_meas) / abs(self.cl_meas)


@dataclass
class Survey:
    comparisons: list[Comparison] = field(default_factory=list)
    points_run: int = 0
    points_total: int = 0
    sample: int | None = None

    def by_table(self, table_num: int) -> list[Comparison]:
        return [c for c in self.comparisons if c.table_num == table_num]

    def errors(self) -> np.ndarray:
        return np.array([c.err for c in self.comparisons
                         if math.isfinite(c.err)])

    def errors_low_mach(self) -> np.ndarray:
        return np.array([c.err for c in self.comparisons
                         if math.isfinite(c.err) and c.m_tip < 0.5])


def section_CL_bem(coll_deg: float, omega_rpm: float, r_over_R: float) -> float:
    """Solve the single-element BEM at the given station; return section CL."""
    R = ROTOR.blade.radius_m
    omega = omega_rpm * math.pi / 30.0
    Omega_R = omega * R
    r = r_over_R * R
    dr = 0.005 * R   # narrow annulus; CL is independent of dr after normalization
    elem = solve_bem_element(
        r=r, dr=dr,
        chord=ROTOR.blade.chord_m, twist_rad=0.0,
        collective_rad=math.radians(coll_deg),
        omega=omega, v_climb=0.0, rho=1.225,
        n_blades=ROTOR.blade.n_blades, radius_m=R,
        polar=POLAR, use_tip_loss=ROTOR.airfoil.tip_loss,
        root_cutout_m=ROTOR.blade.root_cutout_m,
    )
    v_a = elem.lambda_r * Omega_R
    v_t = omega * r * (1.0 + elem.a_prime)
    phi = math.atan2(v_a, v_t)
    alpha = math.radians(coll_deg) - phi
    cl, _ = POLAR.cl_cd(alpha)
    return cl


def load_cl_csv(table_num: int) -> dict[float, float] | None:
    """Return {r/R: CL_measured} from the table's __cl.csv, or None if missing."""
    for path in CSV_DIR.glob(f"page_*_table_{table_num}__cl.csv"):
        with path.open(encoding="ascii") as f:
            reader = csv.reader(f)
            header = next(reader)
            row = next(reader)
            out: dict[float, float] = {}
            for h, v in zip(header, row):
                m = re.match(r"r/R=(\d+\.\d+)", h)
                if m:
                    try:
                        out[float(m.group(1))] = float(v)
                    except ValueError:
                        pass
            return out
    return None


def run_survey(sample: int | None = None) -> Survey:
    """Run the spanwise-CL comparison across (sampled) operating points.

    `sample=None` -> all 32 tables. `sample=N` -> N evenly-spaced tables.
    Returns a `Survey` with per-station comparisons plus aggregate
    helpers. No print output -- callers (`main`, the unit test) format
    as they see fit.
    """
    points = _sample_evenly(OPERATING_POINTS, sample)
    survey = Survey(points_total=len(OPERATING_POINTS), sample=sample)
    for tbl_num, theta, rpm, m_tip in points:
        cl_meas = load_cl_csv(tbl_num)
        if cl_meas is None:
            continue
        survey.points_run += 1
        for r_over_R, cl_m in cl_meas.items():
            cl_b = section_CL_bem(theta, rpm, r_over_R)
            survey.comparisons.append(Comparison(
                table_num=tbl_num, theta_deg=theta, rpm=rpm, m_tip=m_tip,
                r_over_R=r_over_R, cl_meas=cl_m, cl_bem=cl_b))
    return survey


def _print_survey(survey: Survey) -> None:
    if survey.sample is not None and survey.points_run < survey.points_total:
        print(f"** SAMPLED MODE: {survey.points_run} of {survey.points_total} "
              f"operating points -- aggregate stats below are NOT comparable "
              f"to the full-sweep numbers quoted in test_caradonna_spanwise.py")
        print()

    print(f"Rotor: {ROTOR.name}, R={ROTOR.blade.radius_m}, "
          f"chord={ROTOR.blade.chord_m}, sigma={ROTOR.blade.n_blades*ROTOR.blade.chord_m/(math.pi*ROTOR.blade.radius_m):.4f}")
    print()
    print(f"{'tbl':>3} {'theta':>5} {'rpm':>5} {'M_tip':>5}   "
          f"r/R  | {'CL_meas':>8} {'CL_bem':>8} {'err':>6}")
    for c in survey.comparisons:
        print(f" {c.table_num:3d}  {c.theta_deg:5.1f}  {int(c.rpm):5d}  "
              f"{c.m_tip:5.3f}   {c.r_over_R:4.2f} | "
              f"{c.cl_meas:8.4f} {c.cl_bem:8.4f} {c.err:5.1%}")

    arr = survey.errors()
    print()
    print(f"Rows with measured CL: {survey.points_run}/{survey.points_total}")
    print(f"All-point comparisons: {len(arr)}")
    print(f"  mean error       = {arr.mean():.1%}")
    print(f"  median error     = {np.median(arr):.1%}")
    print(f"  RMSE             = {np.sqrt((arr**2).mean()):.1%}")
    print(f"  max error        = {arr.max():.1%}")

    lo = survey.errors_low_mach()
    if lo.size > 0:
        print(f"\nLow-Mach subset (M_tip < 0.5, {lo.size} comparisons):")
        print(f"  mean error       = {lo.mean():.1%}")
        print(f"  median error     = {np.median(lo):.1%}")
        print(f"  RMSE             = {np.sqrt((lo**2).mean()):.1%}")
        print(f"  max error        = {lo.max():.1%}")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else "")
    p.add_argument("--sample", type=int, default=None, metavar="N",
                   help="evenly-sample N of 32 operating points (smoke test); "
                        "omit to run the full sweep")
    args = p.parse_args()
    _print_survey(run_survey(sample=args.sample))
    return 0


if __name__ == "__main__":
    sys.exit(main())
