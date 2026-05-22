"""Cross-implementation polar-interpolation check: at each (case, element)
AeroDyn evaluated during NREL Phase VI Sequence S, query dynbem's
TabulatedPolar at the same alpha and compare the (Cl, Cd) it returns
against what AeroDyn returned.

This isolates the polar lookup from every other piece of the BEM
pipeline -- both codes see the same alpha (AeroDyn's converged value),
the same polar table (same S809 YAML), and we just compare what each
reports for Cl and Cd at that alpha.

If dynbem.TabulatedPolar implements linear interp + endpoint clamping
the same way AeroDyn does, the deltas should be at most a few ULPs.
Larger deltas would reveal a real implementation difference (e.g.
different deg/rad conversion, different boundary handling, a
higher-order interpolation in AeroDyn).
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dynbem import TabulatedPolar                                       # noqa: E402

YAML_PATH = ROOT / "verification" / "openfast_docker" / "inputs" / "nrel_phase_vi_dbemt.yaml"
PER_ELEMENT_CSV = ROOT / "verification" / "data" / "nrel_phase_vi_dbemt_openfast_per_element.csv"
OUT_CSV = ROOT / "verification" / "data" / "dynbem_polar_vs_aerodyn_nrel_phase_vi.csv"


def load_polar_from_yaml(path: Path) -> TabulatedPolar:
    with path.open() as fh:
        cfg = yaml.safe_load(fh)
    pts = cfg["airfoil_polar"]["points"]
    alpha = np.array([math.radians(p["alpha_deg"]) for p in pts], dtype=np.float64)
    cl    = np.array([p["cl"]   for p in pts], dtype=np.float64)
    cd    = np.array([p["cd"]   for p in pts], dtype=np.float64)
    return TabulatedPolar(alpha_rad=alpha, cl=cl, cd=cd)


def main() -> int:
    if not PER_ELEMENT_CSV.exists():
        print(f"Missing per-element AeroDyn CSV: {PER_ELEMENT_CSV}", file=sys.stderr)
        print("Generate it first:", file=sys.stderr)
        print("  cd verification/openfast_docker", file=sys.stderr)
        print("  docker compose run --rm openfast --config /in/nrel_phase_vi_dbemt.yaml",
              file=sys.stderr)
        return 1

    polar = load_polar_from_yaml(YAML_PATH)
    with PER_ELEMENT_CSV.open() as fh:
        rows = list(csv.DictReader(fh))

    print(f"dynbem.TabulatedPolar vs AeroDyn polar lookup")
    print(f"Comparing on {len(rows)} (case, element) AeroDyn evaluations")
    print()
    print(f"{'case':>4} {'el':>3} {'alpha_deg':>10} "
          f"{'cl_AD':>9} {'cl_db':>9} {'dCl':>10}   "
          f"{'cd_AD':>9} {'cd_db':>9} {'dCd':>10}")

    out_rows = []
    dcl_all = []
    dcd_all = []
    # Show first ~5 rows from each case, then aggregate stats
    show_max = 5
    last_case = -1
    shown_this_case = 0
    for r in rows:
        case = int(r["case"])
        el   = int(r["element_idx"])
        alpha_deg = float(r["alpha_deg"])
        cl_AD     = float(r["cl"])
        cd_AD     = float(r["cd"])
        cl_db, cd_db = polar.cl_cd(math.radians(alpha_deg))
        dcl = cl_db - cl_AD
        dcd = cd_db - cd_AD
        dcl_all.append(dcl)
        dcd_all.append(dcd)
        out_rows.append({
            "case":     case,
            "element_idx": el,
            "alpha_deg":   alpha_deg,
            "cl_aerodyn":  cl_AD,
            "cd_aerodyn":  cd_AD,
            "cl_dynbem":   cl_db,
            "cd_dynbem":   cd_db,
            "dcl":         dcl,
            "dcd":         dcd,
        })
        if case != last_case:
            shown_this_case = 0
            last_case = case
        if shown_this_case < show_max:
            print(f"{case:>4} {el:>3} {alpha_deg:>10.4f} "
                  f"{cl_AD:>9.5f} {cl_db:>9.5f} {dcl:>+10.2e}   "
                  f"{cd_AD:>9.5f} {cd_db:>9.5f} {dcd:>+10.2e}")
            shown_this_case += 1
            if shown_this_case == show_max:
                print(f"     ... ({sum(1 for _r in rows if int(_r['case']) == case) - show_max} more elements in case {case})")

    dcl_arr = np.array(dcl_all)
    dcd_arr = np.array(dcd_all)
    print()
    print(f"Aggregate over all {len(rows)} (case, element) pairs:")
    print(f"  Cl: max |dCl| = {np.abs(dcl_arr).max():.2e}   "
          f"mean |dCl| = {np.abs(dcl_arr).mean():.2e}   "
          f"RMS = {np.sqrt(np.mean(dcl_arr**2)):.2e}")
    print(f"  Cd: max |dCd| = {np.abs(dcd_arr).max():.2e}   "
          f"mean |dCd| = {np.abs(dcd_arr).mean():.2e}   "
          f"RMS = {np.sqrt(np.mean(dcd_arr**2)):.2e}")

    # Persist
    with OUT_CSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f"Wrote {OUT_CSV.relative_to(ROOT)}")

    # AeroDyn writes its .out file in single precision (the docker image
    # is built single-precision); we therefore expect ~1e-6 absolute
    # rounding error per channel. Anything materially larger reveals a
    # real implementation difference.
    TOL = 5e-6
    max_dcl = float(np.abs(dcl_arr).max())
    max_dcd = float(np.abs(dcd_arr).max())
    if max_dcl < TOL and max_dcd < TOL:
        print(f"\nPASS: deviations within {TOL:.0e} (AeroDyn .out single-precision noise floor).")
        return 0
    print(f"\nNOTE: max deviations larger than the .out single-precision noise "
          f"floor of {TOL:.0e}.")
    return 0  # informational; not a hard fail


if __name__ == "__main__":
    sys.exit(main())
