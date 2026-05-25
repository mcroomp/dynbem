"""Whole-dataset validation of dynbem.bem at TR 515 operating points.

For each row of Tables I-IV (Wheatley & Hood 1935 PCA-2 wind-tunnel),
the real rotor was in autorotation equilibrium -- so Q_aero == 0 by
construction.  Run the BEM at the prescribed (mu, alpha, Omega, pitch)
and report the residual Q_BEM (normalized).

Two views:
  (1) Q_BEM at prescribed kinematics, no trim -- the dirty baseline.
  (2) Q_BEM at prescribed kinematics with cyclic-trim inner loop --
      the flap-equivalent.

Companion to tests/test_wheatley_autorotation.py, which drives this
module in sampled mode for spot-test bounds.  The 469-row envelope
quoted in that test's docstring comes from running main() with no
sample -- re-run after any BEM change and update the test bounds if
the envelope shifts.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
import time
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

from dynbem.bem import BEMModel
from dynbem import RotorInputs, create_aero
from dynbem.rotor_definition import load as load_rotor
from dynbem.rotor_state import QuasiStaticRotorState
from dynbem.trim import solve_trim_cyclic

CSV_DIR = ROOT / "Research" / "csv" / "Wheatley_Hood_NACA515"
YAML = ROOT / "rotors" / "wheatley_pca2" / "rotor.yaml"
RHO = 1.225

# Each entry: (label, csv_filename, pitch_deg, confidence)
TABLES = [
    ("table_i",   "page_09_table_i.csv",   1.9, "MODERATE"),
    ("table_ii",  "page_09_table_ii.csv",  0.8, "MODERATE"),
    ("table_iii", "page_10_table_iii.csv", 1.9, "HIGH"),
    ("table_iv",  "page_10_table_iv.csv",  2.7, "HIGH"),
]


@dataclass
class Comparison:
    table_label: str
    pitch_deg: float
    mu: float
    alpha_deg: float
    N_rpm: float
    CT_meas: float
    CT_bem: float
    CQ_trim: float
    CQ_notrim: float | None
    tilt_lon: float
    tilt_lat: float

    @property
    def ct_ratio(self) -> float:
        if abs(self.CT_meas) < 1e-9:
            return float("inf")
        return self.CT_bem / self.CT_meas


@dataclass
class Survey:
    comparisons: list[Comparison] = field(default_factory=list)
    sample: int | None = None

    def by_table(self, label: str) -> list[Comparison]:
        return [c for c in self.comparisons if c.table_label == label]

    def cq_trim(self) -> np.ndarray:
        return np.array([c.CQ_trim for c in self.comparisons])

    def cq_notrim(self) -> np.ndarray:
        return np.array([c.CQ_notrim for c in self.comparisons
                         if c.CQ_notrim is not None])


def _load_csv(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open(encoding="ascii") as f:
        for r in csv.DictReader(f):
            try:
                rows.append({
                    "mu":    float(r["mu"]),
                    "alpha": float(r["alpha (deg)"]),
                    "CL":    float(r["CL"]),
                    "CD":    float(r["CD"]),
                    "N":     float(r["N (rpm)"]),
                })
            except (KeyError, ValueError):
                continue
    return rows


def load_table(label: str) -> list[dict[str, float]]:
    """Load one of the four CT tables; returns [] if its CSV is missing."""
    for lbl, name, _, _ in TABLES:
        if lbl == label:
            path = CSV_DIR / name
            return _load_csv(path) if path.exists() else []
    raise KeyError(label)


def load_model() -> BEMModel | None:
    """Build the PCA-2 BEM fixture (cached n_psi_elements=12 for speed).
    Returns None if the rotor.yaml is missing."""
    if not YAML.exists():
        return None
    return create_aero(load_rotor(str(YAML)), "bem", n_psi_elements=12)


def measured_CT(mu: float, alpha_deg: float, CL: float, CD: float) -> float:
    """Wheatley airplane-axes (CL, CD) -> rotor-axis CT.

        T_force  = L * cos(a) + D * sin(a)   where a = disk AoA
        q        = 1/2 * rho * V^2
        V        = OmegaR * mu / cos(a)
        CT_rotor = T / (rho * pi R^2 * (OmegaR)^2)
                 = (CL cos a + CD sin a) * mu^2 / (2 cos^2 a)
    """
    a = math.radians(alpha_deg)
    return (CL * math.cos(a) + CD * math.sin(a)) * mu**2 / (2.0 * math.cos(a)**2)


def _bem_at_point(model: BEMModel, mu: float, alpha_deg: float, N_rpm: float,
                  pitch_deg: float, *, trim: bool):
    """Inner BEM driver: returns (CT, CQ, tilt_lon, tilt_lat) at one row.

    ``trim=True`` runs solve_trim_cyclic first to emulate flap-equivalent
    response (the PCA-2 had freely flapping blades; the cyclic that
    zeros hub moments is its rigid-blade analogue)."""
    R = model.defn.blade.radius_m
    omega = N_rpm * math.pi / 30.0
    a = math.radians(alpha_deg)
    V = omega * R * mu / math.cos(a)
    R_hub = np.array([
        [math.cos(a), 0.0, -math.sin(a)],
        [0.0,         1.0,  0.0        ],
        [math.sin(a), 0.0,  math.cos(a)],
    ])
    v_hub_world = np.zeros(3)
    wind_world  = np.array([V, 0.0, 0.0])

    if trim:
        tr = solve_trim_cyclic(
            model, QuasiStaticRotorState(),
            collective_rad=math.radians(pitch_deg),
            R_hub=R_hub, v_hub_world=v_hub_world, wind_world=wind_world,
            omega_rad_s=omega,
            tilt_min=-math.radians(25.0), tilt_max=math.radians(25.0),
            tolerance_Nm=1.0, max_iterations=20, n_inflow_relax=0,
        )
        tilt_lon, tilt_lat, state = tr.tilt_lon, tr.tilt_lat, tr.final_state
    else:
        tilt_lon = tilt_lat = 0.0
        state = QuasiStaticRotorState()

    inp = RotorInputs(
        collective_rad=math.radians(pitch_deg),
        tilt_lon=tilt_lon, tilt_lat=tilt_lat,
        R_hub=R_hub, v_hub_world=v_hub_world, wind_world=wind_world, t=0.0,
        rho_kg_m3=RHO,
        omega_rad_s=omega,
    )
    res, _ = model.compute_forces(inp, state)
    F_hub = R_hub.T @ res.F_world
    T = -F_hub[2]
    A = math.pi * R**2
    CT = T / (RHO * A * (omega * R)**2)
    CQ = res.Q_spin / (RHO * A * (omega * R)**2 * R)
    return CT, CQ, tilt_lon, tilt_lat


def evaluate_point(model: BEMModel, row: dict[str, float], pitch_deg: float,
                   table_label: str, *, include_notrim: bool = True
                   ) -> Comparison:
    """Run trim + (optional) no-trim BEM at one CSV row; build a Comparison."""
    CT_meas = measured_CT(row["mu"], row["alpha"], row["CL"], row["CD"])
    CT_t, CQ_t, tlon, tlat = _bem_at_point(
        model, row["mu"], row["alpha"], row["N"], pitch_deg, trim=True)
    CQ_nt: float | None
    if include_notrim:
        _, CQ_nt, _, _ = _bem_at_point(
            model, row["mu"], row["alpha"], row["N"], pitch_deg, trim=False)
    else:
        CQ_nt = None
    return Comparison(
        table_label=table_label, pitch_deg=pitch_deg,
        mu=row["mu"], alpha_deg=row["alpha"], N_rpm=row["N"],
        CT_meas=CT_meas, CT_bem=CT_t,
        CQ_trim=CQ_t, CQ_notrim=CQ_nt,
        tilt_lon=tlon, tilt_lat=tlat,
    )


def run_survey(model: BEMModel | None = None, sample: int | None = None,
               *, include_notrim: bool = True) -> Survey:
    """Run the CQ-residual / CT survey across (sampled) TR 515 operating
    points.

    ``sample=None`` runs every available row (~469 total).
    ``sample=N`` evenly-samples N rows per table.
    Returns a `Survey` of Comparison records.  Missing CSVs and missing
    rotor.yaml result in an empty survey (callers can skip)."""
    if model is None:
        model = load_model()
    survey = Survey(sample=sample)
    if model is None:
        return survey
    for label, csv_name, pitch_deg, _ in TABLES:
        path = CSV_DIR / csv_name
        if not path.exists():
            continue
        rows = _sample_evenly(_load_csv(path), sample)
        for r in rows:
            try:
                cmp = evaluate_point(model, r, pitch_deg, label,
                                     include_notrim=include_notrim)
            except Exception:
                continue
            survey.comparisons.append(cmp)
    return survey


def _summarize(label: str, arr: np.ndarray) -> None:
    print(f"  {label:25s} n={len(arr):4d}  "
          f"mean={arr.mean():+.5f}  median={np.median(arr):+.5f}  "
          f"|mean|={abs(arr.mean()):.5f}  "
          f"RMSE={np.sqrt((arr**2).mean()):.5f}  "
          f"max|x|={np.max(np.abs(arr)):.5f}")


def _print_survey(survey: Survey, *, model: BEMModel) -> None:
    print(f"Rotor: {model.defn.name}, "
          f"R={model.defn.blade.radius_m} m, "
          f"sigma={model.defn.blade.n_blades * model.defn.blade.chord_m / (math.pi * model.defn.blade.radius_m):.4f}")
    print()
    for label, csv_name, pitch_deg, _ in TABLES:
        sub = survey.by_table(label)
        if not sub:
            print(f"  skip: {csv_name} (no rows)")
            continue
        cqs_t = np.array([c.CQ_trim for c in sub])
        cqs_nt = np.array([c.CQ_notrim for c in sub if c.CQ_notrim is not None])
        print(f"{csv_name}  (pitch {pitch_deg} deg, {len(sub)} rows):")
        if cqs_nt.size:
            _summarize("CQ no-trim", cqs_nt)
        _summarize("CQ with cyclic trim", cqs_t)
        print()

    print("=" * 70)
    print("ALL TABLES COMBINED:")
    cqs_t = survey.cq_trim()
    cqs_nt = survey.cq_notrim()
    if cqs_nt.size:
        _summarize("CQ no-trim", cqs_nt)
    if cqs_t.size:
        _summarize("CQ trimmed", cqs_t)
    print()
    print("Real rotor at autorotation has CQ == 0 by construction.")
    print("|CQ_BEM| is the residual the BEM produces at the same operating point.")
    print("For reference, hover CQ for a typical rotor is ~0.0003 - 0.0005.")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else "")
    p.add_argument("--sample", type=int, default=None, metavar="N",
                   help="evenly-sample N rows per table (smoke test); "
                        "omit to run the full ~469-row sweep")
    args = p.parse_args()

    model = load_model()
    if model is None:
        print(f"missing rotor fixture: {YAML}")
        return 2

    if args.sample is not None:
        print(f"** SAMPLED MODE: ~{args.sample} rows per table -- aggregate "
              f"stats below are NOT comparable to the 469-row envelope "
              f"quoted in test_wheatley_autorotation.py")
        print()

    t0 = time.time()
    survey = run_survey(model, sample=args.sample, include_notrim=True)
    _print_survey(survey, model=model)
    print(f"\n[walltime {time.time()-t0:.1f}s]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
