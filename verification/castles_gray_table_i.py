"""Validation of Level-1 BEM against NACA TN-2474 Table I.

Table I covers the full flight regime from hover through VRS to WBS/autorotation
for the 6-ft constant-chord, untwisted, 3-blade rotor.

Companion to tests/test_castles_gray.py.  The test imports `load_model`
and `bem_forces` from here as the single source of BEM-driver logic,
and `run_survey` for the Table I aggregate.  Sampled mode keeps the
test fast; the verifier's `main()` runs the full Table I sweep for
re-baselining bounds.

Usage
-----
    .venv/Scripts/python verification/castles_gray_table_i.py
    .venv/Scripts/python verification/castles_gray_table_i.py --sample 3

Regions
-------
  HOVER : lambda2 < 0.10
  DESC  : 0.10 <= lambda2 < 0.30
  VRS   : 0.30 <= lambda2 < 2.00  -- momentum theory expected to fail
  WBS   : lambda2 >= 2.00          -- Glauert/Buhl correction active

Sign convention (NED)
---------------------
  Descent in paper (positive lambda2) maps to v_climb < 0 (air flows upward).
  wind_world[2] = -V_descent for the BEM call.

Delta-CQ definition
-------------------
  DeltaCQ_paper = CQ(theta, V) - CQ_zero_thrust
  CQ_zero_thrust = CQ at theta=0, V=0 (profile drag at idle pitch, same RPM).
  DeltaCQ_model  = CQ_model(theta, V) - CQ_model(theta=0, V=0)
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from dynbem.bem import BEMModel
from dynbem import RotorInputs
import dynbem.rotor_definition as rotor_definition
from dynbem.rotor_state import QuasiStaticRotorState


ROOT = Path(__file__).resolve().parent.parent
ROTOR_YAML = str(ROOT / "rotors" / "castles_gray_6ft" / "rotor.yaml")
RHO = 1.225


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


def load_model() -> BEMModel:
    """Build the Castles-Gray 6-ft rotor BEM model from rotor.yaml."""
    return BEMModel(defn=rotor_definition.load(ROTOR_YAML))


def bem_forces(model: BEMModel, theta_deg: float, omega_rpm: float,
               v_climb_ms: float = 0.0) -> tuple[float, float, float]:
    """Run the BEM at one (theta, RPM, v_climb) point.

    Returns ``(CT, CQ, Q_spin_Nm)``.  v_climb_ms uses NED convention:

      v_climb > 0 : air flows downward through disk (climb / normal inflow)
      v_climb = 0 : hover
      v_climb < 0 : air flows upward through disk (autorotation / WBS)
    """
    omega = omega_rpm * math.pi / 30.0
    R = model.defn.blade.radius_m
    A = math.pi * R ** 2
    inp = RotorInputs(
        collective_rad=math.radians(theta_deg),
        tilt_lon=0.0, tilt_lat=0.0, R_hub=np.eye(3),
        v_hub_world=np.zeros(3),
        wind_world=np.array([0.0, 0.0, v_climb_ms]),
        t=0.0,
    )
    state = QuasiStaticRotorState(omega_rad_s=omega)
    result, _ = model.compute_forces(inp, state)
    T = -result.F_world[2]
    CT = T / (RHO * A * (omega * R) ** 2)
    CQ = result.Q_spin / (RHO * A * (omega * R) ** 2 * R)
    return CT, CQ, result.Q_spin


def region(lambda2: float) -> str:
    if lambda2 < 0.10:
        return "HOVER"
    if lambda2 < 0.30:
        return "DESC "
    if lambda2 < 2.00:
        return "VRS  "
    return "WBS  "


# ---------------------------------------------------------------------------
# Table I data
# Each row: (V_over_OR, theta_075R_deg, delta_CQ_measured, lambda2)
# delta_CQ_measured = None where paper left blank.
# ---------------------------------------------------------------------------

TABLE_I = [
    dict(run=3,  CT=0.004, rpm=1200, rows=[
        (0.0000, 8.83, 0.000282, 0.00),
        (0.0400, 8.74, 0.000249, 0.91),
        (0.0556, 8.83, 0.000251, 1.27),
        (0.0658, 8.83, 0.000266, 1.50),
        (0.0762, 4.87, 0.000009, 1.83),
        (0.0795, 3.99, 0.000039, 1.74),
        (0.0799, 4.10, 0.000009, 1.82),
    ]),
    dict(run=4,  CT=0.004, rpm=1200, rows=[
        (0.0000, 8.66, 0.000285, 0.00),
        (0.0168, 8.40, 0.000270, 0.38),
        (0.0254, 8.49, 0.000262, 0.58),
        (0.0400, 8.32, 0.000243, 0.91),
        (0.0550, 8.37, 0.000239, 1.26),
        (0.0649, 8.54, 0.000256, 1.48),
        (0.0702, 7.24, 0.000177, 1.61),
        (0.0734, 4.79, 0.000046, 1.67),
        (0.0790, 4.10, 0.000018, 1.81),
        (0.0859, 2.87, -0.000024, 1.96),
        (0.0949, 1.23, -0.000112, 2.17),
    ]),
    dict(run=5,  CT=0.004, rpm=1200, rows=[
        (0.0000, 8.86, 0.000292, 0.00),
        (0.0176, 8.70, 0.000279, 0.40),
        (0.0265, 8.53, 0.000270, 0.60),
        (0.0388, 8.40, 0.000245, 0.89),
        (0.0552, 8.53, 0.000249, 1.26),
        (0.0617, 8.48, 0.000256, 1.41),
        (0.0670, 8.46, 0.000250, 1.53),
        (0.0734, 4.60, 0.000039, 1.67),
        (0.0790, 3.99, 0.000011, 1.81),
        (0.0856, 2.58, -0.000039, 1.96),
        (0.0949, 1.18, -0.000116, 2.17),
    ]),
    dict(run=6,  CT=0.002, rpm=1200, rows=[
        (0.0000, 5.33, None, 0.00),
        (0.0190, 4.95, None, 0.61),
        (0.0300, 4.75, None, 0.97),
        (0.0382, 4.92, None, 1.23),
        (0.0447, 5.16, None, 1.44),
        (0.0490, 4.26, None, 1.58),
        (0.0556, 1.44, None, 1.80),
        (0.0596, 0.56, None, 1.92),
    ]),
    dict(run=9,  CT=0.002, rpm=1600, rows=[
        (0.0000, 5.13, 0.000062, 0.00),
        (0.0136, 5.26, 0.000072, 0.43),
        (0.0214, 4.87, 0.000067, 0.69),
        (0.0303, 4.88, 0.000065, 0.98),
        (0.0373, 4.81, 0.000067, 1.20),
        (0.0417, 5.03, 0.000077, 1.35),
        (0.0463, 5.01, 0.000081, 1.50),
        (0.0501, 3.44, 0.000045, 1.62),
        (0.0525, 1.94, 0.000030, 1.69),
        (0.0561, 1.27, -0.000002, 1.81),
        (0.0602, 0.68, -0.000028, 1.94),
        (0.0650, -0.11, -0.000050, 2.10),
        (0.0720, -1.32, -0.000072, 2.32),
        (0.0763, -1.66, -0.000084, 2.46),
    ]),
    dict(run=14, CT=0.005, rpm=1600, rows=[
        (0.0000, 10.06, 0.000414, 0.00),
        (0.0103,  9.99, 0.000409, 0.21),
        (0.0215,  9.81, 0.000400, 0.43),
        (0.0287, 10.02, 0.000396, 0.57),
        (0.0386, 10.13, 0.000397, 0.77),
        (0.0510, 10.36, None,     1.02),
        (0.0548, 10.69, None,     1.10),
        (0.0586, 10.96, None,     1.17),
        (0.0656, 10.92, None,     1.31),
        (0.0722,  9.80, None,     1.44),
        (0.0762,  9.80, None,     1.52),
        (0.0782,  9.79, 0.000377, 1.56),
        (0.0818,  5.07, 0.000091, 1.64),
        (0.0888,  5.49, 0.000086, 1.78),
        (0.0898,  4.89, 0.000051, 1.80),
        (0.0958,  3.52, -0.000016, 1.92),
    ]),
    dict(run=32, CT=0.002, rpm=1200, rows=[
        (0.0000, 5.32, 0.000087, 0.00),
        (0.0048, 5.00, 0.000087, 0.15),
        (0.0161, 5.21, 0.000077, 0.50),
        (0.0298, 4.91, 0.000075, 0.92),
        (0.0376, 5.22, 0.000080, 1.15),
        (0.0442, 5.22, 0.000086, 1.36),
        (0.0484, 5.22, 0.000087, 1.49),
        (0.0552, 1.90, 0.000007, 1.70),
        (0.0592, 1.39, -0.000002, 1.82),
        (0.0641, 0.42, -0.000025, 1.97),
        (0.0676, -0.28, -0.000042, 2.08),
        (0.0721, -0.84, -0.000092, 2.21),
        (0.0771, -1.46, -0.000069, 2.37),  # Flag F
    ]),
    dict(run=34, CT=0.002, rpm=1600, rows=[
        (0.0000, 5.33, 0.000087, 0.00),
        (0.0043, 5.12, 0.000083, 0.36),
        (0.0113, 5.12, 0.000084, 0.61),
        (0.0191, 4.81, 0.000073, 0.91),
        (0.0282, 4.90, 0.000077, 1.13),
        (0.0350, 5.22, 0.000080, 1.27),
        (0.0395, 5.21, 0.000082, 1.43),
        (0.0444, 3.89, 0.000065, 1.55),
        (0.0481, 3.89, 0.000055, 1.63),
        (0.0504, 1.88, 0.000006, 1.72),
        (0.0533, 1.35, -0.000007, 1.84),
        (0.0572, 0.55, -0.000028, 2.01),
        (0.0624, -0.30, -0.000053, 2.22),
        (0.0690, -0.77, -0.000078, 2.39),  # Flag C
    ]),
    dict(run=35, CT=0.004, rpm=1600, rows=[
        (0.0000, 8.68, 0.000271, 0.00),
        (0.0098, 8.64, 0.000270, 0.22),
        (0.0162, 8.50, 0.000261, 0.37),
        (0.0215, 8.50, 0.000261, 0.49),
        (0.0295, 8.45, 0.000248, 0.67),
        (0.0411, 8.30, 0.000234, 0.94),
        (0.0485, 8.30, 0.000229, 1.10),
        (0.0513, 8.30, 0.000229, 1.17),
        (0.0566, 8.30, 0.000239, 1.29),
        (0.0624, 8.30, 0.000243, 1.42),
        (0.0684, 7.22, 0.000175, 1.56),
        (0.0739, 4.36, 0.000048, 1.69),
        (0.0795, 3.67, 0.000020, 1.81),
        (0.0835, 3.29, 0.000001, 1.91),
        (0.0881, 2.46, -0.000045, 2.01),
    ]),
    dict(run=36, CT=0.004, rpm=1200, rows=[
        (0.0000, 8.49, 0.000242, 0.00),   # Flag D
        (0.0066, 8.68, 0.000247, 0.15),
        (0.0163, 8.82, 0.000254, 0.37),
        (0.0238, 8.74, 0.000245, 0.54),
        (0.0353, 8.68, 0.000232, 0.81),
        (0.0459, 8.51, 0.000215, 1.05),
        (0.0525, 8.72, 0.000224, 1.20),
        (0.0590, 8.62, 0.000233, 1.35),
        (0.0652, 8.80, 0.000247, 1.49),
        (0.0684, 8.00, 0.000194, 1.56),
        (0.0715, 6.34, 0.000106, 1.63),
        (0.0764, 5.25, 0.000054, 1.75),
        (0.0831, 4.29, 0.000016, 1.90),
        (0.0914, 2.81, -0.000047, 2.09),
    ]),
    dict(run=38, CT=0.005, rpm=1600, rows=[
        (0.0000, 10.19, None, 0.00),
        (0.0096, 10.18, None, 0.20),
        (0.0161, 10.18, None, 0.33),
        (0.0210, 10.02, None, 0.43),
        (0.0271, 10.02, None, 0.55),
        (0.0364, 10.14, None, 0.74),
        (0.0496, 10.41, None, 1.01),
        (0.0535, 10.22, None, 1.09),
        (0.0569, 10.14, None, 1.16),
        (0.0602, 10.14, None, 1.23),
        (0.0629,  9.98, None, 1.28),
        (0.0662,  9.98, None, 1.35),
        (0.0695,  9.84, None, 1.42),
        (0.0745, 10.22, None, 1.52),
        (0.0765,  9.35, 0.000172, 1.56),
        (0.0801,  6.58, 0.000146, 1.63),
        (0.0819,  6.45, 0.000130, 1.71),
        (0.0883,  5.22, 0.000094, 1.80),
        (0.0940,  4.77, 0.000038, 1.92),
        (0.1032,  3.39, -0.000060, 2.11),
    ]),
]


@dataclass
class Comparison:
    run: int
    region: str
    rpm: int
    theta_deg: float
    v_over_OR: float
    lambda2: float
    CT_nom: float
    CT_pred: float
    dCT_pct: float
    dCQ_meas: float | None
    dCQ_pred: float
    ddCQ_pct: float | None


@dataclass
class Survey:
    comparisons: list[Comparison] = field(default_factory=list)
    sample: int | None = None

    def by_region(self, reg: str) -> list[Comparison]:
        return [c for c in self.comparisons if c.region.strip() == reg.strip()]

    def by_run(self, run: int) -> list[Comparison]:
        return [c for c in self.comparisons if c.run == run]


def run_survey(model: BEMModel | None = None,
               sample: int | None = None) -> Survey:
    """Sweep Table I; sample-N picks evenly-spaced rows from each run.

    sample=None -> every row.  Returns a Survey of Comparison records.
    """
    if model is None:
        model = load_model()
    R = model.defn.blade.radius_m
    survey = Survey(sample=sample)

    cq_base: dict[int, float] = {}
    for rpm in (1200, 1600):
        _, cq0, _ = bem_forces(model, 0.0, float(rpm), 0.0)
        cq_base[rpm] = cq0

    for run_data in TABLE_I:
        run = run_data["run"]
        CT_nom = run_data["CT"]
        rpm = run_data["rpm"]
        omega = rpm * math.pi / 30.0
        for (v_or, theta, dcq_meas, lam2) in _sample_evenly(run_data["rows"], sample):
            v_descent = v_or * omega * R
            v_climb_ms = -v_descent
            CT_pred, CQ_pred, _ = bem_forces(model, theta, float(rpm), v_climb_ms)
            dCQ_pred = CQ_pred - cq_base[rpm]
            dCT_pct = (CT_pred - CT_nom) / CT_nom * 100.0 if CT_nom else float("nan")
            ddCQ_pct: float | None
            if dcq_meas is not None and abs(dcq_meas) > 1e-9:
                ddCQ_pct = (dCQ_pred - dcq_meas) / abs(dcq_meas) * 100.0
            else:
                ddCQ_pct = None
            survey.comparisons.append(Comparison(
                run=run, region=region(lam2), rpm=rpm,
                theta_deg=theta, v_over_OR=v_or, lambda2=lam2,
                CT_nom=CT_nom, CT_pred=CT_pred, dCT_pct=dCT_pct,
                dCQ_meas=dcq_meas, dCQ_pred=dCQ_pred, ddCQ_pct=ddCQ_pct,
            ))
    return survey


def _print_survey(survey: Survey) -> None:
    hdr = (
        f"{'Run':>3}  {'CT':>5}  {'RPM':>4}  "
        f"{'lam2':>5}  {'theta':>6}  "
        f"{'CT_nom':>7}  {'CT_pred':>7}  {'dCT%':>6}  "
        f"{'dCQ_meas':>10}  {'dCQ_pred':>10}  {'ddCQ%':>7}  "
        f"{'region'}"
    )
    sep = "-" * len(hdr)
    print(hdr)
    print(sep)

    stats = {r: {"n": 0, "n_cq": 0, "sumsq_ct": 0.0, "sumsq_cq": 0.0}
             for r in ("HOVER", "DESC ", "VRS  ", "WBS  ")}
    last_run = None
    for c in survey.comparisons:
        if last_run is not None and c.run != last_run:
            print(sep)
        last_run = c.run
        dcq_meas_str = f"{c.dCQ_meas:+.6f}" if c.dCQ_meas is not None else "    --     "
        ddCQ_str = f"{c.ddCQ_pct:+7.1f}" if c.ddCQ_pct is not None else "   --  "
        print(
            f"{c.run:>3}  {c.CT_nom:.3f}  {c.rpm:>4}  "
            f"{c.lambda2:>5.2f}  {c.theta_deg:>6.2f}  "
            f"{c.CT_nom:>7.5f}  {c.CT_pred:>7.5f}  {c.dCT_pct:>+6.1f}  "
            f"{dcq_meas_str:>10}  {c.dCQ_pred:>+10.6f}  {ddCQ_str:>7}  "
            f"{c.region}"
        )
        s = stats[c.region]
        s["n"] += 1
        s["sumsq_ct"] += c.dCT_pct ** 2
        if c.ddCQ_pct is not None and c.dCQ_meas is not None and abs(c.dCQ_meas) > 5e-6:
            s["n_cq"] += 1
            s["sumsq_cq"] += c.ddCQ_pct ** 2
    print(sep)

    print()
    print("RMSE summary (CT error and delta-CQ error by region)")
    print(f"{'Region':<8}  {'n_CT':>5}  {'RMSE_CT%':>9}  {'n_CQ':>5}  {'RMSE_dCQ%':>10}")
    print("-" * 48)
    for reg, s in stats.items():
        rmse_ct = math.sqrt(s["sumsq_ct"] / s["n"]) if s["n"] else float("nan")
        rmse_cq = math.sqrt(s["sumsq_cq"] / s["n_cq"]) if s["n_cq"] else float("nan")
        n_cq_str = str(s["n_cq"]) if s["n_cq"] else "--"
        rmse_cq_str = f"{rmse_cq:9.1f}" if not math.isnan(rmse_cq) else "       --"
        print(f"{reg:<8}  {s['n']:>5}  {rmse_ct:>9.1f}  {n_cq_str:>5}  {rmse_cq_str}")

    print()
    print("VRS note: momentum-theory BEM is not expected to predict VRS accurately.")
    print("WBS note: Glauert/Buhl correction active; model should track data.")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else "")
    p.add_argument("--sample", type=int, default=None, metavar="N",
                   help="evenly-sample N rows per run (smoke test); "
                        "omit to run the full Table I sweep")
    args = p.parse_args()

    if args.sample is not None:
        print(f"** SAMPLED MODE: ~{args.sample} rows per run -- region RMSEs below "
              f"are NOT comparable to the full-sweep numbers used to set "
              f"thresholds in test_castles_gray.py")
        print()

    _print_survey(run_survey(sample=args.sample))
    return 0


if __name__ == "__main__":
    sys.exit(main())
