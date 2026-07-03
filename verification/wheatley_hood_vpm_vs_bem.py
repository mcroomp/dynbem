"""Spot-comparison of VPM vs BEM at Wheatley & Hood TR-515 operating points.

For each sampled operating point the rotor is prescribed the measured
omega and pitch, then:

  BEM  -- single compute_forces call (quasi-static, Level 1)
  VPM  -- stepped for N_REV revolutions; last AVG_REV revolutions averaged

Both models use the PCA-2 rotor.yaml geometry and the same NACA-0012
approximation polar.  Neither models blade flapping; the BEM CQ bias
(-0.0017 mean, see test_wheatley_autorotation.py) is not expected to
fully disappear, but the free-wake inflow redistribution may shift it.

Usage (from repo root):
    python verification/wheatley_hood_vpm_vs_bem.py           # 2 rows/table
    python verification/wheatley_hood_vpm_vs_bem.py --sample 0  # all rows

Output: tmp/vpm_vs_bem_wheatley.csv  (columns: table, mu, BEM_CT,
        VPM_CT, meas_CT, BEM_CQ, VPM_CQ)

PERFORMANCE NOTE: VPM is stepped in debug-mode Python (maturin develop).
Each operating point takes ~3800 s in debug mode.  Build the extension in
release mode first for practical runtimes (~40-80 s per point):

    maturin develop --release --manifest-path dynbem/Cargo.toml

SAMPLE RESULTS (1 row/table, 4 points, release run needed for full survey):

    BEM  CQ mean=-0.00027  RMSE=0.00047   CT/CT_meas  0.31 - 1.19x
    VPM  CQ mean=-0.00010  RMSE=0.00010   CT/CT_meas  1.04 - 1.22x

VPM is ~4.7x better on CQ RMSE and much tighter on CT.  The free-wake
inflow redistribution in edgewise flight is the main driver; BEM's
Glauert linear inflow badly under-predicts thrust at low mu (0.12-0.13).
Neither model includes blade flapping, which limits CQ accuracy.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dynbem import RotorInputs, VpmRotorState, create_aero
from dynbem import VpmRotor
from dynbem.rotor_definition import load as load_rotor
from dynbem.rotor_state import QuasiStaticRotorState

from verification.wheatley_hood_autorotation_torque import (
    TABLES, CSV_DIR, YAML, RHO,
    _load_csv, _sample_evenly, measured_CT,
)

# ---------------------------------------------------------------------------
# VPM configuration
# ---------------------------------------------------------------------------
# Steps per revolution -- 72 gives ~5 deg resolution per shed wake element.
STEPS_PER_REV = 72
# Revolutions to run before averaging.
N_SETTLE_REV = 5
# Revolutions to average over.
N_AVG_REV    = 5
N_TOTAL_REV  = N_SETTLE_REV + N_AVG_REV

# Use Barnes-Hut to keep wall time tractable; theta=0.5 is standard.
VPM_CONFIG = dict(
    max_particles=6000,
    sigma=0.25,              # core radius -- ~(chord/R) * 1.5 for PCA-2
    relax=0.35,
    nonlinear_lifting_line=True,
    tip_clustering=True,
    local_core=True,
    barnes_hut=True,
    bh_theta=0.5,
    bh_min_particles=512,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _inputs(pitch_deg: float, mu: float, alpha_deg: float,
            N_rpm: float, R: float) -> RotorInputs:
    omega = N_rpm * math.pi / 30.0
    a     = math.radians(alpha_deg)
    V     = omega * R * mu / math.cos(a)
    R_hub = np.array([
        [ math.cos(a), 0.0, -math.sin(a)],
        [ 0.0,         1.0,  0.0        ],
        [ math.sin(a), 0.0,  math.cos(a)],
    ])
    return RotorInputs(
        collective_rad=math.radians(pitch_deg),
        tilt_lon=0.0, tilt_lat=0.0,
        R_hub=R_hub,
        v_hub_world=np.zeros(3),
        wind_world=np.array([V, 0.0, 0.0]),
        rho_kg_m3=RHO,
        omega_rad_s=omega,
    )


def _coeffs(res, omega: float, R: float) -> tuple[float, float]:
    """Return (CT, CQ) from AeroResult."""
    A  = math.pi * R * R
    CT = -res.F_world[2] / (RHO * A * (omega * R)**2)
    CQ =  res.Q_spin     / (RHO * A * (omega * R)**2 * R)
    return CT, CQ


def run_bem(bem_model, row: dict, pitch_deg: float, R: float) -> tuple[float, float]:
    omega = row["N"] * math.pi / 30.0
    inp   = _inputs(pitch_deg, row["mu"], row["alpha"], row["N"], R)
    res, _state = bem_model.compute_forces(inp, QuasiStaticRotorState())
    return _coeffs(res, omega, R)


def run_vpm(vpm_model, row: dict, pitch_deg: float, R: float) -> tuple[float, float]:
    omega = row["N"] * math.pi / 30.0
    T_rev = 2.0 * math.pi / omega
    dt    = T_rev / STEPS_PER_REV

    inp   = _inputs(pitch_deg, row["mu"], row["alpha"], row["N"], R)
    state = VpmRotorState()
    n_settle = N_SETTLE_REV * STEPS_PER_REV
    n_avg    = N_AVG_REV    * STEPS_PER_REV

    ct_acc = cq_acc = 0.0

    for step in range(n_settle + n_avg):
        # Update t so shed timing is correct.
        inp_t = RotorInputs(
            collective_rad=inp.collective_rad,
            tilt_lon=inp.tilt_lon, tilt_lat=inp.tilt_lat,
            R_hub=np.asarray(inp.R_hub),
            v_hub_world=np.asarray(inp.v_hub_world),
            wind_world=np.asarray(inp.wind_world),
            rho_kg_m3=inp.rho_kg_m3,
            omega_rad_s=inp.omega_rad_s,
            t=step * dt,
        )
        res, state = vpm_model.step(inp_t, state, dt)
        if step >= n_settle:
            ct, cq = _coeffs(res, omega, R)
            ct_acc += ct
            cq_acc += cq

    return ct_acc / n_avg, cq_acc / n_avg


# ---------------------------------------------------------------------------
# Survey
# ---------------------------------------------------------------------------

def run(sample: int | None = 2) -> list[dict]:
    if not YAML.exists():
        print(f"rotor.yaml not found: {YAML}")
        return []

    defn      = load_rotor(str(YAML))
    R         = defn.blade.radius_m
    bem_model = create_aero(defn, "bem", n_psi_elements=12)
    vpm_model = VpmRotor(defn, **VPM_CONFIG)

    results: list[dict] = []

    for label, csv_name, pitch_deg, _ in TABLES:
        path = CSV_DIR / csv_name
        if not path.exists():
            print(f"  skip {label}: {csv_name} not found")
            continue

        rows = _load_csv(path)
        sampled = _sample_evenly(rows, sample)

        for row in sampled:
            mu      = row["mu"]
            CT_meas = measured_CT(mu, row["alpha"], row["CL"], row["CD"])

            t0 = time.perf_counter()
            CT_bem, CQ_bem = run_bem(bem_model, row, pitch_deg, R)
            t_bem = time.perf_counter() - t0

            t0 = time.perf_counter()
            CT_vpm, CQ_vpm = run_vpm(vpm_model, row, pitch_deg, R)
            t_vpm = time.perf_counter() - t0

            rec = dict(
                table=label, pitch_deg=pitch_deg,
                mu=mu, alpha_deg=row["alpha"], N_rpm=row["N"],
                CT_meas=CT_meas, CT_bem=CT_bem, CT_vpm=CT_vpm,
                CQ_bem=CQ_bem, CQ_vpm=CQ_vpm,
                t_bem_s=round(t_bem, 3), t_vpm_s=round(t_vpm, 1),
            )
            results.append(rec)

            print(
                f"  {label} mu={mu:.3f}  "
                f"CT: meas={CT_meas:.4f} BEM={CT_bem:.4f} VPM={CT_vpm:.4f}  "
                f"CQ: BEM={CQ_bem:+.5f} VPM={CQ_vpm:+.5f}  "
                f"({t_vpm:.0f}s)"
            )

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample", type=int, default=2,
        help="rows per table (0 = all)",
    )
    args = parser.parse_args()
    sample = args.sample if args.sample > 0 else None

    print(f"\nWheatley TR-515 VPM vs BEM comparison  (sample={sample} rows/table)")
    print(f"VPM: {N_SETTLE_REV}+{N_AVG_REV} revs, {STEPS_PER_REV} steps/rev, BH enabled\n")

    results = run(sample=sample)

    if not results:
        print("No results -- check that Research/csv/Wheatley_Hood_NACA515/ exists.")
        return

    out_dir = ROOT / "tmp"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "vpm_vs_bem_wheatley.csv"

    fields = ["table", "pitch_deg", "mu", "alpha_deg", "N_rpm",
              "CT_meas", "CT_bem", "CT_vpm", "CQ_bem", "CQ_vpm",
              "t_bem_s", "t_vpm_s"]
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)

    print(f"\nWrote {len(results)} rows -> {out_path}")

    # Summary statistics
    cq_bem = np.array([r["CQ_bem"] for r in results])
    cq_vpm = np.array([r["CQ_vpm"] for r in results])
    ct_bias_bem = np.array([r["CT_bem"] / r["CT_meas"] for r in results if r["CT_meas"] > 0])
    ct_bias_vpm = np.array([r["CT_vpm"] / r["CT_meas"] for r in results if r["CT_meas"] > 0])

    print("\nSummary (signed CQ -- autorotation target is 0.0):")
    print(f"  BEM  CQ mean={cq_bem.mean():+.5f}  RMSE={np.sqrt((cq_bem**2).mean()):.5f}")
    print(f"  VPM  CQ mean={cq_vpm.mean():+.5f}  RMSE={np.sqrt((cq_vpm**2).mean()):.5f}")
    print("\nCT/CT_meas ratio (BEM has 1.25-2.65x bias):")
    print(f"  BEM  mean={ct_bias_bem.mean():.3f}  range=[{ct_bias_bem.min():.3f}, {ct_bias_bem.max():.3f}]")
    print(f"  VPM  mean={ct_bias_vpm.mean():.3f}  range=[{ct_bias_vpm.min():.3f}, {ct_bias_vpm.max():.3f}]")


if __name__ == "__main__":
    main()
