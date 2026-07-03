"""Step 3: Steady operating-point benchmarks for De Schutter validation.

Runs a compact operating-point sweep and writes CSV outputs for review.
The sweep is deterministic and designed for quick regression checks.

Usage:
    uv run python deschutter/step3_steady_sweep.py
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from dynbem import RotorInputs, create_aero
from dynbem.rotor_definition import load as load_rotor
from dynbem.rotor_state import QuasiStaticRotorState


ROOT = Path(__file__).resolve().parent.parent
ROTOR_YAML = ROOT / "deschutter" / "rotor.yaml"
OUT_DIR = ROOT / "deschutter" / "out"
OUT_CSV = OUT_DIR / "steady_sweep.csv"
RHO = 1.225


@dataclass
class Row:
    u_wind_ms: float
    omega_rpm: float
    collective_deg: float
    thrust_N: float
    q_spin_Nm: float
    ct: float
    cq: float


def run_point(model, u_wind_ms: float, omega_rpm: float, collective_deg: float) -> Row:
    omega = omega_rpm * math.pi / 30.0
    inp = RotorInputs(
        collective_rad=math.radians(collective_deg),
        tilt_lon=0.0,
        tilt_lat=0.0,
        R_hub=np.eye(3),
        v_hub_world=np.zeros(3),
        wind_world=np.array([0.0, 0.0, -u_wind_ms]),
        rho_kg_m3=RHO,
        omega_rad_s=omega,
    )

    result, _ = model.compute_forces(inp, QuasiStaticRotorState())

    r = model.defn.blade.radius_m
    a = math.pi * r * r
    denom_f = RHO * a * (omega * r) ** 2
    denom_q = denom_f * r

    thrust_n = -float(result.F_world[2])
    q_spin_nm = float(result.Q_spin)

    return Row(
        u_wind_ms=u_wind_ms,
        omega_rpm=omega_rpm,
        collective_deg=collective_deg,
        thrust_N=thrust_n,
        q_spin_Nm=q_spin_nm,
        ct=thrust_n / denom_f,
        cq=q_spin_nm / denom_q,
    )


def main() -> int:
    model = create_aero(load_rotor(str(ROTOR_YAML)), model="bem")

    wind_grid = [6.0, 8.0, 10.0, 12.0]
    omega_grid = [220.0, 270.0, 320.0]
    collective_grid = [-8.0, -6.0, -4.0, -2.0, 0.0]

    rows: list[Row] = []
    for u in wind_grid:
        for omega_rpm in omega_grid:
            for collective_deg in collective_grid:
                rows.append(run_point(model, u, omega_rpm, collective_deg))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="", encoding="ascii") as f:
        w = csv.writer(f)
        w.writerow(["u_wind_ms", "omega_rpm", "collective_deg", "thrust_N", "q_spin_Nm", "CT", "CQ"])
        for r in rows:
            w.writerow([
                f"{r.u_wind_ms:.3f}",
                f"{r.omega_rpm:.3f}",
                f"{r.collective_deg:.3f}",
                f"{r.thrust_N:.9e}",
                f"{r.q_spin_Nm:.9e}",
                f"{r.ct:.9e}",
                f"{r.cq:.9e}",
            ])

    # Quick trend sanity: for each (u, omega), CT should increase with collective.
    failures = 0
    for u in wind_grid:
        for omega_rpm in omega_grid:
            series = [r for r in rows if r.u_wind_ms == u and r.omega_rpm == omega_rpm]
            series.sort(key=lambda x: x.collective_deg)
            ct_vals = np.array([r.ct for r in series], dtype=float)
            if np.any(np.diff(ct_vals) < -1e-8):
                failures += 1

    print(f"Wrote {OUT_CSV}")
    if failures == 0:
        print("Step 3 sweep complete. Monotonic CT trend checks passed.")
        return 0

    print(f"Step 3 sweep complete with {failures} trend warning(s).")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
