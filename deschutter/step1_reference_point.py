"""Step 1: De Schutter reference operating point using dynbem.

This script sets up one canonical operating point and prints normalized
outputs. It is intended to verify sign conventions, frame mapping, and
basic magnitude sanity before larger sweeps.

Usage:
    uv run python deschutter/step1_reference_point.py
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from dynbem import RotorInputs, create_aero
from dynbem.rotor_definition import load as load_rotor
from dynbem.rotor_state import QuasiStaticRotorState


ROOT = Path(__file__).resolve().parent.parent
ROTOR_YAML = ROOT / "deschutter" / "rotor.yaml"
RHO = 1.225


def run_reference_point() -> dict[str, float]:
    model = create_aero(load_rotor(str(ROTOR_YAML)), model="bem")

    # Reference point chosen to match the design-scale order of magnitude
    # discussed in De Schutter style studies.
    omega_rpm = 270.0
    omega = omega_rpm * math.pi / 30.0
    u_wind = 10.0

    inputs = RotorInputs(
        collective_rad=math.radians(-6.0),
        tilt_lon=0.0,
        tilt_lat=0.0,
        R_hub=np.eye(3),
        v_hub_world=np.zeros(3),
        wind_world=np.array([0.0, 0.0, -u_wind]),
        t=0.0,
        rho_kg_m3=RHO,
        omega_rad_s=omega,
    )

    result, _ = model.compute_forces(inputs, QuasiStaticRotorState())

    r = model.defn.blade.radius_m
    a = math.pi * r * r
    denom_f = RHO * a * (omega * r) ** 2
    denom_q = denom_f * r

    thrust_n = -float(result.F_world[2])
    q_spin_nm = float(result.Q_spin)
    ct = thrust_n / denom_f
    cq = q_spin_nm / denom_q

    return {
        "omega_rpm": omega_rpm,
        "u_wind_ms": u_wind,
        "collective_deg": -6.0,
        "thrust_N": thrust_n,
        "Q_spin_Nm": q_spin_nm,
        "CT": ct,
        "CQ": cq,
    }


def main() -> int:
    out = run_reference_point()
    print(json.dumps(out, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
