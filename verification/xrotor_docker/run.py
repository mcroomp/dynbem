"""Drive xrotor-python for a rotor described by a YAML.

XROTOR (Mark Drela, MIT, GPL) is one of the few open-source rotor
solvers that handles static-thrust hover cleanly -- its BEM iterates
on a hover-safe variable, not on the wind-turbine `a = v_a/V_inf` that
breaks down at V_inf = 0. We use the DARcorporation/xrotor-python
wrapper (PyPI: xrotor) which exposes a Case-based Python API.

YAML schema matches verification/openfast_docker/inputs/*.yaml as
closely as possible:

    name: <free-form label, used to name the output CSV>
    rotor:
      n_blades: int
      R_hub_m: float        # blade inboard cutout (m from rotor centre)
      R_tip_m: float
    air:
      rho:  float           # optional, defaults 1.225
      mu:   float           # optional, defaults 1.81206e-5
      a:    float           # speed of sound; optional, defaults 340
    blade_stations:
      - {r_m: float, chord_m: float, twist_deg: float}
      - ...
    airfoil_polar:
      Re: float             # informational
      points:
        - {alpha_deg: float, cl: float, cd: float, cm: float}
        - ...               # cm optional (defaults 0)
    operating_points:
      - {U_wind_ms: float, Omega_rpm: float, pitch_deg: float}
      - ...

XROTOR's geometry uses r/R and c/R (normalised by tip radius). We
convert from absolute m before handing off.

For hover (U_wind_ms ~= 0): XROTOR's free-tip potential-formulation
solver handles V_inf = 0 cleanly; we pass vel=1e-3 m/s rather than
zero strictly because the advance-ratio normalisation 1/adv blows up
at exactly zero, but the thrust at this V is indistinguishable from
true hover at typical rotor tip speeds.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np
import yaml

from xrotor import XRotor                                               # noqa: E402
from xrotor.model import Case                                          # noqa: E402


def build_case(cfg: dict, op: dict) -> Case:
    """Translate one (cfg, operating-point) pair into an XRotor Case."""
    r_tip = float(cfg["rotor"]["R_tip_m"])
    r_hub = float(cfg["rotor"]["R_hub_m"])
    n_blds = int(cfg["rotor"]["n_blades"])
    air = cfg.get("air", {})
    rho  = float(air.get("rho", 1.225))
    vso  = float(air.get("a", 340.0))
    rmu  = float(air.get("mu", 1.81206e-5))

    # XRotor expects r/R and c/R, twist already including collective pitch.
    stations = sorted(cfg["blade_stations"], key=lambda s: s["r_m"])
    radii = np.array([s["r_m"] / r_tip for s in stations], dtype=np.float32)
    chord = np.array([s["chord_m"] / r_tip for s in stations], dtype=np.float32)
    pitch_deg = float(op["pitch_deg"])
    # XRotor's "twist" column is the absolute local blade-element pitch
    # (blade twist + collective). For our YAML, twist_deg is the design
    # blade twist; we add the collective on top.
    twist = np.array([s["twist_deg"] + pitch_deg for s in stations],
                     dtype=np.float32)
    ubody = np.zeros_like(radii)

    # Polar -> (alpha_deg, cl, cd, cm) rows, keyed by r/R = 0 (single
    # global polar shared across radii, matching how AeroDyn/CCBlade are
    # being driven in the cross-checks).
    pts = cfg["airfoil_polar"]["points"]
    polar_arr = np.array(
        [[float(p["alpha_deg"]), float(p["cl"]), float(p["cd"]),
          float(p.get("cm", 0.0))] for p in pts],
        dtype=np.float32,
    )

    U = float(op["U_wind_ms"])
    Omega_rpm = float(op["Omega_rpm"])
    omega = Omega_rpm * math.pi / 30.0  # rad/s
    # XRotor's advance ratio is V/(Omega*R). For hover we use a tiny
    # but non-zero velocity to keep the advance-ratio normalisation
    # finite without changing the physics.
    vel = max(U, 1.0e-3)
    adv = vel / (omega * r_tip)

    return Case.from_dict({
        "conditions": {
            "rho": rho, "vso": vso, "rmu": rmu, "alt": 1,
            "vel": vel, "adv": adv,
        },
        "disk": {
            "n_blds": n_blds,
            "blade": {
                "geometry": {
                    "r_hub":  r_hub / r_tip,
                    "r_tip":  1.0,           # normalised; XRotor scales by r_tip * R
                    "r_wake": 0.0,
                    "rake":   0.0,
                    "radii":  radii,
                    "chord":  chord,
                    "twist":  twist,
                    "ubody":  ubody,
                },
                "polars": {
                    0.0: polar_arr,
                },
            },
        },
    })


def run_case(cfg: dict, op: dict) -> dict[str, float]:
    """Run one operating point through XRotor and return a CSV row."""
    xr = XRotor()
    xr.print = False
    xr.case = build_case(cfg, op)
    # Fixed-RPM operation. xr.operate(rpm=...) returns RMS of last iter.
    rms = xr.operate(rpm=float(op["Omega_rpm"]))
    perf = xr.performance
    return {
        "U_wind_ms":  float(op["U_wind_ms"]),
        "Omega_rpm":  float(op["Omega_rpm"]),
        "pitch_deg":  float(op["pitch_deg"]),
        "T_N":        perf.thrust,
        "Q_Nm":       perf.torque,
        "P_W":        perf.power,
        "efficiency": perf.efficiency,
        "rms":        float(rms),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--output-dir", default="/out")
    args = ap.parse_args()

    cfg = yaml.safe_load(Path(args.config).read_text())
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    name = cfg["name"]

    print(f"Cases: {len(cfg['operating_points'])}")
    rows: list[dict[str, float]] = []
    for i, op in enumerate(cfg["operating_points"]):
        print(f"  [{i+1}/{len(cfg['operating_points'])}] "
              f"U={op['U_wind_ms']} Omega={op['Omega_rpm']} pitch={op['pitch_deg']}",
              flush=True)
        try:
            row = run_case(cfg, op)
            rows.append(row)
        except Exception as exc:                                       # noqa: BLE001
            print(f"    ERROR: {exc}", file=sys.stderr)

    out_path = out_dir / f"{name}_xrotor.csv"
    if rows:
        cols = ["U_wind_ms", "Omega_rpm", "pitch_deg",
                "T_N", "Q_Nm", "P_W", "efficiency", "rms"]
        with out_path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        print(f"Wrote {out_path}")
    else:
        print("No rows produced.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
