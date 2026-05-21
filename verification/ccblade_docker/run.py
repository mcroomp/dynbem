"""Generate CCBlade reference results from a YAML config.

The YAML schema is intentionally minimal so a new rotor can be added by
dropping one file into `verification/ccblade_docker/inputs/` and
re-running the container. Schema:

    name: <free-form label, used to name the output CSV>
    rotor:
      n_blades: int
      R_hub_m: float        # radius at the start of the airfoil section
      R_tip_m: float
      precone_deg: float    # optional, defaults 0
      tilt_deg: float       # optional, defaults 0
    air:
      rho: float            # optional, defaults 1.225
      mu: float             # optional, defaults 1.81206e-5
    blade_stations:
      - {r_m: float, chord_m: float, twist_deg: float}
      - ...                 # one row per radial station, sorted by r_m
    airfoil_polar:
      Re: float
      points:
        - {alpha_deg: float, cl: float, cd: float}
        - ...
    operating_points:
      - {U_wind_ms: float, Omega_rpm: float, pitch_deg: float}
      - ...

The output CSV has one row per operating point with columns:
    U_wind_ms, Omega_rpm, pitch_deg, T_N, Q_Nm, P_W, CP, CT, CQ.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import yaml
from wisdem.ccblade.ccblade import CCAirfoil, CCBlade


def build_model(cfg: dict) -> CCBlade:
    stations = sorted(cfg["blade_stations"], key=lambda s: s["r_m"])
    r = np.array([s["r_m"] for s in stations])
    chord = np.array([s["chord_m"] for s in stations])
    twist_deg = np.array([s["twist_deg"] for s in stations])

    polar = cfg["airfoil_polar"]
    alpha = np.array([p["alpha_deg"] for p in polar["points"]])
    cl = np.array([p["cl"] for p in polar["points"]])
    cd = np.array([p["cd"] for p in polar["points"]])
    af = CCAirfoil(alpha, [polar["Re"]], cl, cd)
    af_per_station = [af] * len(r)

    rotor = cfg["rotor"]
    air = cfg.get("air", {})
    return CCBlade(
        r, chord, twist_deg, af_per_station,
        Rhub=rotor["R_hub_m"], Rtip=rotor["R_tip_m"],
        B=rotor["n_blades"],
        rho=air.get("rho", 1.225),
        mu=air.get("mu", 1.81206e-5),
        precone=rotor.get("precone_deg", 0.0),
        tilt=rotor.get("tilt_deg", 0.0),
        yaw=0.0, shearExp=0.0, hubHt=10.0,
        nSector=1,
    )


def run(cfg: dict) -> list[dict]:
    model = build_model(cfg)
    rows: list[dict] = []
    for op in cfg["operating_points"]:
        U = float(op["U_wind_ms"])
        Omega = float(op["Omega_rpm"])
        pitch = float(op["pitch_deg"])
        # evaluate accepts arrays; pass length-1 arrays and unpack.
        out, deriv = model.evaluate(
            np.array([U]), np.array([Omega]), np.array([pitch]),
            coefficients=True,
        )
        rows.append({
            "U_wind_ms": U,
            "Omega_rpm": Omega,
            "pitch_deg": pitch,
            "T_N":  float(out["T"][0]),
            "Q_Nm": float(out["Q"][0]),
            "P_W":  float(out["P"][0]),
            "CP":   float(out["CP"][0]),
            "CT":   float(out["CT"][0]),
            "CQ":   float(out["CQ"][0]),
        })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else "")
    ap.add_argument("--config", required=True,
                    help="path to YAML config (typically /in/<name>.yaml)")
    ap.add_argument("--out-dir", default="/out",
                    help="directory to write the CSV into (default /out)")
    args = ap.parse_args()

    cfg_path = Path(args.config)
    cfg = yaml.safe_load(cfg_path.read_text())
    name = cfg.get("name") or cfg_path.stem
    rows = run(cfg)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{name}.csv"
    with out_path.open("w", encoding="ascii", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out_path}  ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
