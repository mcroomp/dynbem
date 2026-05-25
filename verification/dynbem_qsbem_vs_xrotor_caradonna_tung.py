"""Compare dynbem.QuasiStaticBEM against XROTOR in pure hover on the
Caradonna-Tung rotor (2-blade NACA 0012, R = 1.143 m, 1250 rpm,
pitch 5/8/12 deg).

This is the hover-regime validation we could not do against AeroDyn
(its wind-turbine BEM falls over at V_inf = 0). XROTOR's free-tip
potential formulation handles static thrust natively, so this is the
first peer cross-check of dynbem's hover BEM against another
open-source tool.

The XROTOR reference is produced by:

    cd verification/xrotor_docker
    docker compose run --rm xrotor --config /in/caradonna_tung_hover.yaml

Three numbers per pitch are reported:
  - T_paper:  the Caradonna-Tung published CT scaled to thrust
  - T_xrotor: from the cached CSV above
  - T_dynbem: computed live from dynbem.QuasiStaticBEM in hover

so we can see where each code lands relative to the others.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dynbem import RotorInputs                                          # noqa: E402
from dynbem.bem import BEMModel                                          # noqa: E402
from dynbem.rotor_state import QuasiStaticRotorState                    # noqa: E402
from dynbem.rotor_definition import (                                    # noqa: E402
    BladeGeometry, LinearPolarParameters, RotorDefinition,
)
from dynbem.polar import TabulatedPolar                                 # noqa: E402

DEFAULT_YAML  = ROOT / "verification" / "xrotor_docker" / "inputs" / "caradonna_tung_hover.yaml"
DEFAULT_XCSV  = ROOT / "verification" / "data" / "caradonna_tung_hover_xrotor.csv"

# Published CT from Caradonna-Tung TM-81232 (paper figures 3,4,5):
#   pitch 5 deg  -> CT in 0.00213..0.00218 -> midpoint 0.00215
#   pitch 8 deg  -> CT in 0.00455..0.00473 -> midpoint 0.00465
#   pitch 12 deg -> CT in 0.00792..0.00807 -> midpoint 0.00800
PAPER_CT = {5.0: 0.00215, 8.0: 0.00465, 12.0: 0.00800}
RHO = 1.225


def build_dynbem_model(cfg: dict) -> BEMModel:
    """Translate the YAML into a dynbem QuasiStaticBEM model."""
    r_tip = float(cfg["rotor"]["R_tip_m"])
    stations = sorted(cfg["blade_stations"], key=lambda s: s["r_m"])
    blade = BladeGeometry(
        n_blades       = int(cfg["rotor"]["n_blades"]),
        radius_m       = r_tip,
        root_cutout_m  = float(cfg["rotor"]["R_hub_m"]),
        chord_m        = float(stations[0]["chord_m"]),  # constant chord here
        twist_deg      = 0.0,
        n_elements     = len(stations),
        r_stations_m       = [s["r_m"]       for s in stations],
        chord_stations_m   = [s["chord_m"]   for s in stations],
        twist_stations_deg = [s["twist_deg"] for s in stations],
    )
    airfoil = LinearPolarParameters(
        Re_design        = int(float(cfg["airfoil_polar"]["Re"])),
        CL0              = 0.0,
        CL_alpha_per_rad = 2.0 * math.pi,  # placeholder; we override w/ polar
        CD0              = 0.008,
        alpha_stall_deg  = 12.0,
    )
    pts = cfg["airfoil_polar"]["points"]
    polar = TabulatedPolar(
        alpha_rad=np.array([math.radians(p["alpha_deg"]) for p in pts]),
        cl       =np.array([p["cl"] for p in pts]),
        cd       =np.array([p["cd"] for p in pts]),
    )
    defn = RotorDefinition(blade=blade, airfoil=airfoil)
    return BEMModel(defn, polar, 36)


def dynbem_hover_thrust(model: BEMModel, pitch_deg: float,
                        omega_rpm: float, rho: float) -> tuple[float, float]:
    """Return (T, Q) at hover for the given collective + RPM."""
    omega = omega_rpm * math.pi / 30.0
    inp = RotorInputs(
        collective_rad=math.radians(pitch_deg),
        tilt_lon=0.0, tilt_lat=0.0,
        R_hub=np.eye(3),
        v_hub_world=np.zeros(3),
        wind_world=np.zeros(3),  # hover: no axial wind
        t=0.0,
        rho_kg_m3=rho,
        omega_rad_s=omega,
    )
    state = QuasiStaticRotorState()
    result, _ = model.compute_forces(inp, state)
    T = -float(result.F_world[2])
    Q = float(result.Q_spin)
    return T, Q


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--yaml", type=Path, default=DEFAULT_YAML,
                    help="rotor + polar YAML (same one fed to XROTOR)")
    ap.add_argument("--xrotor-csv", type=Path, default=DEFAULT_XCSV,
                    help="XROTOR results CSV (from verification/data/)")
    args = ap.parse_args()

    if not args.xrotor_csv.exists():
        print(f"Missing XROTOR reference: {args.xrotor_csv}", file=sys.stderr)
        print("Generate it first:", file=sys.stderr)
        print("  cd verification/xrotor_docker", file=sys.stderr)
        print(f"  docker compose run --rm xrotor --config /in/{args.yaml.name}",
              file=sys.stderr)
        return 1

    with args.yaml.open() as fh:
        cfg = yaml.safe_load(fh)
    model = build_dynbem_model(cfg)
    R = model.defn.blade.radius_m
    A = math.pi * R * R

    with args.xrotor_csv.open() as fh:
        xr_rows = {float(r["pitch_deg"]): r for r in csv.DictReader(fh)}

    print(f"Caradonna-Tung hover  --  dynbem.QuasiStaticBEM vs XROTOR vs paper")
    print(f"R = {R:.3f} m, 2 blades NACA 0012, 1250 rpm")
    print()
    print(f"{'pitch':>5}  "
          f"{'T_paper':>9} {'T_xrotor':>9} {'T_dynbem':>9}   "
          f"{'CT_paper':>9} {'CT_xrotor':>10} {'CT_dynbem':>10}   "
          f"{'dynbem/paper':>13} {'dynbem/xrotor':>14}")

    rows_out = []
    for pitch in sorted(xr_rows):
        xr_row = xr_rows[pitch]
        Omega_rpm = float(xr_row["Omega_rpm"])
        omega = Omega_rpm * math.pi / 30.0
        T_xr = float(xr_row["T_N"])
        # dynbem
        T_db, Q_db = dynbem_hover_thrust(model, pitch, Omega_rpm, RHO)
        # paper-derived thrust
        CT_paper = PAPER_CT.get(pitch, float("nan"))
        T_paper  = CT_paper * RHO * A * (omega * R) ** 2
        CT_xr    = T_xr / (RHO * A * (omega * R) ** 2)
        CT_db    = T_db / (RHO * A * (omega * R) ** 2)
        r_db_paper  = T_db / T_paper if T_paper > 0 else float("nan")
        r_db_xrotor = T_db / T_xr    if T_xr    > 0 else float("nan")
        print(f"{pitch:5.1f}  "
              f"{T_paper:9.1f} {T_xr:9.1f} {T_db:9.1f}   "
              f"{CT_paper:9.5f} {CT_xr:10.5f} {CT_db:10.5f}   "
              f"{r_db_paper:>+13.3f} {r_db_xrotor:>+14.3f}")
        rows_out.append({
            "pitch_deg":     pitch,
            "Omega_rpm":     Omega_rpm,
            "T_paper_N":     T_paper,
            "T_xrotor_N":    T_xr,
            "T_dynbem_N":    T_db,
            "CT_paper":      CT_paper,
            "CT_xrotor":     CT_xr,
            "CT_dynbem":     CT_db,
            "Q_dynbem_Nm":   Q_db,
        })

    # Name the output after the YAML stem so coarse-polar vs xfoil-polar
    # runs don't clobber each other.
    out_path = ROOT / "verification" / "data" / f"dynbem_qsbem_vs_xrotor_{args.yaml.stem}.csv"
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys()))
        w.writeheader()
        w.writerows(rows_out)
    print()
    print(f"Wrote {out_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
