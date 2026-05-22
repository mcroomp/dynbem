"""Compare dynbem.OyeBEMModel against OpenFAST AeroDyn DBEMT_Mod=1 on
NREL Phase VI Sequence S.

Both codes implement the same algorithm: Oye's 2-stage annular
dynamic-inflow filter with constant tau1. Differences should reflect
implementation choices (mass-flow term, W_qs solver linearisation,
hover floor, tip-loss formulation), not the underlying physics.

The OpenFAST reference is produced by:

    cd verification/openfast_docker
    docker compose run --rm openfast --config /in/nrel_phase_vi_dbemt.yaml

Convention reconciliation (same as
verification/dynbem_vs_ccblade_nrel_phase_vi.py):
  - dynbem positive collective is pitch-to-stall (helicopter), the
    OpenFAST YAML uses positive pitch-to-feather (turbine).  -> negate
    pitch when handing off to dynbem.
  - Wind blows up through the disk in NED: wind_world=(0, 0, -U_wind).
  - dynbem Q < 0 in energy-extraction mode; AeroDyn's RtAeroMxh is
    negative in the same regime.  abs() both for comparison.
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dynbem import OyeBEMModel, RotorInputs                            # noqa: E402
from dynbem.rotor_definition import load as load_rotor                  # noqa: E402
from dynbem.rotor_state import OyeRotorState                            # noqa: E402

ROTOR_YAML = ROOT / "rotors" / "nrel_phase_vi" / "rotor.yaml"
OPENFAST_CSV = ROOT / "verification" / "data" / "nrel_phase_vi_dbemt_openfast.csv"
RHO = 1.225


def relax_to_steady(
    model: OyeBEMModel, inputs: RotorInputs, omega_rad_s: float,
    settle_s: float = 30.0, dt: float = 0.02,
) -> tuple[float, float]:
    """Drive OyeBEMModel forward in time until the dynamic-inflow state
    settles. Returns (T, Q) at the final converged step.

    Plain explicit Euler on W_int/W (tau1 ~ 0.5-3 s; dt = 0.02 s gives
    50 steps per shortest tau, well-resolved). omega and spin angle are
    held constant.
    """
    n_elements = model.defn.blade.n_elements
    state = OyeRotorState(
        W_int=np.zeros(n_elements),
        W=np.zeros(n_elements),
        omega_rad_s=omega_rad_s,
    )
    n_steps = int(round(settle_s / dt))
    last_T = last_Q = 0.0
    for _ in range(n_steps):
        result, deriv = model.compute_forces(inputs, state)
        # explicit Euler on the inflow states; mechanical stays fixed
        state = OyeRotorState(
            W_int=state.W_int + dt * deriv.W_int,
            W=state.W + dt * deriv.W,
            omega_rad_s=omega_rad_s,
            spin_angle_rad=state.spin_angle_rad,
        )
        last_T = -float(result.F_world[2])
        last_Q = float(result.Q_spin)
    return last_T, last_Q


def main() -> int:
    if not ROTOR_YAML.exists():
        print(f"Missing rotor YAML: {ROTOR_YAML}", file=sys.stderr)
        return 1
    if not OPENFAST_CSV.exists():
        print(f"Missing OpenFAST reference: {OPENFAST_CSV}", file=sys.stderr)
        print("Generate it first:", file=sys.stderr)
        print("  cd verification/openfast_docker", file=sys.stderr)
        print("  docker compose run --rm openfast --config /in/nrel_phase_vi_dbemt.yaml",
              file=sys.stderr)
        return 1

    model = OyeBEMModel(defn=load_rotor(str(ROTOR_YAML)))
    R = model.defn.blade.radius_m

    with OPENFAST_CSV.open() as fh:
        of_rows = list(csv.DictReader(fh))

    print(f"NREL Phase VI Sequence S  --  dynbem.OyeBEMModel vs AeroDyn DBEMT_Mod=1")
    print(f"R = {R:.3f} m, n_elements = {model.defn.blade.n_elements}, "
          f"polar = {model.defn.airfoil.name or '(unnamed)'}")
    print()
    print(f"{'U':>4}  "
          f"{'T_OF':>9} {'T_OYE':>9}  {'dT%':>6}   "
          f"{'P_OF':>9} {'P_OYE':>9}  {'dP%':>6}")

    rows_out = []
    for r in of_rows:
        U_wind  = float(r["U_wind_ms"])
        Omega   = float(r["Omega_rpm"])
        pitch   = float(r["pitch_deg"])
        T_OF    = abs(float(r["RtAeroFxh"]))
        Q_OF    = abs(float(r["RtAeroMxh"]))
        P_OF    = abs(float(r["RtAeroPwr"]))
        omega   = Omega * math.pi / 30.0

        inp = RotorInputs(
            collective_rad=math.radians(-pitch),  # turbine -> helicopter sign
            tilt_lon=0.0, tilt_lat=0.0,
            R_hub=np.eye(3),
            v_hub_world=np.zeros(3),
            wind_world=np.array([0.0, 0.0, -U_wind]),
            t=0.0,
        )
        T_oye, Q_oye = relax_to_steady(model, inp, omega)
        # dynbem returns Q_spin in *helicopter* convention -- the rotor's
        # aerodynamic drag torque.  In wind-turbine mode (energy extraction)
        # Q_spin < 0, meaning the wind drives the rotor; abs() for the
        # magnitude-comparison.
        Q_oye = abs(Q_oye)
        P_oye = Q_oye * omega

        dT = (T_oye - T_OF) / T_OF * 100 if T_OF else 0.0
        dP = (P_oye - P_OF) / P_OF * 100 if P_OF else 0.0
        print(f"{U_wind:4.0f}  "
              f"{T_OF:9.1f} {T_oye:9.1f}  {dT:5.1f}%   "
              f"{P_OF:9.1f} {P_oye:9.1f}  {dP:5.1f}%")
        rows_out.append({
            "U_wind_ms": U_wind, "Omega_rpm": Omega, "pitch_deg": pitch,
            "T_openfast_N": T_OF, "T_oye_N": T_oye,
            "Q_openfast_Nm": Q_OF, "Q_oye_Nm": Q_oye,
            "P_openfast_W": P_OF, "P_oye_W": P_oye,
        })

    # Aggregate stats
    dT_pct = np.array([(r["T_oye_N"] - r["T_openfast_N"]) / r["T_openfast_N"] * 100
                        for r in rows_out if r["T_openfast_N"] > 0])
    dP_pct = np.array([(r["P_oye_W"] - r["P_openfast_W"]) / r["P_openfast_W"] * 100
                        for r in rows_out if r["P_openfast_W"] > 0])
    print()
    print(f"|dT|  mean = {np.mean(np.abs(dT_pct)):5.2f}%   max = {np.max(np.abs(dT_pct)):5.2f}%")
    print(f"|dP|  mean = {np.mean(np.abs(dP_pct)):5.2f}%   max = {np.max(np.abs(dP_pct)):5.2f}%")

    # Persist CSV for downstream use
    out_path = ROOT / "verification" / "data" / "dynbem_oye_vs_openfast_nrel_phase_vi.csv"
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys()))
        w.writeheader()
        w.writerows(rows_out)
    print(f"Wrote {out_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
