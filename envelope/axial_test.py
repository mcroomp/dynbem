"""Axial, zero-gravity baseline for isolating the v_along integrator.

Tether horizontal (elevation=0 deg), gravity disabled.  With horizontal wind
the hub axis aligns with the tether, so the wind passes exactly axially
through the disk and lambda_c = lambda_s = 0 throughout.  Any blow-up
isolates to the v_along control loop rather than Pitt-Peters cyclic
dynamics or geometric tilt.

Usage
-----
    .venv\\Scripts\\python -m envelope.axial_test
    .venv\\Scripts\\python -m envelope.axial_test --v_target -2.0 --t_max 1500
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

from envelope.point_mass import ramp_column_worker


def run_axial_ramp(
    rotor_yaml: str | Path = "rotors/beaupoil_2026/rotor.yaml",
    mass_kg: float = 5.0,
    v_target: float = -0.5,
    wind_speed: float = 10.0,
    t_min: float = 100.0,
    t_max: float = 1500.0,
    ramp_rate: float = 0.5,
    sample_dn: float = 1.0,
    settle_time: float = 20.0,
    dt: float = 0.02,
    omega_init: float = 20.0,
    kp_col: float = 0.01,
    ki_col: float = 0.02,
    col_min: float = -0.25,
    col_max: float = 0.20,
) -> dict:
    """Run the axial baseline ramp and return the worker output dict.

    Uses ramp_column_worker with elevation_deg=0 and gravity_mps2=0 so the
    production code path is exercised — only the geometry/gravity are
    neutralised.  With horizontal wind, this orientation makes wind axial
    through the disk.  Any divergence reflects the integrator + control
    loop only.
    """
    from dynbem.rotor_definition import load as load_rotor

    defn = load_rotor(str(rotor_yaml))
    args = {
        "defn":          defn,
        "mass_kg":       mass_kg,
        "v_target":      v_target,
        "wind_speed":    wind_speed,
        "elevation_deg": 0.0,
        "gravity_mps2":  0.0,
        "t_min":         t_min,
        "t_max":         t_max,
        "sample_dn":     sample_dn,
        "omega_init":    omega_init,
        "settle_time":   settle_time,
        "ramp_rate":     ramp_rate,
        "dt":            dt,
        "kp_col":        kp_col,
        "ki_col":        ki_col,
        "col_min":       col_min,
        "col_max":       col_max,
        "project_root":  str(Path(__file__).parent.parent),
    }
    return ramp_column_worker(args)


def _summarise(res: dict) -> None:
    ts = res["tensions"]
    cols = np.degrees(res["cols"])
    vs = res["v_alongs"]
    om = res["omegas"]
    finite = np.isfinite(cols)
    last_t = float(ts[finite][-1]) if finite.any() else float("nan")
    print(f"Ramp covered {ts[0]:.0f} -> {last_t:.0f} N  ({len(ts)} samples)")
    print()
    print(f"{'T(N)':>8} {'v(m/s)':>9} {'omega':>8} {'col(deg)':>10}")
    step = max(1, len(ts) // 12)
    for i in range(0, len(ts), step):
        print(f"{ts[i]:8.1f} {vs[i]:+9.3f} {om[i]:8.2f} {cols[i]:+10.3f}")
    if len(ts) and (math.isfinite(cols[-1])):
        print(f"{ts[-1]:8.1f} {vs[-1]:+9.3f} {om[-1]:8.2f} {cols[-1]:+10.3f}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--v_target",   type=float, default=-0.5)
    p.add_argument("--wind",       type=float, default=10.0)
    p.add_argument("--mass",       type=float, default=5.0)
    p.add_argument("--t_min",      type=float, default=100.0)
    p.add_argument("--t_max",      type=float, default=1500.0)
    p.add_argument("--ramp_rate",  type=float, default=0.5)
    p.add_argument("--dt",         type=float, default=0.02)
    p.add_argument("--rotor",      default="rotors/beaupoil_2026/rotor.yaml")
    args = p.parse_args()

    res = run_axial_ramp(
        rotor_yaml=args.rotor,
        mass_kg=args.mass,
        v_target=args.v_target,
        wind_speed=args.wind,
        t_min=args.t_min,
        t_max=args.t_max,
        ramp_rate=args.ramp_rate,
        dt=args.dt,
    )
    _summarise(res)
