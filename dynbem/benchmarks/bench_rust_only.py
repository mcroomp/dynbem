"""Stable Rust-only timing.

Runs the compute_forces hot path many times in tight inner loops and
reports the MIN across 7 outer trials (Python time.perf_counter is
fine-grained but background load corrupts means; the min is the most
stable indicator of true single-call cost). Helpful when you want to
see the effect of a build-flag change without Python noise drowning it.
"""
from __future__ import annotations

import math
import time

import numpy as np

import dynbem as rs

AIR = dict(Re_design=500_000, CL0=0.0, CL_alpha_per_rad=5.7,
           CD0=0.011, alpha_stall_deg=12.0)
BLADE = dict(n_blades=2, radius_m=0.914, root_cutout_m=0.10,
             chord_m=0.058, n_elements=15)
N_PSI = 36
N_INNER = 20000
N_OUTER = 7

R_HUB = np.eye(3)
OPS = [
    ("hover",   6.0, (0,0,0),    (0,0,0),     0.0),
    ("forward", 6.0, (10,0,0),   (0,0,0),     0.0),
    ("cyclic",  6.0, (0,0,0),    (0,0,0),     math.radians(2.0)),
]


def _build(model_name):
    defn = rs.RotorDefinition(
        blade=rs.BladeGeometry(**BLADE),
        airfoil=rs.LinearPolarParameters(**AIR),
        autorotation=rs.AutorotationProperties(I_ode_kgm2=0.02),
    )
    return rs.create_aero(defn, model=model_name, n_psi_elements=N_PSI)


def _inputs(op, omega=125.0):
    _, col, vhub, wind, tlon = op
    return rs.RotorInputs(
        collective_rad=math.radians(col),
        tilt_lon=tlon, tilt_lat=0.0,
        R_hub=R_HUB, v_hub_world=np.array(vhub, dtype=float),
        wind_world=np.array(wind, dtype=float),
        omega_rad_s=omega,
    )


def _state(model):
    return model.initial_rotor_state()


def _min_time(model, inputs, state, n_inner, n_outer):
    # warm up
    for _ in range(100):
        model.compute_forces(inputs, state)
    best = float("inf")
    for _ in range(n_outer):
        t0 = time.perf_counter()
        for _ in range(n_inner):
            model.compute_forces(inputs, state)
        t = time.perf_counter() - t0
        if t < best:
            best = t
    return best / n_inner


def main():
    print(f"Rust-only timing -- min of {N_OUTER} trials, {N_INNER} calls/trial")
    print("-" * 60)
    print(f"{'model':<14} {'scenario':<10} {'us/call':>10}")
    print("-" * 60)
    for model_name in ("bem", "pitt_peters", "oye"):
        model = _build(model_name)
        for op in OPS:
            inp = _inputs(op)
            st  = _state(model)
            us  = _min_time(model, inp, st, N_INNER, N_OUTER) * 1e6
            print(f"{model_name:<14} {op[0]:<10} {us:>10.3f}")
    print("-" * 60)


if __name__ == "__main__":
    main()
