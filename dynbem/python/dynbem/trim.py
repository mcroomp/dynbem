"""Cyclic trim and inflow-relax helpers.

The Rust extension exposes `solve_trim_cyclic_py(aero, state, base_inputs, ...)`
and `relax_inflow_py(aero, state, inputs, ...)` -- both take a pre-built
RotorInputs. The legacy Python API instead accepted `collective_rad`,
`R_hub`, `v_hub_world`, `wind_world` as direct kwargs and constructed
the inputs internally. These shims accept either calling convention.
"""
from __future__ import annotations

import numpy as np

from . import RotorInputs
from ._dynbem import (
    TrimResult,
    relax_inflow_py as _relax_inflow_py,
    solve_trim_cyclic_py as _solve_trim_cyclic_py,
)

__all__ = ["TrimResult", "solve_trim_cyclic", "relax_inflow"]


def _build_inputs(
    collective_rad: float,
    tilt_lon: float,
    tilt_lat: float,
    R_hub,
    v_hub_world,
    wind_world,
    rho_kg_m3: float = 1.225,
    motor_torque_Nm: float = 0.0,
    t: float = 0.0,
) -> RotorInputs:
    return RotorInputs(
        collective_rad=collective_rad,
        tilt_lon=tilt_lon,
        tilt_lat=tilt_lat,
        R_hub=np.ascontiguousarray(R_hub, dtype=float),
        v_hub_world=np.ascontiguousarray(v_hub_world, dtype=float),
        wind_world=np.ascontiguousarray(wind_world, dtype=float),
        rho_kg_m3=rho_kg_m3,
        motor_torque_Nm=motor_torque_Nm,
        t=t,
    )


def solve_trim_cyclic(aero, state, base_inputs=None, **kwargs):
    """Solve for the cyclic tilts that drive (Mx_hub, My_hub) to a target.

    Two calling conventions:

    1. New (Rust-native): pass a pre-built ``base_inputs: RotorInputs`` plus
       the solver knobs::

           solve_trim_cyclic(aero, state, inputs, target_moment=(0,0), ...)

    2. Legacy (compat): pass the input fields directly::

           solve_trim_cyclic(aero, state,
                             collective_rad=0.15,
                             R_hub=np.eye(3),
                             v_hub_world=np.zeros(3),
                             wind_world=np.zeros(3),
                             ...)
    """
    if base_inputs is None:
        # legacy form: extract input-building kwargs, build a RotorInputs
        input_keys = {"collective_rad", "R_hub", "v_hub_world", "wind_world",
                      "rho_kg_m3", "motor_torque_Nm", "t"}
        input_kwargs = {k: kwargs.pop(k) for k in list(kwargs) if k in input_keys}
        # The legacy trim solver searches for tilt_lon/tilt_lat, so the
        # inputs RotorInputs starts at the initial tilts (default 0, 0).
        tilt_lon_init = kwargs.get("tilt_lon_init", 0.0)
        tilt_lat_init = kwargs.get("tilt_lat_init", 0.0)
        input_kwargs.setdefault("tilt_lon", tilt_lon_init)
        input_kwargs.setdefault("tilt_lat", tilt_lat_init)
        base_inputs = _build_inputs(**input_kwargs)
    return _solve_trim_cyclic_py(aero, state, base_inputs, **kwargs)


def relax_inflow(aero, state, inputs=None, **kwargs):
    """Advance state to quasi-steady inflow at fixed inputs.

    Two calling conventions (see solve_trim_cyclic docstring):

    1. New: ``relax_inflow(aero, state, inputs, n_steps=200, dt=0.005)``
    2. Legacy: ``relax_inflow(aero, state, collective_rad=..., tilt_lon=...,
                              tilt_lat=..., R_hub=..., v_hub_world=...,
                              wind_world=..., n_steps=200, dt=0.005)``
    """
    if inputs is None:
        input_keys = {"collective_rad", "tilt_lon", "tilt_lat",
                      "R_hub", "v_hub_world", "wind_world",
                      "rho_kg_m3", "motor_torque_Nm", "t"}
        input_kwargs = {k: kwargs.pop(k) for k in list(kwargs) if k in input_keys}
        inputs = _build_inputs(**input_kwargs)
    return _relax_inflow_py(aero, state, inputs, **kwargs)
