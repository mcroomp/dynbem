"""Steady-state cyclic trim solver.

Given an operating point (collective, R_hub, hub velocity, wind, rotor
state), find the (tilt_lon, tilt_lat) cyclic inputs that drive the
in-plane hub moments (M_x, M_y in hub frame) to a user-specified target
(zero by default).

Why this exists
---------------
In forward flight the rotor produces a non-trivial steady-state hub
moment from advancing/retreating blade asymmetry (and from cyclic-inflow
harmonics in models that resolve them).  An attitude controller with
only a proportional gain on attitude error therefore has a steady-state
attitude offset proportional to that disturbance.

The trim solver finds the cyclic offsets that cancel the disturbance so
the outer controller can be pre-loaded with the trim values and then
operate as a small-perturbation regulator around equilibrium — the
textbook helicopter trim/regulator architecture.

API
---
``solve_trim_cyclic`` returns a :class:`TrimResult` containing the trim
cyclic, the final residual moments, the number of Newton iterations, a
convergence flag, and **the final RotorState at the trim point** so the
caller can start a simulation in the trimmed equilibrium without
re-settling.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from . import AeroBase, RotorInputs
from .rotor_state import RotorState


@dataclass
class TrimResult:
    """Outcome of :func:`solve_trim_cyclic`.

    Fields
    ------
    tilt_lon         Converged longitudinal cyclic tilt [rad].
    tilt_lat         Converged lateral cyclic tilt [rad].
    Mx_residual      Hub-frame X moment at the trim, minus the target [N·m].
                     Should satisfy |Mx_residual| < tolerance on success.
    My_residual      Hub-frame Y moment at the trim, minus the target [N·m].
    iterations       Number of Newton iterations performed.
    converged        True iff both residuals < ``tolerance_Nm`` at exit.
    final_state      RotorState at the trim cyclic, relaxed to its
                     quasi-steady inflow.  Use this to start a simulation
                     in the trimmed operating point.
    """
    tilt_lon:    float
    tilt_lat:    float
    Mx_residual: float
    My_residual: float
    iterations:  int
    converged:   bool
    final_state: RotorState


def _semi_implicit_step(
    aero:          AeroBase,
    state:         RotorState,
    derivative:    RotorState,
    inputs:        RotorInputs,
    dt:            float,
    fix_omega_to:  Optional[float],
) -> RotorState:
    """One semi-implicit Euler step: damp dynamic-inflow states by
    ``1/(1 + dt/τ)`` from ``dynbem.inflow_taus``; plain Euler on mechanical
    states (τ = ∞).  Unconditionally stable for stiff inflow ODEs.

    When ``fix_omega_to`` is not None, ω is held at that value
    (convention: ``arr[-2] = omega_rad_s``).
    """
    taus  = aero.inflow_taus(inputs, state)
    arr   = state.to_array()
    darr  = derivative.to_array()
    damp  = np.ones_like(arr)
    finite = np.isfinite(taus)
    damp[finite] = 1.0 / (1.0 + dt / taus[finite])
    new_arr = arr + dt * darr * damp
    if fix_omega_to is not None:
        new_arr[-2] = float(fix_omega_to)
    return state.from_array(new_arr)


def relax_inflow(
    aero:          AeroBase,
    state:         RotorState,
    *,
    collective_rad: float,
    tilt_lon:       float,
    tilt_lat:       float,
    R_hub:          np.ndarray,
    v_hub_world:    np.ndarray,
    wind_world:     np.ndarray,
    n_steps:        int   = 200,
    dt:             float = 0.005,
    fix_omega:      bool  = True,
    t:              float = 0.0,
    rho_kg_m3:      float = 1.225,
) -> RotorState:
    """Advance the rotor state to a quasi-steady inflow at fixed inputs.

    Useful both as a building block for :func:`solve_trim_cyclic` and as
    a standalone helper for callers that just want a settled state at
    some operating point.

    When ``fix_omega`` is True, ω is clamped to its initial value (the
    state.omega_rad_s at entry) so the inflow chases a stationary
    operating point without rotor-torque feedback dragging ω away.

    Returns the relaxed state.
    """
    omega0 = float(state.omega_rad_s) if fix_omega else None
    inputs = RotorInputs(
        collective_rad=collective_rad,
        tilt_lon=tilt_lon, tilt_lat=tilt_lat,
        R_hub=R_hub, v_hub_world=v_hub_world, wind_world=wind_world,
        t=t, rho_kg_m3=rho_kg_m3,
    )
    for _ in range(n_steps):
        _, deriv = aero.compute_forces(inputs, state)
        state = _semi_implicit_step(aero, state, deriv, inputs, dt, omega0)
    return state


def solve_trim_cyclic(
    aero:           AeroBase,
    rotor_state:    RotorState,
    *,
    collective_rad: float,
    R_hub:          np.ndarray,
    v_hub_world:    np.ndarray,
    wind_world:     np.ndarray,
    target_moment:  Tuple[float, float] = (0.0, 0.0),
    tilt_lon_init:  float = 0.0,
    tilt_lat_init:  float = 0.0,
    tilt_min:       float = -math.radians(15.0),
    tilt_max:       float =  math.radians(15.0),
    tolerance_Nm:   float = 0.02,
    max_iterations: int   = 50,
    probe_rad:      float = math.radians(0.5),
    dt_relax:       float = 0.005,
    n_inflow_relax: int   = 100,
    n_settle:       int   = 0,
    fix_omega:      bool  = True,
    t:              float = 0.0,
    rho_kg_m3:      float = 1.225,
) -> TrimResult:
    """Find (tilt_lon, tilt_lat) that drive (M_x_hub, M_y_hub) to ``target_moment``.

    Hub-frame moments are computed as ``R_hub.T @ AeroResult.M_orbital``.
    Cyclic inputs are clipped to ``[tilt_min, tilt_max]``.

    Algorithm
    ---------
    1. (optional) Settle inflow at the initial cyclic for ``n_settle``
       semi-implicit-Euler steps.
    2. Repeat until convergence or ``max_iterations``:
       a. Read (Mx, My) at the current trim.
       b. Probe ∂My/∂tilt_lon and ∂Mx/∂tilt_lat by a finite-difference
          step of ``probe_rad`` and reading the residuals.
       c. Damped Newton (half-step) update on each axis, clipped to
          bounds.
       d. Relax the inflow for ``n_inflow_relax`` steps at the new
          cyclic so the next probe reflects the new operating point.

    Parameters
    ----------
    aero            Aero model (any subclass of AeroBase).
    rotor_state     Initial rotor state.  Not mutated — the final state
                    is returned in ``TrimResult.final_state``.
    collective_rad  Held fixed during trim.
    R_hub           Hub orientation (body→world).
    v_hub_world     Hub velocity in world frame [m/s].
    wind_world      Wind velocity in world frame [m/s].
    target_moment   Target hub-frame ``(M_x, M_y)`` in N·m.  Default
                    ``(0, 0)`` solves the conventional "null the moments"
                    trim; non-zero targets can model an offset CG or
                    intentional disk-plane offset.
    tilt_lon_init,
    tilt_lat_init   Initial cyclic guess.
    tilt_min,
    tilt_max        Hard limits on cyclic [rad].
    tolerance_Nm    |Mx − target_x|, |My − target_y| target [N·m].
    max_iterations  Newton iteration cap.
    probe_rad       Finite-difference probe offset [rad].
    dt_relax        Timestep used to relax inflow.
    n_inflow_relax  Inflow-relax steps after each cyclic update.
    n_settle        Optional pre-settle steps at the initial cyclic
                    (use when the supplied state is far from settled).
    fix_omega       Hold ω at its initial value throughout (recommended
                    for ω-regulated regimes).
    t, rho_kg_m3    Passed through to RotorInputs.

    Returns
    -------
    :class:`TrimResult`.  Even on non-convergence the returned cyclic
    is the best-effort trim and ``final_state`` is the relaxed inflow
    at that cyclic.
    """
    target_x, target_y = float(target_moment[0]), float(target_moment[1])
    omega0 = float(rotor_state.omega_rad_s)
    fix_to = omega0 if fix_omega else None

    tilt_lon = float(np.clip(tilt_lon_init, tilt_min, tilt_max))
    tilt_lat = float(np.clip(tilt_lat_init, tilt_min, tilt_max))
    state    = rotor_state

    def _eval(s: RotorState, tlon: float, tlat: float):
        """Single compute_forces call; returns (Mx_hub − target_x,
        My_hub − target_y, derivative, inputs)."""
        inputs = RotorInputs(
            collective_rad=collective_rad,
            tilt_lon=tlon, tilt_lat=tlat,
            R_hub=R_hub, v_hub_world=v_hub_world, wind_world=wind_world,
            t=t, rho_kg_m3=rho_kg_m3,
        )
        result, deriv = aero.compute_forces(inputs, s)
        M_hub = R_hub.T @ result.M_orbital
        return (float(M_hub[0]) - target_x,
                float(M_hub[1]) - target_y,
                deriv, inputs)

    def _relax(s: RotorState, tlon: float, tlat: float,
               n: int) -> RotorState:
        for _ in range(n):
            _, _, deriv, inputs = _eval(s, tlon, tlat)
            s = _semi_implicit_step(aero, s, deriv, inputs, dt_relax, fix_to)
        return s

    # Optional pre-settle at the initial cyclic.
    if n_settle > 0:
        state = _relax(state, tilt_lon, tilt_lat, n_settle)

    # Settle at the starting trim guess so the first moment read is steady.
    state = _relax(state, tilt_lon, tilt_lat, n_inflow_relax)
    Mx, My, _, _ = _eval(state, tilt_lon, tilt_lat)
    converged    = abs(Mx) < tolerance_Nm and abs(My) < tolerance_Nm
    iterations   = 0

    for iterations in range(1, max_iterations + 1):
        if converged:
            break
        # Probe ∂My/∂tilt_lon
        _, My_probe, _, _ = _eval(state, tilt_lon + probe_rad, tilt_lat)
        dMy_dlon = (My_probe - My) / probe_rad
        # Probe ∂Mx/∂tilt_lat
        Mx_probe, _, _, _ = _eval(state, tilt_lon, tilt_lat + probe_rad)
        dMx_dlat = (Mx_probe - Mx) / probe_rad

        if abs(dMy_dlon) > 1e-6:
            tilt_lon = float(np.clip(
                tilt_lon - 0.5 * My / dMy_dlon, tilt_min, tilt_max,
            ))
        if abs(dMx_dlat) > 1e-6:
            tilt_lat = float(np.clip(
                tilt_lat - 0.5 * Mx / dMx_dlat, tilt_min, tilt_max,
            ))

        state = _relax(state, tilt_lon, tilt_lat, n_inflow_relax)
        Mx, My, _, _ = _eval(state, tilt_lon, tilt_lat)
        converged = abs(Mx) < tolerance_Nm and abs(My) < tolerance_Nm

    return TrimResult(
        tilt_lon=tilt_lon, tilt_lat=tilt_lat,
        Mx_residual=Mx, My_residual=My,
        iterations=iterations, converged=converged,
        final_state=state,
    )
