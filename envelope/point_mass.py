"""Single-point tethered-rotor equilibrium simulation (1-D along-tether).

Coordinate system: NED (X North, Y East, Z Down).

Tether convention
-----------------
  tether_hat = unit vector from hub toward anchor, in Y-Z plane:
    tether_hat = [0, cos(el), sin(el)]
  elevation_deg=0  → horizontal tether pointing East
  elevation_deg=90 → tether pointing straight down

Force on hub from tether: F_tether = tension_n * tether_hat  (toward anchor)

Wind convention (turbine/kite mode)
-------------------------------------
  Upward wind → wind_world = [0, 0, -wind_speed]  (negative Z in NED)

Velocity convention
-------------------
  v_along = dot(vel, tether_hat)
  v_along > 0 : hub moving toward anchor (reel-in)
  v_along < 0 : hub moving away from anchor (reel-out / power stroke)

Motion is 1-D: hub velocity is constrained to the tether direction.
This eliminates the need for cyclic control (which is unimplemented in the
Level-2 aero model anyway) and matches the physical intuition for the map.
"""
from __future__ import annotations

import math
from collections import deque
from typing import Optional

import numpy as np

from aero import PittPetersModel, RotorInputs
from aero.rotor_definition import RotorDefinition
from aero.rotor_state import PittPetersRotorState

G = 9.81  # m/s²


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def tether_hat(elevation_deg: float) -> np.ndarray:
    """Unit vector from hub toward anchor (NED, tether in Y-Z plane)."""
    el = math.radians(elevation_deg)
    return np.array([0.0, math.cos(el), math.sin(el)])


def balance_bz(t_hat: np.ndarray, tension_n: float, mass_kg: float) -> np.ndarray:
    """Hub Z-axis in NED that aligns thrust to oppose tether + gravity.

    At equilibrium, F_world = -T * bz, and bz must point along the resultant
    of tether tension and weight so the thrust exactly cancels both loads.
    """
    f_load = tension_n * t_hat + np.array([0.0, 0.0, mass_kg * G])
    return f_load / float(np.linalg.norm(f_load))


def _build_r_hub(body_z: np.ndarray) -> np.ndarray:
    """Rotation matrix (hub frame → NED) with hub Z-axis = body_z."""
    bz = body_z / float(np.linalg.norm(body_z))
    ref = np.array([1.0, 0.0, 0.0]) if abs(bz[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    bx = np.cross(ref, bz)
    bx /= float(np.linalg.norm(bx))
    by = np.cross(bz, bx)
    return np.column_stack([bx, by, bz])


# ---------------------------------------------------------------------------
# Single-point equilibrium simulator
# ---------------------------------------------------------------------------

def simulate_point(
    model: PittPetersModel,
    mass_kg: float,
    col_rad: float,
    elevation_deg: float,
    tension_n: float,
    wind_speed_ms: float,
    omega_init: float = 100.0,
    v_along_init: float = 0.0,
    v_target: Optional[float] = None,
    dt: float = 0.02,
    t_max: float = 60.0,
    kp_col: float = 0.01,
    ki_col: float = 0.02,
    col_min: float = -0.25,
    col_max: float = 0.20,
    conv_window_s: float = 5.0,
    omega_conv_tol: float = 2.0,
    warm_start: Optional[dict] = None,
) -> dict:
    """Simulate rotor equilibrium at one operating point.

    Parameters
    ----------
    model           PittPetersModel instance (create once, reuse across calls)
    mass_kg         vehicle mass for force balance
    col_rad         initial collective pitch (rad); ignored when warm_start provided
    elevation_deg   tether elevation above horizontal (degrees)
    tension_n       tether tension magnitude (N)
    wind_speed_ms   upward wind speed (m/s); internally converted to [0,0,-V] NED
    omega_init      initial rotor speed (rad/s)
    v_along_init    initial velocity along tether (m/s)
    v_target        if set, PI loop drives v_along toward this value
    warm_start      dict from a previous call's 'final_state' key, to continue

    Returns
    -------
    dict with keys:
        history       list of dicts sampled at ~10 Hz
        eq            equilibrium averages over last conv_window_s (if converged)
        converged     bool
        col_saturated bool
        final_state   dict suitable as warm_start for the next call
    """
    t_hat = tether_hat(elevation_deg)
    wind_world = np.array([0.0, 0.0, -wind_speed_ms])

    # Precompute fixed geometry terms
    bz = balance_bz(t_hat, tension_n, mass_kg)
    R_hub = _build_r_hub(bz)
    f_grav_along = mass_kg * G * float(t_hat[2])   # gravity component along tether

    if warm_start is not None:
        state = PittPetersRotorState(
            lambda_0=float(warm_start.get("lambda_0", 0.0)),
            lambda_c=float(warm_start.get("lambda_c", 0.0)),
            lambda_s=float(warm_start.get("lambda_s", 0.0)),
            omega_rad_s=float(warm_start.get("omega_rad_s", omega_init)),
            spin_angle_rad=float(warm_start.get("spin_angle_rad", 0.0)),
        )
        v_along = float(warm_start.get("v_along", v_along_init))
        col_now = float(warm_start.get("col", col_rad))
        int_vcol = float(warm_start.get("int_vcol", 0.0))
    else:
        state = PittPetersRotorState(omega_rad_s=omega_init)
        v_along = v_along_init
        col_now = col_rad
        int_vcol = 0.0

    n_steps = int(t_max / dt)
    conv_n = max(1, int(conv_window_s / dt))
    log_every = max(1, int(0.1 / dt))  # ~10 Hz

    v_hist = deque(maxlen=conv_n)
    om_hist = deque(maxlen=conv_n)
    history: list[dict] = []
    col_saturated = False

    int_max = (col_max - col_min) / max(ki_col, 1e-9)

    for i in range(n_steps):
        vel = v_along * t_hat

        inputs = RotorInputs(
            collective_rad=col_now,
            tilt_lon=0.0,
            tilt_lat=0.0,
            R_hub=R_hub,
            v_hub_world=vel,
            wind_world=wind_world,
            t=i * dt,
        )

        aero, dstate = model.compute_forces(inputs, state)

        # Integrate rotor state (Euler)
        arr = state.to_array() + dt * dstate.to_array()
        state = state.from_array(arr)
        if state.omega_rad_s < 1.0:
            arr2 = state.to_array()
            arr2[3] = 1.0
            state = state.from_array(arr2)

        # 1-D force integration along tether
        f_thrust_along = float(np.dot(aero.F_world, t_hat))
        f_along = f_thrust_along + f_grav_along + tension_n
        v_along += dt * f_along / mass_kg

        # Collective PI (absolute form with anti-windup)
        if v_target is not None:
            err = v_along - v_target
            int_vcol = max(-int_max, min(int_max, int_vcol + err * dt))
            col_raw = kp_col * err + ki_col * int_vcol
            col_saturated = col_raw < col_min or col_raw > col_max
            col_now = max(col_min, min(col_max, col_raw))

        v_hist.append(v_along)
        om_hist.append(state.omega_rad_s)

        if i % log_every == 0:
            history.append({
                "t": i * dt,
                "v_along": v_along,
                "omega": state.omega_rad_s,
                "col": col_now,
                "lambda_0": state.lambda_0,
            })

    converged = False
    eq: dict = {}
    if len(v_hist) >= conv_n:
        v_range = max(v_hist) - min(v_hist)
        om_range = max(om_hist) - min(om_hist)
        converged = v_range < 0.05 and om_range < omega_conv_tol
        if converged:
            eq = {
                "v_along": sum(v_hist) / len(v_hist),
                "omega": sum(om_hist) / len(om_hist),
                "col": col_now,
            }

    return {
        "history": history,
        "eq": eq,
        "converged": converged,
        "col_saturated": col_saturated,
        "final_state": {
            "lambda_0": state.lambda_0,
            "lambda_c": state.lambda_c,
            "lambda_s": state.lambda_s,
            "omega_rad_s": state.omega_rad_s,
            "spin_angle_rad": state.spin_angle_rad,
            "v_along": v_along,
            "col": col_now,
            "int_vcol": int_vcol,
        },
    }


# ---------------------------------------------------------------------------
# Tension-ramp worker  (module-level so ProcessPoolExecutor can pickle it)
# ---------------------------------------------------------------------------

def ramp_column_worker(args: dict) -> dict:
    """Ramp tether tension continuously for one (v_target, wind, elevation).

    Runs as a separate process. Receives a plain dict (fully picklable).

    Settles at t_min, then ramps to t_max at ramp_rate N/s, recording the
    PI state at every sample_dn N increment along the way.

    args keys
    ---------
    defn            RotorDefinition (frozen dataclass, picklable)
    mass_kg         float
    v_target        float  (m/s, negative = reel-out)
    wind_speed      float  (m/s)
    elevation_deg   float  (degrees)
    t_min           float  (N, start tension)
    t_max           float  (N, end tension)
    sample_dn       float  (N between recorded samples, default 1.0)
    omega_init      float  (rad/s, default 20)
    settle_time     float  (s at t_min before ramp, default 20)
    ramp_rate       float  (N/s, default 0.5)
    dt              float  (integration timestep, default 0.02)
    kp_col, ki_col  floats (PI gains)
    col_min, col_max floats (collective limits rad)

    Returns
    -------
    dict with:
        elevation_deg, wind_speed, v_target  — echoed for assembly
        tensions   np.ndarray  sampled tension values (N)
        cols       np.ndarray  collective at each sample (rad)
        v_alongs   np.ndarray  v_along at each sample (m/s)
        sats       np.ndarray  bool, collective saturated at sample
    """
    import sys as _sys
    import os as _os

    # Ensure the aero package is importable in the spawned process
    project_root = args.get("project_root", "")
    if project_root and project_root not in _sys.path:
        _sys.path.insert(0, project_root)

    from aero import PittPetersModel
    from aero.rotor_state import PittPetersRotorState

    defn = args["defn"]
    model = PittPetersModel(defn=defn)

    mass_kg = float(args["mass_kg"])
    v_target = float(args["v_target"])
    wind_speed = float(args["wind_speed"])
    elevation_deg = float(args["elevation_deg"])
    T_min = float(args.get("t_min", 100.0))
    T_max = float(args.get("t_max", 1000.0))
    sample_dn = float(args.get("sample_dn", 1.0))
    omega_init = float(args.get("omega_init", 20.0))
    settle_time = float(args.get("settle_time", 20.0))
    ramp_rate = float(args.get("ramp_rate", 0.5))
    dt = float(args.get("dt", 0.02))
    kp_col = float(args.get("kp_col", 0.01))
    ki_col = float(args.get("ki_col", 0.02))
    col_min = float(args.get("col_min", -0.25))
    col_max = float(args.get("col_max", 0.20))

    t_hat_arr = tether_hat(elevation_deg)
    wind_world = np.array([0.0, 0.0, -wind_speed])
    f_grav_along = mass_kg * G * float(t_hat_arr[2])
    int_max = (col_max - col_min) / max(ki_col, 1e-9)

    state = PittPetersRotorState(omega_rad_s=omega_init)
    v_along = 0.0
    col_now = 0.0
    int_vcol = 0.0
    col_saturated = False
    aero_clamped = False

    _clip_lo = np.array([-10, -10, -10, 0.5, -np.inf])
    _clip_hi = np.array([ 10,  10,  10, 300,  np.inf])

    def _step(tension_now: float, sim_t: float) -> None:
        nonlocal state, v_along, col_now, int_vcol, col_saturated, aero_clamped

        bz = balance_bz(t_hat_arr, tension_now, mass_kg)
        R_hub = _build_r_hub(bz)
        vel = v_along * t_hat_arr

        inputs = RotorInputs(
            collective_rad=col_now,
            tilt_lon=0.0,
            tilt_lat=0.0,
            R_hub=R_hub,
            v_hub_world=vel,
            wind_world=wind_world,
            t=sim_t,
        )
        aero, dstate = model.compute_forces(inputs, state)

        arr_raw = state.to_array() + dt * dstate.to_array()
        arr_clipped = np.clip(arr_raw, _clip_lo, _clip_hi)
        if not np.array_equal(arr_raw, arr_clipped):
            aero_clamped = True
        state = state.from_array(arr_clipped)

        f_thrust_along = float(np.dot(aero.F_world, t_hat_arr))
        f_along = f_thrust_along + f_grav_along + tension_now
        v_new = v_along + dt * f_along / mass_kg
        if v_new < -30.0 or v_new > 30.0:
            aero_clamped = True
        v_along = max(-30.0, min(30.0, v_new))

        err = v_along - v_target
        int_vcol = max(-int_max, min(int_max, int_vcol + err * dt))
        col_raw = kp_col * err + ki_col * int_vcol
        col_saturated = col_raw < col_min or col_raw > col_max
        col_now = max(col_min, min(col_max, col_raw))

    # Phase 1: settle at T_min
    settle_steps = int(settle_time / dt)
    for i in range(settle_steps):
        _step(T_min, i * dt)

    def _tilt_deg(tension_now: float) -> float:
        """Rotor tilt from vertical (deg): angle between hub axis and -Z (down)."""
        bz = balance_bz(t_hat_arr, tension_now, mass_kg)
        return float(math.degrees(math.acos(max(-1.0, min(1.0, bz[2])))))

    # Record settled state at T_min
    t_samples = [T_min]
    col_samples = [col_now]
    v_samples = [v_along]
    sat_samples = [col_saturated]
    omega_samples = [state.omega_rad_s]
    tilt_samples = [_tilt_deg(T_min)]
    next_sample_t = T_min + sample_dn

    # Phase 2: ramp from T_min to T_max, recording at every sample_dn N
    ramp_steps = int((T_max - T_min) / ramp_rate / dt) + 1
    for i in range(ramp_steps):
        t_now = min(T_min + ramp_rate * (i + 1) * dt, T_max)
        _step(t_now, settle_time + i * dt)

        if aero_clamped:
            break

        while next_sample_t <= t_now + 1e-9 and next_sample_t <= T_max + 1e-9:
            t_samples.append(next_sample_t)
            col_samples.append(col_now)
            v_samples.append(v_along)
            sat_samples.append(col_saturated)
            omega_samples.append(state.omega_rad_s)
            tilt_samples.append(_tilt_deg(next_sample_t))
            next_sample_t += sample_dn

        if t_now >= T_max - 1e-9:
            break

    return {
        "elevation_deg": elevation_deg,
        "wind_speed": wind_speed,
        "v_target": v_target,
        "tensions": np.array(t_samples),
        "cols": np.array(col_samples),
        "v_alongs": np.array(v_samples),
        "sats": np.array(sat_samples, dtype=bool),
        "omegas": np.array(omega_samples),
        "tilts":  np.array(tilt_samples),
    }
