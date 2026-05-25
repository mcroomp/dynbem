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
  Wind blows east-to-west (anchor side toward kite side):
    wind_world = [0, -wind_speed, 0]  (toward -Y in NED)
  Tether geometry already places the anchor east of the hub
  (tether_hat = [0, cos(el), sin(el)] with positive Y), so the kite is
  blown downwind (westward) into the tether and supported by the wind.

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

from dynbem import RotorInputs, create_aero
from dynbem.mechanical import omega_derivative
from dynbem.rotor_definition import RotorDefinition

G = 9.81  # m/s²


# ---------------------------------------------------------------------------
# Semi-implicit Euler for dynamic-inflow states
# ---------------------------------------------------------------------------
#
# The Pitt-Peters / Øye inflow states follow  dλ/dt = (λ_ss − λ)/τ  with
# small τ for small rotors at high V_T.  Plain explicit Euler is unstable
# when dt/τ > 2; we damp the step by (1 + dt/τ)⁻¹ — the implicit-Euler
# limit for fixed λ_ss.  Unconditionally stable; reduces to explicit
# Euler when dt ≪ τ.  τ comes from the aero model's inflow_taus()
# method (returns ∞ for mechanical states, which then integrate as
# plain Euler).

def _clip_state(arr: np.ndarray, state=None) -> np.ndarray:
    """Clip inflow states to +/-10.  Omega is tracked outside the state."""
    return np.clip(arr, -10.0, 10.0)


def _step_state_semi_implicit(model, state, dstate, dt: float,
                              inputs) -> np.ndarray:
    """Semi-implicit Euler on dynamic-inflow states; explicit on the rest.

    Queries ``model.inflow_taus(inputs, state)`` for per-state time
    constants.  Mechanical states (ω, ψ) report τ=∞ and integrate as
    plain explicit Euler.  Model-agnostic — works with any AeroBase
    subclass that overrides inflow_taus correctly.
    """
    taus = np.asarray(model.inflow_taus(inputs, state), dtype=float)
    arr = state.to_array()
    darr = dstate.to_array()

    if taus.size != arr.size:
        raise ValueError(
            f"inflow_taus size mismatch: got {taus.size}, expected {arr.size}"
        )

    # damp = 1/(1+dt/τ); =1 when τ=∞ (plain Euler)
    finite = np.isfinite(taus)
    damp = np.ones_like(arr)
    damp[finite] = 1.0 / (1.0 + dt / taus[finite])
    return arr + dt * darr * damp


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
    model,
    mass_kg: float,
    col_rad: float,
    elevation_deg: float,
    tension_n: float,
    wind_speed_ms: float,
    omega_init: float = 100.0,
    v_along_init: float = 0.0,
    v_target: Optional[float] = None,
    dt: float = 0.005,
    t_max: float = 60.0,
    kp_col: float = 0.01,
    ki_col: float = 0.02,
    col_min: float = -0.25,
    col_max: float = 0.20,
    conv_window_s: float = 5.0,
    omega_conv_tol: float = 2.0,
    warm_start: Optional[dict] = None,
    gravity_mps2: float = G,
) -> dict:
    """Simulate rotor equilibrium at one operating point.

    Parameters
    ----------
    model           AeroBase instance (create once, reuse across calls) —
                    any of PittPetersModel, PittPetersModelJIT, OyeBEMModel, etc.
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
    wind_world = np.array([0.0, -wind_speed_ms, 0.0])

    # Precompute fixed geometry terms
    f_load = tension_n * t_hat + np.array([0.0, 0.0, mass_kg * gravity_mps2])
    bz = f_load / float(np.linalg.norm(f_load))
    R_hub = _build_r_hub(bz)
    f_grav_along = mass_kg * gravity_mps2 * float(t_hat[2])   # gravity component along tether

    # Generic state initialisation — works with any AeroBase subclass.
    _autorot = getattr(model.defn, "autorotation", None)
    I_ode = float((_autorot.I_ode_kgm2 if _autorot is not None else None) or 1.0)
    if warm_start is not None and "state_arr" in warm_start:
        state = model.initial_rotor_state().from_array(np.asarray(warm_start["state_arr"]))
        omega = float(warm_start.get("omega", omega_init))
        v_along = float(warm_start.get("v_along", v_along_init))
        col_now = float(warm_start.get("col", col_rad))
        int_vcol = float(warm_start.get("int_vcol", 0.0))
    else:
        state = model.initial_rotor_state()
        omega = omega_init
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
            omega_rad_s=omega,
            rho_kg_m3=1.225,
        )

        aero, dstate = model.compute_forces(inputs, state)
        arr = _step_state_semi_implicit(model, state, dstate, dt, inputs)

        # Clip dynamic-inflow states (guards against divergence at very
        # short time constants or model start-up transients).
        state = state.from_array(_clip_state(arr))
        omega = max(0.5, min(300.0,
                             omega + dt * omega_derivative(aero.Q_spin, 0.0, I_ode)))

        # 1-D force integration along tether (explicit Euler).  The previous
        # semi-implicit form used a finite-difference probe of ∂F/∂v_along
        # to damp stiff thrust; that's no longer needed at dt = 0.005 s
        # (typical |k_v| ≪ 2m/dt = 2000 N·s/m stability bound) and the probe
        # doubled the BEM cost per step.
        f_thrust_along = float(np.dot(aero.F_world, t_hat))
        f_along = f_thrust_along + f_grav_along + tension_n
        v_along += (dt / mass_kg) * f_along

        # Collective PI (absolute form with anti-windup)
        if v_target is not None:
            err = v_along - v_target
            int_vcol = max(-int_max, min(int_max, int_vcol + err * dt))
            col_raw = kp_col * err + ki_col * int_vcol
            col_saturated = col_raw < col_min or col_raw > col_max
            col_now = max(col_min, min(col_max, col_raw))

        v_hist.append(v_along)
        om_hist.append(omega)

        if i % log_every == 0:
            history.append({
                "t": i * dt,
                "v_along": v_along,
                "omega": omega,
                "col": col_now,
                "lambda_0": getattr(state, "lambda_0", 0.0),
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
            # Model-agnostic state serialisation: pass state_arr + omega back
            # to a later simulate_point() call via warm_start to resume.
            "state_arr": state.to_array(),
            "omega":     omega,
            "v_along":   v_along,
            "col":       col_now,
            "int_vcol":  int_vcol,
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

    # Ensure the aero package is importable in the spawned process
    project_root = args.get("project_root", "")
    if project_root and project_root not in _sys.path:
        _sys.path.insert(0, project_root)

    from dynbem import create_aero

    defn        = args["defn"]
    model_name  = args.get("model", "pitt_peters_jit")
    model       = create_aero(defn, model=model_name)

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
    gravity_mps2 = float(args.get("gravity_mps2", G))

    t_hat_arr = tether_hat(elevation_deg)
    wind_world = np.array([0.0, -wind_speed, 0.0])
    f_grav_along = mass_kg * gravity_mps2 * float(t_hat_arr[2])
    int_max = (col_max - col_min) / max(ki_col, 1e-9)

    I_ode = float(defn.autorotation.I_ode_kgm2 or 1.0)

    # Model-agnostic initial state.
    state = model.initial_rotor_state()
    omega = omega_init
    v_along = 0.0
    col_now = 0.0
    int_vcol = 0.0
    col_saturated = False
    aero_clamped = False

    def _step(tension_now: float, sim_t: float, dt_in: float) -> None:
        """Single explicit step of (state, v_along, PI) at the given dt."""
        nonlocal state, omega, v_along, col_now, int_vcol, col_saturated, aero_clamped

        f_load = tension_now * t_hat_arr + np.array([0.0, 0.0, mass_kg * gravity_mps2])
        bz = f_load / float(np.linalg.norm(f_load))
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
            omega_rad_s=omega,
            rho_kg_m3=1.225,
        )
        aero, dstate = model.compute_forces(inputs, state)

        # Semi-implicit Euler step + clip.
        arr_raw = _step_state_semi_implicit(model, state, dstate, dt_in, inputs)
        arr_clipped = _clip_state(arr_raw)
        if not np.array_equal(arr_raw, arr_clipped):
            aero_clamped = True
        state = state.from_array(arr_clipped)
        omega_new = omega + dt_in * omega_derivative(aero.Q_spin, 0.0, I_ode)
        omega_clipped = max(0.5, min(300.0, omega_new))
        if omega_clipped != omega_new:
            aero_clamped = True
        omega = omega_clipped

        # Explicit Euler on v_along.
        f_thrust_along = float(np.dot(aero.F_world, t_hat_arr))
        f_along = f_thrust_along + f_grav_along + tension_now
        v_new = v_along + (dt_in / mass_kg) * f_along
        if v_new < -30.0 or v_new > 30.0:
            aero_clamped = True
        v_along = max(-30.0, min(30.0, v_new))

        err = v_along - v_target
        int_vcol = max(-int_max, min(int_max, int_vcol + err * dt_in))
        col_raw = kp_col * err + ki_col * int_vcol
        col_saturated = col_raw < col_min or col_raw > col_max
        col_now = max(col_min, min(col_max, col_raw))

    # Phase 1: settle at T_min
    settle_steps = int(settle_time / dt)
    for i in range(settle_steps):
        _step(T_min, i * dt, dt)

    def _tilt_deg(tension_now: float) -> float:
        """Rotor tilt from vertical (deg): angle between hub axis and -Z (down)."""
        bz = balance_bz(t_hat_arr, tension_now, mass_kg)
        return float(math.degrees(math.acos(max(-1.0, min(1.0, bz[2])))))

    # λ_c / λ_s are Pitt-Peters-specific (global harmonics).  Other
    # models (e.g. Øye, which has per-annulus inflow only) report 0 —
    # the heatmap's cyclic widget then just shows a centred dot.
    def _lam_c(st):
        return float(getattr(st, "lambda_c", 0.0))

    def _lam_s(st):
        return float(getattr(st, "lambda_s", 0.0))

    # Record settled state at T_min
    t_samples = [T_min]
    col_samples = [col_now]
    v_samples = [v_along]
    sat_samples = [col_saturated]
    omega_samples = [omega]
    tilt_samples = [_tilt_deg(T_min)]
    lamc_samples = [_lam_c(state)]
    lams_samples = [_lam_s(state)]

    # ----- Phase 2: ramp from T_min to T_max with backtracking -----
    #
    # We walk the tension axis in "segments" of `backtrack_dn` N each.  At
    # the start of every segment we snapshot the full state.  We step
    # through the segment at the current dt; when a segment finishes
    # cleanly (no aero clamp), we keep the samples and advance.  When it
    # ends in a clamp, we restore the segment-start snapshot, halve the
    # local dt, and replay the same segment — up to a recursion cap.  Once
    # past the troublesome segment, dt is restored to the user-set value.
    #
    # This keeps the per-step code simple (one dt, no adaptive bookkeeping)
    # and only pays the smaller-dt cost in the few segments that need it.
    backtrack_dn = 10.0 * sample_dn        # ~10 samples between snapshots
    min_dt = dt / 64.0                     # cost cap
    seg_t_start = T_min
    seg_sim_t_start = settle_time
    local_dt = dt

    def _snapshot():
        return {
            "state":          state,                  # RotorState is a frozen dataclass
            "v_along":        v_along,
            "col_now":        col_now,
            "int_vcol":       int_vcol,
            "col_saturated":  col_saturated,
            # samples taken so far — restored to this length on backtrack
            "n_samples":      len(t_samples),
        }

    def _restore(snap):
        nonlocal state, v_along, col_now, int_vcol, col_saturated
        state         = snap["state"]
        v_along       = snap["v_along"]
        col_now       = snap["col_now"]
        int_vcol      = snap["int_vcol"]
        col_saturated = snap["col_saturated"]
        n = snap["n_samples"]
        del t_samples [n:]; del col_samples[n:]; del v_samples [n:]
        del sat_samples[n:]; del omega_samples[n:]; del tilt_samples[n:]
        del lamc_samples[n:]; del lams_samples[n:]

    snap = _snapshot()
    next_sample_t = T_min + sample_dn
    sim_t = settle_time

    while seg_t_start < T_max - 1e-9:
        seg_t_end = min(seg_t_start + backtrack_dn, T_max)
        # Step from seg_t_start to seg_t_end at local_dt, recording samples.
        ramp_step_dt = ramp_rate * local_dt   # tension increment per step
        n_steps_in_seg = max(1, int(math.ceil((seg_t_end - seg_t_start) / ramp_step_dt)))
        t_now = seg_t_start
        for k in range(n_steps_in_seg):
            t_now = min(seg_t_start + ramp_rate * (k + 1) * local_dt, seg_t_end)
            _step(t_now, sim_t, local_dt)
            sim_t += local_dt

            if aero_clamped:
                break

            while (next_sample_t <= t_now + 1e-9
                   and next_sample_t <= T_max + 1e-9):
                t_samples.append(next_sample_t)
                col_samples.append(col_now)
                v_samples.append(v_along)
                sat_samples.append(col_saturated)
                omega_samples.append(omega)
                tilt_samples.append(_tilt_deg(next_sample_t))
                lamc_samples.append(_lam_c(state))
                lams_samples.append(_lam_s(state))
                next_sample_t += sample_dn

        if not aero_clamped:
            # Segment succeeded — advance, snapshot, try to recover dt.
            seg_t_start = seg_t_end
            seg_sim_t_start = sim_t
            local_dt = min(dt, local_dt * 2.0)   # gentle dt recovery
            snap = _snapshot()
            continue

        # Aero clamped during segment — restore and retry at smaller dt.
        if local_dt <= min_dt:
            # Already at the floor; nothing more to do — record clamp and stop.
            break
        _restore(snap)
        aero_clamped = False
        sim_t = seg_sim_t_start
        local_dt *= 0.5

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
        "lambda_c": np.array(lamc_samples),
        "lambda_s": np.array(lams_samples),
    }
