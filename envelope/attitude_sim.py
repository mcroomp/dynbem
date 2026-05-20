"""Cyclic-pitch attitude regulator around a steady operating point.

Companion simulator to ``envelope/compute_map.py``.  Where the latter
treats the hub orientation as a quasi-static force balance (the hub Z
axis tracks the resultant of tether and gravity instantly), this one
treats hub pitch and roll as **dynamic states** driven by the rotor's
in-plane hub moments and regulated by two cyclic-pitch PIDs.

Modelling scope
---------------
Inputs (set externally; obtained from a `compute_map` sweep point):

  * elevation_deg, tension_n, wind_speed_ms  — operating point
  * collective_eq_rad, omega_init            — converged collective + ω
  * mass_kg, I_airframe_kgm2                 — vehicle mass & inertia

Equilibrium hub orientation (the "set-point" the attitude controller
holds):

    bz_eq = (T · t̂_tether + m · g · ẑ) / |…|         (same as compute_map)

with bx_eq, by_eq picked orthogonally via ``_build_r_hub``.

State variables (in addition to the rotor's inflow + ω + ψ):

    pitch, roll           — RHR rotation angles about (by_eq, bx_eq)
    pitch_rate, roll_rate — d/dt of those
    v_along               — along-tether velocity (as in compute_map)

Rotation parameterisation:

    R_hub = R_eq · R_y_local(pitch) · R_x_local(roll)

so pitch=roll=0 reproduces compute_map's force-balance R_hub_eq.  Local
rotations are intrinsic in the equilibrium frame.

Simplification — perfect yaw compensator
-----------------------------------------
We assume an external yaw mechanism that absorbs the rotor reaction
torque (``Q_spin·hub_axis``) and holds the body frame yaw-fixed.  The
airframe therefore has only two rotational DOF: pitch about ``by_eq``
and roll about ``bx_eq``.

Sign conventions (NED, right-hand rule)
---------------------------------------
  pitch > 0  ⇔  +X side of disk tilts toward +Z (down)   i.e. "nose down"
  roll  > 0  ⇔  +Y side of disk tilts toward +Z (down)   i.e. "roll right"

Plant transfer (from CLAUDE.md "Cyclic pitch convention"):

    tilt_lon > 0  →  M_y_hub < 0  →  d(pitch_rate)/dt = M_y/I < 0
    tilt_lat > 0  →  M_x_hub > 0  →  d(roll_rate)/dt  = M_x/I > 0

Restoring (stable) cyclic feedback therefore uses:

    tilt_lon = +k · pitch + k_d · pitch_rate           (positive gains)
    tilt_lat = -k · roll  - k_d · roll_rate            (negative gains)

The asymmetry comes from the opposite signs of the (tilt → moment)
transfers on the two axes.

Collective control
------------------
The along-tether velocity loop is identical to ``compute_map``'s — PI on
v_along with anti-windup, output clipped to (col_min, col_max).

Limitations
-----------
* Scalar airframe inertia (Ixx = Iyy = I_airframe_kgm2).  Real aircraft
  have a tensor with cross-coupling, but for this controller demo a
  spherical inertia is enough.
* No gyroscopic precession.  A fast-spinning rotor with hub-axis
  perturbations would actually precess (the rotor's angular momentum
  resists tilting).  Adding it requires the rotor inertia and
  rotor-axis projection; for steady-state hold it's negligible (small
  rates) — for fast manoeuvres it matters.
* No tether dynamics.  Tension is treated as a fixed scalar; in
  reality it varies with v_along.
* No off-axis cyclic-inflow harmonics under Øye.  Pitt-Peters reports
  λ_c/λ_s response to cyclic; Øye reports zero.  Both produce the
  hub moments needed for the attitude controller (which is the
  outer-loop quantity); only the transient response shape differs.

The simulator is *deliberately* light — the goal is to demonstrate
two-PID cyclic-attitude control on top of the existing rotor model,
not to be a complete flight simulator.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from dynbem import RotorInputs, create_aero, solve_trim_cyclic
from dynbem.rotor_definition import RotorDefinition
from envelope.point_mass import (
    G,
    _build_r_hub,
    _clip_state,
    _step_state_semi_implicit,
    tether_hat,
)


# ---------------------------------------------------------------------------
# Small PID with anti-windup, used identically for pitch and roll
# ---------------------------------------------------------------------------

@dataclass
class PID:
    """Standard PID with integral anti-windup.

    Signs follow the textbook form ``u = kp·err + ki·∫err + kd·d_err/dt``
    with ``err = setpoint − measurement``.  Callers handle plant-side
    sign conventions externally (e.g. by negating the output for the
    roll axis where ``tilt_lat`` has the opposite sign relationship to
    its measured state).
    """
    kp: float
    ki: float
    kd: float
    out_min: float = -math.inf
    out_max: float = math.inf
    _integral:    float = field(default=0.0, init=False)
    _last_err:    float = field(default=0.0, init=False)
    _initialised: bool  = field(default=False, init=False)

    def reset(self) -> None:
        self._integral    = 0.0
        self._last_err    = 0.0
        self._initialised = False

    def update(self, setpoint: float, measurement: float, dt: float) -> float:
        err = setpoint - measurement
        # First call: no derivative term (avoids spurious kick from initial conditions).
        if not self._initialised:
            d_err = 0.0
            self._initialised = True
        else:
            d_err = (err - self._last_err) / dt
        self._last_err = err
        # Anti-windup: only integrate if the previous output wasn't saturated.
        u_unsat = self.kp * err + self.ki * self._integral + self.kd * d_err
        if self.out_min <= u_unsat <= self.out_max:
            self._integral += err * dt
        u = self.kp * err + self.ki * self._integral + self.kd * d_err
        return max(self.out_min, min(self.out_max, u))


# ---------------------------------------------------------------------------
# Attitude simulator
# ---------------------------------------------------------------------------

def simulate_attitude(
    *,
    # ---- rotor / aero ----
    defn:                RotorDefinition,
    model:               str = "pitt_peters_jit",
    # ---- operating point (typically from a compute_map sweep cell) ----
    elevation_deg:       float,
    tension_n:           float,
    wind_speed_ms:       float,
    collective_eq_rad:   float,
    omega_init:          float,
    # ---- vehicle ----
    mass_kg:             float = 5.0,
    I_airframe_kgm2:     float = 0.5,
    # ---- attitude PID gains (positive = stable per the sign analysis) ----
    kp_att:              float = 4.0,
    ki_att:              float = 0.0,
    kd_att:              float = 1.0,
    tilt_min:            float = -math.radians(15.0),
    tilt_max:            float =  math.radians(15.0),
    # ---- collective PI (matches compute_map defaults) ----
    v_target:            Optional[float] = None,
    fix_omega:           bool = True,            # hold ω at omega_init throughout
    trim_tolerance_Nm:   float = 0.02,           # |Mx|, |My| target for Newton trim
    kp_col:              float = 0.01,
    ki_col:              float = 0.02,
    col_min:             float = -0.25,
    col_max:             float =  0.20,
    # ---- initial perturbation ----
    pitch_init:          float = 0.0,
    roll_init:           float = 0.0,
    pitch_rate_init:     float = 0.0,
    roll_rate_init:      float = 0.0,
    v_along_init:        float = 0.0,
    # ---- sim ----
    dt:                  float = 0.005,
    t_max:               float = 10.0,
    settle_time:         float = 0.0,
    gravity_mps2:        float = G,
) -> dict:
    """Simulate the closed-loop attitude response from a perturbed initial state.

    Returns a dict with arrays sampled at every timestep:

        t              [s]            simulation time
        pitch, roll    [rad]          hub-axis perturbation angles
        pitch_rate, roll_rate [rad/s] their time derivatives
        tilt_lon       [rad]          longitudinal cyclic command
        tilt_lat       [rad]          lateral cyclic command
        col            [rad]          collective command (PI on v_along)
        v_along        [m/s]          along-tether velocity
        omega          [rad/s]        rotor speed
        Mx_hub, My_hub [N·m]          rotor in-plane hub moments
        T_aero         [N]            rotor thrust magnitude

    Use the time histories to inspect convergence rate, overshoot,
    steady-state error, etc.
    """
    aero = create_aero(defn, model=model)

    # -----------------------------------------------------------------------
    # Equilibrium frame from force balance (same construction as compute_map)
    # -----------------------------------------------------------------------
    t_hat   = tether_hat(elevation_deg)
    wind    = np.array([0.0, -wind_speed_ms, 0.0])
    f_load  = tension_n * t_hat + np.array([0.0, 0.0, mass_kg * gravity_mps2])
    bz_eq   = f_load / float(np.linalg.norm(f_load))
    R_hub_eq = _build_r_hub(bz_eq)              # columns = bx_eq, by_eq, bz_eq

    f_grav_along = mass_kg * gravity_mps2 * float(t_hat[2])

    # -----------------------------------------------------------------------
    # PIDs.  Pitch loop uses *positive* feedback on the measured state
    # (tilt_lon > 0 → M_y < 0 → pitch decreases, so positive gain restores).
    # Roll loop uses *negative* feedback (tilt_lat > 0 → M_x > 0 → roll
    # increases, so the gain sign is flipped).  We implement both with the
    # standard textbook PID by sign-flipping the setpoint/measurement
    # before feeding it in.
    #
    # Net effect (PD only, ki=0): tilt_lon = +(kp·pitch + kd·pitch_rate)
    #                             tilt_lat = -(kp·roll  + kd·roll_rate)
    # -----------------------------------------------------------------------
    pid_pitch = PID(kp=kp_att, ki=ki_att, kd=kd_att,
                    out_min=tilt_min, out_max=tilt_max)
    pid_roll  = PID(kp=kp_att, ki=ki_att, kd=kd_att,
                    out_min=tilt_min, out_max=tilt_max)
    int_col   = 0.0   # collective PI integral (separate from PID class for parity with compute_map)
    int_max   = (col_max - col_min) / max(ki_col, 1e-9)

    # -----------------------------------------------------------------------
    # Rotor state init — model-agnostic
    # -----------------------------------------------------------------------
    zero    = aero.initial_rotor_state()
    arr0    = zero.to_array()
    arr0[-2] = omega_init                       # arr[-2] = ω by convention
    rotor_state = zero.from_array(arr0)

    pitch       = 0.0
    roll        = 0.0
    pitch_rate  = 0.0
    roll_rate   = 0.0
    v_along     = v_along_init
    col_now     = collective_eq_rad
    trim_tilt_lon = 0.0   # accumulated during settle, pre-loads PID integrators
    trim_tilt_lat = 0.0

    def _advance_one_step(pitch, roll, pitch_rate, roll_rate,
                          v_along, col_now, rotor_state, int_col, t_now):
        """Compute one timestep.  Returns the updated locals + telemetry."""
        cp, sp = math.cos(pitch), math.sin(pitch)
        cr, sr = math.cos(roll),  math.sin(roll)
        R_y = np.array([[ cp, 0.0,  sp],
                        [0.0, 1.0, 0.0],
                        [-sp, 0.0,  cp]])
        R_x = np.array([[1.0, 0.0, 0.0],
                        [0.0,  cr, -sr],
                        [0.0,  sr,  cr]])
        R_hub = R_hub_eq @ R_y @ R_x
        vel   = v_along * t_hat

        tilt_lon = pid_pitch.update(setpoint=0.0, measurement=-pitch, dt=dt)
        tilt_lat = pid_roll .update(setpoint=0.0, measurement=+roll,  dt=dt)

        inputs = RotorInputs(
            collective_rad=col_now,
            tilt_lon=tilt_lon, tilt_lat=tilt_lat,
            R_hub=R_hub, v_hub_world=vel, wind_world=wind, t=t_now,
        )
        aero_res, drv = aero.compute_forces(inputs, rotor_state)

        M_eq = R_hub_eq.T @ aero_res.M_orbital
        Mx_eq, My_eq = float(M_eq[0]), float(M_eq[1])

        # Attitude integration (Euler — slow modes vs dt)
        d_pitch_rate = My_eq / I_airframe_kgm2
        d_roll_rate  = Mx_eq / I_airframe_kgm2
        pitch_rate += dt * d_pitch_rate
        roll_rate  += dt * d_roll_rate
        pitch      += dt * pitch_rate
        roll       += dt * roll_rate

        # v_along along the (fixed) tether direction
        f_thrust_along = float(np.dot(aero_res.F_world, t_hat))
        f_along = f_thrust_along + f_grav_along + tension_n
        v_along = max(-30.0, min(30.0, v_along + dt * f_along / mass_kg))

        if v_target is not None:
            err = v_along - v_target
            int_col = max(-int_max, min(int_max, int_col + err * dt))
            col_raw = kp_col * err + ki_col * int_col
            col_now = max(col_min, min(col_max, col_raw))

        new_arr = _step_state_semi_implicit(aero, rotor_state, drv, dt, inputs)
        if fix_omega:
            new_arr[-2] = omega_init
        new_arr = _clip_state(new_arr, rotor_state)
        rotor_state = rotor_state.from_array(new_arr)

        return (pitch, roll, pitch_rate, roll_rate, v_along, col_now,
                rotor_state, int_col,
                tilt_lon, tilt_lat, Mx_eq, My_eq, aero_res)

    # -----------------------------------------------------------------------
    # Phase 1a: pre-settle inflow + v_along + collective at zero cyclic.
    # ω is held fixed (= omega_init) so the inflow chases a stationary
    # operating point without rotor-torque feedback dragging ω away.
    # -----------------------------------------------------------------------
    def _settle_step(rotor_state, v_along, col_now, int_col, tlon, tlat):
        vel = v_along * t_hat
        inputs = RotorInputs(
            collective_rad=col_now, tilt_lon=tlon, tilt_lat=tlat,
            R_hub=R_hub_eq, v_hub_world=vel, wind_world=wind, t=0.0,
        )
        res, drv = aero.compute_forces(inputs, rotor_state)
        f_thrust_along = float(np.dot(res.F_world, t_hat))
        f_along = f_thrust_along + f_grav_along + tension_n
        v_along = max(-30.0, min(30.0, v_along + dt * f_along / mass_kg))
        if v_target is not None:
            err = v_along - v_target
            int_col = max(-int_max, min(int_max, int_col + err * dt))
            col_now = max(col_min, min(col_max,
                                       kp_col * err + ki_col * int_col))
        new_arr = _step_state_semi_implicit(aero, rotor_state, drv, dt, inputs)
        new_arr[-2] = omega_init
        new_arr = _clip_state(new_arr, rotor_state)
        rotor_state = rotor_state.from_array(new_arr)
        return rotor_state, v_along, col_now, int_col

    n_inflow_settle = max(1, int(round(settle_time / dt)))
    for _ in range(n_inflow_settle):
        rotor_state, v_along, col_now, int_col = _settle_step(
            rotor_state, v_along, col_now, int_col, 0.0, 0.0,
        )

    # -----------------------------------------------------------------------
    # Phase 1b: trim solver.  The rotor produces non-trivial steady-state
    # hub moments in forward flight (advancing/retreating asymmetry), so
    # the open-loop attitude system is unstable from zero cyclic.
    # ``dynbem.solve_trim_cyclic`` finds the (tilt_lon, tilt_lat) that null
    # those moments at the equilibrium attitude; we then pre-load the
    # PID integrators with the trim values below.
    # -----------------------------------------------------------------------
    vel_eq = v_along * t_hat
    trim   = solve_trim_cyclic(
        aero, rotor_state,
        collective_rad=col_now,
        R_hub=R_hub_eq, v_hub_world=vel_eq, wind_world=wind,
        tilt_min=tilt_min, tilt_max=tilt_max,
        tolerance_Nm=trim_tolerance_Nm,
        dt_relax=dt, n_inflow_relax=100,
        fix_omega=True,
    )
    trim_tilt_lon = trim.tilt_lon
    trim_tilt_lat = trim.tilt_lat
    Mx_trim       = trim.Mx_residual
    My_trim       = trim.My_residual
    rotor_state   = trim.final_state

    # Pre-load PID integrators with the trim cyclic so the controller
    # starts in "trim hold" rather than ramping up from zero.  The PID
    # output formula u = kp·err + ki·∫err + kd·d_err — at err = 0 the
    # output equals ki·∫err, so ∫err = trim/ki yields output = trim.
    if pid_pitch.ki != 0.0:
        pid_pitch._integral = trim_tilt_lon / pid_pitch.ki
    if pid_roll.ki != 0.0:
        # Roll loop uses sign-flipped measurement (we feed +roll), so the
        # PID output is the *negation* of the desired tilt_lat.  We want
        # u = trim_tilt_lat at err = 0, so ∫err = trim/ki.
        pid_roll._integral  = trim_tilt_lat / pid_roll.ki

    # Apply user-requested perturbation on top of the trimmed state.
    pitch      = pitch_init
    roll       = roll_init
    pitch_rate = pitch_rate_init
    roll_rate  = roll_rate_init

    # -----------------------------------------------------------------------
    # Pre-allocate history (one row per step)
    # -----------------------------------------------------------------------
    n_steps = int(round(t_max / dt))
    h = {
        k: np.empty(n_steps + 1) for k in [
            "t", "pitch", "roll", "pitch_rate", "roll_rate",
            "tilt_lon", "tilt_lat", "col", "v_along", "omega",
            "Mx_hub", "My_hub", "T_aero",
        ]
    }

    # -----------------------------------------------------------------------
    # Phase 2: record telemetry while the controller settles the
    # perturbed state back toward equilibrium.
    # -----------------------------------------------------------------------
    for i in range(n_steps + 1):
        # Record the *current* state at the start of step i.
        h["t"][i]          = i * dt
        h["pitch"][i]      = pitch
        h["roll"][i]       = roll
        h["pitch_rate"][i] = pitch_rate
        h["roll_rate"][i]  = roll_rate
        h["col"][i]        = col_now
        h["v_along"][i]    = v_along
        h["omega"][i]      = rotor_state.omega_rad_s
        if i == n_steps:
            # No update on the last index — leave tilt/M/T at the last
            # computed values (filled below for i < n_steps).
            h["tilt_lon"][i] = h["tilt_lon"][i-1] if i > 0 else 0.0
            h["tilt_lat"][i] = h["tilt_lat"][i-1] if i > 0 else 0.0
            h["Mx_hub"][i]   = h["Mx_hub"][i-1]   if i > 0 else 0.0
            h["My_hub"][i]   = h["My_hub"][i-1]   if i > 0 else 0.0
            h["T_aero"][i]   = h["T_aero"][i-1]   if i > 0 else 0.0
            break

        (pitch, roll, pitch_rate, roll_rate, v_along, col_now,
         rotor_state, int_col,
         tilt_lon, tilt_lat, Mx_eq, My_eq, aero_res) = _advance_one_step(
            pitch, roll, pitch_rate, roll_rate, v_along, col_now,
            rotor_state, int_col, t_now=settle_time + i * dt,
        )
        h["tilt_lon"][i] = tilt_lon
        h["tilt_lat"][i] = tilt_lat
        h["Mx_hub"][i]   = Mx_eq
        h["My_hub"][i]   = My_eq
        h["T_aero"][i]   = -float(aero_res.F_world @ R_hub_eq[:, 2])

    return {
        "history":      h,
        "R_hub_eq":     R_hub_eq,
        "bz_eq":        bz_eq,
        "operating_point": {
            "elevation_deg":     elevation_deg,
            "tension_n":         tension_n,
            "wind_speed_ms":     wind_speed_ms,
            "collective_eq_rad": collective_eq_rad,
            "omega_init":        omega_init,
            "mass_kg":           mass_kg,
            "I_airframe_kgm2":   I_airframe_kgm2,
        },
        "trim": {
            "tilt_lon":  trim_tilt_lon,
            "tilt_lat":  trim_tilt_lat,
            "Mx_resid":  Mx_trim,
            "My_resid":  My_trim,
        },
        "final_state": {
            "pitch":      pitch,
            "roll":       roll,
            "pitch_rate": pitch_rate,
            "roll_rate":  roll_rate,
            "v_along":    v_along,
            "omega":      rotor_state.omega_rad_s,
        },
    }
