# dynbem API reference

All public names are importable directly from `dynbem`.

```python
import dynbem
from dynbem import (
    create_aero, RotorInputs, AeroResult,
    solve_trim_cyclic, relax_inflow, TrimResult,
    omega_derivative, euler_step_omega,
    QuasiStaticRotorState, PittPetersRotorState, OyeRotorState,
    LinearPolar, TabulatedPolar,
)
```

---

## Factory

### `create_aero(defn, model, *, n_psi_elements=36, polar=None)`

Build an aero model from a `RotorDefinition`.

| `model` string | Type returned | Dynamic inflow |
|---|---|---|
| `"quasi_static"` / `"bem"` | `QuasiStaticBEM` | none (steady BEM) |
| `"pitt_peters"` / `"pitt_peters_jit"` | `PittPetersModel` | Pitt-Peters 3-state |
| `"oye"` / `"oye_bem"` | `OyeBEMModel` | Oye 2-stage annular |

`n_psi_elements` controls the number of azimuth stations in the
forward-flight psi-loop (default 36; ignored in axial flight).
`polar` overrides auto-construction from `defn.airfoil`.

### Direct model constructors

The three model classes can also be instantiated directly with the same
defaults. Use this when you need a specific type for `isinstance()` checks
or want to skip the `model` string dispatch:

```python
QuasiStaticBEM(defn, polar=None, n_psi_elements=36)
PittPetersModel(defn, polar=None, n_psi_elements=36)
OyeBEMModel(defn, polar=None, n_psi_elements=36, coupling_k=0.6)
```

In all cases `polar=None` means the polar is auto-built from `defn.airfoil`
(a `TabulatedPolar` if `polar_csv` is set, otherwise a `LinearPolar`).
Pass an explicit polar only when overriding the one implied by the rotor
definition.

`OyeBEMModel.coupling_k` (default `0.6`) is the empirical Oye coupling
factor `k` between the intermediate and quasi-steady inflow stages
(OpenFAST DBEMT default).

---

## Core loop

### `model.initial_rotor_state()`

Return a zero-initialised state of the correct type for this model:
`QuasiStaticRotorState`, `PittPetersRotorState`, or `OyeRotorState`.

### `model.compute_forces(inputs: RotorInputs, state) -> (AeroResult, state_derivative)`

Run one aerodynamic evaluation. Returns:

- `AeroResult` — forces, moments, and shaft torque for the current step.
- `state_derivative` — inflow state rates `d/dt`. Integrate this with
  the same `dt` used for `omega` to advance the inflow state:

```python
result, dstate = model.compute_forces(inputs, state)

# Advance inflow (explicit Euler; use semi-implicit for stiff cases)
arr = state.to_array() + dt * dstate.to_array()
state = state.from_array(arr)

# Advance rotor speed (caller owns the mechanical ODE)
omega += dt * omega_derivative(result.Q_spin, motor_torque_Nm, I_kgm2)
inputs = RotorInputs(..., omega_rad_s=omega, ...)
```

---

## `RotorInputs`

Constructed by keyword — all fields required.

| Field | Type | Description |
|---|---|---|
| `collective_rad` | float | Collective blade pitch [rad] |
| `tilt_lon` | float | Longitudinal swashplate tilt [rad]. Positive = nose-down |
| `tilt_lat` | float | Lateral swashplate tilt [rad]. Positive = roll right |
| `R_hub` | ndarray (3x3) | Rotation matrix: hub frame -> NED world |
| `v_hub_world` | ndarray (3,) | Hub velocity in NED world frame [m/s] |
| `wind_world` | ndarray (3,) | Ambient wind velocity in NED world frame [m/s] |
| `omega_rad_s` | float | Rotor speed [rad/s]. Caller advances this each step |
| `rho_kg_m3` | float | Air density [kg/m^3] (default 1.225) |
| `t` | float | Simulation time [s] (used for logging / future unsteady models) |

---

## `AeroResult`

Returned by `compute_forces`.

| Field | Type | Description |
|---|---|---|
| `F_world` | ndarray (3,) | Total rotor force in NED world frame [N]. `F_world[2] < 0` for upward lift |
| `M_orbital` | ndarray (3,) | Hub moments in NED world frame [N*m]. Roll about X, pitch about Y |
| `Q_spin` | float | Aerodynamic reaction torque on the shaft [N*m]. Positive opposes rotation. Feed directly to `omega_derivative` |
| `M_spin` | ndarray (3,) | Full spin moment vector in NED world frame [N*m] |

---

## Rotor states

All three states expose:

- `state.to_array() -> ndarray` — serialise inflow components to a 1-D array.
- `state.from_array(arr) -> state` — deserialise (returns a new instance).

### `QuasiStaticRotorState`

No inflow DOFs. `to_array()` returns an empty array.

### `PittPetersRotorState(lambda_0, lambda_c, lambda_s)`

| Attribute | Description |
|---|---|
| `lambda_0` | Uniform induced inflow ratio `v_i / (Omega*R)` |
| `lambda_c` | Cosine (longitudinal) inflow harmonic |
| `lambda_s` | Sine (lateral) inflow harmonic |

`to_array()` order: `[lambda_0, lambda_c, lambda_s]`.

### `OyeRotorState`

Per-annulus two-stage filter state. Constructed via
`model.initial_rotor_state()` (requires knowing `n_elements` from the
model). `to_array()` concatenates `[W_int..., W...]` (length
`2 * n_elements`).

---

## Mechanical ODE helpers

### `omega_derivative(Q_aero, motor_torque_Nm, I_ode_kgm2) -> float`

Return `d(omega)/dt = (motor_torque - Q_aero) / I` [rad/s^2].
Pass `AeroResult.Q_spin` directly for `Q_aero`.

### `euler_step_omega(omega, spin_angle, Q_aero, motor_torque_Nm, I_ode_kgm2, dt) -> (omega_new, spin_angle_new)`

Forward-Euler step for both `omega` and the spin angle.

---

## Trim solver

### `solve_trim_cyclic(aero, state, base_inputs, *, ...) -> TrimResult`

Find the `(tilt_lon, tilt_lat)` pair that drives in-plane hub moments
to `target_moment` (default `(0, 0)`) within `tolerance_Nm`.

`base_inputs` is a `RotorInputs` with the operating-point fields set;
`tilt_lon` and `tilt_lat` on it are the initial guess (both default 0).

| Keyword | Default | Description |
|---|---|---|
| `target_moment` | `(0.0, 0.0)` | Target `(Mx_hub, My_hub)` [N*m] |
| `tilt_lon_init` | `0.0` | Override initial tilt_lon guess [rad] |
| `tilt_lat_init` | `0.0` | Override initial tilt_lat guess [rad] |
| `tilt_min` | `-0.262` (-15 deg) | Lower cyclic bound [rad] |
| `tilt_max` | `+0.262` (+15 deg) | Upper cyclic bound [rad] |
| `tolerance_Nm` | `0.02` | Convergence threshold [N*m] |
| `max_iterations` | `50` | Newton iteration limit |
| `probe_rad` | `0.00873` (0.5 deg) | Finite-difference probe step [rad] |
| `dt_relax` | `0.005` | Timestep for inflow relaxation within each iteration [s] |
| `n_inflow_relax` | `100` | Inflow relax steps per iteration |
| `n_settle` | `0` | Extra settle steps after convergence |

### `relax_inflow(aero, state, inputs, *, n_steps=200, dt=0.005) -> state`

Advance the inflow state to quasi-steady equilibrium at fixed `inputs`.
Returns the settled state. Useful for initialising before a trim solve
or a time-domain run.

### `TrimResult`

| Attribute | Type | Description |
|---|---|---|
| `tilt_lon` | float | Trim longitudinal cyclic [rad] |
| `tilt_lat` | float | Trim lateral cyclic [rad] |
| `Mx_residual` | float | Hub roll moment at trim [N*m] |
| `My_residual` | float | Hub pitch moment at trim [N*m] |
| `converged` | bool | Whether tolerance was met |
| `iterations` | int | Newton iterations used |
| `final_state` | state | Settled inflow state at trim cyclic |

---

## Polars

For normal use you do not need to construct a polar directly — model
constructors and `create_aero` auto-build the polar from `defn.airfoil`
when `polar=None` (the default). The polar types are public for two
narrower purposes: (a) querying `cl_cd(alpha)` for debugging / plotting
without running a full BEM call, and (b) constructing a custom polar to
pass as the `polar=` override.

### `LinearPolar(CL0, CL_alpha_per_rad, CD0, alpha_stall_rad)`

Flat-plate linear lift with constant drag. Lift clips to `CL_max =
CL0 + CL_alpha * alpha_stall` beyond stall.

### `TabulatedPolar(alpha_rad, cl, cd)`

Cubic-spline interpolation over sampled `(alpha, CL, CD)` arrays.
Load from CSV with `dynbem.load_tabulated_polar(path)`.

Both types expose:

- `polar.cl_cd(alpha_rad: float) -> (cl, cd)` — evaluate the polar at
  a single angle of attack.

### `build_polar(airfoil) -> LinearPolar | TabulatedPolar`

Convert an `AirfoilProperties` to the matching polar type
(uses `polar_csv` if set, otherwise builds a `LinearPolar` from
`CL0` / `CL_alpha_per_rad` / `CD0` / `alpha_stall_deg`).

---

## Rotor definition

### `dynbem.rotor_definition.load(path) -> RotorDefinition`

Parse a rotor YAML file. Example files under `rotors/`.

```python
from dynbem.rotor_definition import load
defn = load("rotors/beaupoil_2026/rotor.yaml")
```

Key attributes on `RotorDefinition`:

| Attribute | Description |
|---|---|
| `defn.blade` | `BladeGeometry` (radius, chord, twist, n_elements) |
| `defn.airfoil` | `AirfoilProperties` (polar parameters or CSV path) |
| `defn.control` | `ControlProperties` or `None` (swashplate gain + phase) |

---

## Minimal hover example

```python
import numpy as np
import dynbem
from dynbem.rotor_definition import load

defn  = load("rotors/beaupoil_2026/rotor.yaml")
model = dynbem.create_aero(defn, model="pitt_peters")
state = model.initial_rotor_state()

omega = 28.0           # rad/s
dt    = 0.005          # s
I     = 0.08           # kg*m^2 (rotor inertia)

inputs = dynbem.RotorInputs(
    collective_rad=0.14,
    tilt_lon=0.0, tilt_lat=0.0,
    R_hub=np.eye(3),
    v_hub_world=np.zeros(3),
    wind_world=np.zeros(3),
    omega_rad_s=omega,
    rho_kg_m3=1.225,
    t=0.0,
)

for step in range(500):
    result, dstate = model.compute_forces(inputs, state)

    # advance inflow
    arr   = state.to_array() + dt * dstate.to_array()
    state = state.from_array(arr)

    # advance rotor speed (free-spin, no motor torque)
    omega   += dt * dynbem.omega_derivative(result.Q_spin, 0.0, I)
    inputs   = dynbem.RotorInputs(
        collective_rad=0.14, tilt_lon=0.0, tilt_lat=0.0,
        R_hub=np.eye(3), v_hub_world=np.zeros(3), wind_world=np.zeros(3),
        omega_rad_s=omega, rho_kg_m3=1.225, t=(step + 1) * dt,
    )

print(f"Thrust: {-result.F_world[2]:.1f} N")
print(f"omega:  {omega:.2f} rad/s")
```

## Cyclic trim example

```python
import numpy as np
import dynbem
from dynbem.rotor_definition import load

defn  = load("rotors/beaupoil_2026/rotor.yaml")
model = dynbem.create_aero(defn, model="pitt_peters")
state = model.initial_rotor_state()

base  = dynbem.RotorInputs(
    collective_rad=0.14, tilt_lon=0.0, tilt_lat=0.0,
    R_hub=np.eye(3),
    v_hub_world=np.zeros(3),
    wind_world=np.array([10.0, 0.0, 0.0]),   # 10 m/s headwind
    omega_rad_s=28.0, rho_kg_m3=1.225, t=0.0,
)

tr = dynbem.solve_trim_cyclic(model, state, base, tolerance_Nm=0.05)
if tr.converged:
    print(f"tilt_lon = {np.degrees(tr.tilt_lon):.2f} deg")
    print(f"tilt_lat = {np.degrees(tr.tilt_lat):.2f} deg")
    print(f"Mx={tr.Mx_residual:.3f} N*m  My={tr.My_residual:.3f} N*m")
```
