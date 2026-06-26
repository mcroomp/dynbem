# dynbem_rs

Pure-Rust BEM / Pitt-Peters / Oye dynamic-inflow rotor aerodynamics library.
No pyo3, no numpy, no file I/O -- just the math.

This crate is the computational core behind the
[dynbem](https://github.com/mcroomp/dynbem) Python package, which wraps it
via PyO3 + maturin and provides the public Python API.

## Models

| Struct | Inflow model |
|---|---|
| `QuasiStaticBEM<P>` | Quasi-static BEM (steady annular momentum) |
| `PittPetersModel<P>` | Pitt-Peters 3-state dynamic inflow |
| `OyeBEMModel<P>` | Oye 2-stage annular dynamic inflow (DBEMT equivalent) |

All three are generic over `P: Polar`. Pick `LinearPolar` for a flat-plate
lift curve or `TabulatedPolar` for interpolated alpha/CL/CD data; implement
the `Polar` trait to supply your own.

All three implement the `AeroModel` trait:

```rust
use dynbem_rs::aero_model::IntegrationMethod;

// Built-in state step with selectable inflow integrator.
fn step(&self, inputs: &RotorInputs, state: &Self::State, dt: f64, method: IntegrationMethod)
    -> (AeroResult, Self::State);

fn initial_state(&self) -> Self::State;

// If you want to integrate inflow yourself, use these
fn compute_forces(&self, inputs: &RotorInputs, state: &Self::State)
    -> (AeroResult, Self::State);

fn inflow_taus(&self, inputs: &RotorInputs, state: &Self::State)
    -> Vec<f64>;

```

Integration methods:

- `IntegrationMethod::SemiImplicitEuler`
- `IntegrationMethod::ExplicitEuler`
- `IntegrationMethod::ExponentialRelaxation`

## Polar types

| Type | When to use |
|---|---|
| `LinearPolar` | Flat-plate / linear lift curve; constructed from CL0, CL_alpha, CD0, alpha_stall |
| `TabulatedPolar` | Interpolated from alpha/CL/CD arrays |

`LinearPolar::from_properties(props: &LinearPolarParameters)` builds a
`LinearPolar` directly from the airfoil block of a `RotorDefinition`.

## Custom polar types

Implement the `Polar` trait to supply your own polar:

```rust
pub trait Polar {
    fn cl_cd(&self, alpha_rad: f64) -> (f64, f64);
}
```

Then pass it to any model constructor:

```rust
let model = PittPetersModel::build(defn, 36, my_polar);
```

## Coordinate system

NED frame throughout. Rotor spins CCW when viewed from above (American /
Bell / Sikorsky convention). See the repository README for full sign
conventions, cyclic pitch convention, and Pitt-Peters L-matrix derivation.

## Usage

Add to `Cargo.toml`:

```toml
[dependencies]
dynbem_rs = "0.3"
```

```rust
use dynbem_rs::{
    aero_io::{Mat3, RotorInputs, Vec3},
    aero_model::AeroModel,
    polar::LinearPolar,
    pitt_peters::PittPetersModel,
    rotor_definition::
        LinearPolarParameters, BladeGeometry, RotorDefinition,
    },
};
use std::f64::consts::PI;

let airfoil = LinearPolarParameters {
    CL0: 0.0, CL_alpha_per_rad: 2.0 * PI, CD0: 0.01,
    alpha_stall_deg: 15.0,
};
let defn = RotorDefinition {
    blade: BladeGeometry::uniform(2, 0.15, 0.015, 0.025, 0.0, 12),  // tip_loss defaults to true
    airfoil,
    control: None,
    name: "my_rotor".into(),
    description: String::new(),
};

let polar = LinearPolar::from_properties(&defn.airfoil);
let model = PittPetersModel::build(defn, 36, polar);
let state = model.initial_state();

let inputs = RotorInputs {
    collective_rad: 0.14,
    tilt_lon: 0.0,
    tilt_lat: 0.0,
    R_hub: Mat3::eye(),
    v_hub_world: Vec3::zero(),
    wind_world: Vec3::zero(),
    omega_rad_s: 628.0,
    rho_kg_m3: 1.225,
};
let dt = 0.01;
let (result, _state_next) = model.step(
    &inputs,
    &state,
    dt,
    dynbem_rs::aero_model::IntegrationMethod::SemiImplicitEuler,
);
println!("Thrust: {:.1} N", -result.F_world.0[2]);
```

## Inputs: `RotorInputs`

All quantities are in SI units. The frame is NED (North-East-Down) throughout.

| Field | Type | Description |
|---|---|---|
| `collective_rad` | `f64` | Blade collective pitch angle (rad). Positive increases pitch and thrust. In servo-flap mode this is the flap command, not direct blade pitch. |
| `tilt_lon` | `f64` | Longitudinal swashplate tilt (rad). Positive tilts the disk nose-down (forward stick). |
| `tilt_lat` | `f64` | Lateral swashplate tilt (rad). Positive rolls the disk right. |
| `R_hub` | `Mat3` | Rotation matrix from hub frame to world (NED) frame. `R_hub * [0,0,1]` gives the hub spin axis in world coordinates. Pass `Mat3::eye()` for a level rotor with +Z spin axis pointing down. |
| `v_hub_world` | `Vec3` | Hub-centre velocity in world (NED) frame, m/s. Positive +Z is downward. |
| `wind_world` | `Vec3` | Ambient wind velocity in world (NED) frame, m/s. Only the relative velocity `wind_world - v_hub_world` matters -- equal and opposite values produce identical aerodynamics. |
| `omega_rad_s` | `f64` | Rotor angular velocity (rad/s). Positive for the correct spin direction (CCW from above, American convention). The caller advances this externally: `omega += dt * (motor_torque - Q_spin) / I_kgm2`. |
| `rho_kg_m3` | `f64` | Air density (kg/m^3). Use 1.225 for ISA sea level. |
| `t` | `f64` | *(Python API only)* Simulation time (s) carried as a convenience field; not present in the Rust `RotorInputs` struct and not read by any model. |

### Wind vs velocity

`v_hub_world` and `wind_world` cover the same physical concept: what matters is the
relative velocity of the air with respect to the hub. The two inputs are combined
immediately as `v_rel = wind_world - v_hub_world` inside `kinematics()` and the
individual values are never used again. Providing both separately is a convenience
so that a flight simulator can pass hub velocity and ambient wind independently
without pre-subtracting.

## Outputs: `AeroResult`

All quantities are in SI units in the world (NED) frame.

| Field | Type | Description |
|---|---|---|
| `F_world` | `Vec3` | Net aerodynamic force on the rotor hub, N, in NED world frame. In hover with the rotor level, thrust appears as a negative Z component (i.e. `F_world.0[2] < 0` lifts the aircraft). |
| `M_hub_world` | `Vec3` | Aerodynamic hub moment (pitching/rolling) about the hub centre, N·m, in NED world frame. Corresponds to `M_x`/`M_y` in standard helicopter notation. This is what the airframe sees after blade flap reduction (if `FlapProperties` is set). In NED, positive Mx is roll-right, positive My is pitch-up. |
| `Q_spin` | `f64` | Aerodynamic reaction torque opposing rotation, N·m. Always positive for a powered rotor (torque opposes spin). Use in the angular acceleration equation: `omega += dt * (motor_torque - Q_spin) / I_kgm2`. |
| `M_spin` | `Vec3` | `Q_spin` expressed as a vector along the rotor spin axis in world frame, N·m. Direction is `R_hub * [0,0,1]` (the hub +Z axis in world coordinates). Convenience for rigid-body integrators that accept a 3-vector torque input. |

### State Stepping

Use `step(...)` with an explicit integration method.

`AeroModel::step` runs `compute_forces` then advances the inflow state with
the selected method in one call:

```rust
// Returns forces at state(t) and the integrated state at t+dt.
use dynbem_rs::aero_model::IntegrationMethod;

let (result, new_state) = model.step(
    &inputs,
    &state,
    dt,
    IntegrationMethod::SemiImplicitEuler,
);
```

Available methods:

- `IntegrationMethod::SemiImplicitEuler` (recommended default for dynamic inflow)
- `IntegrationMethod::ExplicitEuler` (useful for baseline matching and some regressions)
- `IntegrationMethod::ExponentialRelaxation` (exact for frozen first-order lag over the step)

For each inflow DOF `i` the update is:

```
new_lam[i] = (lam[i] + dt * dlam[i]) / (1 + dt / tau[i])
```

where `tau[i]` comes from `inflow_taus()` when using `SemiImplicitEuler`.
For the scalar first-order lag form,
this update is unconditionally stable in `dt`; in the full coupled nonlinear
system it is typically much more robust than explicit Euler and strongly damps
stiff inflow modes. When `tau[i]` is infinite (quasi-static BEM) it reduces
to plain explicit Euler.

### Advanced: Integrate Yourself

Use this when you need higher-order accuracy (RK4), adaptive timestepping, or
want to couple the inflow ODE into a larger system integrator.

`compute_forces` returns `(AeroResult, Self::State)`. For dynamic-inflow models
(Pitt-Peters, Oye) the second return value is the **time-derivative** of the
inflow state, not the new state directly.

```rust
// Single Euler step -- simplest custom integration.
let (result, d_state) = model.compute_forces(&inputs, &state);
let inflow = state.get_inflow();
let d_inflow = d_state.get_inflow();
let new_inflow: Vec<f64> = inflow.iter()
    .zip(d_inflow.iter())
    .map(|(lam, dlam)| lam + dt * dlam)
    .collect();
let mut new_state = state.clone();
new_state.set_inflow(new_inflow);
```

`inflow_taus()` returns the time constant for each DOF; use these to
set sub-step size or as the damping term in a semi-implicit scheme.

For the quasi-static BEM model `d_state` is always zero and both options
are equivalent.

## License

MIT
