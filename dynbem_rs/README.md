# dynbem_rs

Pure-Rust BEM / Pitt-Peters / Oye / VPM rotor aerodynamics library.
No pyo3, no numpy, no file I/O -- just the math.

This crate is the computational core behind the
[dynbem](https://github.com/mcroomp/dynbem) Python package, which wraps it
via PyO3 + maturin and provides the public Python API.

## Models

Algebraic dynamic-inflow / BEM models (one converged evaluation per call):

| Struct | Inflow model |
|---|---|
| `QuasiStaticBEM<P>` | Quasi-static BEM (steady annular momentum) |
| `PittPetersModel<P>` | Pitt-Peters 3-state dynamic inflow |
| `OyeBEMModel<P>` | Oye 2-stage annular dynamic inflow (DBEMT equivalent) |

Wake-resolving model (free-wake, marched one azimuth step per call):

| Struct | Model |
|---|---|
| `VpmRotor<P>` | Vortex particle method free-wake, forward flight (see [`../docs/VPM_DESIGN.md`](../docs/VPM_DESIGN.md)) |

All are generic over `P: Polar`. Pick `LinearPolar` for a flat-plate
lift curve or `TabulatedPolar` for interpolated alpha/CL/CD data; implement
the `Polar` trait to supply your own.

All implement the `AeroModel` trait:

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

The three algebraic models return a converged answer from a single
`compute_forces` call and integrate a small scalar inflow state. `VpmRotor`
is different: it carries the whole free wake as its state, so `compute_forces`
panics (a free wake has no single-shot evaluation) and `AeroModel::step`
advances the wake by exactly one azimuth increment per call (`dt` is the
convection step; the caller drives the loop). Use `VpmRotor::march` to settle
a periodic wake and `VpmRotor::step_one` for a single advance. See
[`../docs/VPM_DESIGN.md`](../docs/VPM_DESIGN.md) for the full design.

Integration methods:

- `IntegrationMethod::SemiImplicitEuler`
- `IntegrationMethod::ExplicitEuler`
- `IntegrationMethod::ExponentialRelaxation`

## Blade pitch actuation

`RotorDefinition::pitch_actuation` (`PitchActuation` enum) selects how the
swashplate commands reach the blade:

- `DirectMechanical` (default): the swashplate sets blade pitch directly.
- `ServoFlap(ServoFlapActuation)`: a trailing-edge servo-flap drives a free
  feathering DOF (feathering + damper architecture); the solved feathering
  angle replaces the direct swashplate pitch path. The quasi-static harmonic
  solve is in `servoflap.rs`; the time-domain feathering ODE is in
  `vpm_rotor.rs`. All `ServoFlapActuation` constants are physical and
  measurable (pitch inertia, bearing damper, AC offset aero spring, blade
  camber moment) -- see [`../AGENTS.md`](../AGENTS.md) for the sign
  conventions and the two servo-flap architectures.

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
dynbem_rs = "0.5"
```

```rust
use dynbem_rs::{
    aero_io::{Mat3, RotorInputs, Vec3},
    aero_model::AeroModel,
    polar::LinearPolar,
    pitt_peters::PittPetersModel,
    rotor_definition::{
        BladeGeometry, LinearPolarParameters, PitchActuation, RotorDefinition,
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
    pitch_actuation: PitchActuation::DirectMechanical,
    flap: None,
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

## Developer notes

### Hard rules

1. **No `pyo3` / `numpy` imports.** Python-facing helpers belong in
   [`../dynbem/`](../dynbem/), not here.
2. **No file IO outside `rotor_definition.rs`.** The YAML fields in
   `RotorDefinition` are parsed by the Python layer (`dynbem.rotor_definition`
   via `yaml.safe_load`); the Rust struct is populated by the PyO3 glue.
   Math modules (`bem_common`, `pitt_peters`, `oye`, `servoflap`, `polar`,
   `cyclic`, `trim`, `common`) must stay free of `std::fs` / `serde` so
   they remain embeddable and decoupled from file-format concerns.
3. **Public API stability matters.** The `dynbem/` glue crate depends on
   every public field, struct, and function here. Renaming or moving
   things requires matching edits there. Prefer additive changes.
4. **Sign conventions are NOT documented in this crate.** They live in
   [`../AGENTS.md`](../AGENTS.md) and [`../docs/BEM_COMMON.md`](../docs/BEM_COMMON.md).
   Refer to those; do not duplicate.

### Module map

    src/
    +-- lib.rs                public module declarations
    +-- aero_io.rs            Vec3, Mat3, RotorInputs, AeroResult
    +-- aero_model.rs         AeroModel trait + RotorStateExt trait + IntegrationMethod
    +-- bem_common.rs         RadialGrid, PolarTable, kinematics(), element_force(),
    |                         run_psi_loop<K: PsiKernel>, assemble_result(),
    |                         apply_flap_reduction()
    +-- common.rs             numerical floors (EPS_*), vrs_lambda1 VRS polynomial
    +-- cyclic.rs             swashplate -> theta_1c, theta_1s mapping
    +-- oye.rs                OyeBEMModel (annular 2-stage filter)
    +-- pitt_peters.rs        PittPetersModel (3-state L-matrix ODE)
    +-- polar.rs              LinearPolar, TabulatedPolar, Polar trait
    +-- quasi_static_bem.rs   QuasiStaticBEM + Ning windmill Brent solver
    +-- rotor_definition.rs   all RotorDefinition types (blade, airfoil, control,
    |                         servo-flap, flap, inertia)
    +-- servoflap.rs          quasi-static servo-flap-driven feathering model
    +-- trim.rs               solve_trim_cyclic<M>, relax_inflow<M> (generic over model)
    +-- vpm.rs                ParticleField + regularized Biot-Savart engine
    |                         (private; used only by vpm_rotor)
    +-- vpm_rotor.rs          VpmRotor free-wake forward-flight coupling
                              (see docs/VPM_DESIGN.md)

    bin/
    +-- rotor_profile.rs      per-step timing across all models; VPM direct vs
    |                         Barnes-Hut at matched N, sequential vs parallel
    |                         (`--long` sweep, `--seq` / `--par`)
    +-- bh_profile.rs         Barnes-Hut velocity-eval microbenchmark
    +-- profile_kernels.rs    low-level kernel timing

State types (`QuasiStaticRotorState`, `PittPetersRotorState`, `OyeRotorState`,
`VpmRotorState`) are defined in each model's own module. `RotorStateExt`
(the serialization trait) is declared in `aero_model.rs`.

### Adding a new aero model

1. Add `src/foo.rs` with the model struct and `impl AeroModel for FooModel`.
   Implement `AeroModel::step`, `compute_forces`, `initial_state`,
   `inflow_taus`.
2. Add `FooRotorState` to `foo.rs` and the `RotorStateExt` impl there.
   `RotorStateExt` serializes **inflow states only** via
   `get_inflow()` / `set_inflow(Vec<f64>)`. There are no mechanical
   fields; `omega_rad_s` is passed through `RotorInputs` on every call.
3. Add `pub mod foo;` to `lib.rs`.
4. Add a wrapper newtype in [`../dynbem/src/wrappers.rs`](../dynbem/src/wrappers.rs)
   (mark `subclass = true` if Python should auto-build the polar) and an
   `AeroAny` variant in [`../dynbem/src/trim_py.rs`](../dynbem/src/trim_py.rs).
5. Wire the model name into `create_aero()` in
   [`../dynbem/python/dynbem/factory.py`](../dynbem/python/dynbem/factory.py).

### Hot-path conventions

- Once-per-call kinematics prelude lives in `bem_common::kinematics`.
  Result assembly is `bem_common::assemble_result`. Per-element BEM
  integrand is `bem_common::element_force` (`#[inline(always)]`). All
  three BEM models call these -- do not duplicate.
- The psi x r sweep is one generic function `bem_common::run_psi_loop<K: PsiKernel>`.
  Pitt-Peters and Oye each implement `PsiKernel` for their own `lam_local`
  formula. Monomorphization + `#[inline(always)]` on the trait methods give
  the same codegen as a hand-rolled loop with no runtime dispatch.
- Plain `for` loops over `&[f64]`, no SIMD intrinsics. LLVM autovectorizes
  the if-converted bodies.
- `Vec3`/`Mat3` are `Copy` newtypes around plain f64 arrays. Operators
  lower to the same scalar FMA chains as hand-rolled index arithmetic.

### Numerical floors

All in `common.rs`:

| Constant | Value | Purpose |
|---|---|---|
| `EPS_DENOM` | 1e-9 | generic division / ratio safety |
| `EPS_OMEGA_R` | 1e-6 | rotor-not-spinning threshold |
| `MIN_LOSS_FACTOR` | 1e-4 | Prandtl tip+hub loss floor |
| `MASS_FLOW_HOVER_FLOOR_FRAC` | 1e-2 | mass-flow floor at hover / zero thrust |
| `VRS_DESCENT_THRESHOLD` | 1e-3 | VRS detection guard against hover chattering |
| `MU_T_FLOOR` | 0.05 | Pitt-Peters L-matrix denominator floor |

Tuned to keep hover / climb / descent / VRS / autorotation all stable in one
code path. Do not change without running the full suite (`uv run pytest tests/ -q`).

## License

MIT
