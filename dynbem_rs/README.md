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
fn compute_forces(&self, inputs: &RotorInputs, state: &Self::State)
    -> (AeroResult, Self::State);

fn initial_state(&self) -> Self::State;

fn inflow_taus(&self, inputs: &RotorInputs, state: &Self::State)
    -> Vec<f64>;
```

## Polar types

| Type | When to use |
|---|---|
| `LinearPolar` | Flat-plate / linear lift curve; constructed from CL0, CL_alpha, CD0, alpha_stall |
| `TabulatedPolar` | Interpolated from alpha/CL/CD arrays |

`LinearPolar::from_properties(props: &AirfoilProperties)` builds a
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
dynbem_rs = "0.2"
```

```rust
use dynbem_rs::{
    aero_io::{Mat3, RotorInputs, Vec3},
    aero_model::AeroModel,
    polar::LinearPolar,
    pitt_peters::PittPetersModel,
    rotor_definition::{
        AirfoilProperties, BladeGeometry, RotorDefinition,
    },
};
use std::f64::consts::PI;

let airfoil = AirfoilProperties {
    CL0: 0.0, CL_alpha_per_rad: 2.0 * PI, CD0: 0.01,
    alpha_stall_deg: 15.0, tip_loss: true,
};
let defn = RotorDefinition {
    blade: BladeGeometry::uniform(2, 0.15, 0.015, 0.025, 0.0, 12),
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
    t: 0.0,
};
let (result, _dstate) = model.compute_forces(&inputs, &state);
println!("Thrust: {:.1} N", -result.F_world.0[2]);
```

## License

MIT
