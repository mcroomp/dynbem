# steady_flight_simple

A Bevy 3D sandbox for the four Rust aero models. It renders a rotor, a tether,
and (with VPM) wake particles. A graphical desktop and graphics driver are
required; audio and gamepad dependencies are disabled.

Source layout:

- `src/main.rs`: app entrypoint and wiring
- `src/simulation.rs`: flight dynamics, control, and aero model stepping
- `src/rendering.rs`: Bevy scene setup and visualization systems

Run from the repository root in release mode (VPM is expensive in debug):

```text
cargo run --release -p steady_flight_simple -- bem
cargo run --release -p steady_flight_simple -- oye
cargo run --release -p steady_flight_simple -- pitt_peters
cargo run --release -p steady_flight_simple -- vpm
```

The default model is `vpm`; unknown model names are errors. Close the window
to stop. There are no interactive flight controls.

## Model scope

The hardcoded four-blade rotor has a 2.5 m radius and fixed 20.148 rad/s speed.
The 5 kg body starts on a 60 m radius at an elevation near the 35 deg PID
setpoint. The tether is a constant-magnitude force toward the origin, not a
length constraint. Translation and attitude use simplified damping and
acceleration/rate limits. This is a visualization/control sandbox, not an
empirically validated steady-flight or tether model.

Simulation uses fixed 0.02 s steps, at most three per frame. Excess backlog
is discarded, so slow machines show slow-motion rather than accumulating
unbounded work. Simulation, rotor transforms, camera, and wake rendering run
in that order. NED simulation coordinates are mapped to Bevy's Y-up frame;
blades and wake use the same CCW-from-above azimuth convention.

Headless smoke and coordinate-mapping tests:

```text
cargo test --release -p steady_flight_simple
```
