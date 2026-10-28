# steady_flight_simple

Standalone steady-flight sandbox executable moved out of `dynbem_rs/bin`.
This version opens a real-time 3D window using Bevy and renders the rotor and tether.

Source layout:

- `src/main.rs`: app entrypoint and wiring
- `src/simulation.rs`: flight dynamics, control, and aero model stepping
- `src/rendering.rs`: Bevy scene setup and visualization systems

Run:

cargo run -p steady_flight_simple -- bem
cargo run -p steady_flight_simple -- oye
cargo run -p steady_flight_simple -- pitt_peters
cargo run -p steady_flight_simple -- vpm

Close the window to stop the program.
