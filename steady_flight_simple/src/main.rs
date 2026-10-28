//! Real-time steady-flight visual sandbox.
//!
//! Main only wires modules and starts the Bevy app.

mod rendering;
mod simulation;

use bevy::prelude::*;
use simulation::{FlightSim, ModelChoice, SimRunner};
use std::env;

fn main() {
    let choice = env::args()
        .nth(1)
        .and_then(|s| ModelChoice::from_arg(&s))
        .unwrap_or(ModelChoice::Vpm);

    println!(
        "Launching real-time viewer model={} (close window to exit)",
        choice.as_str()
    );

    App::new()
        .insert_resource(SimRunner {
            sim: FlightSim::new(choice),
            accumulator_s: 0.0,
        })
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: format!("steady_flight_simple - {} (real-time)", choice.as_str()),
                resolution: (1280.0, 800.0).into(),
                ..default()
            }),
            ..default()
        }))
        .add_systems(Startup, rendering::setup_scene)
        .add_systems(
            Update,
            (
                simulation::update_simulation,
                rendering::update_visuals,
                rendering::update_particles,
                rendering::follow_camera,
            ),
        )
        .run();
}
