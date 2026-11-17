//! Real-time steady-flight visual sandbox.
//!
//! Main only wires modules and starts the Bevy app.

mod rendering;
mod simulation;

use bevy::prelude::*;
use simulation::{FlightSim, ModelChoice, SimRunner};
use std::env;

fn main() -> Result<(), String> {
    let choice = match env::args().nth(1) {
        Some(arg) => ModelChoice::from_arg(&arg)
            .ok_or_else(|| format!("Unknown model {arg:?}; use bem, oye, pitt_peters, or vpm"))?,
        None => ModelChoice::Vpm,
    };

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
                resolution: (1280.0_f32, 800.0).into(),
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
                rendering::follow_camera,
                rendering::update_particles,
            )
                .chain(),
        )
        .run();
    Ok(())
}
