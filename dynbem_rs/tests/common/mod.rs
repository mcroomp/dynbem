// Shared fixtures for integration tests.
#![allow(dead_code)]

use dynbem_rs::polar::LinearPolar;
use dynbem_rs::quasi_static_bem::QuasiStaticBEM;
use dynbem_rs::rotor_definition::{
    BladeGeometry, ControlProperties, LinearPolarParameters, PitchActuation, RotorDefinition,
};

/// Beaupoil 2026 / RAWES rotor definition (matches rotors/beaupoil_2026/rotor.yaml).
pub fn beaupoil_rotor() -> RotorDefinition {
    RotorDefinition {
        blade: BladeGeometry {
            n_blades: 4,
            radius_m: 2.5,
            root_cutout_m: 0.5,
            chord_m: 0.20,
            twist_deg: 0.0,
            n_elements: 10,
            tip_loss: true,
            r_stations_m: Vec::new(),
            chord_stations_m: Vec::new(),
            twist_stations_deg: Vec::new(),
        },
        airfoil: LinearPolarParameters {
            CL0: 0.393,
            CL_alpha_per_rad: 5.79,
            CD0: 0.0079,
            alpha_stall_deg: 13.0,
        },
        control: Some(ControlProperties {
            swashplate_pitch_gain_rad: 0.3,
            swashplate_phase_deg: Some(0.0),
        }),
        pitch_actuation: PitchActuation::DirectMechanical,
        flap: None,
        name: "beaupoil_2026".to_string(),
        description: String::new(),
    }
}

pub fn qs_model() -> QuasiStaticBEM<LinearPolar> {
    let defn = beaupoil_rotor();
    let polar = LinearPolar::from_properties(&defn.airfoil);
    QuasiStaticBEM::build(defn, 36, polar)
}
