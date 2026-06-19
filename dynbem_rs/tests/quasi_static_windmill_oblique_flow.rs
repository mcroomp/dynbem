mod common;

use dynbem_rs::aero_io::{Mat3, RotorInputs, Vec3};
use dynbem_rs::aero_model::AeroModel;
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::quasi_static_bem::QuasiStaticBEM;

fn rawes_row122_inputs() -> RotorInputs {
    RotorInputs {
        collective_rad: -0.18146020320832462,
        tilt_lon: 0.012614825392536453,
        tilt_lat: 0.035447368174067954,
        R_hub: Mat3([
            [-0.007232662001129507, -0.9995813722829041, 0.02801372495419454],
            [0.684832857996051, -0.025365018889292906, -0.728258589019448],
            [0.7286642884136498, 0.013917471093409295, 0.684729624547885],
        ]),
        v_hub_world: Vec3::new(-0.5280840184646872, -0.17333482520461627, -0.4766214883089747),
        wind_world: Vec3::new(0.0, 10.0, 0.0),
        t: 51.102500000000254,
        rho_kg_m3: 1.225,
        omega_rad_s: 37.02311435435481,
    }
}

#[test]
fn rawes_row122_bem_force_should_remain_opposite_body_z() {
    let defn = common::beaupoil_rotor();
    let polar = LinearPolar::from_properties(&defn.airfoil);
    let model = QuasiStaticBEM::build(defn, 36, polar);
    let inputs = rawes_row122_inputs();
    let state = model.initial_state();

    let (result, _dstate) = model.compute_forces(&inputs, &state);
    let body_z = Vec3::new(inputs.R_hub.0[0][2], inputs.R_hub.0[1][2], inputs.R_hub.0[2][2]);
    let rel_wind_body = inputs.R_hub.transpose() * (inputs.wind_world - inputs.v_hub_world);
    let minus_f_dot_body_z = -result.F_world.dot(body_z);

    eprintln!("RAWES row 122 BEM reproduction");
    eprintln!("  collective_rad = {:+.9}", inputs.collective_rad);
    eprintln!("  tilt_lon       = {:+.9}", inputs.tilt_lon);
    eprintln!("  tilt_lat       = {:+.9}", inputs.tilt_lat);
    eprintln!("  omega_rad_s    = {:+.9}", inputs.omega_rad_s);
    eprintln!("  v_hub_world    = {:?}", inputs.v_hub_world.0);
    eprintln!("  wind_world     = {:?}", inputs.wind_world.0);
    eprintln!("  rel_wind_body  = {:?}", rel_wind_body.0);
    eprintln!("  body_z_world   = {:?}", body_z.0);
    eprintln!("  F_world        = {:?}", result.F_world.0);
    eprintln!("  -F dot body_z  = {:+.9}", minus_f_dot_body_z);
    eprintln!("  Q_spin         = {:+.9}", result.Q_spin);

    assert!(
        minus_f_dot_body_z > 0.0,
        "BEM returned force along +body_z: -F dot body_z = {minus_f_dot_body_z:+.9} N"
    );
}