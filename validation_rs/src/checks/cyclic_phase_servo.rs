// Cyclic phase shift: direct-mechanical vs servo-flap (Kaman-type) actuation,
// using the Beaupoil 2026 rotor geometry and airfoil.
//
// Physics summary
// ---------------
// Direct-mechanical: the swashplate sets blade pitch directly. tilt_lon > 0
// peaks blade pitch at the tail (psi = pi), producing maximum thrust at the
// tail and a nose-down pitching moment: My < 0, Mx ~ 0.
//
// Servo-flap: the swashplate drives a trailing-edge flap, which in turn
// drives blade feathering through an aerodynamic pitching moment. At 1/rev,
// the feathering DOF acts as a resonant spring-damper system that introduces
// a ~90 deg phase lag in the direction of rotation. For a CCW-from-above
// rotor, the same tilt_lon command therefore peaks blade pitch ~90 deg later
// at psi = 3*pi/2 (the retreating / right side), producing maximum thrust on
// the right and a roll-right moment: Mx > 0, My ~ 0.
//
// Assertions
// ----------
// Direct:    |My| > |Mx|  and  My < 0   (pitching dominates, nose-down)
// ServoFlap: |Mx| > |My|  and  Mx > 0   (rolling dominates, roll-right)
//
// Rotor: Beaupoil 2026 (4-blade, R=2.5 m, chord=0.20 m, SG6040 airfoil).
// Hover at collective=5 deg, omega=20 rad/s, rho=1.225.
// VPM: fast_test preset with sigma scaled to 1.5*chord = 0.30 m.

use crate::helpers::*;
use crate::report::Report;
use dynbem_rs::cyclic::ControlGains;
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::rotor_definition::{
    BladeGeometry, LinearPolarParameters, PitchActuation, RotorDefinition, ServoFlapActuation,
    ServoFlapGeometry,
};
use dynbem_rs::vpm::{FlightCondition, VpmRotor, VpmRotorConfig};

const B_OMEGA: f64 = 20.0; // rad/s, near measured equilibrium at V_wind=10 m/s
const B_R: f64 = 2.5;
const B_ROOT: f64 = 0.5;
const B_CHORD: f64 = 0.20;

fn beaupoil_defn(actuation: PitchActuation) -> RotorDefinition {
    RotorDefinition {
        blade: BladeGeometry {
            n_blades: 4,
            radius_m: B_R,
            root_cutout_m: B_ROOT,
            chord_m: B_CHORD,
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
        control: None, // gain=1, phase=0 -- isolates actuator effect
        pitch_actuation: actuation,
        flap: None,
        name: "beaupoil_2026_test".to_string(),
        description: String::new(),
    }
}

fn beaupoil_servo() -> ServoFlapActuation {
    // Feathering + damper architecture (Kaman ideal: axis at AC, no aero spring).
    // The bearing damper alone sets the cyclic phase lag. The cyclic test drives
    // collective=0, so no DC restoring is needed; ac_offset=0 and blade_Cm_AC=0.
    ServoFlapActuation {
        I_theta_kgm2: 0.05,
        damper_Nms_per_rad: 2.0,
        ac_offset_m: 0.0, // Kaman ideal: feathering axis at AC
        blade_Cm_AC: 0.0,
        flap: ServoFlapGeometry {
            C_M_delta_per_rad: -1.0,
            r_inner_m: 1.0,       // mid-span flap start
            r_outer_m: 0.9 * B_R, // near-tip flap end
        },
    }
}

fn beaupoil_fast_rotor(defn: &RotorDefinition) -> VpmRotor<LinearPolar> {
    // The Beaupoil has a much larger chord/radius than the default test rotor.
    // VPM stability requires sigma/delta_s > ~0.5 where delta_s is the arc
    // spacing between consecutive shed particles at the tip (= R * dpsi).
    //
    // With dpsi = PI/6 (30 deg/step, 12 steps/rev):
    //   delta_s = 2.5 * PI/6 = 1.31 m
    //   sigma = 0.45 m  ->  sigma/delta_s = 0.34  (marginally stable)
    // max_particles = 2500 retains ~2 revolutions of wake:
    //   4 blades x 21 particles/step x 12 steps/rev x 2.5 rev = 2520
    //
    // NAN_DEBUG=1: enable the scalar nan-asserting induction path so any bad
    // particle pair is reported with full diagnostics before the NaN propagates.
    let nan_debug = std::env::var("NAN_DEBUG")
        .map(|v| v == "1")
        .unwrap_or(false);
    let cfg = VpmRotorConfig {
        max_particles: 2500,
        sigma: 0.45,
        nonlinear_lifting_line: false,
        use_rayon: !nan_debug, // single-threaded when asserting
        use_scalar_nan_check: nan_debug,
        ..VpmRotorConfig::fast_test()
    };
    VpmRotor::new(defn, polar_for(&defn.airfoil), ControlGains::default(), cfg)
}

/// Step-by-step march that panics on the first NaN result, printing which step
/// produced it.  Set NAN_DEBUG=1 in the environment to enable this path
/// (otherwise it falls through to the normal batch march).
fn march_nan_debug(
    rotor: &VpmRotor<LinearPolar>,
    fc: &FlightCondition,
    dt: f64,
    n_steps: usize,
    label: &str,
) -> dynbem_rs::vpm::VpmRotorResult {
    use dynbem_rs::vpm::VpmRotorState;
    let mut state: Option<VpmRotorState> = None;
    let avg_window = (n_steps / 2).max(1);
    let mut results = Vec::with_capacity(avg_window);

    for step in 0..n_steps {
        let (res, next_state) =
            rotor.step_one(fc, state.as_ref().unwrap_or(&VpmRotorState::default()), dt);

        if res.thrust.is_nan() || res.mx_hub.is_nan() || res.my_hub.is_nan() {
            panic!(
                "[{}] NaN at step {}: thrust={} mx={} my={} n_particles={}",
                label, step, res.thrust, res.mx_hub, res.my_hub, res.n_particles
            );
        }

        if step >= n_steps - avg_window {
            results.push(res);
        }
        state = Some(next_state);
    }

    // Average over the window (mirrors march() behaviour)
    let n = results.len() as f64;
    dynbem_rs::vpm::VpmRotorResult {
        thrust: results.iter().map(|r| r.thrust).sum::<f64>() / n,
        torque: results.iter().map(|r| r.torque).sum::<f64>() / n,
        mx_hub: results.iter().map(|r| r.mx_hub).sum::<f64>() / n,
        my_hub: results.iter().map(|r| r.my_hub).sum::<f64>() / n,
        n_particles: results.last().map(|r| r.n_particles).unwrap_or(0),
        wake_centroid: results.last().map(|r| r.wake_centroid).unwrap_or([0.0; 3]),
    }
}

pub fn check_cyclic_phase_servo(r: &mut Report) {
    r.begin_module(
        "cyclic_phase_servo",
        "Cyclic phase: direct-mech My dominates, servo-flap Mx dominates (90-deg lag, Beaupoil rotor)",
    );

    let tilt_lon_rad = 3.0_f64.to_radians();
    // dpsi = PI/6 (30 deg/step, 12 steps/rev): gives sigma/delta_s ~ 0.34
    // which keeps the Beaupoil wake stable. Run 12 revolutions (144 steps),
    // averaged over the last 6 (72 steps).
    let dpsi = std::f64::consts::PI / 6.0;
    let dt = dpsi / B_OMEGA;
    let n_steps = 12 * 12; // 12 revolutions x 12 steps/rev

    let fc = FlightCondition {
        collective_rad: 5.0_f64.to_radians(),
        tilt_lon: tilt_lon_rad,
        tilt_lat: 0.0,
        v_hub: [0.0, 0.0, 0.0],
        omega_rad_s: B_OMEGA,
        rho: RHO,
    };

    let nan_debug = std::env::var("NAN_DEBUG")
        .map(|v| v == "1")
        .unwrap_or(false);

    // --- direct-mechanical: pitching moment (My) should dominate ---
    let defn_direct = beaupoil_defn(PitchActuation::DirectMechanical);
    let rotor_direct = beaupoil_fast_rotor(&defn_direct);
    let res_d = if nan_debug {
        march_nan_debug(&rotor_direct, &fc, dt, n_steps, "direct")
    } else {
        rotor_direct.march(&fc, None, dt, n_steps).0
    };

    // For servo-flap, drive the flap with cyclic only (collective_rad=0).
    // With no aero spring (ac_offset=0) a constant collective flap command has
    // no DC restoring and the feathering integrates without bound, driving the
    // blade into deep stall. Setting collective=0 isolates the cyclic (phase)
    // response and keeps theta_f within the linear polar range.
    let fc_sf = FlightCondition {
        collective_rad: 0.0,
        ..fc
    };

    r.info("direct", "mx_hub_Nm", res_d.mx_hub, f64::NAN);
    r.info("direct", "my_hub_Nm", res_d.my_hub, f64::NAN);
    r.assert_bool(
        "direct",
        "my_lt_0_nose_down",
        res_d.my_hub,
        0.0,
        res_d.my_hub < 0.0,
        &format!(
            "direct: tilt_lon>0 should give My<0 (nose-down), got My={:.3} Nm",
            res_d.my_hub
        ),
    );
    r.assert_bool(
        "direct",
        "pitch_dominates_roll",
        res_d.my_hub.abs(),
        res_d.mx_hub.abs(),
        res_d.my_hub.abs() > res_d.mx_hub.abs(),
        &format!(
            "direct: pitching should dominate: |My|={:.3} > |Mx|={:.3} Nm",
            res_d.my_hub.abs(),
            res_d.mx_hub.abs()
        ),
    );

    // --- servo-flap: 90-deg phase lag shifts response to rolling (Mx) ---
    let defn_sf = beaupoil_defn(PitchActuation::ServoFlap(beaupoil_servo()));
    let rotor_sf = beaupoil_fast_rotor(&defn_sf);
    let (res_s, state_s) = if nan_debug {
        let r = march_nan_debug(&rotor_sf, &fc_sf, dt, n_steps, "servo_flap");
        (r, None)
    } else {
        let (r, s) = rotor_sf.march(&fc_sf, None, dt, n_steps);
        (r, Some(s))
    };

    // Diagnostic: print settled theta_f (feathering angle) for each blade.
    // theta_f(psi_b) for CCW rotation: positive delta_theta_1s -> peaks at psi=pi/2
    // (LEFT side, -Y), which should give Mx > 0 (roll-right) from Mx = sum(r*dT*sin(psi)).
    if let Some(ref s) = state_s {
        if let Some(ref tf) = s.theta_f {
            eprintln!("servo_flap theta_f (rad): {:?}", tf);
            eprintln!("  psi at averaging end = {:.3} rad", s.psi);
        }
    }

    // Cross-check sign convention: direct tilt_lat>0 should give Mx>0 (roll-right).
    {
        let fc_lat = FlightCondition {
            tilt_lon: 0.0,
            tilt_lat: tilt_lon_rad, // same magnitude, lateral
            ..fc
        };
        let res_lat = rotor_direct.march(&fc_lat, None, dt, n_steps).0;
        eprintln!(
            "direct tilt_lat>0: Mx={:.3} My={:.3} (Mx>0 expected for roll-right)",
            res_lat.mx_hub, res_lat.my_hub
        );
    }

    r.info("servo_flap", "mx_hub_Nm", res_s.mx_hub, f64::NAN);
    r.info("servo_flap", "my_hub_Nm", res_s.my_hub, f64::NAN);
    r.assert_bool(
        "servo_flap",
        "mx_gt_0_roll_right",
        res_s.mx_hub,
        0.0,
        res_s.mx_hub > 0.0,
        &format!(
            "servo-flap: 90-deg lag should give Mx>0 (roll-right), got Mx={:.3} Nm",
            res_s.mx_hub
        ),
    );
    r.assert_bool(
        "servo_flap",
        "roll_dominates_pitch",
        res_s.mx_hub.abs(),
        res_s.my_hub.abs(),
        res_s.mx_hub.abs() > res_s.my_hub.abs(),
        &format!(
            "servo-flap: rolling should dominate: |Mx|={:.3} > |My|={:.3} Nm",
            res_s.mx_hub.abs(),
            res_s.my_hub.abs()
        ),
    );
}
