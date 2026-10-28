use bevy::prelude::{Res, ResMut, Resource, Time};
use dynbem_rs::aero_io::{AeroResult, Mat3, RotorInputs, Vec3 as DVec3};
use dynbem_rs::aero_model::{AeroModel, IntegrationMethod};
use dynbem_rs::cyclic::ControlGains;
use dynbem_rs::oye::{OyeBEMModel, OyeRotorState};
use dynbem_rs::pitt_peters::{PittPetersModel, PittPetersRotorState};
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::quasi_static_bem::{QuasiStaticBEM, QuasiStaticRotorState};
use dynbem_rs::rotor_definition::{
    BladeGeometry, ControlProperties, LinearPolarParameters, PitchActuation, RotorDefinition,
};
use dynbem_rs::vpm::{VpmRotor, VpmRotorConfig, VpmRotorState};

const RHO: f64 = 1.225;
const G: f64 = 9.81;

// Hardcoded requested setup
const MASS_KG: f64 = 5.0;
const TETHER_LENGTH_M: f64 = 60.0;
const OMEGA_RAD_S: f64 = 20.148;
const TETHER_FORCE_N: f64 = MASS_KG * G;

pub const SIM_DT: f64 = 0.02;

const ELEV_SETPOINT_DEG: f64 = 35.0;

const COLLECTIVE_MIN_RAD: f64 = -0.35;
const COLLECTIVE_MAX_RAD: f64 = 0.65;
const CYCLIC_MAX_RAD: f64 = 0.35;
const MAX_ACCEL_MPS2: f64 = 12.0;
const MAX_ANG_RATE_RAD_S: f64 = 0.35;

#[derive(Clone, Copy)]
struct Pid {
    kp: f64,
    ki: f64,
    kd: f64,
    integ: f64,
    prev_err: f64,
    integ_min: f64,
    integ_max: f64,
}

impl Pid {
    fn new(kp: f64, ki: f64, kd: f64, integ_min: f64, integ_max: f64) -> Self {
        Self {
            kp,
            ki,
            kd,
            integ: 0.0,
            prev_err: 0.0,
            integ_min,
            integ_max,
        }
    }

    fn update(&mut self, err: f64, dt: f64) -> f64 {
        self.integ = (self.integ + err * dt).clamp(self.integ_min, self.integ_max);
        let derr = (err - self.prev_err) / dt;
        self.prev_err = err;
        self.kp * err + self.ki * self.integ + self.kd * derr
    }
}

#[derive(Clone, Copy, Debug)]
pub enum ModelChoice {
    Bem,
    PittPeters,
    Oye,
    Vpm,
}

impl ModelChoice {
    pub fn from_arg(s: &str) -> Option<Self> {
        match s.to_ascii_lowercase().as_str() {
            "bem" => Some(Self::Bem),
            "pitt_peters" | "pitt-peters" | "pp" => Some(Self::PittPeters),
            "oye" => Some(Self::Oye),
            "vpm" => Some(Self::Vpm),
            _ => None,
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Bem => "bem",
            Self::PittPeters => "pitt_peters",
            Self::Oye => "oye",
            Self::Vpm => "vpm",
        }
    }
}

#[derive(Resource)]
pub struct SimRunner {
    pub sim: FlightSim,
    pub accumulator_s: f64,
}

enum AeroRunner {
    Bem {
        model: QuasiStaticBEM<LinearPolar>,
        state: QuasiStaticRotorState,
    },
    PittPeters {
        model: PittPetersModel<LinearPolar>,
        state: PittPetersRotorState,
    },
    Oye {
        model: OyeBEMModel<LinearPolar>,
        state: OyeRotorState,
    },
    Vpm {
        model: VpmRotor<LinearPolar>,
        state: VpmRotorState,
    },
}

impl AeroRunner {
    fn new(choice: ModelChoice, defn: &RotorDefinition, polar: &LinearPolar) -> Self {
        let n_psi = 36;
        match choice {
            ModelChoice::Bem => {
                let model = QuasiStaticBEM::build(defn.clone(), n_psi, polar.clone());
                let state = model.initial_state();
                Self::Bem { model, state }
            }
            ModelChoice::PittPeters => {
                let model = PittPetersModel::build(defn.clone(), n_psi, polar.clone());
                let state = model.initial_state();
                Self::PittPeters { model, state }
            }
            ModelChoice::Oye => {
                let model = OyeBEMModel::build(defn.clone(), n_psi, polar.clone());
                let state = model.initial_state();
                Self::Oye { model, state }
            }
            ModelChoice::Vpm => {
                let config = VpmRotorConfig {
                    max_particles: 6000,
                    sigma: 0.30,
                    barnes_hut: true,
                    bh_theta: 0.5,
                    bh_min_particles: 200,
                    ..VpmRotorConfig::default()
                };
                let model = VpmRotor::new(defn, polar.clone(), ControlGains::default(), config);
                let state = model.initial_state();
                Self::Vpm { model, state }
            }
        }
    }

    fn step(&mut self, inputs: &RotorInputs, dt: f64) -> AeroResult {
        match self {
            Self::Bem { model, state } => {
                let (res, next) =
                    model.step(inputs, state, dt, IntegrationMethod::SemiImplicitEuler);
                *state = next;
                res
            }
            Self::PittPeters { model, state } => {
                let (res, next) =
                    model.step(inputs, state, dt, IntegrationMethod::SemiImplicitEuler);
                *state = next;
                res
            }
            Self::Oye { model, state } => {
                let (res, next) =
                    model.step(inputs, state, dt, IntegrationMethod::SemiImplicitEuler);
                *state = next;
                res
            }
            Self::Vpm { model, state } => {
                let (res, next) =
                    model.step(inputs, state, dt, IntegrationMethod::SemiImplicitEuler);
                *state = next;
                res
            }
        }
    }

    /// Free-wake particle positions in the hub frame (metres) paired with their
    /// vortex strength magnitude |alpha| (m^3/s), or empty for the non-VPM
    /// models which carry no explicit particle cloud.
    fn wake_particles_hub(&self) -> Vec<([f32; 3], f32)> {
        match self {
            Self::Vpm { state, .. } => match &state.wake {
                Some(field) => field
                    .particles()
                    .map(|(pos, strength, _sigma)| {
                        let mag = (strength[0] * strength[0]
                            + strength[1] * strength[1]
                            + strength[2] * strength[2])
                            .sqrt();
                        (pos, mag)
                    })
                    .collect(),
                None => Vec::new(),
            },
            _ => Vec::new(),
        }
    }
}

fn beaupoil_rotor_definition() -> RotorDefinition {
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
        description: "Hardcoded Beaupoil rotor for steady-flight binary".to_string(),
    }
}

fn rot_x(roll: f64) -> Mat3 {
    let (s, c) = roll.sin_cos();
    Mat3([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])
}

fn rot_y(pitch: f64) -> Mat3 {
    let (s, c) = pitch.sin_cos();
    Mat3([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])
}

fn mat_mul(a: Mat3, b: Mat3) -> Mat3 {
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = a.0[i][0] * b.0[0][j] + a.0[i][1] * b.0[1][j] + a.0[i][2] * b.0[2][j];
        }
    }
    Mat3(out)
}

pub fn hub_rotation(roll: f64, pitch: f64) -> Mat3 {
    // Y then X: pitch then roll for a simple level-referenced orientation.
    mat_mul(rot_y(pitch), rot_x(roll))
}

fn vec_norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn desired_roll_pitch_from_force_balance(f_ext: [f64; 3]) -> (f64, f64) {
    // Rotor thrust points opposite hub +Z axis, so to balance external force
    // (gravity + tether), we align hub +Z with the external-force direction.
    let mag = vec_norm(f_ext).max(1e-9);
    let axis = [f_ext[0] / mag, f_ext[1] / mag, f_ext[2] / mag];

    // For R = Ry(pitch)*Rx(roll), hub +Z axis in world is:
    // [sin(pitch)cos(roll), -sin(roll), cos(pitch)cos(roll)]
    let roll = (-axis[1]).asin();
    let pitch = axis[0].atan2(axis[2]);
    (roll, pitch)
}

pub struct FlightSim {
    aero: AeroRunner,
    pos: [f64; 3],
    vel: [f64; 3],
    roll: f64,
    pitch: f64,
    p_rate: f64,
    q_rate: f64,
    pid_elev: Pid,
    pid_roll: Pid,
    pid_pitch: Pid,
    rotor_phase: f64,
    collective: f64,
    tilt_lon: f64,
    tilt_lat: f64,
}

impl FlightSim {
    pub fn new(choice: ModelChoice) -> Self {
        let defn = beaupoil_rotor_definition();
        let polar = LinearPolar::from_properties(&defn.airfoil);
        let aero = AeroRunner::new(choice, &defn, &polar);
        let target_elev = ELEV_SETPOINT_DEG.to_radians();
        Self {
            aero,
            pos: [
                TETHER_LENGTH_M * (target_elev - 0.08).cos(),
                0.0,
                -TETHER_LENGTH_M * (target_elev - 0.08).sin(),
            ],
            vel: [0.0, 0.0, 0.0],
            roll: 4.0_f64.to_radians(),
            pitch: -3.0_f64.to_radians(),
            p_rate: 0.0,
            q_rate: 0.0,
            pid_elev: Pid::new(0.35, 0.08, 0.12, -0.5, 0.5),
            pid_roll: Pid::new(1.2, 0.05, 0.25, -0.3, 0.3),
            pid_pitch: Pid::new(1.2, 0.05, 0.25, -0.3, 0.3),
            rotor_phase: 0.0,
            collective: 0.0,
            tilt_lon: 0.0,
            tilt_lat: 0.0,
        }
    }

    pub fn position_ned(&self) -> [f64; 3] {
        self.pos
    }

    pub fn roll_rad(&self) -> f64 {
        self.roll
    }

    pub fn pitch_rad(&self) -> f64 {
        self.pitch
    }

    pub fn rotor_phase_rad(&self) -> f64 {
        self.rotor_phase
    }

    /// Free-wake particle positions in the hub frame (metres) paired with vortex
    /// strength magnitude |alpha|. Empty unless the active aero model is VPM.
    pub fn wake_particles_hub(&self) -> Vec<([f32; 3], f32)> {
        self.aero.wake_particles_hub()
    }

    pub fn step(&mut self, dt: f64) {
        let target_elev = ELEV_SETPOINT_DEG.to_radians();
        let collective_trim = -0.03;

        let horiz = (self.pos[0] * self.pos[0] + self.pos[1] * self.pos[1])
            .sqrt()
            .max(1e-9);
        let elev = (-self.pos[2]).atan2(horiz);
        let elev_err = target_elev - elev;

        let col_cmd = self.pid_elev.update(elev_err, dt);
        self.collective = (collective_trim + col_cmd).clamp(COLLECTIVE_MIN_RAD, COLLECTIVE_MAX_RAD);

        let r_mag = vec_norm(self.pos).max(1e-9);
        let tether_dir = [
            -self.pos[0] / r_mag,
            -self.pos[1] / r_mag,
            -self.pos[2] / r_mag,
        ];
        let f_tether = [
            TETHER_FORCE_N * tether_dir[0],
            TETHER_FORCE_N * tether_dir[1],
            TETHER_FORCE_N * tether_dir[2],
        ];
        let f_grav = [0.0, 0.0, MASS_KG * G];
        let f_ext = [
            f_grav[0] + f_tether[0],
            f_grav[1] + f_tether[1],
            f_grav[2] + f_tether[2],
        ];
        let (roll_target, pitch_target) = desired_roll_pitch_from_force_balance(f_ext);

        let roll_err = roll_target - self.roll;
        let pitch_err = pitch_target - self.pitch;
        self.tilt_lon =
            (-self.pid_pitch.update(pitch_err, dt)).clamp(-CYCLIC_MAX_RAD, CYCLIC_MAX_RAD);
        self.tilt_lat = self
            .pid_roll
            .update(roll_err, dt)
            .clamp(-CYCLIC_MAX_RAD, CYCLIC_MAX_RAD);

        let inputs = RotorInputs {
            collective_rad: self.collective,
            tilt_lon: self.tilt_lon,
            tilt_lat: self.tilt_lat,
            R_hub: hub_rotation(self.roll, self.pitch),
            v_hub_world: DVec3::new(self.vel[0], self.vel[1], self.vel[2]),
            wind_world: DVec3::zero(),
            rho_kg_m3: RHO,
            omega_rad_s: OMEGA_RAD_S,
        };

        let res = self.aero.step(&inputs, dt);
        let f_aero = res.F_world.0;
        let f_total = [
            f_aero[0] + f_ext[0],
            f_aero[1] + f_ext[1],
            f_aero[2] + f_ext[2],
        ];

        let translational_damping = 0.60;
        let f_damp = [
            -translational_damping * self.vel[0],
            -translational_damping * self.vel[1],
            -translational_damping * self.vel[2],
        ];

        let acc = [
            (f_total[0] + f_damp[0]) / MASS_KG,
            (f_total[1] + f_damp[1]) / MASS_KG,
            (f_total[2] + f_damp[2]) / MASS_KG,
        ];

        let acc_mag = vec_norm(acc);
        let acc_limited = if acc_mag > MAX_ACCEL_MPS2 {
            let s = MAX_ACCEL_MPS2 / acc_mag;
            [acc[0] * s, acc[1] * s, acc[2] * s]
        } else {
            acc
        };

        self.vel[0] += acc_limited[0] * dt;
        self.vel[1] += acc_limited[1] * dt;
        self.vel[2] += acc_limited[2] * dt;

        self.pos[0] += self.vel[0] * dt;
        self.pos[1] += self.vel[1] * dt;
        self.pos[2] += self.vel[2] * dt;

        let mx = res.M_hub_world.0[0];
        let my = res.M_hub_world.0[1];
        let i_roll = 8.0;
        let i_pitch = 8.0;
        let ang_damp = 1.5;
        let p_dot = (mx - ang_damp * self.p_rate) / i_roll;
        let q_dot = (my - ang_damp * self.q_rate) / i_pitch;
        self.p_rate += p_dot * dt;
        self.q_rate += q_dot * dt;
        self.p_rate = self.p_rate.clamp(-MAX_ANG_RATE_RAD_S, MAX_ANG_RATE_RAD_S);
        self.q_rate = self.q_rate.clamp(-MAX_ANG_RATE_RAD_S, MAX_ANG_RATE_RAD_S);
        self.roll += self.p_rate * dt;
        self.pitch += self.q_rate * dt;

        self.rotor_phase += OMEGA_RAD_S * dt;
    }
}

pub fn update_simulation(time: Res<Time>, mut runner: ResMut<SimRunner>) {
    // Cap how much wall time we try to catch up in one frame. VPM steps are
    // expensive, so without this a slow frame would queue several sub-steps and
    // spiral (each frame falling further behind). Dropping excess accumulated
    // time keeps the loop stable at the cost of slight slow-motion under load.
    const MAX_STEPS_PER_FRAME: u32 = 3;

    runner.accumulator_s += time.delta_seconds_f64();
    let mut steps = 0;
    while runner.accumulator_s >= SIM_DT && steps < MAX_STEPS_PER_FRAME {
        runner.sim.step(SIM_DT);
        runner.accumulator_s -= SIM_DT;
        steps += 1;
    }
    if runner.accumulator_s > SIM_DT {
        // Discard the backlog we could not process this frame.
        runner.accumulator_s = 0.0;
    }
}
