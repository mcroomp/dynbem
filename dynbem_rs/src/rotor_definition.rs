// Rotor definition: blade geometry, airfoil, control.
// Pure data structs -- only the fields actually used in Rust math.
// YAML loading and metadata (inertia, KamanFlap, polar_csv, Re_design,
// etc.) stay Python-side. See ../../../AGENTS.md.

use std::f64::consts::PI;

#[inline]
fn lerp_clamped(x: f64, xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len();
    if n == 0 {
        panic!("empty interpolation table");
    }
    if n == 1 || x <= xs[0] {
        return ys[0];
    }
    if x >= xs[n - 1] {
        return ys[n - 1];
    }
    let i = xs.partition_point(|&v| v <= x);
    let i = i.max(1).min(n - 1);
    let t = (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
    ys[i - 1] + t * (ys[i] - ys[i - 1])
}

#[derive(Clone, Debug)]
pub struct BladeGeometry {
    pub n_blades: usize,
    pub radius_m: f64,
    pub root_cutout_m: f64,
    pub chord_m: f64,
    pub twist_deg: f64,
    pub n_elements: usize,
    pub tip_loss: bool,
    pub r_stations_m: Vec<f64>,
    pub chord_stations_m: Vec<f64>,
    pub twist_stations_deg: Vec<f64>,
}

impl BladeGeometry {
    /// Construct a uniform blade (constant chord and twist, no radial stations).
    pub fn uniform(
        n_blades: usize,
        radius_m: f64,
        root_cutout_m: f64,
        chord_m: f64,
        twist_deg: f64,
        n_elements: usize,
    ) -> Self {
        Self {
            n_blades,
            radius_m,
            root_cutout_m,
            chord_m,
            twist_deg,
            n_elements,
            tip_loss: true,
            r_stations_m: vec![],
            chord_stations_m: vec![],
            twist_stations_deg: vec![],
        }
    }

    pub fn span_m(&self) -> f64 {
        self.radius_m - self.root_cutout_m
    }
    pub fn r_cp_m(&self) -> f64 {
        self.root_cutout_m + (2.0 / 3.0) * self.span_m()
    }
    pub fn disk_area_m2(&self) -> f64 {
        PI * (self.radius_m * self.radius_m - self.root_cutout_m * self.root_cutout_m)
    }
    pub fn solidity(&self) -> f64 {
        (self.n_blades as f64) * self.chord_m / (PI * self.radius_m)
    }
    pub fn has_radial_stations(&self) -> bool {
        self.r_stations_m.len() >= 2
            && self.chord_stations_m.len() == self.r_stations_m.len()
            && self.twist_stations_deg.len() == self.r_stations_m.len()
    }
    pub fn chord_at(&self, r: f64) -> f64 {
        if !self.has_radial_stations() {
            self.chord_m
        } else {
            lerp_clamped(r, &self.r_stations_m, &self.chord_stations_m)
        }
    }
    pub fn twist_at(&self, r: f64) -> f64 {
        if !self.has_radial_stations() {
            self.twist_deg
        } else {
            lerp_clamped(r, &self.r_stations_m, &self.twist_stations_deg)
        }
    }
}

#[derive(Clone, Debug)]
#[allow(non_snake_case)]
pub struct LinearPolarParameters {
    pub CL0: f64,
    pub CL_alpha_per_rad: f64,
    pub CD0: f64,
    pub alpha_stall_deg: f64,
}

#[derive(Clone, Debug)]
#[allow(non_snake_case)]
pub struct ControlProperties {
    pub swashplate_pitch_gain_rad: f64,
    pub swashplate_phase_deg: Option<f64>,
}

/// Blade feathering (pitch-bearing) DOF driven by a trailing-edge servo-flap.
///
/// The servo-flap exerts a pitching moment about the feathering axis.  Because
/// there is no centrifugal spring (the feathering axis is parallel to the span),
/// the blade feathers freely, damped only by the mechanical damper at the pitch
/// bearing.  The servo-flap moment drives the feathering angle theta(psi).
/// In this (servo-flap) actuation mode both collective and cyclic are
/// interpreted as flap commands and the feathering response REPLACES the direct
/// swashplate pitch path across the full span (see `servoflap.rs`).
///
/// EOM (psi-domain, 1/rev harmonic balance, no aerodynamic spring):
///   I_theta * theta'' + C_theta * theta' = M_servo(psi)
///
/// For the Kaman design the feathering axis is placed at the aerodynamic centre
/// so the aerodynamic restoring moment is zero; set ac_offset_m=0 for this.
#[allow(non_snake_case)]
#[derive(Clone, Debug)]
pub struct ServoFlapActuation {
    /// Blade pitch moment of inertia about the feathering axis [kg*m^2].
    pub I_theta_kgm2: f64,
    /// Rotary damper coefficient at the pitch bearing [N*m*s/rad].
    pub damper_Nms_per_rad: f64,
    /// Distance from feathering axis to aerodynamic centre [m], positive
    /// when AC is forward of the feathering axis (divergent if negative).
    /// 0.0 = feathering axis exactly at AC (Kaman ideal).
    pub ac_offset_m: f64,
    /// Servo-flap geometry that drives the feathering DOF.
    pub flap: ServoFlapGeometry,
}

/// Geometry of the servo-flap (trailing-edge elevon) driven by the swashplate.
///
/// In servo-flap mode the swashplate collective AND cyclic drive flap deflection
/// delta_f(psi) = delta_f0 + delta_f1c*cos(psi) + delta_f1s*sin(psi); delta_f
/// produces an aerodynamic pitching moment about the feathering axis.
#[allow(non_snake_case)]
#[derive(Clone, Debug)]
pub struct ServoFlapGeometry {
    /// Pitching moment coefficient per unit flap deflection [rad^-1].
    /// Thin-airfoil estimate at AC:
    /// C_M_delta = -(1/pi) * (sin(theta_h) - 0.5*sin(2*theta_h)),
    /// cos(theta_h) = 2*(flap_chord_fraction) - 1.
    /// Negative = nose-down moment for downward (positive) flap deflection.
    pub C_M_delta_per_rad: f64,
    /// Inboard edge of servo-flap along blade span [m] (from shaft centre).
    pub r_inner_m: f64,
    /// Outboard edge of servo-flap along blade span [m].
    pub r_outer_m: f64,
}

/// Quasi-static blade flapping properties (hingeless / equivalent-hinge model).
///
/// Models the blade's out-of-plane flexibility as an equivalent spring-hinge.
/// The blade absorbs most of the aerodynamic pitching/rolling moment; only the
/// fraction determined by the flap frequency ratio reaches the hub (airframe).
///
/// For a centrally-hinged blade with zero stiffness: nu_beta = 1.0 exactly,
/// and the hub moment is zero (all absorbed by flapping).
/// For a rigid blade: nu_beta -> infinity, hub moment = full aero moment.
/// For a flexible blade: nu_beta ~ 1.05-1.15 typically.
#[allow(non_snake_case)]
#[derive(Clone, Debug)]
pub struct FlapProperties {
    /// Blade flap moment of inertia about the (virtual) flap hinge [kg*m^2].
    /// For a uniform blade: I_b = m_blade * R^2 / 3.
    pub I_blade_flap_kgm2: f64,
    /// Non-rotating natural frequency of the blade in flap [rad/s].
    /// Related to root bending stiffness: K_beta = I_b * omega_NR^2.
    /// Set to 0.0 for a freely-hinged blade (no spring).
    pub omega_nr_rad_s: f64,
}

impl FlapProperties {
    /// Rotating flap frequency ratio squared: nu_beta^2 = 1 + (omega_NR / Omega)^2.
    /// The "1" comes from centrifugal stiffening at the rotating speed Omega.
    #[inline]
    pub fn nu_beta_sq(&self, omega: f64) -> f64 {
        if omega.abs() < 1e-6 {
            // Non-rotating: no centrifugal stiffening, just the spring.
            // Return large number (rigid limit) to avoid division issues.
            return 1e6;
        }
        1.0 + (self.omega_nr_rad_s / omega).powi(2)
    }

    /// Fraction of the aerodynamic hub moment that passes to the airframe.
    ///
    /// factor = (nu_beta^2 - 1) / nu_beta^2
    ///
    /// - Freely hinged (omega_NR=0): nu^2=1, factor=0 (no moment transfer)
    /// - Rigid blade (omega_NR>>Omega): factor->1 (full moment transfer)
    /// - Typical hingeless: factor ~ 0.05-0.15
    #[inline]
    pub fn hub_moment_factor(&self, omega: f64) -> f64 {
        let nu2 = self.nu_beta_sq(omega);
        (nu2 - 1.0) / nu2
    }
}

/// How blade pitch is actuated.
///
/// - `DirectMechanical`: the swashplate sets blade pitch directly (collective +
///   cyclic map straight to blade-element pitch). This is the default.
/// - `ServoFlap`: a trailing-edge servo-flap drives a passive feathering DOF;
///   the feathering response replaces the direct swashplate pitch path.
#[derive(Clone, Debug)]
pub enum PitchActuation {
    /// Swashplate -> blade pitch directly (rigid pitch).
    DirectMechanical,
    /// Servo-flap moment -> passive feathering DOF.
    ServoFlap(ServoFlapActuation),
}

impl Default for PitchActuation {
    fn default() -> Self {
        PitchActuation::DirectMechanical
    }
}

#[derive(Clone, Debug)]
pub struct RotorDefinition {
    pub blade: BladeGeometry,
    pub airfoil: LinearPolarParameters,
    pub control: Option<ControlProperties>,
    /// Blade pitch actuation mode. Default = `DirectMechanical` (rigid pitch).
    pub pitch_actuation: PitchActuation,
    /// Quasi-static blade flapping. When present, hub moments are reduced by
    /// the flap frequency ratio (blade absorbs most of the moment via deflection).
    pub flap: Option<FlapProperties>,
    pub name: String,
    pub description: String,
}

impl RotorDefinition {
    pub fn span_m(&self) -> f64 {
        self.blade.span_m()
    }
    pub fn r_cp_m(&self) -> f64 {
        self.blade.r_cp_m()
    }
    pub fn disk_area_m2(&self) -> f64 {
        self.blade.disk_area_m2()
    }
    pub fn solidity(&self) -> f64 {
        self.blade.solidity()
    }

    pub fn control_gains(&self) -> crate::cyclic::ControlGains {
        match &self.control {
            None => crate::cyclic::ControlGains::default(),
            Some(c) => {
                let phase = c.swashplate_phase_deg.unwrap_or(0.0).to_radians();
                crate::cyclic::ControlGains {
                    gain: c.swashplate_pitch_gain_rad,
                    phase_rad: phase,
                }
            }
        }
    }

    /// True when blade pitch is actuated by a servo-flap (vs direct mechanical).
    pub fn is_servoflap(&self) -> bool {
        matches!(self.pitch_actuation, PitchActuation::ServoFlap(_))
    }
}
