// Rotor definition: blade geometry, airfoil, control, inertia, autorotation.
// Pure data structs. YAML loading stays Python-side. See ../../../CLAUDE.md.

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

#[derive(Clone, Debug, Default)]
#[allow(non_snake_case)]
pub struct KamanFlap {
    pub enabled: bool,
    pub chord_fraction: Option<f64>,
    pub span_start_m: Option<f64>,
    pub span_end_m: Option<f64>,
    pub tau: Option<f64>,
    pub CM_gamma_per_rad: Option<f64>,
    pub swashplate_load_fraction: Option<f64>,
    pub notes: String,
}

#[derive(Clone, Debug)]
pub struct BladeGeometry {
    pub n_blades: usize,
    pub radius_m: f64,
    pub root_cutout_m: f64,
    pub chord_m: f64,
    pub twist_deg: f64,
    pub n_elements: usize,
    pub r_stations_m: Vec<f64>,
    pub chord_stations_m: Vec<f64>,
    pub twist_stations_deg: Vec<f64>,
}

impl BladeGeometry {
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
pub struct AirfoilProperties {
    pub Re_design: i64,
    pub CL0: f64,
    pub CL_alpha_per_rad: f64,
    pub CD0: f64,
    pub alpha_stall_deg: f64,
    pub tip_loss: bool,
    pub name: String,
    pub source: String,
    pub polar_csv: Option<String>,
    pub CD_structural: f64,
    pub Re_operating: Option<i64>,
}

#[derive(Clone, Debug, Default)]
#[allow(non_snake_case)]
pub struct InertiaProperties {
    pub mass_kg: Option<f64>,
    pub I_body_kgm2: Vec<f64>,
    pub I_spin_kgm2: Option<f64>,
    pub blade_mass_kg: Option<f64>,
    pub stationary_assembly_mass_kg: Option<f64>,
    pub spinning_hub_shell_mass_kg: Option<f64>,
    pub I_blade_flap_kgm2: Option<f64>,
}

#[derive(Clone, Debug)]
#[allow(non_snake_case)]
pub struct ControlProperties {
    pub swashplate_pitch_gain_rad: f64,
    pub axle_attachment_length_m: Option<f64>,
    pub K_cyc: Option<f64>,
    pub swashplate_phase_deg: Option<f64>,
    pub servo_slew_rate_deg_s: Option<f64>,
    pub servo_travel_deg: Option<f64>,
    pub kaman_flap: KamanFlap,
}

#[derive(Clone, Debug, Default)]
#[allow(non_snake_case)]
pub struct AutorotationProperties {
    pub I_ode_kgm2: Option<f64>,
    pub omega_min_rad_s: Option<f64>,
    pub omega_eq_rad_s: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct RotorDefinition {
    pub blade: BladeGeometry,
    pub airfoil: AirfoilProperties,
    pub control: Option<ControlProperties>,
    pub inertia: InertiaProperties,
    pub autorotation: AutorotationProperties,
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
}
