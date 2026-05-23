// YAML loader for RotorDefinition. Pure Rust (no pyo3 / no numpy);
// the only file IO in dynbem_rs lives here. Mirrors the schema the
// legacy Python loader (dynbem.rotor_definition.load) accepts, so
// every rotor.yaml under rotors/* round-trips through both loaders.

use std::fs;
use std::path::{Path, PathBuf};

use serde::Deserialize;

use crate::rotor_definition::{
    AirfoilProperties, AutorotationProperties, BladeGeometry, ControlProperties,
    InertiaProperties, KamanFlap, RotorDefinition,
};

#[derive(Debug)]
pub enum YamlLoadError {
    Io(std::io::Error),
    Parse(serde_yaml::Error),
    /// A required field was missing (path = dotted YAML path,
    /// e.g. "rotor.n_blades").
    MissingField(String),
    /// A value had the wrong shape (e.g. mass_kg dict with a
    /// non-numeric value).
    InvalidValue { field: String, reason: String },
}

impl std::fmt::Display for YamlLoadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            YamlLoadError::Io(e) => write!(f, "I/O error: {}", e),
            YamlLoadError::Parse(e) => write!(f, "YAML parse error: {}", e),
            YamlLoadError::MissingField(p) => write!(f, "missing required field: {}", p),
            YamlLoadError::InvalidValue { field, reason } => {
                write!(f, "invalid value for {}: {}", field, reason)
            }
        }
    }
}

impl std::error::Error for YamlLoadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            YamlLoadError::Io(e) => Some(e),
            YamlLoadError::Parse(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for YamlLoadError {
    fn from(e: std::io::Error) -> Self {
        YamlLoadError::Io(e)
    }
}

impl From<serde_yaml::Error> for YamlLoadError {
    fn from(e: serde_yaml::Error) -> Self {
        YamlLoadError::Parse(e)
    }
}

// ---------------------------------------------------------------------------
// YAML mirror structs. These match the on-disk schema exactly; the
// Rust-side RotorDefinition is built from them below. Keeping the
// mirror separate from the public structs keeps RotorDefinition free
// of serde derives (and free of file-format concerns).
//
// Intentionally no `deny_unknown_fields`: the legacy Python loader was
// permissive (dict.get on every key, silently ignoring extras), so
// existing YAMLs with notes / commented-out sections / future fields
// keep working when round-tripped through Rust.
// ---------------------------------------------------------------------------

/// `mass_kg` in the YAML may be a number OR a dict whose values are
/// summed (the legacy loader interpreted dicts as a sum of named
/// components, so existing YAMLs that use that shape keep working).
#[derive(Debug, Deserialize)]
#[serde(untagged)]
enum MaybeSum {
    Number(f64),
    Sum(std::collections::BTreeMap<String, f64>),
}

impl MaybeSum {
    fn to_f64(&self) -> f64 {
        match self {
            MaybeSum::Number(v) => *v,
            MaybeSum::Sum(m) => m.values().sum(),
        }
    }
}

#[derive(Debug, Default, Deserialize)]
struct RotorYaml {
    #[serde(default)]
    name: Option<String>,
    #[serde(default)]
    description: Option<String>,
    rotor: BladeYaml,
    airfoil: AirfoilYaml,
    #[serde(default)]
    inertia: InertiaYaml,
    #[serde(default)]
    control: Option<ControlYaml>,
    #[serde(default)]
    autorotation: AutorotationYaml,
}

#[derive(Debug, Default, Deserialize)]
struct BladeYaml {
    n_blades: Option<usize>,
    radius_m: Option<f64>,
    root_cutout_m: Option<f64>,
    chord_m: Option<f64>,
    #[serde(default)]
    twist_deg: Option<f64>,
    #[serde(default)]
    n_elements: Option<usize>,
    #[serde(default)]
    r_stations_m: Vec<f64>,
    #[serde(default)]
    chord_stations_m: Vec<f64>,
    #[serde(default)]
    twist_stations_deg: Vec<f64>,
}

#[allow(non_snake_case)]
#[derive(Debug, Default, Deserialize)]
struct AirfoilYaml {
    Re_design: Option<i64>,
    CL0: Option<f64>,
    CL_alpha_per_rad: Option<f64>,
    CD0: Option<f64>,
    alpha_stall_deg: Option<f64>,
    #[serde(default)]
    tip_loss: Option<bool>,
    #[serde(default)]
    designation: Option<String>,
    #[serde(default)]
    source: Option<String>,
    #[serde(default)]
    polar_csv: Option<String>,
    #[serde(default)]
    CD_structural: Option<f64>,
    #[serde(default)]
    Re_operating: Option<i64>,
}

#[allow(non_snake_case)]
#[derive(Debug, Default, Deserialize)]
struct InertiaYaml {
    #[serde(default)]
    mass_kg: Option<MaybeSum>,
    #[serde(default)]
    I_body_kgm2: Vec<f64>,
    #[serde(default)]
    I_spin_kgm2: Option<MaybeSum>,
    #[serde(default)]
    blade_mass_kg: Option<MaybeSum>,
    #[serde(default)]
    stationary_assembly_mass_kg: Option<MaybeSum>,
    #[serde(default)]
    spinning_hub_shell_mass_kg: Option<MaybeSum>,
    #[serde(default)]
    I_blade_flap_kgm2: Option<MaybeSum>,
}

#[allow(non_snake_case)]
#[derive(Debug, Default, Deserialize)]
struct ControlYaml {
    swashplate_pitch_gain_rad: Option<f64>,
    #[serde(default)]
    axle_attachment_length_m: Option<f64>,
    #[serde(default)]
    K_cyc: Option<f64>,
    #[serde(default)]
    swashplate_phase_deg: Option<f64>,
    #[serde(default)]
    servo_slew_rate_deg_s: Option<f64>,
    #[serde(default)]
    servo_travel_deg: Option<f64>,
    #[serde(default)]
    kaman_flap: Option<KamanFlapYaml>,
}

#[allow(non_snake_case)]
#[derive(Debug, Default, Deserialize)]
struct KamanFlapYaml {
    #[serde(default)]
    enabled: Option<bool>,
    #[serde(default)]
    chord_fraction: Option<f64>,
    #[serde(default)]
    span_start_m: Option<f64>,
    #[serde(default)]
    span_end_m: Option<f64>,
    #[serde(default)]
    tau: Option<f64>,
    #[serde(default)]
    CM_gamma_per_rad: Option<f64>,
    #[serde(default)]
    swashplate_load_fraction: Option<f64>,
    #[serde(default)]
    notes: Option<String>,
}

#[allow(non_snake_case)]
#[derive(Debug, Default, Deserialize)]
struct AutorotationYaml {
    #[serde(default)]
    I_ode_kgm2: Option<f64>,
    #[serde(default)]
    omega_min_rad_s: Option<f64>,
    #[serde(default)]
    omega_eq_rad_s: Option<f64>,
}

// ---------------------------------------------------------------------------
// Public API.
// ---------------------------------------------------------------------------

/// Parse a RotorDefinition from a YAML string.
///
/// `base_dir`, if `Some`, is used to resolve a relative `airfoil.polar_csv`
/// path; if `None`, the field is stored verbatim (callers who already
/// have an absolute path or who want raw input can pass `None`).
pub fn from_yaml_str(text: &str, base_dir: Option<&Path>) -> Result<RotorDefinition, YamlLoadError> {
    let parsed: RotorYaml = serde_yaml::from_str(text)?;
    build(parsed, base_dir)
}

/// Read a YAML file from disk and parse it into a RotorDefinition.
/// `airfoil.polar_csv`, if present and relative, is resolved against
/// the YAML file's parent directory (matching the legacy loader).
pub fn from_yaml_file(path: impl AsRef<Path>) -> Result<RotorDefinition, YamlLoadError> {
    let path = path.as_ref();
    let text = fs::read_to_string(path)?;
    let base = path.parent();
    from_yaml_str(&text, base)
}

impl RotorDefinition {
    /// Convenience wrapper around [`from_yaml_str`].
    pub fn from_yaml_str(
        text: &str,
        base_dir: Option<&Path>,
    ) -> Result<RotorDefinition, YamlLoadError> {
        from_yaml_str(text, base_dir)
    }

    /// Convenience wrapper around [`from_yaml_file`].
    pub fn from_yaml_file(path: impl AsRef<Path>) -> Result<RotorDefinition, YamlLoadError> {
        from_yaml_file(path)
    }
}

// ---------------------------------------------------------------------------
// Build helpers.
// ---------------------------------------------------------------------------

fn req<T>(value: Option<T>, field: &str) -> Result<T, YamlLoadError> {
    value.ok_or_else(|| YamlLoadError::MissingField(field.to_string()))
}

fn resolve_polar_csv(value: Option<String>, base_dir: Option<&Path>) -> Option<String> {
    let v = value?;
    let p = PathBuf::from(&v);
    if p.is_absolute() || base_dir.is_none() {
        return Some(v);
    }
    Some(base_dir.unwrap().join(p).to_string_lossy().into_owned())
}

fn resolve_mass_kg(inertia: &InertiaYaml, n_blades: usize) -> Option<f64> {
    if let Some(explicit) = &inertia.mass_kg {
        return Some(explicit.to_f64());
    }
    // Derive from components: blade_mass * N_blades + stationary + shell.
    let blade = inertia.blade_mass_kg.as_ref().map(|x| x.to_f64());
    let stat = inertia.stationary_assembly_mass_kg.as_ref().map(|x| x.to_f64());
    let shell = inertia.spinning_hub_shell_mass_kg.as_ref().map(|x| x.to_f64());
    match (blade, stat, shell) {
        (Some(b), Some(s), Some(h)) => Some(b * n_blades as f64 + s + h),
        _ => None,
    }
}

fn build(y: RotorYaml, base_dir: Option<&Path>) -> Result<RotorDefinition, YamlLoadError> {
    let r = &y.rotor;
    let n_blades = req(r.n_blades, "rotor.n_blades")?;
    let radius_m = req(r.radius_m, "rotor.radius_m")?;
    let root_cutout_m = req(r.root_cutout_m, "rotor.root_cutout_m")?;
    let chord_m = req(r.chord_m, "rotor.chord_m")?;

    let blade = BladeGeometry {
        n_blades,
        radius_m,
        root_cutout_m,
        chord_m,
        twist_deg: r.twist_deg.unwrap_or(0.0),
        n_elements: r.n_elements.unwrap_or(10),
        r_stations_m: r.r_stations_m.clone(),
        chord_stations_m: r.chord_stations_m.clone(),
        twist_stations_deg: r.twist_stations_deg.clone(),
    };

    let a = &y.airfoil;
    let airfoil = AirfoilProperties {
        Re_design: req(a.Re_design, "airfoil.Re_design")?,
        CL0: req(a.CL0, "airfoil.CL0")?,
        CL_alpha_per_rad: req(a.CL_alpha_per_rad, "airfoil.CL_alpha_per_rad")?,
        CD0: req(a.CD0, "airfoil.CD0")?,
        alpha_stall_deg: req(a.alpha_stall_deg, "airfoil.alpha_stall_deg")?,
        tip_loss: a.tip_loss.unwrap_or(true),
        name: a.designation.clone().unwrap_or_default(),
        source: a.source.clone().unwrap_or_default(),
        polar_csv: resolve_polar_csv(a.polar_csv.clone(), base_dir),
        CD_structural: a.CD_structural.unwrap_or(0.0),
        Re_operating: a.Re_operating,
    };

    let control = match y.control {
        None => None,
        Some(c) => Some(ControlProperties {
            swashplate_pitch_gain_rad: req(
                c.swashplate_pitch_gain_rad,
                "control.swashplate_pitch_gain_rad",
            )?,
            axle_attachment_length_m: c.axle_attachment_length_m,
            K_cyc: c.K_cyc,
            swashplate_phase_deg: c.swashplate_phase_deg,
            servo_slew_rate_deg_s: c.servo_slew_rate_deg_s,
            servo_travel_deg: c.servo_travel_deg,
            kaman_flap: match c.kaman_flap {
                None => KamanFlap::default(),
                Some(k) => KamanFlap {
                    enabled: k.enabled.unwrap_or(false),
                    chord_fraction: k.chord_fraction,
                    span_start_m: k.span_start_m,
                    span_end_m: k.span_end_m,
                    tau: k.tau,
                    CM_gamma_per_rad: k.CM_gamma_per_rad,
                    swashplate_load_fraction: k.swashplate_load_fraction,
                    notes: k.notes.unwrap_or_default(),
                },
            },
        }),
    };

    let inertia = InertiaProperties {
        mass_kg: resolve_mass_kg(&y.inertia, n_blades),
        I_body_kgm2: y.inertia.I_body_kgm2.clone(),
        I_spin_kgm2: y.inertia.I_spin_kgm2.as_ref().map(|x| x.to_f64()),
        blade_mass_kg: y.inertia.blade_mass_kg.as_ref().map(|x| x.to_f64()),
        stationary_assembly_mass_kg: y.inertia
            .stationary_assembly_mass_kg
            .as_ref()
            .map(|x| x.to_f64()),
        spinning_hub_shell_mass_kg: y.inertia
            .spinning_hub_shell_mass_kg
            .as_ref()
            .map(|x| x.to_f64()),
        I_blade_flap_kgm2: y.inertia.I_blade_flap_kgm2.as_ref().map(|x| x.to_f64()),
    };

    let autorotation = AutorotationProperties {
        I_ode_kgm2: y.autorotation.I_ode_kgm2,
        omega_min_rad_s: y.autorotation.omega_min_rad_s,
        omega_eq_rad_s: y.autorotation.omega_eq_rad_s,
    };

    Ok(RotorDefinition {
        blade,
        airfoil,
        control,
        inertia,
        autorotation,
        name: y.name.unwrap_or_default(),
        description: y.description.unwrap_or_default(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const CASTLES_GRAY: &str = r#"
name: Castles-Gray-6ft
description: minimal fixture
rotor:
  n_blades: 3
  radius_m: 0.914
  root_cutout_m: 0.155
  chord_m: 0.0479
  twist_deg: 0.0
  n_elements: 30
airfoil:
  designation: NACA 0015
  source: XFOIL
  Re_design: 256000
  CL0: 0.0
  CL_alpha_per_rad: 5.90
  CD0: 0.01046
  alpha_stall_deg: 15.5
  tip_loss: true
  polar_csv: naca0015_ncrit5_re200k.csv
autorotation:
  I_ode_kgm2: 1.0
"#;

    #[test]
    fn parses_castles_gray_minimal() {
        let defn = from_yaml_str(CASTLES_GRAY, None).expect("parse");
        assert_eq!(defn.blade.n_blades, 3);
        assert!((defn.blade.radius_m - 0.914).abs() < 1e-12);
        assert!((defn.blade.chord_m - 0.0479).abs() < 1e-12);
        assert_eq!(defn.blade.n_elements, 30);
        assert!(!defn.blade.has_radial_stations());
        assert_eq!(defn.airfoil.name, "NACA 0015");
        assert_eq!(defn.airfoil.polar_csv.as_deref(), Some("naca0015_ncrit5_re200k.csv"));
        assert!(defn.control.is_none());
        assert_eq!(defn.autorotation.I_ode_kgm2, Some(1.0));
    }

    #[test]
    fn polar_csv_resolves_against_base_dir() {
        let base = Path::new("/tmp/rotors/foo");
        let defn = from_yaml_str(CASTLES_GRAY, Some(base)).expect("parse");
        let csv = defn.airfoil.polar_csv.unwrap();
        // Resolved path must contain the base dir + the relative file name.
        assert!(csv.contains("naca0015_ncrit5_re200k.csv"));
        assert!(csv.contains("foo"));
    }

    const BEAUPOIL: &str = r#"
name: beaupoil_2026
rotor:
  n_blades: 4
  radius_m: 2.5
  root_cutout_m: 0.5
  chord_m: 0.20
  twist_deg: 0.0
  n_elements: 10
airfoil:
  designation: SG6040
  Re_design: 127000
  CL0: 0.393
  CL_alpha_per_rad: 5.79
  CD0: 0.0079
  alpha_stall_deg: 13.0
  tip_loss: true
inertia:
  mass_kg: 5.00
  I_body_kgm2: [5.0, 5.0, 10.0]
  I_spin_kgm2: 3.94
  blade_mass_kg: 0.35
  stationary_assembly_mass_kg: 0.99
  spinning_hub_shell_mass_kg: 2.61
control:
  swashplate_pitch_gain_rad: 0.3
  axle_attachment_length_m:  0.3
  K_cyc:                     0.4
  swashplate_phase_deg:      0.0
  servo_slew_rate_deg_s:     545.0
  servo_travel_deg:          100.0
  kaman_flap:
    enabled:                  true
    chord_fraction:           0.25
    span_start_m:             1.2
    span_end_m:               2.5
    tau:                      0.45
    CM_gamma_per_rad:        -0.35
    swashplate_load_fraction: 0.1
autorotation:
  I_ode_kgm2: 10.0
  omega_min_rad_s: 0.5
  omega_eq_rad_s: 20.148
"#;

    #[test]
    fn parses_beaupoil_full() {
        let defn = from_yaml_str(BEAUPOIL, None).expect("parse");
        assert_eq!(defn.blade.n_blades, 4);
        let c = defn.control.as_ref().expect("control present");
        assert!((c.swashplate_pitch_gain_rad - 0.3).abs() < 1e-12);
        let k = &c.kaman_flap;
        assert!(k.enabled);
        assert_eq!(k.chord_fraction, Some(0.25));
        assert!((defn.inertia.mass_kg.unwrap() - 5.00).abs() < 1e-12);
        assert_eq!(defn.inertia.I_body_kgm2.len(), 3);
        assert_eq!(defn.autorotation.omega_eq_rad_s, Some(20.148));
    }

    #[test]
    fn mass_kg_derived_from_components_when_omitted() {
        let yaml = r#"
rotor:
  n_blades: 4
  radius_m: 2.5
  root_cutout_m: 0.5
  chord_m: 0.20
airfoil:
  Re_design: 1
  CL0: 0.0
  CL_alpha_per_rad: 6.0
  CD0: 0.01
  alpha_stall_deg: 12.0
inertia:
  blade_mass_kg: 0.35
  stationary_assembly_mass_kg: 0.99
  spinning_hub_shell_mass_kg: 2.61
"#;
        let defn = from_yaml_str(yaml, None).expect("parse");
        let expected = 0.35 * 4.0 + 0.99 + 2.61;
        assert!((defn.inertia.mass_kg.unwrap() - expected).abs() < 1e-12);
    }

    #[test]
    fn mass_kg_dict_sums_values() {
        let yaml = r#"
rotor:
  n_blades: 2
  radius_m: 1.0
  root_cutout_m: 0.1
  chord_m: 0.05
airfoil:
  Re_design: 1
  CL0: 0.0
  CL_alpha_per_rad: 6.0
  CD0: 0.01
  alpha_stall_deg: 12.0
inertia:
  mass_kg:
    foo: 1.5
    bar: 2.0
    baz: 0.25
"#;
        let defn = from_yaml_str(yaml, None).expect("parse");
        assert!((defn.inertia.mass_kg.unwrap() - 3.75).abs() < 1e-12);
    }

    #[test]
    fn missing_required_field_errors() {
        let yaml = r#"
rotor:
  n_blades: 2
  radius_m: 1.0
  root_cutout_m: 0.1
  # chord_m omitted
airfoil:
  Re_design: 1
  CL0: 0.0
  CL_alpha_per_rad: 6.0
  CD0: 0.01
  alpha_stall_deg: 12.0
"#;
        let err = from_yaml_str(yaml, None).unwrap_err();
        match err {
            YamlLoadError::MissingField(p) => assert_eq!(p, "rotor.chord_m"),
            other => panic!("expected MissingField, got {:?}", other),
        }
    }

    #[test]
    fn radial_stations_parse() {
        let yaml = r#"
rotor:
  n_blades: 2
  radius_m: 5.0
  root_cutout_m: 1.0
  chord_m: 0.5
  r_stations_m: [1.0, 2.0, 3.0]
  chord_stations_m: [0.5, 0.4, 0.3]
  twist_stations_deg: [-22.0, -10.0, -1.0]
airfoil:
  Re_design: 1
  CL0: 0.0
  CL_alpha_per_rad: 6.0
  CD0: 0.01
  alpha_stall_deg: 12.0
"#;
        let defn = from_yaml_str(yaml, None).expect("parse");
        assert!(defn.blade.has_radial_stations());
        assert_eq!(defn.blade.r_stations_m.len(), 3);
        assert!((defn.blade.chord_at(2.0) - 0.4).abs() < 1e-12);
    }
}
