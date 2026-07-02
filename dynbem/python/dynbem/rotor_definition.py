"""dynbem.rotor_definition submodule.

Defines pure-Python classes for all rotor definition types, including
metadata fields not needed by the Rust math core (Re_design, polar_csv,
inertia, kaman_flap, etc.).  Each Python class internally holds a ``_rust``
attribute -- a lean _dynbem Rust wrapper populated with only the math fields
-- which is what model constructors receive.

YAML loading is handled here in pure Python via yaml.safe_load.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

from ._dynbem import (
    LinearPolarParameters as _RustLinearPolarParameters,
    BladeGeometry as _RustBladeGeometry,
    ControlProperties as _RustControlProperties,
    ServoFlapActuation as _RustServoFlapActuation,
    ServoFlapGeometry as _RustServoFlapGeometry,
    FlapProperties as _RustFlapProperties,
    RotorDefinition as _RustRotorDefinition,
)


# ---------------------------------------------------------------------------
# BladeGeometry: Python subclass with optional station-override arrays.
# All physics dimensions are required; station arrays default to None
# (treated as empty lists = uniform blade).
# ---------------------------------------------------------------------------

class BladeGeometry(_RustBladeGeometry):
    def __new__(
        cls, n_blades, radius_m, root_cutout_m, chord_m,
        twist_deg=0.0, n_elements=20,
        r_stations_m=None, chord_stations_m=None, twist_stations_deg=None,
        tip_loss=True,
    ):
        return _RustBladeGeometry.__new__(
            cls, n_blades, radius_m, root_cutout_m, chord_m,
            twist_deg, n_elements,
            r_stations_m if r_stations_m is not None else [],
            chord_stations_m if chord_stations_m is not None else [],
            twist_stations_deg if twist_stations_deg is not None else [],
            tip_loss,
        )

    def __init__(self, *args, **kwargs):
        pass

__all__ = [
    "LinearPolarParameters",
    "AutorotationProperties",
    "BladeGeometry",
    "ControlProperties",
    "FlapProperties",
    "ServoFlapActuation",
    "ServoFlapGeometry",
    "InertiaProperties",
    "KamanFlap",
    "RotorDefinition",
    "ValidationIssue",
    "load",
    "loads",
    "default",
]


# ---------------------------------------------------------------------------
# Pure-Python metadata classes (no Rust backing needed)
# ---------------------------------------------------------------------------

class KamanFlap:
    """Kaman servo-flap geometry data (Python-only, not used in Rust math)."""

    def __init__(
        self,
        chord_fraction=None,
        span_start_m=None,
        span_end_m=None,
        tau=None,
        CM_gamma_per_rad=None,
        swashplate_load_fraction=None,
        notes="",
    ):
        self.chord_fraction = chord_fraction
        self.span_start_m = span_start_m
        self.span_end_m = span_end_m
        self.tau = tau
        self.CM_gamma_per_rad = CM_gamma_per_rad
        self.swashplate_load_fraction = swashplate_load_fraction
        self.notes = notes


class ServoFlapGeometry:
    """Servo-flap (trailing-edge elevon) geometry.

    In servo-flap mode the swashplate collective AND cyclic drive flap
    deflection delta_f(psi). The flap exerts an aerodynamic pitching moment
    about the feathering axis. Used inside ServoFlapActuation.flap.
    """

    def __init__(
        self,
        C_M_delta_per_rad,
        r_inner_m,
        r_outer_m,
    ):
        self.C_M_delta_per_rad = C_M_delta_per_rad
        self.r_inner_m = r_inner_m
        self.r_outer_m = r_outer_m

    def _to_rust(self):
        return _RustServoFlapGeometry(
            C_M_delta_per_rad=self.C_M_delta_per_rad,
            r_inner_m=self.r_inner_m,
            r_outer_m=self.r_outer_m,
        )


class ServoFlapActuation:
    """Servo-flap pitch actuation: a passive blade feathering DOF.

    A trailing-edge servo-flap exerts a pitching moment about the feathering
    axis.  The blade feathers freely, damped by the bearing damper.  Both
    collective and cyclic are interpreted as flap commands and the resulting
    delta_theta(psi) REPLACES the direct swashplate pitch path in the psi-loop.

    Fields:
      I_theta_kgm2        -- pitch moment of inertia about feathering axis [kg*m^2]
      damper_Nms_per_rad  -- bearing rotary damper coefficient [N*m*s/rad]
      flap                -- ServoFlapGeometry (required)
      ac_offset_m         -- AC distance forward of feathering axis [m] (0 = Kaman ideal)
      control_stiffness_Nm_per_rad -- pushrod/linkage torsional stiffness about
                            the feathering axis [N*m/rad]. Restoring moment that
                            makes the feathering DOF well-posed (feathering has no
                            centrifugal stiffening). Consumed by the time-domain
                            VPM feathering DOF; the BEM feathering solve ignores it.
    """

    def __init__(
        self,
        I_theta_kgm2,
        damper_Nms_per_rad,
        flap,
        ac_offset_m=0.0,
        control_stiffness_Nm_per_rad=0.0,
    ):
        self.I_theta_kgm2 = I_theta_kgm2
        self.damper_Nms_per_rad = damper_Nms_per_rad
        self.flap = flap
        self.ac_offset_m = ac_offset_m
        self.control_stiffness_Nm_per_rad = control_stiffness_Nm_per_rad

    def _to_rust(self):
        return _RustServoFlapActuation(
            I_theta_kgm2=float(self.I_theta_kgm2),
            damper_Nms_per_rad=float(self.damper_Nms_per_rad),
            flap=self.flap._to_rust(),
            ac_offset_m=float(self.ac_offset_m),
            control_stiffness_Nm_per_rad=float(self.control_stiffness_Nm_per_rad),
        )


class FlapProperties:
    """Quasi-static blade flapping properties (hingeless / equivalent-hinge).

    Models the blade's out-of-plane flexibility as an equivalent spring-hinge.
    The blade absorbs most of the aerodynamic pitching/rolling moment; only the
    fraction determined by the flap frequency ratio reaches the hub (airframe).

    Fields:
      I_blade_flap_kgm2  -- blade flap inertia about (virtual) hinge [kg*m^2]
      omega_nr_rad_s     -- non-rotating flap natural frequency [rad/s]
                            (from root bending stiffness: K_beta = I_b * omega_NR^2)
                            Set to 0.0 for a freely-hinged blade (no spring).
    """

    def __init__(self, I_blade_flap_kgm2, omega_nr_rad_s=0.0):
        self.I_blade_flap_kgm2 = I_blade_flap_kgm2
        self.omega_nr_rad_s = omega_nr_rad_s

    def _to_rust(self):
        return _RustFlapProperties(
            I_blade_flap_kgm2=float(self.I_blade_flap_kgm2),
            omega_nr_rad_s=float(self.omega_nr_rad_s),
        )

    def hub_moment_factor(self, omega_rad_s):
        """Fraction of aero hub moment that passes to airframe at given omega."""
        return self._to_rust().hub_moment_factor(omega_rad_s)


class InertiaProperties:
    """Rotor inertia data (Python-only, not used in Rust math)."""

    def __init__(
        self,
        mass_kg=None,
        I_body_kgm2=None,
        I_spin_kgm2=None,
        blade_mass_kg=None,
        stationary_assembly_mass_kg=None,
        spinning_hub_shell_mass_kg=None,
        I_blade_flap_kgm2=None,
    ):
        self.mass_kg = mass_kg
        self.I_body_kgm2 = list(I_body_kgm2) if I_body_kgm2 is not None else []
        self.I_spin_kgm2 = I_spin_kgm2
        self.blade_mass_kg = blade_mass_kg
        self.stationary_assembly_mass_kg = stationary_assembly_mass_kg
        self.spinning_hub_shell_mass_kg = spinning_hub_shell_mass_kg
        self.I_blade_flap_kgm2 = I_blade_flap_kgm2


class LinearPolarParameters:
    """Linear polar parameters for the BEM solver."""

    def __init__(
        self,
        CL0,
        CL_alpha_per_rad,
        CD0,
        alpha_stall_deg,
        Re_design=None,
        name="",
        source="",
        polar_csv=None,
        CD_structural=0.0,
        Re_operating=None,
    ):
        self.Re_design = Re_design
        self.CL0 = CL0
        self.CL_alpha_per_rad = CL_alpha_per_rad
        self.CD0 = CD0
        self.alpha_stall_deg = alpha_stall_deg
        self.name = name
        self.source = source
        self.polar_csv = polar_csv
        self.CD_structural = CD_structural
        self.Re_operating = Re_operating


class ControlProperties:
    """Control / swashplate properties."""

    def __init__(
        self,
        swashplate_pitch_gain_rad,
        axle_attachment_length_m=None,
        K_cyc=None,
        swashplate_phase_deg=None,
        servo_slew_rate_deg_s=None,
        servo_travel_deg=None,
        kaman_flap=None,
    ):
        self.swashplate_pitch_gain_rad = swashplate_pitch_gain_rad
        self.axle_attachment_length_m = axle_attachment_length_m
        self.K_cyc = K_cyc
        self.swashplate_phase_deg = swashplate_phase_deg
        self.servo_slew_rate_deg_s = servo_slew_rate_deg_s
        self.servo_travel_deg = servo_travel_deg
        self.kaman_flap = kaman_flap


class AutorotationProperties:
    """Autorotation metadata (Python-only; I_ode_kgm2 used by mechanical utilities)."""

    def __init__(
        self,
        I_ode_kgm2=None,
        omega_min_rad_s=None,
        omega_eq_rad_s=None,
    ):
        self.I_ode_kgm2 = I_ode_kgm2
        self.omega_min_rad_s = omega_min_rad_s
        self.omega_eq_rad_s = omega_eq_rad_s


class RotorDefinition:
    """Rotor definition (all fields in Python).

    Pass to ``create_aero`` or the model constructors.  The Python model
    wrappers build the Rust representation on-demand via ``_to_rust_defn``.
    """

    def __init__(
        self,
        blade,
        airfoil,
        control=None,
        inertia=None,
        autorotation=None,
        servoflap=None,
        flap=None,
        name="",
        description="",
    ):
        auto = autorotation if autorotation is not None else AutorotationProperties()
        self.blade = blade
        self.airfoil = airfoil
        self.control = control
        self.inertia = inertia if inertia is not None else InertiaProperties()
        self.autorotation = auto
        self.servoflap = servoflap
        self.flap = flap
        self.name = name
        self.description = description

    # Convenience geometry properties that delegate to blade.
    @property
    def span_m(self):
        return self.blade.span_m

    @property
    def r_cp_m(self):
        return self.blade.r_cp_m

    @property
    def disk_area_m2(self):
        return self.blade.disk_area_m2

    @property
    def solidity(self):
        return self.blade.solidity


# ---------------------------------------------------------------------------
# ValidationIssue + validate() methods (legacy Python-side helper).
# ---------------------------------------------------------------------------

@dataclass
class ValidationIssue:
    level: str
    field: str
    message: str

    def __str__(self) -> str:
        return f"[{self.level}] {self.field}: {self.message}"


def _validate_blade(self) -> List[ValidationIssue]:
    issues: List[ValidationIssue] = []
    if self.radius_m <= 0:
        issues.append(ValidationIssue("ERROR", "blade.radius_m", "must be positive"))
    if self.chord_m <= 0:
        issues.append(ValidationIssue("ERROR", "blade.chord_m", "must be positive"))
    if self.span_m <= 0:
        issues.append(ValidationIssue("ERROR", "blade.span_m",
                                       "root_cutout_m must be less than radius_m"))
    if self.n_blades < 1:
        issues.append(ValidationIssue("ERROR", "blade.n_blades", "must be >= 1"))
    return issues


def _validate_airfoil(self) -> List[ValidationIssue]:
    issues: List[ValidationIssue] = []
    if self.CD0 < 0:
        issues.append(ValidationIssue("ERROR", "airfoil.CD0", "must be >= 0"))
    return issues


def _validate_definition(self) -> List[ValidationIssue]:
    return self.blade.validate() + self.airfoil.validate()


# Attach validate() to BladeGeometry (Rust pyclass) and Python classes.
BladeGeometry.validate = _validate_blade
LinearPolarParameters.validate = _validate_airfoil
RotorDefinition.validate = _validate_definition


def _to_rust_control(control):
    """Build a lean ``_dynbem.ControlProperties`` from a Python wrapper.

    Returns ``None`` for ``None`` input; passes through a Rust object unchanged.
    """
    if control is None or isinstance(control, _RustControlProperties):
        return control
    return _RustControlProperties(
        swashplate_pitch_gain_rad=control.swashplate_pitch_gain_rad,
        swashplate_phase_deg=control.swashplate_phase_deg,
    )


def _to_rust_defn(defn):
    if isinstance(defn, _RustRotorDefinition):
        return defn
    airfoil = defn.airfoil
    control = defn.control
    servoflap = getattr(defn, "servoflap", None)
    flap = getattr(defn, "flap", None)
    return _RustRotorDefinition(
        blade=defn.blade,
        airfoil=_RustLinearPolarParameters(
            CL0=airfoil.CL0,
            CL_alpha_per_rad=airfoil.CL_alpha_per_rad,
            CD0=airfoil.CD0,
            alpha_stall_deg=airfoil.alpha_stall_deg,
        ),
        control=_RustControlProperties(
            swashplate_pitch_gain_rad=control.swashplate_pitch_gain_rad,
            swashplate_phase_deg=control.swashplate_phase_deg,
        ) if control is not None else None,
        servoflap=servoflap._to_rust() if servoflap is not None else None,
        name=defn.name,
        description=defn.description,
        flap=flap._to_rust() if flap is not None else None,
    )


def load(path: str) -> RotorDefinition:
    """Load a RotorDefinition from a YAML file.

    ``airfoil.polar_csv``, if relative, is resolved against the YAML
    file's parent directory.
    """
    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(file_path)
    text = file_path.read_text(encoding="utf-8")
    return loads(text, base_dir=str(file_path.parent))


# ---------------------------------------------------------------------------
# YAML helpers
# ---------------------------------------------------------------------------

def _maybe_sum(value: Any) -> Optional[float]:
    """Accept a scalar or a dict whose values are summed (legacy schema)."""
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, dict):
        return float(sum(value.values()))
    raise ValueError(f"expected number or dict, got {type(value).__name__}")


def _resolve_polar_csv(value: Optional[str], base_dir: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    p = Path(value)
    if p.is_absolute() or base_dir is None:
        return value
    return str(Path(base_dir) / p)


def _build_kaman_flap(k: Dict[str, Any]) -> KamanFlap:
    return KamanFlap(
        chord_fraction=k.get("chord_fraction"),
        span_start_m=k.get("span_start_m"),
        span_end_m=k.get("span_end_m"),
        tau=k.get("tau"),
        CM_gamma_per_rad=k.get("CM_gamma_per_rad"),
        swashplate_load_fraction=k.get("swashplate_load_fraction"),
        notes=k.get("notes", "") or "",
    )


def _build_from_dict(doc: Dict[str, Any], base_dir: Optional[str]) -> RotorDefinition:
    r = doc.get("rotor") or {}
    a = doc.get("airfoil") or {}
    i = doc.get("inertia") or {}
    c_raw = doc.get("control")
    ar = doc.get("autorotation") or {}

    def _req(mapping: Dict[str, Any], key: str, section: str) -> Any:
        v = mapping.get(key)
        if v is None:
            raise ValueError(f"missing required field: {section}.{key}")
        return v

    n_blades = int(_req(r, "n_blades", "rotor"))

    blade = BladeGeometry(
        n_blades=n_blades,
        radius_m=float(_req(r, "radius_m", "rotor")),
        root_cutout_m=float(_req(r, "root_cutout_m", "rotor")),
        chord_m=float(_req(r, "chord_m", "rotor")),
        twist_deg=float(r.get("twist_deg") or 0.0),
        n_elements=int(r.get("n_elements") or 10),
        r_stations_m=list(r.get("r_stations_m") or []),
        chord_stations_m=list(r.get("chord_stations_m") or []),
        twist_stations_deg=list(r.get("twist_stations_deg") or []),
        tip_loss=bool(r.get("tip_loss", True)),
    )

    airfoil = LinearPolarParameters(
        Re_design=int(_req(a, "Re_design", "airfoil")),
        CL0=float(_req(a, "CL0", "airfoil")),
        CL_alpha_per_rad=float(_req(a, "CL_alpha_per_rad", "airfoil")),
        CD0=float(_req(a, "CD0", "airfoil")),
        alpha_stall_deg=float(_req(a, "alpha_stall_deg", "airfoil")),
        name=a.get("designation") or "",
        source=a.get("source") or "",
        polar_csv=_resolve_polar_csv(a.get("polar_csv"), base_dir),
        CD_structural=float(a.get("CD_structural") or 0.0),
        Re_operating=int(a["Re_operating"]) if a.get("Re_operating") is not None else None,
    )

    # Inertia: resolve mass from components if explicit mass_kg is absent.
    explicit_mass = _maybe_sum(i.get("mass_kg"))
    if explicit_mass is None:
        blade_mass = _maybe_sum(i.get("blade_mass_kg"))
        stat_mass = _maybe_sum(i.get("stationary_assembly_mass_kg"))
        shell_mass = _maybe_sum(i.get("spinning_hub_shell_mass_kg"))
        if blade_mass is not None and stat_mass is not None and shell_mass is not None:
            explicit_mass = blade_mass * n_blades + stat_mass + shell_mass

    inertia = InertiaProperties(
        mass_kg=explicit_mass,
        I_body_kgm2=list(i.get("I_body_kgm2") or []),
        I_spin_kgm2=_maybe_sum(i.get("I_spin_kgm2")),
        blade_mass_kg=_maybe_sum(i.get("blade_mass_kg")),
        stationary_assembly_mass_kg=_maybe_sum(i.get("stationary_assembly_mass_kg")),
        spinning_hub_shell_mass_kg=_maybe_sum(i.get("spinning_hub_shell_mass_kg")),
        I_blade_flap_kgm2=_maybe_sum(i.get("I_blade_flap_kgm2")),
    )

    control = None
    if c_raw is not None:
        kf_raw = c_raw.get("kaman_flap")
        control = ControlProperties(
            swashplate_pitch_gain_rad=float(_req(c_raw, "swashplate_pitch_gain_rad", "control")),
            axle_attachment_length_m=c_raw.get("axle_attachment_length_m"),
            K_cyc=c_raw.get("K_cyc"),
            swashplate_phase_deg=c_raw.get("swashplate_phase_deg"),
            servo_slew_rate_deg_s=c_raw.get("servo_slew_rate_deg_s"),
            servo_travel_deg=c_raw.get("servo_travel_deg"),
            kaman_flap=_build_kaman_flap(kf_raw) if kf_raw else None,
        )

    autorotation = AutorotationProperties(
        I_ode_kgm2=ar.get("I_ode_kgm2"),
        omega_min_rad_s=ar.get("omega_min_rad_s"),
        omega_eq_rad_s=ar.get("omega_eq_rad_s"),
    )

    # Pitch actuation: DirectMechanical (default) or ServoFlap.
    # YAML mirrors the Rust `PitchActuation` enum:
    #   pitch_actuation:
    #     servoflap:
    #       I_theta_kgm2: ...
    #       damper_Nms_per_rad: ...
    #       ac_offset_m: ...
    #       flap: { C_M_delta_per_rad: ..., r_inner_m: ..., r_outer_m: ... }
    # Absent `pitch_actuation` (or absent `servoflap`) => direct mechanical.
    servoflap = None
    pa_raw = doc.get("pitch_actuation")
    if pa_raw is not None:
        sf_raw = pa_raw.get("servoflap")
        if sf_raw is not None:
            flap_raw = sf_raw.get("flap")
            if flap_raw is None:
                raise ValueError(
                    "missing required field: pitch_actuation.servoflap.flap"
                )
            flap = ServoFlapGeometry(
                C_M_delta_per_rad=float(_req(flap_raw, "C_M_delta_per_rad",
                                             "pitch_actuation.servoflap.flap")),
                r_inner_m=float(_req(flap_raw, "r_inner_m",
                                     "pitch_actuation.servoflap.flap")),
                r_outer_m=float(_req(flap_raw, "r_outer_m",
                                     "pitch_actuation.servoflap.flap")),
            )
            servoflap = ServoFlapActuation(
                I_theta_kgm2=float(_req(sf_raw, "I_theta_kgm2",
                                        "pitch_actuation.servoflap")),
                damper_Nms_per_rad=float(sf_raw.get("damper_Nms_per_rad") or 0.0),
                flap=flap,
                ac_offset_m=float(sf_raw.get("ac_offset_m") or 0.0),
                control_stiffness_Nm_per_rad=float(
                    sf_raw.get("control_stiffness_Nm_per_rad") or 0.0
                ),
            )

    # Blade flapping properties (quasi-static hub moment reduction).
    # YAML schema:
    #   flap:
    #     I_blade_flap_kgm2: ...
    #     omega_nr_rad_s: ...        (optional, default 0.0 = pure hinge)
    flap_props = None
    flap_raw = doc.get("flap")
    if flap_raw is not None:
        flap_props = FlapProperties(
            I_blade_flap_kgm2=float(_req(flap_raw, "I_blade_flap_kgm2", "flap")),
            omega_nr_rad_s=float(flap_raw.get("omega_nr_rad_s") or 0.0),
        )

    return RotorDefinition(
        blade=blade,
        airfoil=airfoil,
        control=control,
        inertia=inertia,
        autorotation=autorotation,
        servoflap=servoflap,
        flap=flap_props,
        name=doc.get("name") or "",
        description=doc.get("description") or "",
    )


def loads(text: str, base_dir: Optional[str] = None) -> RotorDefinition:
    """Parse a RotorDefinition from a YAML string.

    ``base_dir``, if provided, is the directory used to resolve a
    relative ``airfoil.polar_csv`` path; otherwise the value is stored
    verbatim.
    """
    doc = yaml.safe_load(text)
    return _build_from_dict(doc, base_dir)


def default() -> RotorDefinition:
    raise NotImplementedError("No project default rotor is defined.")
