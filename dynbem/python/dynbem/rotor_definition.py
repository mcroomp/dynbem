"""dynbem.rotor_definition submodule.

Defines pure-Python classes for all rotor definition types, including
metadata fields not needed by the Rust math core (Re_design, polar_csv,
inertia, etc.). Each Python class internally holds a ``_rust``
attribute -- a lean _dynbem Rust wrapper populated with only the math fields
-- which is what model constructors receive.

YAML loading is handled here in pure Python via yaml.safe_load.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

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
        cls: type[BladeGeometry],
        n_blades: int,
        radius_m: float,
        root_cutout_m: float,
        chord_m: float,
        twist_deg: float = 0.0,
        n_elements: int = 20,
        r_stations_m: list[float] | None = None,
        chord_stations_m: list[float] | None = None,
        twist_stations_deg: list[float] | None = None,
        tip_loss: bool = True,
    ) -> BladeGeometry:
        return cast(
            BladeGeometry,
            _RustBladeGeometry.__new__(
                cls,
                n_blades,
                radius_m,
                root_cutout_m,
                chord_m,
                twist_deg,
                n_elements,
                r_stations_m if r_stations_m is not None else [],
                chord_stations_m if chord_stations_m is not None else [],
                twist_stations_deg if twist_stations_deg is not None else [],
                tip_loss,
            ),
        )

    def __init__(
        self,
        n_blades: int,
        radius_m: float,
        root_cutout_m: float,
        chord_m: float,
        twist_deg: float = 0.0,
        n_elements: int = 20,
        r_stations_m: list[float] | None = None,
        chord_stations_m: list[float] | None = None,
        twist_stations_deg: list[float] | None = None,
        tip_loss: bool = True,
    ) -> None:
        pass

    def validate(self) -> list[ValidationIssue]:
        return _validate_blade(self)

__all__ = [
    "LinearPolarParameters",
    "AutorotationProperties",
    "BladeGeometry",
    "ControlProperties",
    "FlapProperties",
    "ServoFlapActuation",
    "ServoFlapGeometry",
    "InertiaProperties",
    "RotorDefinition",
    "ValidationIssue",
    "load",
    "loads",
    "default",
]


class ServoFlapGeometry:
    """Servo-flap (trailing-edge elevon) geometry.

    In servo-flap mode the swashplate collective AND cyclic drive flap
    deflection delta_f(psi). The flap exerts an aerodynamic pitching moment
    about the feathering axis. Used inside ServoFlapActuation.flap.
    """

    def __init__(
        self,
        C_M_delta_per_rad: float,
        r_inner_m: float,
        r_outer_m: float,
    ) -> None:
        self.C_M_delta_per_rad = C_M_delta_per_rad
        self.r_inner_m = r_inner_m
        self.r_outer_m = r_outer_m

class ServoFlapActuation:
    """Servo-flap pitch actuation: feathering + damper architecture.

    This models the feathering + damper servo-flap type: the blade rides a pitch
    bearing and feathers freely as a rigid body, restrained by the bearing damper
    and the aerodynamic spring from any AC offset. A trailing-edge servo-flap
    exerts a pitching moment about the feathering axis. Both collective and cyclic
    are interpreted as flap commands and the resulting delta_theta(psi) REPLACES
    the direct swashplate pitch path in the psi-loop.

    (The alternative torsional-twist architecture -- a rigid-pitch blade
    elastically twisted against its spar stiffness -- is not modelled yet.)

    All fields are measurable physical constants:
      I_theta_kgm2        -- pitch moment of inertia about feathering axis [kg*m^2]
                            (bifilar / swing test)
      damper_Nms_per_rad  -- bearing rotary damper coefficient [N*m*s/rad]
                            (bearing drag test). The only dissipation; in the
                            damper-dominated limit it sets the ~90 deg cyclic lag.
      flap                -- ServoFlapGeometry (required)
      ac_offset_m         -- distance from feathering axis to AC (~25% chord) [m],
                            positive = AC aft of axis (stable aero spring), 0 =
                            axis at AC (Kaman ideal, no spring). Ruler on the
                            section, or a calibrated-moment RPM sweep.
      blade_Cm_AC         -- blade section zero-lift pitching moment coefficient [-]
                            about the AC (negative for cambered airfoils). Sets the
                            physical DC trim. Measurable: run flap-undeflected and
                            read the settled pitch. 0.0 = symmetric or unknown.

    When the blade pitches up it makes more lift; because that lift acts at the
    AC a distance ac_offset_m aft of the pitch axis, it creates a nose-down
    restoring torque proportional to the pitch change (an aerodynamic spring),
    while blade_Cm_AC adds the airfoil's constant built-in twisting moment that
    offsets where this spring settles.
    """

    def __init__(
        self,
        I_theta_kgm2: float,
        damper_Nms_per_rad: float,
        flap: ServoFlapGeometry,
        ac_offset_m: float = 0.0,
        blade_Cm_AC: float = 0.0,
    ) -> None:
        self.I_theta_kgm2 = I_theta_kgm2
        self.damper_Nms_per_rad = damper_Nms_per_rad
        self.flap = flap
        self.ac_offset_m = ac_offset_m
        self.blade_Cm_AC = blade_Cm_AC

    def _to_rust(self) -> _RustServoFlapActuation:
        return _RustServoFlapActuation(
            I_theta_kgm2=float(self.I_theta_kgm2),
            damper_Nms_per_rad=float(self.damper_Nms_per_rad),
            flap=_RustServoFlapGeometry(
                C_M_delta_per_rad=self.flap.C_M_delta_per_rad,
                r_inner_m=self.flap.r_inner_m,
                r_outer_m=self.flap.r_outer_m,
            ),
            ac_offset_m=float(self.ac_offset_m),
            blade_Cm_AC=float(self.blade_Cm_AC),
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

    def __init__(
        self, I_blade_flap_kgm2: float, omega_nr_rad_s: float = 0.0
    ) -> None:
        self.I_blade_flap_kgm2 = I_blade_flap_kgm2
        self.omega_nr_rad_s = omega_nr_rad_s

class InertiaProperties:
    """Rotor inertia data (Python-only, not used in Rust math)."""

    def __init__(
        self,
        mass_kg: float | None = None,
        I_body_kgm2: list[float] | None = None,
        I_spin_kgm2: float | None = None,
        blade_mass_kg: float | None = None,
        stationary_assembly_mass_kg: float | None = None,
        spinning_hub_shell_mass_kg: float | None = None,
        I_blade_flap_kgm2: float | None = None,
    ) -> None:
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
        CL0: float,
        CL_alpha_per_rad: float,
        CD0: float,
        alpha_stall_deg: float,
        Re_design: int | None = None,
        name: str = "",
        source: str = "",
        polar_csv: str | None = None,
        CD_structural: float = 0.0,
        Re_operating: int | None = None,
    ) -> None:
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

    def validate(self) -> list[ValidationIssue]:
        return _validate_airfoil(self)


class ControlProperties:
    """Control / swashplate properties."""

    def __init__(
        self,
        swashplate_pitch_gain_rad: float,
        axle_attachment_length_m: float | None = None,
        K_cyc: float | None = None,
        swashplate_phase_deg: float | None = None,
        servo_slew_rate_deg_s: float | None = None,
        servo_travel_deg: float | None = None,
    ) -> None:
        self.swashplate_pitch_gain_rad = swashplate_pitch_gain_rad
        self.axle_attachment_length_m = axle_attachment_length_m
        self.K_cyc = K_cyc
        self.swashplate_phase_deg = swashplate_phase_deg
        self.servo_slew_rate_deg_s = servo_slew_rate_deg_s
        self.servo_travel_deg = servo_travel_deg


class AutorotationProperties:
    """Autorotation metadata (Python-only; I_ode_kgm2 used by mechanical utilities)."""

    def __init__(
        self,
        I_ode_kgm2: float | None = None,
        omega_min_rad_s: float | None = None,
        omega_eq_rad_s: float | None = None,
    ) -> None:
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
        blade: BladeGeometry,
        airfoil: LinearPolarParameters,
        control: ControlProperties | None = None,
        inertia: InertiaProperties | None = None,
        autorotation: AutorotationProperties | None = None,
        servoflap: ServoFlapActuation | None = None,
        flap: FlapProperties | None = None,
        name: str = "",
        description: str = "",
    ) -> None:
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
    def span_m(self) -> float:
        return self.blade.span_m

    @property
    def r_cp_m(self) -> float:
        return self.blade.r_cp_m

    @property
    def disk_area_m2(self) -> float:
        return self.blade.disk_area_m2

    @property
    def solidity(self) -> float:
        return self.blade.solidity

    def validate(self) -> list[ValidationIssue]:
        return _validate_definition(self)


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


def _validate_blade(self: BladeGeometry) -> list[ValidationIssue]:
    issues: list[ValidationIssue] = []
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


def _validate_airfoil(self: LinearPolarParameters) -> list[ValidationIssue]:
    issues: list[ValidationIssue] = []
    if self.CD0 < 0:
        issues.append(ValidationIssue("ERROR", "airfoil.CD0", "must be >= 0"))
    return issues


def _validate_definition(self: RotorDefinition) -> list[ValidationIssue]:
    return self.blade.validate() + self.airfoil.validate()


def _to_rust_control(  # pyright: ignore[reportUnusedFunction]
    control: ControlProperties | _RustControlProperties | None,
) -> _RustControlProperties | None:
    """Build a lean ``_dynbem.ControlProperties`` from a Python wrapper.

    Returns ``None`` for ``None`` input; passes through a Rust object unchanged.
    """
    if control is None or isinstance(control, _RustControlProperties):
        return control
    return _RustControlProperties(
        swashplate_pitch_gain_rad=control.swashplate_pitch_gain_rad,
        swashplate_phase_deg=control.swashplate_phase_deg,
    )


def _to_rust_defn(  # pyright: ignore[reportUnusedFunction]
    defn: RotorDefinition | _RustRotorDefinition,
) -> _RustRotorDefinition:
    if isinstance(defn, _RustRotorDefinition):
        return defn
    airfoil = defn.airfoil
    control = defn.control
    servoflap = defn.servoflap
    flap = defn.flap
    servoflap_rust = (
        servoflap._to_rust()  # pyright: ignore[reportPrivateUsage]
        if servoflap is not None
        else None
    )
    flap_rust = (
        _RustFlapProperties(
            I_blade_flap_kgm2=float(flap.I_blade_flap_kgm2),
            omega_nr_rad_s=float(flap.omega_nr_rad_s),
        )
        if flap is not None
        else None
    )
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
        servoflap=servoflap_rust,
        name=defn.name,
        description=defn.description,
        flap=flap_rust,
    )


def load(path: str | Path) -> RotorDefinition:
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

def _maybe_sum(value: Any) -> float | None:
    """Accept a scalar or a dict whose values are summed (legacy schema)."""
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, dict):
        return float(sum(cast(dict[str, Any], value).values()))
    raise ValueError(f"expected number or dict, got {type(value).__name__}")


def _resolve_polar_csv(value: str | None, base_dir: str | None) -> str | None:
    if value is None:
        return None
    p = Path(value)
    if p.is_absolute() or base_dir is None:
        return value
    return str(Path(base_dir) / p)


def _mapping(value: Any, section: str) -> dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise TypeError(f"expected a mapping for {section}")
    if not all(isinstance(key, str) for key in cast(dict[object, Any], value)):
        raise ValueError(f"expected string keys in {section}")
    return cast(dict[str, Any], value)


def _build_from_dict(doc: dict[str, Any], base_dir: str | None) -> RotorDefinition:
    r = _mapping(doc.get("rotor"), "rotor")
    a = _mapping(doc.get("airfoil"), "airfoil")
    i = _mapping(doc.get("inertia"), "inertia")
    control_value = doc.get("control")
    c_raw = _mapping(control_value, "control") if control_value is not None else None
    ar = _mapping(doc.get("autorotation"), "autorotation")

    def _req(mapping: dict[str, Any], key: str, section: str) -> Any:
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
        if "kaman_flap" in c_raw:
            raise ValueError(
                "control.kaman_flap is no longer supported; configure active "
                "servo-flap dynamics with pitch_actuation.servoflap"
            )
        control = ControlProperties(
            swashplate_pitch_gain_rad=float(_req(c_raw, "swashplate_pitch_gain_rad", "control")),
            axle_attachment_length_m=c_raw.get("axle_attachment_length_m"),
            K_cyc=c_raw.get("K_cyc"),
            swashplate_phase_deg=c_raw.get("swashplate_phase_deg"),
            servo_slew_rate_deg_s=c_raw.get("servo_slew_rate_deg_s"),
            servo_travel_deg=c_raw.get("servo_travel_deg"),
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
    pa_value = doc.get("pitch_actuation")
    if pa_value is not None:
        pa_raw = _mapping(pa_value, "pitch_actuation")
        sf_value = pa_raw.get("servoflap")
        sf_raw = _mapping(sf_value, "pitch_actuation.servoflap") if sf_value is not None else None
        if sf_raw is not None:
            flap_value = sf_raw.get("flap")
            flap_raw = _mapping(flap_value, "pitch_actuation.servoflap.flap") if flap_value is not None else None
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
                blade_Cm_AC=float(sf_raw.get("blade_Cm_AC") or 0.0),
            )

    # Blade flapping properties (quasi-static hub moment reduction).
    # YAML schema:
    #   flap:
    #     I_blade_flap_kgm2: ...
    #     omega_nr_rad_s: ...        (optional, default 0.0 = pure hinge)
    flap_props = None
    flap_value = doc.get("flap")
    if flap_value is not None:
        flap_raw = _mapping(flap_value, "flap")
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


def loads(text: str, base_dir: str | None = None) -> RotorDefinition:
    """Parse a RotorDefinition from a YAML string.

    ``base_dir``, if provided, is the directory used to resolve a
    relative ``airfoil.polar_csv`` path; otherwise the value is stored
    verbatim.
    """
    doc = yaml.safe_load(text)
    if not isinstance(doc, dict):
        raise ValueError("rotor definition must be a YAML mapping")
    if not all(isinstance(key, str) for key in cast(dict[object, Any], doc)):
        raise ValueError("rotor definition keys must be strings")
    return _build_from_dict(cast(dict[str, Any], doc), base_dir)


def default() -> RotorDefinition:
    raise NotImplementedError("No project default rotor is defined.")
