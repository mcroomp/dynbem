"""dynbem.rotor_definition submodule (compat shim).

Re-exports the rotor definition pyclasses from the Rust extension and
ports the YAML loader + ValidationIssue layer from the legacy
pure-Python implementation. The loader and validation stay in Python
because they depend on PyYAML / dataclasses and don't belong inside a
Rust crate.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from . import (
    AirfoilProperties,
    AutorotationProperties,
    BladeGeometry,
    ControlProperties,
    InertiaProperties,
    KamanFlap,
    RotorDefinition,
)

__all__ = [
    "AirfoilProperties",
    "AutorotationProperties",
    "BladeGeometry",
    "ControlProperties",
    "InertiaProperties",
    "KamanFlap",
    "RotorDefinition",
    "ValidationIssue",
    "load",
    "default",
]


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


# Attach the validate() methods to the Rust pyclasses. Pyo3 pyclass types
# allow class-level attribute assignment by default; if a future pyo3
# version sets the type as frozen, switch these to free functions and
# update the tests to call validate_blade(b) / validate_airfoil(a).
BladeGeometry.validate = _validate_blade
AirfoilProperties.validate = _validate_airfoil
RotorDefinition.validate = _validate_definition


def _require_yaml():
    try:
        import yaml  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "Loading rotor YAML requires PyYAML. Install with: "
            "pip install PyYAML"
        ) from exc


def _maybe_float(value) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, dict):
        return float(sum(float(v) for v in value.values()))
    return float(value)


def _maybe_int(value) -> Optional[int]:
    return None if value is None else int(value)


def _resolve_polar_csv(value, yaml_path: Path) -> Optional[str]:
    if value is None:
        return None
    p = Path(str(value))
    if not p.is_absolute():
        p = yaml_path.parent / p
    return str(p)


def _resolve_mass_kg(inertia_section: dict, n_blades: int) -> Optional[float]:
    explicit = inertia_section.get("mass_kg")
    if explicit is not None:
        return float(explicit)
    blade = _maybe_float(inertia_section.get("blade_mass_kg"))
    stat = _maybe_float(inertia_section.get("stationary_assembly_mass_kg"))
    shell = _maybe_float(inertia_section.get("spinning_hub_shell_mass_kg"))
    if blade is not None and stat is not None and shell is not None:
        return blade * n_blades + stat + shell
    return None


def load(path: str) -> RotorDefinition:
    """Load a RotorDefinition from a YAML file.

    The YAML schema matches the legacy dynbem layout: top-level sections
    `rotor` (blade geometry), `airfoil`, `inertia`, `control`,
    `autorotation`, plus optional `name` / `description`.
    """
    _require_yaml()
    import yaml

    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(file_path)
    with file_path.open(encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}

    rotor = data.get("rotor", {})
    airfoil = data.get("airfoil", {})
    inertia = data.get("inertia", {})
    control = data.get("control", {})
    autorotation = data.get("autorotation", {})

    return RotorDefinition(
        name=str(data.get("name", "")),
        description=str(data.get("description", "")),
        blade=BladeGeometry(
            n_blades=int(rotor["n_blades"]),
            radius_m=float(rotor["radius_m"]),
            root_cutout_m=float(rotor["root_cutout_m"]),
            chord_m=float(rotor["chord_m"]),
            twist_deg=float(rotor.get("twist_deg", 0.0)),
            n_elements=int(rotor.get("n_elements", 10)),
            r_stations_m=[float(v) for v in rotor.get("r_stations_m", [])],
            chord_stations_m=[float(v) for v in rotor.get("chord_stations_m", [])],
            twist_stations_deg=[float(v) for v in rotor.get("twist_stations_deg", [])],
        ),
        airfoil=AirfoilProperties(
            Re_design=int(airfoil["Re_design"]),
            CL0=float(airfoil["CL0"]),
            CL_alpha_per_rad=float(airfoil["CL_alpha_per_rad"]),
            CD0=float(airfoil["CD0"]),
            alpha_stall_deg=float(airfoil["alpha_stall_deg"]),
            tip_loss=bool(airfoil.get("tip_loss", True)),
            name=str(airfoil.get("designation", "")),
            source=str(airfoil.get("source", "")),
            polar_csv=_resolve_polar_csv(airfoil.get("polar_csv"), file_path),
            CD_structural=float(airfoil.get("CD_structural", 0.0)),
            Re_operating=_maybe_int(airfoil.get("Re_operating")),
        ),
        control=ControlProperties(
            swashplate_pitch_gain_rad=float(control["swashplate_pitch_gain_rad"]),
            axle_attachment_length_m=_maybe_float(control.get("axle_attachment_length_m")),
            K_cyc=_maybe_float(control.get("K_cyc")),
            swashplate_phase_deg=_maybe_float(control.get("swashplate_phase_deg")),
            servo_slew_rate_deg_s=_maybe_float(control.get("servo_slew_rate_deg_s")),
            servo_travel_deg=_maybe_float(control.get("servo_travel_deg")),
            kaman_flap=KamanFlap(**control.get("kaman_flap", {})),
        ) if control else None,
        inertia=InertiaProperties(
            mass_kg=_resolve_mass_kg(inertia, n_blades=int(rotor["n_blades"])),
            I_body_kgm2=[float(v) for v in inertia.get("I_body_kgm2", [])],
            I_spin_kgm2=_maybe_float(inertia.get("I_spin_kgm2")),
            blade_mass_kg=_maybe_float(inertia.get("blade_mass_kg")),
            stationary_assembly_mass_kg=_maybe_float(inertia.get("stationary_assembly_mass_kg")),
            spinning_hub_shell_mass_kg=_maybe_float(inertia.get("spinning_hub_shell_mass_kg")),
            I_blade_flap_kgm2=_maybe_float(inertia.get("I_blade_flap_kgm2")),
        ),
        autorotation=AutorotationProperties(
            I_ode_kgm2=_maybe_float(autorotation.get("I_ode_kgm2")),
            omega_min_rad_s=_maybe_float(autorotation.get("omega_min_rad_s")),
            omega_eq_rad_s=_maybe_float(autorotation.get("omega_eq_rad_s")),
        ),
    )


def default() -> RotorDefinition:
    raise NotImplementedError("No project default rotor is defined.")
