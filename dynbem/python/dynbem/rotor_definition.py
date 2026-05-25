"""dynbem.rotor_definition submodule (compat shim).

Re-exports the rotor definition pyclasses from the Rust extension.
YAML loading is handled in pure Python here; the Rust core holds only
the struct definitions. The ``ValidationIssue`` dataclass + ``.validate()``
shims remain for backwards compatibility with the legacy pure-Python dynbem.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

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
    "loads",
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
    )

    airfoil = AirfoilProperties(
        Re_design=int(_req(a, "Re_design", "airfoil")),
        CL0=float(_req(a, "CL0", "airfoil")),
        CL_alpha_per_rad=float(_req(a, "CL_alpha_per_rad", "airfoil")),
        CD0=float(_req(a, "CD0", "airfoil")),
        alpha_stall_deg=float(_req(a, "alpha_stall_deg", "airfoil")),
        tip_loss=bool(a.get("tip_loss", True)),
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
            swashplate_pitch_gain_rad=float(
                _req(c_raw, "swashplate_pitch_gain_rad", "control")
            ),
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

    return RotorDefinition(
        blade=blade,
        airfoil=airfoil,
        control=control,
        inertia=inertia,
        autorotation=autorotation,
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
