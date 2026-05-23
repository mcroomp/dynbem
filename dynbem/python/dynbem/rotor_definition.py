"""dynbem.rotor_definition submodule (compat shim).

Re-exports the rotor definition pyclasses from the Rust extension and
delegates YAML loading to ``dynbem_rs::rotor_yaml`` via the
``_dynbem.load_rotor_yaml`` pyfunction, so pure-Rust callers and
Python callers share a single parser + schema. The Python-side
``ValidationIssue`` dataclass + ``.validate()`` shims remain here for
backwards compatibility with the legacy pure-Python dynbem.
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
from ._dynbem import load_rotor_yaml as _load_rotor_yaml
from ._dynbem import loads_rotor_yaml as _loads_rotor_yaml

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

    Delegates to ``dynbem_rs::rotor_yaml::from_yaml_file`` so Rust and
    Python share a single parser. ``airfoil.polar_csv``, if relative,
    is resolved against the YAML file's parent directory.

    The YAML schema matches the legacy dynbem layout: top-level
    sections ``rotor`` (blade geometry), ``airfoil``, ``inertia``,
    ``control``, ``autorotation``, plus optional ``name`` / ``description``.
    """
    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(file_path)
    return _load_rotor_yaml(str(file_path))


def loads(text: str, base_dir: Optional[str] = None) -> RotorDefinition:
    """Parse a RotorDefinition from a YAML string.

    ``base_dir``, if provided, is the directory used to resolve a
    relative ``airfoil.polar_csv`` path; otherwise the value is stored
    verbatim. Mirrors ``dynbem_rs::rotor_yaml::from_yaml_str``.
    """
    return _loads_rotor_yaml(text, base_dir)


def default() -> RotorDefinition:
    raise NotImplementedError("No project default rotor is defined.")
