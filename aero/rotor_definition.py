"""Rotor definition interfaces and generic YAML loading utilities."""

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None


@dataclass
class ValidationIssue:
    level: str
    field: str
    message: str

    def __str__(self) -> str:
        return f"[{self.level}] {self.field}: {self.message}"


@dataclass(frozen=True)
class KamanFlap:
    enabled: bool = False
    chord_fraction: Optional[float] = None
    span_start_m: Optional[float] = None
    span_end_m: Optional[float] = None
    tau: Optional[float] = None
    CM_gamma_per_rad: Optional[float] = None
    swashplate_load_fraction: Optional[float] = None
    notes: str = ""


@dataclass(frozen=True)
class BladeGeometry:
    n_blades: int
    radius_m: float
    root_cutout_m: float
    chord_m: float
    twist_deg: float = 0.0
    taper_ratio: float = 1.0
    n_elements: int = 10

    @property
    def span_m(self) -> float:
        return self.radius_m - self.root_cutout_m

    @property
    def r_cp_m(self) -> float:
        return self.root_cutout_m + (2.0 / 3.0) * self.span_m

    @property
    def S_w_m2(self) -> float:
        return self.n_blades * self.chord_m * self.span_m

    @property
    def disk_area_m2(self) -> float:
        return math.pi * (self.radius_m**2 - self.root_cutout_m**2)

    @property
    def aspect_ratio(self) -> float:
        return self.span_m / self.chord_m

    @property
    def solidity(self) -> float:
        return self.n_blades * self.chord_m / (math.pi * self.radius_m)

    def validate(self) -> list[ValidationIssue]:
        issues: list[ValidationIssue] = []
        if self.radius_m <= 0:
            issues.append(ValidationIssue("ERROR", "blade.radius_m", "must be positive"))
        if self.chord_m <= 0:
            issues.append(ValidationIssue("ERROR", "blade.chord_m", "must be positive"))
        if self.span_m <= 0:
            issues.append(ValidationIssue("ERROR", "blade.span_m", "root_cutout_m must be less than radius_m"))
        if self.n_blades < 1:
            issues.append(ValidationIssue("ERROR", "blade.n_blades", "must be >= 1"))
        return issues


@dataclass(frozen=True)
class AirfoilProperties:
    Re_design: int
    CL0: float
    CL_alpha_per_rad: float
    CD0: float
    alpha_stall_deg: float
    tip_loss: bool = True
    name: str = ""
    source: str = ""
    polar_csv: Optional[str] = None
    CD_structural: float = 0.0
    Re_operating: Optional[int] = None

    def validate(self) -> list[ValidationIssue]:
        issues: list[ValidationIssue] = []
        if self.CD0 < 0:
            issues.append(ValidationIssue("ERROR", "airfoil.CD0", "must be >= 0"))
        return issues


@dataclass(frozen=True)
class InertiaProperties:
    mass_kg: Optional[float] = None
    I_body_kgm2: tuple[float, ...] = field(default_factory=tuple)
    I_spin_kgm2: Optional[float] = None
    blade_mass_kg: Optional[float] = None
    stationary_assembly_mass_kg: Optional[float] = None
    spinning_hub_shell_mass_kg: Optional[float] = None
    I_blade_flap_kgm2: Optional[float] = None


@dataclass(frozen=True)
class ControlProperties:
    swashplate_pitch_gain_rad: float
    axle_attachment_length_m: Optional[float] = None
    K_cyc: Optional[float] = None
    swashplate_phase_deg: Optional[float] = None
    servo_slew_rate_deg_s: Optional[float] = None
    servo_travel_deg: Optional[float] = None
    kaman_flap: KamanFlap = field(default_factory=KamanFlap)


@dataclass(frozen=True)
class AutorotationProperties:
    I_ode_kgm2: Optional[float] = None
    omega_min_rad_s: Optional[float] = None
    omega_eq_rad_s: Optional[float] = None


@dataclass(frozen=True)
class RotorDefinition:
    blade: BladeGeometry
    airfoil: AirfoilProperties
    control: Optional[ControlProperties] = None
    inertia: InertiaProperties = field(default_factory=InertiaProperties)
    autorotation: AutorotationProperties = field(default_factory=AutorotationProperties)
    name: str = ""
    description: str = ""

    # Forwarding properties for convenience
    @property
    def span_m(self) -> float:
        return self.blade.span_m

    @property
    def r_cp_m(self) -> float:
        return self.blade.r_cp_m

    @property
    def S_w_m2(self) -> float:
        return self.blade.S_w_m2

    @property
    def disk_area_m2(self) -> float:
        return self.blade.disk_area_m2

    @property
    def aspect_ratio(self) -> float:
        return self.blade.aspect_ratio

    @property
    def solidity(self) -> float:
        return self.blade.solidity

    def validate(self) -> list[ValidationIssue]:
        return self.blade.validate() + self.airfoil.validate()


def _require_yaml() -> None:
    if yaml is None:
        raise ImportError("PyYAML is required to load rotor definition files.")


def load(path: str) -> RotorDefinition:
    _require_yaml()

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
            taper_ratio=float(rotor.get("taper_ratio", 1.0)),
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
            mass_kg=_maybe_float(inertia.get("mass_kg")),
            I_body_kgm2=tuple(float(v) for v in inertia.get("I_body_kgm2", [])),
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


def _maybe_float(value: object) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, dict):
        return float(sum(float(v) for v in value.values()))
    return float(value)


def _maybe_int(value: object) -> Optional[int]:
    return None if value is None else int(value)


def _resolve_polar_csv(value: object, yaml_path: Path) -> Optional[str]:
    if value is None:
        return None
    p = Path(str(value))
    if not p.is_absolute():
        p = yaml_path.parent / p
    return str(p)
