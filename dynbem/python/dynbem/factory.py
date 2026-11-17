"""Factory for the three aero models.

Mirrors dynbem.create_aero -- one entry point that builds the right
model + polar from a RotorDefinition. The polar is auto-built from
LinearPolarParameters (LinearPolar from CL0/CL_alpha/CD0/alpha_stall, or
TabulatedPolar from polar_csv) so callers don't have to construct one
manually.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import TYPE_CHECKING, TypeAlias, TypedDict

import numpy as np

from ._dynbem import (
    LinearPolar,
    LinearPolarParameters as RustLinearPolarParameters,
    TabulatedPolar,
    RotorDefinition as RustRotorDefinition,
)
from .rotor_definition import LinearPolarParameters, RotorDefinition

if TYPE_CHECKING:
    from ._dynbem import (
        _OyeBEMModelLinear,  # pyright: ignore[reportPrivateUsage]
        _OyeBEMModelTabulated,  # pyright: ignore[reportPrivateUsage]
        _PittPetersModelLinear,  # pyright: ignore[reportPrivateUsage]
        _PittPetersModelTabulated,  # pyright: ignore[reportPrivateUsage]
        _QuasiStaticBEMLinear,  # pyright: ignore[reportPrivateUsage]
        _QuasiStaticBEMTabulated,  # pyright: ignore[reportPrivateUsage]
        _VpmRotorLinear,  # pyright: ignore[reportPrivateUsage]
        _VpmRotorTabulated,  # pyright: ignore[reportPrivateUsage]
    )

    QuasiStaticModel: TypeAlias = _QuasiStaticBEMLinear | _QuasiStaticBEMTabulated
    PittPetersModel: TypeAlias = _PittPetersModelLinear | _PittPetersModelTabulated
    OyeBEMModel: TypeAlias = _OyeBEMModelLinear | _OyeBEMModelTabulated
    VpmRotor: TypeAlias = _VpmRotorLinear | _VpmRotorTabulated
    AeroModel: TypeAlias = (
        QuasiStaticModel | PittPetersModel | OyeBEMModel | VpmRotor
    )


class VpmConfig(TypedDict, total=False):
    max_particles: int
    sigma: float
    relax: float
    nonlinear_lifting_line: bool
    tip_clustering: bool
    local_core: bool
    barnes_hut: bool
    bh_theta: float
    bh_min_particles: int
    merge_wake: bool
    merge_every: int
    merge_min_particles: int
    merge_kappa: float
    merge_chi_min: float
    merge_region_dist: float
    core_spread_nu: float
    strength_decay_tau_rev: float


def build_polar(
    airfoil: LinearPolarParameters | RustLinearPolarParameters,
) -> LinearPolar | TabulatedPolar:
    """Build a polar from a LinearPolarParameters.

    If polar_csv is provided, load it (airfoiltools.com format -- 9
    metadata rows then Alpha,Cl,Cd, with any non-numeric leading rows
    skipped). Otherwise build a LinearPolar from CL0/CL_alpha/CD0/stall.
    """
    if isinstance(airfoil, LinearPolarParameters) and airfoil.polar_csv is not None:
        return load_tabulated_polar(airfoil.polar_csv)
    return LinearPolar(
        CL0=airfoil.CL0,
        CL_alpha_per_rad=airfoil.CL_alpha_per_rad,
        CD0=airfoil.CD0,
        alpha_stall_rad=math.radians(airfoil.alpha_stall_deg),
    )


def load_tabulated_polar(path: str | Path) -> TabulatedPolar:
    """Load TabulatedPolar from an Alpha,Cl,Cd CSV (airfoiltools.com format)."""
    alphas: list[float] = []
    cls: list[float] = []
    cds: list[float] = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            try:
                alphas.append(math.radians(float(parts[0])))
                cls.append(float(parts[1]))
                cds.append(float(parts[2]))
            except (ValueError, IndexError):
                continue  # skip metadata/header rows
    if not alphas:
        raise ValueError(f"No numeric data found in polar CSV: {path}")
    return TabulatedPolar(
        alpha_rad=np.array(alphas, dtype=float),
        cl=np.array(cls, dtype=float),
        cd=np.array(cds, dtype=float),
    )


def _build_polar_from_defn(  # pyright: ignore[reportUnusedFunction]
    defn: RotorDefinition | RustRotorDefinition,
    polar: LinearPolar | TabulatedPolar | None,
) -> LinearPolar | TabulatedPolar:
    """Return ``polar`` if given, else auto-build from ``defn.airfoil``."""
    if polar is not None:
        return polar
    return build_polar(defn.airfoil)


def create_aero(
    defn: RotorDefinition | RustRotorDefinition,
    model: str,
    *,
    n_psi_elements: int = 36,
    polar: LinearPolar | TabulatedPolar | None = None,
) -> AeroModel:
    """Build an aero model for a RotorDefinition.

    model
      "quasi_static" / "bem"   QuasiStaticBEM (Level 1, no dynamic inflow)
      "pitt_peters"            PittPetersModel (Level 2, Pitt-Peters L-matrix)
      "oye" / "oye_bem"        OyeBEMModel (Level 2, Oye 2-stage annular inflow)
      "vpm" / "vpm_rotor" / "free_wake"
                               VpmRotor (Level 3, forward-flight free-wake VPM;
                               no single-shot compute_forces -- advance with
                               step(inputs, state, dt))
    """
    if polar is None:
        polar = build_polar(defn.airfoil)
    # Lazy import: the public QuasiStaticBEM/PittPetersModel/OyeBEMModel are
    # Python subclasses defined in dynbem/__init__.py after this module is
    # imported, so we resolve them at call time rather than module-load time.
    from . import OyeBEMModel, PittPetersModel, QuasiStaticBEM  # noqa: WPS433
    if model in ("quasi_static", "bem"):
        return QuasiStaticBEM(defn, polar, n_psi_elements)
    if model in ("pitt_peters",):
        return PittPetersModel(defn, polar, n_psi_elements)
    if model in ("oye", "oye_bem"):
        return OyeBEMModel(defn, polar, n_psi_elements, 0.6)
    if model in ("vpm", "vpm_rotor", "free_wake"):
        from . import VpmRotor  # noqa: WPS433
        return VpmRotor(defn, polar)
    raise ValueError(
        f"Unknown aero model {model!r}. "
        "Choose 'quasi_static' (alias 'bem'), 'pitt_peters', 'oye', or 'vpm'."
    )


__all__ = ["create_aero", "build_polar", "load_tabulated_polar"]
