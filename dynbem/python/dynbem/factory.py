"""Factory for the three aero models.

Mirrors dynbem.create_aero -- one entry point that builds the right
model + polar from a RotorDefinition. The polar is auto-built from
AirfoilProperties (LinearPolar from CL0/CL_alpha/CD0/alpha_stall, or
TabulatedPolar from polar_csv) so callers don't have to construct one
manually.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Optional, Union

import numpy as np

from ._dynbem import (
    LinearPolar,
    TabulatedPolar,
    RotorDefinition,
)


def build_polar(airfoil) -> Union[LinearPolar, TabulatedPolar]:
    """Build a polar from an AirfoilProperties.

    If polar_csv is provided, load it (airfoiltools.com format -- 9
    metadata rows then Alpha,Cl,Cd, with any non-numeric leading rows
    skipped). Otherwise build a LinearPolar from CL0/CL_alpha/CD0/stall.
    """
    if airfoil.polar_csv is not None:
        return load_tabulated_polar(airfoil.polar_csv)
    return LinearPolar(
        CL0=airfoil.CL0,
        CL_alpha_per_rad=airfoil.CL_alpha_per_rad,
        CD0=airfoil.CD0,
        alpha_stall_rad=math.radians(airfoil.alpha_stall_deg),
    )


def load_tabulated_polar(path: Union[str, Path]) -> TabulatedPolar:
    """Load TabulatedPolar from an Alpha,Cl,Cd CSV (airfoiltools.com format)."""
    alphas, cls, cds = [], [], []
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


def create_aero(
    defn: RotorDefinition,
    model: str = "pitt_peters",
    *,
    n_psi_elements: int = 36,
    polar: Optional[Union[LinearPolar, TabulatedPolar]] = None,
):
    """Build an aero model for a RotorDefinition.

    model
      "quasi_static" / "bem"   QuasiStaticBEM (Level 1, no dynamic inflow)
      "pitt_peters"            PittPetersModel (Level 2, Pitt-Peters L-matrix)
      "pitt_peters_jit"        alias for pitt_peters (Rust is already compiled)
      "oye" / "oye_bem"        OyeBEMModel (Level 2, Oye 2-stage annular inflow)
    """
    if polar is None:
        polar = build_polar(defn.airfoil)
    # Lazy import: the public QuasiStaticBEM/PittPetersModel/OyeBEMModel are
    # Python subclasses defined in dynbem/__init__.py after this module is
    # imported, so we resolve them at call time rather than module-load time.
    from . import OyeBEMModel, PittPetersModel, QuasiStaticBEM  # noqa: WPS433
    if model in ("quasi_static", "bem"):
        return QuasiStaticBEM(defn, polar, n_psi_elements=n_psi_elements)
    if model in ("pitt_peters", "pitt_peters_jit", "jit", "pitt_peters_numpy"):
        return PittPetersModel(defn, polar, n_psi_elements=n_psi_elements)
    if model in ("oye", "oye_bem"):
        return OyeBEMModel(defn, polar, n_psi_elements=n_psi_elements)
    raise ValueError(
        f"Unknown aero model {model!r}. "
        "Choose 'quasi_static' (alias 'bem'), 'pitt_peters', or 'oye'."
    )


# silence unused-import warning for csv (kept in case future extensions need it)
_ = csv

__all__ = ["create_aero", "build_polar", "load_tabulated_polar"]
