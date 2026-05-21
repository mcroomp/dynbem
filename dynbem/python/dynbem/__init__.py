"""dynbem - Rust BEM / Pitt-Peters / Oye dynamic-inflow rotor models.

Maturin builds the compiled extension `_dynbem`; this shim re-exports
the public surface. (Formerly published as `dynbem_rs`; the original
pure-Python implementation lives on as `dynbem_old` for reference.)
"""
from abc import ABC

from ._dynbem import (  # noqa: F401
    vrs_lambda1,
    cyclic_coeffs,
    prandtl_tip_loss,
    prandtl_hub_loss,
    LinearPolar,
    TabulatedPolar,
    KamanFlap,
    BladeGeometry,
    AirfoilProperties,
    InertiaProperties,
    ControlProperties,
    AutorotationProperties,
    RotorDefinition,
    QuasiStaticRotorState,
    PittPetersRotorState,
    OyeRotorState,
    RotorInputs,
    AeroResult,
    TrimResult,
)
from ._dynbem import BEMModel as _BEMModel  # noqa: F401
from ._dynbem import PittPetersModel as _PittPetersModel  # noqa: F401
from ._dynbem import OyeBEMModel as _OyeBEMModel  # noqa: F401
from .factory import create_aero, build_polar, load_tabulated_polar  # noqa: F401
from .trim import solve_trim_cyclic, relax_inflow  # noqa: F401


# ---------------------------------------------------------------------------
# Python subclasses of the Rust model pyclasses that auto-build a polar
# from the rotor's AirfoilProperties when none is given. This restores the
# legacy `BEMModel(defn=...)` ergonomics, including loading the tabulated
# polar from the YAML's polar_csv (which Rust's auto-default LinearPolar
# can't do because it doesn't read files).
# ---------------------------------------------------------------------------

class BEMModel(_BEMModel):
    def __new__(cls, defn, polar=None, n_psi_elements=36):
        if polar is None:
            polar = build_polar(defn.airfoil)
        return _BEMModel.__new__(cls, defn, polar, n_psi_elements)


class PittPetersModel(_PittPetersModel):
    def __new__(cls, defn, polar=None, n_psi_elements=36):
        if polar is None:
            polar = build_polar(defn.airfoil)
        return _PittPetersModel.__new__(cls, defn, polar, n_psi_elements)


class OyeBEMModel(_OyeBEMModel):
    def __new__(cls, defn, polar=None, n_psi_elements=36, coupling_k=0.6):
        if polar is None:
            polar = build_polar(defn.airfoil)
        return _OyeBEMModel.__new__(
            cls, defn, polar, n_psi_elements, coupling_k,
        )


# ---------------------------------------------------------------------------
# Virtual ABCs for isinstance() compatibility with the legacy Python API.
#
# The Rust pyclasses don't inherit from Python ABCs, so we declare two
# marker base classes and register the concrete types with them.
# ---------------------------------------------------------------------------

class RotorState(ABC):
    """Virtual base for rotor state vectors. The three concrete types
    (QuasiStaticRotorState, PittPetersRotorState, OyeRotorState) are
    registered with this ABC at import time."""


RotorState.register(QuasiStaticRotorState)
RotorState.register(PittPetersRotorState)
RotorState.register(OyeRotorState)


class AeroBase(ABC):
    """Virtual base for aero models. The three concrete model classes
    (BEMModel, PittPetersModel, OyeBEMModel) are registered at import."""


AeroBase.register(BEMModel)
AeroBase.register(PittPetersModel)
AeroBase.register(OyeBEMModel)


# Legacy alias: dynbem_old exposed `AirfoilPolar` as a marker for the
# two concrete polar types. Provide it as a tuple usable in isinstance().
AirfoilPolar = (LinearPolar, TabulatedPolar)


__all__ = [
    # functions
    "vrs_lambda1", "cyclic_coeffs",
    "prandtl_tip_loss", "prandtl_hub_loss",
    "create_aero", "build_polar", "load_tabulated_polar",
    "solve_trim_cyclic", "relax_inflow",
    # types
    "LinearPolar", "TabulatedPolar", "AirfoilPolar",
    "KamanFlap", "BladeGeometry", "AirfoilProperties",
    "InertiaProperties", "ControlProperties", "AutorotationProperties",
    "RotorDefinition",
    "QuasiStaticRotorState", "PittPetersRotorState", "OyeRotorState",
    "RotorInputs", "AeroResult",
    "BEMModel", "PittPetersModel", "OyeBEMModel",
    "TrimResult",
    # virtual ABCs
    "AeroBase", "RotorState",
]
