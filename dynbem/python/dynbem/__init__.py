"""dynbem - Rust BEM / Pitt-Peters / Oye dynamic-inflow rotor models.

Maturin builds the compiled extension `_dynbem`; this shim re-exports
the public surface. (Formerly published as `dynbem_rs`; the original
pure-Python implementation lives on as `dynbem_old` for reference.)
"""
from abc import ABC

from ._dynbem import (  # noqa: F401
    vrs_lambda1,
    cyclic_coeffs as _cyclic_coeffs_rust,
    prandtl_tip_loss,
    prandtl_hub_loss,
    LinearPolar,
    TabulatedPolar,
    QuasiStaticRotorState,
    PittPetersRotorState,
    OyeRotorState,
    RotorInputs,
    AeroResult,
    TrimResult,
)
from .rotor_definition import (  # noqa: F401
    BladeGeometry,
    KamanFlap,
    InertiaProperties,
    AirfoilProperties,
    ControlProperties,
    AutorotationProperties,
    RotorDefinition,
)
from ._dynbem import QuasiStaticBEM as _QuasiStaticBEM  # noqa: F401
from ._dynbem import PittPetersModel as _PittPetersModel  # noqa: F401
from ._dynbem import OyeBEMModel as _OyeBEMModel  # noqa: F401
from .factory import create_aero, build_polar, load_tabulated_polar  # noqa: F401
from .trim import solve_trim_cyclic, relax_inflow  # noqa: F401
from .mechanical import omega_derivative, euler_step_omega  # noqa: F401


def cyclic_coeffs(tilt_lon, tilt_lat, control=None):  # noqa: F401
    """Compute (theta_1c, theta_1s) cyclic coefficients.

    Accepts either a _dynbem.ControlProperties (lean Rust class) or a
    dynbem.ControlProperties Python wrapper; in the latter case the internal
    ``._rust`` handle is forwarded to the Rust implementation.
    """
    rust_ctrl = getattr(control, "_rust", control)
    return _cyclic_coeffs_rust(tilt_lon, tilt_lat, rust_ctrl)


# ---------------------------------------------------------------------------
# Python subclasses of the Rust model pyclasses that auto-build a polar
# from the rotor's AirfoilProperties when none is given and extract the
# lean ._rust RotorDefinition from the Python wrapper before passing to Rust.
# ---------------------------------------------------------------------------

def _rust_defn(defn):
    """Extract the lean _dynbem.RotorDefinition from a Python wrapper."""
    return getattr(defn, "_rust", defn)


class QuasiStaticBEM(_QuasiStaticBEM):
    def __new__(cls, defn, polar, n_psi_elements):
        return _QuasiStaticBEM.__new__(cls, _rust_defn(defn), polar, n_psi_elements)

    def __init__(self, defn, polar, n_psi_elements):
        self._defn = defn

    @property
    def defn(self):
        return self._defn


# Backwards-compat alias. The model used to be named `BEMModel`, but
# "BEM" is the family that PittPetersModel and OyeBEMModel also belong
# to -- the distinguishing feature of this one is quasi-static inflow.
# Existing callers and YAML configs that reference `BEMModel` keep working.
BEMModel = QuasiStaticBEM


class PittPetersModel(_PittPetersModel):
    def __new__(cls, defn, polar, n_psi_elements):
        return _PittPetersModel.__new__(cls, _rust_defn(defn), polar, n_psi_elements)

    def __init__(self, defn, polar, n_psi_elements):
        self._defn = defn

    @property
    def defn(self):
        return self._defn


class OyeBEMModel(_OyeBEMModel):
    def __new__(cls, defn, polar, n_psi_elements, coupling_k):
        return _OyeBEMModel.__new__(
            cls, _rust_defn(defn), polar, n_psi_elements, coupling_k,
        )

    def __init__(self, defn, polar, n_psi_elements, coupling_k):
        self._defn = defn

    @property
    def defn(self):
        return self._defn


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
    (QuasiStaticBEM, PittPetersModel, OyeBEMModel) are registered at import."""


AeroBase.register(QuasiStaticBEM)
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
    "QuasiStaticBEM", "BEMModel", "PittPetersModel", "OyeBEMModel",
    "TrimResult",
    # virtual ABCs
    "AeroBase", "RotorState",
]
