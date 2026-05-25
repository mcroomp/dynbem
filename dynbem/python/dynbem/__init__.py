"""dynbem - Rust BEM / Pitt-Peters / Oye dynamic-inflow rotor models.

Maturin builds the compiled extension `_dynbem`; this shim re-exports
the public surface. (Formerly published as `dynbem_rs`; the original
pure-Python implementation lives on as `dynbem_old` for reference.)
"""

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
from ._dynbem import _QuasiStaticBEMLinear, _QuasiStaticBEMTabulated  # noqa: F401
from ._dynbem import _PittPetersModelLinear, _PittPetersModelTabulated  # noqa: F401
from ._dynbem import _OyeBEMModelLinear, _OyeBEMModelTabulated  # noqa: F401
from .rotor_definition import (  # noqa: F401
    BladeGeometry,
    KamanFlap,
    InertiaProperties,
    LinearPolarParameters,
    ControlProperties,
    AutorotationProperties,
    RotorDefinition,
    _to_rust_defn,
    _to_rust_control,
)
from .factory import create_aero, build_polar, load_tabulated_polar, _build_polar_from_defn  # noqa: F401
from .trim import solve_trim_cyclic, relax_inflow  # noqa: F401
from .mechanical import omega_derivative, euler_step_omega  # noqa: F401


def cyclic_coeffs(tilt_lon, tilt_lat, control=None):  # noqa: F401
    """Compute (theta_1c, theta_1s) cyclic coefficients."""
    return _cyclic_coeffs_rust(tilt_lon, tilt_lat, _to_rust_control(control))


# ---------------------------------------------------------------------------
# Factory functions — build the concrete Rust model depending on polar type.
# These replace the old Python subclasses (which needed __new__ tricks to
# wrap PyO3 pyclasses). Duck-typed: isinstance checks work via the Rust
# pyclasses themselves (_QuasiStaticBEMLinear / _QuasiStaticBEMTabulated etc.)
# ---------------------------------------------------------------------------


def QuasiStaticBEM(defn, polar=None, n_psi_elements=36):  # noqa: N802
    polar_rs = _build_polar_from_defn(defn, polar)
    defn_rs = _to_rust_defn(defn)
    if isinstance(polar_rs, LinearPolar):
        return _QuasiStaticBEMLinear(defn_rs, polar_rs, n_psi_elements)
    return _QuasiStaticBEMTabulated(defn_rs, polar_rs, n_psi_elements)


# Backwards-compat alias (was BEMModel before the quasi-static qualifier was added).
BEMModel = QuasiStaticBEM


def PittPetersModel(defn, polar=None, n_psi_elements=36):  # noqa: N802
    polar_rs = _build_polar_from_defn(defn, polar)
    defn_rs = _to_rust_defn(defn)
    if isinstance(polar_rs, LinearPolar):
        return _PittPetersModelLinear(defn_rs, polar_rs, n_psi_elements)
    return _PittPetersModelTabulated(defn_rs, polar_rs, n_psi_elements)


def OyeBEMModel(defn, polar=None, n_psi_elements=36, coupling_k=0.6):  # noqa: N802
    polar_rs = _build_polar_from_defn(defn, polar)
    defn_rs = _to_rust_defn(defn)
    if isinstance(polar_rs, LinearPolar):
        return _OyeBEMModelLinear(defn_rs, polar_rs, n_psi_elements, coupling_k)
    return _OyeBEMModelTabulated(defn_rs, polar_rs, n_psi_elements, coupling_k)



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
    "KamanFlap", "BladeGeometry", "LinearPolarParameters",
    "InertiaProperties", "ControlProperties", "AutorotationProperties",
    "RotorDefinition",
    "QuasiStaticRotorState", "PittPetersRotorState", "OyeRotorState",
    "RotorInputs", "AeroResult",
    "QuasiStaticBEM", "BEMModel", "PittPetersModel", "OyeBEMModel",
    "TrimResult",
]

