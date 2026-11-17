"""dynbem - Rust BEM / Pitt-Peters / Oye dynamic-inflow rotor models.

Maturin builds the compiled extension `_dynbem`; this shim re-exports
the public surface. (Formerly published as `dynbem_rs`; supersedes the
legacy pure-Python implementation.)
"""

from __future__ import annotations

from typing import TYPE_CHECKING, TypeAlias, Unpack

from ._dynbem import (  # noqa: F401
    cyclic_coeffs as _cyclic_coeffs_rust,
    LinearPolar,
    TabulatedPolar,
    QuasiStaticRotorState,
    PittPetersRotorState,
    OyeRotorState,
    VpmRotorState,
    RotorInputs,
    AeroResult,
    TrimResult,
    ControlProperties as _RustControlProperties,
)
from ._dynbem import _QuasiStaticBEMLinear, _QuasiStaticBEMTabulated  # noqa: F401  # pyright: ignore[reportPrivateUsage]
from ._dynbem import _PittPetersModelLinear, _PittPetersModelTabulated  # noqa: F401  # pyright: ignore[reportPrivateUsage]
from ._dynbem import _OyeBEMModelLinear, _OyeBEMModelTabulated  # noqa: F401  # pyright: ignore[reportPrivateUsage]
from ._dynbem import _VpmRotorLinear, _VpmRotorTabulated  # noqa: F401  # pyright: ignore[reportPrivateUsage]
from .rotor_definition import (  # noqa: F401
    BladeGeometry,
    InertiaProperties,
    LinearPolarParameters,
    ControlProperties,
    AutorotationProperties,
    RotorDefinition,
)
from .rotor_definition import _to_rust_control, _to_rust_defn  # pyright: ignore[reportPrivateUsage]
from .factory import build_polar, create_aero, load_tabulated_polar
from .factory import _build_polar_from_defn  # pyright: ignore[reportPrivateUsage]
from .trim import solve_trim_cyclic, relax_inflow  # noqa: F401
from .mechanical import omega_derivative, step_omega  # noqa: F401

if TYPE_CHECKING:
    from .factory import (
        OyeBEMModel as _OyeBEMModel,
        PittPetersModel as _PittPetersModel,
        QuasiStaticModel as _QuasiStaticModel,
        VpmRotor as _VpmRotor,
    )
    from .factory import VpmConfig as _VpmConfig
    from ._dynbem import RotorDefinition as _RustRotorDefinition

Polar: TypeAlias = LinearPolar | TabulatedPolar


def cyclic_coeffs(
    tilt_lon: float,
    tilt_lat: float,
    control: ControlProperties | _RustControlProperties | None = None,
) -> tuple[float, float]:
    """Compute (theta_1c, theta_1s) cyclic coefficients."""
    return _cyclic_coeffs_rust(tilt_lon, tilt_lat, _to_rust_control(control))


# ---------------------------------------------------------------------------
# Factory functions — build the concrete Rust model depending on polar type.
# These replace the old Python subclasses (which needed __new__ tricks to
# wrap PyO3 pyclasses). Duck-typed: isinstance checks work via the Rust
# pyclasses themselves (_QuasiStaticBEMLinear / _QuasiStaticBEMTabulated etc.)
# ---------------------------------------------------------------------------


def QuasiStaticBEM(
    defn: RotorDefinition | _RustRotorDefinition,
    polar: Polar | None = None,
    n_psi_elements: int = 36,
) -> _QuasiStaticModel:  # noqa: N802
    polar_rs = _build_polar_from_defn(defn, polar)
    defn_rs = _to_rust_defn(defn)
    if isinstance(polar_rs, LinearPolar):
        return _QuasiStaticBEMLinear(defn_rs, polar_rs, n_psi_elements)
    return _QuasiStaticBEMTabulated(defn_rs, polar_rs, n_psi_elements)


def PittPetersModel(
    defn: RotorDefinition | _RustRotorDefinition,
    polar: Polar | None = None,
    n_psi_elements: int = 36,
) -> _PittPetersModel:  # noqa: N802
    polar_rs = _build_polar_from_defn(defn, polar)
    defn_rs = _to_rust_defn(defn)
    if isinstance(polar_rs, LinearPolar):
        return _PittPetersModelLinear(defn_rs, polar_rs, n_psi_elements)
    return _PittPetersModelTabulated(defn_rs, polar_rs, n_psi_elements)


def OyeBEMModel(
    defn: RotorDefinition | _RustRotorDefinition,
    polar: Polar | None = None,
    n_psi_elements: int = 36,
    coupling_k: float = 0.6,
) -> _OyeBEMModel:  # noqa: N802
    polar_rs = _build_polar_from_defn(defn, polar)
    defn_rs = _to_rust_defn(defn)
    if isinstance(polar_rs, LinearPolar):
        return _OyeBEMModelLinear(defn_rs, polar_rs, n_psi_elements, coupling_k)
    return _OyeBEMModelTabulated(defn_rs, polar_rs, n_psi_elements, coupling_k)


def VpmRotor(
    defn: RotorDefinition | _RustRotorDefinition,
    polar: Polar | None = None,
    **config: Unpack[_VpmConfig],
) -> _VpmRotor:  # noqa: N802
    """Forward-flight free-wake VPM rotor (Level 3).

    Like the BEM-family models, the Rust core is generic over the polar type,
    so this factory dispatches to _VpmRotorLinear or _VpmRotorTabulated based on
    the polar instance. Has no single-shot compute_forces: advance it in time
    with step(inputs, state, dt) -- one sub-step of duration dt per call.

    Config kwargs: max_particles (default 4800), sigma, relax,
    nonlinear_lifting_line, tip_clustering, local_core, barnes_hut,
    bh_theta, bh_min_particles.
    """
    polar_rs = _build_polar_from_defn(defn, polar)
    defn_rs = _to_rust_defn(defn)
    if isinstance(polar_rs, LinearPolar):
        return _VpmRotorLinear(defn_rs, polar_rs, **config)
    return _VpmRotorTabulated(defn_rs, polar_rs, **config)

__all__ = [
    # functions
    "cyclic_coeffs",
    "create_aero", "build_polar", "load_tabulated_polar",
    "solve_trim_cyclic", "relax_inflow",
    "omega_derivative", "step_omega",
    # types
    "LinearPolar", "TabulatedPolar",
    "BladeGeometry", "LinearPolarParameters",
    "InertiaProperties", "ControlProperties", "AutorotationProperties",
    "RotorDefinition",
    "QuasiStaticRotorState", "PittPetersRotorState", "OyeRotorState",
    "VpmRotorState",
    "RotorInputs", "AeroResult",
    "QuasiStaticBEM", "PittPetersModel", "OyeBEMModel",
    "VpmRotor",
    "TrimResult",
]

