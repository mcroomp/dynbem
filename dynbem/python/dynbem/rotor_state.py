"""dynbem.rotor_state submodule (compat shim).

Re-exports the three concrete state types from the Rust extension and
declares a virtual base class `RotorState` so that
``isinstance(state, RotorState)`` keeps working in legacy code.

The Rust types are #[pyclass]es that don't inherit from a Python ABC,
so we register them explicitly here.
"""
from abc import ABC

from . import OyeRotorState, PittPetersRotorState, QuasiStaticRotorState


class RotorState(ABC):
    """Virtual base for rotor state vectors. The three concrete types
    (QuasiStaticRotorState, PittPetersRotorState, OyeRotorState) are
    registered with this ABC at import time so isinstance() works."""


RotorState.register(QuasiStaticRotorState)
RotorState.register(PittPetersRotorState)
RotorState.register(OyeRotorState)


__all__ = [
    "RotorState",
    "QuasiStaticRotorState",
    "PittPetersRotorState",
    "OyeRotorState",
]
