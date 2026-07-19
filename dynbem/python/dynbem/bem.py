"""dynbem.bem submodule (compat shim).

Re-exports QuasiStaticBEM (legacy alias BEMModel) so legacy dotted-path
imports continue to work.
"""
from . import BEMModel, QuasiStaticBEM

__all__ = [
    "QuasiStaticBEM",
    "BEMModel",
]
