"""dynbem.bem submodule (compat shim).

Re-exports QuasiStaticBEM (legacy alias BEMModel) and the Prandtl loss
helpers so legacy dotted-path imports continue to work.
"""
from . import BEMModel, QuasiStaticBEM, prandtl_hub_loss, prandtl_tip_loss

__all__ = [
    "QuasiStaticBEM",
    "BEMModel",
    "prandtl_hub_loss",
    "prandtl_tip_loss",
]
