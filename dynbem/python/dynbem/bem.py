"""dynbem.bem submodule (compat shim).

Re-exports BEMModel, the Prandtl loss helpers, and the per-annulus
solve_bem_element (with BEMElementResult) so legacy dotted-path imports
continue to work.
"""
from . import BEMModel, prandtl_hub_loss, prandtl_tip_loss
from ._dynbem import BEMElementResult, solve_bem_element

__all__ = [
    "BEMModel",
    "BEMElementResult",
    "prandtl_hub_loss",
    "prandtl_tip_loss",
    "solve_bem_element",
]
