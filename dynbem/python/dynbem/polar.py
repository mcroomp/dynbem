"""dynbem.polar submodule (compat shim)."""
from . import LinearPolar, TabulatedPolar

# Legacy alias: the original pure-Python package exposed an AirfoilPolar
# base type. The Rust
# port treats LinearPolar and TabulatedPolar as the two concrete polar
# classes (any duck-typed object with .cl_cd(alpha) -> (cl, cd) works);
# AirfoilPolar is kept here as a marker for isinstance() checks.
AirfoilPolar = (LinearPolar, TabulatedPolar)

__all__ = ["LinearPolar", "TabulatedPolar", "AirfoilPolar"]
