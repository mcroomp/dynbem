"""Airfoil polar abstractions."""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path

import numpy as np


class AirfoilPolar(ABC):
    """Maps angle of attack to (CL, CD) for a 2D airfoil section."""

    @abstractmethod
    def cl_cd(self, alpha_rad: float) -> tuple[float, float]: ...


@dataclass(frozen=True)
class LinearPolar(AirfoilPolar):
    """Linear lift curve with constant drag, clipped at stall.

    Below stall:   CL = CL0 + CL_alpha * alpha,  CD = CD0
    At/above stall: CL is capped at the stall value; CD grows linearly past it.
    """

    CL0: float
    CL_alpha_per_rad: float
    CD0: float
    alpha_stall_rad: float

    def cl_cd(self, alpha_rad: float) -> tuple[float, float]:
        if abs(alpha_rad) < self.alpha_stall_rad:
            return self.CL0 + self.CL_alpha_per_rad * alpha_rad, self.CD0
        cl = math.copysign(
            self.CL0 + self.CL_alpha_per_rad * self.alpha_stall_rad, alpha_rad
        )
        cd = self.CD0 + (abs(alpha_rad) - self.alpha_stall_rad)
        return cl, cd

    @classmethod
    def from_properties(cls, props: "AirfoilProperties") -> "LinearPolar":
        return cls(
            CL0=props.CL0,
            CL_alpha_per_rad=props.CL_alpha_per_rad,
            CD0=props.CD0,
            alpha_stall_rad=math.radians(props.alpha_stall_deg),
        )


class TabulatedPolar(AirfoilPolar):
    """Interpolated polar from a tabulated (alpha, CL, CD) dataset.

    Reads the airfoiltools.com CSV format (9 metadata lines, then a header
    line, then rows of: Alpha, Cl, Cd, ...).  Any non-numeric leading lines
    are skipped automatically.

    Outside the tabulated alpha range the nearest endpoint values are used
    (clamp, not extrapolate) to avoid divergence in the BEM iteration.
    """

    def __init__(self, alpha_rad: np.ndarray, cl: np.ndarray, cd: np.ndarray) -> None:
        self._alpha = alpha_rad
        self._cl = cl
        self._cd = cd

    def cl_cd(self, alpha_rad: float) -> tuple[float, float]:
        cl = float(np.interp(alpha_rad, self._alpha, self._cl))
        cd = float(np.interp(alpha_rad, self._alpha, self._cd))
        return cl, cd

    @classmethod
    def from_csv(cls, path: str | Path) -> "TabulatedPolar":
        """Load from an airfoiltools.com polar CSV (or any Alpha,Cl,Cd CSV)."""
        alphas, cls_, cds = [], [], []
        with open(path, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split(",")
                try:
                    alphas.append(math.radians(float(parts[0])))
                    cls_.append(float(parts[1]))
                    cds.append(float(parts[2]))
                except (ValueError, IndexError):
                    continue  # skip header / metadata rows
        if not alphas:
            raise ValueError(f"No numeric data found in polar CSV: {path}")
        return cls(np.array(alphas), np.array(cls_), np.array(cds))
