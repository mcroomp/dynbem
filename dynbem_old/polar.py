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

    def cl_cd_arr(self, alpha_rad: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Batched cl_cd. Default falls back to scalar loop; override for speed."""
        a = np.asarray(alpha_rad, dtype=float).ravel()
        cl = np.empty_like(a)
        cd = np.empty_like(a)
        for i, x in enumerate(a):
            cl[i], cd[i] = self.cl_cd(float(x))
        return cl.reshape(np.shape(alpha_rad)), cd.reshape(np.shape(alpha_rad))


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

    def cl_cd_arr(self, alpha_rad: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Vectorized cl_cd over an array of angles of attack."""
        a = np.asarray(alpha_rad, dtype=float)
        cl_lin = self.CL0 + self.CL_alpha_per_rad * a
        cd = np.full_like(a, self.CD0)
        stall = np.abs(a) >= self.alpha_stall_rad
        if stall.any():
            cl_stall_mag = self.CL0 + self.CL_alpha_per_rad * self.alpha_stall_rad
            cl_stall = np.copysign(cl_stall_mag, a)
            cl = np.where(stall, cl_stall, cl_lin)
            cd = np.where(stall, self.CD0 + (np.abs(a) - self.alpha_stall_rad), cd)
            return cl, cd
        return cl_lin, cd

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
        self._alpha = np.ascontiguousarray(alpha_rad, dtype=float)
        self._cl = np.ascontiguousarray(cl, dtype=float)
        self._cd = np.ascontiguousarray(cd, dtype=float)

    def cl_cd(self, alpha_rad: float) -> tuple[float, float]:
        cl = float(np.interp(alpha_rad, self._alpha, self._cl))
        cd = float(np.interp(alpha_rad, self._alpha, self._cd))
        return cl, cd

    def cl_cd_arr(self, alpha_rad: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Vectorized cl_cd via np.interp (its C binary search beats a
        hand-rolled index lookup at the array sizes we use, ~10 elements)."""
        cl = np.interp(alpha_rad, self._alpha, self._cl)
        cd = np.interp(alpha_rad, self._alpha, self._cd)
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
