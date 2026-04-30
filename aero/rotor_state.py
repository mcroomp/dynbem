"""Rotor dynamic state abstractions for ODE integration."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np


class RotorState(ABC):
    """All rotor-level states owned and integrated by the external ODE solver.

    Every concrete subclass must expose omega_rad_s and spin_angle_rad as
    plain fields alongside whatever inflow states the model requires.

    When returned from compute_forces() the array represents dstate/dt
    (derivatives). The caller integrates:

        next_state = state.from_array(
            state.to_array() + dt * derivative.to_array()
        )
    """

    @property
    @abstractmethod
    def n_states(self) -> int:
        """Total number of scalar states in this model."""
        ...

    @abstractmethod
    def to_array(self) -> np.ndarray:
        """Pack state into a 1-D float64 array for the ODE integrator."""
        ...

    @abstractmethod
    def from_array(self, arr: np.ndarray) -> "RotorState":
        """Return a new instance reconstructed from an integrated array slice."""
        ...


@dataclass
class QuasiStaticRotorState(RotorState):
    """Quasi-static inflow — BEM converges to steady state each step.

    Carries only the two mechanical states. Acts as the no-op inflow
    model so that quasi-static and dynamic implementations satisfy the
    same interface.

    Fields
    ------
    omega_rad_s     rotor angular velocity, rad/s
    spin_angle_rad  rotor azimuth angle, rad
    """

    omega_rad_s:    float = 0.0
    spin_angle_rad: float = 0.0

    @property
    def n_states(self) -> int:
        return 2

    def to_array(self) -> np.ndarray:
        return np.array([self.omega_rad_s, self.spin_angle_rad], dtype=np.float64)

    def from_array(self, arr: np.ndarray) -> "QuasiStaticRotorState":
        if arr.shape != (2,):
            raise ValueError(f"QuasiStaticRotorState expects 2 states, got {arr.shape}")
        return QuasiStaticRotorState(float(arr[0]), float(arr[1]))


@dataclass
class PittPetersRotorState(RotorState):
    """Pitt-Peters three-state dynamic inflow plus mechanical states.

    Inflow states (indices 0-2)
    ---------------------------
    lambda_0  uniform inflow through the disk
    lambda_c  fore-aft tilt of the inflow distribution (cosine harmonic)
    lambda_s  lateral tilt of the inflow distribution (sine harmonic)

    Mechanical states (indices 3-4)
    --------------------------------
    omega_rad_s     rotor angular velocity, rad/s
    spin_angle_rad  rotor azimuth angle, rad
    """

    lambda_0:       float = 0.0
    lambda_c:       float = 0.0
    lambda_s:       float = 0.0
    omega_rad_s:    float = 0.0
    spin_angle_rad: float = 0.0

    @property
    def n_states(self) -> int:
        return 5

    def to_array(self) -> np.ndarray:
        return np.array(
            [self.lambda_0, self.lambda_c, self.lambda_s,
             self.omega_rad_s, self.spin_angle_rad],
            dtype=np.float64,
        )

    def from_array(self, arr: np.ndarray) -> "PittPetersRotorState":
        if arr.shape != (5,):
            raise ValueError(f"PittPetersRotorState expects 5 states, got {arr.shape}")
        return PittPetersRotorState(
            float(arr[0]), float(arr[1]), float(arr[2]),
            float(arr[3]), float(arr[4]),
        )
