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

    def to_array(self) -> np.ndarray:
        return np.array([self.omega_rad_s, self.spin_angle_rad], dtype=np.float64)

    def from_array(self, arr: np.ndarray) -> "QuasiStaticRotorState":
        if arr.shape != (2,):
            raise ValueError(f"QuasiStaticRotorState expects 2 states, got {arr.shape}")
        return QuasiStaticRotorState(float(arr[0]), float(arr[1]))


@dataclass
class OyeRotorState(RotorState):
    """Øye 2-stage dynamic inflow plus mechanical states.

    Each radial annulus i has its own pair of filter states (W_int_i, W_i)
    representing induced inflow ratios v_i/(Ω·R).  W is what the blade
    actually sees; W_int is the intermediate filter stage between the
    quasi-steady target W_qs and the dynamic W.  Both arrays have length
    ``n_elements`` (set by the rotor's BladeGeometry).

    Inflow ODE (per annulus, Øye 1990; OpenFAST AD Theory §6.3.4):

        τ₁ · dW_int/dt + W_int = W_qs + k · τ₁ · dW_qs/dt
        τ₂ · dW/dt     + W     = W_int

    with k=0.6 empirical coupling and τ₁(a), τ₂(r, a) per OpenFAST.

    Why this instead of Pitt-Peters: each annulus is INDEPENDENT — no
    global L-matrix coupling, no λ_c/λ_s harmonic states, no
    BEM-driven feedback through hub moments.  Numerically much friendlier
    at high advance ratios / descent regimes.

    Array layout (length 2·n_elements + 2):
        [W_int_0..W_int_{n-1}, W_0..W_{n-1}, omega_rad_s, spin_angle_rad]

    The trailing two states match every other rotor state for the
    envelope integrator's convention (arr[-2]=ω, arr[-1]=ψ).
    """

    W_int:          np.ndarray
    W:              np.ndarray
    omega_rad_s:    float = 0.0
    spin_angle_rad: float = 0.0

    def to_array(self) -> np.ndarray:
        return np.concatenate([
            np.asarray(self.W_int, dtype=np.float64),
            np.asarray(self.W,     dtype=np.float64),
            np.array([self.omega_rad_s, self.spin_angle_rad], dtype=np.float64),
        ])

    def from_array(self, arr: np.ndarray) -> "OyeRotorState":
        n_total = arr.shape[0]
        if n_total < 4 or (n_total - 2) % 2 != 0:
            raise ValueError(
                f"OyeRotorState array length {n_total} invalid; expected "
                f"2*n_elements + 2 for some n_elements >= 1."
            )
        n = (n_total - 2) // 2
        return OyeRotorState(
            W_int=np.asarray(arr[:n], dtype=np.float64),
            W=np.asarray(arr[n:2*n], dtype=np.float64),
            omega_rad_s=float(arr[-2]),
            spin_angle_rad=float(arr[-1]),
        )

    @classmethod
    def zeros(cls, n_elements: int, omega_rad_s: float = 0.0
              ) -> "OyeRotorState":
        """Convenience: zero inflow at given rotor speed."""
        return cls(
            W_int=np.zeros(n_elements, dtype=np.float64),
            W=np.zeros(n_elements, dtype=np.float64),
            omega_rad_s=omega_rad_s,
        )


@dataclass
class PittPetersRotorState(RotorState):
    """Pitt-Peters three-state dynamic inflow plus mechanical states.

    Inflow harmonic decomposition (see CLAUDE.md "Rotor rotation direction"
    and "Pitt-Peters inflow ODE"):

        λ(r, ψ) = λ_0 + (r/R) · (λ_c · cos(ψ) + λ_s · sin(ψ))

    where ψ = 0 is at hub +X (nose) and ψ increases CCW from above
    (American convention). All inflow states are non-dimensional
    ratios v/(Ω·R).

    Inflow states (indices 0-2)
    ---------------------------
    lambda_0  uniform inflow through the disk (positive = downwash)
    lambda_c  cosine harmonic of inflow. Negative for forward flight
              along +X (Glauert wake-skew: more inflow at the back of
              disk = ψ = π = tail). Driven by C_M_hub (pitch moment)
              and by C_T cross-coupling in forward flight.
    lambda_s  sine harmonic of inflow. Positive when more inflow is
              on the −Y (left/West) side of disk, e.g. under roll-right
              cyclic input. Driven by C_L_hub (roll moment).

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
