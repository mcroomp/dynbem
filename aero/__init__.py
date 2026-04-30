"""Aerodynamic interfaces and factory hooks."""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
import numpy as np

from .rotor_state import PittPetersRotorState, QuasiStaticRotorState, RotorState
from .polar import AirfoilPolar, LinearPolar


@dataclass
class AeroResult:
    """Return type for compute_forces()."""
    F_world:   np.ndarray   # [3]
    M_orbital: np.ndarray   # [3]
    Q_spin:    float
    M_spin:    np.ndarray   # [3]

    def _as_array(self) -> np.ndarray:
        return np.concatenate([self.F_world, self.M_orbital + self.M_spin])

    def __getitem__(self, key):
        return self._as_array()[key]

    def __len__(self) -> int:
        return 6

    def __array__(self, dtype=None, copy=None):
        arr = self._as_array()
        return arr if dtype is None else arr.astype(dtype)

    @property
    def wrench(self) -> np.ndarray:
        return self._as_array()


@dataclass
class RotorInputs:
    """External driving conditions imposed on the rotor each timestep.

    These are quantities the vehicle or environment sets — nothing here
    has a derivative owned by the aero model. Rotor mechanical and inflow
    states (omega, spin_angle, lambda) live in RotorState instead.

    Fields
    ------
    collective_rad   collective pitch angle, rad
    tilt_lon         longitudinal cyclic tilt, rad
    tilt_lat         lateral cyclic tilt, rad
    R_hub            [3x3] rotation matrix, hub frame to world frame
    v_hub_world      [3] hub velocity in world frame, m/s
    wind_world       [3] wind velocity in world frame, m/s
    t                simulation time, s
    rho_kg_m3        air density, kg/m³ (default ISA sea level)
    motor_torque_Nm  shaft torque applied by motor/generator, N·m
                     positive = driving rotor, negative = braking
                     zero (default) = pure autorotation
    """

    collective_rad:  float
    tilt_lon:        float
    tilt_lat:        float
    R_hub:           np.ndarray
    v_hub_world:     np.ndarray
    wind_world:      np.ndarray
    t:               float
    rho_kg_m3:       float = field(default=1.225)
    motor_torque_Nm: float = field(default=0.0)


class AeroBase(ABC):
    """Abstract interface for aero models."""

    def initial_rotor_state(self) -> RotorState:
        """Return the zero rotor state for this model.

        Override in subclasses that use dynamic inflow (e.g. return
        PittPetersRotorState() for a Pitt-Peters model).  The integrator
        calls this once at initialisation to allocate the right state type.
        """
        return QuasiStaticRotorState()

    @abstractmethod
    def compute_forces(
        self,
        inputs: RotorInputs,
        state:  RotorState,
    ) -> "tuple[AeroResult, RotorState]":
        """Compute aerodynamic forces and rotor state derivatives.

        Parameters
        ----------
        inputs  External driving conditions for this timestep.
        state   Current rotor state (inflow + omega + spin_angle).

        Returns
        -------
        result      Forces and moments.
        derivative  dstate/dt as a RotorState — the caller integrates.
        """
        ...

    @abstractmethod
    def to_dict(self) -> dict: ...

    @classmethod
    def from_definition(cls, defn: "RotorDefinition") -> "AeroBase":
        return cls(defn)


from . import rotor_definition  # noqa: F401, E402
from .bem import BEMModel, prandtl_tip_loss, solve_bem_element  # noqa: F401, E402
from .pitt_peters import PittPetersModel, vrs_lambda1, prescribed_element_forces  # noqa: F401, E402
from .rotor_definition import (  # noqa: F401, E402
    AirfoilProperties,
    AutorotationProperties,
    BladeGeometry,
    ControlProperties,
    InertiaProperties,
    KamanFlap,
    RotorDefinition,
    ValidationIssue,
)

__all__ = [
    "AeroResult",
    "AeroBase",
    "BEMModel",
    "PittPetersModel",
    "vrs_lambda1",
    "prescribed_element_forces",
    "RotorInputs",
    "RotorState",
    "PittPetersRotorState",
    "QuasiStaticRotorState",
    "AirfoilPolar",
    "LinearPolar",
    "create_aero",
    "rotor_definition",
    "AirfoilProperties",
    "AutorotationProperties",
    "BladeGeometry",
    "ControlProperties",
    "InertiaProperties",
    "RotorDefinition",
    "KamanFlap",
    "ValidationIssue",
]


def create_aero(defn=None, model: str = "custom") -> AeroBase:
    _ = (defn, model)
    raise NotImplementedError("No aero model implementation is registered in this package.")
