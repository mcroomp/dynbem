"""Aerodynamic interfaces and factory hooks."""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
import numpy as np

from .rotor_state import PittPetersRotorState, QuasiStaticRotorState, RotorState
from .polar import AirfoilPolar, LinearPolar


@dataclass
class AeroResult:
    """Return type for compute_forces().

    Fields
    ------
    F_world   [3] Net aerodynamic force on the rotor in world (NED) frame, N.
              Equal to ``-T_total * hub_axis_ned``; for a level rotor with
              positive thrust, ``F_world[2] < 0`` (upward = −Z).
    M_orbital [3] In-plane hub moments from non-axisymmetric thrust, world
              frame. Accumulated per-element as ``r · dT · [sin ψ, cos ψ, 0]``
              in hub frame and rotated via ``R_hub``. Non-zero in forward
              flight (advancing/retreating velocity asymmetry) and under
              cyclic input. Zero in axisymmetric hover. See CLAUDE.md.
    Q_spin    Shaft drag torque magnitude, N·m. Positive in powered hover
              (aero drag opposes rotor motion). Drives the rotor speed ODE
              via ``d_omega = (-Q_spin + motor_torque) / I_ode``.
    M_spin    [3] Reaction torque on the airframe from the rotor system
              (shaft + motor stator), world frame. For our CCW-from-above
              convention this is ``+Q_spin * hub_axis_ned`` — airframe is
              pushed to spin CW from above (American helicopter yaw-right
              tendency without tail rotor).
    """
    F_world:   np.ndarray   # [3]
    M_orbital: np.ndarray   # [3]
    Q_spin:    float
    M_spin:    np.ndarray   # [3]


@dataclass
class RotorInputs:
    """External driving conditions imposed on the rotor each timestep.

    These are quantities the vehicle or environment sets — nothing here
    has a derivative owned by the aero model. Rotor mechanical and inflow
    states (omega, spin_angle, lambda) live in RotorState instead.

    Fields
    ------
    collective_rad   collective pitch angle, rad
    tilt_lon         longitudinal swashplate tilt, rad. Helicopter-standard
                     sign: positive → nose-down (forward stick). Mapped to
                     blade pitch via aero.cyclic.cyclic_coeffs() using the
                     rotor's ControlProperties (gain, phase). With default
                     gain=1, phase=0 (control=None) this acts as the
                     direct θ_1c blade-pitch amplitude with helicopter signs.
    tilt_lat         lateral swashplate tilt, rad. Helicopter-standard sign:
                     positive → roll right.
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

    def inflow_taus(
        self,
        inputs: RotorInputs,
        state:  RotorState,
    ) -> np.ndarray:
        """Time constants per state component for semi-implicit integration.

        Returns an array of the same length as ``state.to_array()``.  Each
        element is the time constant τ of that state's first-order lag
        (used by the envelope integrator's semi-implicit damping
        ``damp = 1/(1 + dt/τ)``).  Use ``np.inf`` for states that should
        be integrated as plain explicit Euler (mechanical ω, ψ, and any
        quasi-static states with no dynamics).

        Default implementation: all-infinity (plain explicit Euler for
        every state).  Models with stiff dynamic inflow override this so
        the envelope can damp the explicit step.
        """
        return np.full(state.to_array().shape, np.inf)



from . import rotor_definition  # noqa: F401, E402
from .bem import (  # noqa: F401, E402
    BEMModel,
    prandtl_hub_loss,
    prandtl_tip_loss,
    solve_bem_element,
)
from .pitt_peters import PittPetersModel, vrs_lambda1  # noqa: F401, E402
from .oye import OyeBEMModel  # noqa: F401, E402
from .rotor_state import OyeRotorState  # noqa: F401, E402
from .trim import TrimResult, relax_inflow, solve_trim_cyclic  # noqa: F401, E402
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
    "OyeBEMModel",
    "vrs_lambda1",
    "RotorInputs",
    "RotorState",
    "PittPetersRotorState",
    "OyeRotorState",
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
    "TrimResult",
    "solve_trim_cyclic",
    "relax_inflow",
]


def create_aero(defn: "RotorDefinition", model: str = "pitt_peters_jit") -> AeroBase:
    """Factory for the aero models in this package.

    model
      "bem"               BEMModel (Level 1, quasi-static inflow)
      "pitt_peters"       PittPetersModel (Level 2, Pitt-Peters L-matrix, numpy)
      "pitt_peters_jit"   PittPetersModelJIT (Level 2, JIT-compiled — default)
      "oye"               OyeBEMModel (Level 2, Øye 2-stage annular inflow,
                          JIT-compiled).  Stable alternative to Pitt-Peters
                          at high advance ratios / descent + edgewise wind.
    """
    if model == "bem":
        from .bem import BEMModel
        return BEMModel(defn=defn)
    if model in ("pitt_peters", "pitt_peters_numpy"):
        from .pitt_peters import PittPetersModel
        return PittPetersModel(defn=defn)
    if model in ("pitt_peters_jit", "jit"):
        from .pitt_peters_jit import PittPetersModelJIT
        return PittPetersModelJIT(defn=defn)
    if model in ("oye", "oye_bem"):
        from .oye import OyeBEMModel
        return OyeBEMModel(defn=defn)
    raise ValueError(
        f"Unknown aero model {model!r}. "
        f"Choose 'bem', 'pitt_peters', 'pitt_peters_jit', or 'oye'."
    )
