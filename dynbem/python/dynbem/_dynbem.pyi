from __future__ import annotations

from typing import TypeAlias, TypedDict, Unpack

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


class LinearPolar:
    CL0: float
    CL_alpha_per_rad: float
    CD0: float
    alpha_stall_rad: float

    def __init__(
        self, CL0: float, CL_alpha_per_rad: float, CD0: float, alpha_stall_rad: float
    ) -> None: ...
    @staticmethod
    def from_properties(airfoil: LinearPolarParameters) -> LinearPolar: ...
    def cl_cd(self, alpha_rad: float) -> tuple[float, float]: ...
    def cl_cd_arr(self, alpha_rad: FloatArray) -> tuple[FloatArray, FloatArray]: ...


class TabulatedPolar:
    def __init__(self, alpha_rad: FloatArray, cl: FloatArray, cd: FloatArray) -> None: ...
    def cl_cd(self, alpha_rad: float) -> tuple[float, float]: ...
    def cl_cd_arr(self, alpha_rad: FloatArray) -> tuple[FloatArray, FloatArray]: ...


class BladeGeometry:
    n_blades: int
    radius_m: float
    root_cutout_m: float
    chord_m: float
    twist_deg: float
    n_elements: int
    tip_loss: bool
    span_m: float
    r_cp_m: float
    disk_area_m2: float
    solidity: float
    has_radial_stations: bool

    def __new__(
        cls: type[BladeGeometry],
        n_blades: int,
        radius_m: float,
        root_cutout_m: float,
        chord_m: float,
        twist_deg: float,
        n_elements: int,
        r_stations_m: list[float],
        chord_stations_m: list[float],
        twist_stations_deg: list[float],
        tip_loss: bool = True,
    ) -> BladeGeometry: ...
    def __init__(
        self,
        n_blades: int,
        radius_m: float,
        root_cutout_m: float,
        chord_m: float,
        twist_deg: float,
        n_elements: int,
        r_stations_m: list[float],
        chord_stations_m: list[float],
        twist_stations_deg: list[float],
        tip_loss: bool = True,
    ) -> None: ...
    def chord_at(self, r: float) -> float: ...
    def twist_at(self, r: float) -> float: ...


class LinearPolarParameters:
    CL0: float
    CL_alpha_per_rad: float
    CD0: float
    alpha_stall_deg: float

    def __init__(
        self, CL0: float, CL_alpha_per_rad: float, CD0: float, alpha_stall_deg: float
    ) -> None: ...


class ControlProperties:
    swashplate_pitch_gain_rad: float
    swashplate_phase_deg: float | None

    def __init__(
        self, swashplate_pitch_gain_rad: float, swashplate_phase_deg: float | None
    ) -> None: ...


class ServoFlapGeometry:
    C_M_delta_per_rad: float
    r_inner_m: float
    r_outer_m: float

    def __init__(
        self, C_M_delta_per_rad: float, r_inner_m: float, r_outer_m: float
    ) -> None: ...


class ServoFlapActuation:
    I_theta_kgm2: float
    damper_Nms_per_rad: float
    ac_offset_m: float
    blade_Cm_AC: float
    flap: ServoFlapGeometry

    def __init__(
        self,
        I_theta_kgm2: float,
        damper_Nms_per_rad: float,
        flap: ServoFlapGeometry,
        ac_offset_m: float = 0.0,
        blade_Cm_AC: float = 0.0,
    ) -> None: ...


class FlapProperties:
    I_blade_flap_kgm2: float
    omega_nr_rad_s: float

    def __init__(self, I_blade_flap_kgm2: float, omega_nr_rad_s: float) -> None: ...


class RotorDefinition:
    blade: BladeGeometry
    airfoil: LinearPolarParameters
    name: str
    description: str
    control: ControlProperties | None
    servoflap: ServoFlapActuation | None
    flap: FlapProperties | None

    def __init__(
        self,
        blade: BladeGeometry,
        airfoil: LinearPolarParameters,
        control: ControlProperties | None,
        name: str,
        description: str,
        servoflap: ServoFlapActuation | None = None,
        flap: FlapProperties | None = None,
    ) -> None: ...


class QuasiStaticRotorState:
    def __init__(self) -> None: ...
    def to_array(self) -> FloatArray: ...
    def from_array(self, arr: FloatArray) -> QuasiStaticRotorState: ...


class PittPetersRotorState:
    lambda_0: float
    lambda_c: float
    lambda_s: float

    def __init__(self, lambda_0: float, lambda_c: float, lambda_s: float) -> None: ...
    def to_array(self) -> FloatArray: ...
    def from_array(self, arr: FloatArray) -> PittPetersRotorState: ...


class OyeRotorState:
    W_int: FloatArray
    W: FloatArray

    def __init__(self, W_int: FloatArray, W: FloatArray) -> None: ...
    @staticmethod
    def zeros(n_elements: int) -> OyeRotorState: ...
    def to_array(self) -> FloatArray: ...
    def from_array(self, arr: FloatArray) -> OyeRotorState: ...


class VpmRotorState:
    n_particles: int

    def __init__(self) -> None: ...
    @staticmethod
    def zeros(n_elements: int = 0) -> VpmRotorState: ...
    def to_array(self) -> FloatArray: ...
    def from_array(self, arr: FloatArray) -> VpmRotorState: ...


class RotorInputs:
    collective_rad: float
    tilt_lon: float
    tilt_lat: float
    R_hub: NDArray[np.float64]
    v_hub_world: FloatArray
    wind_world: FloatArray
    omega_rad_s: float
    rho_kg_m3: float

    def __init__(
        self,
        collective_rad: float,
        tilt_lon: float,
        tilt_lat: float,
        R_hub: NDArray[np.float64],
        v_hub_world: FloatArray,
        wind_world: FloatArray,
        omega_rad_s: float,
        rho_kg_m3: float,
    ) -> None: ...


class AeroResult:
    F_world: FloatArray
    m_hub_world: FloatArray
    Q_spin: float
    M_spin: FloatArray


class VpmConfig(TypedDict, total=False):
    max_particles: int
    sigma: float
    relax: float
    nonlinear_lifting_line: bool
    tip_clustering: bool
    local_core: bool
    barnes_hut: bool
    bh_theta: float
    bh_min_particles: int
    merge_wake: bool
    merge_every: int
    merge_min_particles: int
    merge_kappa: float
    merge_chi_min: float
    merge_region_dist: float
    core_spread_nu: float
    strength_decay_tau_rev: float


class _QuasiStaticBEMLinear:
    def __init__(
        self, defn: RotorDefinition, polar: LinearPolar, n_psi_elements: int
    ) -> None: ...
    @property
    def defn(self) -> RotorDefinition: ...
    @property
    def n_psi_elements(self) -> int: ...
    def initial_rotor_state(self) -> QuasiStaticRotorState: ...
    def compute_forces(
        self, inputs: RotorInputs, state: QuasiStaticRotorState
    ) -> tuple[AeroResult, QuasiStaticRotorState]: ...
    def step(
        self,
        inputs: RotorInputs,
        state: QuasiStaticRotorState,
        dt: float,
        integration_method: str | None = None,
    ) -> tuple[AeroResult, QuasiStaticRotorState]: ...
    def inflow_taus(
        self, inputs: RotorInputs, state: QuasiStaticRotorState
    ) -> FloatArray: ...


class _QuasiStaticBEMTabulated:
    def __init__(
        self, defn: RotorDefinition, polar: TabulatedPolar, n_psi_elements: int
    ) -> None: ...
    @property
    def defn(self) -> RotorDefinition: ...
    @property
    def n_psi_elements(self) -> int: ...
    def initial_rotor_state(self) -> QuasiStaticRotorState: ...
    def compute_forces(
        self, inputs: RotorInputs, state: QuasiStaticRotorState
    ) -> tuple[AeroResult, QuasiStaticRotorState]: ...
    def step(
        self,
        inputs: RotorInputs,
        state: QuasiStaticRotorState,
        dt: float,
        integration_method: str | None = None,
    ) -> tuple[AeroResult, QuasiStaticRotorState]: ...
    def inflow_taus(
        self, inputs: RotorInputs, state: QuasiStaticRotorState
    ) -> FloatArray: ...


class _PittPetersModelLinear:
    def __init__(
        self, defn: RotorDefinition, polar: LinearPolar, n_psi_elements: int
    ) -> None: ...
    @property
    def defn(self) -> RotorDefinition: ...
    @property
    def n_psi_elements(self) -> int: ...
    def initial_rotor_state(self) -> PittPetersRotorState: ...
    def compute_forces(
        self, inputs: RotorInputs, state: PittPetersRotorState
    ) -> tuple[AeroResult, PittPetersRotorState]: ...
    def step(
        self,
        inputs: RotorInputs,
        state: PittPetersRotorState,
        dt: float,
        integration_method: str | None = None,
    ) -> tuple[AeroResult, PittPetersRotorState]: ...
    def inflow_taus(
        self, inputs: RotorInputs, state: PittPetersRotorState
    ) -> FloatArray: ...


class _PittPetersModelTabulated:
    def __init__(
        self, defn: RotorDefinition, polar: TabulatedPolar, n_psi_elements: int
    ) -> None: ...
    @property
    def defn(self) -> RotorDefinition: ...
    @property
    def n_psi_elements(self) -> int: ...
    def initial_rotor_state(self) -> PittPetersRotorState: ...
    def compute_forces(
        self, inputs: RotorInputs, state: PittPetersRotorState
    ) -> tuple[AeroResult, PittPetersRotorState]: ...
    def step(
        self,
        inputs: RotorInputs,
        state: PittPetersRotorState,
        dt: float,
        integration_method: str | None = None,
    ) -> tuple[AeroResult, PittPetersRotorState]: ...
    def inflow_taus(
        self, inputs: RotorInputs, state: PittPetersRotorState
    ) -> FloatArray: ...


class _OyeBEMModelLinear:
    def __init__(
        self, defn: RotorDefinition, polar: LinearPolar, n_psi_elements: int, coupling_k: float
    ) -> None: ...
    @property
    def defn(self) -> RotorDefinition: ...
    @property
    def n_psi_elements(self) -> int: ...
    @property
    def coupling_k(self) -> float: ...
    def initial_rotor_state(self) -> OyeRotorState: ...
    def compute_forces(
        self, inputs: RotorInputs, state: OyeRotorState
    ) -> tuple[AeroResult, OyeRotorState]: ...
    def step(
        self,
        inputs: RotorInputs,
        state: OyeRotorState,
        dt: float,
        integration_method: str | None = None,
    ) -> tuple[AeroResult, OyeRotorState]: ...
    def inflow_taus(self, inputs: RotorInputs, state: OyeRotorState) -> FloatArray: ...


class _OyeBEMModelTabulated:
    def __init__(
        self,
        defn: RotorDefinition,
        polar: TabulatedPolar,
        n_psi_elements: int,
        coupling_k: float,
    ) -> None: ...
    @property
    def defn(self) -> RotorDefinition: ...
    @property
    def n_psi_elements(self) -> int: ...
    @property
    def coupling_k(self) -> float: ...
    def initial_rotor_state(self) -> OyeRotorState: ...
    def compute_forces(
        self, inputs: RotorInputs, state: OyeRotorState
    ) -> tuple[AeroResult, OyeRotorState]: ...
    def step(
        self,
        inputs: RotorInputs,
        state: OyeRotorState,
        dt: float,
        integration_method: str | None = None,
    ) -> tuple[AeroResult, OyeRotorState]: ...
    def inflow_taus(self, inputs: RotorInputs, state: OyeRotorState) -> FloatArray: ...


class _VpmRotorLinear:
    def __init__(
        self,
        defn: RotorDefinition,
        polar: LinearPolar,
        **config: Unpack[VpmConfig],
    ) -> None: ...
    def initial_rotor_state(self) -> VpmRotorState: ...
    def compute_forces(
        self, inputs: RotorInputs, state: VpmRotorState
    ) -> tuple[AeroResult, VpmRotorState]: ...
    def step(
        self,
        inputs: RotorInputs,
        state: VpmRotorState,
        dt: float,
        integration_method: str | None = None,
    ) -> tuple[AeroResult, VpmRotorState]: ...


class _VpmRotorTabulated:
    def __init__(
        self,
        defn: RotorDefinition,
        polar: TabulatedPolar,
        **config: Unpack[VpmConfig],
    ) -> None: ...
    def initial_rotor_state(self) -> VpmRotorState: ...
    def compute_forces(
        self, inputs: RotorInputs, state: VpmRotorState
    ) -> tuple[AeroResult, VpmRotorState]: ...
    def step(
        self,
        inputs: RotorInputs,
        state: VpmRotorState,
        dt: float,
        integration_method: str | None = None,
    ) -> tuple[AeroResult, VpmRotorState]: ...


class TrimResult:
    tilt_lon: float
    tilt_lat: float
    Mx_residual: float
    My_residual: float
    iterations: int
    converged: bool
    final_state: QuasiStaticRotorState | PittPetersRotorState | OyeRotorState


QuasiStaticModel: TypeAlias = _QuasiStaticBEMLinear | _QuasiStaticBEMTabulated
PittPetersModel: TypeAlias = _PittPetersModelLinear | _PittPetersModelTabulated
OyeBEMModel: TypeAlias = _OyeBEMModelLinear | _OyeBEMModelTabulated
VpmRotor: TypeAlias = _VpmRotorLinear | _VpmRotorTabulated
AeroModel: TypeAlias = QuasiStaticModel | PittPetersModel | OyeBEMModel | VpmRotor


def cyclic_coeffs(
    tilt_lon: float, tilt_lat: float, control: ControlProperties | None = None
) -> tuple[float, float]: ...

def omega_derivative(
    omega: float,
    q_aero: float,
    motor_torque_nm: float,
    i_ode_kgm2: float,
    bearing_friction_nm: float = 0.0,
) -> float: ...

def step_omega(
    omega: float,
    q_aero: float,
    motor_torque_nm: float,
    i_ode_kgm2: float,
    dt: float,
    bearing_friction_nm: float = 0.0,
) -> float: ...

def solve_trim_cyclic_py(
    aero: object,
    state: QuasiStaticRotorState | PittPetersRotorState | OyeRotorState,
    base_inputs: RotorInputs,
    *,
    target_moment: tuple[float, float] = (0.0, 0.0),
    tilt_lon_init: float = 0.0,
    tilt_lat_init: float = 0.0,
    tilt_min: float = -0.2617993877991494,
    tilt_max: float = 0.2617993877991494,
    tolerance_Nm: float = 0.02,
    max_iterations: int = 50,
    probe_rad: float = 0.008726646259971648,
    dt_relax: float = 0.005,
    n_inflow_relax: int = 100,
    n_settle: int = 0,
) -> TrimResult: ...

def relax_inflow_py(
    aero: object,
    state: QuasiStaticRotorState | PittPetersRotorState | OyeRotorState,
    inputs: RotorInputs,
    n_steps: int = 200,
    dt: float = 0.005,
) -> QuasiStaticRotorState | PittPetersRotorState | OyeRotorState: ...
