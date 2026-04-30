"""Level 2: Pitt-Peters 3-state dynamic inflow with VRS empirical correction.

Inflow distribution (azimuth-averaged for axial flight):
    λ(r, ψ) = λ_0 + (r/R)(λ_c cosψ + λ_s sinψ)

Inflow ODE — first-order lags from Peters & HaQuang (1988):
    dλ_0/dt = (λ_0_ss − λ_0) / τ_0     τ_0  = 8R / (3π V_T)
    dλ_c/dt = (λ_c_ss − λ_c) / τ_cs    τ_cs = 16R / (45π V_T)
    dλ_s/dt = (λ_s_ss − λ_s) / τ_cs

Steady-state target for λ_0:
    hover / climb / WBS  →  momentum theory: λ_0_ss = T / (2ρA V_T ΩR)
    VRS  (0 < λ₂ < 2)   →  Leishman (2000) empirical polynomial

Cyclic targets λ_c_ss = λ_s_ss = 0 for symmetric axial flight.
Forward-flight azimuth integration is deferred to Level 3.

References
----------
Pitt & Peters (1981), Vertica 5(1), 21-34.
Peters & HaQuang (1988), JAHS 33(4), 64-68.
Leishman (2000), §12.4, §12.7.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from . import AeroBase, AeroResult, RotorInputs
from .bem import prandtl_tip_loss
from .rotor_definition import RotorDefinition
from .rotor_state import PittPetersRotorState, RotorState
from .polar import AirfoilPolar, LinearPolar, TabulatedPolar


# ---------------------------------------------------------------------------
# VRS empirical polynomial  (Leishman 2000, §12.7)
# ---------------------------------------------------------------------------
# λ_1/V_h = 1 + C[0]·λ₂ + C[1]·λ₂² + C[2]·λ₂³ + C[3]·λ₂⁴
# where λ₂ = V_descent / V_h > 0.  Valid for 0 ≤ λ₂ ≤ 2.
# Fit to Castles-Gray (NACA TN-2474) and Coleman (1945) measured data.
_VRS_C = (1.125, -1.372, 1.718, -0.655)


def vrs_lambda1(lambda2: float) -> float:
    """Normalized induced velocity λ₁ = v_i/V_h from Leishman VRS polynomial.

    lambda2  V_descent / V_h, must be in [0, 2]
    Returns  v_i / V_h  (= 1.0 at λ₂=0 hover; ≈1.0 at λ₂=2 WBS boundary)
    """
    k = lambda2
    return 1.0 + _VRS_C[0]*k + _VRS_C[1]*k**2 + _VRS_C[2]*k**3 + _VRS_C[3]*k**4


# ---------------------------------------------------------------------------
# Prescribed-inflow blade element
# ---------------------------------------------------------------------------

def prescribed_element_forces(
    r: float,
    dr: float,
    chord: float,
    twist_rad: float,
    collective_rad: float,
    omega: float,
    lambda_r: float,
    rho: float,
    n_blades: int,
    radius_m: float,
    polar: AirfoilPolar,
    use_tip_loss: bool,
) -> tuple[float, float]:
    """Blade element with prescribed inflow ratio lambda_r = v_a / (ΩR).

    Tangential induction (a′) is neglected — valid when thrust loading is
    moderate (a′ ≪ 1 in hover).

    Returns (dT [N], dQ [N·m]).
    """
    Omega_R = omega * radius_m
    if Omega_R < 1e-6:
        return 0.0, 0.0

    v_a = lambda_r * Omega_R
    v_t = omega * r
    if v_t < 1e-9:
        return 0.0, 0.0

    phi = math.atan2(v_a, v_t)
    alpha = (collective_rad + twist_rad) - phi
    cl, cd = polar.cl_cd(alpha)

    x = r / radius_m
    F = prandtl_tip_loss(n_blades, x, phi) if use_tip_loss else 1.0
    F = max(F, 1e-4)

    cn = cl * math.cos(phi) - cd * math.sin(phi)
    ct = cl * math.sin(phi) + cd * math.cos(phi)

    q = 0.5 * rho * (v_a**2 + v_t**2) * chord * dr * n_blades
    return q * cn, q * ct * r


# ---------------------------------------------------------------------------
# Pitt-Peters model
# ---------------------------------------------------------------------------

@dataclass
class PittPetersModel(AeroBase):
    """Level 2 BEM: Pitt-Peters 3-state dynamic inflow + VRS correction (NED).

    Accepts PittPetersRotorState — [λ_0, λ_c, λ_s, ω, ψ].
    Inflow states are dimensionless (v / ΩR).

    Initialization note: PittPetersRotorState() starts with λ_0 = 0 (no
    induced flow).  The ODE will converge to hover equilibrium within
    ~5 time constants τ_0 = 8R/(3π V_h).  For a 0.914 m rotor at 1200 rpm
    that is roughly 0.08 s — fast enough to ignore in most simulations.
    """

    defn: RotorDefinition
    _polar: AirfoilPolar = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if self.defn.airfoil.polar_csv is not None:
            self._polar = TabulatedPolar.from_csv(self.defn.airfoil.polar_csv)
        else:
            self._polar = LinearPolar.from_properties(self.defn.airfoil)

    def initial_rotor_state(self) -> PittPetersRotorState:
        return PittPetersRotorState()

    def compute_forces(
        self,
        inputs: RotorInputs,
        state: RotorState,
    ) -> tuple[AeroResult, RotorState]:
        assert isinstance(state, PittPetersRotorState)

        blade   = self.defn.blade
        omega   = state.omega_rad_s
        rho     = inputs.rho_kg_m3
        R       = blade.radius_m
        A       = math.pi * R**2
        Omega_R = omega * R

        hub_axis = inputs.R_hub @ np.array([0.0, 0.0, 1.0])
        v_rel    = inputs.wind_world - inputs.v_hub_world
        v_climb  = float(np.dot(v_rel, hub_axis))  # NED: < 0 for descent

        lam0  = state.lambda_0
        lam_c = state.lambda_c
        lam_s = state.lambda_s

        # ------------------------------------------------------------------
        # Blade element loop — prescribed uniform inflow lam0
        # (λ_c, λ_s contributions average to zero over azimuth)
        # ------------------------------------------------------------------
        n    = blade.n_elements
        r0   = blade.root_cutout_m
        dr   = (R - r0) / n
        r_mid = np.linspace(r0 + 0.5 * dr, R - 0.5 * dr, n)
        twist = math.radians(blade.twist_deg)

        # Total inflow = induced (lam0) + freestream (v_climb/ΩR).
        # Blade elements must see the combined axial velocity so that WBS
        # produces net-upward flow (lambda_total < 0) and autorotating torque.
        lambda_climb = v_climb / Omega_R if Omega_R > 1e-6 else 0.0
        lambda_total = lam0 + lambda_climb

        T_total = Q_total = 0.0
        if Omega_R > 1e-6:
            for rv in r_mid:
                dT, dQ = prescribed_element_forces(
                    r=float(rv), dr=dr, chord=blade.chord_m,
                    twist_rad=twist, collective_rad=inputs.collective_rad,
                    omega=omega, lambda_r=lambda_total, rho=rho,
                    n_blades=blade.n_blades, radius_m=R,
                    polar=self._polar,
                    use_tip_loss=self.defn.airfoil.tip_loss,
                )
                T_total += dT
                Q_total += dQ

        # ------------------------------------------------------------------
        # Pitt-Peters inflow ODE
        # ------------------------------------------------------------------
        # V_h: hover-equivalent inflow speed from current thrust
        T_pos = max(T_total, 0.0)
        V_h   = math.sqrt(T_pos / (2.0 * rho * A)) if T_pos > 1e-6 else 0.0

        # v0: current induced inflow in m/s; V_T: total axial flow through disk
        v0  = lam0 * Omega_R
        V_T = max(abs(v_climb + v0), 1e-2 * max(Omega_R, 1.0))

        # λ₂ = V_descent / V_h  (> 0 during descent)
        V_c   = max(-v_climb, 0.0)
        lam2  = (V_c / V_h) if V_h > 1e-3 else 0.0

        # Steady-state uniform inflow target
        if v_climb < -1e-3 and 0.0 < lam2 < 2.0:
            # VRS: Leishman polynomial replaces momentum theory
            lam0_ss = (vrs_lambda1(lam2) * V_h / Omega_R) if Omega_R > 1e-6 else 0.0
        else:
            # Hover / climb / WBS: momentum theory T = 2ρA v_i V_T
            lam0_ss = T_total / (2.0 * rho * A * V_T * Omega_R) if Omega_R > 1e-6 else 0.0

        # Apparent-mass time constants (Peters & HaQuang 1988)
        tau_0  = (8.0 * R) / (3.0 * math.pi * V_T)
        tau_cs = (16.0 * R) / (45.0 * math.pi * V_T)

        d_lam0  = (lam0_ss - lam0)  / tau_0
        d_lam_c = (0.0    - lam_c)  / tau_cs   # cyclic decays to zero (axial flight)
        d_lam_s = (0.0    - lam_s)  / tau_cs

        # ------------------------------------------------------------------
        # Mechanical states
        # ------------------------------------------------------------------
        I_ode = (
            self.defn.autorotation.I_ode_kgm2
            if self.defn.autorotation.I_ode_kgm2 is not None
            else 1.0
        )
        d_omega      = (-Q_total + inputs.motor_torque_Nm) / I_ode
        d_spin_angle = omega

        # ------------------------------------------------------------------
        # Assemble outputs
        # ------------------------------------------------------------------
        F_world = -T_total * hub_axis
        result  = AeroResult(
            F_world=F_world,
            M_orbital=np.zeros(3),
            Q_spin=Q_total,
            M_spin=Q_total * hub_axis,
        )
        derivative = PittPetersRotorState(
            lambda_0=d_lam0,
            lambda_c=d_lam_c,
            lambda_s=d_lam_s,
            omega_rad_s=d_omega,
            spin_angle_rad=d_spin_angle,
        )
        return result, derivative

    def to_dict(self) -> dict:
        return {"model": "PittPeters_Level2", "n_elements": self.defn.blade.n_elements}

    @classmethod
    def from_definition(cls, defn: RotorDefinition) -> "PittPetersModel":
        return cls(defn=defn)
