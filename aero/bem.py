"""Level 1 BEM solver — multi-element, NED frame.

Coordinate system: NED (North-East-Down).
  - Rotor hub axis points in the +Z (down) direction for a level hover.
  - Thrust opposes gravity: F_world[2] < 0 (upward force = negative Z).
  - v_climb > 0: axial freestream flows downward through disk (helicopter climb).
  - v_climb < 0: axial freestream flows upward through disk (autorotation/turbine).
  - v_climb = 0: hover — rotor generates its own induced inflow.

Inflow iteration uses the total inflow ratio λ_r = v_a/(Ω·R) rather than the
wind-turbine induction-factor form (which degenerates at v_climb = 0 / hover).

Momentum-BEM derivation
-----------------------
Momentum annulus (with Prandtl tip-loss F):
    dCT/dx = 4·F·x·λ_r·(λ_r − λ_c)

Blade element (N blades, chord c, local solidity σ_r = N·c/(2π·r)):
    dCT/dx = σ_r·x·cn·(λ_r² + x²)

Cancel x, define k = σ_r·cn/(4·F):
    k·(λ_r² + x²) = λ_r·(λ_r − λ_c)   ← the iteration equation
    (k−1)·λ_r² + λ_c·λ_r + k·x² = 0   ← quadratic form

Integration contract
--------------------
compute_forces() returns (AeroResult, derivative) where derivative holds dstate/dt.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import NamedTuple

import numpy as np

from . import AeroBase, AeroResult, RotorInputs
from .rotor_definition import RotorDefinition
from .rotor_state import QuasiStaticRotorState, RotorState
from .polar import AirfoilPolar, LinearPolar, TabulatedPolar


# ---------------------------------------------------------------------------
# Prandtl tip-loss (public — used in component tests)
# ---------------------------------------------------------------------------

def prandtl_tip_loss(n_blades: int, x: float, phi_rad: float) -> float:
    """Prandtl tip-loss factor F at normalised radius x = r/R.

    Returns F ∈ (0, 1].  F → 1 far from tip; F → 0 at the tip for small phi.

    n_blades  number of blades
    x         r/R, normalised radius ∈ (0, 1]
    phi_rad   flow angle from rotor plane (rad); may be negative in turbine mode
    """
    if abs(phi_rad) < 1e-9 or x >= 1.0:
        return 1.0
    f = (n_blades / 2.0) * (1.0 - x) / (x * abs(math.sin(phi_rad)))
    return (2.0 / math.pi) * math.acos(min(1.0, math.exp(-f)))


# ---------------------------------------------------------------------------
# BEM element (public — exposed for unit testing)
# ---------------------------------------------------------------------------

_MAX_BEM_ITER = 60
_BEM_TOL = 1e-7


class BEMElementResult(NamedTuple):
    """Converged state and forces for one blade-element annulus."""
    lambda_r: float   # total axial inflow ratio v_a/(Ω·R)
    a_prime: float    # tangential induction factor
    dT: float         # thrust contribution (N) — positive = upward for level rotor
    dQ: float         # reaction torque contribution (N·m) — positive = drag on rotor
    momentum_residual: float  # |4F·λ_r·(λ_r−λ_c) − σ_r·cn·(λ_r²+x²)| at convergence


def solve_bem_element(
    r: float,
    dr: float,
    chord: float,
    twist_rad: float,
    collective_rad: float,
    omega: float,
    v_climb: float,
    rho: float,
    n_blades: int,
    radius_m: float,
    polar: AirfoilPolar,
    use_tip_loss: bool,
) -> BEMElementResult:
    """Solve one blade-element annulus; return converged state and forces.

    v_climb  External axial freestream in the hub-axis direction (NED):
               > 0  air from above (helicopter climb)
               = 0  hover
               < 0  air from below (autorotation / flying wind turbine)

    dT > 0: thrust opposing inflow (upward for level rotor = −Z in NED).
    dQ > 0: reaction torque opposing rotor spin.
    """
    x = r / radius_m
    Omega_R = omega * radius_m
    if Omega_R < 1e-6:
        return BEMElementResult(0.0, 0.0, 0.0, 0.0, 0.0)

    sigma_r = n_blades * chord / (2.0 * math.pi * r)
    theta = collective_rad + twist_rad
    lambda_c = v_climb / Omega_R

    if lambda_c >= 0.0:
        lambda_r = max(lambda_c + 0.03, 0.02)
    else:
        lambda_r = min(lambda_c * 0.85, -0.02)

    a_prime = 0.0
    cn_final = 0.0
    F_final = 1.0

    for _ in range(_MAX_BEM_ITER):
        v_a = lambda_r * Omega_R
        v_t = omega * r * (1.0 + a_prime)
        if v_t < 1e-9:
            break

        phi = math.atan2(v_a, v_t)
        alpha = theta - phi
        cl, cd = polar.cl_cd(alpha)

        F = prandtl_tip_loss(n_blades, x, phi) if use_tip_loss else 1.0
        F = max(F, 1e-4)
        F_final = F

        cn = cl * math.cos(phi) - cd * math.sin(phi)
        ct = cl * math.sin(phi) + cd * math.cos(phi)
        cn_final = cn

        # Momentum-BEM quadratic: k = σ_r·cn/(4·F)
        k = sigma_r * cn / (4.0 * F)

        if abs(k - 1.0) > 1e-6:
            disc = max(0.0, lambda_c**2 - 4.0 * (k - 1.0) * k * x**2)
            sq = math.sqrt(disc)
            denom = 2.0 * (k - 1.0)
            r1 = (-lambda_c + sq) / denom
            r2 = (-lambda_c - sq) / denom
            if lambda_c >= 0.0:
                lambda_r_new = r2 if r2 > 0.0 else r1
            else:
                lambda_r_new = r1 if r1 < 0.0 else r2
        else:
            if abs(lambda_c) > 1e-8:
                lambda_r_new = -k * x**2 / lambda_c
            else:
                lambda_r_new = x * math.sqrt(max(0.0, k))

        lambda_r_new = max(-2.0, min(2.0, lambda_r_new))

        sc = math.sin(phi) * math.cos(phi)
        if abs(sc) > 1e-8 and abs(ct) > 1e-10:
            ap_denom = 4.0 * F * sc / (sigma_r * ct) - 1.0
            a_prime_new = (1.0 / ap_denom) if abs(ap_denom) > 1e-8 else 0.0
            a_prime_new = max(-0.5, min(0.5, a_prime_new))
        else:
            a_prime_new = 0.0

        converged = (
            abs(lambda_r_new - lambda_r) < _BEM_TOL
            and abs(a_prime_new - a_prime) < _BEM_TOL
        )
        lambda_r = 0.5 * lambda_r + 0.5 * lambda_r_new
        a_prime = 0.5 * a_prime + 0.5 * a_prime_new
        if converged:
            break

    # Recompute final forces at converged state
    v_a_f = lambda_r * Omega_R
    v_t_f = omega * r * (1.0 + a_prime)
    v_rel = math.sqrt(v_a_f**2 + v_t_f**2)
    phi_f = math.atan2(v_a_f, v_t_f)
    alpha_f = theta - phi_f
    cl_f, cd_f = polar.cl_cd(alpha_f)
    cn_f = cl_f * math.cos(phi_f) - cd_f * math.sin(phi_f)
    ct_f = cl_f * math.sin(phi_f) + cd_f * math.cos(phi_f)

    q_dyn = 0.5 * rho * v_rel**2 * chord * dr * n_blades
    dT = q_dyn * cn_f
    dQ = q_dyn * ct_f * r

    # Momentum-BEM balance residual — should be near zero at convergence.
    # 4·F·λ_r·(λ_r−λ_c) == σ_r·cn·(λ_r²+x²)
    momentum_side = 4.0 * F_final * lambda_r * (lambda_r - lambda_c)
    blade_side = sigma_r * cn_f * (lambda_r**2 + x**2)
    residual = abs(momentum_side - blade_side)

    return BEMElementResult(lambda_r, a_prime, dT, dQ, residual)


# ---------------------------------------------------------------------------
# BEM model
# ---------------------------------------------------------------------------

@dataclass
class BEMModel(AeroBase):
    """Multi-element BEM rotor model (Level 1).

    NED frame throughout.  Returns QuasiStaticRotorState derivatives.
    """

    defn: RotorDefinition
    _polar: AirfoilPolar = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if self.defn.airfoil.polar_csv is not None:
            self._polar = TabulatedPolar.from_csv(self.defn.airfoil.polar_csv)
        else:
            self._polar = LinearPolar.from_properties(self.defn.airfoil)

    def initial_rotor_state(self) -> QuasiStaticRotorState:
        return QuasiStaticRotorState()

    def compute_forces(
        self,
        inputs: RotorInputs,
        state: RotorState,
    ) -> tuple[AeroResult, RotorState]:
        assert isinstance(state, QuasiStaticRotorState)

        blade = self.defn.blade
        airfoil = self.defn.airfoil
        omega = state.omega_rad_s

        hub_axis_ned = inputs.R_hub @ np.array([0.0, 0.0, 1.0])
        v_rel_world = inputs.wind_world - inputs.v_hub_world
        v_climb = float(np.dot(v_rel_world, hub_axis_ned))

        n = blade.n_elements
        r_root, r_tip = blade.root_cutout_m, blade.radius_m
        dr = (r_tip - r_root) / n
        r_stations = np.linspace(r_root + 0.5 * dr, r_tip - 0.5 * dr, n)
        twist_rad = math.radians(blade.twist_deg)

        T_total = 0.0
        Q_total = 0.0
        for r in r_stations:
            elem = solve_bem_element(
                r=float(r), dr=dr, chord=blade.chord_m,
                twist_rad=twist_rad, collective_rad=inputs.collective_rad,
                omega=omega, v_climb=v_climb, rho=inputs.rho_kg_m3,
                n_blades=blade.n_blades, radius_m=blade.radius_m,
                polar=self._polar, use_tip_loss=airfoil.tip_loss,
            )
            T_total += elem.dT
            Q_total += elem.dQ

        F_world = -T_total * hub_axis_ned
        M_orbital = np.zeros(3)
        M_spin = Q_total * hub_axis_ned

        I_ode = (
            self.defn.autorotation.I_ode_kgm2
            if self.defn.autorotation.I_ode_kgm2 is not None
            else 1.0
        )
        d_omega = (-Q_total + inputs.motor_torque_Nm) / I_ode
        d_spin_angle = omega

        result = AeroResult(
            F_world=F_world,
            M_orbital=M_orbital,
            Q_spin=Q_total,
            M_spin=M_spin,
        )
        derivative = QuasiStaticRotorState(
            omega_rad_s=d_omega,
            spin_angle_rad=d_spin_angle,
        )
        return result, derivative

    def to_dict(self) -> dict:
        return {"model": "BEM_Level1", "n_elements": self.defn.blade.n_elements}

    @classmethod
    def from_definition(cls, defn: RotorDefinition) -> "BEMModel":
        return cls(defn=defn)
