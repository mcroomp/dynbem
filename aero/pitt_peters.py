"""Level 2: Pitt-Peters 3-state dynamic inflow with VRS empirical correction.

Inflow distribution (uniform + linear tilt):
    λ(r, ψ) = λ_0 + (r/R)(λ_c cosψ + λ_s sinψ)

Inflow ODE — first-order lags from Peters & HaQuang (1988):
    dλ_0/dt = (λ_0_ss − λ_0) / τ_0     τ_0  = 8R / (3π V_T)
    dλ_c/dt = (λ_c_ss − λ_c) / τ_cs    τ_cs = 16R / (45π V_T)
    dλ_s/dt = (λ_s_ss − λ_s) / τ_cs

Steady-state target for λ_0:
    hover / climb / WBS  →  momentum theory: λ_0_ss = T / (2ρA V_T ΩR)
    VRS  (0 < λ₂ < 2)   →  Leishman (2000) empirical polynomial

Cyclic steady-state targets (Glauert skewed-wake model):
    axial flight  → λ_c_ss = λ_s_ss = 0
    forward flight (µ > 0.01) → λ_c_ss = −µ_x · tan(χ/2)
                                 λ_s_ss = −µ_y · tan(χ/2)
    where tan(χ/2) = µ / (√(µ² + λ_total²) + |λ_total|)
    and µ_x, µ_y are advance ratio components in the hub frame.

Forward-flight blade element loop averages over n_psi_elements azimuth
stations, using the full λ(r,ψ) inflow distribution at each station.
V_T uses total disk velocity √(v_edge² + (v_climb + v_i)²) in forward flight.

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
    v_t_extra: float = 0.0,
) -> tuple[float, float]:
    """Blade element with prescribed inflow ratio lambda_r = v_a / (ΩR).

    Tangential induction (a′) is neglected — valid when thrust loading is
    moderate (a′ ≪ 1 in hover).

    v_t_extra  In-plane wind tangential contribution at this azimuth (m/s).
               Positive adds to blade tangential velocity (advancing side).
               Zero (default) gives the standard axial-flight result.

    Returns (dT [N], dQ [N·m]).
    """
    Omega_R = omega * radius_m
    if Omega_R < 1e-6:
        return 0.0, 0.0

    v_a = lambda_r * Omega_R
    v_t = omega * r + v_t_extra
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

    In forward flight (µ > 0.01) the blade element loop integrates over
    n_psi_elements azimuth stations using the full λ(r,ψ) distribution,
    and cyclic steady-state targets are set by the Glauert skewed-wake model.

    Initialization note: PittPetersRotorState() starts with λ_0 = 0 (no
    induced flow).  The ODE will converge to hover equilibrium within
    ~5 time constants τ_0 = 8R/(3π V_h).  For a 0.914 m rotor at 1200 rpm
    that is roughly 0.08 s — fast enough to ignore in most simulations.
    """

    defn: RotorDefinition
    n_psi_elements: int = 36
    _polar: AirfoilPolar = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if self.defn.airfoil.polar_csv is not None:
            self._polar = TabulatedPolar.from_csv(self.defn.airfoil.polar_csv)
        else:
            self._polar = LinearPolar.from_properties(self.defn.airfoil)

        # Cache fixed radial geometry (n_elements, root_cutout, radius, etc.
        # don't change per call).  This shaves a linspace + per-step Python
        # overhead from compute_forces.
        blade = self.defn.blade
        R = blade.radius_m
        n = blade.n_elements
        r0 = blade.root_cutout_m
        self._dr = (R - r0) / n
        self._r_mid = np.linspace(r0 + 0.5 * self._dr, R - 0.5 * self._dr, n)
        self._x_mid = self._r_mid / R
        self._twist_rad = math.radians(blade.twist_deg)

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

        # In-plane (edgewise) wind component
        v_inplane = v_rel - v_climb * hub_axis
        v_edge    = float(np.linalg.norm(v_inplane))
        mu        = v_edge / max(Omega_R, 1e-6)

        # Advance ratio components in hub frame (for cyclic targets)
        v_inplane_hub = inputs.R_hub.T @ v_inplane
        mu_x = float(v_inplane_hub[0]) / max(Omega_R, 1e-6)
        mu_y = float(v_inplane_hub[1]) / max(Omega_R, 1e-6)

        lam0  = state.lambda_0
        lam_c = state.lambda_c
        lam_s = state.lambda_s

        # Total axial inflow (induced + freestream) seen by the blade disk
        lambda_climb = v_climb / Omega_R if Omega_R > 1e-6 else 0.0
        lambda_total = lam0 + lambda_climb

        # ------------------------------------------------------------------
        # Blade element forces (vectorized over r; loop over ψ in fwd flight)
        # ------------------------------------------------------------------
        r_arr = self._r_mid               # (n,)
        x_arr = self._x_mid               # (n,)
        dr    = self._dr
        twist = self._twist_rad
        chord = blade.chord_m
        n_b   = blade.n_blades
        col   = inputs.collective_rad

        T_total = Q_total = 0.0

        if Omega_R > 1e-6:
            if mu > 0.01 and omega > 1.0:
                # Forward flight: ψ-loop, vectorized over r at each station.
                n_psi = self.n_psi_elements
                T_acc = 0.0
                Q_acc = 0.0
                v_t_base = omega * r_arr     # (n,)
                for i_psi in range(n_psi):
                    psi = 2.0 * math.pi * i_psi / n_psi
                    cos_psi = math.cos(psi)
                    sin_psi = math.sin(psi)
                    # CCW-from-above (American convention), ψ=0 at +X.
                    # t_hat = [-sin(ψ), -cos(ψ), 0] in hub frame.
                    t_hat_ned = inputs.R_hub @ np.array([-sin_psi, -cos_psi, 0.0])
                    # v_t_extra = -v_inplane · t_hat (so advancing side gets +V).
                    v_t_extra = -float(np.dot(v_inplane, t_hat_ned))

                    lam_local = lambda_total + x_arr * (lam_c * cos_psi + lam_s * sin_psi)
                    v_a = lam_local * Omega_R                       # (n,)
                    v_t = v_t_base + v_t_extra                       # (n,)
                    valid = v_t > 0.0                                # reverse-flow mask
                    if not valid.any():
                        continue
                    phi   = np.arctan2(v_a, np.where(valid, v_t, 1.0))
                    alpha = (col + twist) - phi
                    cl, cd = self._polar.cl_cd_arr(alpha)
                    cn = cl * np.cos(phi) - cd * np.sin(phi)
                    ct = cl * np.sin(phi) + cd * np.cos(phi)
                    q  = 0.5 * rho * (v_a*v_a + v_t*v_t) * chord * dr * n_b
                    dT = np.where(valid, q * cn, 0.0)
                    dQ = np.where(valid, q * ct * r_arr, 0.0)
                    T_acc += float(dT.sum())
                    Q_acc += float(dQ.sum())

                T_total = T_acc / n_psi
                Q_total = Q_acc / n_psi

            else:
                # Axial flight: uniform inflow over the disk.  Fully vectorized.
                v_a = lambda_total * Omega_R                         # scalar
                v_t = omega * r_arr                                  # (n,)
                phi   = np.arctan2(v_a, v_t)
                alpha = (col + twist) - phi
                cl, cd = self._polar.cl_cd_arr(alpha)
                cn = cl * np.cos(phi) - cd * np.sin(phi)
                ct = cl * np.sin(phi) + cd * np.cos(phi)
                q  = 0.5 * rho * (v_a*v_a + v_t*v_t) * chord * dr * n_b
                T_total = float((q * cn).sum())
                Q_total = float((q * ct * r_arr).sum())

        # ------------------------------------------------------------------
        # Pitt-Peters inflow ODE
        # ------------------------------------------------------------------
        T_pos = max(T_total, 0.0)
        V_h   = math.sqrt(T_pos / (2.0 * rho * A)) if T_pos > 1e-6 else 0.0

        v0 = lam0 * Omega_R
        # V_T: total flow speed through disk — includes in-plane component
        # in forward flight.  Floor prevents τ → ∞ in VRS.
        V_T = max(
            math.sqrt(v_edge**2 + (v_climb + v0)**2),
            1e-2 * max(Omega_R, 1.0),
        )

        # λ₂ = V_descent / V_h  (> 0 during descent)
        V_c  = max(-v_climb, 0.0)
        lam2 = (V_c / V_h) if V_h > 1e-3 else 0.0

        # Steady-state uniform inflow
        if v_climb < -1e-3 and 0.0 < lam2 < 2.0:
            lam0_ss = (vrs_lambda1(lam2) * V_h / Omega_R) if Omega_R > 1e-6 else 0.0
        else:
            lam0_ss = T_total / (2.0 * rho * A * V_T * Omega_R) if Omega_R > 1e-6 else 0.0

        # Cyclic steady-state targets — Glauert skewed-wake model.
        # CCW-from-above, ψ=0 at +X, r_hat = [cos(ψ), -sin(ψ), 0].
        # Max inflow at back of disk gives:
        #   λ_c_ss = +µ_x · tan(χ/2),  λ_s_ss = −µ_y · tan(χ/2)
        # Asymmetric sign reflects the y-flip in r_hat (see CLAUDE.md).
        # µ here is wind-relative (v_inplane_hub / Ω_R), opposite-sign of
        # the vehicle-advance-ratio convention used in most texts.
        mu_sq = mu_x**2 + mu_y**2
        lam_sq = lambda_total**2
        denom = math.sqrt(mu_sq + lam_sq) + max(abs(lambda_total), 1e-6)
        tan_half_chi = math.sqrt(mu_sq) / denom if mu_sq > 1e-8 else 0.0
        lam_c_ss = +mu_x * tan_half_chi
        lam_s_ss = -mu_y * tan_half_chi

        # Apparent-mass time constants (Peters & HaQuang 1988)
        tau_0  = (8.0 * R) / (3.0 * math.pi * V_T)
        tau_cs = (16.0 * R) / (45.0 * math.pi * V_T)

        d_lam0  = (lam0_ss  - lam0)  / tau_0
        d_lam_c = (lam_c_ss - lam_c) / tau_cs
        d_lam_s = (lam_s_ss - lam_s) / tau_cs

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
