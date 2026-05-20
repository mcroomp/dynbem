"""Level 2: Pitt-Peters 3-state dynamic inflow with VRS empirical correction.

Inflow distribution (uniform + linear tilt):
    λ(r, ψ) = λ_0 + (r/R)·(λ_c·cos ψ + λ_s·sin ψ)

Inflow ODE — first-order lags from Peters & HaQuang (1988):
    dλ_0/dt = (λ_0_ss − λ_0) / τ_0     τ_0  = 8R / (3π·V_T)
    dλ_c/dt = (λ_c_ss − λ_c) / τ_cs    τ_cs = 16R / (45π·V_T)
    dλ_s/dt = (λ_s_ss − λ_s) / τ_cs

Steady-state targets — Pitt-Peters L matrix (Peters 2009 Eq 10) with
X = tan(χ/2), translated to our ψ=0-at-+X convention:
    λ_0_ss = C_T/(2·µ_T)       + (15π·X/64) · C_M_hub / µ_T
    λ_c_ss = −(15π·X/64) · C_T + 4·cos(χ) / (1+cos χ)  · C_M_hub  / µ_T
    λ_s_ss =                     4         / (1+cos χ) · C_L_hub  / µ_T

The −(15π·X/64)·C_T term in λ_c_ss is the cross-coupling that produces
Glauert wake-skew naturally from thrust forcing.  In the VRS regime
(0 < V_descent/V_h < 2) λ_0_ss is overridden by the Leishman (2000)
empirical polynomial; cross-coupling is then skipped (momentum theory
invalid in recirculating wake).

Forward-flight blade element loop averages over n_psi_elements azimuth
stations using the full λ(r,ψ) inflow distribution at each station.
V_T uses total disk velocity √(v_edge² + (v_climb + v_i)²) in forward
flight.

For numerical stability at high advance ratios and descent + edgewise
wind regimes, the BEM-driven feedback through the L matrix can be
stiff — see ``dynbem/oye.py`` for the alternative annulus-local
formulation that avoids the global L coupling.

See CLAUDE.md "Pitt-Peters inflow ODE" for the canonical L-matrix
formulation and sign conventions.

References
----------
Peters, D.A. (2009).  How Dynamic Inflow Survives in the Competitive
  World of Rotorcraft Aerodynamics: The Alexander Nikolsky Honorary
  Lecture.  JAHS 54(1):011001.
Pitt & Peters (1981), Vertica 5(1), 21-34.
Peters & HaQuang (1988), JAHS 33(4), 64-68.
Leishman (2000), §12.4, §12.7.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from . import AeroBase, AeroResult, RotorInputs
from ._bem_common import radial_grid, vrs_lambda1
from .bem import prandtl_tip_loss
from .cyclic import cyclic_coeffs
from .rotor_definition import RotorDefinition
from .rotor_state import PittPetersRotorState, RotorState
from .polar import AirfoilPolar, LinearPolar, TabulatedPolar


# Re-export so external callers (tests, docs) can still do
# ``from dynbem.pitt_peters import vrs_lambda1``.
__all__ = ["PittPetersModel", "vrs_lambda1", "prescribed_element_forces"]


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

    The ψ-loop integrates over n_psi_elements azimuth stations whenever
    µ > 0.01, cyclic input is present, or the cyclic inflow state is
    nonzero — using the full λ(r,ψ) distribution.  Steady-state cyclic
    targets come from the Pitt-Peters L matrix (Peters 2009 Eq 10).

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
        # Cache fixed radial geometry — shared with OyeBEMModel and
        # PittPetersModelJIT via dynbem/_bem_common.py.
        (self._dr, self._r_mid, self._x_mid, self._x_hub_unused,
         self._twist_rad) = radial_grid(self.defn.blade)

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
        # In-plane hub-frame moments from non-axisymmetric thrust
        # (forward flight, cyclic inflow distribution, and/or cyclic pitch).
        Mx_hub = 0.0
        My_hub = 0.0

        # Cyclic pitch coefficients: θ_cyclic(ψ) = θ_1c·cos(ψ) + θ_1s·sin(ψ).
        theta_1c, theta_1s = cyclic_coeffs(
            inputs.tilt_lon, inputs.tilt_lat, self.defn.control
        )
        has_cyclic = abs(theta_1c) + abs(theta_1s) > 1e-12

        if Omega_R > 1e-6:
            if (mu > 0.01 or has_cyclic or abs(lam_c) + abs(lam_s) > 1e-12) and omega > 1.0:
                # ψ-loop, vectorized over r at each station.  Triggered by
                # forward flight, nonzero cyclic, or nonzero λ_c/λ_s state.
                n_psi = self.n_psi_elements
                T_acc = 0.0
                Q_acc = 0.0
                Mx_acc = 0.0
                My_acc = 0.0
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

                    col_psi = col + theta_1c * cos_psi + theta_1s * sin_psi

                    lam_local = lambda_total + x_arr * (lam_c * cos_psi + lam_s * sin_psi)
                    v_a = lam_local * Omega_R                       # (n,)
                    v_t = v_t_base + v_t_extra                       # (n,)
                    valid = v_t > 0.0                                # reverse-flow mask
                    if not valid.any():
                        continue
                    phi   = np.arctan2(v_a, np.where(valid, v_t, 1.0))
                    alpha = (col_psi + twist) - phi
                    cl, cd = self._polar.cl_cd_arr(alpha)
                    cn = cl * np.cos(phi) - cd * np.sin(phi)
                    ct = cl * np.sin(phi) + cd * np.cos(phi)
                    q  = 0.5 * rho * (v_a*v_a + v_t*v_t) * chord * dr * n_b
                    dT = np.where(valid, q * cn, 0.0)
                    dQ = np.where(valid, q * ct * r_arr, 0.0)
                    dT_sum = float(dT.sum())
                    T_acc  += dT_sum
                    Q_acc  += float(dQ.sum())
                    # Per-element moment: r·dT·[sin(ψ), cos(ψ), 0] in hub frame.
                    rdT_sum = float((dT * r_arr).sum())
                    Mx_acc += rdT_sum * sin_psi
                    My_acc += rdT_sum * cos_psi

                T_total = T_acc / n_psi
                Q_total = Q_acc / n_psi
                Mx_hub  = Mx_acc / n_psi
                My_hub  = My_acc / n_psi

            else:
                # Axial flight, no cyclic, no cyclic inflow: axisymmetric.
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

        # ------------------------------------------------------------------
        # Pitt-Peters L-matrix steady-state targets (hub axes).
        # See Peters (2009) Eq 10; signs translated to our ψ=0-at-+X
        # convention.  At hover χ=0, so L_off=0 and L_cc=L_ss=2 — cyclic
        # forcing (C_M_hub, C_L_hub) drives λ_c_ss, λ_s_ss through the
        # diagonal; in forward flight the off-diagonal −L_off·C_T term
        # reproduces Glauert wake-skew from thrust forcing.
        #
        # NOT in wind axes: treats (C_M_hub, C_L_hub) as wind-axis directly,
        # exact for axial/longitudinal flight; approximate for µ_y ≠ 0.
        # See CLAUDE.md "Wind-axis rotation" for the reverted attempt.
        # ------------------------------------------------------------------
        mu_T_eff = max(V_T / Omega_R if Omega_R > 1e-6 else 0.0, 0.05)
        mu_inplane = v_edge / max(Omega_R, 1e-6)
        chi = math.atan2(mu_inplane, abs(lambda_total) + 1e-6)
        cos_chi = math.cos(chi)
        tan_half_chi = math.tan(0.5 * chi)
        L_off = (15.0 * math.pi / 64.0) * tan_half_chi
        L_cc  = 4.0 * cos_chi / (1.0 + cos_chi)
        L_ss  = 4.0 / (1.0 + cos_chi)

        norm = rho * A * Omega_R * R * V_T
        C_L_hub = Mx_hub / norm if norm > 1e-9 else 0.0
        C_M_hub = My_hub / norm if norm > 1e-9 else 0.0

        if v_climb < -1e-3 and 0.0 < lam2 < 2.0:
            lam0_ss = (vrs_lambda1(lam2) * V_h / Omega_R) if Omega_R > 1e-6 else 0.0
        else:
            lam0_ss = (
                T_total / (2.0 * rho * A * V_T * Omega_R)
                + L_off * C_M_hub / mu_T_eff
            ) if Omega_R > 1e-6 else 0.0

        C_T = T_total / (rho * A * Omega_R * Omega_R) if Omega_R > 1e-6 else 0.0
        lam_c_ss = (-L_off * C_T + L_cc * C_M_hub) / mu_T_eff
        lam_s_ss = ( L_ss * C_L_hub) / mu_T_eff

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
        F_world   = -T_total * hub_axis
        M_orbital = inputs.R_hub @ np.array([Mx_hub, My_hub, 0.0])
        result  = AeroResult(
            F_world=F_world,
            M_orbital=M_orbital,
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

    def inflow_taus(self, inputs, state):
        """Time constants per state component: [τ_0, τ_cs, τ_cs, ∞, ∞].

        Mirrors the apparent-mass formulas used inside compute_forces:
            τ_0  = 8R/(3π·V_T),  τ_cs = 16R/(45π·V_T)
        with V_T = √(v_edge² + (v_climb + v_0)²) and a hover floor.  Used
        by the envelope integrator's semi-implicit damping.
        """
        omega   = state.omega_rad_s
        R       = self.defn.blade.radius_m
        Omega_R = omega * R
        hub_axis = inputs.R_hub @ np.array([0.0, 0.0, 1.0])
        v_rel    = inputs.wind_world - inputs.v_hub_world
        v_climb  = float(np.dot(v_rel, hub_axis))
        v_edge   = float(np.linalg.norm(v_rel - v_climb * hub_axis))
        v0       = state.lambda_0 * Omega_R
        V_T = max(math.sqrt(v_edge*v_edge + (v_climb + v0)**2),
                  1e-2 * max(Omega_R, 1.0))
        tau_0  = (8.0 * R) / (3.0 * math.pi * V_T)
        tau_cs = (16.0 * R) / (45.0 * math.pi * V_T)
        return np.array([tau_0, tau_cs, tau_cs, np.inf, np.inf])

    def to_dict(self) -> dict:
        return {"model": "PittPeters_Level2", "n_elements": self.defn.blade.n_elements}

    @classmethod
    def from_definition(cls, defn: RotorDefinition) -> "PittPetersModel":
        return cls(defn=defn)
