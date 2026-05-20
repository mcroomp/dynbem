"""Level 2 BEM with Øye 2-stage annular dynamic inflow.

Annulus-local alternative to Pitt-Peters.  Each radial annulus has its
own pair of first-order lag filters chasing the quasi-steady momentum
target; **no global L-matrix coupling** between annuli, no λ_c/λ_s
harmonic states.  Numerically much friendlier than Pitt-Peters at high
advance ratios and in descent + edgewise wind, where the
Pitt-Peters BEM-feedback (λ_c → BEM → C_L_hub → λ_s_ss) is stiff.

Inflow ODE (per annulus i, Øye 1990; OpenFAST AD Theory §6.3.4):

    τ₁ · dW_int/dt + W_int = W_qs + k · τ₁ · dW_qs/dt
    τ₂ · dW/dt     + W     = W_int

with coupling constant k = 0.6 (empirical, Øye original) and

    τ₁ = 1.1 / (1 - 1.3·min(a, 0.5)) · R / V_inf
    τ₂(r) = (0.39 - 0.26·(r/R)²) · τ₁

The blade sees W (not W_int, not W_qs).  This module's DBEMT_Mod
equivalent is OpenFAST's Mod=1: dW_qs/dt is treated as zero across the
outer integration step (constant τ assumption), which is exact for
operating-envelope sweeps and a small approximation for fast transients.

ψ-loop
------
Identical structure to PittPetersModelJIT._fwd_forces, except the local
inflow at each (r, ψ) is simply ``λ_climb + W[r]`` — no harmonic
decomposition.  Cyclic-frequency variation in blade loading still
emerges naturally because cyclic pitch and v_t_extra make local AoA
vary with ψ, and the hub moments Mx/My are still accumulated from the
azimuthally-varying dT.  But the *inflow itself* doesn't have explicit
cosine/sine harmonic states.

References
----------
Øye, S. (1990).  A simple vortex model.  IEA Symposium.
Snel, H. & Schepers, J.G. (1995).  Joint investigation of dynamic
  inflow effects and implementation of an engineering method.  ECN.
OpenFAST AeroDyn Theory v3.5, §6.3.4 (DBEMT).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numba import njit

from . import AeroBase, AeroResult, RotorInputs
from ._bem_common import (
    _interp_polar,
    build_polar_arrays,
    radial_grid,
    vrs_lambda1,
)
from .cyclic import cyclic_coeffs
from .rotor_definition import RotorDefinition
from .rotor_state import OyeRotorState, RotorState
from .polar import AirfoilPolar, LinearPolar, TabulatedPolar


# Øye coupling constant.  0.6 is the value Øye originally proposed and
# the OpenFAST default; it damps the coupling between W_qs and W_int.
_OYE_K = 0.6


# ---------------------------------------------------------------------------
# JIT'd ψ-loop kernel
# ---------------------------------------------------------------------------

@njit(cache=True, fastmath=True)
def _oye_psi_loop(r_arr, dr, chord, twist, col, omega, Omega_R,
                  lambda_climb, W, rho, n_b, n_psi,
                  v_in_hub_x, v_in_hub_y,
                  theta_1c, theta_1s,
                  alpha_tab, cl_tab, cd_tab,
                  out_dT_avg):
    """ψ-loop with per-annulus inflow W[r].

    Mirrors PittPetersModelJIT._fwd_forces but the local inflow comes
    from the annulus filter W[r] instead of harmonic decomposition.

    Out-parameter (pre-allocated by caller):
        out_dT_avg[i]  azimuth-averaged dT per annulus i (N)

    Returns (T_total, Q_total, Mx_hub, My_hub).  Each is azimuth-averaged
    (= 1/n_psi sum, matching the existing convention).
    """
    n_r = r_arr.shape[0]
    for i in range(n_r):
        out_dT_avg[i] = 0.0

    T_acc  = 0.0
    Q_acc  = 0.0
    Mx_acc = 0.0
    My_acc = 0.0
    inv_n_psi = 1.0 / n_psi

    for ipsi in range(n_psi):
        psi = 2.0 * math.pi * ipsi / n_psi
        cos_psi = math.cos(psi)
        sin_psi = math.sin(psi)
        v_t_extra = v_in_hub_x * sin_psi + v_in_hub_y * cos_psi
        col_psi = col + theta_1c * cos_psi + theta_1s * sin_psi
        rdT_sum = 0.0
        for i in range(n_r):
            r = r_arr[i]
            v_t = omega * r + v_t_extra
            if v_t <= 0.0:
                continue
            # Local inflow: uniform across annulus, no harmonics.
            lam_local = lambda_climb + W[i]
            v_a = lam_local * Omega_R
            phi = math.atan2(v_a, v_t)
            alpha = col_psi + twist - phi
            cl, cd = _interp_polar(alpha, alpha_tab, cl_tab, cd_tab)
            cos_p = math.cos(phi)
            sin_p = math.sin(phi)
            cn = cl * cos_p - cd * sin_p
            ct = cl * sin_p + cd * cos_p
            # No Prandtl tip loss in this JIT kernel — matches the
            # convention in pitt_peters_jit.  W_qs solver below uses
            # F=1 to stay consistent.
            q  = 0.5 * rho * (v_a*v_a + v_t*v_t) * chord * dr * n_b
            dT_elem = q * cn
            T_acc  += dT_elem
            Q_acc  += q * ct * r
            rdT_sum += r * dT_elem
            out_dT_avg[i] += dT_elem * inv_n_psi
        Mx_acc += rdT_sum * sin_psi
        My_acc += rdT_sum * cos_psi

    return T_acc * inv_n_psi, Q_acc * inv_n_psi, Mx_acc * inv_n_psi, My_acc * inv_n_psi


# ---------------------------------------------------------------------------
# Quasi-steady W_qs from per-annulus momentum balance
# ---------------------------------------------------------------------------

def _solve_W_qs(dT_avg: np.ndarray, x_arr: np.ndarray,
                dr: float, R: float, Omega_R: float, mu_T: float,
                rho: float) -> np.ndarray:
    """Quasi-steady induced inflow ratio per annulus, Glauert form.

    Annulus momentum with Glauert's resultant-velocity mass flow:

        dT = 4·π·r·dr·ρ·V_resultant·v_i

    where V_resultant = √(v_edge² + (v_climb + v_i)²) ≈ Ω·R·μ_T uses
    the rotor-mean μ_T (computed externally from the converged state).
    Linearised in W_qs, this gives the closed-form annulus solution

        W_qs[i] = dCT/dx[i] / (4·x[i]·μ_T)

    with dCT/dx = dT·R / (ρ·π·R²·Ω_R²·dr).  This matches the steady-
    state target Pitt-Peters uses (lam0_ss = T/(2ρA·V_T·Ω_R) in
    aggregate form) and stays stable in forward flight where the pure
    axial-momentum form  4·x·λ_r·W = dCT/dx  produces an enormous W
    when λ_r is small (the autorotation / VRS attractor).

    No Prandtl loss factor: matches the JIT ψ-loop convention.

    Returns W_qs array, same length as the radial grid.
    """
    if Omega_R < 1e-6:
        return np.zeros_like(dT_avg)
    A = math.pi * R * R
    rho_norm = rho * A * Omega_R * Omega_R * dr / R
    dCdx = dT_avg / max(rho_norm, 1e-30)
    W_qs = dCdx / (4.0 * np.maximum(x_arr, 1e-6) * max(mu_T, 0.05))
    # Sanity clamp — physical W is small; clip avoids NaN if a transient
    # state produces wild dT.
    return np.clip(W_qs, -1.0, 1.0)


def _oye_taus(R: float, x_arr: np.ndarray, V_inf: float, a_avg: float
              ) -> tuple[np.ndarray, np.ndarray]:
    """Øye time constants per annulus.

    τ₁ = 1.1 / (1 − 1.3·min(a_avg, 0.5)) · R / V_inf      (rotor-mean)
    τ₂(r) = (0.39 − 0.26·(r/R)²) · τ₁                      (radius-dependent)

    The 0.5 clamp on a_avg keeps τ₁ finite as induction approaches the
    actuator-disk limit a=0.5.  Matches OpenFAST DBEMT_Mod=1.
    """
    a_clamped = max(0.0, min(0.5, a_avg))
    tau1 = 1.1 / (1.0 - 1.3 * a_clamped) * R / max(V_inf, 1e-2)
    tau2_arr = (0.39 - 0.26 * x_arr * x_arr) * tau1
    tau1_arr = np.full_like(x_arr, tau1)
    return tau1_arr, tau2_arr


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

@dataclass
class OyeBEMModel(AeroBase):
    """Level 2 BEM with Øye 2-stage annular dynamic inflow.

    Stable alternative to PittPetersModel for high-advance-ratio and
    descent + edgewise operating regimes where Pitt-Peters' global
    L-matrix coupling can be numerically stiff.  Trade-off: doesn't
    capture cyclic-driven λ_c/λ_s harmonics explicitly — the inflow
    response to hub moments is purely through per-annulus loading
    asymmetry (which the ψ-loop captures), not through a coupled L
    matrix.
    """

    defn: RotorDefinition
    n_psi_elements: int = 36
    coupling_k:     float = _OYE_K
    _polar: AirfoilPolar = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if self.defn.airfoil.polar_csv is not None:
            self._polar = TabulatedPolar.from_csv(self.defn.airfoil.polar_csv)
        else:
            self._polar = LinearPolar.from_properties(self.defn.airfoil)
        # Shared with PittPetersModel(JIT) — see dynbem/_bem_common.py.
        (self._dr, self._r_mid, self._x_mid, self._x_hub,
         self._twist_rad) = radial_grid(self.defn.blade)
        self._alpha_tab, self._cl_tab, self._cd_tab = build_polar_arrays(self._polar)
        self._n_elements = self.defn.blade.n_elements
        self._use_tip_loss = bool(self.defn.airfoil.tip_loss)

        # Scratch buffer for the JIT kernel — avoid per-call allocation.
        self._dT_avg_buf = np.zeros(self._n_elements, dtype=np.float64)

    # -----------------------------------------------------------------------

    def initial_rotor_state(self) -> OyeRotorState:
        return OyeRotorState.zeros(self._n_elements)

    def compute_forces(
        self,
        inputs: RotorInputs,
        state:  RotorState,
    ) -> tuple[AeroResult, RotorState]:
        assert isinstance(state, OyeRotorState)
        blade   = self.defn.blade
        omega   = state.omega_rad_s
        rho     = inputs.rho_kg_m3
        R       = blade.radius_m

        Omega_R = omega * R
        hub_axis = inputs.R_hub @ np.array([0.0, 0.0, 1.0])
        v_rel    = inputs.wind_world - inputs.v_hub_world
        v_climb  = float(np.dot(v_rel, hub_axis))
        v_inplane = v_rel - v_climb * hub_axis
        v_edge   = float(np.linalg.norm(v_inplane))
        mu       = v_edge / max(Omega_R, 1e-6)
        v_inplane_hub = inputs.R_hub.T @ v_inplane

        # Cyclic pitch coefficients.
        theta_1c, theta_1s = cyclic_coeffs(
            inputs.tilt_lon, inputs.tilt_lat, self.defn.control
        )
        has_cyclic = abs(theta_1c) + abs(theta_1s) > 1e-12

        # -------------------------------------------------------------------
        # Blade element forces — JIT ψ-loop with per-annulus W
        # -------------------------------------------------------------------
        T_total = Q_total = 0.0
        Mx_hub = My_hub = 0.0
        lambda_climb = v_climb / Omega_R if Omega_R > 1e-6 else 0.0

        if Omega_R > 1e-6 and omega > 1.0:
            T_total, Q_total, Mx_hub, My_hub = _oye_psi_loop(
                self._r_mid, self._dr,
                blade.chord_m, self._twist_rad,
                inputs.collective_rad, omega, Omega_R,
                lambda_climb, state.W, rho, blade.n_blades,
                self.n_psi_elements,
                float(v_inplane_hub[0]), float(v_inplane_hub[1]),
                theta_1c, theta_1s,
                self._alpha_tab, self._cl_tab, self._cd_tab,
                self._dT_avg_buf,
            )

        # -------------------------------------------------------------------
        # W_qs from per-annulus momentum balance (Glauert form), VRS override
        # -------------------------------------------------------------------
        # Rotor-mean μ_T = V_T / Ω·R, where V_T is the resultant velocity
        # through the disk.  Matches Pitt-Peters' aggregate denominator.
        v0_mean = float(np.mean(state.W)) * Omega_R
        V_T = max(math.sqrt(v_edge*v_edge + (v_climb + v0_mean)**2),
                  1e-2 * max(Omega_R, 1.0))
        mu_T = V_T / max(Omega_R, 1e-6)

        # Vortex Ring State: momentum theory diverges when 0 < V_descent/V_h < 2.
        # Use the Leishman empirical polynomial (same as PittPetersModel) and
        # apply it uniformly across annuli — Øye doesn't resolve VRS locally.
        rho_A = rho * math.pi * R * R
        T_pos = max(T_total, 0.0)
        V_h   = math.sqrt(T_pos / (2.0 * rho_A)) if T_pos > 1e-6 else 0.0
        V_c   = max(-v_climb, 0.0)
        lam2  = (V_c / V_h) if V_h > 1e-3 else 0.0
        if v_climb < -1e-3 and 0.0 < lam2 < 2.0 and Omega_R > 1e-6:
            W_uniform = vrs_lambda1(lam2) * V_h / Omega_R
            W_qs = np.full(self._n_elements, W_uniform)
        else:
            W_qs = _solve_W_qs(
                self._dT_avg_buf, self._x_mid,
                self._dr, R, Omega_R, mu_T, rho,
            )

        # -------------------------------------------------------------------
        # Øye 2-stage filter ODE
        # -------------------------------------------------------------------
        # τ₁, τ₂ use the same V_T (=mass-flow speed at the disk) and rotor-
        # mean a as the W_qs solver above.
        a_avg = float(np.mean(state.W) * Omega_R / V_T) if V_T > 1e-3 else 0.0
        tau1, tau2 = _oye_taus(R, self._x_mid, V_T, a_avg)

        # Mod=1: treat dW_qs/dt = 0 across the outer step (OpenFAST default
        # for envelope-style sweeps; the integrator's small dt makes this
        # exact for steady operating points).
        dW_int_dt = (W_qs - state.W_int) / tau1
        dW_dt     = (state.W_int - state.W) / tau2

        # -------------------------------------------------------------------
        # Mechanical states + assemble outputs
        # -------------------------------------------------------------------
        I_ode = (
            self.defn.autorotation.I_ode_kgm2
            if self.defn.autorotation.I_ode_kgm2 is not None
            else 1.0
        )
        d_omega = (-Q_total + inputs.motor_torque_Nm) / I_ode
        d_psi   = omega

        F_world   = -T_total * hub_axis
        M_orbital = inputs.R_hub @ np.array([Mx_hub, My_hub, 0.0])
        result = AeroResult(
            F_world=F_world,
            M_orbital=M_orbital,
            Q_spin=Q_total,
            M_spin=Q_total * hub_axis,
        )
        derivative = OyeRotorState(
            W_int=dW_int_dt,
            W=dW_dt,
            omega_rad_s=d_omega,
            spin_angle_rad=d_psi,
        )
        # Silence unused-variable lints — `mu`, `has_cyclic` are kept for
        # potential downstream gating but aren't currently used (Øye runs
        # the same ψ-loop in all regimes).
        _ = (mu, has_cyclic)
        return result, derivative

    # -----------------------------------------------------------------------

    def inflow_taus(self, inputs: RotorInputs, state: RotorState) -> np.ndarray:
        """Time constants for [W_int(n), W(n), ω, ψ].

        W_int uses τ₁, W uses τ₂(r); mechanical states get ∞.  Used by
        the envelope integrator's semi-implicit damping.
        """
        assert isinstance(state, OyeRotorState)
        omega = state.omega_rad_s
        R     = self.defn.blade.radius_m
        Omega_R = omega * R
        if Omega_R < 1e-6:
            n = self._n_elements
            return np.concatenate([
                np.full(n, np.inf), np.full(n, np.inf),
                [np.inf, np.inf],
            ])

        hub_axis = inputs.R_hub @ np.array([0.0, 0.0, 1.0])
        v_rel    = inputs.wind_world - inputs.v_hub_world
        v_climb  = float(np.dot(v_rel, hub_axis))
        v_edge   = float(np.linalg.norm(v_rel - v_climb * hub_axis))
        V_inf = max(math.sqrt(v_edge*v_edge + (v_climb + np.mean(state.W)*Omega_R)**2),
                    1e-2 * max(Omega_R, 1.0))
        a_avg = float(np.mean(state.W) * Omega_R / V_inf) if V_inf > 1e-3 else 0.0
        tau1, tau2 = _oye_taus(R, self._x_mid, V_inf, a_avg)
        return np.concatenate([tau1, tau2, [np.inf, np.inf]])

    def to_dict(self) -> dict:
        return {
            "model": "Oye_BEM_Level2",
            "n_elements": self.defn.blade.n_elements,
            "coupling_k": self.coupling_k,
        }

    @classmethod
    def from_definition(cls, defn: RotorDefinition) -> "OyeBEMModel":
        return cls(defn=defn)
