"""Numba-JIT'd Pitt-Peters variant — same physics as PittPetersModel.

Mirrors ``PittPetersModel.compute_forces`` exactly except the blade-element
radial loop (axial branch) and the (psi x r) loop (forward-flight branch)
are compiled with ``@njit``.  The Pitt-Peters ODE update, VRS polynomial,
and inflow geometry are identical.  Use ``tests/test_pitt_peters_jit.py``
to validate numerical agreement against the reference implementation.

Why a separate module: keeps the well-validated PittPetersModel intact as
the reference, while letting the JIT version evolve freely.  When the JIT
version is stable, the original can be deprecated.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numba import njit

from . import AeroBase, AeroResult, RotorInputs
from .cyclic import cyclic_coeffs
from .pitt_peters import vrs_lambda1
from .rotor_definition import RotorDefinition
from .rotor_state import PittPetersRotorState, RotorState
from .polar import AirfoilPolar, LinearPolar, TabulatedPolar


# ---------------------------------------------------------------------------
# JIT'd kernels
# ---------------------------------------------------------------------------

@njit(cache=True, fastmath=True)
def _interp_polar(alpha, alpha_tab, cl_tab, cd_tab):
    """Scalar binary-search interpolation (cl, cd) at angle alpha (rad).

    Clamps to endpoints outside the tabulated range, matching np.interp.
    """
    n = alpha_tab.shape[0]
    if alpha <= alpha_tab[0]:
        return cl_tab[0], cd_tab[0]
    if alpha >= alpha_tab[n - 1]:
        return cl_tab[n - 1], cd_tab[n - 1]
    lo = 0
    hi = n - 1
    while hi - lo > 1:
        mid = (lo + hi) >> 1
        if alpha_tab[mid] <= alpha:
            lo = mid
        else:
            hi = mid
    a_lo = alpha_tab[lo]
    a_hi = alpha_tab[hi]
    t = (alpha - a_lo) / (a_hi - a_lo)
    return (
        cl_tab[lo] + t * (cl_tab[hi] - cl_tab[lo]),
        cd_tab[lo] + t * (cd_tab[hi] - cd_tab[lo]),
    )


@njit(cache=True, fastmath=True)
def _axial_forces(r_arr, dr, chord, twist, col, omega, Omega_R,
                  lambda_total, rho, n_b,
                  alpha_tab, cl_tab, cd_tab):
    """Sum thrust and torque over radial elements, axial flight (mu = 0)."""
    T = 0.0
    Q = 0.0
    v_a = lambda_total * Omega_R
    for i in range(r_arr.shape[0]):
        r = r_arr[i]
        v_t = omega * r
        if v_t < 1e-9:
            continue
        phi = math.atan2(v_a, v_t)
        alpha = col + twist - phi
        cl, cd = _interp_polar(alpha, alpha_tab, cl_tab, cd_tab)
        cos_p = math.cos(phi)
        sin_p = math.sin(phi)
        cn = cl * cos_p - cd * sin_p
        ct = cl * sin_p + cd * cos_p
        q = 0.5 * rho * (v_a * v_a + v_t * v_t) * chord * dr * n_b
        T += q * cn
        Q += q * ct * r
    return T, Q


@njit(cache=True, fastmath=True)
def _fwd_forces(r_arr, x_arr, dr, chord, twist, col, omega, Omega_R,
                lambda_total, lam_c, lam_s, rho, n_b, n_psi,
                v_in_hub_x, v_in_hub_y,
                theta_1c, theta_1s,
                alpha_tab, cl_tab, cd_tab):
    """Sum thrust, torque, and in-plane hub moments over (psi, r).

    v_in_hub_x, v_in_hub_y are the in-plane wind components in hub frame.
    theta_1c, theta_1s are the cyclic pitch coefficients
    (theta(psi) = collective + theta_1c*cos(psi) + theta_1s*sin(psi)).

    CCW-from-above (American convention), psi=0 at +X.
    Tangential direction in hub frame: t_hat = [-sin(psi), -cos(psi), 0].
    v_t_extra = -v_inplane . t_hat = +v_in_hub_x*sin(psi) + v_in_hub_y*cos(psi).
    Per-element moment: dM_hub = r*dT*[sin(psi), cos(psi), 0].

    Returns (T, Q, Mx_hub, My_hub), each azimuth-averaged.
    """
    T_acc = 0.0
    Q_acc = 0.0
    Mx_acc = 0.0
    My_acc = 0.0
    for ipsi in range(n_psi):
        psi = 2.0 * math.pi * ipsi / n_psi
        cos_psi = math.cos(psi)
        sin_psi = math.sin(psi)
        v_t_extra = v_in_hub_x * sin_psi + v_in_hub_y * cos_psi
        col_psi = col + theta_1c * cos_psi + theta_1s * sin_psi
        rdT_sum = 0.0
        for i in range(r_arr.shape[0]):
            r = r_arr[i]
            v_t = omega * r + v_t_extra
            if v_t <= 0.0:
                continue
            x = x_arr[i]
            lam_local = lambda_total + x * (lam_c * cos_psi + lam_s * sin_psi)
            v_a = lam_local * Omega_R
            phi = math.atan2(v_a, v_t)
            alpha = col_psi + twist - phi
            cl, cd = _interp_polar(alpha, alpha_tab, cl_tab, cd_tab)
            cos_p = math.cos(phi)
            sin_p = math.sin(phi)
            cn = cl * cos_p - cd * sin_p
            ct = cl * sin_p + cd * cos_p
            q = 0.5 * rho * (v_a * v_a + v_t * v_t) * chord * dr * n_b
            dT_elem = q * cn
            T_acc += dT_elem
            Q_acc += q * ct * r
            rdT_sum += r * dT_elem
        Mx_acc += rdT_sum * sin_psi
        My_acc += rdT_sum * cos_psi
    return T_acc / n_psi, Q_acc / n_psi, Mx_acc / n_psi, My_acc / n_psi


# ---------------------------------------------------------------------------
# Model class
# ---------------------------------------------------------------------------

@dataclass
class PittPetersModelJIT(AeroBase):
    """Pitt-Peters Level 2 with JIT-compiled blade-element loop.

    API and physics identical to PittPetersModel; only the radial integration
    is compiled.  Validated by tests/test_pitt_peters_jit.py.
    """

    defn: RotorDefinition
    n_psi_elements: int = 36
    _polar: AirfoilPolar = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if self.defn.airfoil.polar_csv is not None:
            self._polar = TabulatedPolar.from_csv(self.defn.airfoil.polar_csv)
        else:
            self._polar = LinearPolar.from_properties(self.defn.airfoil)

        blade = self.defn.blade
        R = blade.radius_m
        n = blade.n_elements
        r0 = blade.root_cutout_m
        self._dr = (R - r0) / n
        self._r_mid = np.ascontiguousarray(
            np.linspace(r0 + 0.5 * self._dr, R - 0.5 * self._dr, n)
        )
        self._x_mid = np.ascontiguousarray(self._r_mid / R)
        self._twist_rad = math.radians(blade.twist_deg)

        # Polar arrays passed into JIT kernels.  For tabulated polars we use
        # the existing arrays; for analytical polars we sample onto a uniform
        # grid (no behaviour change since LinearPolar is piecewise-affine).
        if isinstance(self._polar, TabulatedPolar):
            self._alpha_tab = np.ascontiguousarray(self._polar._alpha)
            self._cl_tab = np.ascontiguousarray(self._polar._cl)
            self._cd_tab = np.ascontiguousarray(self._polar._cd)
        else:
            n_p = 4001
            a = np.linspace(-math.pi / 2, math.pi / 2, n_p)
            cl = np.empty(n_p)
            cd = np.empty(n_p)
            for i in range(n_p):
                cl_i, cd_i = self._polar.cl_cd(float(a[i]))
                cl[i] = cl_i
                cd[i] = cd_i
            self._alpha_tab = np.ascontiguousarray(a)
            self._cl_tab = np.ascontiguousarray(cl)
            self._cd_tab = np.ascontiguousarray(cd)

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
        v_climb  = float(np.dot(v_rel, hub_axis))

        v_inplane = v_rel - v_climb * hub_axis
        v_edge    = float(np.linalg.norm(v_inplane))
        mu        = v_edge / max(Omega_R, 1e-6)

        v_inplane_hub = inputs.R_hub.T @ v_inplane
        mu_x = float(v_inplane_hub[0]) / max(Omega_R, 1e-6)
        mu_y = float(v_inplane_hub[1]) / max(Omega_R, 1e-6)

        lam0  = state.lambda_0
        lam_c = state.lambda_c
        lam_s = state.lambda_s

        lambda_climb = v_climb / Omega_R if Omega_R > 1e-6 else 0.0
        lambda_total = lam0 + lambda_climb

        # Cyclic pitch: θ(ψ) = collective + θ_1c·cos(ψ) + θ_1s·sin(ψ).
        theta_1c, theta_1s = cyclic_coeffs(
            inputs.tilt_lon, inputs.tilt_lat, self.defn.control
        )
        has_cyclic = abs(theta_1c) + abs(theta_1s) > 1e-12

        # ------------------------------------------------------------------
        # Blade element forces — JIT kernels
        # ------------------------------------------------------------------
        T_total = Q_total = 0.0
        Mx_hub = My_hub = 0.0
        if Omega_R > 1e-6:
            if (mu > 0.01 or has_cyclic or abs(lam_c) + abs(lam_s) > 1e-12) and omega > 1.0:
                T_total, Q_total, Mx_hub, My_hub = _fwd_forces(
                    self._r_mid, self._x_mid, self._dr, blade.chord_m,
                    self._twist_rad, inputs.collective_rad, omega, Omega_R,
                    lambda_total, lam_c, lam_s, rho, blade.n_blades,
                    self.n_psi_elements,
                    float(v_inplane_hub[0]), float(v_inplane_hub[1]),
                    theta_1c, theta_1s,
                    self._alpha_tab, self._cl_tab, self._cd_tab,
                )
            else:
                T_total, Q_total = _axial_forces(
                    self._r_mid, self._dr, blade.chord_m, self._twist_rad,
                    inputs.collective_rad, omega, Omega_R, lambda_total,
                    rho, blade.n_blades,
                    self._alpha_tab, self._cl_tab, self._cd_tab,
                )

        # ------------------------------------------------------------------
        # Pitt-Peters inflow ODE — copied verbatim from PittPetersModel
        # ------------------------------------------------------------------
        T_pos = max(T_total, 0.0)
        V_h   = math.sqrt(T_pos / (2.0 * rho * A)) if T_pos > 1e-6 else 0.0

        v0  = lam0 * Omega_R
        V_T = max(
            math.sqrt(v_edge**2 + (v_climb + v0)**2),
            1e-2 * max(Omega_R, 1.0),
        )

        V_c  = max(-v_climb, 0.0)
        lam2 = (V_c / V_h) if V_h > 1e-3 else 0.0

        if v_climb < -1e-3 and 0.0 < lam2 < 2.0:
            lam0_ss = (vrs_lambda1(lam2) * V_h / Omega_R) if Omega_R > 1e-6 else 0.0
        else:
            lam0_ss = T_total / (2.0 * rho * A * V_T * Omega_R) if Omega_R > 1e-6 else 0.0

        mu_sq = mu_x**2 + mu_y**2
        lam_sq = lambda_total**2
        denom = math.sqrt(mu_sq + lam_sq) + max(abs(lambda_total), 1e-6)
        tan_half_chi = math.sqrt(mu_sq) / denom if mu_sq > 1e-8 else 0.0
        # CCW-from-above, ψ=0 at +X: λ_c_ss = +mu_x·tan(χ/2),
        # λ_s_ss = -mu_y·tan(χ/2). See PittPetersModel for full derivation.
        lam_c_ss = +mu_x * tan_half_chi
        lam_s_ss = -mu_y * tan_half_chi

        tau_0  = (8.0 * R) / (3.0 * math.pi * V_T)
        tau_cs = (16.0 * R) / (45.0 * math.pi * V_T)

        d_lam0  = (lam0_ss  - lam0)  / tau_0
        d_lam_c = (lam_c_ss - lam_c) / tau_cs
        d_lam_s = (lam_s_ss - lam_s) / tau_cs

        I_ode = (
            self.defn.autorotation.I_ode_kgm2
            if self.defn.autorotation.I_ode_kgm2 is not None
            else 1.0
        )
        d_omega      = (-Q_total + inputs.motor_torque_Nm) / I_ode
        d_spin_angle = omega

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

    def to_dict(self) -> dict:
        return {"model": "PittPeters_Level2_JIT", "n_elements": self.defn.blade.n_elements}

    @classmethod
    def from_definition(cls, defn: RotorDefinition) -> "PittPetersModelJIT":
        return cls(defn=defn)
