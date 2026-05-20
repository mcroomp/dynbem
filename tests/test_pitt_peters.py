"""Validation of the Level-2 Pitt-Peters model (aero/pitt_peters.py).

Three checks:
1. Hover CT converges to within 5% of the Level-1 BEM quasi-static result.
   Both models share the same blade geometry and airfoil, so they should
   agree in steady hover.  Tight tolerance catches integration bugs.

2. VRS no-blow-up: at a shallow descent (lambda2 ~ 0.4) the Level-1 BEM
   instantly jumps to 2.4x CT; the Pitt-Peters VRS correction should hold
   CT below 1.6x at steady state.

3. WBS autorotation: in deep WBS (lambda2 > 2.5), CQ must be negative
   (rotor extracting energy from the descent wind — autorotation/turbine mode).
"""

import math
from pathlib import Path
import numpy as np
import pytest

from dynbem.pitt_peters import PittPetersModel, vrs_lambda1
from dynbem.bem import BEMModel
from dynbem import RotorInputs
import dynbem.rotor_definition as rotor_definition
from dynbem.rotor_state import PittPetersRotorState, QuasiStaticRotorState

_ROTOR_YAML = str(
    Path(__file__).parent.parent / "rotors" / "castles_gray_6ft" / "rotor.yaml"
)
_RHO = 1.225


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def defn():
    return rotor_definition.load(_ROTOR_YAML)


@pytest.fixture(scope="module")
def pp_model(defn):
    return PittPetersModel(defn=defn)


@pytest.fixture(scope="module")
def bem_model(defn):
    return BEMModel(defn=defn)


# ---------------------------------------------------------------------------
# Integration helper
# ---------------------------------------------------------------------------

def _euler_to_steady(model, theta_deg: float, rpm: float, v_climb_ms: float,
                     n_steps: int = 6000, dt: float = 0.001,
                     lam0_init: float = 0.0) -> tuple[float, float, float]:
    """Euler-integrate the Pitt-Peters inflow ODE to steady state.

    Returns (CT, CQ, lam0_final) at the end of integration.
    """
    omega = rpm * math.pi / 30.0
    R = model.defn.blade.radius_m
    A = math.pi * R**2

    inp = RotorInputs(
        collective_rad=math.radians(theta_deg),
        tilt_lon=0.0, tilt_lat=0.0, R_hub=np.eye(3),
        v_hub_world=np.zeros(3),
        wind_world=np.array([0.0, 0.0, v_climb_ms]),
        t=0.0,
    )
    lam0 = lam0_init
    for _ in range(n_steps):
        state = PittPetersRotorState(lambda_0=lam0, omega_rad_s=omega)
        res, drv = model.compute_forces(inp, state)
        lam0 += drv.lambda_0 * dt

    T = -res.F_world[2]
    CT = T / (_RHO * A * (omega * R) ** 2)
    CQ = res.Q_spin / (_RHO * A * (omega * R) ** 2 * R)
    return CT, CQ, lam0


def _bem_ct(model: BEMModel, theta_deg: float, rpm: float, v_climb_ms: float) -> float:
    omega = rpm * math.pi / 30.0
    R = model.defn.blade.radius_m
    A = math.pi * R**2
    inp = RotorInputs(
        collective_rad=math.radians(theta_deg),
        tilt_lon=0.0, tilt_lat=0.0, R_hub=np.eye(3),
        v_hub_world=np.zeros(3),
        wind_world=np.array([0.0, 0.0, v_climb_ms]),
        t=0.0,
    )
    res, _ = model.compute_forces(inp, QuasiStaticRotorState(omega_rad_s=omega))
    return -res.F_world[2] / (_RHO * A * (omega * R) ** 2)


# ---------------------------------------------------------------------------
# VRS polynomial unit tests
# ---------------------------------------------------------------------------

class TestVRSPolynomial:
    """Boundary-value checks on the Leishman VRS polynomial."""

    def test_hover_boundary(self):
        """At lambda2=0 (hover), lambda1/V_h == 1.0 exactly."""
        assert vrs_lambda1(0.0) == pytest.approx(1.0, abs=1e-10)

    def test_wbs_boundary(self):
        """At lambda2=2 (WBS entry), lambda1/V_h returns to ~1.0 (+/-5%)."""
        val = vrs_lambda1(2.0)
        assert 0.95 <= val <= 1.10, f"vrs_lambda1(2.0) = {val:.4f}, expected near 1.0"

    def test_vrs_peak_above_hover(self):
        """In mid-VRS (lambda2 ~1), induced velocity exceeds hover value."""
        assert vrs_lambda1(1.0) > 1.2, (
            f"vrs_lambda1(1.0) = {vrs_lambda1(1.0):.4f}, expected > 1.2"
        )


# ---------------------------------------------------------------------------
# Scenario 1 — Hover convergence
# ---------------------------------------------------------------------------

class TestPittPetersHover:
    """Pitt-Peters at steady hover should agree with Level-1 BEM within 5%.

    Both models use the same blade geometry, so any gap is a bug.
    """

    @pytest.mark.parametrize("theta_deg,rpm", [
        (8.86, 1200),   # Run 5, CT_nom=0.004
        (5.13, 1600),   # Run 9, CT_nom=0.002
    ])
    def test_hover_ct_matches_bem(self, pp_model, bem_model, theta_deg, rpm):
        """Pitt-Peters steady hover CT is within 5% of Level-1 BEM CT."""
        ct_bem = _bem_ct(bem_model, theta_deg, rpm, 0.0)
        ct_pp, _, _lam0 = _euler_to_steady(pp_model, theta_deg, rpm, 0.0)
        err = abs(ct_pp - ct_bem) / ct_bem
        assert err < 0.05, (
            f"theta={theta_deg} rpm={rpm}: "
            f"CT_PP={ct_pp:.5f}, CT_BEM={ct_bem:.5f}, err={err:.1%}"
        )


# ---------------------------------------------------------------------------
# Scenario 2 — VRS no CT blow-up
# ---------------------------------------------------------------------------

class TestPittPetersVRS:
    """Pitt-Peters must not produce the Level-1 BEM VRS CT explosion.

    Level-1 BEM: CT jumps to ~2.4x nominal at V/OR = 0.02 (lambda2 ~ 0.4).
    Pitt-Peters + VRS correction: CT stays below 1.6x at lambda2 ~ 0.4.
    """

    def test_shallow_vrs_ct_not_blown_up(self, pp_model, defn):
        """CT at shallow VRS (lambda2 ~ 0.4) is < 1.6x nominal."""
        theta_deg, rpm = 8.86, 1200
        CT_nom = 0.004
        omega = rpm * math.pi / 30.0
        R = defn.blade.radius_m
        A = math.pi * R**2

        # Find hover equilibrium first, then apply descent
        hover_ct_bem = _bem_ct(BEMModel(defn=defn), theta_deg, rpm, 0.0)
        _ct, _, lam0_hover = _euler_to_steady(pp_model, theta_deg, rpm, 0.0)

        v_or_target = 0.020  # V/OR = 0.020, lambda2 ~ 0.4 in Level-1 blow-up zone
        v_descent = v_or_target * omega * R
        ct_pp, _, _l = _euler_to_steady(pp_model, theta_deg, rpm, -v_descent,
                                        lam0_init=lam0_hover)
        ratio = ct_pp / CT_nom
        assert ratio < 1.6, (
            f"Shallow VRS CT = {ct_pp:.5f} ({ratio:.2f}x nominal); "
            f"Pitt-Peters should suppress VRS blow-up to < 1.6x. "
            f"Level-1 BEM produces 2.4x at the same operating point."
        )


# ---------------------------------------------------------------------------
# Scenario 3 — WBS autorotation (CQ < 0 in deep descent)
# ---------------------------------------------------------------------------

class TestPittPetersWBS:
    """In deep WBS the rotor must be in autorotation mode (CQ < 0)."""

    def test_deep_wbs_cq_negative(self, pp_model, defn):
        """CQ < 0 in deep WBS (lambda2 ~ 2.5) with hover collective."""
        theta_deg, rpm = 8.86, 1200
        omega = rpm * math.pi / 30.0
        R = defn.blade.radius_m

        # V/OR = 0.200 gives lambda2 ~ 2.9 in the scan
        v_or = 0.200
        v_descent = v_or * omega * R
        _, cq, _l = _euler_to_steady(pp_model, theta_deg, rpm, -v_descent)
        assert cq < 0, (
            f"Deep WBS CQ = {cq:.6f} should be negative (autorotation). "
            f"Check the total-inflow computation in PittPetersModel.compute_forces "
            f"(blade sees λ_total = λ_0 + v_climb/ΩR)."
        )

    def test_inflow_lag_lam0_evolves_gradually(self, pp_model, defn):
        """Pitt-Peters lam0 responds with a finite time constant, not instantly.

        At hover equilibrium dλ_0/dt ≈ 0.  After a step collective increase,
        dλ_0/dt must be nonzero (ODE is driving lam0 toward the new target),
        and after one time constant τ₀ = 8R/(3π V_T) ≈ 0.14 s the state must
        have moved but not yet reached the new steady state.  This distinguishes
        Pitt-Peters from a quasi-static model where lam0 would jump instantly.
        """
        theta_deg, rpm = 8.86, 1200
        omega = rpm * math.pi / 30.0
        R = defn.blade.radius_m
        A = math.pi * R**2

        # Find hover equilibrium lam0
        lam0_eq = 0.0
        inp_hover = RotorInputs(
            collective_rad=math.radians(theta_deg),
            tilt_lon=0.0, tilt_lat=0.0, R_hub=np.eye(3),
            v_hub_world=np.zeros(3), wind_world=np.zeros(3), t=0.0,
        )
        for _ in range(5000):
            state = PittPetersRotorState(lambda_0=lam0_eq, omega_rad_s=omega)
            res, drv = pp_model.compute_forces(inp_hover, state)
            lam0_eq += drv.lambda_0 * 0.001
        lam0_hover = lam0_eq

        # At hover equilibrium: dλ_0/dt should be near zero
        state_hover = PittPetersRotorState(lambda_0=lam0_hover, omega_rad_s=omega)
        _, drv_hover = pp_model.compute_forces(inp_hover, state_hover)
        assert abs(drv_hover.lambda_0) < 0.01, (
            f"At hover equilibrium dλ_0/dt = {drv_hover.lambda_0:.4f}/s, expected ≈0"
        )

        # Step up collective by 2°; new target lam0 is higher
        inp_step = RotorInputs(
            collective_rad=math.radians(theta_deg + 2.0),
            tilt_lon=0.0, tilt_lat=0.0, R_hub=np.eye(3),
            v_hub_world=np.zeros(3), wind_world=np.zeros(3), t=0.0,
        )
        # Immediately after step: lam0 is still at hover value, ODE is active
        state_step = PittPetersRotorState(lambda_0=lam0_hover, omega_rad_s=omega)
        _, drv_step = pp_model.compute_forces(inp_step, state_step)
        assert drv_step.lambda_0 > 0.1, (
            f"dλ_0/dt = {drv_step.lambda_0:.4f}/s after +2° collective step; "
            f"should be positive and large (driving lam0 upward)"
        )

        # After 20 ms (short relative to tau_0 ~ 0.14 s) lam0 should have moved
        # noticeably but not converged — demonstrates the inflow lag is present.
        dt = 0.001
        n_short = 20  # 20 ms
        lam0 = lam0_hover
        for _ in range(n_short):
            state = PittPetersRotorState(lambda_0=lam0, omega_rad_s=omega)
            _, drv = pp_model.compute_forces(inp_step, state)
            lam0 += drv.lambda_0 * dt
        lam0_at_20ms = lam0

        _, _cq, lam0_ss_new = _euler_to_steady(pp_model, theta_deg + 2.0, rpm, 0.0)

        # After 20 ms: lam0 has started moving (lag is not zero)
        assert lam0_at_20ms > lam0_hover, "lam0 did not increase after collective step"
        # After 20 ms: lam0 has NOT yet reached new steady state (lag is nonzero)
        frac = (lam0_at_20ms - lam0_hover) / (lam0_ss_new - lam0_hover)
        assert frac < 0.90, (
            f"lam0 moved {frac:.0%} of the way to steady state in only 20 ms; "
            f"expected lag to slow convergence (Pitt-Peters tau_0 ~ 0.14 s). "
            f"lam0: {lam0_hover:.5f} → {lam0_at_20ms:.5f} → {lam0_ss_new:.5f}"
        )


# Note: a wind-axis L-matrix rotation was implemented and reverted (it
# produced rotational covariance for `µ_y ≠ 0` but destabilised the tethered-
# rotor envelope via λ_c → BEM → C_L_hub → λ_s feedback).  See pitt_peters.py.

