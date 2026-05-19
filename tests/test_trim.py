"""Unit tests for aero.solve_trim_cyclic.

The trim solver finds the (tilt_lon, tilt_lat) cyclic that nulls the
in-plane hub moments at a given operating point.  These tests verify:

  * In axisymmetric hover the trim is the zero cyclic (no asymmetry to
    cancel).
  * In forward flight the trim residual is below tolerance.
  * The reported residual matches a direct ``compute_forces`` evaluation
    at the trim cyclic.
  * The solver handles both Pitt-Peters and Øye models without
    model-specific knobs.
  * Trim cyclic clips inside the requested ``tilt_min/tilt_max``
    bounds.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from aero import RotorInputs, create_aero, relax_inflow, solve_trim_cyclic
from aero.rotor_definition import load as load_rotor


_ROTOR_YAML = str(
    Path(__file__).parent.parent / "rotors" / "beaupoil_2026" / "rotor.yaml"
)
_OMEGA      = 28.0           # rad/s
_COLLECTIVE = math.radians(-9.0)
_TOL_NM     = 0.05            # newtonmeter target tolerance


@pytest.fixture(scope="module")
def defn():
    return load_rotor(_ROTOR_YAML)


def _level_R():
    """Level FRD hub: body_z points DOWN through the disk."""
    return np.eye(3)


def _moments_at(aero, state, *, collective, tilt_lon, tilt_lat,
                R_hub, v_hub_world, wind_world):
    inputs = RotorInputs(
        collective_rad=collective,
        tilt_lon=tilt_lon, tilt_lat=tilt_lat,
        R_hub=R_hub, v_hub_world=v_hub_world, wind_world=wind_world, t=0.0,
    )
    res, _ = aero.compute_forces(inputs, state)
    M_hub = R_hub.T @ res.M_orbital
    return float(M_hub[0]), float(M_hub[1])


@pytest.mark.parametrize("model_name", ["pitt_peters_jit", "oye"])
class TestTrimSolver:

    def test_hover_trim_is_near_zero(self, defn, model_name):
        """Axisymmetric hover (no wind) ⇒ no asymmetry ⇒ trim cyclic ≈ 0."""
        aero  = create_aero(defn, model=model_name)
        state = aero.initial_rotor_state()
        state.omega_rad_s = _OMEGA

        result = solve_trim_cyclic(
            aero, state,
            collective_rad=_COLLECTIVE,
            R_hub=_level_R(),
            v_hub_world=np.zeros(3), wind_world=np.zeros(3),
            tolerance_Nm=_TOL_NM,
        )
        assert result.converged, (
            f"{model_name}: hover trim did not converge in "
            f"{result.iterations} iters (Mx={result.Mx_residual:.4f}, "
            f"My={result.My_residual:.4f})"
        )
        assert abs(result.tilt_lon) < math.radians(0.5)
        assert abs(result.tilt_lat) < math.radians(0.5)

    def test_forward_flight_trim_residual_below_tolerance(self, defn, model_name):
        """10 m/s East wind across level disk ⇒ trim residual < tolerance."""
        aero  = create_aero(defn, model=model_name)
        state = aero.initial_rotor_state()
        state.omega_rad_s = _OMEGA

        result = solve_trim_cyclic(
            aero, state,
            collective_rad=_COLLECTIVE,
            R_hub=_level_R(),
            v_hub_world=np.zeros(3),
            wind_world=np.array([0.0, 10.0, 0.0]),
            tolerance_Nm=_TOL_NM,
        )
        assert result.converged, (
            f"{model_name}: trim did not converge: "
            f"Mx={result.Mx_residual:.4f}, My={result.My_residual:.4f}"
        )
        assert abs(result.Mx_residual) < _TOL_NM
        assert abs(result.My_residual) < _TOL_NM

    def test_trim_residual_matches_direct_evaluation(self, defn, model_name):
        """The residual reported by the solver equals the moment at the
        returned trim cyclic, sampled directly from compute_forces using
        the returned ``final_state`` (which holds the relaxed inflow)."""
        aero  = create_aero(defn, model=model_name)
        state = aero.initial_rotor_state()
        state.omega_rad_s = _OMEGA
        wind = np.array([0.0, 10.0, 0.0])

        result = solve_trim_cyclic(
            aero, state,
            collective_rad=_COLLECTIVE,
            R_hub=_level_R(),
            v_hub_world=np.zeros(3), wind_world=wind,
            tolerance_Nm=_TOL_NM,
        )
        Mx, My = _moments_at(
            aero, result.final_state,
            collective=_COLLECTIVE,
            tilt_lon=result.tilt_lon, tilt_lat=result.tilt_lat,
            R_hub=_level_R(), v_hub_world=np.zeros(3), wind_world=wind,
        )
        assert abs(Mx - result.Mx_residual) < 1e-6
        assert abs(My - result.My_residual) < 1e-6

    def test_trim_to_nonzero_target_moment(self, defn, model_name):
        """target_moment = (0, M_target) should produce a trim where
        My_hub ≈ M_target (within tolerance)."""
        aero  = create_aero(defn, model=model_name)
        state = aero.initial_rotor_state()
        state.omega_rad_s = _OMEGA
        wind = np.array([0.0, 10.0, 0.0])
        M_target = 5.0   # N·m on body-Y, small enough to be feasible

        result = solve_trim_cyclic(
            aero, state,
            collective_rad=_COLLECTIVE,
            R_hub=_level_R(),
            v_hub_world=np.zeros(3), wind_world=wind,
            target_moment=(0.0, M_target),
            tolerance_Nm=_TOL_NM,
        )
        assert result.converged
        Mx, My = _moments_at(
            aero, result.final_state,
            collective=_COLLECTIVE,
            tilt_lon=result.tilt_lon, tilt_lat=result.tilt_lat,
            R_hub=_level_R(), v_hub_world=np.zeros(3), wind_world=wind,
        )
        assert abs(Mx)            < _TOL_NM, f"Mx={Mx:.4f} should be near 0"
        assert abs(My - M_target) < _TOL_NM, f"My={My:.4f} should be near {M_target}"

    def test_relax_inflow_settles_to_steady_state(self, defn, model_name):
        """``relax_inflow`` reaches a fixed point on the inflow states.

        We check inflow only — ω is held fixed by ``fix_omega`` and the
        spin angle ψ accumulates monotonically (rotor keeps spinning),
        so excluding both is the right convergence diagnostic.
        """
        aero = create_aero(defn, model=model_name)
        s0   = aero.initial_rotor_state()
        s0.omega_rad_s = _OMEGA

        kw = dict(
            collective_rad=_COLLECTIVE, tilt_lon=0.0, tilt_lat=0.0,
            R_hub=_level_R(), v_hub_world=np.zeros(3),
            wind_world=np.array([0.0, 10.0, 0.0]),
            n_steps=500, dt=0.005, fix_omega=True,
        )
        s1 = relax_inflow(aero, s0, **kw)
        s2 = relax_inflow(aero, s1, **kw)
        # arr[-2] = ω (held fixed), arr[-1] = ψ (monotonic).
        # Inflow states are at arr[:-2].
        inflow_delta = float(np.linalg.norm(s2.to_array()[:-2] - s1.to_array()[:-2]))
        assert inflow_delta < 1e-4, (
            f"{model_name}: inflow not settled (Δ_inflow={inflow_delta:.4e})"
        )

    def test_solver_cancels_baseline_disturbance(self, defn, model_name):
        """At zero cyclic the wind-driven baseline moment is large; at the
        solver's trim cyclic the same moment is small."""
        aero  = create_aero(defn, model=model_name)
        state = aero.initial_rotor_state()
        state.omega_rad_s = _OMEGA
        wind = np.array([0.0, 10.0, 0.0])

        # Settle inflow at zero cyclic, then read baseline moment.
        for _ in range(200):
            inputs = RotorInputs(
                collective_rad=_COLLECTIVE, tilt_lon=0.0, tilt_lat=0.0,
                R_hub=_level_R(), v_hub_world=np.zeros(3), wind_world=wind, t=0.0,
            )
            _, deriv = aero.compute_forces(inputs, state)
            arr = state.to_array() + 0.005 * deriv.to_array()
            arr[-2] = _OMEGA   # hold ω fixed
            state = state.from_array(arr)
        Mx0, My0 = _moments_at(
            aero, state,
            collective=_COLLECTIVE, tilt_lon=0.0, tilt_lat=0.0,
            R_hub=_level_R(), v_hub_world=np.zeros(3), wind_world=wind,
        )
        baseline_mag = math.hypot(Mx0, My0)

        # Now run the solver from this state.
        result = solve_trim_cyclic(
            aero, state,
            collective_rad=_COLLECTIVE,
            R_hub=_level_R(),
            v_hub_world=np.zeros(3), wind_world=wind,
            tolerance_Nm=_TOL_NM,
        )
        trim_mag = math.hypot(result.Mx_residual, result.My_residual)

        # Sanity: baseline disturbance is large enough that the test is meaningful.
        assert baseline_mag > 10.0, (
            f"{model_name}: baseline disturbance too small to test "
            f"({baseline_mag:.2f} N*m) — pick a stronger wind/operating point"
        )
        # The solver should knock the residual down by ≥ 100×.
        assert trim_mag < baseline_mag / 100.0, (
            f"{model_name}: solver did not cancel disturbance "
            f"(baseline {baseline_mag:.2f}, trim {trim_mag:.4f} N*m)"
        )


def test_trim_clips_to_bounds(defn):
    """Trim cyclic is clipped to the requested [tilt_min, tilt_max]."""
    aero  = create_aero(defn, model="oye")
    state = aero.initial_rotor_state()
    state.omega_rad_s = _OMEGA

    tight = math.radians(1.0)
    result = solve_trim_cyclic(
        aero, state,
        collective_rad=_COLLECTIVE,
        R_hub=np.eye(3),
        v_hub_world=np.zeros(3),
        wind_world=np.array([0.0, 10.0, 0.0]),
        tilt_min=-tight, tilt_max=tight,
        tolerance_Nm=0.01,
        max_iterations=20,
    )
    assert -tight <= result.tilt_lon <= tight
    assert -tight <= result.tilt_lat <= tight
