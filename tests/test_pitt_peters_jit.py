"""Validate PittPetersModelJIT against PittPetersModel on a grid of inputs.

Both implementations should produce numerically identical outputs (within
floating-point tolerance) over the operating envelope we care about: axial
flight (hover, climb, descent, autorotation) and forward flight (mu > 0.01).
"""
from pathlib import Path

import math
import numpy as np
import pytest

import aero.rotor_definition as rotor_definition
from aero import RotorInputs
from aero.pitt_peters import PittPetersModel
from aero.pitt_peters_jit import PittPetersModelJIT
from aero.rotor_state import PittPetersRotorState


@pytest.fixture(scope="module")
def defn():
    return rotor_definition.load(
        str(Path(__file__).parent.parent / "rotors" / "beaupoil_2026" / "rotor.yaml")
    )


@pytest.fixture(scope="module")
def models(defn):
    return PittPetersModel(defn=defn), PittPetersModelJIT(defn=defn)


def _make_state(omega, lam0=0.0, lam_c=0.0, lam_s=0.0):
    return PittPetersRotorState(
        lambda_0=lam0, lambda_c=lam_c, lambda_s=lam_s,
        omega_rad_s=omega, spin_angle_rad=0.0,
    )


def _make_inputs(col, v_hub_world=(0.0, 0.0, 0.0), wind=(0.0, 0.0, -10.0),
                 R_hub=None):
    if R_hub is None:
        R_hub = np.eye(3)
    return RotorInputs(
        collective_rad=col,
        tilt_lon=0.0,
        tilt_lat=0.0,
        R_hub=R_hub,
        v_hub_world=np.array(v_hub_world, dtype=float),
        wind_world=np.array(wind, dtype=float),
        t=0.0,
    )


def _assert_match(ref_result, jit_result, ref_dstate, jit_dstate,
                  rtol=1e-9, atol=1e-9):
    np.testing.assert_allclose(jit_result.F_world, ref_result.F_world,
                               rtol=rtol, atol=atol)
    np.testing.assert_allclose(jit_result.Q_spin, ref_result.Q_spin,
                               rtol=rtol, atol=atol)
    np.testing.assert_allclose(jit_result.M_spin, ref_result.M_spin,
                               rtol=rtol, atol=atol)
    np.testing.assert_allclose(jit_dstate.to_array(), ref_dstate.to_array(),
                               rtol=rtol, atol=atol)


# ---------------------------------------------------------------------------
# Axial: hover, climb, descent, autorotation (covers VRS region)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("col_deg", [-10.0, -5.0, 0.0, 5.0, 10.0])
@pytest.mark.parametrize("omega", [20.0, 60.0, 100.0])
@pytest.mark.parametrize("v_climb", [-10.0, -3.0, 0.0, 2.0, 8.0])
@pytest.mark.parametrize("lam0", [-0.05, 0.0, 0.05])
def test_axial_match(models, col_deg, omega, v_climb, lam0):
    ref, jit = models
    state = _make_state(omega, lam0=lam0)
    inputs = _make_inputs(
        col=math.radians(col_deg),
        v_hub_world=(0.0, 0.0, 0.0),
        wind=(0.0, 0.0, v_climb),
    )
    r_ref, d_ref = ref.compute_forces(inputs, state)
    r_jit, d_jit = jit.compute_forces(inputs, state)
    _assert_match(r_ref, r_jit, d_ref, d_jit, rtol=1e-8, atol=1e-9)


# ---------------------------------------------------------------------------
# Forward flight: mu > 0.01, non-zero cyclic inflow states
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("col_deg", [-5.0, 0.0, 8.0])
@pytest.mark.parametrize("omega", [40.0, 100.0])
@pytest.mark.parametrize("v_edge", [3.0, 10.0])
@pytest.mark.parametrize("lam_c", [-0.05, 0.0, 0.05])
def test_forward_flight_match(models, col_deg, omega, v_edge, lam_c):
    ref, jit = models
    state = _make_state(omega, lam0=0.02, lam_c=lam_c, lam_s=-lam_c)
    # Edgewise wind in +x direction, plus a bit of axial
    inputs = _make_inputs(
        col=math.radians(col_deg),
        v_hub_world=(0.0, 0.0, 0.0),
        wind=(v_edge, 0.0, -5.0),
    )
    r_ref, d_ref = ref.compute_forces(inputs, state)
    r_jit, d_jit = jit.compute_forces(inputs, state)
    _assert_match(r_ref, r_jit, d_ref, d_jit, rtol=1e-7, atol=1e-9)


# ---------------------------------------------------------------------------
# Tilted hub (non-identity R_hub) — guards against frame mixups
# ---------------------------------------------------------------------------

def test_tilted_hub_match(models):
    ref, jit = models
    # 30 deg tilt about x-axis
    c, s = math.cos(math.radians(30.0)), math.sin(math.radians(30.0))
    R = np.array([[1, 0, 0], [0, c, -s], [0, s, c]], dtype=float)
    state = _make_state(omega=80.0, lam0=0.02)
    inputs = _make_inputs(col=math.radians(-3.0), wind=(0.0, 0.0, -10.0), R_hub=R)
    r_ref, d_ref = ref.compute_forces(inputs, state)
    r_jit, d_jit = jit.compute_forces(inputs, state)
    _assert_match(r_ref, r_jit, d_ref, d_jit, rtol=1e-7, atol=1e-9)


# ---------------------------------------------------------------------------
# Edge cases: zero omega, very small omega
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("omega", [0.0, 0.5])
def test_low_omega(models, omega):
    ref, jit = models
    state = _make_state(omega=omega)
    inputs = _make_inputs(col=0.0)
    r_ref, d_ref = ref.compute_forces(inputs, state)
    r_jit, d_jit = jit.compute_forces(inputs, state)
    _assert_match(r_ref, r_jit, d_ref, d_jit, rtol=1e-7, atol=1e-9)
