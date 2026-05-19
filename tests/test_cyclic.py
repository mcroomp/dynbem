"""Cyclic-control wiring: helicopter-standard signs on all three models.

Convention (see CLAUDE.md "Rotor rotation direction"):
  CCW from above, ψ=0 at +X (hub-frame nose).
  tilt_lon > 0  ⇒  nose-down (forward stick)  ⇒  M_y < 0
  tilt_lat > 0  ⇒  roll right                 ⇒  M_x > 0

These tests run at hover so any moment is purely from cyclic — no Glauert
wake-skew contamination.  The level rotor (R_hub = I) makes hub-frame and
world-frame moments identical.
"""

import math
from pathlib import Path

import numpy as np
import pytest

from aero import RotorInputs, create_aero
from aero.cyclic import cyclic_coeffs
from aero.rotor_definition import load as load_rotor
from aero.rotor_state import PittPetersRotorState, QuasiStaticRotorState


_ROTOR_YAML = str(
    Path(__file__).parent.parent / "rotors" / "castles_gray_6ft" / "rotor.yaml"
)
_RPM = 1200.0
_OMEGA = _RPM * math.pi / 30.0
_COLLECTIVE = math.radians(8.0)


# ---------------------------------------------------------------------------
# Helper builders
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def defn():
    return load_rotor(_ROTOR_YAML)


def _make_inputs(tilt_lon: float, tilt_lat: float) -> RotorInputs:
    return RotorInputs(
        collective_rad=_COLLECTIVE,
        tilt_lon=tilt_lon, tilt_lat=tilt_lat,
        R_hub=np.eye(3),
        v_hub_world=np.zeros(3),
        wind_world=np.zeros(3),
        t=0.0,
    )


def _initial_state(model):
    s = model.initial_rotor_state()
    # Spin the rotor up by replacing omega; other fields keep defaults.
    if isinstance(s, PittPetersRotorState):
        return PittPetersRotorState(omega_rad_s=_OMEGA)
    return QuasiStaticRotorState(omega_rad_s=_OMEGA)


# ---------------------------------------------------------------------------
# cyclic_coeffs helper
# ---------------------------------------------------------------------------

def test_cyclic_coeffs_no_control_defaults_to_direct_amplitudes():
    # gain=1, phase=0 ⇒ θ_1c = -tilt_lon, θ_1s = +tilt_lat.
    c, s = cyclic_coeffs(tilt_lon=0.05, tilt_lat=0.0, control=None)
    assert c == pytest.approx(-0.05)
    assert s == pytest.approx(0.0)

    c, s = cyclic_coeffs(tilt_lon=0.0, tilt_lat=0.03, control=None)
    assert c == pytest.approx(0.0)
    assert s == pytest.approx(0.03)


def test_cyclic_coeffs_phase_rotates_cos_to_sin():
    # phi = 90° advances the cos coefficient to sin.
    from aero.rotor_definition import ControlProperties

    ctrl = ControlProperties(
        swashplate_pitch_gain_rad=1.0,
        swashplate_phase_deg=90.0,
    )
    c, s = cyclic_coeffs(tilt_lon=0.05, tilt_lat=0.0, control=ctrl)
    # -tilt_lon·cos(90°) - tilt_lat·sin(90°) = 0
    assert c == pytest.approx(0.0, abs=1e-12)
    # -tilt_lon·sin(90°) + tilt_lat·cos(90°) = -0.05
    assert s == pytest.approx(-0.05)


# ---------------------------------------------------------------------------
# Model-level wiring — helicopter-standard signs
# ---------------------------------------------------------------------------

@pytest.fixture(params=["bem", "pitt_peters", "pitt_peters_jit"])
def model(defn, request):
    return create_aero(defn, model=request.param)


def test_zero_cyclic_zero_in_plane_moment(model):
    """Pure hover with no cyclic: M_orbital ≈ 0 by symmetry."""
    res, _ = model.compute_forces(_make_inputs(0.0, 0.0), _initial_state(model))
    # Thrust must be substantial — sanity that we're in a real operating point.
    assert -res.F_world[2] > 50.0
    assert np.linalg.norm(res.M_orbital) < 1e-3


def test_tilt_lon_positive_gives_nose_down_moment(model):
    """tilt_lon > 0 ⇒ M_orbital_y < 0 (nose-down)."""
    res, _ = model.compute_forces(
        _make_inputs(tilt_lon=math.radians(2.0), tilt_lat=0.0),
        _initial_state(model),
    )
    # M_y should be clearly negative; the orthogonal M_x should be near zero.
    assert res.M_orbital[1] < -1.0
    assert abs(res.M_orbital[0]) < 0.1 * abs(res.M_orbital[1])


def test_tilt_lat_positive_gives_roll_right_moment(model):
    """tilt_lat > 0 ⇒ M_orbital_x > 0 (roll right)."""
    res, _ = model.compute_forces(
        _make_inputs(tilt_lon=0.0, tilt_lat=math.radians(2.0)),
        _initial_state(model),
    )
    assert res.M_orbital[0] > 1.0
    assert abs(res.M_orbital[1]) < 0.1 * abs(res.M_orbital[0])


def test_cyclic_sign_symmetry(model):
    """Reversing the cyclic input reverses the moment."""
    pos, _ = model.compute_forces(
        _make_inputs(tilt_lon=math.radians(1.5), tilt_lat=0.0),
        _initial_state(model),
    )
    neg, _ = model.compute_forces(
        _make_inputs(tilt_lon=-math.radians(1.5), tilt_lat=0.0),
        _initial_state(model),
    )
    np.testing.assert_allclose(pos.M_orbital, -neg.M_orbital, atol=1e-6, rtol=1e-3)


def test_thrust_roughly_invariant_to_small_cyclic(model):
    """Small cyclic should mostly redistribute thrust, not change the total."""
    base, _ = model.compute_forces(_make_inputs(0.0, 0.0), _initial_state(model))
    cyc,  _ = model.compute_forces(
        _make_inputs(tilt_lon=math.radians(1.0), tilt_lat=math.radians(1.0)),
        _initial_state(model),
    )
    T_base = -base.F_world[2]
    T_cyc  = -cyc.F_world[2]
    # Thrust changes by <5% for ~1° cyclic at hover.
    assert abs(T_cyc - T_base) / T_base < 0.05


# ---------------------------------------------------------------------------
# Pitt-Peters dynamic-inflow coupling to cyclic — diagonal momentum balance.
# These run only on the Pitt-Peters models (not Level-1 BEM).
# ---------------------------------------------------------------------------

def _euler_integrate_pp(model, inputs, n_steps=8000, dt=0.0005):
    """Time-integrate a Pitt-Peters model to (approximate) steady state.

    Uses a small fixed-step Euler — long enough for cyclic inflow harmonics
    (τ_cs ≈ tens of ms at hover) to converge.  Returns the final state and
    the final compute_forces result.
    """
    from aero.rotor_state import PittPetersRotorState

    state = PittPetersRotorState(omega_rad_s=_OMEGA)
    res = None
    for _ in range(n_steps):
        res, drv = model.compute_forces(inputs, state)
        state = PittPetersRotorState(
            lambda_0=state.lambda_0 + drv.lambda_0 * dt,
            lambda_c=state.lambda_c + drv.lambda_c * dt,
            lambda_s=state.lambda_s + drv.lambda_s * dt,
            omega_rad_s=_OMEGA,  # hold omega fixed for this test
        )
    return state, res


@pytest.fixture(params=["pitt_peters", "pitt_peters_jit"])
def pp_model(defn, request):
    return create_aero(defn, model=request.param)


def test_hover_cyclic_drives_lambda_c_negative_for_nose_down(pp_model):
    """tilt_lon > 0 (nose-down) ⇒ M_y < 0 ⇒ λ_c_ss < 0 ⇒ converged λ_c < 0.

    Physically: peak thrust at ψ=π (tail) needs peak inflow at ψ=π to balance
    by momentum theory.  λ(r,ψ) = λ_0 + (r/R)·λ_c·cos(ψ); cos(π) = -1, so peak
    inflow at ψ=π requires λ_c < 0.
    """
    state, _ = _euler_integrate_pp(
        pp_model, _make_inputs(tilt_lon=math.radians(2.0), tilt_lat=0.0)
    )
    assert state.lambda_c < -1e-4
    assert abs(state.lambda_s) < 0.1 * abs(state.lambda_c)


def test_hover_cyclic_drives_lambda_s_positive_for_roll_right(pp_model):
    """tilt_lat > 0 (roll right) ⇒ M_x > 0 ⇒ λ_s_ss > 0 ⇒ converged λ_s > 0."""
    state, _ = _euler_integrate_pp(
        pp_model, _make_inputs(tilt_lon=0.0, tilt_lat=math.radians(2.0))
    )
    assert state.lambda_s > 1e-4
    assert abs(state.lambda_c) < 0.1 * abs(state.lambda_s)


def test_hover_no_cyclic_keeps_lambda_cs_zero(pp_model):
    """Pure hover with zero cyclic: axisymmetric — λ_c, λ_s stay at zero."""
    state, _ = _euler_integrate_pp(
        pp_model, _make_inputs(0.0, 0.0)
    )
    assert abs(state.lambda_c) < 1e-6
    assert abs(state.lambda_s) < 1e-6
    # And λ_0 should converge to a sensible positive hover induced inflow.
    assert state.lambda_0 > 0.01


def test_cyclic_inflow_reduces_hub_moment(pp_model):
    """Steady-state inflow should partially cancel the cyclic-driven moment.

    Compare M_orbital before/after letting Pitt-Peters reach steady state with
    a cyclic input held constant: the converged inflow tilt reduces local AoA
    where pitch is increased, so the *integrated* moment is smaller than the
    transient (pre-inflow) value.
    """
    inputs = _make_inputs(tilt_lon=math.radians(2.0), tilt_lat=0.0)
    state0 = PittPetersRotorState(omega_rad_s=_OMEGA)
    res_initial, _ = pp_model.compute_forces(inputs, state0)
    state_ss, res_ss = _euler_integrate_pp(pp_model, inputs)
    # M_y should be reduced in magnitude (inflow opposing the asymmetry).
    assert 0 < abs(res_ss.M_orbital[1]) < abs(res_initial.M_orbital[1])
    # Sign preserved — still nose-down.
    assert res_ss.M_orbital[1] < 0


def test_forward_flight_glauert_emerges_from_momentum_balance(pp_model):
    """Forward flight, no cyclic: λ_c, λ_s should still develop nonzero values
    purely from the advancing/retreating velocity asymmetry — this is the
    Glauert wake-skew effect emerging from the diagonal momentum balance.
    """
    inputs = RotorInputs(
        collective_rad=_COLLECTIVE,
        tilt_lon=0.0, tilt_lat=0.0,
        R_hub=np.eye(3),
        v_hub_world=np.array([10.0, 0.0, 0.0]),  # 10 m/s in +X
        wind_world=np.zeros(3),
        t=0.0,
    )
    state, _ = _euler_integrate_pp(pp_model, inputs)
    # Some longitudinal inflow harmonic should appear from velocity asymmetry.
    # Sign: vehicle moves +X ⇒ relative wind in -X ⇒ wake skews -X ⇒
    # more inflow at -X (tail = ψ=π) ⇒ -λ_c > 0 ⇒ λ_c < 0.
    assert state.lambda_c < -1e-5
