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
