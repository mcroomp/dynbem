"""Python API smoke tests for trim wrappers.

Detailed trim-solver behavior tests live in Rust integration tests:
`dynbem_rs/tests/trim_solver.rs`.
This module verifies Python binding plumbing only.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from dynbem import RotorInputs, create_aero, relax_inflow, solve_trim_cyclic
from dynbem.rotor_definition import load as load_rotor


_ROTOR_YAML = str(
    Path(__file__).parent.parent / "rotors" / "beaupoil_2026" / "rotor.yaml"
)
_OMEGA = 28.0
_COLLECTIVE = math.radians(-9.0)


@pytest.fixture(scope="module")
def defn():
    return load_rotor(_ROTOR_YAML)


def _inputs() -> RotorInputs:
    return RotorInputs(
        collective_rad=_COLLECTIVE,
        tilt_lon=0.0,
        tilt_lat=0.0,
        R_hub=np.eye(3),
        v_hub_world=np.zeros(3),
        wind_world=np.array([0.0, 10.0, 0.0]),
        rho_kg_m3=1.225,
        omega_rad_s=_OMEGA,
        t=0.0,
    )


@pytest.mark.parametrize("model_name", ["pitt_peters", "oye"])
def test_solve_trim_cyclic_api_smoke(defn, model_name):
    aero = create_aero(defn, model=model_name)
    state = aero.initial_rotor_state()

    result = solve_trim_cyclic(
        aero,
        state,
        _inputs(),
        tolerance_Nm=1.7 if model_name == "pitt_peters" else 0.05,
    )

    assert isinstance(result.converged, bool)
    assert result.iterations >= 0
    assert math.isfinite(result.tilt_lon)
    assert math.isfinite(result.tilt_lat)
    assert hasattr(result.final_state, "to_array")


@pytest.mark.parametrize("model_name", ["pitt_peters", "oye"])
def test_relax_inflow_api_smoke(defn, model_name):
    aero = create_aero(defn, model=model_name)
    s0 = aero.initial_rotor_state()
    s1 = relax_inflow(aero, s0, _inputs(), n_steps=100, dt=0.005)

    a0 = np.asarray(s0.to_array(), dtype=float)
    a1 = np.asarray(s1.to_array(), dtype=float)
    assert a0.shape == a1.shape
    assert np.all(np.isfinite(a1))


def test_trim_clips_to_bounds_api(defn):
    aero = create_aero(defn, model="oye")
    state = aero.initial_rotor_state()
    tight = math.radians(1.0)

    result = solve_trim_cyclic(
        aero,
        state,
        _inputs(),
        tilt_min=-tight,
        tilt_max=tight,
        tolerance_Nm=0.01,
        max_iterations=20,
    )

    assert -tight <= result.tilt_lon <= tight
    assert -tight <= result.tilt_lat <= tight
