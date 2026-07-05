import math

import numpy as np

from dynbem import AeroResult
from dynbem.rotor_definition import (
    AutorotationProperties,
    BladeGeometry,
    LinearPolarParameters,
    RotorDefinition,
)
from dynbem.rotor_state import OyeRotorState, PittPetersRotorState, QuasiStaticRotorState
from tests.helpers import hover_inputs, make_bem, make_oye, make_pitt_peters


def _test_defn() -> RotorDefinition:
    blade = BladeGeometry(
        n_blades=2,
        radius_m=1.143,
        root_cutout_m=0.1,
        chord_m=0.1905,
        twist_deg=0.0,
        n_elements=20,
    )
    airfoil = LinearPolarParameters(
        Re_design=1_000_000,
        CL0=0.0,
        CL_alpha_per_rad=2 * math.pi,
        CD0=0.008,
        alpha_stall_deg=15.0,
    )
    return RotorDefinition(
        blade=blade,
        airfoil=airfoil,
        autorotation=AutorotationProperties(I_ode_kgm2=1.0),
        name="step-api-test-rotor",
    )


def test_step_api_available_and_returns_types_for_all_models():
    defn = _test_defn()
    models = [
        (make_bem(defn), QuasiStaticRotorState),
        (make_pitt_peters(defn), PittPetersRotorState),
        (make_oye(defn), OyeRotorState),
    ]

    inp = hover_inputs(collective_deg=8.0, omega_rad_s=1250 * math.pi / 30.0)
    dt = 0.01

    for model, state_type in models:
        state = model.initial_rotor_state()
        assert hasattr(model, "step")

        result, next_state = model.step(inp, state, dt)

        assert isinstance(result, AeroResult)
        assert isinstance(next_state, state_type)


def test_step_matches_compute_forces_for_quasi_static_state_shape():
    defn = _test_defn()
    model = make_bem(defn)
    inp = hover_inputs(collective_deg=8.0, omega_rad_s=1250 * math.pi / 30.0)
    state = model.initial_rotor_state()

    result_c, deriv = model.compute_forces(inp, state)
    result_s, next_state = model.step(inp, state, dt=0.02)

    assert isinstance(result_c, AeroResult)
    assert isinstance(result_s, AeroResult)
    assert isinstance(deriv, QuasiStaticRotorState)
    assert isinstance(next_state, QuasiStaticRotorState)
    np.testing.assert_array_equal(deriv.to_array(), np.array([], dtype=float))
    np.testing.assert_array_equal(next_state.to_array(), np.array([], dtype=float))
