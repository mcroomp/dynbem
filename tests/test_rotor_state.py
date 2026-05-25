import numpy as np
import pytest

from dynbem.rotor_state import PittPetersRotorState, QuasiStaticRotorState, RotorState


class TestQuasiStaticRotorState:
    def test_to_array_shape_and_dtype(self):
        arr = QuasiStaticRotorState().to_array()
        assert arr.shape == (0,)
        assert arr.dtype == np.float64

    def test_roundtrip(self):
        s = QuasiStaticRotorState()
        r = s.from_array(s.to_array())
        assert isinstance(r, QuasiStaticRotorState)

    def test_from_array_wrong_shape_raises(self):
        with pytest.raises((ValueError, TypeError)):
            QuasiStaticRotorState().from_array(np.array([1.0]))

    def test_is_rotor_state(self):
        assert isinstance(QuasiStaticRotorState(), RotorState)


class TestPittPetersRotorState:
    def test_default_fields_zero(self):
        s = PittPetersRotorState(0.0, 0.0, 0.0)
        assert s.lambda_0 == 0.0
        assert s.lambda_c == 0.0
        assert s.lambda_s == 0.0

    def test_to_array_shape_and_dtype(self):
        arr = PittPetersRotorState(0.0, 0.0, 0.0).to_array()
        assert arr.shape == (3,)
        assert arr.dtype == np.float64

    def test_to_array_values(self):
        s = PittPetersRotorState(
            lambda_0=0.1, lambda_c=0.2, lambda_s=0.3,
        )
        np.testing.assert_array_equal(s.to_array(), [0.1, 0.2, 0.3])

    def test_inflow_states_at_indices_0_to_2(self):
        s = PittPetersRotorState(lambda_0=0.1, lambda_c=0.2, lambda_s=0.3)
        arr = s.to_array()
        assert arr[0] == pytest.approx(0.1)
        assert arr[1] == pytest.approx(0.2)
        assert arr[2] == pytest.approx(0.3)

    def test_roundtrip(self):
        s = PittPetersRotorState(
            lambda_0=0.1, lambda_c=-0.05, lambda_s=0.02,
        )
        r = s.from_array(s.to_array())
        assert r.lambda_0 == pytest.approx(0.1)
        assert r.lambda_c == pytest.approx(-0.05)
        assert r.lambda_s == pytest.approx(0.02)

    def test_from_array_wrong_shape_raises(self):
        with pytest.raises((ValueError, TypeError)):
            PittPetersRotorState(0.0, 0.0, 0.0).from_array(np.array([0.1, 0.2]))

    def test_from_array_returns_new_instance(self):
        s = PittPetersRotorState(0.1, 0.0, 0.0)
        r = s.from_array(np.array([0.5, 0.0, 0.0]))
        assert r.lambda_0 == pytest.approx(0.5)
        assert s.lambda_0 == pytest.approx(0.1)  # original unchanged

    def test_is_rotor_state(self):
        assert isinstance(PittPetersRotorState(0.0, 0.0, 0.0), RotorState)
