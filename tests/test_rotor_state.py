import numpy as np
import pytest

from aero.rotor_state import PittPetersRotorState, QuasiStaticRotorState, RotorState


class TestQuasiStaticRotorState:
    def test_n_states(self):
        assert QuasiStaticRotorState().n_states == 2

    def test_default_fields_zero(self):
        s = QuasiStaticRotorState()
        assert s.omega_rad_s == 0.0
        assert s.spin_angle_rad == 0.0

    def test_to_array_shape_and_dtype(self):
        arr = QuasiStaticRotorState().to_array()
        assert arr.shape == (2,)
        assert arr.dtype == np.float64

    def test_to_array_values(self):
        s = QuasiStaticRotorState(omega_rad_s=10.0, spin_angle_rad=1.5)
        np.testing.assert_array_equal(s.to_array(), [10.0, 1.5])

    def test_roundtrip(self):
        s = QuasiStaticRotorState(omega_rad_s=10.0, spin_angle_rad=1.5)
        r = s.from_array(s.to_array())
        assert r.omega_rad_s == pytest.approx(10.0)
        assert r.spin_angle_rad == pytest.approx(1.5)

    def test_from_array_wrong_shape_raises(self):
        with pytest.raises(ValueError):
            QuasiStaticRotorState().from_array(np.array([1.0]))

    def test_from_array_returns_new_instance(self):
        s = QuasiStaticRotorState(omega_rad_s=5.0, spin_angle_rad=0.5)
        r = s.from_array(np.array([20.0, 2.0]))
        assert r.omega_rad_s == pytest.approx(20.0)
        assert s.omega_rad_s == pytest.approx(5.0)  # original unchanged

    def test_is_rotor_state(self):
        assert isinstance(QuasiStaticRotorState(), RotorState)


class TestPittPetersRotorState:
    def test_n_states(self):
        assert PittPetersRotorState().n_states == 5

    def test_default_fields_zero(self):
        s = PittPetersRotorState()
        assert s.lambda_0 == 0.0
        assert s.lambda_c == 0.0
        assert s.lambda_s == 0.0
        assert s.omega_rad_s == 0.0
        assert s.spin_angle_rad == 0.0

    def test_to_array_shape_and_dtype(self):
        arr = PittPetersRotorState().to_array()
        assert arr.shape == (5,)
        assert arr.dtype == np.float64

    def test_to_array_values(self):
        s = PittPetersRotorState(
            lambda_0=0.1, lambda_c=0.2, lambda_s=0.3,
            omega_rad_s=50.0, spin_angle_rad=1.0,
        )
        np.testing.assert_array_equal(s.to_array(), [0.1, 0.2, 0.3, 50.0, 1.0])

    def test_inflow_states_at_indices_0_to_2(self):
        s = PittPetersRotorState(lambda_0=0.1, lambda_c=0.2, lambda_s=0.3)
        arr = s.to_array()
        assert arr[0] == pytest.approx(0.1)
        assert arr[1] == pytest.approx(0.2)
        assert arr[2] == pytest.approx(0.3)

    def test_mechanical_states_at_indices_3_and_4(self):
        s = PittPetersRotorState(omega_rad_s=50.0, spin_angle_rad=1.0)
        arr = s.to_array()
        assert arr[3] == pytest.approx(50.0)
        assert arr[4] == pytest.approx(1.0)

    def test_roundtrip(self):
        s = PittPetersRotorState(
            lambda_0=0.1, lambda_c=-0.05, lambda_s=0.02,
            omega_rad_s=47.0, spin_angle_rad=3.1,
        )
        r = s.from_array(s.to_array())
        assert r.lambda_0 == pytest.approx(0.1)
        assert r.lambda_c == pytest.approx(-0.05)
        assert r.lambda_s == pytest.approx(0.02)
        assert r.omega_rad_s == pytest.approx(47.0)
        assert r.spin_angle_rad == pytest.approx(3.1)

    def test_from_array_wrong_shape_raises(self):
        with pytest.raises(ValueError):
            PittPetersRotorState().from_array(np.array([0.1, 0.2, 0.3]))

    def test_from_array_returns_new_instance(self):
        s = PittPetersRotorState(lambda_0=0.1, omega_rad_s=10.0)
        r = s.from_array(np.array([0.5, 0.0, 0.0, 20.0, 0.0]))
        assert r.lambda_0 == pytest.approx(0.5)
        assert r.omega_rad_s == pytest.approx(20.0)
        assert s.lambda_0 == pytest.approx(0.1)  # original unchanged
        assert s.omega_rad_s == pytest.approx(10.0)

    def test_is_rotor_state(self):
        assert isinstance(PittPetersRotorState(), RotorState)
