import math

import numpy as np
import pytest

from dynbem.polar import LinearPolar, TabulatedPolar


@pytest.fixture
def polar():
    return LinearPolar(CL0=0.0, CL_alpha_per_rad=2 * math.pi, CD0=0.01, alpha_stall_rad=math.radians(15))


class TestLinearPolar:
    def test_zero_alpha(self, polar):
        cl, cd = polar.cl_cd(0.0)
        assert cl == pytest.approx(0.0)
        assert cd == pytest.approx(0.01)

    def test_positive_alpha_attached(self, polar):
        alpha = math.radians(5)
        cl, cd = polar.cl_cd(alpha)
        assert cl == pytest.approx(2 * math.pi * alpha, rel=1e-6)
        assert cd == pytest.approx(0.01)

    def test_negative_alpha_attached(self, polar):
        alpha = math.radians(-5)
        cl, cd = polar.cl_cd(alpha)
        assert cl == pytest.approx(2 * math.pi * alpha, rel=1e-6)
        assert cd == pytest.approx(0.01)

    def test_at_stall_boundary(self, polar):
        alpha = math.radians(15)
        cl_just_below, _ = polar.cl_cd(math.radians(14.99))
        cl_at, _ = polar.cl_cd(alpha)
        # at stall the linear and stall models agree (boundary is inclusive of stall branch)
        assert cl_at == pytest.approx(cl_just_below, rel=1e-2)

    def test_above_stall_cl_capped(self, polar):
        cl_at_stall, _ = polar.cl_cd(polar.alpha_stall_rad)
        cl_further, _ = polar.cl_cd(polar.alpha_stall_rad * 2.0)
        assert cl_further == pytest.approx(cl_at_stall, rel=1e-6)

    def test_above_stall_cd_increases(self, polar):
        _, cd_attached = polar.cl_cd(math.radians(5))
        _, cd_stalled = polar.cl_cd(math.radians(20))
        assert cd_stalled > cd_attached

    def test_stall_symmetry_negative(self, polar):
        cl_pos, _ = polar.cl_cd(math.radians(20))
        cl_neg, _ = polar.cl_cd(math.radians(-20))
        assert cl_neg == pytest.approx(-cl_pos, rel=1e-6)

    def test_from_properties(self):
        from dynbem.rotor_definition import AirfoilProperties
        props = AirfoilProperties(
            Re_design=500_000,
            CL0=0.1,
            CL_alpha_per_rad=5.5,
            CD0=0.015,
            alpha_stall_deg=14.0,
            tip_loss=True,
        )
        p = LinearPolar.from_properties(props)
        assert p.CL0 == pytest.approx(0.1)
        assert p.CL_alpha_per_rad == pytest.approx(5.5)
        assert p.CD0 == pytest.approx(0.015)
        assert p.alpha_stall_rad == pytest.approx(math.radians(14.0))


class TestTabulatedPolar:
    """TabulatedPolar must implement numpy.interp's semantics exactly:
    linear between knots, clamp at endpoints, no periodic wrap. The full
    S809 cross-check lives in verification/dynbem_polar_interp_check.py
    and verification/dynbem_polar_vs_aerodyn_nrel_phase_vi.py; this is the
    minimal in-tree spec-compliance test."""

    @pytest.fixture
    def small_polar(self):
        alpha = np.array([-0.2, -0.1, 0.0, 0.1, 0.2], dtype=np.float64)
        cl    = np.array([-1.0, -0.4, 0.05, 0.55, 1.05], dtype=np.float64)
        cd    = np.array([0.05, 0.02, 0.012, 0.018, 0.04], dtype=np.float64)
        return alpha, cl, cd, TabulatedPolar(alpha_rad=alpha, cl=cl, cd=cd)

    def test_exact_at_knots(self, small_polar):
        alpha, cl, cd, p = small_polar
        for i in range(len(alpha)):
            got_cl, got_cd = p.cl_cd(float(alpha[i]))
            assert got_cl == pytest.approx(cl[i], abs=0.0, rel=1e-15)
            assert got_cd == pytest.approx(cd[i], abs=0.0, rel=1e-15)

    def test_linear_between_knots(self, small_polar):
        alpha, cl, cd, p = small_polar
        # Mid-point between knots 1 and 2 (-0.05 rad) -> linear average
        a = (alpha[1] + alpha[2]) / 2
        got_cl, got_cd = p.cl_cd(a)
        expected_cl = (cl[1] + cl[2]) / 2
        expected_cd = (cd[1] + cd[2]) / 2
        assert got_cl == pytest.approx(expected_cl, abs=1e-14)
        assert got_cd == pytest.approx(expected_cd, abs=1e-14)

    def test_clamp_below(self, small_polar):
        alpha, cl, cd, p = small_polar
        got_cl, got_cd = p.cl_cd(alpha[0] - 10.0)
        assert got_cl == pytest.approx(cl[0])
        assert got_cd == pytest.approx(cd[0])

    def test_clamp_above(self, small_polar):
        alpha, cl, cd, p = small_polar
        got_cl, got_cd = p.cl_cd(alpha[-1] + 10.0)
        assert got_cl == pytest.approx(cl[-1])
        assert got_cd == pytest.approx(cd[-1])

    def test_matches_numpy_interp(self, small_polar):
        """The whole point of this implementation: bit-for-bit np.interp."""
        alpha, cl, cd, p = small_polar
        rng = np.random.default_rng(seed=7)
        test = rng.uniform(alpha[0] - 0.05, alpha[-1] + 0.05, size=200)
        got_cl = np.array([p.cl_cd(a)[0] for a in test])
        got_cd = np.array([p.cl_cd(a)[1] for a in test])
        ref_cl = np.interp(test, alpha, cl)
        ref_cd = np.interp(test, alpha, cd)
        # f64 division reformulation: a few ULPs is acceptable
        assert np.abs(got_cl - ref_cl).max() < 1e-13
        assert np.abs(got_cd - ref_cd).max() < 1e-13
