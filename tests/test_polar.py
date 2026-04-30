import math
import pytest

from aero.polar import LinearPolar


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
        from aero.rotor_definition import AirfoilProperties
        props = AirfoilProperties(
            Re_design=500_000,
            CL0=0.1,
            CL_alpha_per_rad=5.5,
            CD0=0.015,
            alpha_stall_deg=14.0,
        )
        p = LinearPolar.from_properties(props)
        assert p.CL0 == pytest.approx(0.1)
        assert p.CL_alpha_per_rad == pytest.approx(5.5)
        assert p.CD0 == pytest.approx(0.015)
        assert p.alpha_stall_rad == pytest.approx(math.radians(14.0))
