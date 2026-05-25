import pytest

from dynbem.rotor_definition import (
    AirfoilProperties,
    AutorotationProperties,
    BladeGeometry,
    ControlProperties,
    InertiaProperties,
    RotorDefinition,
)
from tests.helpers import make_airfoil, make_blade, make_control


@pytest.fixture
def blade():
    return make_blade(n_blades=3, radius_m=5.0, root_cutout_m=0.5, chord_m=0.3)


@pytest.fixture
def airfoil():
    return make_airfoil(
        CL0=0.0,
        CL_alpha_per_rad=5.7,
        CD0=0.01,
        alpha_stall_deg=15.0,
        Re_design=500_000,
    )


@pytest.fixture
def control():
    return make_control(swashplate_pitch_gain_rad=0.1)


@pytest.fixture
def defn(blade, airfoil, control):
    return RotorDefinition(blade=blade, airfoil=airfoil, control=control)


class TestBladeGeometry:
    def test_span(self, blade):
        assert blade.span_m == pytest.approx(4.5)

    def test_r_cp(self, blade):
        assert blade.r_cp_m == pytest.approx(0.5 + (2 / 3) * 4.5)

    def test_solidity(self, blade):
        import math
        assert blade.solidity == pytest.approx(3 * 0.3 / (math.pi * 5.0))

    def test_disk_area(self, blade):
        import math
        assert blade.disk_area_m2 == pytest.approx(math.pi * (5.0**2 - 0.5**2))

    def test_validate_ok(self, blade):
        assert blade.validate() == []

    def test_validate_negative_span(self):
        b = make_blade(n_blades=3, radius_m=0.4, root_cutout_m=0.5, chord_m=0.3)
        issues = b.validate()
        assert any(i.level == "ERROR" and "span" in i.field for i in issues)

    def test_validate_zero_blades(self):
        b = make_blade(n_blades=0, radius_m=5.0, root_cutout_m=0.5, chord_m=0.3)
        issues = b.validate()
        assert any(i.level == "ERROR" and "n_blades" in i.field for i in issues)


class TestAirfoilProperties:
    def test_validate_ok(self, airfoil):
        assert airfoil.validate() == []

    def test_validate_negative_cd0(self):
        a = make_airfoil(
            CL0=0.0, CL_alpha_per_rad=5.7, CD0=-0.01,
            alpha_stall_deg=15.0, Re_design=500_000,
        )
        issues = a.validate()
        assert any("CD0" in i.field for i in issues)


class TestRotorDefinition:
    def test_forwarding_span(self, defn):
        assert defn.span_m == defn.blade.span_m

    def test_forwarding_solidity(self, defn):
        assert defn.solidity == defn.blade.solidity

    def test_validate_ok(self, defn):
        assert defn.validate() == []

    def test_validate_aggregates_blade_and_airfoil(self):
        bad_blade = make_blade(n_blades=0, radius_m=5.0, root_cutout_m=0.5, chord_m=0.3)
        bad_airfoil = make_airfoil(
            CL0=0.0, CL_alpha_per_rad=5.7, CD0=-0.01, alpha_stall_deg=15.0,
            Re_design=500_000,
        )
        d = RotorDefinition(blade=bad_blade, airfoil=bad_airfoil)
        issues = d.validate()
        assert len(issues) >= 2

    def test_inertia_defaults(self, defn):
        assert defn.inertia.mass_kg is None
        assert list(defn.inertia.I_body_kgm2) == []

    def test_autorotation_defaults(self, defn):
        assert defn.autorotation.omega_eq_rad_s is None

    def test_frozen(self, defn):
        with pytest.raises((TypeError, AttributeError)):
            defn.name = "modified"
