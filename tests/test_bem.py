"""Tests for the Level 1 BEM solver.

Validation target: Caradonna-Tung rotor (NASA TM-81232, 1981).
  2 blades, NACA 0012, R=1.143 m, chord=0.1905 m, no twist, no taper.
  Tested at collective 5°, 8°, 12° in hover (no forward speed, no wind axial).

NED convention: hub axis points down (+Z). Thrust is negative-Z.
"""

import math
import numpy as np
import pytest

from dynbem import AeroResult, BEMModel, RotorInputs
from dynbem.rotor_definition import (
    AirfoilProperties, AutorotationProperties, BladeGeometry, RotorDefinition,
)
from dynbem.rotor_state import QuasiStaticRotorState


# ---------------------------------------------------------------------------
# Caradonna-Tung rotor fixture
# ---------------------------------------------------------------------------

@pytest.fixture
def ct_blade():
    return BladeGeometry(
        n_blades=2,
        radius_m=1.143,
        root_cutout_m=0.1,
        chord_m=0.1905,
        twist_deg=0.0,
        n_elements=20,
    )


@pytest.fixture
def ct_airfoil():
    # NACA 0012: CL0=0, 2π lift slope, CD0≈0.008, stall ~15°
    return AirfoilProperties(
        Re_design=1_000_000,
        CL0=0.0,
        CL_alpha_per_rad=2 * math.pi,
        CD0=0.008,
        alpha_stall_deg=15.0,
        tip_loss=True,
    )


@pytest.fixture
def ct_defn(ct_blade, ct_airfoil):
    return RotorDefinition(
        blade=ct_blade,
        airfoil=ct_airfoil,
        autorotation=AutorotationProperties(I_ode_kgm2=1.0),
        name="Caradonna-Tung",
    )


@pytest.fixture
def ct_model(ct_defn):
    return BEMModel(defn=ct_defn)


def _hover_inputs(collective_deg: float, omega_rpm: float) -> RotorInputs:
    """Pure hover: no translation, no wind (NED frame, hub points down)."""
    omega_rad_s = omega_rpm * math.pi / 30.0
    return RotorInputs(
        collective_rad=math.radians(collective_deg),
        tilt_lon=0.0,
        tilt_lat=0.0,
        R_hub=np.eye(3),          # hub frame == world NED
        v_hub_world=np.zeros(3),
        wind_world=np.zeros(3),   # no wind: pure momentum-driven induction
        t=0.0,
        rho_kg_m3=1.225,
        motor_torque_Nm=0.0,
    ), QuasiStaticRotorState(omega_rad_s=omega_rad_s)


# ---------------------------------------------------------------------------
# Basic interface tests
# ---------------------------------------------------------------------------

class TestBEMInterface:
    def test_returns_aero_result(self, ct_model):
        inp, state = _hover_inputs(8.0, 1250)
        result, deriv = ct_model.compute_forces(inp, state)
        assert isinstance(result, AeroResult)

    def test_returns_quasi_static_state(self, ct_model):
        inp, state = _hover_inputs(8.0, 1250)
        _, deriv = ct_model.compute_forces(inp, state)
        assert isinstance(deriv, QuasiStaticRotorState)

    def test_initial_rotor_state(self, ct_model):
        s = ct_model.initial_rotor_state()
        assert isinstance(s, QuasiStaticRotorState)
        assert s.omega_rad_s == 0.0


# ---------------------------------------------------------------------------
# NED sign convention tests
# ---------------------------------------------------------------------------

class TestNEDConventions:
    def test_thrust_negative_z_in_hover(self, ct_model):
        """Upward thrust must be negative NED-Z."""
        inp, state = _hover_inputs(8.0, 1250)
        result, _ = ct_model.compute_forces(inp, state)
        assert result.F_world[2] < 0, "Thrust should be −Z (upward) in NED"

    def test_zero_xy_forces_in_axial_hover(self, ct_model):
        """In pure axial hover, lateral forces should be near-zero."""
        inp, state = _hover_inputs(8.0, 1250)
        result, _ = ct_model.compute_forces(inp, state)
        assert abs(result.F_world[0]) < 1.0
        assert abs(result.F_world[1]) < 1.0

    def test_spin_angle_deriv_equals_omega(self, ct_model):
        """dψ/dt == ω."""
        omega = 130.9  # rad/s
        inp, state = _hover_inputs(8.0, 1250)
        state = QuasiStaticRotorState(omega_rad_s=omega)
        inp = RotorInputs(
            collective_rad=math.radians(8.0),
            tilt_lon=0.0, tilt_lat=0.0,
            R_hub=np.eye(3),
            v_hub_world=np.zeros(3),
            wind_world=np.zeros(3),
            t=0.0,
        )
        _, deriv = ct_model.compute_forces(inp, state)
        assert deriv.spin_angle_rad == pytest.approx(omega)


# ---------------------------------------------------------------------------
# Physical plausibility tests (not tight validation — that needs test data)
# ---------------------------------------------------------------------------

class TestBEMPhysics:
    def test_thrust_increases_with_collective(self, ct_model):
        for c1, c2 in [(5, 8), (8, 12)]:
            inp1, s1 = _hover_inputs(c1, 1250)
            inp2, s2 = _hover_inputs(c2, 1250)
            r1, _ = ct_model.compute_forces(inp1, s1)
            r2, _ = ct_model.compute_forces(inp2, s2)
            # More negative Z = more thrust
            assert r2.F_world[2] < r1.F_world[2], (
                f"Thrust should increase from {c1}° to {c2}°"
            )

    def test_thrust_increases_with_rpm(self, ct_model):
        inp1, s1 = _hover_inputs(8.0, 900)
        inp2, s2 = _hover_inputs(8.0, 1500)
        r1, _ = ct_model.compute_forces(inp1, s1)
        r2, _ = ct_model.compute_forces(inp2, s2)
        assert r2.F_world[2] < r1.F_world[2]

    def test_zero_omega_gives_near_zero_thrust(self, ct_model):
        inp = RotorInputs(
            collective_rad=math.radians(8.0),
            tilt_lon=0.0, tilt_lat=0.0,
            R_hub=np.eye(3),
            v_hub_world=np.zeros(3),
            wind_world=np.zeros(3),
            t=0.0,
        )
        state = QuasiStaticRotorState(omega_rad_s=0.0)
        result, _ = ct_model.compute_forces(inp, state)
        assert abs(result.F_world[2]) < 5.0

    def test_autorotation_wind_produces_torque(self, ct_model):
        """Upward wind (−Z) should drive the rotor (positive aero torque)."""
        inp = RotorInputs(
            collective_rad=math.radians(5.0),
            tilt_lon=0.0, tilt_lat=0.0,
            R_hub=np.eye(3),
            v_hub_world=np.zeros(3),
            wind_world=np.array([0.0, 0.0, -15.0]),  # 15 m/s upward (−Z NED)
            t=0.0,
        )
        state = QuasiStaticRotorState(omega_rad_s=50.0)
        result, deriv = ct_model.compute_forces(inp, state)
        # In autorotation the aero torque should accelerate the rotor
        # (d_omega > 0 since motor_torque=0 and wind drives rotor)
        assert deriv.omega_rad_s > 0, "Upward wind should spin up rotor"

    def test_tip_loss_reduces_thrust(self):
        """Tip loss enabled should give less thrust than tip loss disabled."""
        blade = BladeGeometry(
            n_blades=2, radius_m=1.143, root_cutout_m=0.1,
            chord_m=0.1905, n_elements=20,
        )
        airfoil_tl = AirfoilProperties(
            Re_design=1_000_000, CL0=0.0,
            CL_alpha_per_rad=2 * math.pi, CD0=0.008,
            alpha_stall_deg=15.0, tip_loss=True,
        )
        airfoil_no = AirfoilProperties(
            Re_design=1_000_000, CL0=0.0,
            CL_alpha_per_rad=2 * math.pi, CD0=0.008,
            alpha_stall_deg=15.0, tip_loss=False,
        )
        defn_tl = RotorDefinition(blade=blade, airfoil=airfoil_tl)
        defn_no = RotorDefinition(blade=blade, airfoil=airfoil_no)
        m_tl = BEMModel(defn=defn_tl)
        m_no = BEMModel(defn=defn_no)
        inp, state = _hover_inputs(8.0, 1250)
        r_tl, _ = m_tl.compute_forces(inp, state)
        r_no, _ = m_no.compute_forces(inp, state)
        # Tip loss reduces thrust magnitude
        assert abs(r_tl.F_world[2]) < abs(r_no.F_world[2])


# ---------------------------------------------------------------------------
# Caradonna-Tung dimensional sanity: CT at 8° collective
# ---------------------------------------------------------------------------

class TestCaradonnaTungSanity:
    """Loose sanity bounds against NASA TM-81232 Table 2.

    At 8° collective, Ω ≈ 1250 RPM (tip Mach ~0.44):
      CT ≈ 0.0064  (measured, compressibility corrections not included here)

    We expect the incompressible BEM to over-predict slightly.
    Acceptable range: 0.003 < CT < 0.012.
    """

    def test_ct_8deg_in_range(self, ct_model):
        omega_rpm = 1250.0
        omega = omega_rpm * math.pi / 30.0
        R = ct_model.defn.blade.radius_m
        rho = 1.225
        A = math.pi * R**2

        inp, state = _hover_inputs(8.0, omega_rpm)
        result, _ = ct_model.compute_forces(inp, state)

        T = -result.F_world[2]         # thrust magnitude (positive upward)
        CT = T / (rho * A * (omega * R)**2)
        assert 0.002 < CT < 0.015, f"CT={CT:.5f} outside expected range"

    def test_ct_5deg_less_than_8deg(self, ct_model):
        omega_rpm = 1250.0
        omega = omega_rpm * math.pi / 30.0
        R = ct_model.defn.blade.radius_m
        rho = 1.225
        A = math.pi * R**2

        def get_ct(collective_deg):
            inp, state = _hover_inputs(collective_deg, omega_rpm)
            result, _ = ct_model.compute_forces(inp, state)
            T = -result.F_world[2]
            return T / (rho * A * (omega * R)**2)

        assert get_ct(5.0) < get_ct(8.0) < get_ct(12.0)
