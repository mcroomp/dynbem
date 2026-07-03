"""Tests for quasi-static blade flapping (hub moment reduction).

Covers:
  1. FlapProperties hub_moment_factor correctness at known omega.
  2. With flap configured, hub moments are reduced vs rigid blade.
  3. Thrust is unchanged (flapping only affects moments, not axial force).
  4. Zero omega_nr (freely hinged) gives near-zero hub moment.
  5. Large omega_nr (stiff blade) gives nearly unchanged hub moment.
  6. YAML round-trip for flap section.
  7. All three models (pitt_peters, oye, bem) support flap.
"""

import math
import pytest
import numpy as np
import dynbem
from dynbem.rotor_definition import (
    load, FlapProperties, RotorDefinition,
    BladeGeometry, LinearPolarParameters, ControlProperties,
)
from dynbem._dynbem import RotorInputs, FlapProperties as _RustFlapProperties


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_defn(flap=None):
    blade = BladeGeometry(
        n_blades=2, radius_m=1.2, root_cutout_m=0.1,
        chord_m=0.15, twist_deg=0.0, n_elements=6,
    )
    airfoil = LinearPolarParameters(
        Re_design=500000, CL0=0.0, CL_alpha_per_rad=2 * math.pi,
        CD0=0.01, alpha_stall_deg=14.0,
    )
    control = ControlProperties(swashplate_pitch_gain_rad=1.0)
    return RotorDefinition(
        blade=blade, airfoil=airfoil, control=control,
        flap=flap, name="test_flap", description="",
    )


def make_inputs(omega=25.0, collective=0.1, tilt_lat=0.0, v_wind=5.0):
    return RotorInputs(
        omega_rad_s=omega, collective_rad=collective,
        tilt_lon=0.0, tilt_lat=tilt_lat, rho_kg_m3=1.225,
        v_hub_world=np.zeros(3), wind_world=np.array([v_wind, 0.0, 0.0]),
        R_hub=np.eye(3), t=0.0,
    )


def thrust_N(result):
    return -float(result.F_world[2])


# ---------------------------------------------------------------------------
# 1. FlapProperties.hub_moment_factor correctness
# ---------------------------------------------------------------------------

def test_hub_moment_factor_math():
    """Verify the reduction factor formula directly."""
    fp = _RustFlapProperties(I_blade_flap_kgm2=0.01, omega_nr_rad_s=10.0)
    omega = 30.0
    # nu_beta^2 = 1 + (10/30)^2 = 1 + 1/9 = 10/9
    # factor = (10/9 - 1) / (10/9) = (1/9) / (10/9) = 1/10
    expected = 1.0 / 10.0
    assert fp.hub_moment_factor(omega) == pytest.approx(expected, rel=1e-10)


def test_hub_moment_factor_zero_stiffness():
    """Freely-hinged blade: factor should be 0 (no moment transfer)."""
    fp = _RustFlapProperties(I_blade_flap_kgm2=0.01, omega_nr_rad_s=0.0)
    assert fp.hub_moment_factor(30.0) == pytest.approx(0.0, abs=1e-15)


def test_hub_moment_factor_stiff():
    """Very stiff blade: factor approaches 1."""
    fp = _RustFlapProperties(I_blade_flap_kgm2=0.01, omega_nr_rad_s=1000.0)
    # nu^2 = 1 + (1000/30)^2 ~ 1112; factor = 1111/1112 ~ 0.9991
    f = fp.hub_moment_factor(30.0)
    assert f > 0.99
    assert f < 1.0


# ---------------------------------------------------------------------------
# 2. Hub moments are reduced with flap vs rigid
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("model_name", ["pitt_peters", "oye", "bem"])
def test_flap_reduces_hub_moment(model_name):
    """With flap, Mx and My should be smaller than rigid blade."""
    # omega_nr = 8 rad/s, omega = 25 rad/s
    # nu^2 = 1 + (8/25)^2 = 1.1024, factor = 0.1024/1.1024 ~ 0.093
    flap = FlapProperties(I_blade_flap_kgm2=0.01, omega_nr_rad_s=8.0)

    defn_rigid = make_defn(flap=None)
    defn_flap = make_defn(flap=flap)

    model_rigid = dynbem.create_aero(defn_rigid, model=model_name)
    model_flap = dynbem.create_aero(defn_flap, model=model_name)

    inputs = make_inputs(tilt_lat=0.05)

    res_rigid, _ = model_rigid.compute_forces(inputs, model_rigid.initial_rotor_state())
    res_flap, _ = model_flap.compute_forces(inputs, model_flap.initial_rotor_state())

    mx_rigid = abs(float(res_rigid.m_hub_world[0]))
    mx_flap = abs(float(res_flap.m_hub_world[0]))

    # Flap should reduce moments significantly
    assert mx_flap < mx_rigid * 0.5, (
        f"flap moment {mx_flap:.4f} should be < 50% of rigid {mx_rigid:.4f}"
    )


# ---------------------------------------------------------------------------
# 3. Thrust unchanged by flapping
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("model_name", ["pitt_peters", "oye", "bem"])
def test_flap_does_not_change_thrust(model_name):
    """Flapping should not affect total thrust."""
    flap = FlapProperties(I_blade_flap_kgm2=0.01, omega_nr_rad_s=8.0)

    defn_rigid = make_defn(flap=None)
    defn_flap = make_defn(flap=flap)

    model_rigid = dynbem.create_aero(defn_rigid, model=model_name)
    model_flap = dynbem.create_aero(defn_flap, model=model_name)

    inputs = make_inputs(tilt_lat=0.05)

    res_rigid, _ = model_rigid.compute_forces(inputs, model_rigid.initial_rotor_state())
    res_flap, _ = model_flap.compute_forces(inputs, model_flap.initial_rotor_state())

    t_rigid = thrust_N(res_rigid)
    t_flap = thrust_N(res_flap)

    assert t_flap == pytest.approx(t_rigid, rel=1e-10)


# ---------------------------------------------------------------------------
# 4. Freely-hinged blade gives near-zero hub moment
# ---------------------------------------------------------------------------

def test_free_hinge_zero_moment():
    """omega_nr=0 means all moment absorbed by flapping."""
    flap = FlapProperties(I_blade_flap_kgm2=0.01, omega_nr_rad_s=0.0)
    defn = make_defn(flap=flap)
    model = dynbem.create_aero(defn, model="pitt_peters")
    inputs = make_inputs(tilt_lat=0.1)

    res, _ = model.compute_forces(inputs, model.initial_rotor_state())
    mx = abs(float(res.M_orbital[0]))
    my = abs(float(res.M_orbital[1]))

    assert mx < 1e-12
    assert my < 1e-12


# ---------------------------------------------------------------------------
# 5. Very stiff blade gives nearly unchanged moment
# ---------------------------------------------------------------------------

def test_stiff_blade_unchanged():
    """Large omega_nr -> factor ~1, moments unchanged."""
    flap = FlapProperties(I_blade_flap_kgm2=0.01, omega_nr_rad_s=10000.0)
    defn_rigid = make_defn(flap=None)
    defn_flap = make_defn(flap=flap)

    model_rigid = dynbem.create_aero(defn_rigid, model="pitt_peters")
    model_flap = dynbem.create_aero(defn_flap, model="pitt_peters")

    inputs = make_inputs(tilt_lat=0.05)

    res_rigid, _ = model_rigid.compute_forces(inputs, model_rigid.initial_rotor_state())
    res_flap, _ = model_flap.compute_forces(inputs, model_flap.initial_rotor_state())

    mx_rigid = float(res_rigid.M_orbital[0])
    mx_flap = float(res_flap.M_orbital[0])

    assert mx_flap == pytest.approx(mx_rigid, rel=1e-4)


# ---------------------------------------------------------------------------
# 6. YAML round-trip
# ---------------------------------------------------------------------------

def test_yaml_roundtrip():
    """Flap section loads from YAML correctly."""
    from dynbem.rotor_definition import loads

    yaml_text = """\
name: flap_test
rotor:
  n_blades: 2
  radius_m: 1.2
  root_cutout_m: 0.1
  chord_m: 0.15
  n_elements: 6
airfoil:
  Re_design: 500000
  CL0: 0.0
  CL_alpha_per_rad: 6.28
  CD0: 0.01
  alpha_stall_deg: 14.0
control:
  swashplate_pitch_gain_rad: 1.0
flap:
  I_blade_flap_kgm2: 0.012
  omega_nr_rad_s: 7.5
"""
    defn = loads(yaml_text)
    assert defn.flap is not None
    assert defn.flap.I_blade_flap_kgm2 == pytest.approx(0.012)
    assert defn.flap.omega_nr_rad_s == pytest.approx(7.5)


def test_yaml_no_flap_section():
    """Absent flap section -> flap is None (rigid blade)."""
    from dynbem.rotor_definition import loads

    yaml_text = """\
name: no_flap
rotor:
  n_blades: 2
  radius_m: 1.2
  root_cutout_m: 0.1
  chord_m: 0.15
airfoil:
  Re_design: 500000
  CL0: 0.0
  CL_alpha_per_rad: 6.28
  CD0: 0.01
  alpha_stall_deg: 14.0
control:
  swashplate_pitch_gain_rad: 1.0
"""
    defn = loads(yaml_text)
    assert defn.flap is None


# ---------------------------------------------------------------------------
# 7. Quantitative check of reduction factor in model output
# ---------------------------------------------------------------------------

def test_moment_reduction_matches_factor():
    """The ratio of flap/rigid moments should match hub_moment_factor."""
    omega = 25.0
    omega_nr = 8.0
    # factor = (omega_nr/omega)^2 / (1 + (omega_nr/omega)^2)
    nu2 = 1.0 + (omega_nr / omega) ** 2
    expected_factor = (nu2 - 1.0) / nu2

    flap = FlapProperties(I_blade_flap_kgm2=0.01, omega_nr_rad_s=omega_nr)
    defn_rigid = make_defn(flap=None)
    defn_flap = make_defn(flap=flap)

    model_rigid = dynbem.create_aero(defn_rigid, model="pitt_peters")
    model_flap = dynbem.create_aero(defn_flap, model="pitt_peters")

    inputs = make_inputs(omega=omega, tilt_lat=0.05)

    res_rigid, _ = model_rigid.compute_forces(inputs, model_rigid.initial_rotor_state())
    res_flap, _ = model_flap.compute_forces(inputs, model_flap.initial_rotor_state())

    mx_rigid = float(res_rigid.M_orbital[0])
    mx_flap = float(res_flap.M_orbital[0])

    if abs(mx_rigid) > 1e-10:
        actual_ratio = mx_flap / mx_rigid
        assert actual_ratio == pytest.approx(expected_factor, rel=1e-6)
