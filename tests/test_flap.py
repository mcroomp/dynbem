"""Tests for the servo-flap pitch actuation model.

Covers:
  1. Direct-mechanical blade (no servo-flap) is unchanged -- all models run.
  2. ServoFlapActuation round-trips through rotor YAML loader (test_flap_rotor).
  3. ServoFlapActuation round-trips through _to_rust.
  4. Servo-flap command (tilt_lat) changes hub moment vs. direct mechanical.
  5. Less damping -> more feathering authority.
  6. OyeBEMModel also runs with servo-flap actuation configured.
  7. Raw feathering lag: tilt_lat couples into pitch axis without compensation.
  8. Servo mode disables the direct cyclic pitch path.
  9. Full beaupoil servoflaps rotor runs compute_forces without error.
"""

import math
import pytest
import numpy as np
import dynbem
from dynbem.rotor_definition import (
    load, ServoFlapActuation, ServoFlapGeometry, RotorDefinition,
    BladeGeometry, LinearPolarParameters, ControlProperties,
)
from dynbem._dynbem import (
    ServoFlapActuation as _RustServoFlapActuation,
    RotorInputs,
)
from pathlib import Path

FEATHERING_ROTOR_YAML = Path(__file__).parent.parent / "rotors" / "test_flap_rotor" / "rotor.yaml"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_minimal_defn(servoflap=False, damper=5.0):
    blade = BladeGeometry(
        n_blades=2, radius_m=1.2, root_cutout_m=0.1,
        chord_m=0.15, twist_deg=0.0, n_elements=6,
    )
    airfoil = LinearPolarParameters(
        Re_design=500000, CL0=0.0, CL_alpha_per_rad=2 * math.pi,
        CD0=0.01, alpha_stall_deg=14.0,
    )
    control = ControlProperties(swashplate_pitch_gain_rad=1.0)
    servoflap_act = None
    if servoflap:
        servoflap_act = ServoFlapActuation(
            I_theta_kgm2=0.001, damper_Nms_per_rad=damper,
            flap=ServoFlapGeometry(
                C_M_delta_per_rad=-1.5, r_inner_m=0.5, r_outer_m=1.2,
            ),
            ac_offset_m=0.0,
        )
    return RotorDefinition(
        blade=blade, airfoil=airfoil, control=control,
        servoflap=servoflap_act, name="test", description="",
    )


def make_inputs(omega=25.0, collective=0.1, tilt_lon=0.0, tilt_lat=0.0, v_wind=0.0):
    return RotorInputs(
        omega_rad_s=omega, collective_rad=collective,
        tilt_lon=tilt_lon, tilt_lat=tilt_lat, rho_kg_m3=1.225,
        v_hub_world=np.zeros(3), wind_world=np.array([v_wind, 0.0, 0.0]),
        R_hub=np.eye(3), t=0.0,
    )


def thrust_N(result):
    return -float(result.F_world[2])


# ---------------------------------------------------------------------------
# 1. Direct-mechanical blade is unchanged
# ---------------------------------------------------------------------------

def test_rigid_blade_pitt_peters_runs():
    defn = make_minimal_defn(servoflap=False)
    model = dynbem.create_aero(defn, model="pitt_peters")
    state = model.initial_rotor_state()
    result, _ = model.compute_forces(make_inputs(), state)
    assert math.isfinite(thrust_N(result))
    assert thrust_N(result) > 0


# ---------------------------------------------------------------------------
# 2. YAML loader round-trip
# ---------------------------------------------------------------------------

def test_feathering_rotor_yaml_loads():
    defn = load(str(FEATHERING_ROTOR_YAML))
    assert defn.servoflap is not None, "test_flap_rotor should have pitch_actuation.servoflap: section"
    act = defn.servoflap
    assert act.I_theta_kgm2 == pytest.approx(0.0012)
    assert act.damper_Nms_per_rad == pytest.approx(0.5)
    assert act.ac_offset_m == pytest.approx(0.0)
    assert act.flap is not None
    assert act.flap.C_M_delta_per_rad == pytest.approx(-1.5)
    assert act.flap.r_inner_m == pytest.approx(1.2)
    assert act.flap.r_outer_m == pytest.approx(2.5)


# ---------------------------------------------------------------------------
# 3. ServoFlapActuation round-trips through _to_rust
# ---------------------------------------------------------------------------

def test_feathering_rust_roundtrip():
    act = ServoFlapActuation(
        I_theta_kgm2=0.002, damper_Nms_per_rad=3.0,
        flap=ServoFlapGeometry(C_M_delta_per_rad=-1.2, r_inner_m=0.4, r_outer_m=1.1),
        ac_offset_m=0.005,
    )
    rust = act._to_rust()
    assert isinstance(rust, _RustServoFlapActuation)
    assert rust.I_theta_kgm2 == pytest.approx(0.002)
    assert rust.damper_Nms_per_rad == pytest.approx(3.0)
    assert rust.ac_offset_m == pytest.approx(0.005)
    assert rust.flap.C_M_delta_per_rad == pytest.approx(-1.2)


# ---------------------------------------------------------------------------
# 4. Servo command changes hub moment vs. direct mechanical
# ---------------------------------------------------------------------------

def test_servo_command_changes_hub_moment():
    defn_feat = make_minimal_defn(servoflap=True)
    defn_rigid = make_minimal_defn(servoflap=False)
    model_feat = dynbem.create_aero(defn_feat, model="pitt_peters")
    model_rigid = dynbem.create_aero(defn_rigid, model="pitt_peters")
    inputs = make_inputs(tilt_lat=0.1, v_wind=3.0)
    res_feat, _ = model_feat.compute_forces(inputs, model_feat.initial_rotor_state())
    res_rigid, _ = model_rigid.compute_forces(inputs, model_rigid.initial_rotor_state())
    mx_feat = float(res_feat.m_hub_world[0])
    mx_rigid = float(res_rigid.m_hub_world[0])
    assert mx_feat != pytest.approx(mx_rigid, abs=1e-3), (
        f"servo-flap model should differ from rigid: feat={mx_feat:.4f} rigid={mx_rigid:.4f}"
    )


# ---------------------------------------------------------------------------
# 5. Less damping -> more feathering authority
# ---------------------------------------------------------------------------

def test_authority_increases_with_less_damping():
    inputs = make_inputs(tilt_lat=0.1, v_wind=0.0)
    defn_r = make_minimal_defn(servoflap=False)
    model_r = dynbem.create_aero(defn_r, model="pitt_peters")
    res_r, _ = model_r.compute_forces(inputs, model_r.initial_rotor_state())
    my_rigid = float(res_r.m_hub_world[1])

    results = {}
    for damper, label in [(1.0, "lo"), (10.0, "hi")]:
        defn = make_minimal_defn(servoflap=True, damper=damper)
        model = dynbem.create_aero(defn, model="pitt_peters")
        res, _ = model.compute_forces(inputs, model.initial_rotor_state())
        results[label] = abs(float(res.m_hub_world[1]) - my_rigid)

    assert results["lo"] > results["hi"], (
        f"lo_damper dev={results['lo']:.4f} should > hi_damper dev={results['hi']:.4f}"
    )


# ---------------------------------------------------------------------------
# 6. OyeBEMModel runs with servo-flap actuation
# ---------------------------------------------------------------------------

def test_oye_with_feathering_runs():
    defn = make_minimal_defn(servoflap=True)
    model = dynbem.create_aero(defn, model="oye")
    result, _ = model.compute_forces(make_inputs(v_wind=3.0), model.initial_rotor_state())
    assert math.isfinite(thrust_N(result))


# ---------------------------------------------------------------------------
# 7. Raw feathering lag: tilt_lat couples into pitch axis without compensation
# ---------------------------------------------------------------------------

def test_phase_lag_lateral_command():
    """With no internal phase compensation, 1/rev feathering lag is visible.

    A lateral command produces a notable pitch-axis response relative to
    the direct-mechanical case.
    """
    defn_feat = make_minimal_defn(servoflap=True)
    defn_rigid = make_minimal_defn(servoflap=False)
    inputs = make_inputs(tilt_lat=0.1, v_wind=0.0)
    model_feat = dynbem.create_aero(defn_feat, model="pitt_peters")
    model_rigid = dynbem.create_aero(defn_rigid, model="pitt_peters")
    res_feat, _ = model_feat.compute_forces(inputs, model_feat.initial_rotor_state())
    res_rigid, _ = model_rigid.compute_forces(inputs, model_rigid.initial_rotor_state())
    assert math.isfinite(thrust_N(res_feat))
    assert math.isfinite(thrust_N(res_rigid))

    my_feat = float(res_feat.M_orbital[1])
    my_rigid = float(res_rigid.M_orbital[1])
    assert my_feat != pytest.approx(my_rigid, abs=1e-6)


def test_servo_mode_disables_direct_cyclic_pitch_path():
    """In servo-flap mode, cyclic must act through flap moments only.

    If C_M_delta is zero, cyclic flap commands should have negligible effect.
    """
    defn = make_minimal_defn(servoflap=True)
    defn.servoflap.flap.C_M_delta_per_rad = 0.0

    model = dynbem.create_aero(defn, model="pitt_peters")
    state = model.initial_rotor_state()
    res0, _ = model.compute_forces(make_inputs(tilt_lat=0.0, v_wind=3.0), state)
    res1, _ = model.compute_forces(make_inputs(tilt_lat=0.1, v_wind=3.0), state)

    mx0 = float(res0.M_orbital[0])
    mx1 = float(res1.M_orbital[0])
    assert mx1 == pytest.approx(mx0, abs=1e-6), (
        f"cyclic effect should be near zero when flap moment gain is zero: {mx0:.6f} vs {mx1:.6f}"
    )


# ---------------------------------------------------------------------------
# 9. Full beaupoil servo-flap rotor runs compute_forces
# ---------------------------------------------------------------------------

def test_beaupoil_feathering_compute_forces():
    defn = load(str(FEATHERING_ROTOR_YAML))
    model = dynbem.create_aero(defn, model="pitt_peters")
    result, _ = model.compute_forces(
        make_inputs(omega=20.0, collective=0.05, v_wind=5.0, tilt_lat=0.05),
        model.initial_rotor_state(),
    )
    assert math.isfinite(thrust_N(result))
