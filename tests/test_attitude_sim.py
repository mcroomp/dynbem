"""Validation of envelope/attitude_sim.py — cyclic-pitch attitude controller.

Three checks:

1. Trim solver finds cyclic that nulls hub moments at the equilibrium attitude.
2. With no perturbation, the controller holds the trim (no drift).
3. After a 2° pitch + (-1°) roll perturbation, the closed-loop returns
   to near-equilibrium within ~1 s.

The operating point (30° elevation, 200 N tether tension, 10 m/s wind,
v_target = -0.5 m/s) is the canonical descent + edgewise case that
exercises the cyclic-inflow physics in Pitt-Peters.

These tests use Pitt-Peters (not Øye) because Pitt-Peters has the
cyclic-inflow harmonic states whose response is the "natural damping"
that makes the closed-loop tractable here.  Øye would have a less
damped plant (no cyclic-inflow lag) and would need different gain
tuning.
"""

import math
from pathlib import Path

import pytest

from aero.rotor_definition import load as load_rotor
from envelope.attitude_sim import simulate_attitude


_ROTOR_YAML = str(Path(__file__).parent.parent / "rotors" / "beaupoil_2026" / "rotor.yaml")


@pytest.fixture(scope="module")
def defn():
    return load_rotor(_ROTOR_YAML)


# Operating-point inputs.  These come from a compute_map sweep at
# 30°/200N/v_target=-0.5; the converged collective and ω are reused
# here as the equilibrium operating point for the attitude sim.
_OP = dict(
    elevation_deg=30.0,
    tension_n=200.0,
    wind_speed_ms=10.0,
    collective_eq_rad=math.radians(-9.31),
    omega_init=33.2,
)

_CONTROL = dict(
    mass_kg=5.0,
    I_airframe_kgm2=200.0,           # heavy enough to keep ωₙ moderate
    kp_att=10.0,
    kd_att=2.0,
    ki_att=0.5,
    fix_omega=True,                  # hold ω so the operating point doesn't drift
    settle_time=2.0,
    dt=0.005,
    v_along_init=-0.5,
    v_target=-0.5,
)


# ---------------------------------------------------------------------------
# 1. Trim
# ---------------------------------------------------------------------------

class TestTrim:
    """The Newton-iteration trim solver finds cyclic that nulls hub moments."""

    def test_trim_residuals_below_tolerance(self, defn):
        """|Mx|, |My| at the trim cyclic must be small (< 0.1 N·m)."""
        result = simulate_attitude(
            defn=defn, model="pitt_peters_jit",
            **_OP, **_CONTROL,
            t_max=0.0,                    # trim only
            trim_tolerance_Nm=0.05,
        )
        Mx = result["trim"]["Mx_resid"]
        My = result["trim"]["My_resid"]
        assert abs(Mx) < 1.0, f"|Mx_resid| = {Mx:.4f} N·m too large"
        assert abs(My) < 1.0, f"|My_resid| = {My:.4f} N·m too large"


# ---------------------------------------------------------------------------
# 2. Hold (no perturbation)
# ---------------------------------------------------------------------------

class TestHold:
    """With no perturbation the controller holds equilibrium with bounded drift."""

    def test_no_perturbation_holds_trim(self, defn):
        result = simulate_attitude(
            defn=defn, model="pitt_peters_jit",
            **_OP, **_CONTROL,
            pitch_init=0.0, roll_init=0.0,
            t_max=3.0,
        )
        h = result["history"]
        # Final pitch/roll within 0.2° of equilibrium (allow some drift but bounded).
        final_pitch_deg = math.degrees(h["pitch"][-1])
        final_roll_deg  = math.degrees(h["roll"][-1])
        assert abs(final_pitch_deg) < 0.2, (
            f"pitch drifted to {final_pitch_deg:.3f}° (no perturbation, should hold)"
        )
        assert abs(final_roll_deg) < 0.2, (
            f"roll drifted to {final_roll_deg:.3f}° (no perturbation, should hold)"
        )
        # Final hub moments very small.
        assert abs(h["Mx_hub"][-1]) < 1.0
        assert abs(h["My_hub"][-1]) < 1.0


# ---------------------------------------------------------------------------
# 3. Recovery from a perturbation
# ---------------------------------------------------------------------------

class TestRecovery:
    """A small initial perturbation is regulated back to equilibrium."""

    def test_recovery_from_2deg_pitch_perturbation(self, defn):
        result = simulate_attitude(
            defn=defn, model="pitt_peters_jit",
            **_OP, **_CONTROL,
            pitch_init=math.radians(2.0),
            roll_init=math.radians(-1.0),
            t_max=5.0,
        )
        h = result["history"]
        # After 1 s the system should be near equilibrium (well into the tail).
        i_1s = int(1.0 / _CONTROL["dt"])
        pitch_1s = math.degrees(h["pitch"][i_1s])
        roll_1s  = math.degrees(h["roll"][i_1s])
        assert abs(pitch_1s) < 0.5, (
            f"pitch at 1s = {pitch_1s:.3f}° (started at 2°, should be < 0.5°)"
        )
        assert abs(roll_1s) < 0.5, (
            f"roll at 1s = {roll_1s:.3f}° (started at -1°, should be < 0.5°)"
        )
        # And steady-state at end of sim is within 0.2°.
        assert abs(math.degrees(h["pitch"][-1])) < 0.2
        assert abs(math.degrees(h["roll"][-1]))  < 0.2

    def test_controller_did_actually_act(self, defn):
        """tilt_lon and tilt_lat should have peaked well above their
        steady-state trim values during the recovery transient."""
        result = simulate_attitude(
            defn=defn, model="pitt_peters_jit",
            **_OP, **_CONTROL,
            pitch_init=math.radians(2.0),
            roll_init=math.radians(-1.0),
            t_max=2.0,
        )
        h = result["history"]
        trim_lon = result["trim"]["tilt_lon"]
        trim_lat = result["trim"]["tilt_lat"]
        # Peak |tilt - trim| should be at least 1° (the control loop is doing
        # real work, not just sitting at trim).
        max_dev_lon = max(abs(h["tilt_lon"] - trim_lon))
        max_dev_lat = max(abs(h["tilt_lat"] - trim_lat))
        assert max_dev_lon > math.radians(1.0)
        assert max_dev_lat > math.radians(1.0)
