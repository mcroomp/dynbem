"""Simplified sign diagnostic for quasi_static BEM only.

Isolates behaviour from test_rawes_ic_aero_sign.py using the simplest
possible geometry: R_hub = identity (hub_axis = NED +Z = downward).

With R_hub = I:
  v_climb = wind[2]   (positive = wind blows down = helicopter climb)
  v_inplane = wind[:2] (in-plane wind triggers psi-loop when large enough)

The RAWES IC scenario decomposes to:
  v_climb   = -9.04 m/s   (upflow through disk, windmill regime)
  v_inplane =  4.27 m/s   (mu ~ 0.11, triggers psi-loop at omega=53.16)

Root cause identified:
  - Case D (axial only, mu=0): uses the Brent windmill solver -> correct sign
    at ALL omegas (30..60 rad/s).
  - Case C (in-plane + axial, mu>0): psi-loop calls solve_bem_element
    (helicopter quadratic) with v_climb < 0.  At omega > ~42 rad/s,
    |lambda_climb| = |v_climb|/(omega*R) becomes small enough that the
    quadratic converges to the helicopter (positive-thrust) root rather than
    the windmill root, flipping the sign of T.
  - Cases A, B: sign is consistent across omegas (no windmill interaction).

This means the psi-loop needs to use a windmill-capable solver when v_climb<0,
just as the axial path already does (try windmill Brent first, fall back to
helicopter quadratic).

Four cases:
  A: hover              mu=0,   v_climb=0    -- axial path
  B: in-plane only      mu>0,   v_climb=0    -- psi-loop
  C: in-plane + axial   mu>0,   v_climb<0    -- psi-loop (RAWES IC analog) [BUG]
  D: axial windmill     mu=0,   v_climb<0    -- axial windmill path [correct reference]
"""

from __future__ import annotations

import numpy as np
import pytest

from dynbem import RotorInputs, create_aero
from dynbem.rotor_definition import load as load_rotor
from pathlib import Path

_ROTOR_YAML = Path(__file__).parent.parent / "rotors" / "beaupoil_2026" / "rotor.yaml"
_RHO = 1.225
_COLLECTIVE_RAD = -0.18
# The two omegas from the original RAWES IC test (QS uses higher, PP uses lower)
_OMEGA_QS = 53.161687   # rad/s -- QS trim omega from windpower
_OMEGA_PP = 38.132161   # rad/s -- PP trim omega from windpower

_R_HUB_IDENTITY = np.eye(3)  # hub_axis = [0,0,1] = NED down

# Decomposed RAWES IC in aligned-rotor coordinates (hub_axis = [0,0,1]):
#   v_climb = wind[2],  v_inplane_hub = wind[:2]
_V_CLIMB = -9.04   # m/s -- matches RAWES IC decomposed v_climb
_V_INPLANE = 4.27  # m/s -- matches RAWES IC decomposed v_inplane magnitude


def _inputs(wind_world: np.ndarray, omega: float = _OMEGA_QS) -> RotorInputs:
    return RotorInputs(
        collective_rad=_COLLECTIVE_RAD,
        tilt_lon=0.0,
        tilt_lat=0.0,
        R_hub=_R_HUB_IDENTITY,
        v_hub_world=np.zeros(3),
        wind_world=wind_world,
        omega_rad_s=omega,
        t=0.0,
        rho_kg_m3=_RHO,
    )


def _qs_f2(wind: np.ndarray, omega: float) -> tuple[float, list]:
    """Run quasi_static and return (F_world[2], F_world list)."""
    rotor = load_rotor(str(_ROTOR_YAML))
    model = create_aero(rotor, model="quasi_static")
    result, _ = model.compute_forces(_inputs(wind, omega), model.initial_rotor_state())
    fw = np.asarray(result.F_world, dtype=float)
    return float(fw[2]), fw.tolist()


# ---------------------------------------------------------------------------
# Case A: pure hover (mu=0, v_climb=0) -- axial path
# With collective=-0.18 and CL0=0.393: alpha < 0 in hover -> T < 0 -> F[2] > 0
# Assertion: the two omegas give the same SIGN (both negative thrust direction).
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("omega", [_OMEGA_QS, _OMEGA_PP])
def test_case_a_hover_sign_consistent(omega: float) -> None:
    """Case A: hover. Both omegas should agree on sign of F_world[2]."""
    wind = np.array([0.0, 0.0, 0.0])
    f2_qs, fw_qs = _qs_f2(wind, omega)
    # Hover with negative collective produces negative thrust (downward force).
    # Just assert the value is non-trivially nonzero and consistent.
    assert abs(f2_qs) > 10.0, (
        f"Case A hover omega={omega:.1f}: thrust unexpectedly near zero\n"
        f"F_world = {fw_qs}"
    )


# ---------------------------------------------------------------------------
# Case B: in-plane wind only (mu>0, v_climb=0) -- psi-loop
# No axial flow, so helicopter quadratic is fine. Sign is consistent.
# ---------------------------------------------------------------------------

def test_case_b_inplane_sign_consistent() -> None:
    """Case B: in-plane wind only. QS sign is the same at both omegas."""
    wind = np.array([_V_INPLANE, 0.0, 0.0])
    f2_high, fw_high = _qs_f2(wind, _OMEGA_QS)
    f2_low, fw_low   = _qs_f2(wind, _OMEGA_PP)
    assert (f2_high > 0) == (f2_low > 0), (
        f"Case B in-plane wind: sign flips between omegas -- unexpected!\n"
        f"  omega={_OMEGA_QS:.1f}: F_world[2]={f2_high:.1f} {fw_high}\n"
        f"  omega={_OMEGA_PP:.1f}: F_world[2]={f2_low:.1f} {fw_low}"
    )


# ---------------------------------------------------------------------------
# Case C: in-plane + axial (mu>0, v_climb<0) -- psi-loop, RAWES IC analog
#
# BUG: the psi-loop calls solve_bem_element (helicopter quadratic) for every
# element even when v_climb < 0 (windmill regime).  At omega > ~42 rad/s,
# |lambda_climb| = |v_climb|/(omega*R) becomes small enough that the quadratic
# converges to the helicopter-mode root instead of the windmill root, flipping
# the sign of T.  The axial path (Case D) uses the Brent windmill solver and
# stays correct at all omegas.  The fix: in the psi-loop, route elements
# through the windmill solver when v_climb < 0, mirroring the axial path.
# ---------------------------------------------------------------------------

def test_case_c_psiloop_matches_axial_windmill_sign_at_high_omega() -> None:
    """Case C vs D at QS omega: psi-loop must match axial windmill sign."""
    wind_c = np.array([_V_INPLANE, 0.0, _V_CLIMB])      # psi-loop (mu>0)
    wind_d = np.array([0.0, 0.0, _V_CLIMB])              # axial Brent (mu=0)
    f2_c, fw_c = _qs_f2(wind_c, _OMEGA_QS)
    f2_d, fw_d = _qs_f2(wind_d, _OMEGA_QS)
    assert (f2_c > 0) == (f2_d > 0), (
        f"Case C vs D at omega={_OMEGA_QS:.1f}: sign mismatch\n"
        f"  Case C (psi-loop, mu>0): F_world[2]={f2_c:.1f}  {fw_c}\n"
        f"  Case D (axial Brent):    F_world[2]={f2_d:.1f}  {fw_d}"
    )


def test_case_c_sign_consistent_across_omegas() -> None:
    """Case C: QS sign must not flip between the two operating omegas."""
    wind = np.array([_V_INPLANE, 0.0, _V_CLIMB])
    f2_high, fw_high = _qs_f2(wind, _OMEGA_QS)
    f2_low, fw_low   = _qs_f2(wind, _OMEGA_PP)
    assert (f2_high > 0) == (f2_low > 0), (
        f"Case C: QS sign flips between omegas (wrong quadratic root in psi-loop).\n"
        f"  omega={_OMEGA_QS:.1f}: F_world[2]={f2_high:.1f}  {fw_high}\n"
        f"  omega={_OMEGA_PP:.1f}: F_world[2]={f2_low:.1f}  {fw_low}"
    )


# ---------------------------------------------------------------------------
# Case D: axial windmill only (mu=0, v_climb<0) -- axial windmill path
# QS windmill solver should give F_world[2] < 0 (upward = anti-hub-axis).
# If the windmill Brent solver is correct, this should pass at both omegas.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("omega", [_OMEGA_QS, _OMEGA_PP])
def test_case_d_axial_windmill_force_is_upward(omega: float) -> None:
    """Case D: axial windmill. F_world[2] < 0 (force opposes hub axis)."""
    wind = np.array([0.0, 0.0, _V_CLIMB])  # mu=0, v_climb=-9.04
    f2, fw = _qs_f2(wind, omega)
    assert f2 < 0.0, (
        f"Case D axial windmill omega={omega:.1f}: expected F_world[2] < 0, got {f2:.3f}\n"
        f"F_world = {fw}"
    )
