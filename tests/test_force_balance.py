"""Force-balance and frame-transform validation for the rotor models.

Coverage (gaps left by test_castles_gray.py, test_cyclic.py, test_pitt_peters.py):

1. **TestHubFrameTransforms** — Force / moment directions under non-identity
   R_hub (pitched and rolled hubs). Catches frame-transform regressions.

2. **TestHoverForceBalance** — Hover collective sweep against vehicle weight,
   anchored to Castles-Gray Table I Runs 3/4/5 (CT=0.004 @ θ_0.75R≈8.78°,
   1200 rpm). Verifies that F_world cancels gravity to within a few %, that
   M_orbital is essentially zero in axisymmetric hover, and that M_spin
   aligns with hub axis with the correct sign.

3. **TestNonHoverForceBalance** — Vehicle translating in +X under gravity,
   with and without ambient wind. Pitt-Peters settles to a non-trivial
   steady state; F_world still opposes hub axis as expected.

4. **TestSwashplateMapping** — End-to-end through cyclic_coeffs:
   helicopter-standard signs verified at the AeroResult interface, plus
   the swashplate phase parameter `swashplate_phase_deg` rotates the
   disk-tilt response by the expected amount.

Rotor: Castles-Gray 6-ft (Research/Castles_TN2474/, rotors/castles_gray_6ft/).
"""

import math
from pathlib import Path

import numpy as np
import pytest

from dynbem import RotorInputs, create_aero
from dynbem.rotor_definition import (
    ControlProperties,
    RotorDefinition,
    load as load_rotor,
)
from dynbem.rotor_state import PittPetersRotorState
from tests.helpers import pp_state

_ROTOR_YAML = str(
    Path(__file__).parent.parent / "rotors" / "castles_gray_6ft" / "rotor.yaml"
)
_RHO = 1.225
_GRAVITY = 9.81


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def cg_defn() -> RotorDefinition:
    return load_rotor(_ROTOR_YAML)


@pytest.fixture(scope="module")
def cg_pp(cg_defn):
    return create_aero(cg_defn, model="pitt_peters_jit")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_inputs(*,
                 collective_deg: float,
                 tilt_lon_deg: float = 0.0,
                 tilt_lat_deg: float = 0.0,
                 R_hub: np.ndarray | None = None,
                 v_hub_world=(0.0, 0.0, 0.0),
                 wind_world=(0.0, 0.0, 0.0),
                 omega_rpm: float = 1200.0) -> RotorInputs:
    if R_hub is None:
        R_hub = np.eye(3)
    omega_rad_s = omega_rpm * math.pi / 30.0
    return RotorInputs(
        collective_rad=math.radians(collective_deg),
        tilt_lon=math.radians(tilt_lon_deg),
        tilt_lat=math.radians(tilt_lat_deg),
        R_hub=R_hub,
        v_hub_world=np.asarray(v_hub_world, dtype=float),
        wind_world=np.asarray(wind_world, dtype=float),
        t=0.0,
        rho_kg_m3=1.225,
        omega_rad_s=omega_rad_s,
    )


def _R_pitch(theta_deg: float) -> np.ndarray:
    """Rotation that pitches the hub axis forward (about +Y body axis).

    Hub Z axis (originally +Z_world) rotates by +theta about +Y, so the
    hub axis in world becomes [sin θ, 0, cos θ] — points partly +X, partly +Z.
    Physically: nose-down rotor disk tilt (forward flight stick).
    """
    t = math.radians(theta_deg)
    c, s = math.cos(t), math.sin(t)
    return np.array([[ c, 0, s],
                     [ 0, 1, 0],
                     [-s, 0, c]], dtype=float)


def _R_roll(phi_deg: float) -> np.ndarray:
    """Rotation that rolls the hub axis to the right (about +X body axis).

    Hub axis in world becomes [0, -sin φ, cos φ] — points partly -Y (left),
    partly +Z. Physically: roll-left disk tilt for positive phi.

    (Our naming convention: a +φ roll about +X sends +Z toward +Y in the
    standard right-hand rule, but for an upward-thrust rotor the tilt of the
    *axis* sends thrust toward -Y. Read the components, not the name.)
    """
    p = math.radians(phi_deg)
    c, s = math.cos(p), math.sin(p)
    return np.array([[1, 0,  0],
                     [0, c, -s],
                     [0, s,  c]], dtype=float)


def _euler_to_steady(model, inputs: RotorInputs, *, rpm: float,
                     n_steps: int = 6000, dt: float = 0.001):
    """Integrate Pitt-Peters dynamic inflow to (approximate) steady state.

    Holds omega fixed at the given rpm; only lambda_0, lambda_c, lambda_s
    evolve. Returns (final_state, final_result).
    """
    omega = rpm * math.pi / 30.0
    state = pp_state()
    res = None
    for _ in range(n_steps):
        res, drv = model.compute_forces(inputs, state)
        state = PittPetersRotorState(
            lambda_0=state.lambda_0 + drv.lambda_0 * dt,
            lambda_c=state.lambda_c + drv.lambda_c * dt,
            lambda_s=state.lambda_s + drv.lambda_s * dt,
        )
    return state, res


def _disk_area(defn: RotorDefinition) -> float:
    R = defn.blade.radius_m
    return math.pi * R * R


# ===========================================================================
# 1. Hub-frame transforms — force/moment directions
# ===========================================================================

class TestHubFrameTransforms:
    """Verify F_world, M_orbital, M_spin transform with R_hub correctly."""

    def test_level_hub_thrust_along_minus_z(self, cg_pp):
        """R_hub = I, hover: F_world should be purely along -Z."""
        state, res = _euler_to_steady(
            cg_pp, _make_inputs(collective_deg=8.0), rpm=1200.0,
        )
        T = -res.F_world[2]
        assert T > 50.0, "Should have non-trivial hover thrust"
        # In-plane components must be effectively zero in axisymmetric hover.
        assert abs(res.F_world[0]) < 1e-6 * T
        assert abs(res.F_world[1]) < 1e-6 * T

    def test_pitched_hub_thrust_tilts_forward(self, cg_pp):
        """Pitching the hub forward by 10° tilts F_world by 10°.

        With R_hub = R_pitch(10°), the hub-axis vector becomes
        [sin 10°, 0, cos 10°]. F_world = -T * hub_axis, so F_world should
        have x-component ≈ -T·sin 10° and z-component ≈ -T·cos 10°.
        """
        R = _R_pitch(10.0)
        _, res = _euler_to_steady(
            cg_pp, _make_inputs(collective_deg=8.0, R_hub=R), rpm=1200.0,
        )
        T_mag = np.linalg.norm(res.F_world)
        s, c = math.sin(math.radians(10.0)), math.cos(math.radians(10.0))
        assert res.F_world[0] == pytest.approx(-T_mag * s, rel=1e-3, abs=0.5)
        assert abs(res.F_world[1]) < 1e-6 * T_mag
        assert res.F_world[2] == pytest.approx(-T_mag * c, rel=1e-3)

    def test_rolled_hub_thrust_tilts_sideways(self, cg_pp):
        """Rolling the hub by 10° gives F_world a Y component.

        With R_hub = R_roll(10°), hub axis = R @ [0,0,1] = [0, -sin 10°, cos 10°].
        F_world = -T * hub_axis = [0, +T·sin 10°, -T·cos 10°].
        """
        R = _R_roll(10.0)
        _, res = _euler_to_steady(
            cg_pp, _make_inputs(collective_deg=8.0, R_hub=R), rpm=1200.0,
        )
        T_mag = np.linalg.norm(res.F_world)
        s, c = math.sin(math.radians(10.0)), math.cos(math.radians(10.0))
        assert abs(res.F_world[0]) < 1e-6 * T_mag
        assert res.F_world[1] == pytest.approx(+T_mag * s, rel=1e-3, abs=0.5)
        assert res.F_world[2] == pytest.approx(-T_mag * c, rel=1e-3)

    def test_m_spin_always_aligned_with_hub_axis(self, cg_pp):
        """M_spin must point along hub-axis-in-world for any R_hub."""
        for label, R in [("level", np.eye(3)),
                         ("pitched", _R_pitch(15.0)),
                         ("rolled", _R_roll(20.0))]:
            _, res = _euler_to_steady(
                cg_pp, _make_inputs(collective_deg=8.0, R_hub=R), rpm=1200.0,
            )
            hub_axis = R @ np.array([0.0, 0.0, 1.0])
            # M_spin / Q_spin should equal hub_axis exactly.
            m_norm = np.linalg.norm(res.M_spin)
            if m_norm > 1e-9:
                np.testing.assert_allclose(
                    res.M_spin / m_norm, hub_axis, atol=1e-9,
                    err_msg=f"M_spin direction wrong for {label} hub",
                )
            # And M_spin = Q_spin * hub_axis (signed equality).
            np.testing.assert_allclose(
                res.M_spin, res.Q_spin * hub_axis, atol=1e-9,
                err_msg=f"M_spin sign/magnitude wrong for {label} hub",
            )

    def test_m_orbital_is_in_plane(self, cg_pp):
        """In-plane hub moments must be perpendicular to hub axis.

        M_orbital is built from r·dT·[sin ψ, cos ψ, 0] in hub frame; the
        z-component in hub frame is zero by construction. Rotating to world
        gives M_orbital · hub_axis_world = 0 exactly.
        """
        # Use a tilted hub and a cyclic input so M_orbital is non-trivial.
        R = _R_pitch(10.0)
        _, res = _euler_to_steady(
            cg_pp,
            _make_inputs(collective_deg=8.0, tilt_lon_deg=2.0, R_hub=R),
            rpm=1200.0,
        )
        hub_axis = R @ np.array([0.0, 0.0, 1.0])
        m_dot = float(np.dot(res.M_orbital, hub_axis))
        assert abs(m_dot) < 1e-9, (
            f"M_orbital has out-of-plane component {m_dot:.3e}; "
            f"should be perpendicular to hub axis."
        )


# ===========================================================================
# 2. Hover force balance — anchored to Castles-Gray Table I
# ===========================================================================

class TestHoverForceBalance:
    """Hover thrust balances gravity to within physical tolerance."""

    def test_hover_thrust_matches_castles_gray_run345(self, cg_pp, cg_defn):
        """θ_0.75R = 8.78° (Runs 3/4/5 mean), 1200 rpm → CT ≈ 0.004.

        The Castles-Gray rotor is untwisted so θ_0.75R == collective. Runs
        3/4/5 report CT = 0.004 at θ ∈ [8.66°, 8.86°] (mean 8.78°). Pitt-Peters
        with uniform inflow on this rotor typically gives CT within ~10% of
        the measured value (depending on the polar tail).
        """
        rpm = 1200.0
        omega = rpm * math.pi / 30.0
        R = cg_defn.blade.radius_m
        A = _disk_area(cg_defn)
        _, res = _euler_to_steady(
            cg_pp, _make_inputs(collective_deg=8.78), rpm=rpm,
        )
        T = -res.F_world[2]
        CT = T / (_RHO * A * (omega * R) ** 2)
        # 15% tolerance — Pitt-Peters with uniform inflow on a small rotor
        # against experimental CT.  The bundled XFOIL polar (NCrit=5,
        # Re=200k) gives Cl_α ≈ 6.0–6.5/rad in the operating alpha range vs
        # the test paper's measured 5.90/rad (Table VIII), which inflates
        # CT by ~13% at this collective.  See
        # Research/Castles_TN2474/naca0015_polar.md.
        assert CT == pytest.approx(0.004, rel=0.15), (
            f"CT @ θ=8.78°, 1200 rpm = {CT:.5f}; "
            f"Castles-Gray Runs 3/4/5 report CT≈0.004"
        )

    @pytest.mark.parametrize("weight_N", [25.0, 40.0, 60.0])
    def test_collective_solves_for_weight(self, cg_pp, cg_defn, weight_N):
        """Find collective that lifts the given weight; check F_world[2] = -W.

        Uses bisection on collective_deg over [0°, 16°] until thrust
        matches weight within 0.5%.
        """
        rpm = 1200.0
        lo, hi = 0.0, 16.0
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            _, res = _euler_to_steady(
                cg_pp, _make_inputs(collective_deg=mid), rpm=rpm,
            )
            T = -res.F_world[2]
            if T < weight_N:
                lo = mid
            else:
                hi = mid
            if abs(T - weight_N) / weight_N < 0.005:
                break
        # Verify the final solution
        _, res = _euler_to_steady(
            cg_pp, _make_inputs(collective_deg=mid), rpm=rpm,
        )
        assert -res.F_world[2] == pytest.approx(weight_N, rel=0.01)
        # Side checks at the equilibrium point:
        # F_world has no in-plane component (axisymmetric hover).
        assert abs(res.F_world[0]) < 1e-6 * weight_N
        assert abs(res.F_world[1]) < 1e-6 * weight_N
        # M_orbital essentially zero (axisymmetric hover).
        assert np.linalg.norm(res.M_orbital) < 1e-3
        # M_spin aligned with +Z hub axis (level rotor) — Q_spin > 0 in
        # powered hover, so M_spin[2] > 0.
        assert res.Q_spin > 0
        assert res.M_spin[2] == pytest.approx(res.Q_spin, rel=1e-12)

    def test_zero_collective_zero_thrust(self, cg_pp):
        """Symmetric airfoil + zero collective + zero cyclic ⇒ ~zero thrust.

        Castles-Gray rotor uses NACA 0015 (symmetric). The Pitt-Peters
        steady state at θ=0 should give thrust within a fraction of a
        Newton (limited by the polar's tabulated cl(0) precision).
        """
        _, res = _euler_to_steady(
            cg_pp, _make_inputs(collective_deg=0.0), rpm=1200.0,
        )
        T = -res.F_world[2]
        assert abs(T) < 1.0, f"|T| at θ=0 = {abs(T):.4f} N, expected ~0"


# ===========================================================================
# 3. Non-hover force balance — translating flight, wind, gravity
# ===========================================================================

class TestNonHoverForceBalance:
    """Translating flight and ambient wind produce sensible force tilts."""

    def test_climb_needs_more_collective_than_hover(self, cg_pp, cg_defn):
        """Same weight, vertical climb requires higher collective than hover."""
        rpm = 1200.0
        weight = 35.0  # N

        def thrust(coll_deg, v_climb):
            _, res = _euler_to_steady(
                cg_pp,
                _make_inputs(collective_deg=coll_deg,
                             v_hub_world=(0.0, 0.0, -v_climb)),
                rpm=rpm,
            )
            return -res.F_world[2]

        def solve_coll(v_climb):
            lo, hi = 0.0, 16.0
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                if thrust(mid, v_climb) < weight:
                    lo = mid
                else:
                    hi = mid
            return mid

        coll_hover = solve_coll(0.0)
        coll_climb = solve_coll(3.0)  # 3 m/s climb
        assert coll_climb > coll_hover + 0.3, (
            f"Climb collective {coll_climb:.3f}° not noticeably above hover "
            f"{coll_hover:.3f}° (expected ≥+0.3°)."
        )

    def test_headwind_into_hovering_rotor_tilts_thrust_axis_loaded(self, cg_pp):
        """Edgewise wind into a level rotor: F_world still mostly -Z, but
        induced moments appear via the BEM ψ-loop velocity asymmetry.

        Set wind in +X (north) into a level hovering rotor. F_world remains
        primarily downward (rotor produces lift regardless), but the rotor
        now sees forward-flight asymmetry → Pitt-Peters develops nonzero
        λ_c (Glauert) and λ_s (velocity-asymmetry rolling moment).
        """
        # 10 m/s wind from the north, vehicle stationary
        state, res = _euler_to_steady(
            cg_pp,
            _make_inputs(collective_deg=8.0,
                         wind_world=(10.0, 0.0, 0.0),
                         v_hub_world=(0.0, 0.0, 0.0)),
            rpm=1200.0,
        )
        T = -res.F_world[2]
        # Thrust still well above zero
        assert T > 30.0
        # Inflow asymmetry: at least one of λ_c/λ_s must be non-trivial.
        assert abs(state.lambda_c) + abs(state.lambda_s) > 1e-3, (
            f"Expected non-zero cyclic inflow under edgewise wind; got "
            f"λ_c={state.lambda_c:.4e}, λ_s={state.lambda_s:.4e}"
        )

    def test_wind_relative_to_hub_is_what_matters(self, cg_pp):
        """v_rel = wind − v_hub: equal-and-opposite vehicle motion and wind
        should give the same answer as the static case.

        Compare: (a) hover, wind=0, v_hub=0 vs (b) hover, wind=(V,0,0),
        v_hub=(V,0,0). The relative wind through the disk is the same in
        both — both should give identical thrust and inflow state.
        """
        state_a, res_a = _euler_to_steady(
            cg_pp, _make_inputs(collective_deg=8.0), rpm=1200.0,
        )
        state_b, res_b = _euler_to_steady(
            cg_pp,
            _make_inputs(collective_deg=8.0,
                         v_hub_world=(5.0, 0.0, 0.0),
                         wind_world=(5.0, 0.0, 0.0)),
            rpm=1200.0,
        )
        # The relative wind seen by the disk is identical, so results should
        # match to integration precision.
        np.testing.assert_allclose(res_a.F_world, res_b.F_world, atol=1e-6)
        assert res_a.Q_spin == pytest.approx(res_b.Q_spin, rel=1e-9)
        assert state_a.lambda_0 == pytest.approx(state_b.lambda_0, rel=1e-9)


# ===========================================================================
# 4. Swashplate mapping — sign + phase end-to-end through cyclic_coeffs
# ===========================================================================

def _defn_with_phase(base_defn: RotorDefinition, phase_deg: float) -> RotorDefinition:
    """Return a copy of base_defn with the given swashplate phase angle.

    Gain is unity so tilt_lon, tilt_lat map directly to blade pitch amplitudes
    (with the helicopter-standard sign convention from cyclic_coeffs).
    """
    ctrl = ControlProperties(
        swashplate_pitch_gain_rad=1.0,
        swashplate_phase_deg=phase_deg,
    )
    # The Rust pyclass isn't a @dataclass, so dataclasses.replace() doesn't
    # work -- construct a new RotorDefinition explicitly with the same
    # parts and the swapped control.
    return RotorDefinition(
        blade=base_defn.blade,
        airfoil=base_defn.airfoil,
        control=ctrl,
        inertia=base_defn.inertia,
        autorotation=base_defn.autorotation,
        name=base_defn.name,
        description=base_defn.description,
    )


class TestSwashplateMapping:
    """End-to-end cyclic input → AeroResult force/moment direction."""

    def test_tilt_lon_positive_produces_nose_down_moment(self, cg_pp):
        """tilt_lon > 0 ⇒ M_y_world < 0 (nose-down pitch moment, NED).

        Already covered in test_cyclic.py at the basic-sign level. Here we
        also check ORDER OF MAGNITUDE: at 2° cyclic on a 50N-thrust rotor
        with R≈0.91 m, the resulting pitch moment should be a few N·m.
        """
        _, res = _euler_to_steady(
            cg_pp,
            _make_inputs(collective_deg=8.0, tilt_lon_deg=2.0),
            rpm=1200.0,
        )
        T = -res.F_world[2]
        assert res.M_orbital[1] < -0.05, (
            f"M_y = {res.M_orbital[1]:.4e} N·m, expected nose-down (<0)"
        )
        # Order-of-magnitude check: |M_y| should be on the order of
        # (cyclic angle in rad) * T * R/2 — i.e. a few N·m for 2° on this rotor.
        scale = math.radians(2.0) * T * cg_pp.defn.blade.radius_m
        assert 0.01 * scale < abs(res.M_orbital[1]) < 1.5 * scale

    def test_tilt_lat_positive_produces_roll_right_moment(self, cg_pp):
        """tilt_lat > 0 ⇒ M_x_world > 0 (roll right, NED)."""
        _, res = _euler_to_steady(
            cg_pp,
            _make_inputs(collective_deg=8.0, tilt_lat_deg=2.0),
            rpm=1200.0,
        )
        assert res.M_orbital[0] > 0.05

    def test_cyclic_is_linear_for_small_inputs(self, cg_pp):
        """Doubling small cyclic input doubles the resulting moment."""
        _, res1 = _euler_to_steady(
            cg_pp, _make_inputs(collective_deg=8.0, tilt_lon_deg=0.5),
            rpm=1200.0,
        )
        _, res2 = _euler_to_steady(
            cg_pp, _make_inputs(collective_deg=8.0, tilt_lon_deg=1.0),
            rpm=1200.0,
        )
        # 5% tolerance — second-order BEM nonlinearity at small cyclic.
        assert res2.M_orbital[1] == pytest.approx(2.0 * res1.M_orbital[1], rel=0.05)

    def test_swashplate_phase_rotates_response_by_phi(self, cg_defn):
        """A 90° swashplate phase advance should make `tilt_lon` produce a
        roll moment (where it normally produces a pitch moment).

        At φ=0: tilt_lon > 0 ⇒ θ_1c < 0, θ_1s = 0  ⇒ peak pitch at ψ=π
        (tail) ⇒ M_y < 0 (nose-down).

        At φ=90°: from cyclic_coeffs,
          θ_1c = -tilt_lon·cos(π/2) - tilt_lat·sin(π/2) = -tilt_lat
          θ_1s = -tilt_lon·sin(π/2) + tilt_lat·cos(π/2) = -tilt_lon
        So tilt_lon now drives θ_1s < 0 (peak pitch at ψ=3π/2 = right side)
        ⇒ M_x < 0 (roll LEFT).
        """
        cg_phi0  = create_aero(_defn_with_phase(cg_defn,  0.0), model="pitt_peters_jit")
        cg_phi90 = create_aero(_defn_with_phase(cg_defn, 90.0), model="pitt_peters_jit")

        _, res0  = _euler_to_steady(
            cg_phi0,  _make_inputs(collective_deg=8.0, tilt_lon_deg=2.0),
            rpm=1200.0,
        )
        _, res90 = _euler_to_steady(
            cg_phi90, _make_inputs(collective_deg=8.0, tilt_lon_deg=2.0),
            rpm=1200.0,
        )

        # φ=0 baseline: tilt_lon > 0 ⇒ nose-down (M_y < 0), M_x ≈ 0
        assert res0.M_orbital[1] < -0.05
        assert abs(res0.M_orbital[0]) < 0.1 * abs(res0.M_orbital[1])

        # φ=90° response: moment swings to the roll axis with the sign flip
        # derived above (roll-LEFT).
        assert res90.M_orbital[0] < -0.05, (
            f"φ=90°: expected M_x < 0 (roll left), got {res90.M_orbital[0]:.4e}"
        )
        assert abs(res90.M_orbital[1]) < 0.1 * abs(res90.M_orbital[0])

        # And the magnitudes of pitch-moment-at-φ=0 and roll-moment-at-φ=90°
        # should match (it's the same physical asymmetry, rotated).
        assert abs(res90.M_orbital[0]) == pytest.approx(
            abs(res0.M_orbital[1]), rel=0.1
        )
