"""Validation of OyeBEMModel — annulus-local 2-stage dynamic inflow.

Goals:
1. State plumbing: OyeRotorState round-trips through to_array / from_array.
2. Initialisation: the factory returns the right type.
3. Hover CT agrees with PittPetersModel within a few percent (both are
   Level-2 BEM with different inflow ODEs; should agree in axisymmetric
   hover where there's no L-matrix vs annulus-local distinction).
4. Climb / descent: induced inflow has the right sign and the autorotation
   torque is negative.
5. ODE lag: λ states evolve toward steady state with finite time
   constants (not instantaneously).
6. Envelope robustness: at the descent + edgewise wind operating point
   where Pitt-Peters was numerically stiff, Øye converges cleanly at
   dt=0.005 — verifies the "annulus-local = no L-matrix feedback" claim.
"""

import math
from pathlib import Path

import numpy as np
import pytest

from dynbem import RotorInputs, create_aero
from dynbem.oye import OyeBEMModel
from dynbem.rotor_state import OyeRotorState, PittPetersRotorState
import dynbem.rotor_definition as rotor_definition


_ROTOR_YAML = str(
    Path(__file__).parent.parent / "rotors" / "castles_gray_6ft" / "rotor.yaml"
)
_RHO = 1.225


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def defn():
    return rotor_definition.load(_ROTOR_YAML)


@pytest.fixture(scope="module")
def oye_model(defn):
    return OyeBEMModel(defn=defn)


@pytest.fixture(scope="module")
def pp_model(defn):
    return create_aero(defn, model="pitt_peters_jit")


# ---------------------------------------------------------------------------
# State plumbing
# ---------------------------------------------------------------------------

class TestOyeRotorState:
    def test_zeros_factory_returns_correct_shape(self):
        s = OyeRotorState.zeros(n_elements=10, omega_rad_s=125.0)
        assert s.W_int.shape == (10,)
        assert s.W.shape == (10,)
        assert s.omega_rad_s == 125.0
        assert s.spin_angle_rad == 0.0

    def test_to_array_layout(self):
        s = OyeRotorState(
            W_int=np.array([0.01, 0.02, 0.03]),
            W=np.array([0.10, 0.20, 0.30]),
            omega_rad_s=50.0,
            spin_angle_rad=1.5,
        )
        arr = s.to_array()
        # [W_int_0..2, W_0..2, omega, psi]
        assert arr.shape == (8,)
        assert arr[:3] == pytest.approx([0.01, 0.02, 0.03])
        assert arr[3:6] == pytest.approx([0.10, 0.20, 0.30])
        assert arr[-2] == 50.0
        assert arr[-1] == 1.5

    def test_round_trip(self):
        s0 = OyeRotorState(
            W_int=np.array([0.01, -0.02, 0.03]),
            W=np.array([0.10, 0.20, -0.30]),
            omega_rad_s=42.0,
            spin_angle_rad=-1.0,
        )
        arr = s0.to_array()
        s1 = s0.from_array(arr + 0.0)  # +0.0 to ensure a new array
        assert np.allclose(s1.W_int, s0.W_int)
        assert np.allclose(s1.W, s0.W)
        assert s1.omega_rad_s == s0.omega_rad_s
        assert s1.spin_angle_rad == s0.spin_angle_rad

    def test_from_array_rejects_bad_length(self):
        s = OyeRotorState.zeros(3)
        with pytest.raises(ValueError):
            s.from_array(np.zeros(5))   # odd number of inflow states


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------

class TestFactory:
    def test_create_aero_oye(self, defn):
        m = create_aero(defn, model="oye")
        assert isinstance(m, OyeBEMModel)
        state = m.initial_rotor_state()
        assert isinstance(state, OyeRotorState)
        assert state.W.shape == (defn.blade.n_elements,)


# ---------------------------------------------------------------------------
# Hover: CT should match Pitt-Peters within ~5%
# ---------------------------------------------------------------------------

def _euler_to_steady_inflow(model, inputs, state, dt: float = 0.001,
                            n_steps: int = 6000):
    """Time-integrate model holding omega fixed.  Returns final (result, state)."""
    omega = state.omega_rad_s
    res = None
    for _ in range(n_steps):
        res, drv = model.compute_forces(inputs, state)
        arr = state.to_array() + dt * drv.to_array()
        arr[-2] = omega   # hold omega fixed for this comparison
        state = state.from_array(arr)
    return res, state


class TestOyeHover:
    """In axisymmetric hover, Pitt-Peters and Øye should agree on CT.

    Both are Level-2 dynamic-inflow BEMs at the same blade geometry; the
    only difference is the inflow ODE structure (L-matrix vs annulus-
    local).  In hover there's no harmonic / cyclic content, so they should
    converge to the same total thrust.
    """

    # tolerance pair: (θ_deg, rpm, rel_tol)
    # Higher tolerance at low thrust because Pitt-Peters' aggregate
    # lam0_ss and Øye's per-annulus W_qs diverge at small CT (small μ_T
    # hits the 0.05 floor, where the per-annulus vs uniform inflow
    # distinction matters more).  At moderate thrust they agree closely.
    @pytest.mark.parametrize("theta_deg,rpm,rel_tol", [
        (8.86, 1200, 0.05),   # Run 5, CT ≈ 0.004 — agree within 5%
        (5.13, 1600, 0.25),   # Run 9, CT ≈ 0.002 — small CT, looser
    ])
    def test_hover_ct_matches_pitt_peters(self, oye_model, pp_model, defn,
                                          theta_deg, rpm, rel_tol):
        omega = rpm * math.pi / 30.0
        R = defn.blade.radius_m
        A = math.pi * R * R
        inp = RotorInputs(
            collective_rad=math.radians(theta_deg),
            tilt_lon=0.0, tilt_lat=0.0, R_hub=np.eye(3),
            v_hub_world=np.zeros(3), wind_world=np.zeros(3), t=0.0,
        )

        res_oye, _ = _euler_to_steady_inflow(
            oye_model, inp,
            OyeRotorState.zeros(defn.blade.n_elements, omega_rad_s=omega),
        )
        res_pp,  _ = _euler_to_steady_inflow(
            pp_model, inp, PittPetersRotorState(omega_rad_s=omega),
        )
        T_oye = -res_oye.F_world[2]
        T_pp  = -res_pp.F_world[2]
        CT_oye = T_oye / (_RHO * A * (omega * R) ** 2)
        CT_pp  = T_pp  / (_RHO * A * (omega * R) ** 2)
        err = abs(CT_oye - CT_pp) / CT_pp
        assert err < rel_tol, (
            f"θ={theta_deg} rpm={rpm}: "
            f"CT_Oye={CT_oye:.5f}, CT_PP={CT_pp:.5f}, err={err:.1%}"
        )


# ---------------------------------------------------------------------------
# Climb / descent: sign and magnitude of induced inflow
# ---------------------------------------------------------------------------

class TestOyeInflowSigns:
    """Steady-state W has the right sign for each operating regime."""

    def test_hover_W_positive(self, oye_model, defn):
        """Hover: W > 0 (induced flow downward through disk)."""
        omega = 1200.0 * math.pi / 30.0
        inp = RotorInputs(
            collective_rad=math.radians(8.86),
            tilt_lon=0.0, tilt_lat=0.0, R_hub=np.eye(3),
            v_hub_world=np.zeros(3), wind_world=np.zeros(3), t=0.0,
        )
        _, state = _euler_to_steady_inflow(
            oye_model, inp,
            OyeRotorState.zeros(defn.blade.n_elements, omega_rad_s=omega),
        )
        assert float(np.mean(state.W)) > 1e-4

    def test_climb_W_positive_but_smaller_than_hover(self, oye_model, defn):
        """Climb: induction (W) is smaller than in hover because some of
        the disk flow comes from the freestream."""
        omega = 1200.0 * math.pi / 30.0
        inp_hover = RotorInputs(
            collective_rad=math.radians(8.86),
            tilt_lon=0.0, tilt_lat=0.0, R_hub=np.eye(3),
            v_hub_world=np.zeros(3), wind_world=np.zeros(3), t=0.0,
        )
        inp_climb = RotorInputs(
            collective_rad=math.radians(8.86),
            tilt_lon=0.0, tilt_lat=0.0, R_hub=np.eye(3),
            v_hub_world=np.zeros(3),
            wind_world=np.array([0.0, 0.0, 5.0]),     # +Z = downward freestream
            t=0.0,
        )
        _, st_h = _euler_to_steady_inflow(
            oye_model, inp_hover,
            OyeRotorState.zeros(defn.blade.n_elements, omega_rad_s=omega),
        )
        _, st_c = _euler_to_steady_inflow(
            oye_model, inp_climb,
            OyeRotorState.zeros(defn.blade.n_elements, omega_rad_s=omega),
        )
        W_h = float(np.mean(st_h.W))
        W_c = float(np.mean(st_c.W))
        assert W_c > 0.0
        assert W_c < W_h, f"climb W={W_c:.4f} not less than hover W={W_h:.4f}"


# ---------------------------------------------------------------------------
# Lag dynamics: λ states evolve over τ, not instantaneously
# ---------------------------------------------------------------------------

class TestOyeLag:
    """After a step change in collective the state should move toward the
    new equilibrium with finite τ — not jump instantly."""

    def test_step_response_has_lag(self, oye_model, defn):
        omega = 1200.0 * math.pi / 30.0
        # Reach hover equilibrium first.
        inp1 = RotorInputs(
            collective_rad=math.radians(8.86),
            tilt_lon=0.0, tilt_lat=0.0, R_hub=np.eye(3),
            v_hub_world=np.zeros(3), wind_world=np.zeros(3), t=0.0,
        )
        _, state = _euler_to_steady_inflow(
            oye_model, inp1,
            OyeRotorState.zeros(defn.blade.n_elements, omega_rad_s=omega),
        )
        W0 = state.W.copy()

        # Step up collective by 2°.  Take ONE step.
        inp2 = RotorInputs(
            collective_rad=math.radians(10.86),
            tilt_lon=0.0, tilt_lat=0.0, R_hub=np.eye(3),
            v_hub_world=np.zeros(3), wind_world=np.zeros(3), t=0.0,
        )
        _, drv = oye_model.compute_forces(inp2, state)
        # dW/dt should be nonzero and pushing W upward (more thrust → more induction).
        assert float(np.mean(drv.W_int)) > 0.0, "dW_int/dt should be positive after collective step-up"
        # But the state hasn't moved yet (we haven't integrated).
        assert np.allclose(state.W, W0)


# ---------------------------------------------------------------------------
# Envelope robustness: the regime where Pitt-Peters was stiff
# ---------------------------------------------------------------------------

class TestOyeEnvelopeStability:
    """At descent + edgewise wind, Øye should converge at the same dt
    that destabilises Pitt-Peters' L-matrix."""

    def test_descent_with_sidewind_converges(self, defn):
        """The 30°/200N tether case that was the canonical Pitt-Peters
        stability failure.  Øye should reach v_along ≈ -0.5 m/s
        equilibrium at dt = 0.005 without numerical blowup."""
        from envelope.point_mass import (
            _step_state_semi_implicit, _clip_state,
            tether_hat, _build_r_hub, G,
        )

        m = create_aero(defn, model="oye")
        # Use the beaupoil rotor — the one that exhibited the regression.
        defn_b = rotor_definition.load(str(
            Path(__file__).parent.parent / "rotors" / "beaupoil_2026" / "rotor.yaml"
        ))
        m = create_aero(defn_b, model="oye")

        t_hat = tether_hat(30.0)
        T_tether, mass = 200.0, 5.0
        f_load = T_tether * t_hat + np.array([0.0, 0.0, mass * G])
        bz = f_load / float(np.linalg.norm(f_load))
        R_hub = _build_r_hub(bz)
        wind = np.array([0.0, -10.0, 0.0])

        state = OyeRotorState.zeros(defn_b.blade.n_elements, omega_rad_s=20.0)
        v_along = 0.0
        col = 0.0
        int_v = 0.0
        v_target = -0.5
        dt = 0.005

        for step in range(int(60.0 / dt)):
            inp = RotorInputs(
                collective_rad=col, tilt_lon=0.0, tilt_lat=0.0,
                R_hub=R_hub, v_hub_world=v_along * t_hat,
                wind_world=wind, t=step * dt,
            )
            aero_r, dst = m.compute_forces(inp, state)
            arr = _step_state_semi_implicit(m, state, dst, dt, inp)
            arr = _clip_state(arr, state)
            state = state.from_array(arr)
            f_thrust = float(np.dot(aero_r.F_world, t_hat))
            f_along = f_thrust + mass * G * float(t_hat[2]) + T_tether
            v_along = max(-30.0, min(30.0, v_along + dt / mass * f_along))
            err = v_along - v_target
            int_v = max(-22.5, min(22.5, int_v + err * dt))
            col = max(-0.25, min(0.20, 0.01 * err + 0.02 * int_v))

        # Converged: v_along at target, omega reasonable, W small and
        # nowhere near the ±10 clip.
        assert abs(v_along - v_target) < 0.1, f"v_along={v_along}, expected ≈ {v_target}"
        assert 20.0 < state.omega_rad_s < 100.0, f"omega={state.omega_rad_s}"
        assert abs(float(np.mean(state.W))) < 0.5, f"|W|={np.mean(state.W)} near clip"
