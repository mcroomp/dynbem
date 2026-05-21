"""Cross-check: dynbem.bem vs CCBlade on real rotors.

Two cross-check pairs, both driven by their verification module in
sampled mode (one BEM-driver path, fast tests, full-sweep script
elsewhere for re-baselining bounds):

NREL Phase VI (TestDynbemVsCCBladePhaseVI)
------------------------------------------
The right cross-check: a real twisted/tapered wind turbine in the
operating envelope CCBlade was calibrated for.  S809 airfoil, 2 blades,
R = 5.029 m, 72 RPM, +3 deg tip pitch, V_wind 5..25 m/s.

dynbem.bem dispatches dynamically: when v_climb < 0 (wind blowing
axially through the disk) it tries the windmill BEM iteration first
(a / (1-a) = sigma_r Cn / (4 F sin^2 phi)).  If that converges to a
valid windmill state (0 < a < 0.5, AoA below stall) the result is
used; otherwise the existing helicopter momentum quadratic takes
over.  No fixture flag -- the regime is read from the flow state.

Full-sweep numbers from running the verifier with no --sample (21
operating points):

    CT err: median 2.0 %, mean 1.9 %, RMSE 1.9 %, max 2.3 %
    CQ err: median 2.4 %, mean 2.4 %, RMSE 2.4 %, max 2.9 %

Both BEMs now agree within ~3 % at every operating point.  The
windmill solver uses Brent's-method-on-phi (Ning 2014) to find the
physical inflow angle robustly, with Buhl's quadratic taking over
from classical momentum theory at a = 0.4 (smooth transition by
construction).

Beaupoil RAWES rotor (TestDynbemVsCCBladeBeaupoil)
---------------------------------------------------
4-blade untwisted helicopter-style rotor (R = 2.5 m, chord 0.2 m,
SG6040 airfoil) driven by axial wind across 25 operating points
spanning V_wind = 5..16 m/s and Omega = 100..300 rpm.  Initially
this was a smoke test because the old fixed-point windmill solver
collapsed into the stalled basin of attraction on a cambered
airfoil at zero pitch.  With Brent's-method-on-phi (Ning 2014) on
the windmill side now finding the clean root, agreement is in line
with Phase VI:

    CT err: median 1.4 %, mean 1.4 %, RMSE 1.4 %, max 1.8 %
    CQ err: median 2.7 %, mean 3.3 %, RMSE 3.7 %, max 11.1 %

The 11.1 % CQ max sits at the autorotation crossing
(V_wind=5 m/s, Omega=192.4 rpm) where |CQ_cc| = 2e-5: the absolute
error is ~2e-6 N*m/normalisation, but relative error blows up at
the near-zero denominator.  The CQ envelope test filters those
near-zero points to avoid spurious failures from that effect.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from verification.dynbem_vs_ccblade_beaupoil import (  # noqa: E402
    CCBLADE_CSV as BEAUPOIL_CSV,
    ROTOR_YAML as BEAUPOIL_YAML,
    run_survey as run_beaupoil_survey,
)
from verification.dynbem_vs_ccblade_nrel_phase_vi import (  # noqa: E402
    CCBLADE_CSV as PHASE_VI_CSV,
    ROTOR_YAML as PHASE_VI_YAML,
    run_survey as run_phase_vi_survey,
)


# ---------------------------------------------------------------------------
# NREL Phase VI -- the real cross-check
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def phase_vi_survey():
    if not PHASE_VI_CSV.exists():
        pytest.skip(
            f"missing CCBlade reference CSV: {PHASE_VI_CSV}.  "
            f"Build it via verification/ccblade_docker/ (see README).")
    if not PHASE_VI_YAML.exists():
        pytest.skip(f"missing rotor fixture: {PHASE_VI_YAML}")
    # sample=None -> all 21 rows; the sweep is cheap (~1s).
    return run_phase_vi_survey(sample=None)


class TestDynbemVsCCBladePhaseVI:
    def test_thrust_sign_and_finiteness(self, phase_vi_survey):
        """Every CCBlade row + dynbem prediction is finite and positive
        (both BEMs agree the rotor reacts upwind into the airflow)."""
        for c in phase_vi_survey.comparisons:
            for name, v in (("CT_db", c.CT_db), ("CQ_db", c.CQ_db),
                            ("T_db_N", c.T_db_N), ("Q_db_Nm", c.Q_db_Nm)):
                assert math.isfinite(v), (
                    f"non-finite {name}={v} at U={c.U_wind_ms}")
            assert c.T_cc_N > 0 and c.T_db_N > 0
            assert c.Q_cc_Nm > 0 and c.Q_db_Nm > 0, (
                f"both BEMs should predict turbine extraction at "
                f"U={c.U_wind_ms}: CCBlade Q={c.Q_cc_Nm:+.1f}, "
                f"dynbem Q={c.Q_db_Nm:+.1f}")

    def test_ct_envelope_within_3pct(self, phase_vi_survey):
        """Every per-point CT error stays within +/-3 % across the
        whole sweep.  Current full-sweep figure is max 2.3 %."""
        offenders = [c for c in phase_vi_survey.comparisons if c.ct_err >= 0.03]
        assert not offenders, (
            "CT err > 3 %:\n" + "\n".join(
                f"  V={c.U_wind_ms} CT_cc={c.CT_cc:.4f} CT_db={c.CT_db:.4f} "
                f"err={c.ct_err:.1%}" for c in offenders))

    def test_cq_envelope_within_4pct(self, phase_vi_survey):
        """Per-point CQ error stays within +/-4 % everywhere.
        Current full-sweep figure is max 2.9 %."""
        offenders = [c for c in phase_vi_survey.comparisons if c.cq_err >= 0.04]
        assert not offenders, (
            "CQ err > 4 %:\n" + "\n".join(
                f"  V={c.U_wind_ms} CQ_cc={c.CQ_cc:+.5f} CQ_db={c.CQ_db:+.5f} "
                f"err={c.cq_err:.1%}" for c in offenders))


# ---------------------------------------------------------------------------
# Beaupoil RAWES rotor -- known-open investigation, smoke test only
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def beaupoil_survey():
    if not BEAUPOIL_CSV.exists():
        pytest.skip(f"missing CCBlade reference CSV: {BEAUPOIL_CSV}")
    if not BEAUPOIL_YAML.exists():
        pytest.skip(f"missing rotor fixture: {BEAUPOIL_YAML}")
    # 25 operating points; full sweep is fast (well under 1s).
    return run_beaupoil_survey(sample=None)


class TestDynbemVsCCBladeBeaupoil:
    def test_results_are_finite_and_signs_agree(self, beaupoil_survey):
        """No NaN/Inf; both BEMs predict positive thrust at every point."""
        for c in beaupoil_survey.comparisons:
            for name, v in (("CT_db", c.CT_db), ("CQ_db", c.CQ_db),
                            ("T_db_N", c.T_db_N), ("Q_db_Nm", c.Q_db_Nm)):
                assert math.isfinite(v), (
                    f"non-finite {name}={v} at U={c.U_wind_ms}")
            assert c.T_cc_N > 0 and c.T_db_N > 0, (
                f"thrust sign mismatch at U={c.U_wind_ms} Omega={c.Omega_rpm}")

    def test_ct_envelope_within_3pct(self, beaupoil_survey):
        """Every per-point CT error stays within +/-3 %.  Current
        full-sweep figure is max 1.8 %."""
        offenders = [c for c in beaupoil_survey.comparisons if c.ct_err >= 0.03]
        assert not offenders, (
            "CT err > 3 %:\n" + "\n".join(
                f"  V={c.U_wind_ms} Omega={c.Omega_rpm} "
                f"CT_cc={c.CT_cc:.4f} CT_db={c.CT_db:.4f} "
                f"err={c.ct_err:.1%}" for c in offenders))

    def test_cq_envelope_excluding_near_zero(self, beaupoil_survey):
        """Per-point CQ error stays within +/-15 % at points where
        |CQ_cc| >= 1e-4 (i.e. not at the autorotation crossing where
        denominators go to zero).  Current full-sweep figure is max
        ~5 % away from the crossing."""
        offenders = [c for c in beaupoil_survey.comparisons
                     if abs(c.CQ_cc) >= 1e-4 and c.cq_err >= 0.15]
        assert not offenders, (
            "CQ err > 15 % (excluding |CQ_cc| < 1e-4 near-zero rows):\n" +
            "\n".join(
                f"  V={c.U_wind_ms} Omega={c.Omega_rpm} "
                f"CQ_cc={c.CQ_cc:+.5f} CQ_db={c.CQ_db:+.5f} "
                f"err={c.cq_err:.1%}" for c in offenders))
