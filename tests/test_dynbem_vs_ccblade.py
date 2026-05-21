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
Kept around as a known-open investigation, not a passing target:
Beaupoil is a helicopter-style rotor (untwisted, constant chord,
cambered SG6040, zero pitch) being driven by wind -- outside the
design envelope of standard turbine BEMs.  dynbem puts every blade
station into stall at the design point; CCBlade does not.  Smoke
test only (finiteness + thrust-sign + factor-of-10 CT envelope).
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

_SAMPLE_N = 5


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
    return run_beaupoil_survey(sample=_SAMPLE_N)


class TestDynbemVsCCBladeBeaupoil:
    def test_sample_populated(self, beaupoil_survey):
        assert len(beaupoil_survey.comparisons) >= 3

    def test_results_are_finite(self, beaupoil_survey):
        for c in beaupoil_survey.comparisons:
            for name, v in (("CT_db", c.CT_db), ("CQ_db", c.CQ_db),
                            ("T_db_N", c.T_db_N), ("Q_db_Nm", c.Q_db_Nm)):
                assert math.isfinite(v), (
                    f"non-finite {name}={v} at U={c.U_wind_ms}")

    def test_thrust_sign_agrees(self, beaupoil_survey):
        for c in beaupoil_survey.comparisons:
            assert c.T_cc_N > 0 and c.T_db_N > 0, (
                f"thrust sign mismatch at U={c.U_wind_ms} Omega={c.Omega_rpm}")

    def test_ct_within_order_of_magnitude(self, beaupoil_survey):
        """Beaupoil cross-check is mismatched-by-design (helicopter rotor
        forced into wind-turbine flow).  Loose factor-of-10 bound only."""
        for c in beaupoil_survey.comparisons:
            ratio = c.CT_db / c.CT_cc if c.CT_cc != 0 else float("inf")
            assert 0.1 < ratio < 10.0, (
                f"CT ratio {ratio:.2f} outside [0.1, 10] at "
                f"U={c.U_wind_ms} Omega={c.Omega_rpm}")
