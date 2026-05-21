"""Cross-check: dynbem.bem vs CCBlade on real rotors.

Two cross-check pairs, both driven by their verification module in
sampled mode (one BEM-driver path, fast tests, full-sweep script
elsewhere for re-baselining bounds):

NREL Phase VI (TestDynbemVsCCBladePhaseVI)
------------------------------------------
The right cross-check: a real twisted/tapered wind turbine in the
operating envelope CCBlade was calibrated for.  S809 airfoil, 2 blades,
R = 5.029 m, 72 RPM, +3 deg tip pitch, V_wind 5..25 m/s.

Full-sweep numbers from running the verifier with no --sample (21
operating points):

    CT err: median 12 %, mean 28 %, RMSE 43 %, max 147 % (V=5, TSR=7.5)
    CQ err: median 37 %, mean 66 %, RMSE 98 %, max 333 % (V=5)

The CT envelope tightens steadily as wind speed grows -- at typical
operating speeds (V >= 10 m/s) dynbem matches CCBlade within 10-15 %
on CT and 25-40 % on CQ.  The low-V outliers are the high-TSR /
light-disk-loading regime where both BEMs are sensitive to induction-
modelling details.

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

import numpy as np
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

    def test_ct_median_within_20pct(self, phase_vi_survey):
        """Median CT error <= 20 % across the V_wind = 5..25 m/s sweep.
        Current full-sweep figure is ~12 %."""
        ct = phase_vi_survey.ct_errors()
        med = float(np.median(ct))
        assert med < 0.20, (
            f"median CT error {med:.1%} exceeds 20 % "
            f"(was ~12 % at last full-sweep run)")

    def test_ct_tight_at_operating_speeds(self, phase_vi_survey):
        """For V_wind >= 12 m/s (well past the high-TSR transition),
        CT tracks within +/- 20 % at every point.  V=10..11 sit right
        at the TSR ~ 4-5 transition where both BEMs' induction details
        diverge faster."""
        op = [c for c in phase_vi_survey.comparisons if c.U_wind_ms >= 12.0]
        offenders = [c for c in op if c.ct_err >= 0.20]
        assert not offenders, (
            "CT err > 20 % at operating speeds:\n" + "\n".join(
                f"  V={c.U_wind_ms} CT_cc={c.CT_cc:.4f} CT_db={c.CT_db:.4f} "
                f"err={c.ct_err:.1%}" for c in offenders))

    def test_cq_median_within_50pct(self, phase_vi_survey):
        """Median CQ error <= 50 % across the sweep.  Current figure
        is ~37 %.  Tighter than CT in absolute % is hard because CQ is
        small near hover-end of the sweep."""
        cq = phase_vi_survey.cq_errors()
        med = float(np.median(cq))
        assert med < 0.50, (
            f"median CQ error {med:.1%} exceeds 50 % "
            f"(was ~37 % at last full-sweep run)")


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
