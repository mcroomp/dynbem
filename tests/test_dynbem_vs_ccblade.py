"""Cross-check: dynbem.bem vs CCBlade on the Beaupoil RAWES rotor.

The verification module
[`verification/dynbem_vs_ccblade_beaupoil.py`](../verification/dynbem_vs_ccblade_beaupoil.py)
drives the full comparison; this test runs it in sampled mode.

Status -- KNOWN-OPEN INVESTIGATION
----------------------------------
At present dynbem and CCBlade disagree substantially on the Beaupoil
rotor (after sign-convention reconciliation):

    median CT error ~40 %, max ~500 %
    median |CQ| error ~250 %, max ~19000 %  (very-low-disk-loading rows)

The disagreement is largest at low U_wind / high Omega (high TSR,
light disk loading) and tightest at the high U_wind / low Omega
corner (heavy disk loading).  Candidate root causes being investigated
include the Prandtl tip-loss formulation at light loading, polar
interpretation near zero AoA, and 4-blade solidity arithmetic.

Until the disagreement is narrowed down, this test does NOT assert
a tight numerical envelope -- doing so would just bake in whichever
side is currently wrong.  Instead it exercises the integration
end-to-end:

  - the CCBlade cached CSV exists and parses,
  - the rotor.yaml fixture loads,
  - dynbem runs cleanly at every sampled operating point,
  - signs are sane (both BEMs predict positive thrust),
  - results are finite.

When the disagreement narrows after debugging, replace the loose
checks with a real bound (median CT err < X, etc.) sourced from
running the verifier with no sample.
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
    CCBLADE_CSV, ROTOR_YAML, run_survey,
)

_SAMPLE_N = 5


@pytest.fixture(scope="module")
def survey():
    if not CCBLADE_CSV.exists():
        pytest.skip(
            f"missing CCBlade reference CSV: {CCBLADE_CSV}.  "
            f"Build it via verification/ccblade_docker/ (see README)."
        )
    if not ROTOR_YAML.exists():
        pytest.skip(f"missing rotor fixture: {ROTOR_YAML}")
    return run_survey(sample=_SAMPLE_N)


class TestDynbemVsCCBladeBeaupoil:
    def test_sample_populated(self, survey):
        """Sanity: the sampled survey produced records."""
        assert len(survey.comparisons) >= 3, (
            f"expected >=3 comparisons at sample={_SAMPLE_N}, "
            f"got {len(survey.comparisons)}")

    def test_results_are_finite(self, survey):
        """No NaN / Inf -- catches blow-ups in the BEM solver."""
        for c in survey.comparisons:
            for name, v in (
                ("CT_db", c.CT_db), ("CQ_db", c.CQ_db),
                ("T_db_N", c.T_db_N), ("Q_db_Nm", c.Q_db_Nm),
            ):
                assert math.isfinite(v), (
                    f"non-finite {name}={v} at U={c.U_wind_ms} "
                    f"Omega={c.Omega_rpm} pitch={c.pitch_deg}")

    def test_thrust_sign_agrees(self, survey):
        """Both BEMs predict a positive thrust (rotor reacts upwind into
        the wind) on Beaupoil at every operating point in the sweep."""
        for c in survey.comparisons:
            assert c.T_cc_N > 0 and c.T_db_N > 0, (
                f"thrust sign mismatch at U={c.U_wind_ms} Omega={c.Omega_rpm}: "
                f"CCBlade T={c.T_cc_N:.1f}N, dynbem T={c.T_db_N:.1f}N")

    def test_ct_within_order_of_magnitude(self, survey):
        """Very loose: dynbem CT must be within a factor of 10 of CCBlade.
        A factor-of-100 disagreement would indicate a fundamental
        misconfiguration (sign flipped, wrong solidity, etc.); a
        factor of 2-5 is the current state of the investigation."""
        for c in survey.comparisons:
            ratio = c.CT_db / c.CT_cc if c.CT_cc != 0 else float("inf")
            assert 0.1 < ratio < 10.0, (
                f"CT ratio {ratio:.2f} outside [0.1, 10] at U={c.U_wind_ms} "
                f"Omega={c.Omega_rpm}: CCBlade={c.CT_cc:.4f}, dynbem={c.CT_db:.4f}")
