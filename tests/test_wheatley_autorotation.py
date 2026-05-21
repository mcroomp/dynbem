"""Validation of dynbem.bem against Wheatley & Hood (1935) NACA TR 515 --
PCA-2 autogiro full-scale wind-tunnel tests.

This test drives the same comparison loop that
verification/wheatley_hood_autorotation_torque.py uses for its 469-row
whole-dataset sweep -- just with sample=2 per table (8 rows total) so
it fits in the pytest budget.  Keep the BEM-call / trim / measured-CT
logic in the verification module; this file only asserts on the
returned Comparison records.

This is the **only forward-flight autorotation dataset** in the
empirical validation set; Castles-Gray TN-2474 covers vertical descent
only.  See EMPIRICAL_VALIDATION.md section 4 for citation, expected
variance, and the relationship to NACA TR 487 and Harris CR-2008-215370.

Whole-dataset survey numbers (from running the verifier with no sample,
469 rows across Tables I-IV):

    no trim:  CQ mean = -0.00169, RMSE 0.00171, max|CQ| = 0.00220
    trimmed:  CQ mean = -0.00113, RMSE 0.00116, max|CQ| = 0.00176

The sign is consistently negative across all 469 rows: the BEM thinks
the air drives the rotor harder than the real autorotation balance.
Likely cause: the rigid-blade BEM cannot reproduce flapping-induced
reductions in local angle of attack on the advancing side, so
integrated induced power comes out higher than measured.

CT bias (BEM / measured): ratio ranges from ~1.25 (high mu) to ~2.65
(low mu).  Tightening would need the real PCA-2 polar
(Harris Appendix 11.9) and ideally flap dynamics.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from verification.wheatley_hood_autorotation_torque import (  # noqa: E402
    Comparison, evaluate_point, load_model, load_table, run_survey,
)


# Per-table sample size used for the aggregate fixture.  2 evenly-spaced
# rows per table x 4 tables = 8 comparisons covering the full mu range
# at all three pitch settings.  At ~0.8s per trim+notrim pair this fits
# comfortably under the 10s pytest timeout (compute the survey once,
# module-scope).
_SAMPLE_N = 2

# Order-of-magnitude bound on CT_bem / CT_meas; matches the
# 1.25x - 2.65x bias observed in the whole-dataset survey with margin
# for the simplified fixture (constant chord, no twist, NACA 0015 polar).
_BIAS_MIN = 0.50
_BIAS_MAX = 4.00

# Envelope across the 469-row trimmed survey; sampled rows must stay
# inside it.
_CQ_TRIM_ENVELOPE = 0.003


@pytest.fixture(scope="module")
def pca2_model():
    model = load_model()
    if model is None:
        pytest.skip("PCA-2 fixture not found")
    return model


@pytest.fixture(scope="module")
def survey(pca2_model):
    """Sampled survey: 2 evenly-spaced rows per table, trim only."""
    s = run_survey(pca2_model, sample=_SAMPLE_N, include_notrim=False)
    if not s.comparisons:
        pytest.skip("no Wheatley CSV data available")
    return s


class TestWheatleyCT:
    def test_sample_loaded(self, survey):
        """Sanity: sampled survey actually populated comparisons."""
        labels = {c.table_label for c in survey.comparisons}
        assert "table_iii" in labels and "table_iv" in labels, (
            f"sampled survey should cover HIGH-confidence tables III + IV; "
            f"got {sorted(labels)}")

    def test_ct_ratio_within_bias_band(self, survey):
        """Every sampled BEM CT must land within the documented
        order-of-magnitude band relative to the airplane-axes
        measurement transformed to rotor axis."""
        offenders: list[Comparison] = []
        for c in survey.comparisons:
            ratio = c.ct_ratio
            if not math.isfinite(ratio):
                continue
            if not (_BIAS_MIN < ratio < _BIAS_MAX):
                offenders.append(c)
        assert not offenders, (
            f"CT_bem/CT_meas outside [{_BIAS_MIN}, {_BIAS_MAX}]:\n" +
            "\n".join(
                f"  {c.table_label} mu={c.mu:.3f} alpha={c.alpha_deg:.1f} "
                f"N={c.N_rpm:.1f}: ratio={c.ct_ratio:.2f}"
                for c in offenders))


class TestAutorotationEquilibrium:
    """Each TR 515 data point IS at autorotation by construction
    (the wind tunnel had no rotor drive), so for the real rotor
    Q_aero == 0 exactly.  Sampled BEM rows must fall inside the
    469-row trimmed envelope and carry the dataset-wide negative
    sign bias."""

    def test_cq_residual_bounded(self, survey):
        offenders = [c for c in survey.comparisons
                     if abs(c.CQ_trim) >= _CQ_TRIM_ENVELOPE]
        assert not offenders, (
            f"|CQ_BEM| >= {_CQ_TRIM_ENVELOPE}:\n" +
            "\n".join(
                f"  {c.table_label} mu={c.mu:.3f} alpha={c.alpha_deg:.1f} "
                f"N={c.N_rpm:.1f}: CQ={c.CQ_trim:+.5f}"
                for c in offenders))

    def test_cq_sign_matches_dataset_bias(self, survey):
        """Across all 469 rows the BEM gives CQ < 0.  Sampled rows must
        carry the same sign so regressions that flip Q convention or
        change the bias direction get caught."""
        offenders = [c for c in survey.comparisons if c.CQ_trim >= 0.0]
        assert not offenders, (
            "expected CQ_BEM < 0 (aero over-driving rotor):\n" +
            "\n".join(
                f"  {c.table_label} mu={c.mu:.3f} alpha={c.alpha_deg:.1f}: "
                f"CQ={c.CQ_trim:+.5f}"
                for c in offenders))


class TestTrends:
    """Trend assertions need specific row pairs, not aggregate samples.
    Call evaluate_point directly with the chosen rows.  All BEM-driving
    logic still lives in the verification module."""

    def test_ct_decreases_with_mu_at_fixed_omega(self, pca2_model):
        """Table III at fixed N ~ 98 rpm: CT_meas decreases with mu; BEM
        must reproduce the same monotonic trend even if absolute values
        differ."""
        rows = load_table("table_iii")
        if not rows:
            pytest.skip("table III csv missing")
        rows_98 = sorted([r for r in rows if 97 < r["N"] < 100],
                         key=lambda x: x["mu"])
        assert len(rows_98) >= 4
        low = evaluate_point(pca2_model, rows_98[0],  1.9, "table_iii",
                             include_notrim=False)
        high = evaluate_point(pca2_model, rows_98[-1], 1.9, "table_iii",
                              include_notrim=False)
        assert low.CT_bem > high.CT_bem, (
            f"BEM CT should decrease from mu={low.mu:.2f} to "
            f"mu={high.mu:.2f}: got {low.CT_bem:.5f} -> {high.CT_bem:.5f}")

    def test_higher_pitch_gives_higher_ct(self, pca2_model):
        """At matched mu, CT(pitch=2.7) > CT(pitch=1.9)."""
        rows3 = load_table("table_iii")
        rows4 = load_table("table_iv")
        if not rows3 or not rows4:
            pytest.skip("table III or IV csv missing")
        target_mu = 0.30
        r3 = min(rows3, key=lambda r: abs(r["mu"] - target_mu))
        r4 = min(rows4, key=lambda r: abs(r["mu"] - target_mu))
        c3 = evaluate_point(pca2_model, r3, 1.9, "table_iii",
                            include_notrim=False)
        c4 = evaluate_point(pca2_model, r4, 2.7, "table_iv",
                            include_notrim=False)
        assert c4.CT_bem > c3.CT_bem, (
            f"BEM at pitch 2.7 deg, mu={c4.mu:.3f} should give higher CT "
            f"than pitch 1.9 deg, mu={c3.mu:.3f}: "
            f"got CT4={c4.CT_bem:.5f}, CT3={c3.CT_bem:.5f}")
