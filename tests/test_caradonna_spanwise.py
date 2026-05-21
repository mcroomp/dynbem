"""Spanwise sectional CL vs Caradonna & Tung (1981) NASA TM-81232.

This test drives the same comparison loop that
verification/caradonna_tung_spanwise_cl.py uses for its whole-dataset
sweep -- just with `sample=` set small so it fits in the pytest budget.
Keep the BEM-call logic in the verification module; this file only
asserts on the returned aggregate.

Whole-dataset survey numbers (from running the verifier with no
sample, 151 comparisons across 32 operating points):

    All operating points          median error 30.8%  RMSE 82.7%
    Low-Mach subset (M_tip < 0.5) median error 31.5%  RMSE 89.6%

The BEM systematically over-predicts section CL by 25-50%, consistent
with the same inviscid-incompressible BEM bias seen on Castles-Gray
hover CT. Outliers (>200%) come from documented OCR damage in some
Caradonna-Tung tables (see CaradonnaTung/CLAUDE.md). At the rotor tip
(r/R = 0.96) the Prandtl tip-loss correction brings BEM within ~10-20%
of measurement, so that's the cleanest single station.

Sample size N=8 (every 4th operating point) takes well under 1s and
exercises the full theta/rpm/M_tip range. Per-point bound is
+/-150% to absorb the OCR-damaged outliers that the wider sample
sees; comparisons with |CL_meas| <= 0.01 are skipped (the verifier
itself filters these out of the median/RMSE stats).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

# Pull the verification module in as the single source of truth for
# how a section-CL comparison is run.  No BEM call logic here.
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from verification.caradonna_tung_spanwise_cl import run_survey  # noqa: E402


# Eight evenly-spaced points across the 32-table dataset -- captures
# the M_tip = 0.226 - 0.890 range and all three theta steps.
_SAMPLE_N = 8


@pytest.fixture(scope="module")
def survey():
    return run_survey(sample=_SAMPLE_N)


class TestCaradonnaTungSpanwise:
    def test_some_data_loaded(self, survey):
        """Sanity: the sampled survey actually ran and loaded CL data."""
        assert survey.points_run >= 5, (
            f"sample of {_SAMPLE_N} should yield >=5 tables with CL data, "
            f"got {survey.points_run}")
        assert len(survey.comparisons) >= 20, (
            f"expected >=20 per-station comparisons, got {len(survey.comparisons)}")

    def test_section_cl_within_band(self, survey):
        """Every per-station BEM CL is within +/-150% of measurement
        (comparisons with |CL_meas| <= 0.01 are skipped: no signal to
        normalise against, matches the survey's own filtering).

        Catches factor-of-2 BEM bugs without demanding what the
        Level-1 polar can't deliver across the documented OCR-damaged
        outliers in some tables.
        """
        offenders = [c for c in survey.comparisons
                     if math.isfinite(c.err) and c.err >= 1.50]
        assert not offenders, (
            "per-station err > 150%:\n" + "\n".join(
                f"  table {c.table_num} r/R={c.r_over_R:.2f}: "
                f"BEM={c.cl_bem:.4f} meas={c.cl_meas:.4f} err={c.err:.1%}"
                for c in offenders))

    def test_tip_station_cl_tight(self, survey):
        """At r/R = 0.96 the Prandtl tip-loss correction is most
        influential; the BEM is within +/-50% there across the survey.
        Skip near-zero-CL points (no signal)."""
        tip = [c for c in survey.comparisons
               if abs(c.r_over_R - 0.96) < 1e-6 and math.isfinite(c.err)]
        assert tip, "sample should contain at least one tip-station comparison with signal"
        offenders = [c for c in tip if c.err >= 0.50]
        assert not offenders, (
            "tip-station err > 50%:\n" + "\n".join(
                f"  table {c.table_num}: BEM={c.cl_bem:.4f} "
                f"meas={c.cl_meas:.4f} err={c.err:.1%}"
                for c in offenders))

    def test_aggregate_median_bounded(self, survey):
        """Sampled-sweep median absolute error stays below 80%.

        Full-sweep median is ~31% (docstring). The 80% bound gives
        headroom for sample-to-sample variance without making the test
        a no-op."""
        arr = survey.errors()
        med = float(sorted(arr)[len(arr) // 2])
        assert med < 0.80, f"median error {med:.1%} exceeds 80%"
