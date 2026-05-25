"""Smoke test for the flight-envelope sweep (envelope/compute_map.py).

Runs a tiny 1-elevation x 1-v_target x 1-wind grid with a short settle
time and coarse tension range so it completes quickly under pytest.
Asserts:
  - No job crashes (cols_arr row is not all-NaN).
  - Collective values land inside the PI limits.
  - omega_arr values are in a physically reasonable band.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from envelope.compute_map import compute_grid  # noqa: E402

_SAMPLED_PARAMS = {
    "v_targets":   [0.5],
    "winds":       [10.0],
    "elevations":  [60.0],
    "t_min":       150.0,
    "t_max":       400.0,
    "sample_dn":   50.0,        # 6 tension samples -- fast
    "mass_kg":     5.0,
    "omega_init":  20.0,
    "settle_time": 5.0,         # shortened from 20 s
    "ramp_rate":   2.0,         # faster ramp
    "dt":          0.02,        # coarser timestep acceptable for smoke test
    "kp_col":      0.01,
    "ki_col":      0.02,
    "col_min":     -0.25,
    "col_max":      0.20,
    "n_workers":   1,           # single worker avoids multiprocessing overhead
    "model":       "pitt_peters",
}


def test_envelope_smoke():
    """Envelope runs without errors and produces physically plausible output."""
    data = compute_grid(_SAMPLED_PARAMS)

    cols    = data["cols_arr"]       # shape (1, 1, 1, n_samples)
    omegas  = data["omegas_arr"]
    sats    = data["sats_arr"]

    # At least one non-NaN sample per curve.
    assert np.any(np.isfinite(cols)), "all collective samples are NaN -- job crashed"

    finite_cols = cols[np.isfinite(cols)]
    assert np.all(finite_cols >= -0.30), f"collective below floor: {finite_cols.min():.3f}"
    assert np.all(finite_cols <=  0.25), f"collective above ceiling: {finite_cols.max():.3f}"

    finite_omegas = omegas[np.isfinite(omegas)]
    assert np.all(finite_omegas > 0), "omega not positive"
    assert np.all(finite_omegas < 350), f"omega unrealistically high: {finite_omegas.max():.1f}"
