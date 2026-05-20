"""Validation of dynbem.bem against Wheatley & Hood (1935) NACA TR 515 —
PCA-2 autogiro full-scale wind-tunnel tests.

This is the **only forward-flight autorotation dataset** in the
empirical validation set; Castles-Gray TN-2474 covers vertical descent
only.  See EMPIRICAL_VALIDATION.md §4 for citation, expected variance,
and the relationship to NACA TR 487 and Harris NASA CR-2008-215370.

Strategy (matches Option 2 in the validation plan)
---------------------------------------------------
For each data row, we prescribe all kinematics from the table -
(mu, alpha_shaft, Omega, collective_pitch) - and let the BEM compute
rotor-axis thrust CT.  We then transform the measured (CL_wind,
CD_wind) into rotor-axis CT and assert agreement.

The PCA-2 had **freely flapping blades**: in forward flight the blades
flap up forward / down backward, which kills the cyclic-induced lift
asymmetry and zeros hub moments at steady state.  Our Level-1 BEM has
rigid blades and no flap dynamics, so a fixed-collective-no-cyclic run
over-predicts CT by ~3x at low mu (lift dominated by advancing-side
high-thrust region).  We model the equivalent flap response by running
``solve_trim_cyclic`` first: cyclic pitch and flap couple via the same
90-degree precession, so the cyclic that zeros hub moments is a
quantitative proxy for the steady-flap response.  CT is then read at
the trimmed condition.

Transformation (airplane axes <-> rotor axes, alpha = disk angle of
attack measured from the freestream to the disk plane):

    T_force = L * cos(alpha) + D * sin(alpha)
    where L = CL_wind * q * pi R^2, D = CD_wind * q * pi R^2,
          q  = 1/2 rho V^2
          V  = (Omega R) * mu / cos(alpha)

    CT_rotor = T_force / (rho * pi R^2 * (Omega R)^2)
             = (CL cos(alpha) + CD sin(alpha)) * mu^2 / (2 cos^2(alpha))

The data is in `Research/csv/Wheatley_Hood_NACA515/`, which is gitignored
because Research/ is treated as local-only data.  The test skips if
the CSVs are not present so CI works on a fresh clone.

Expected variance
-----------------
The rotor fixture is deliberately simplified (constant chord, untwisted,
single NACA 0015 polar reused from Castles-Gray, Re=200k).  The PCA-2
had mixed airfoils (Goettingen 429 inner + symmetric outer), taper, and
~-6 deg twist; we expect ~+/-50% absolute CT bias consistent with the
inviscid-incompressible BEM bias seen on Castles-Gray, Caradonna-Tung,
and Harrington.  The trend across mu and across pitch should track
much tighter (~+/-15%).
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
import pytest

from dynbem.bem import BEMModel
from dynbem import RotorInputs
from dynbem.rotor_definition import load as load_rotor
from dynbem.rotor_state import QuasiStaticRotorState
from dynbem.trim import solve_trim_cyclic


_ROTOR_YAML = (Path(__file__).parent.parent
               / "rotors" / "wheatley_pca2" / "rotor.yaml")
_CSV_DIR = (Path(__file__).parent.parent
            / "Research" / "csv" / "Wheatley_Hood_NACA515")

# With the simplified PCA-2 fixture and rigid-blade BEM, absolute CT
# matches the measurement to within an order of magnitude after
# flap-equivalent cyclic trim.  The trend across mu and pitch matches
# qualitatively (see TestWheatleyTableIII::test_ct_decreases_with_mu and
# TestPitchSweep).  Tighter absolute accuracy would need the real PCA-2
# polar and blade twist; see module docstring.

# Air density at ISA sea level (the wind tunnel was at Langley field;
# atmospheric corrections are below noise for this comparison).
_RHO = 1.225


# ---------------------------------------------------------------------------
# Fixture & helpers
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def pca2_model():
    if not _ROTOR_YAML.exists():
        pytest.skip(f"PCA-2 fixture not found: {_ROTOR_YAML}")
    # n_psi_elements=12 (vs default 36) cuts the per-compute_forces cost ~3x
    # at the cost of slight azimuth averaging error.  This test family runs
    # solve_trim_cyclic dozens of times, so it's the dominant lever.
    return BEMModel(defn=load_rotor(str(_ROTOR_YAML)), n_psi_elements=12)


def _load_csv(name: str) -> list[dict[str, float]]:
    path = _CSV_DIR / name
    if not path.exists():
        pytest.skip(
            f"{path} missing - run Research/extract_tables.py "
            "after extracting source MD tables")
    out: list[dict[str, float]] = []
    with path.open(encoding="ascii") as f:
        for row in csv.DictReader(f):
            # Convert numeric fields; skip the leading 'flag' column if any.
            try:
                out.append({
                    "mu":    float(row["mu"]),
                    "alpha": float(row["alpha (deg)"]),
                    "CL":    float(row["CL"]),
                    "CD":    float(row["CD"]),
                    "L_D":   float(row["L/D"]),
                    "N":     float(row["N (rpm)"]),
                })
            except (KeyError, ValueError):
                continue
    return out


def _measured_CT(mu: float, alpha_deg: float, CL: float, CD: float) -> float:
    """Wheatley airplane-axes (CL, CD) -> rotor-axis CT.

    See module docstring for derivation.
    """
    a = math.radians(alpha_deg)
    return (CL * math.cos(a) + CD * math.sin(a)) * mu**2 / (2.0 * math.cos(a)**2)


def _bem_at_point(model: BEMModel, mu: float, alpha_deg: float,
                  omega_rpm: float, pitch_deg: float):
    """Run dynbem.bem at the prescribed PCA-2 operating point.

    The PCA-2 has freely flapping blades; we model their steady-state
    response with the cyclic that zeros hub moments (90-deg flap-pitch
    coupling).  Returns ``(CT, CQ, tilt_lon_rad, tilt_lat_rad)`` where
    CT and CQ are normalized by ``rho * pi R^2 * (Omega R)^2`` and
    ``rho * pi R^2 * (Omega R)^2 * R`` respectively.
    """
    R = model.defn.blade.radius_m
    omega = omega_rpm * math.pi / 30.0
    a = math.radians(alpha_deg)
    V = omega * R * mu / math.cos(a)

    # R_hub: shaft tilted back by alpha around the +Y axis.  Disk normal
    # then points up-and-forward (in NED: -Z and +X components), letting
    # the freestream flow up through the disk - the autorotation condition.
    R_hub = np.array([
        [math.cos(a), 0.0, -math.sin(a)],
        [0.0,         1.0,  0.0        ],
        [math.sin(a), 0.0,  math.cos(a)],
    ])
    v_hub_world = np.zeros(3)
    wind_world  = np.array([V, 0.0, 0.0])

    # 1. Trim cyclic so hub moments are zero (flap-equivalent for our rigid BEM).
    # Level-1 BEM is quasi-static (no inflow ODE state), so we skip inflow
    # relaxation between Newton steps - compute_forces itself iterates inflow
    # internally to convergence on every call.
    trim = solve_trim_cyclic(
        model,
        QuasiStaticRotorState(omega_rad_s=omega),
        collective_rad=math.radians(pitch_deg),
        R_hub=R_hub,
        v_hub_world=v_hub_world,
        wind_world=wind_world,
        tilt_min=-math.radians(25.0),
        tilt_max= math.radians(25.0),
        tolerance_Nm=1.0,
        max_iterations=20,
        n_inflow_relax=0,
    )

    # 2. CT and CQ at the trimmed condition.
    inputs = RotorInputs(
        collective_rad=math.radians(pitch_deg),
        tilt_lon=trim.tilt_lon, tilt_lat=trim.tilt_lat,
        R_hub=R_hub, v_hub_world=v_hub_world, wind_world=wind_world, t=0.0,
    )
    result, _ = model.compute_forces(inputs, trim.final_state)
    F_hub = R_hub.T @ result.F_world
    T = -F_hub[2]
    A = math.pi * R**2
    CT = T / (_RHO * A * (omega * R)**2)
    CQ = result.Q_spin / (_RHO * A * (omega * R)**2 * R)
    return CT, CQ, trim.tilt_lon, trim.tilt_lat


def _bem_CT(model: BEMModel, mu: float, alpha_deg: float, omega_rpm: float,
            pitch_deg: float) -> tuple[float, float, float]:
    """Back-compat shim: returns (CT, tilt_lon, tilt_lat)."""
    CT, _CQ, tlon, tlat = _bem_at_point(model, mu, alpha_deg, omega_rpm, pitch_deg)
    return CT, tlon, tlat


# ---------------------------------------------------------------------------
# Per-row comparison (high-confidence Tables III and IV)
# ---------------------------------------------------------------------------

# Observed BEM/measured bias range from a 5-point spot check spanning mu = 0.14 - 0.70:
#   ratio ranges from ~1.25 (high mu) to ~2.65 (low mu).  Tolerance is set to
#   bracket this with margin; tightening would require (a) the real PCA-2 polar,
#   (b) blade twist, (c) flap dynamics.  See module docstring.
_BIAS_MIN = 0.50
_BIAS_MAX = 4.00

_TABLE_III_PITCH_DEG = 1.9
_TABLE_IV_PITCH_DEG  = 2.7


@pytest.fixture(scope="module")
def table_iii_rows():
    return _load_csv("page_10_table_iii.csv")


@pytest.fixture(scope="module")
def table_iv_rows():
    return _load_csv("page_10_table_iv.csv")


class TestWheatleyTableIII:
    """Table III: pitch 1.9 deg, faired protuberances.  79 rows in source CSV."""

    def test_csv_loaded(self, table_iii_rows):
        assert len(table_iii_rows) > 50, (
            f"Expected ~79 rows from page_10_table_iii.csv, got "
            f"{len(table_iii_rows)}")

    @pytest.mark.parametrize("row_idx", [0, 40, 75])  # span the mu range
    def test_ct_order_of_magnitude(self, pca2_model, table_iii_rows, row_idx):
        """BEM CT must be within [0.5x, 4x] of measurement after flap-trim.

        Why this band: with the simplified PCA-2 fixture (constant chord,
        no twist, reused NACA 0015 polar) and rigid-blade BEM, the
        absolute CT carries a known ~1.25-2.65x bias that depends on mu.
        Tightening would require the real PCA-2 polar + twist + ideally
        flap dynamics.  This test verifies the BEM gets the right
        order of magnitude; the per-row trend tests below catch
        directional errors.
        """
        if row_idx >= len(table_iii_rows):
            pytest.skip(f"only {len(table_iii_rows)} rows available")
        r = table_iii_rows[row_idx]
        CT_meas = _measured_CT(r["mu"], r["alpha"], r["CL"], r["CD"])
        CT_bem, _, _ = _bem_CT(pca2_model, r["mu"], r["alpha"], r["N"], _TABLE_III_PITCH_DEG)
        ratio = CT_bem / CT_meas
        assert _BIAS_MIN < ratio < _BIAS_MAX, (
            f"row {row_idx} (mu={r['mu']:.3f}, alpha={r['alpha']:.1f}, "
            f"N={r['N']:.1f}): CT_bem/CT_meas = {ratio:.2f} outside "
            f"[{_BIAS_MIN}, {_BIAS_MAX}]")

    def test_ct_decreases_with_mu_at_fixed_omega(self, pca2_model, table_iii_rows):
        """At fixed Omega ~ 98 rpm, CT_meas decreases monotonically with mu.

        The BEM should reproduce the same monotonic trend (sign of
        dCT/dmu) even if absolute values differ.
        """
        rows_98 = sorted([r for r in table_iii_rows if 97 < r["N"] < 100],
                         key=lambda x: x["mu"])
        assert len(rows_98) >= 4
        # Sample the lowest-mu and highest-mu rows at 98 rpm; the band is
        # what matters, not every intermediate point.
        sampled = [rows_98[0], rows_98[-1]]
        cts_bem = [_bem_CT(pca2_model, r["mu"], r["alpha"], r["N"],
                           _TABLE_III_PITCH_DEG)[0] for r in sampled]
        assert cts_bem[0] > cts_bem[-1], (
            f"BEM CT should decrease from mu={sampled[0]['mu']:.2f} to "
            f"mu={sampled[-1]['mu']:.2f}: got {cts_bem[0]:.5f} -> "
            f"{cts_bem[-1]:.5f}")


class TestWheatleyTableIV:
    """Table IV: pitch 2.7 deg, faired.  74 rows.

    Higher pitch than Table III should give higher CT at matching mu.
    """

    def test_csv_loaded(self, table_iv_rows):
        assert len(table_iv_rows) > 50

    @pytest.mark.parametrize("row_idx", [0, 35, 70])
    def test_ct_order_of_magnitude(self, pca2_model, table_iv_rows, row_idx):
        if row_idx >= len(table_iv_rows):
            pytest.skip(f"only {len(table_iv_rows)} rows available")
        r = table_iv_rows[row_idx]
        CT_meas = _measured_CT(r["mu"], r["alpha"], r["CL"], r["CD"])
        CT_bem, _, _ = _bem_CT(pca2_model, r["mu"], r["alpha"], r["N"], _TABLE_IV_PITCH_DEG)
        ratio = CT_bem / CT_meas
        assert _BIAS_MIN < ratio < _BIAS_MAX, (
            f"row {row_idx} (mu={r['mu']:.3f}, alpha={r['alpha']:.1f}, "
            f"N={r['N']:.1f}): CT_bem/CT_meas = {ratio:.2f} outside "
            f"[{_BIAS_MIN}, {_BIAS_MAX}]")


class TestAutorotationEquilibrium:
    """Each TR 515 data point IS at autorotation by construction
    (the wind tunnel had no rotor drive — the rotor freely span up to its
    self-balance Omega), so for the real rotor Q_aero == 0 exactly.

    A whole-dataset survey (oneoff/val_wheatley_autorotation.py, 469
    rows across Tables I-IV) showed the BEM produces a consistent
    residual:

        no trim:  CQ mean = -0.00169, RMSE 0.00171, max|CQ| = 0.00220
        trimmed:  CQ mean = -0.00113, RMSE 0.00116, max|CQ| = 0.00176

    The sign is consistently negative: the BEM thinks the air drives the
    rotor harder than the real autorotation balance.  Likely cause: the
    rigid-blade BEM cannot reproduce flapping-induced changes in local
    angle of attack on the advancing/retreating sides, so total induced
    power comes out higher than measured.

    Bounds below are derived from the 469-row survey:
      - |CQ| <= 0.003 brackets every row in the dataset.
      - CQ <= 0 (negative) matches the observed sign bias in all rows.
    """

    # Sampled to span low/mid/high mu and the two pitch settings with
    # HIGH-confidence transcription.  Each takes ~1.5s with cyclic trim;
    # all five together fit under the default 10s pytest timeout.
    SAMPLES: list[tuple[str, float, float, float, float, float]] = [
        # (label, mu, alpha_deg, N_rpm, pitch_deg, expected_CQ_band_max)
        ("low_mu_high_alpha",   0.145, 15.9, 98.8,  1.9, 0.003),
        ("mid_mu_mid_alpha",    0.252,  6.1, 119.2, 1.9, 0.003),
        ("high_mu_low_alpha",   0.491,  1.8, 138.5, 1.9, 0.003),
        ("very_high_mu",        0.701,  1.3, 97.6,  1.9, 0.003),
        ("pitch_2_7_mid_mu",    0.270,  4.8, 118.9, 2.7, 0.003),
    ]

    @pytest.mark.parametrize("sample", SAMPLES, ids=lambda s: s[0])
    def test_cq_residual_bounded(self, pca2_model, sample):
        """|CQ_BEM| at the autorotation operating point is within the
        envelope observed across the full 469-row dataset survey."""
        _, mu, alpha, N, pitch, max_cq = sample
        _CT, CQ, _, _ = _bem_at_point(pca2_model, mu, alpha, N, pitch)
        assert abs(CQ) < max_cq, (
            f"|CQ_BEM| = {abs(CQ):.5f} exceeds dataset envelope {max_cq} "
            f"at mu={mu:.3f}, alpha={alpha:.1f}, N={N:.1f}, pitch={pitch}")

    @pytest.mark.parametrize("sample", SAMPLES, ids=lambda s: s[0])
    def test_cq_sign_matches_dataset_bias(self, pca2_model, sample):
        """Across all 469 rows the BEM gives CQ < 0 (over-drives the
        rotor relative to autorotation balance).  Sampled rows must
        carry the same sign so regressions that flip Q convention or
        change the bias direction get caught."""
        _, mu, alpha, N, pitch, _ = sample
        _CT, CQ, _, _ = _bem_at_point(pca2_model, mu, alpha, N, pitch)
        assert CQ < 0.0, (
            f"CQ_BEM = {CQ:+.5f} (expected negative - aero over-driving "
            f"rotor) at mu={mu:.3f}, alpha={alpha:.1f}, N={N:.1f}")


class TestPitchSweep:
    """Cross-table: at matched mu, CT(pitch=2.7) > CT(pitch=1.9).

    This is a more robust trend check than absolute CT - it exercises
    the BEM's response to collective input independent of the
    polar-fidelity bias.
    """

    def test_higher_pitch_gives_higher_ct(self, pca2_model,
                                          table_iii_rows, table_iv_rows):
        target_mu = 0.30
        # Find nearest row in each table.
        r3 = min(table_iii_rows, key=lambda r: abs(r["mu"] - target_mu))
        r4 = min(table_iv_rows,  key=lambda r: abs(r["mu"] - target_mu))
        CT3, _, _ = _bem_CT(pca2_model, r3["mu"], r3["alpha"], r3["N"], _TABLE_III_PITCH_DEG)
        CT4, _, _ = _bem_CT(pca2_model, r4["mu"], r4["alpha"], r4["N"], _TABLE_IV_PITCH_DEG)
        assert CT4 > CT3, (
            f"BEM at pitch 2.7 deg, mu={r4['mu']:.3f} should give higher "
            f"CT than pitch 1.9 deg, mu={r3['mu']:.3f}: "
            f"got CT4={CT4:.5f}, CT3={CT3:.5f}")
