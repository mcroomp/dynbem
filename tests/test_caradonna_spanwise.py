"""Spanwise sectional CL vs Caradonna & Tung (1981) NASA TM-81232.

The whole-dataset survey lives in oneoff/val_caradonna_spanwise.py
(152 comparisons across 32 operating points).  Summary findings:

    All operating points          median error 27.6%  RMSE 63.7%
    Low-Mach subset (M_tip < 0.5) median error 31.2%  RMSE 73.3%

The BEM systematically over-predicts section CL by 25-50%, consistent
with the same inviscid-incompressible BEM bias seen on Castles-Gray
hover CT.  Outliers (>200%) come from documented OCR damage in the
Caradonna-Tung tables (see CaradonnaTung/CLAUDE.md page index).  At
the rotor tip (r/R = 0.96) the Prandtl tip-loss correction brings the
BEM closer to measurement (~10-15% error), so it's the cleanest single
station for sanity-checking the BEM-vs-data agreement.

The tests below sample three operating points from the cleanest
low-Mach subset:
    Table 17 (theta=8 deg, Omega=1250 rpm, M_tip=0.439) - C-T-marked
        "primary validation case"
    Table 27 (theta=12 deg, Omega=650 rpm,  M_tip=0.226) - lowest M_tip
    Table 28 (theta=12 deg, Omega=1250 rpm, M_tip=0.433)

At each, we assert:
  (1) BEM CL agrees with measured to within +/-100% per station
      (loose - catches factor-of-2 bugs without demanding what the
      Level-1 polar can't deliver),
  (2) BEM tip CL (r/R = 0.96) is within +/-50% (tighter band where
      tip-loss makes the comparison cleanest), and
  (3) the BEM CL profile is monotone-ish from r/R = 0.50 to the peak
      then drops at the tip (qualitative spanwise shape).
"""
from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import pytest

from dynbem.bem import solve_bem_element
from dynbem.rotor_definition import (
    AirfoilProperties, AutorotationProperties, BladeGeometry, RotorDefinition)
from dynbem.polar import LinearPolar


_CSV_DIR = (Path(__file__).parent.parent
            / "Research" / "csv" / "CaradonnaTung")

_ROTOR = RotorDefinition(
    blade=BladeGeometry(
        n_blades=2, radius_m=1.143, root_cutout_m=0.1,
        chord_m=0.1905, twist_deg=0.0, n_elements=30),
    airfoil=AirfoilProperties(
        Re_design=1_000_000, CL0=0.0,
        CL_alpha_per_rad=2 * math.pi, CD0=0.008,
        alpha_stall_deg=15.0, tip_loss=True),
    autorotation=AutorotationProperties(I_ode_kgm2=1.0),
    name="Caradonna-Tung",
)
_POLAR = LinearPolar(_ROTOR.airfoil.CL0, _ROTOR.airfoil.CL_alpha_per_rad,
                     _ROTOR.airfoil.CD0,
                     math.radians(_ROTOR.airfoil.alpha_stall_deg))


# (label, table_num, theta_deg, omega_rpm, M_tip_for_annotation)
_SAMPLES: list[tuple[str, int, float, float, float]] = [
    ("table_17_theta8_1250rpm",  17,  8.0, 1250, 0.439),
    ("table_27_theta12_650rpm",  27, 12.0,  650, 0.226),
    ("table_28_theta12_1250rpm", 28, 12.0, 1250, 0.433),
]


def _section_CL_bem(coll_deg: float, omega_rpm: float, r_over_R: float) -> float:
    R = _ROTOR.blade.radius_m
    omega = omega_rpm * math.pi / 30.0
    Omega_R = omega * R
    r = r_over_R * R
    elem = solve_bem_element(
        r=r, dr=0.005 * R,
        chord=_ROTOR.blade.chord_m, twist_rad=0.0,
        collective_rad=math.radians(coll_deg),
        omega=omega, v_climb=0.0, rho=1.225,
        n_blades=_ROTOR.blade.n_blades, radius_m=R,
        polar=_POLAR, use_tip_loss=_ROTOR.airfoil.tip_loss,
        root_cutout_m=_ROTOR.blade.root_cutout_m,
    )
    v_a = elem.lambda_r * Omega_R
    v_t = omega * r * (1.0 + elem.a_prime)
    phi = math.atan2(v_a, v_t)
    alpha = math.radians(coll_deg) - phi
    cl, _ = _POLAR.cl_cd(alpha)
    return cl


def _load_measured_CL(table_num: int) -> dict[float, float] | None:
    for path in _CSV_DIR.glob(f"page_*_table_{table_num}__cl.csv"):
        with path.open(encoding="ascii") as f:
            reader = csv.reader(f)
            header = next(reader)
            row = next(reader)
            out: dict[float, float] = {}
            for h, v in zip(header, row):
                m = re.match(r"r/R=(\d+\.\d+)", h)
                if m:
                    try:
                        out[float(m.group(1))] = float(v)
                    except ValueError:
                        pass
            return out
    return None


class TestCaradonnaTungSpanwise:
    @pytest.mark.parametrize("sample", _SAMPLES, ids=lambda s: s[0])
    def test_section_cl_within_band(self, sample):
        """At low-Mach C-T operating points, BEM section CL agrees with
        measurement within +/-100% per radial station."""
        _, tbl, theta, omega, _M = sample
        meas = _load_measured_CL(tbl)
        if meas is None:
            pytest.skip(f"no CL csv for C-T table {tbl}")
        for r_over_R, cl_m in meas.items():
            cl_b = _section_CL_bem(theta, omega, r_over_R)
            assert abs(cl_b - cl_m) / abs(cl_m) < 1.0, (
                f"table {tbl}, r/R={r_over_R:.2f}: "
                f"BEM CL={cl_b:.4f}, measured CL={cl_m:.4f}, "
                f"err={abs(cl_b-cl_m)/abs(cl_m):.1%}")

    @pytest.mark.parametrize("sample", _SAMPLES, ids=lambda s: s[0])
    def test_tip_section_cl_tight(self, sample):
        """At the rotor tip (r/R = 0.96), the Prandtl tip-loss correction
        is most influential; the BEM is within +/-50% there in the survey."""
        _, tbl, theta, omega, _M = sample
        meas = _load_measured_CL(tbl)
        if meas is None:
            pytest.skip(f"no CL csv for C-T table {tbl}")
        cl_m = meas.get(0.96)
        if cl_m is None:
            pytest.skip(f"table {tbl} has no r/R=0.96 station")
        cl_b = _section_CL_bem(theta, omega, 0.96)
        assert abs(cl_b - cl_m) / abs(cl_m) < 0.50, (
            f"table {tbl}, r/R=0.96: BEM CL={cl_b:.4f}, "
            f"measured CL={cl_m:.4f}, err={abs(cl_b-cl_m)/abs(cl_m):.1%}")

    @pytest.mark.parametrize("sample", _SAMPLES, ids=lambda s: s[0])
    def test_spanwise_pattern_has_tip_drop(self, sample):
        """The classic BEM spanwise loading rises from root, peaks
        inboard of the tip, and drops at the tip due to Prandtl tip
        loss.  Verify BEM reproduces this qualitative pattern."""
        _, tbl, theta, omega, _M = sample
        stations = [0.50, 0.68, 0.80, 0.89, 0.96]
        cls = [_section_CL_bem(theta, omega, r) for r in stations]
        # Inner-mid increase: CL at r/R=0.80 > CL at r/R=0.50
        assert cls[2] > cls[0], (
            f"table {tbl}: expected CL(r/R=0.80) > CL(r/R=0.50): "
            f"got {cls[2]:.4f} <= {cls[0]:.4f}")
        # Tip drop: CL at tip (0.96) < peak (0.89 or 0.80)
        peak = max(cls[2:4])
        assert cls[4] < peak, (
            f"table {tbl}: expected tip CL(0.96)={cls[4]:.4f} < "
            f"peak={peak:.4f} (Prandtl tip loss)")

    def test_section_cl_csv_loaded(self):
        """Sanity: at least one CL csv exists where expected."""
        for _, tbl, *_ in _SAMPLES:
            assert _load_measured_CL(tbl) is not None, (
                f"missing CL CSV for C-T table {tbl} - "
                f"run Research/extract_tables.py")
