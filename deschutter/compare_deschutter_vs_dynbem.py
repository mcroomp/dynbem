"""Compare De Schutter-style closed-form reference vs dynbem BEM.

This script builds a compact operating-point grid and computes:
  - dynbem CT/CQ from the same setup used in stage 3
  - a De Schutter-style representative-radius closed-form estimate

Outputs:
  - CSV with per-point values and errors
  - printed aggregate error metrics (median, RMSE, max)

Usage:
    c:\repos\aero\.venv\Scripts\python.exe deschutter\compare_deschutter_vs_dynbem.py
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from dynbem import RotorInputs, create_aero
from dynbem.rotor_definition import load as load_rotor
from dynbem.rotor_state import QuasiStaticRotorState


ROOT = Path(__file__).resolve().parent.parent
ROTOR_YAML = ROOT / "deschutter" / "rotor.yaml"
OUT_DIR = ROOT / "deschutter" / "out"
OUT_CSV = OUT_DIR / "deschutter_vs_dynbem.csv"
RHO = 1.225


@dataclass
class Comparison:
    u_wind_ms: float
    omega_rpm: float
    collective_deg: float
    ct_dynbem: float
    cq_dynbem: float
    ct_ref: float
    cq_ref: float

    @property
    def ct_rel_err(self) -> float:
        if abs(self.ct_ref) < 1e-12:
            return float("nan")
        return abs(self.ct_dynbem - self.ct_ref) / abs(self.ct_ref)

    @property
    def cq_rel_err(self) -> float:
        if abs(self.cq_ref) < 1e-12:
            return float("nan")
        return abs(self.cq_dynbem - self.cq_ref) / abs(self.cq_ref)

    @property
    def abs_cq_rel_err(self) -> float:
        denom = abs(self.cq_ref)
        if denom < 1e-12:
            return float("nan")
        return abs(abs(self.cq_dynbem) - abs(self.cq_ref)) / denom


def dynbem_ct_cq(model, u_wind_ms: float, omega_rpm: float, collective_deg: float) -> tuple[float, float]:
    omega = omega_rpm * math.pi / 30.0
    inp = RotorInputs(
        collective_rad=math.radians(collective_deg),
        tilt_lon=0.0,
        tilt_lat=0.0,
        R_hub=np.eye(3),
        v_hub_world=np.zeros(3),
        wind_world=np.array([0.0, 0.0, -u_wind_ms]),
        rho_kg_m3=RHO,
        omega_rad_s=omega,
    )
    result, _ = model.compute_forces(inp, QuasiStaticRotorState())

    r = model.defn.blade.radius_m
    a = math.pi * r * r
    denom_f = RHO * a * (omega * r) ** 2
    denom_q = denom_f * r

    thrust_n = -float(result.F_world[2])
    q_spin_nm = float(result.Q_spin)
    return thrust_n / denom_f, q_spin_nm / denom_q


def deschutter_reference_ct_cq(model, u_wind_ms: float, omega_rpm: float, collective_deg: float) -> tuple[float, float]:
    """Representative-radius De Schutter-style estimate.

    This is intentionally a reduced model for A/B comparison, not a full
    re-implementation. It uses Eq.25-style CL/CD and Eq.30/31-style
    decomposition at one representative radius.
    """
    blade = model.defn.blade
    af = model.defn.airfoil

    n_blades = blade.n_blades
    r_tip = blade.radius_m
    r_root = blade.root_cutout_m
    chord = blade.chord_m
    span = r_tip - r_root
    area_blade = chord * span
    area_disk = math.pi * r_tip * r_tip

    # Two-thirds span from root in the De Schutter summary.
    r_cp = r_root + (2.0 / 3.0) * span

    omega = omega_rpm * math.pi / 30.0
    v_tan = omega * r_cp
    v_ax = u_wind_ms

    # Inflow angle and AoA proxy (alpha = beta + phi in summary doc).
    phi = math.atan2(v_ax, max(v_tan, 1e-9))
    beta = math.radians(collective_deg)
    alpha = beta + phi

    # Keep within the linearized De Schutter validity envelope.
    alpha_lim = math.radians(getattr(af, "alpha_stall_deg", 15.0))
    alpha = float(np.clip(alpha, -alpha_lim, alpha_lim))

    cl_alpha = af.CL_alpha_per_rad
    cl = af.CL0 + cl_alpha * alpha
    cd = af.CD0 + (cl * cl) / (math.pi * 12.0 * 0.8)
    cd_t = 0.021
    cd_total = cd + cd_t

    v_app = math.hypot(v_tan, v_ax)
    q_dyn = 0.5 * RHO * v_app * v_app

    # Per blade forces from Eq.30/31 style decomposition.
    lift = q_dyn * area_blade * cl
    drag = q_dyn * area_blade * cd_total

    thrust_per_blade = lift * math.cos(phi) - drag * math.sin(phi)
    inplane_per_blade = lift * math.sin(phi) + drag * math.cos(phi)

    thrust_total = n_blades * thrust_per_blade

    # dynbem CQ sign is helicopter convention; use negative sign so that
    # wind-driven (turbine-like) torque maps negative, matching dynbem.
    q_spin = -n_blades * r_cp * inplane_per_blade

    denom_f = RHO * area_disk * (omega * r_tip) ** 2
    denom_q = denom_f * r_tip

    ct = thrust_total / denom_f
    cq = q_spin / denom_q
    return ct, cq


def summarize(vals: np.ndarray) -> tuple[float, float, float]:
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return float("nan"), float("nan"), float("nan")
    med = float(np.median(vals))
    rmse = float(np.sqrt(np.mean(vals * vals)))
    maxv = float(np.max(vals))
    return med, rmse, maxv


def main() -> int:
    model = create_aero(load_rotor(str(ROTOR_YAML)), model="bem")

    wind_grid = [6.0, 8.0, 10.0, 12.0]
    omega_grid = [220.0, 270.0, 320.0]
    collective_grid = [-8.0, -6.0, -4.0, -2.0, 0.0]

    comps: list[Comparison] = []
    for u in wind_grid:
        for omega_rpm in omega_grid:
            for collective_deg in collective_grid:
                ct_d, cq_d = dynbem_ct_cq(model, u, omega_rpm, collective_deg)
                ct_r, cq_r = deschutter_reference_ct_cq(model, u, omega_rpm, collective_deg)
                comps.append(
                    Comparison(
                        u_wind_ms=u,
                        omega_rpm=omega_rpm,
                        collective_deg=collective_deg,
                        ct_dynbem=ct_d,
                        cq_dynbem=cq_d,
                        ct_ref=ct_r,
                        cq_ref=cq_r,
                    )
                )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="", encoding="ascii") as f:
        w = csv.writer(f)
        w.writerow([
            "u_wind_ms",
            "omega_rpm",
            "collective_deg",
            "ct_dynbem",
            "cq_dynbem",
            "ct_ref",
            "cq_ref",
            "ct_rel_err",
            "cq_rel_err",
            "abs_cq_rel_err",
        ])
        for c in comps:
            w.writerow([
                f"{c.u_wind_ms:.3f}",
                f"{c.omega_rpm:.3f}",
                f"{c.collective_deg:.3f}",
                f"{c.ct_dynbem:.9e}",
                f"{c.cq_dynbem:.9e}",
                f"{c.ct_ref:.9e}",
                f"{c.cq_ref:.9e}",
                f"{c.ct_rel_err:.9e}",
                f"{c.cq_rel_err:.9e}",
                f"{c.abs_cq_rel_err:.9e}",
            ])

    ct_err = np.array([c.ct_rel_err for c in comps], dtype=float)
    cq_err = np.array([c.cq_rel_err for c in comps], dtype=float)
    abs_cq_err = np.array([c.abs_cq_rel_err for c in comps], dtype=float)

    ct_med, ct_rmse, ct_max = summarize(ct_err)
    cq_med, cq_rmse, cq_max = summarize(cq_err)
    acq_med, acq_rmse, acq_max = summarize(abs_cq_err)

    print(f"Wrote {OUT_CSV}")
    print("Relative error summary vs De Schutter-style reference:")
    print(f"  CT       median={ct_med:.3f}  rmse={ct_rmse:.3f}  max={ct_max:.3f}")
    print(f"  CQ       median={cq_med:.3f}  rmse={cq_rmse:.3f}  max={cq_max:.3f}")
    print(f"  |CQ|     median={acq_med:.3f}  rmse={acq_rmse:.3f}  max={acq_max:.3f}")
    print("Note: this compares dynbem against a reduced representative-radius")
    print("reference, not a full De Schutter DAE + induction solve.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
