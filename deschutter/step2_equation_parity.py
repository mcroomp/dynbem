"""Step 2: Equation-level parity checks for De Schutter core formulas.

This script validates the closed-form coefficients used in the 2018 model
summary against known numeric reference values.

Usage:
    uv run python deschutter/step2_equation_parity.py
"""
from __future__ import annotations

import math


def cl_alpha_eq25(aspect_ratio: float) -> float:
    return (2.0 * math.pi) / (1.0 + 2.0 / aspect_ratio)


def cl_cd_eq25(alpha_rad: float, aspect_ratio: float, oswald_eff: float, cd0: float) -> tuple[float, float]:
    cl_alpha = cl_alpha_eq25(aspect_ratio)
    cl = cl_alpha * alpha_rad
    cd = cd0 + (cl * cl) / (math.pi * aspect_ratio * oswald_eff)
    return cl, cd


def cd_t_eq29(cyl_cd: float, cable_d_m: float, cable_l_m: float, wing_area_m2: float) -> float:
    return cyl_cd * cable_d_m * cable_l_m / wing_area_m2


def assert_close(name: str, value: float, expected: float, tol: float) -> None:
    if abs(value - expected) > tol:
        raise AssertionError(f"{name} failed: value={value:.10f}, expected={expected:.10f}, tol={tol:.3e}")


def main() -> int:
    # Eq.25 parameters from De Schutter style summary.
    ar = 12.0
    oswald = 0.8
    cd0 = 0.01

    # 1) Lift slope parity.
    cl_alpha = cl_alpha_eq25(ar)
    assert_close("CL_alpha", cl_alpha, 12.0 * math.pi / 7.0, 1e-12)

    # 2) Example point alpha=5 deg parity.
    cl_5, cd_5 = cl_cd_eq25(math.radians(5.0), ar, oswald, cd0)
    assert_close("CL(alpha=5deg)", cl_5, 0.4699799208, 1e-6)
    assert_close("CD(alpha=5deg)", cd_5, 0.0173193645, 1e-6)

    # 3) Eq.29 structural drag coefficient parity.
    # Using d_cable=1.5 mm, L_cable=2.6 m, S_w=0.1875 m^2.
    cd_t = cd_t_eq29(1.0, 0.0015, 2.6, 0.1875)
    assert_close("CD_T", cd_t, 0.0208, 1e-4)

    print("Step 2 parity checks passed.")
    print(f"CL_alpha={cl_alpha:.9f}  CL@5deg={cl_5:.9f}  CD@5deg={cd_5:.9f}  CD_T={cd_t:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
