"""Off-plane cyclic-trim sweep (quasi-static aero).

Purpose
-------
A tethered rotor that drifts off the vertical wind plane has no direct
"wind-direction sensor".  This sweep tests whether the *trimmed cyclic*
itself can serve as that sensor: hold the rotor at a constant tether
elevation (so altitude is preserved) but rotate it about the vertical by
an azimuth angle phi out of the wind plane, then ask the trim solver for
the cyclic (tilt_lon, tilt_lat) that nulls the steady hub moments at each
phi.  If one cyclic component crosses zero monotonically at phi = 0, it is
an implicit off-plane sensor and a slow restoring law can null it.

Why quasi-static
----------------
QuasiStaticBEM takes hub moments straight from the azimuth-resolved
blade-element integral with the true tangential apparent wind
``v_t_extra(psi)``, so the directional (off-plane) signal is trustworthy
for oblique inflow.  Pitt-Peters applies a wind-axis approximation to the
L-matrix that is only exact for axial / pure-longitudinal flight (see
AGENTS.md), so it is not used here.

Geometry (NED, X North, Y East, Z Down)
----------------------------------------
On-plane tether (hub -> anchor) lives in the Y-Z plane:

    tether_hat(phi=0) = [0, cos(el), sin(el)]

Swinging off-plane by azimuth phi about the vertical, *keeping the
vertical (sin el) component fixed* so altitude is unchanged:

    tether_hat(phi) = [cos(el) * sin(phi),
                       cos(el) * cos(phi),
                       sin(el)]

The hub frame is pinned to the *tether* (not the wind): the longitudinal
axis lies in the tether's vertical plane and the lateral axis is
horizontal and perpendicular to it (``tether_aligned_r_hub``).  This is
the only frame a real vehicle can build -- GPS gives the anchor and the
vehicle position, hence the tether vector; there is no wind sensor.

Findings (quasi-static, beaupoil_2026 rotor, el=30 deg)
-------------------------------------------------------
The off-plane sensor is the longitudinal cyclic ``tilt_lon`` (NOT
``tilt_lat`` -- an earlier guess that was wrong; the tether-pinned frame
puts the signal on the longitudinal axis).

WHY THE NULL IS EXACT (symmetry argument -- the load-bearing result):
On-plane the tether-wind geometry is mirror-symmetric about the tether
vertical plane (the bx-bz plane).  Under that mirror the longitudinal hub
moment My_hub flips sign, so symmetry FORCES My_hub = 0 on-plane, hence
``tilt_lon = 0`` at phi = 0 -- regardless of wind speed, tension, or
collective.  Off-plane breaks the mirror: My_hub ~ sin(phi), so the
longitudinal cyclic that nulls it is an odd, monotonic function of phi
with a guaranteed zero at phi = 0.  That is a calibration-free direction
sensor: sign(tilt_lon) = -sign(phi) (at 5 m/s), so a slow restoring law
``turn_rate_cmd = -k * tilt_lon`` drives the kite back onto the plane and
only the LOOP GAIN (not the zero) needs scheduling.

MATHEMATICAL RELATIONSHIP:
    trim:  My_aero(phi) + (dMy/d tilt_lon) * tilt_lon = 0
           => tilt_lon(phi) = - My_aero(phi) / (dMy/d tilt_lon)
    In the tether frame the horizontal wind splits as
           mu_x ~ V cos(phi)  (longitudinal),  mu_y ~ V sin(phi) (lateral)
    Symmetry makes My_aero odd in mu_y, so My_aero ~ K(mu_x) * mu_y, and
    with cyclic authority g = dMy/d tilt_lon ~ const (set by collective,
    Omega):
           tilt_lon(phi) ~ -(K V / g) sin(phi) ~ c1 * phi  (small phi)

EMPIRICAL FIT over +/-20 deg (least squares vs the swept curve):
    5 m/s : near-plane slope c1 ~ 0.90 deg/deg; tilt_lon ~ -41 deg sin(phi)
            (mild compression at the edges; sin(phi)cos(phi) fits best,
            consistent with My_aero ~ mu_x mu_y ~ sin(phi)cos(phi)).
    7 m/s : essentially linear, c1 ~ 0.43 deg/deg.
   10 m/s : UNRELIABLE -- the 2-D Newton trim diverges once tilt_lat rails
            (see caveats); do not read magnitudes from it.

Inversion for a controller:
    near plane :  phi ~ tilt_lon / c1
    full range :  phi ~ arcsin(tilt_lon / A)   (A ~ -41 deg at 5 m/s)

CAVEATS:
  * The zero of tilt_lon is exact (symmetry); the GAIN c1 / amplitude A is
    operating-point dependent and here DROPS with wind (0.90 -> 0.43
    deg/deg, 5 -> 7 m/s) because the sweep holds collective fixed instead
    of thrust-trimming it, so the cyclic authority g is not constant.
    Schedule c1 on wind/airspeed/collective for an angle estimate; for a
    sign-only restoring law it does not matter.
  * tilt_lat is RAILED at +/-15 deg across the sweep and the trim reports
    conv = n.  That is NOT a sensor problem: My_hub (the sensor axis) is
    nulled cleanly everywhere.  The rail is the rotor's real inherent
    lateral flapping moment Mx_hub in edgewise flow (~238 N*m on-plane,
    roughly even in phi, not symmetry-protected), which exceeds +/-15 deg
    of lateral cyclic authority at this operating point.  It is a
    flight-trim / control-authority issue (can the kite hold roll?),
    separate from sensing.
  * The 10 m/s breakdown is the 2-D Newton interacting badly with the
    railed lateral axis.  A 1-D longitudinal-only trim (target My_hub = 0,
    leave tilt_lat = 0) would make all wind speeds as clean as 5 m/s and
    give a trustworthy c1(V); not yet implemented here.

In this frame the trimmed ``tilt_lon`` is the off-plane sensor: ~0
on-plane, growing ~ sin(phi) off-plane with sign = -sign(phi).

Run
---
    # with uv on PATH:
    uv run python -m envelope.offplane_sweep --winds 5 7 10
    # or call the uv-managed venv directly (Windows PowerShell):
    & ".venv\Scripts\python.exe" -m envelope.offplane_sweep --winds 5 7 10

Output: an ASCII table per wind speed, an .npz next to the script, and
(optionally, with --plot) a PNG if matplotlib is available.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from dynbem import RotorInputs, create_aero, solve_trim_cyclic
from dynbem.rotor_definition import RotorDefinition
from dynbem.rotor_definition import load as load_rotor

from envelope.point_mass import G

_ROTOR_YAML = str(Path(__file__).parent.parent / "rotors" / "beaupoil_2026" / "rotor.yaml")


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def offplane_tether_hat(elevation_deg: float, phi_deg: float) -> np.ndarray:
    """Hub -> anchor unit vector swung off the wind plane by azimuth phi.

    The vertical (Z, sin el) component is held fixed so the rotor keeps
    altitude; only the horizontal projection rotates about the vertical.
    phi = 0 reproduces the in-plane tether ``[0, cos el, sin el]``.
    """
    el = math.radians(elevation_deg)
    phi = math.radians(phi_deg)
    ce, se = math.cos(el), math.sin(el)
    return np.array([ce * math.sin(phi), ce * math.cos(phi), se])


def tether_aligned_r_hub(bz: np.ndarray, t_hat: np.ndarray) -> np.ndarray:
    """Hub frame (hub -> NED) pinned to the tether, not to global North.

    Unlike ``point_mass._build_r_hub`` (which resolves the free azimuth
    DOF with an arbitrary world vector), this pins the hub azimuth to the
    *tether's vertical plane* -- the one quantity a real vehicle can build
    from GPS (it knows the anchor and its own position).  There is NO wind
    sensor, so the frame must never reference the wind.

        bx (longitudinal)  lies in the tether's vertical plane
        by (lateral)       is horizontal, perpendicular to that plane
        bz (thrust)        from the tether+gravity force balance

    In this frame the wind (fixed in the world) is purely longitudinal
    on-plane (phi=0 -> tilt_lat = 0) and acquires a lateral component
    ~ V*sin(phi) off-plane, so the trimmed ``tilt_lat`` is the off-plane
    sensor: zero at phi=0, sign = sign(phi).
    """
    bz = bz / float(np.linalg.norm(bz))
    t_horiz = np.array([t_hat[0], t_hat[1], 0.0])
    nrm = float(np.linalg.norm(t_horiz))
    if nrm < 1e-9:
        # Tether straight down: azimuth undefined; fall back to North.
        t_horiz = np.array([1.0, 0.0, 0.0])
    else:
        t_horiz /= nrm
    bx = t_horiz - float(np.dot(t_horiz, bz)) * bz   # longitudinal, in tether plane
    bx /= float(np.linalg.norm(bx))
    by = np.cross(bz, bx)                            # lateral, perp to tether plane
    by /= float(np.linalg.norm(by))
    return np.column_stack([bx, by, bz])


# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------

def sweep_offplane(
    defn: RotorDefinition,
    *,
    elevation_deg: float = 30.0,
    tension_n: float = 300.0,
    wind_speed_ms: float = 10.0,
    mass_kg: float = 5.0,
    collective_rad: float = 0.10,
    omega_rad_s: float = 40.0,
    phi_deg: Optional[Sequence[float]] = None,
    rho_kg_m3: float = 1.225,
    tilt_limit_deg: float = 15.0,
    tolerance_Nm: float = 0.02,
    model: str = "quasi_static",
) -> dict:
    """Trim cyclic across off-plane azimuth phi at fixed elevation.

    Holds elevation (altitude), tension, wind, collective and omega fixed;
    rotates the tether off the wind plane by each phi and records the
    cyclic (tilt_lon, tilt_lat) that nulls the steady hub moments.

    Returns a dict of equal-length arrays:
        phi_deg, tilt_lon, tilt_lat   [deg, rad, rad]
        Mx_res, My_res                [N*m] trim residuals
        converged                     [bool]
    """
    if phi_deg is None:
        phi_deg = np.arange(-40.0, 40.0 + 1e-9, 5.0)
    phi_deg = np.asarray(phi_deg, dtype=float)

    aero = create_aero(defn, model=model)
    wind_world = np.array([0.0, -wind_speed_ms, 0.0])
    tilt_lim = math.radians(tilt_limit_deg)

    n = phi_deg.size
    tilt_lon = np.full(n, np.nan)
    tilt_lat = np.full(n, np.nan)
    mx_res = np.full(n, np.nan)
    my_res = np.full(n, np.nan)
    converged = np.zeros(n, dtype=bool)

    # Warm-start each trim from the previous phi's solution for continuity
    # and fewer Newton iterations along the curve.
    lon_init, lat_init = 0.0, 0.0

    for i, phi in enumerate(phi_deg):
        t_hat = offplane_tether_hat(elevation_deg, float(phi))
        f_load = tension_n * t_hat + np.array([0.0, 0.0, mass_kg * G])
        bz = f_load / float(np.linalg.norm(f_load))
        R_hub = tether_aligned_r_hub(bz, t_hat)

        base_inputs = RotorInputs(
            collective_rad=collective_rad,
            tilt_lon=0.0,
            tilt_lat=0.0,
            R_hub=R_hub,
            v_hub_world=np.zeros(3),  # station-keeping; moments from wind + spin
            wind_world=wind_world,
            t=0.0,
            omega_rad_s=omega_rad_s,
            rho_kg_m3=rho_kg_m3,
        )

        trim = solve_trim_cyclic(
            aero,
            aero.initial_rotor_state(),
            base_inputs,
            tilt_lon_init=lon_init,
            tilt_lat_init=lat_init,
            tilt_min=-tilt_lim,
            tilt_max=tilt_lim,
            tolerance_Nm=tolerance_Nm,
        )

        tilt_lon[i] = trim.tilt_lon
        tilt_lat[i] = trim.tilt_lat
        mx_res[i] = trim.Mx_residual
        my_res[i] = trim.My_residual
        converged[i] = trim.converged

        if trim.converged:
            lon_init, lat_init = trim.tilt_lon, trim.tilt_lat

    return {
        "phi_deg": phi_deg,
        "tilt_lon": tilt_lon,
        "tilt_lat": tilt_lat,
        "Mx_res": mx_res,
        "My_res": my_res,
        "converged": converged,
    }


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def _print_table(wind_speed_ms: float, res: dict) -> None:
    print()
    print(f"wind = {wind_speed_ms:.1f} m/s")
    print(f"  {'phi[deg]':>9}  {'tilt_lon[deg]':>13}  {'tilt_lat[deg]':>13}  "
          f"{'Mx_res[Nm]':>11}  {'My_res[Nm]':>11}  conv")
    for i in range(res["phi_deg"].size):
        print(f"  {res['phi_deg'][i]:>9.1f}  "
              f"{math.degrees(res['tilt_lon'][i]):>13.4f}  "
              f"{math.degrees(res['tilt_lat'][i]):>13.4f}  "
              f"{res['Mx_res'][i]:>11.4f}  {res['My_res'][i]:>11.4f}  "
              f"{'y' if res['converged'][i] else 'n'}")


def _slope_at_zero(phi_deg: np.ndarray, signal: np.ndarray) -> float:
    """Central-difference slope d(signal)/d(phi) [per deg] near phi = 0."""
    order = np.argsort(phi_deg)
    p = phi_deg[order]
    s = signal[order]
    i0 = int(np.argmin(np.abs(p)))
    if 0 < i0 < p.size - 1:
        return float((s[i0 + 1] - s[i0 - 1]) / (p[i0 + 1] - p[i0 - 1]))
    return float("nan")


def _maybe_plot(results: list[tuple[float, dict]], out_png: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - plotting is optional
        print(f"[plot skipped: matplotlib unavailable: {exc}]")
        return

    fig, (ax_lat, ax_lon) = plt.subplots(1, 2, figsize=(11, 4.2), sharex=True)
    for wind, res in results:
        phi = res["phi_deg"]
        ax_lat.plot(phi, np.degrees(res["tilt_lat"]), marker="o",
                    label=f"{wind:.0f} m/s")
        ax_lon.plot(phi, np.degrees(res["tilt_lon"]), marker="o",
                    label=f"{wind:.0f} m/s")
    for ax, title in ((ax_lat, "tilt_lat (lateral cyclic)"),
                      (ax_lon, "tilt_lon (longitudinal cyclic)")):
        ax.axhline(0.0, color="k", lw=0.6)
        ax.axvline(0.0, color="k", lw=0.6)
        ax.set_xlabel("off-plane azimuth phi [deg]")
        ax.set_ylabel("trimmed cyclic [deg]")
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.legend(title="wind")
    fig.suptitle("Off-plane cyclic trim (quasi-static)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    print(f"saved plot -> {out_png}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rotor", default=_ROTOR_YAML, help="rotor YAML path")
    ap.add_argument("--elevation", type=float, default=30.0,
                    help="tether elevation above horizontal [deg]")
    ap.add_argument("--tension", type=float, default=300.0, help="tether tension [N]")
    ap.add_argument("--mass", type=float, default=5.0, help="vehicle mass [kg]")
    ap.add_argument("--collective", type=float, default=0.10, help="collective pitch [rad]")
    ap.add_argument("--omega", type=float, default=40.0, help="rotor speed [rad/s]")
    ap.add_argument("--winds", type=float, nargs="+", default=[10.0],
                    help="wind speeds to sweep [m/s] (one curve each)")
    ap.add_argument("--phi-min", type=float, default=-40.0, help="min off-plane azimuth [deg]")
    ap.add_argument("--phi-max", type=float, default=40.0, help="max off-plane azimuth [deg]")
    ap.add_argument("--phi-step", type=float, default=5.0, help="off-plane azimuth step [deg]")
    ap.add_argument("--plot", action="store_true", help="save a PNG (needs matplotlib)")
    ap.add_argument("--npz", default=str(Path(__file__).parent / "offplane_sweep.npz"),
                    help="output .npz path")
    args = ap.parse_args(argv)

    defn = load_rotor(args.rotor)
    phi_deg = np.arange(args.phi_min, args.phi_max + 1e-9, args.phi_step)

    results: list[tuple[float, dict]] = []
    save: dict[str, np.ndarray] = {"phi_deg": phi_deg, "winds": np.asarray(args.winds)}

    for wind in args.winds:
        res = sweep_offplane(
            defn,
            elevation_deg=args.elevation,
            tension_n=args.tension,
            wind_speed_ms=wind,
            mass_kg=args.mass,
            collective_rad=args.collective,
            omega_rad_s=args.omega,
            phi_deg=phi_deg,
        )
        results.append((wind, res))
        _print_table(wind, res)
        tag = f"w{wind:g}"
        save[f"tilt_lon_{tag}"] = res["tilt_lon"]
        save[f"tilt_lat_{tag}"] = res["tilt_lat"]

        lat_slope = _slope_at_zero(res["phi_deg"], res["tilt_lat"])
        lon_slope = _slope_at_zero(res["phi_deg"], res["tilt_lon"])
        i0 = int(np.argmin(np.abs(res["phi_deg"])))
        print(f"  near phi=0: tilt_lat = {math.degrees(res['tilt_lat'][i0]):+.4f} deg, "
              f"d(tilt_lat)/d(phi) = {math.degrees(lat_slope):+.5f} deg/deg, "
              f"d(tilt_lon)/d(phi) = {math.degrees(lon_slope):+.5f} deg/deg")

    out_npz = Path(args.npz)
    np.savez(out_npz, **save)
    print(f"\nsaved data -> {out_npz}")

    if args.plot:
        _maybe_plot(results, out_npz.with_suffix(".png"))


if __name__ == "__main__":
    main()
