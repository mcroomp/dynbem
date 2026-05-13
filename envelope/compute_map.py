"""Stable-regime envelope map for a tethered rotor.

Sweeps (v_target × wind × elevation) combinations, ramping tether tension
continuously from t_min to t_max and recording collective angle at every
sample_dn N along the ramp.  Output is a smooth curve of collective vs
tension for each (elevation, wind, v_target) combination.

Usage (CLI)
-----------
    python -m envelope.compute_map --quick --save quick.npz
    python -m envelope.compute_map --full  --save full.npz --plot c:/temp/plots
    python -m envelope.compute_map --load  full.npz --plot c:/temp/plots

Output arrays (shape: [n_vtargets, n_winds, n_elevations, n_samples])
-----------------------------------------------------------------------
    cols_arr    equilibrium collective (rad) at each sampled tension
    v_alongs_arr  actual v_along (m/s) at each sample
    sats_arr    bool — was collective clamped at that sample?
    tensions_arr  1-D tension axis (N), shared across all curves
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np

# ---------------------------------------------------------------------------
# Grid presets
# ---------------------------------------------------------------------------

_ROTOR_YAML = str(
    Path(__file__).parent.parent / "rotors" / "beaupoil_2026" / "rotor.yaml"
)

QUICK_GRID: dict[str, Any] = {
    "v_targets":   [-0.5, 0.5],
    "winds":       [10.0],
    "elevations":  [30.0, 40.0, 50.0, 60.0, 70.0, 80.0],
    "t_min":       100.0,
    "t_max":       1000.0,
    "sample_dn":   1.0,          # N between recorded samples
    "mass_kg":     5.0,
    "omega_init":  20.0,
    "settle_time": 20.0,
    "ramp_rate":   0.5,
    "dt":          0.02,
    "kp_col":      0.01,
    "ki_col":      0.02,
    "col_min":     -0.25,
    "col_max":      0.20,
    "n_workers":   4,
    "rotor_yaml":  _ROTOR_YAML,
}

FULL_GRID: dict[str, Any] = {
    "v_targets":   [-0.5, 0.5],
    "winds":       [10.0],
    "elevations":  [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0],
    "t_min":       100.0,
    "t_max":       1000.0,
    "sample_dn":   1.0,
    "mass_kg":     5.0,
    "omega_init":  20.0,
    "settle_time": 20.0,
    "ramp_rate":   0.5,
    "dt":          0.02,
    "kp_col":      0.01,
    "ki_col":      0.02,
    "col_min":     -0.25,
    "col_max":      0.20,
    "n_workers":   os.cpu_count() or 4,
    "rotor_yaml":  _ROTOR_YAML,
}


# ---------------------------------------------------------------------------
# Grid computation
# ---------------------------------------------------------------------------

def compute_grid(params: dict[str, Any] | None = None, **kwargs: Any) -> dict[str, Any]:
    """Compute the envelope curves.

    Returns a dict with 4-D arrays [n_vt, n_winds, n_el, n_samples] and
    a shared 1-D tensions_arr axis.
    """
    cfg = dict(FULL_GRID)
    if params is not None:
        cfg.update(params)
    cfg.update(kwargs)

    from aero.rotor_definition import load as _load_rotor
    defn = _load_rotor(cfg["rotor_yaml"])

    v_targets:  list[float] = list(cfg["v_targets"])
    winds:      list[float] = list(cfg["winds"])
    elevations: list[float] = list(cfg["elevations"])
    t_min:      float = float(cfg["t_min"])
    t_max:      float = float(cfg["t_max"])
    sample_dn:  float = float(cfg["sample_dn"])
    mass_kg:    float = float(cfg["mass_kg"])
    n_workers:  int   = int(cfg.get("n_workers", 4))

    nv  = len(v_targets)
    nw  = len(winds)
    na  = len(elevations)
    n_samples = int(round((t_max - t_min) / sample_dn)) + 1
    tensions_arr = np.linspace(t_min, t_max, n_samples)

    shape = (nv, nw, na, n_samples)
    cols_arr     = np.full(shape, np.nan)
    v_alongs_arr = np.full(shape, np.nan)
    omegas_arr   = np.full(shape, np.nan)
    tilts_arr    = np.full(shape, np.nan)
    sats_arr     = np.zeros(shape, dtype=bool)

    project_root = str(Path(__file__).parent.parent)
    worker_jobs: list[tuple[tuple[int, int, int], dict]] = []

    for vi, vt in enumerate(v_targets):
        for wi, ws in enumerate(winds):
            for ai, el in enumerate(elevations):
                worker_jobs.append(((vi, wi, ai), {
                    "defn":          defn,
                    "mass_kg":       mass_kg,
                    "v_target":      vt,
                    "wind_speed":    ws,
                    "elevation_deg": el,
                    "t_min":         t_min,
                    "t_max":         t_max,
                    "sample_dn":     sample_dn,
                    "omega_init":    float(cfg.get("omega_init", 20.0)),
                    "settle_time":   float(cfg.get("settle_time", 20.0)),
                    "ramp_rate":     float(cfg.get("ramp_rate", 0.5)),
                    "dt":            float(cfg.get("dt", 0.02)),
                    "kp_col":        float(cfg.get("kp_col", 0.01)),
                    "ki_col":        float(cfg.get("ki_col", 0.02)),
                    "col_min":       float(cfg.get("col_min", -0.25)),
                    "col_max":       float(cfg.get("col_max", 0.20)),
                    "project_root":  project_root,
                }))

    total = len(worker_jobs)
    print(f"Grid: {nv} v_targets x {nw} winds x {na} elevations  "
          f"({n_samples} samples per curve, {t_min:.0f}-{t_max:.0f} N)")
    print(f"Launching {total} jobs across {n_workers} workers ...")

    from envelope.point_mass import ramp_column_worker

    t0 = time.time()
    done = 0

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        future_to_idx = {
            pool.submit(ramp_column_worker, wargs): idx
            for idx, wargs in worker_jobs
        }
        for future in as_completed(future_to_idx):
            vi, wi, ai = future_to_idx[future]
            done += 1
            try:
                res = future.result()
            except Exception as exc:
                print(f"  [{done}/{total}] el={elevations[ai]:.0f}deg "
                      f"w={winds[wi]} vt={v_targets[vi]:+.1f} -- ERROR: {exc}")
                continue

            # Align returned arrays onto the shared tension axis
            ret_t = res["tensions"]
            ret_c = res["cols"]
            ret_v = res["v_alongs"]
            ret_s = res["sats"]
            ret_o = res["omegas"]
            ret_tl = res["tilts"]
            for si, t_s in enumerate(tensions_arr):
                idx_near = int(np.argmin(np.abs(ret_t - t_s)))
                if abs(ret_t[idx_near] - t_s) < sample_dn:
                    cols_arr    [vi, wi, ai, si] = ret_c[idx_near]
                    v_alongs_arr[vi, wi, ai, si] = ret_v[idx_near]
                    sats_arr    [vi, wi, ai, si] = ret_s[idx_near]
                    omegas_arr  [vi, wi, ai, si] = ret_o[idx_near]
                    tilts_arr   [vi, wi, ai, si] = ret_tl[idx_near]

            elapsed = time.time() - t0
            rate = done / elapsed
            eta = (total - done) / rate if rate > 0 else 0
            print(f"  [{done}/{total}] el={elevations[ai]:.0f}deg "
                  f"w={winds[wi]:.0f} vt={v_targets[vi]:+.1f}  ETA {eta:.0f}s")

    print(f"Done in {time.time() - t0:.1f}s.")

    return {
        "cols_arr":     cols_arr,
        "v_alongs_arr": v_alongs_arr,
        "omegas_arr":   omegas_arr,
        "tilts_arr":    tilts_arr,
        "sats_arr":     sats_arr,
        "tensions_arr": tensions_arr,
        "v_targets":    np.array(v_targets),
        "winds":        np.array(winds),
        "elevations":   np.array(elevations),
        "mass_kg":      mass_kg,
        "rotor_name":   defn.name,
    }


# ---------------------------------------------------------------------------
# Save / load
# ---------------------------------------------------------------------------

def save_grid(data: dict[str, Any], path: str) -> None:
    np.savez_compressed(
        path,
        cols_arr=data["cols_arr"],
        v_alongs_arr=data["v_alongs_arr"],
        omegas_arr=data["omegas_arr"],
        tilts_arr=data["tilts_arr"],
        sats_arr=data["sats_arr"],
        tensions_arr=data["tensions_arr"],
        v_targets=data["v_targets"],
        winds=data["winds"],
        elevations=data["elevations"],
        mass_kg=np.array([data["mass_kg"]]),
        rotor_name=np.array([data["rotor_name"]]),
    )
    print(f"Saved -> {path}")


def load_grid(path: str) -> dict[str, Any]:
    raw = np.load(path, allow_pickle=True)
    return {
        "cols_arr":     raw["cols_arr"],
        "v_alongs_arr": raw["v_alongs_arr"],
        "omegas_arr":   raw["omegas_arr"],
        "tilts_arr":    raw["tilts_arr"],
        "sats_arr":     raw["sats_arr"],
        "tensions_arr": raw["tensions_arr"],
        "v_targets":    raw["v_targets"],
        "winds":        raw["winds"],
        "elevations":   raw["elevations"],
        "mass_kg":      float(raw["mass_kg"][0]),
        "rotor_name":   str(raw["rotor_name"][0]),
    }


# ---------------------------------------------------------------------------
# PNG plots  — one figure per wind speed, curves per elevation
# ---------------------------------------------------------------------------

def plot_curves(data: dict[str, Any], out_dir: str = ".") -> list[str]:
    """Save one PNG per wind speed: collective (deg) vs tension (N).

    Each PNG has one subplot per v_target.  Elevation is encoded by colour
    (viridis).  Saturated regions are shown as a thicker translucent overlay.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
    except ImportError:
        print("matplotlib not available -- skipping PNG output")
        return []

    tensions    = data["tensions_arr"]
    cols_deg    = np.degrees(data["cols_arr"])   # [nv, nw, na, ns]
    sats        = data["sats_arr"]
    v_targets   = data["v_targets"]
    winds       = data["winds"]
    elevations  = data["elevations"]
    col_min_deg = np.degrees(float(data.get("col_min", -0.25)))
    col_max_deg = np.degrees(float(data.get("col_max",  0.20)))
    rotor_name  = data.get("rotor_name", "")
    mass_kg     = data.get("mass_kg", "")

    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    na = len(elevations)
    cmap = matplotlib.colormaps.get_cmap("viridis").resampled(na)
    colors = [cmap(i) for i in range(na)]

    norm = mcolors.Normalize(vmin=elevations[0], vmax=elevations[-1])
    sm = cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])

    saved: list[str] = []
    nv = len(v_targets)

    for wi, ws in enumerate(winds):
        fig, axes = plt.subplots(1, nv, figsize=(7 * nv, 5), sharey=True)
        if nv == 1:
            axes = [axes]

        fig.suptitle(
            f"Equilibrium collective vs tether tension  --  {rotor_name}  "
            f"mass={mass_kg} kg  wind={ws:.0f} m/s",
            fontsize=11,
        )

        for vi, vt in enumerate(v_targets):
            ax = axes[vi]
            label = "reel-out" if vt < 0 else "reel-in"

            for ei, el in enumerate(elevations):
                c_deg = cols_deg[vi, wi, ei]      # (n_samples,)
                sat   = sats    [vi, wi, ei]

                ax.plot(tensions, c_deg, color=colors[ei], linewidth=1.5,
                        label=f"{el:.0f} deg")

                # Overlay saturated segments as thick translucent band
                if sat.any():
                    c_sat = np.where(sat, c_deg, np.nan)
                    ax.plot(tensions, c_sat, color=colors[ei],
                            linewidth=5, alpha=0.35, solid_capstyle="round")

            # Collective limits
            ax.axhline(col_min_deg, color="0.4", linewidth=0.8,
                       linestyle="--", label=f"col_min ({col_min_deg:.1f} deg)")
            ax.axhline(col_max_deg, color="0.4", linewidth=0.8,
                       linestyle=":",  label=f"col_max ({col_max_deg:.1f} deg)")

            ax.set_xlabel("Tether tension (N)")
            ax.set_title(f"v_target = {vt:+.1f} m/s ({label})")
            ax.grid(True, alpha=0.3)

        axes[0].set_ylabel("Equilibrium collective (deg)")

        cb = fig.colorbar(sm, ax=axes, label="Elevation (deg)", shrink=0.8)

        fig.tight_layout()
        fname = out_path / f"envelope_wind{ws:.0f}ms.png"
        fig.savefig(fname, dpi=150)
        plt.close(fig)
        print(f"Saved -> {fname}")
        saved.append(str(fname))

    return saved


# ---------------------------------------------------------------------------
# Heatmap plot  — sampled from the continuous curves every heatmap_dn N
# ---------------------------------------------------------------------------

def plot_grid(data: dict[str, Any], out_dir: str = ".",
              heatmap_dn: float = 10.0,
              t_max_plot: float | None = None,
              quantity: str = "col") -> list[str]:
    """Save one PNG per wind speed as a heatmap.

    quantity  "col"  — equilibrium collective (deg), diverging RdBu_r
              "rpm"  — rotor speed (RPM), sequential plasma
              "tilt" — rotor tilt from vertical (deg), sequential viridis
    Cells beyond the first aero-clamp or PI-saturation are left blank.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
    except ImportError:
        print("matplotlib not available -- skipping PNG output")
        return []

    tensions_full = data["tensions_arr"]
    sats_full     = data["sats_arr"]
    v_targets     = data["v_targets"]
    winds         = data["winds"]
    elevations    = data["elevations"]
    rotor_name    = data.get("rotor_name", "")
    mass_kg       = data.get("mass_kg", "")

    # Select quantity to display
    if quantity == "rpm":
        vals_full = data["omegas_arr"] * (60.0 / (2.0 * np.pi))
        qty_label = "Rotor speed (RPM)"
        fmt = "{:.0f}"
    elif quantity == "tilt":
        vals_full = data["tilts_arr"]
        qty_label = "Rotor tilt from vertical (deg)"
        fmt = "{:.1f}"
    else:
        vals_full = np.degrees(data["cols_arr"])
        qty_label = "Collective (deg)"
        fmt = "{:+.1f}"

    t_min = tensions_full[0]
    t_max = tensions_full[-1] if t_max_plot is None else min(t_max_plot, tensions_full[-1])
    t_grid = np.arange(t_min, t_max + heatmap_dn * 0.5, heatmap_dn)

    def _idx(t: float) -> int:
        return int(np.argmin(np.abs(tensions_full - t)))

    t_indices = [_idx(t) for t in t_grid]
    vals_hm = vals_full[:, :, :, t_indices]   # [nv, nw, na, nt_grid]
    sats_hm = sats_full[:, :, :, t_indices]

    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    finite = vals_hm[np.isfinite(vals_hm)]
    if quantity == "col":
        vmax = max(float(np.nanmax(np.abs(finite))), 1.0)
        norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
        cmap = "RdBu_r"
    else:
        vmin = float(np.nanmin(finite)) if len(finite) else 0.0
        vmax = float(np.nanmax(finite)) if len(finite) else 1.0
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        cmap = "plasma" if quantity == "rpm" else "viridis"

    na = len(elevations)
    nt = len(t_grid)
    cell_w = max(0.5, 5.0 / na)
    cell_h = max(0.2, 4.0 / nt)
    fontsize = max(5, min(8, int(min(cell_w, cell_h) * 13)))

    el_edges = _cell_edges(np.array(elevations, dtype=float))
    t_edges  = _cell_edges(t_grid)

    saved: list[str] = []
    nv = len(v_targets)

    for wi, ws in enumerate(winds):
        fig_w = (na * cell_w + 2.0) * nv
        fig_h = nt * cell_h + 1.8
        fig, axes = plt.subplots(1, nv, figsize=(fig_w, fig_h), squeeze=False)
        fig.suptitle(
            f"{qty_label} -- {rotor_name}  mass={mass_kg} kg  wind={ws:.0f} m/s",
            fontsize=11,
        )

        # Tilt overlay: shown as a second line in each cell for collective plot
        tilt_hm = data["tilts_arr"][:, :, :, t_indices] if quantity == "col" else None

        for vi, vt in enumerate(v_targets):
            ax = axes[0, vi]
            label = "reel-out" if vt < 0 else "reel-in"

            z = vals_hm[vi, wi, :, :].copy()    # (na, nt_grid)
            sat = sats_hm[vi, wi, :, :]

            # Blank out everything from the first saturated sample onward
            for ai in range(na):
                first_sat = np.argmax(sat[ai])
                if sat[ai, first_sat]:
                    z[ai, first_sat:] = np.nan

            mesh = ax.pcolormesh(
                el_edges, t_edges, z.T,
                cmap=cmap, norm=norm, shading="flat",
            )

            el_ctrs = np.array(elevations, dtype=float)
            for ai in range(na):
                for ti in range(nt):
                    xc = el_ctrs[ai]
                    yc = t_grid[ti]
                    val = z[ai, ti]
                    if not np.isfinite(val):
                        continue
                    brightness = (val - norm.vmin) / max(norm.vmax - norm.vmin, 1e-9)
                    tc = "white" if brightness > 0.6 else "black"
                    if tilt_hm is not None:
                        tilt_val = tilt_hm[vi, wi, ai, ti]
                        cell_label = f"{fmt.format(val)}\n{tilt_val:.0f}°"
                    else:
                        cell_label = fmt.format(val)
                    ax.text(xc, yc, cell_label, ha="center", va="center",
                            fontsize=fontsize, color=tc, fontweight="bold",
                            linespacing=1.1)

            ax.set_xlabel("Tether elevation (deg)")
            ax.set_ylabel("Tether tension (N)")
            ax.set_title(f"v_target={vt:+.1f} m/s ({label})")
            ax.set_xlim(el_edges[0], el_edges[-1])
            ax.set_ylim(t_edges[0], t_edges[-1])
            ax.set_xticks(el_ctrs)
            ax.set_yticks(t_grid)
            fig.colorbar(mesh, ax=ax, label=qty_label)

        fig.tight_layout()
        fname = out_path / f"envelope_{quantity}_wind{ws:.0f}ms.png"
        fig.savefig(fname, dpi=150)
        plt.close(fig)
        print(f"Saved -> {fname}")
        saved.append(str(fname))

    return saved


def _cell_edges(centres: np.ndarray) -> np.ndarray:
    c = np.asarray(centres, dtype=float)
    mids = 0.5 * (c[:-1] + c[1:])
    return np.concatenate([[c[0] - (mids[0] - c[0])], mids, [c[-1] + (c[-1] - mids[-1])]])


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tethered-rotor envelope map")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--quick", action="store_true", help="Coarse grid (fast preview)")
    mode.add_argument("--full",  action="store_true", help="Full grid")
    mode.add_argument("--load",  metavar="FILE",      help="Load existing .npz and re-plot")
    parser.add_argument("--save",    metavar="FILE", help="Save result to .npz")
    parser.add_argument("--plot",    metavar="DIR",  help="Save PNG plots to DIR")
    parser.add_argument("--rotor",   metavar="YAML", help="Override rotor YAML path")
    parser.add_argument("--mass",    type=float,     help="Override vehicle mass (kg)")
    parser.add_argument("--workers", type=int,       help="Number of parallel workers")
    parser.add_argument("--tmax",     type=float,     help="Max tension to plot (N)")
    parser.add_argument("--quantity", default="col",
                        choices=["col", "rpm", "tilt"],
                        help="Quantity to plot: col (collective), rpm, tilt")
    args = parser.parse_args()

    if args.load:
        data = load_grid(args.load)
        plot_dir = args.plot or str(Path(args.load).parent)
        plot_grid(data, plot_dir, t_max_plot=args.tmax, quantity=args.quantity)
        sys.exit(0)

    params = dict(QUICK_GRID if args.quick else FULL_GRID)
    if args.rotor:
        params["rotor_yaml"] = args.rotor
    if args.mass:
        params["mass_kg"] = args.mass
    if args.workers:
        params["n_workers"] = args.workers

    data = compute_grid(params)

    if args.save:
        save_grid(data, args.save)

    plot_dir = args.plot or (str(Path(args.save).parent) if args.save else ".")
    plot_grid(data, plot_dir, t_max_plot=args.tmax, quantity=args.quantity)
