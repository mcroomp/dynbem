"""Isolate where De Schutter-style reference diverges from dynbem.

Reads deschutter/out/deschutter_vs_dynbem.csv and emits:
- grouped error summaries by wind, rpm, collective
- top-N worst points for CT and |CQ|
- sign mismatch rates and correlations

Outputs are printed and also saved to:
  deschutter/out/divergence_report.txt

Usage:
    c:\repos\aero\.venv\Scripts\python.exe deschutter\isolate_divergence.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
IN_CSV = ROOT / "deschutter" / "out" / "deschutter_vs_dynbem.csv"
OUT_TXT = ROOT / "deschutter" / "out" / "divergence_report.txt"


def _table_group(df: pd.DataFrame, group_col: str, err_col: str) -> pd.DataFrame:
    g = (
        df.groupby(group_col, as_index=False)[err_col]
        .agg(median="median", mean="mean", max="max")
        .sort_values("mean", ascending=False)
    )
    return g


def _fmt(df: pd.DataFrame, float_cols: list[str]) -> str:
    d = df.copy()
    for c in float_cols:
        if c in d.columns:
            d[c] = d[c].map(lambda x: f"{x:.6f}")
    return d.to_string(index=False)


def main() -> int:
    if not IN_CSV.exists():
        raise FileNotFoundError(f"Missing input CSV: {IN_CSV}")

    df = pd.read_csv(IN_CSV)

    top_ct = (
        df.sort_values("ct_rel_err", ascending=False)
        .head(10)[["u_wind_ms", "omega_rpm", "collective_deg", "ct_dynbem", "ct_ref", "ct_rel_err"]]
    )
    top_cq = (
        df.sort_values("abs_cq_rel_err", ascending=False)
        .head(10)[["u_wind_ms", "omega_rpm", "collective_deg", "cq_dynbem", "cq_ref", "abs_cq_rel_err"]]
    )

    ct_by_wind = _table_group(df, "u_wind_ms", "ct_rel_err")
    ct_by_rpm = _table_group(df, "omega_rpm", "ct_rel_err")
    ct_by_collective = _table_group(df, "collective_deg", "ct_rel_err")

    acq_by_wind = _table_group(df, "u_wind_ms", "abs_cq_rel_err")
    acq_by_rpm = _table_group(df, "omega_rpm", "abs_cq_rel_err")
    acq_by_collective = _table_group(df, "collective_deg", "abs_cq_rel_err")

    sign_ct = float((np.sign(df["ct_dynbem"]) != np.sign(df["ct_ref"])).mean())
    sign_cq = float((np.sign(df["cq_dynbem"]) != np.sign(df["cq_ref"])).mean())

    corr_ct = float(df["ct_dynbem"].corr(df["ct_ref"]))
    corr_cq = float(df["cq_dynbem"].corr(df["cq_ref"]))
    corr_abs_cq = float(np.corrcoef(np.abs(df["cq_dynbem"]), np.abs(df["cq_ref"]))[0, 1])

    lines: list[str] = []
    lines.append("Divergence isolation report")
    lines.append(f"rows={len(df)}")
    lines.append("")

    lines.append("Top 10 CT relative error points")
    lines.append(_fmt(top_ct, ["u_wind_ms", "omega_rpm", "collective_deg", "ct_dynbem", "ct_ref", "ct_rel_err"]))
    lines.append("")

    lines.append("Top 10 |CQ| relative error points")
    lines.append(_fmt(top_cq, ["u_wind_ms", "omega_rpm", "collective_deg", "cq_dynbem", "cq_ref", "abs_cq_rel_err"]))
    lines.append("")

    lines.append("CT error by wind")
    lines.append(_fmt(ct_by_wind, ["u_wind_ms", "median", "mean", "max"]))
    lines.append("")

    lines.append("CT error by rpm")
    lines.append(_fmt(ct_by_rpm, ["omega_rpm", "median", "mean", "max"]))
    lines.append("")

    lines.append("CT error by collective")
    lines.append(_fmt(ct_by_collective, ["collective_deg", "median", "mean", "max"]))
    lines.append("")

    lines.append("|CQ| error by wind")
    lines.append(_fmt(acq_by_wind, ["u_wind_ms", "median", "mean", "max"]))
    lines.append("")

    lines.append("|CQ| error by rpm")
    lines.append(_fmt(acq_by_rpm, ["omega_rpm", "median", "mean", "max"]))
    lines.append("")

    lines.append("|CQ| error by collective")
    lines.append(_fmt(acq_by_collective, ["collective_deg", "median", "mean", "max"]))
    lines.append("")

    lines.append("Global diagnostics")
    lines.append(f"CT sign mismatch rate={sign_ct:.6f}")
    lines.append(f"CQ sign mismatch rate={sign_cq:.6f}")
    lines.append(f"CT correlation={corr_ct:.6f}")
    lines.append(f"CQ correlation={corr_cq:.6f}")
    lines.append(f"|CQ| correlation={corr_abs_cq:.6f}")

    report = "\n".join(lines)
    OUT_TXT.parent.mkdir(parents=True, exist_ok=True)
    OUT_TXT.write_text(report, encoding="ascii")

    print(report)
    print()
    print(f"Wrote {OUT_TXT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
