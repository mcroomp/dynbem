"""Deskew scanned pages by aligning long ruled lines to 0/90 degrees.

Detects long near-horizontal and near-vertical line segments via the
probabilistic Hough transform, computes the median angular deviation
from the nearest axis, and rotates the image to cancel it.

Usage:
    python deskew.py page_10.png [page_11.png ...]            # in-place, .bak written
    python deskew.py --out-dir straight/ page_*.png           # write copies
    python deskew.py --dry-run page_10.png                    # report angle only
    python deskew.py --no-backup page_10.png                  # overwrite, no backup
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import cv2
import numpy as np


def estimate_skew_deg(
    img: np.ndarray,
    min_line_frac: float = 0.25,
    angle_tol_deg: float = 8.0,
) -> tuple[float, int]:
    """Return (skew_deg, n_lines_used). Positive skew = image rotated CCW from upright.

    `min_line_frac` is the minimum line length as a fraction of min(H, W).
    `angle_tol_deg` is how far from 0/90 a segment may be and still count.
    """
    if img.ndim == 3:
        gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    else:
        gray = img

    h, w = gray.shape
    min_len = int(min_line_frac * min(h, w))

    edges = cv2.Canny(gray, 50, 150, apertureSize=3)
    segments = cv2.HoughLinesP(
        edges,
        rho=1,
        theta=np.pi / 1800,
        threshold=120,
        minLineLength=min_len,
        maxLineGap=10,
    )
    if segments is None:
        return 0.0, 0

    deviations: list[float] = []
    weights: list[float] = []
    for x1, y1, x2, y2 in segments[:, 0, :]:
        dx, dy = x2 - x1, y2 - y1
        length = float(np.hypot(dx, dy))
        ang = np.degrees(np.arctan2(dy, dx))
        # Fold to nearest axis: map angle to deviation from 0 or 90.
        # Bring into (-90, 90]:
        ang = ((ang + 90.0) % 180.0) - 90.0
        # Distance to nearest of {-90, 0, 90}:
        if abs(ang) <= angle_tol_deg:
            dev = ang  # near horizontal
        elif abs(abs(ang) - 90.0) <= angle_tol_deg:
            dev = ang - 90.0 if ang > 0 else ang + 90.0  # near vertical
        else:
            continue
        deviations.append(dev)
        weights.append(length)

    if not deviations:
        return 0.0, 0

    devs = np.array(deviations)
    wts = np.array(weights)
    # Weighted median.
    order = np.argsort(devs)
    devs_s, wts_s = devs[order], wts[order]
    cum = np.cumsum(wts_s)
    cutoff = cum[-1] / 2.0
    skew = float(devs_s[np.searchsorted(cum, cutoff)])
    return skew, len(deviations)


def rotate_image(img: np.ndarray, angle_deg: float) -> np.ndarray:
    """Rotate `img` by `angle_deg` (CCW) about its centre, expanding the canvas
    and filling new pixels with white."""
    h, w = img.shape[:2]
    cx, cy = w / 2.0, h / 2.0
    M = cv2.getRotationMatrix2D((cx, cy), angle_deg, 1.0)
    cos = abs(M[0, 0])
    sin = abs(M[0, 1])
    new_w = int(h * sin + w * cos)
    new_h = int(h * cos + w * sin)
    M[0, 2] += new_w / 2.0 - cx
    M[1, 2] += new_h / 2.0 - cy
    border = (255, 255, 255) if img.ndim == 3 else 255
    return cv2.warpAffine(
        img,
        M,
        (new_w, new_h),
        flags=cv2.INTER_CUBIC,
        borderMode=cv2.BORDER_CONSTANT,
        borderValue=border,
    )


def deskew_file(
    src: Path,
    dst: Path,
    *,
    backup: bool,
    dry_run: bool,
    min_skew_deg: float = 0.05,
) -> None:
    img = cv2.imdecode(np.fromfile(str(src), dtype=np.uint8), cv2.IMREAD_UNCHANGED)
    if img is None:
        print(f"  {src.name}: could not read", file=sys.stderr)
        return

    skew, n = estimate_skew_deg(img)
    print(f"  {src.name}: skew = {skew:+.3f} deg from {n} lines", end="")

    if dry_run:
        print(" (dry-run)")
        return
    if abs(skew) < min_skew_deg:
        print(" (below threshold, skipped)")
        return

    rotated = rotate_image(img, skew)

    if backup and src == dst and not src.with_suffix(src.suffix + ".bak").exists():
        shutil.copy2(src, src.with_suffix(src.suffix + ".bak"))

    dst.parent.mkdir(parents=True, exist_ok=True)
    ok, buf = cv2.imencode(dst.suffix, rotated)
    if not ok:
        print(" -> encode failed", file=sys.stderr)
        return
    buf.tofile(str(dst))
    print(f" -> wrote {dst.name}")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    p.add_argument("inputs", nargs="+", type=Path, help="image files to deskew")
    p.add_argument("--out-dir", type=Path, default=None,
                   help="write straightened copies here instead of in place")
    p.add_argument("--no-backup", action="store_true",
                   help="when writing in place, do not save a .bak")
    p.add_argument("--dry-run", action="store_true",
                   help="report detected skew but write nothing")
    p.add_argument("--min-skew", type=float, default=0.05,
                   help="below this many degrees the image is left alone")
    args = p.parse_args(argv)

    print(f"Deskewing {len(args.inputs)} file(s):")
    for src in args.inputs:
        if not src.is_file():
            print(f"  {src}: not a file", file=sys.stderr)
            continue
        dst = (args.out_dir / src.name) if args.out_dir else src
        deskew_file(
            src,
            dst,
            backup=not args.no_backup,
            dry_run=args.dry_run,
            min_skew_deg=args.min_skew,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
