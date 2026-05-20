"""Walk Research/ and convert every markdown table into a CSV file.

Use this whenever a `.md` table extraction is updated to keep the
parallel `.csv` siblings in sync.  Idempotent — safe to re-run.

Usage
-----
    python Research/extract_tables.py            # process all of Research/
    python Research/extract_tables.py CaradonnaTung Castles_TN2474
                                                 # process only those subfolders

CSV layout
----------
All CSVs land under `Research/csv/`, mirroring the source paper folder
structure.  For `Research/<Paper>/foo.md` containing markdown tables:
- 1 table   → `Research/csv/<Paper>/foo.csv`
- N tables  → `Research/csv/<Paper>/foo__<slug>.csv` per table, where
              `<slug>` comes from the nearest preceding `##`/`###`
              heading (or `table_01`, ... if no heading was found).

Rules for what counts as a markdown table
-----------------------------------------
- Header row: any line starting with `|` and containing at least two `|`.
- Separator row: the immediately following line, also starting with `|`,
  containing only `|`, `:`, `-`, and whitespace.
- Data rows: every subsequent line starting with `|` until we hit a blank
  line or a line not starting with `|`.

Header cell text is normalised lightly for CSV (collapse internal
whitespace, strip leading/trailing whitespace) but otherwise preserved
verbatim — Unicode greek letters, subscripts, units in parentheses, etc.
all stay.

Output
------
Each CSV is **ASCII**, newline-delimited, comma-separated, with a header
row and a trailing newline.  Greek letters, subscripts, em-dashes and
other Unicode in the markdown headers/cells are transliterated to
plain-ASCII equivalents (e.g. ``mu``, ``alpha``, ``r/R=0.75``, ``-``).
The full table in `_ASCII_MAP` documents the substitutions.  Empty
markdown cells (typically the OCR "|     |" alignment columns) become
empty CSV cells.

Skipping
--------
The following filenames are skipped because they are prose, not data:
- `CLAUDE.md` (any level)
- `README.md`
- `EMPIRICAL_VALIDATION.md`
- `*.summary.md`, `*.abstract.md`
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_ROOT = ROOT / "csv"
SKIP_NAMES = {"CLAUDE.md", "README.md", "EMPIRICAL_VALIDATION.md",
              "summary.md", "abstract.md", "paper.md", "extraction.md"}

# Unicode -> ASCII substitutions for CSV output.  The MD source keeps the
# original symbols; only the CSV writer uses this map.  Edit here if a
# paper uses a character not yet covered.
_ASCII_MAP = {
    # Greek letters
    "α": "alpha", "β": "beta",  "γ": "gamma", "δ": "delta",
    "ε": "eps",   "ζ": "zeta",  "η": "eta",   "θ": "theta",
    "λ": "lambda","μ": "mu",    "ν": "nu",    "ξ": "xi",
    "π": "pi",    "ρ": "rho",   "σ": "sigma", "τ": "tau",
    "φ": "phi",   "χ": "chi",   "ψ": "psi",   "ω": "omega",
    "Δ": "Delta", "Σ": "Sigma", "Ω": "Omega", "Θ": "Theta",
    "Π": "Pi",    "Φ": "Phi",   "Ψ": "Psi",
    # Subscript digits
    "₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4",
    "₅": "5", "₆": "6", "₇": "7", "₈": "8", "₉": "9",
    # Subscript symbols
    "₊": "+", "₋": "-", "₌": "=", "₍": "(", "₎": ")",
    # Superscript digits / symbols
    "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4",
    "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9",
    "⁺": "+", "⁻": "-",
    # Dashes / quotes / spaces
    "‐": "-", "‑": "-", "‒": "-", "–": "-", "—": "-",
    "―": "-", "−": "-",
    "‘": "'", "’": "'", "“": '"', "”": '"',
    " ": " ", " ": " ", " ": " ", " ": " ", "​": "",
    # Math / units
    "°": "deg", "×": "x", "÷": "/", "±": "+-",
    "≈": "~", "≠": "!=", "≤": "<=", "≥": ">=",
    "√": "sqrt", "∞": "inf", "∂": "d", "∇": "del",
    "∫": "int", "∑": "sum", "∏": "prod",
    # Fractions
    "¼": "1/4", "½": "1/2", "¾": "3/4",
    "⅓": "1/3", "⅔": "2/3",
    # Misc
    "·": ".", "•": "*", "…": "...",
    "»": ">>", "«": "<<", "→": "->", "←": "<-",
    "↑": "^",  "↓": "v",
    # Check / cross marks (often used in cross-validation tables)
    "✓": "OK", "✔": "OK", "✗": "X", "✘": "X",
    "⚠": "!", "❌": "X", "✅": "OK",
    # Zero-width / variation selectors (emoji presentation, BOMs, joiners)
    "️": "", "︎": "", "​": "", "‌": "", "‍": "",
    "﻿": "",
}


def to_ascii(text: str) -> str:
    """Apply _ASCII_MAP, then strip anything still non-ASCII (with a warning)."""
    for k, v in _ASCII_MAP.items():
        if k in text:
            text = text.replace(k, v)
    if not text.isascii():
        bad = sorted({c for c in text if not c.isascii()})
        print(f"warning: dropping non-ASCII chars {bad!r} from {text!r}",
              file=sys.stderr)
        text = "".join(c if c.isascii() else "?" for c in text)
    return text


def slugify(text: str, max_len: int = 60) -> str:
    """Lowercase, replace anything non-alphanum with single '_'."""
    text = text.strip().lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    text = text.strip("_")
    return text[:max_len] or "table"


def is_separator(line: str) -> bool:
    """True if line looks like '|---|:--:|---|' etc."""
    stripped = line.strip()
    if not stripped.startswith("|"):
        return False
    body = stripped.strip("|")
    return bool(body) and all(c in "-:| \t" for c in body)


def split_row(line: str) -> list[str]:
    """Split '| a | b | c |' into ['a', 'b', 'c'] (drop leading/trailing empties)."""
    line = line.strip()
    if line.startswith("|"):
        line = line[1:]
    if line.endswith("|"):
        line = line[:-1]
    return [cell.strip() for cell in line.split("|")]


def find_tables(md_text: str) -> list[tuple[str, list[list[str]]]]:
    """Return list of (heading_slug, table_rows) for each markdown table.

    `heading_slug` is derived from the nearest preceding `## ...` or
    `### ...` heading; empty string if no heading was found.
    `table_rows[0]` is the header; the rest are data.
    """
    lines = md_text.splitlines()
    tables: list[tuple[str, list[list[str]]]] = []
    current_heading = ""
    i = 0
    while i < len(lines):
        line = lines[i]
        m = re.match(r"^(#{2,4})\s+(.*?)\s*$", line)
        if m:
            current_heading = slugify(m.group(2))
            i += 1
            continue
        # Detect table start: header row followed by separator row.
        if (line.strip().startswith("|") and i + 1 < len(lines)
                and is_separator(lines[i + 1])):
            header = split_row(line)
            i += 2
            rows: list[list[str]] = [header]
            while i < len(lines) and lines[i].strip().startswith("|"):
                if is_separator(lines[i]):
                    i += 1
                    continue
                row = split_row(lines[i])
                # pad / truncate to header length
                if len(row) < len(header):
                    row += [""] * (len(header) - len(row))
                elif len(row) > len(header):
                    row = row[:len(header)]
                rows.append(row)
                i += 1
            # drop tables with no data
            if len(rows) > 1:
                tables.append((current_heading, rows))
            continue
        i += 1
    return tables


def csv_target(md_path: Path, suffix: str) -> Path:
    """Map Research/<Paper>/foo.md  →  Research/csv/<Paper>/foo<suffix>.csv."""
    rel = md_path.relative_to(ROOT)            # e.g. CaradonnaTung/page_10_table_1.md
    target = CSV_ROOT / rel.with_suffix("")    # e.g. Research/csv/CaradonnaTung/page_10_table_1
    target = target.with_name(target.name + suffix + ".csv")
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def emit_csvs(md_path: Path, tables: list[tuple[str, list[list[str]]]]) -> list[Path]:
    """Write one or more CSV files under Research/csv/.  Returns list of paths."""
    written: list[Path] = []
    if not tables:
        return written
    if len(tables) == 1:
        out = csv_target(md_path, "")
        with out.open("w", encoding="ascii", newline="") as f:
            csv.writer(f).writerows([[to_ascii(c) for c in r] for r in tables[0][1]])
        written.append(out)
        return written
    # Multiple tables: disambiguate by heading slug or counter.
    used_slugs: set[str] = set()
    for idx, (slug, rows) in enumerate(tables, start=1):
        candidate = slug or f"table_{idx:02d}"
        n = candidate
        k = 2
        while n in used_slugs:
            n = f"{candidate}_{k}"
            k += 1
        used_slugs.add(n)
        out = csv_target(md_path, f"__{n}")
        with out.open("w", encoding="ascii", newline="") as f:
            csv.writer(f).writerows([[to_ascii(c) for c in r] for r in rows])
        written.append(out)
    return written


def should_skip(path: Path) -> bool:
    if path.name in SKIP_NAMES:
        return True
    if path.name.endswith((".summary.md", ".abstract.md")):
        return True
    return False


def main(argv: list[str]) -> int:
    targets: list[Path]
    if argv:
        targets = [ROOT / a for a in argv]
    else:
        targets = [ROOT]
    md_files: list[Path] = []
    for t in targets:
        if not t.exists():
            print(f"skip: {t} does not exist", file=sys.stderr)
            continue
        md_files.extend(sorted(t.rglob("*.md")))

    total_tables = 0
    total_csvs = 0
    n_md_with_tables = 0
    for md in md_files:
        if should_skip(md):
            continue
        text = md.read_text(encoding="utf-8")
        tables = find_tables(text)
        if not tables:
            continue
        n_md_with_tables += 1
        written = emit_csvs(md, tables)
        total_tables += len(tables)
        total_csvs += len(written)
        rel = md.relative_to(ROOT)
        if len(written) == 1:
            print(f"{rel}  ->  {written[0].relative_to(ROOT)}")
        else:
            print(f"{rel}  ->  {len(written)} CSVs in {written[0].parent.relative_to(ROOT)}/")
            for w in written:
                print(f"     {w.name}")
    print()
    print(f"processed {n_md_with_tables} md files, "
          f"{total_tables} tables, wrote {total_csvs} CSVs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
