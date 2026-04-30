#!/usr/bin/env python3
"""Validate extracted Cp table markdown files for physical plausibility.

Checks
------
  x/c_range       x/c not in [0, 1]                         FAIL
  x/c_monotonic   x/c decreases within a radial station      WARN
  x/c_duplicate   same x/c appears twice in a station        WARN
  cp_magnitude    |Cp| > CP_WARN (5.0) or > CP_FAIL (15.0)  WARN / FAIL
  cl_range        CL < -0.05 or > 1.0                        FAIL
  cl_expected     CL outside rough θc-based range            WARN
  stray_sep       separator row found inside a data table    WARN
  bad_chars       non-ASCII / lookalike characters in a cell FAIL
  unreadable      cell value is '?'                          INFO

Usage
-----
    python validate_tables.py            # all page_*_table_*.md here
    python validate_tables.py file.md    # one or more specific files
    python validate_tables.py -v ...     # include INFO lines
"""

import re
import sys
import unicodedata
from pathlib import Path
from dataclasses import dataclass
from typing import Optional

CP_WARN = 5.0
CP_FAIL = 15.0

# Characters that look like ASCII punctuation/digits but are not.
# Maps code point → human-readable name for the error message.
LOOKALIKES: dict[int, str] = {
    0x2212: 'MINUS SIGN (U+2212)',        # −  looks like -
    0x2013: 'EN DASH (U+2013)',           # –
    0x2014: 'EM DASH (U+2014)',           # —
    0x2010: 'HYPHEN (U+2010)',            # ‐
    0x2011: 'NON-BREAKING HYPHEN',        # ‑
    0x00AD: 'SOFT HYPHEN (U+00AD)',       # ­
    0xFE63: 'SMALL HYPHEN-MINUS',         # ﹣
    0xFF0D: 'FULLWIDTH HYPHEN-MINUS',     # －
    0x00B7: 'MIDDLE DOT (U+00B7)',        # ·  looks like decimal point
    0xFF0E: 'FULLWIDTH FULL STOP',        # ．
    **{0xFF10 + i: f'FULLWIDTH DIGIT {i}' for i in range(10)},   # ０–９
    **{0x2070 + i: f'SUPERSCRIPT {i}' for i in range(10)},       # ⁰–⁹ (sparse)
    0x00B9: 'SUPERSCRIPT ONE',
    0x00B2: 'SUPERSCRIPT TWO',
    0x00B3: 'SUPERSCRIPT THREE',
    0x207A: 'SUPERSCRIPT PLUS',
    0x207B: 'SUPERSCRIPT MINUS',
    0x0045: None,  # plain E — fine, just here as a reminder NOT to flag it
}
# Remove the sentinel
del LOOKALIKES[0x0045]


def check_cell_chars(value: str, context: str) -> list['Issue']:
    """Return FAIL issues for any non-ASCII or lookalike character in a numeric cell."""
    issues = []
    if not value or value in ('?', ''):
        return issues
    seen: set[int] = set()
    for ch in value:
        cp = ord(ch)
        if cp <= 127:
            continue
        if cp in seen:
            continue
        seen.add(cp)
        if cp in LOOKALIKES:
            name = LOOKALIKES[cp]
            issues.append(Issue('FAIL', context, None,
                                f'lookalike character {name} (U+{cp:04X}) '
                                f'in cell {value!r} — replace with ASCII'))
        else:
            cat = unicodedata.category(ch)
            name = unicodedata.name(ch, f'U+{cp:04X}')
            issues.append(Issue('FAIL', context, None,
                                f'non-ASCII character {name!r} (category {cat}) '
                                f'in cell {value!r}'))
    return issues

STATIONS = ['r/R=0.50', 'r/R=0.68', 'r/R=0.80', 'r/R=0.89', 'r/R=0.96']

# Expected CL range per θc (lo, hi) — generous bounds, just catches gross errors
CL_EXPECTED = {
    0:  (0.00, 0.05),
    2:  (0.01, 0.15),
    4:  (0.02, 0.25),
    5:  (0.05, 0.30),
    7:  (0.10, 0.50),
    8:  (0.15, 0.55),
    10: (0.20, 0.70),
    12: (0.25, 0.80),
}


@dataclass
class Issue:
    level: str          # FAIL / WARN / INFO
    section: str
    station: Optional[str]
    detail: str

    def __str__(self):
        loc = self.section
        if self.station:
            loc += f' / {self.station}'
        return f'  {self.level:<4}  {loc}: {self.detail}'


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def parse_float(s: str) -> Optional[float]:
    s = s.strip()
    if not s or s == '?':
        return None
    try:
        return float(s)
    except ValueError:
        return None


def parse_row(line: str) -> list[str]:
    s = line.strip()
    if s.startswith('|'):
        s = s[1:]
    if s.endswith('|'):
        s = s[:-1]
    return [c.strip() for c in s.split('|')]


def is_separator(cells: list[str]) -> bool:
    return bool(cells) and all(re.fullmatch(r'-+', c) or c == '' for c in cells)


def parse_header(text: str) -> tuple[Optional[int], Optional[float]]:
    """Return (theta_deg, m_tip) from the title line."""
    m = re.search(r'θc=(\d+)°', text)
    theta = int(m.group(1)) if m else None
    m = re.search(r'M_tip=([\d.]+)', text)
    mtip = float(m.group(1)) if m else None
    return theta, mtip


def extract_section(lines: list[str], keyword: str, stop_keywords: list[str]) -> list[str]:
    """Return lines belonging to the named section."""
    start = None
    for i, l in enumerate(lines):
        if keyword.lower() in l.lower():
            start = i
            break
    if start is None:
        return []
    result = []
    for l in lines[start + 1:]:
        if any(k.lower() in l.lower() for k in stop_keywords):
            break
        result.append(l)
    return result


def parse_table(lines: list[str]) -> tuple[list[str], list[list[str]], list[int]]:
    """
    Parse a markdown table from a block of lines.
    Returns (headers, data_rows, stray_sep_indices).
    stray_sep_indices are 0-based indices into `lines` of separator rows
    that appear after the first (expected) header separator.
    """
    table_lines = [(i, l) for i, l in enumerate(lines) if l.strip().startswith('|')]
    if not table_lines:
        return [], [], []

    headers = None
    header_sep_done = False
    data_rows = []
    stray_seps = []

    for raw_i, (i, line) in enumerate(table_lines):
        cells = parse_row(line)
        if headers is None:
            headers = cells
        elif not header_sep_done:
            if is_separator(cells):
                header_sep_done = True
            # else: treat as data (no separator row — unusual)
        else:
            if is_separator(cells):
                stray_seps.append(i)
            else:
                data_rows.append(cells)

    return headers or [], data_rows, stray_seps


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_station_column(rows: list[list[str]], xc_col: int, cp_col: int,
                             station: str, section: str) -> list[Issue]:
    issues = []
    prev_xc: Optional[float] = None
    seen_xc: set[float] = set()

    for cells in rows:
        xc_raw = cells[xc_col] if xc_col < len(cells) else ''
        cp_raw = cells[cp_col] if cp_col < len(cells) else ''

        if not xc_raw and not cp_raw:
            continue  # blank padding row for this station — fine

        # --- character checks ---
        ctx = f'{section} / {station}'
        issues.extend(check_cell_chars(xc_raw, ctx))
        issues.extend(check_cell_chars(cp_raw, ctx))

        # --- x/c checks ---
        xc = parse_float(xc_raw)
        if xc_raw and xc_raw != '?':
            if xc is None:
                issues.append(Issue('WARN', section, station,
                                    f'x/c unparseable: {xc_raw!r}'))
            else:
                if xc < 0.0 or xc > 1.0:
                    issues.append(Issue('FAIL', section, station,
                                        f'x/c={xc} outside [0, 1]'))
                elif xc in seen_xc:
                    issues.append(Issue('WARN', section, station,
                                        f'duplicate x/c={xc}'))
                elif prev_xc is not None and xc < prev_xc:
                    issues.append(Issue('WARN', section, station,
                                        f'x/c={xc} non-monotonic (previous {prev_xc:.4f})'))
                seen_xc.add(xc)
                if xc not in seen_xc or (prev_xc is None or xc > prev_xc):
                    prev_xc = xc

        # --- Cp checks ---
        if cp_raw == '?':
            issues.append(Issue('INFO', section, station,
                                f'unreadable value at x/c={xc_raw}'))
        elif cp_raw:
            cp = parse_float(cp_raw)
            if cp is None:
                issues.append(Issue('WARN', section, station,
                                    f'Cp unparseable: {cp_raw!r} at x/c={xc_raw}'))
            elif abs(cp) > CP_FAIL:
                issues.append(Issue('FAIL', section, station,
                                    f'|Cp|={abs(cp):.3g} at x/c={xc_raw} — '
                                    f'exceeds {CP_FAIL}, almost certainly scan error'))
            elif abs(cp) > CP_WARN:
                issues.append(Issue('WARN', section, station,
                                    f'|Cp|={abs(cp):.3g} at x/c={xc_raw} — '
                                    f'suspicious (>{CP_WARN}); check scan'))

    return issues


def validate_file(path: Path, verbose: bool = False) -> list[Issue]:
    text = path.read_text(encoding='utf-8')
    lines = text.splitlines()
    issues: list[Issue] = []

    theta, mtip = parse_header(text)

    # --- Upper and lower surface ---
    for section_key, section_label in (('upper surface', 'Upper surface'),
                                        ('lower surface', 'Lower surface')):
        stop = ['lower surface', '**cl', '## cl'] if 'upper' in section_key \
               else ['**cl', '## cl']
        sec_lines = extract_section(lines, section_key, stop)
        _, data_rows, stray_seps = parse_table(sec_lines)

        for si in stray_seps:
            issues.append(Issue('WARN', section_label, None,
                                f'stray separator row at section line {si}'))

        for s_idx, station in enumerate(STATIONS):
            xc_col = s_idx * 2
            cp_col = s_idx * 2 + 1
            issues.extend(
                validate_station_column(data_rows, xc_col, cp_col, station, section_label)
            )

    # --- CL row ---
    cl_lines = extract_section(lines, '**cl', [])
    if not cl_lines:
        cl_lines = extract_section(lines, '## cl', [])
    _, cl_data, _ = parse_table(cl_lines)
    if cl_data:
        cells = cl_data[0]
        lo, hi = CL_EXPECTED.get(theta, (0.0, 1.0)) if theta is not None else (0.0, 1.0)
        for s_idx, station in enumerate(STATIONS):
            raw = cells[s_idx] if s_idx < len(cells) else ''
            issues.extend(check_cell_chars(raw, f'CL / {station}'))
            if raw == '?':
                issues.append(Issue('INFO', 'CL', station, 'unreadable (?)'))
                continue
            val = parse_float(raw)
            if val is None:
                continue
            if val < -0.05:
                issues.append(Issue('FAIL', 'CL', station,
                                    f'CL={val:.4f} — negative for positive θc'))
            elif val > 1.0:
                issues.append(Issue('FAIL', 'CL', station,
                                    f'CL={val:.4f} — exceeds 1.0'))
            elif not (lo <= val <= hi):
                issues.append(Issue('WARN', 'CL', station,
                                    f'CL={val:.4f} outside expected [{lo:.3f}, {hi:.3f}] '
                                    f'for θc={theta}°'))

    return issues


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

    verbose = '-v' in sys.argv
    args = [a for a in sys.argv[1:] if not a.startswith('-')]

    if args:
        targets = [Path(a) for a in args]
    else:
        here = Path(__file__).parent
        targets = sorted(here.glob('page_*_table_*.md'))

    if not targets:
        print('No files found.')
        return

    total_fails = 0
    total_warns = 0
    total_info = 0

    for path in targets:
        issues = validate_file(path, verbose)
        fails  = [i for i in issues if i.level == 'FAIL']
        warns  = [i for i in issues if i.level == 'WARN']
        infos  = [i for i in issues if i.level == 'INFO']
        total_fails += len(fails)
        total_warns += len(warns)
        total_info  += len(infos)

        if not fails and not warns and (not infos or not verbose):
            n_q = len(infos)
            suffix = f'  ({n_q} unreadable)' if n_q else ''
            print(f'  OK    {path.name}{suffix}')
            continue

        status = 'FAIL' if fails else 'WARN'
        print(f'\n[{status}] {path.name}')
        for iss in issues:
            if iss.level == 'INFO' and not verbose:
                continue
            print(iss)
        if infos and not verbose:
            print(f'  INFO  {len(infos)} unreadable value(s) — run with -v to list')

    print(f'\n{"=" * 60}')
    print(f'Files checked : {len(targets)}')
    print(f'FAIL          : {total_fails}')
    print(f'WARN          : {total_warns}')
    print(f'INFO (unreads): {total_info}')
    if total_fails == 0 and total_warns == 0:
        print('All clean.')


if __name__ == '__main__':
    main()
