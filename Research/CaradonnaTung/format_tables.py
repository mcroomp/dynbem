#!/usr/bin/env python3
"""Reformat markdown tables so columns align in plain-text view.

Usage:
    python format_tables.py            # formats all page_*_table_*.md in this directory
    python format_tables.py file.md    # formats a single file
"""

import re
import sys
from pathlib import Path


def parse_row(line):
    s = line.strip()
    if s.startswith('|'):
        s = s[1:]
    if s.endswith('|'):
        s = s[:-1]
    return [c.strip() for c in s.split('|')]


def is_separator(cells):
    return all(re.fullmatch(r'-+', c) or c == '' for c in cells)


def format_table(lines):
    rows = [parse_row(l) for l in lines]
    n_cols = max(len(r) for r in rows)
    for r in rows:
        while len(r) < n_cols:
            r.append('')

    widths = [3] * n_cols
    for r in rows:
        if not is_separator(r):
            for i, cell in enumerate(r):
                widths[i] = max(widths[i], len(cell))

    out = []
    for r in rows:
        if is_separator(r):
            cells = ['-' * widths[i] for i in range(n_cols)]
        else:
            cells = [r[i].ljust(widths[i]) for i in range(n_cols)]
        out.append('| ' + ' | '.join(cells) + ' |')
    return out


def process_file(path):
    lines = path.read_text(encoding='utf-8').splitlines()
    out = []
    i = 0
    while i < len(lines):
        if lines[i].strip().startswith('|'):
            block = []
            while i < len(lines) and lines[i].strip().startswith('|'):
                block.append(lines[i])
                i += 1
            out.extend(format_table(block))
        else:
            out.append(lines[i])
            i += 1
    path.write_text('\n'.join(out) + '\n', encoding='utf-8')


def main():
    if len(sys.argv) > 1:
        targets = [Path(a) for a in sys.argv[1:]]
    else:
        here = Path(__file__).parent
        targets = sorted(here.glob('page_*_table_*.md'))

    if not targets:
        print('No files found.')
        return

    for f in targets:
        print(f'  {f.name}')
        process_file(f)
    print(f'Done — {len(targets)} file(s) formatted.')


if __name__ == '__main__':
    main()
