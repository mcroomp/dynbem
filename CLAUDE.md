# CLAUDE.md — AI assistant instructions

Human-facing docs (install, usage, coordinate conventions, design notes,
implementation roadmap, research sources) live in [README.md](README.md).
**Read it first** — most of what you'd want to know about this codebase is
there, and this file does not repeat it.

This file holds only directives that are specifically for you (the AI
assistant) and would be noise in the README.

## Workflow

- **Python**: use the venv at `.venv\`. Activate with
  `.venv\Scripts\activate` (Windows), or invoke directly via
  `.venv\Scripts\python` / `.venv\Scripts\pytest`. Don't install packages
  globally or create a new venv.
- **Coordinate frame**: NED everywhere. See README "Coordinate system" —
  the "coordinate trap" section especially matters when you adapt
  equations from a paper, because most rotor literature uses a different
  frame and the sign flips are easy to miss.
- **Sign conventions**: before changing any inflow / thrust / torque
  sign, re-read the README "BEM solver design" and "Pitt-Peters design
  notes" sections. The signs are load-bearing and were tuned to make
  hover, climb, descent, VRS, and autorotation all work in one code
  path.

## When extending the aero models

- New levels (e.g. Peters-He) plug in behind the `AeroBase` interface
  in `aero/__init__.py`. Don't break existing call sites — keep
  `compute_forces(inputs, state) -> (AeroResult, RotorState)`.
- Validation data lives under `Research/`. When adding a new model
  level, add a `tests/test_<model>.py` and, if appropriate, a
  `val_step*.py` script that compares against a specific paper's data.
- Don't store derived results inside `Research/` — that directory is
  for source-paper extractions only.

## Subfolder CLAUDE.md files

- `Research/CLAUDE.md` — extraction conventions for paper sources.
- `Research/CaradonnaTung/CLAUDE.md` — Caradonna-Tung page index, CT
  tables, validation notes.

Defer to those when working inside the respective directories.
