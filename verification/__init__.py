"""Whole-dataset validation scripts. Each module exposes a
`run_*` / `run_survey` function that unit tests in `tests/` call with
a small `sample=` argument; the module's own `main()` runs the full
sweep for re-baselining bounds. See CLAUDE.md for the policy."""
