# De Schutter Validation Stages (1-3)

This folder contains a staged validation scaffold for applying the
De Schutter 2018 model logic to dynbem aero outputs.

## Files

- rotor.yaml
  Reference rotor fixture for staged checks.

- step1_reference_point.py
  Stage 1: one canonical operating-point run to validate setup and
  sign/magnitude sanity.

- step2_equation_parity.py
  Stage 2: equation-level parity checks for key closed-form terms
  (Eq.25 lift/drag and Eq.29 structural drag coefficient).

- step3_steady_sweep.py
  Stage 3: deterministic steady operating-point benchmark sweep with
  CSV output at deschutter/out/steady_sweep.csv.

- compare_deschutter_vs_dynbem.py
  Point-by-point A/B comparison between a De Schutter-style reduced
  closed-form reference and dynbem across the same grid. Writes
  deschutter/out/deschutter_vs_dynbem.csv and prints aggregate errors.

- isolate_divergence.py
  Pandas-based analysis of deschutter_vs_dynbem.csv that isolates where
  errors cluster by wind, rpm, and collective, and writes
  deschutter/out/divergence_report.txt.

## Run

From repository root:

    uv run python deschutter/step1_reference_point.py
    uv run python deschutter/step2_equation_parity.py
    uv run python deschutter/step3_steady_sweep.py
    uv run python deschutter/compare_deschutter_vs_dynbem.py
    uv run python deschutter/isolate_divergence.py

## Notes

- The fixture is intentionally simple and linearized for early-stage
  validation before adding pumping-cycle dynamics.
- Stage 3 includes a monotonic CT-vs-collective trend sanity check.
