# ccblade_docker

Containerised CCBlade (NREL's open-source pure-BEM solver, distributed
inside WISDEM) for generating reference BEM curves that
`verification/dynbem_vs_ccblade_*.py` compares against.

The point of the container is isolation: WISDEM pulls ~100 transitive
dependencies (OpenMDAO, MoorPy, netCDF, pandas, polars, Jupyter,
orbit-nrel, ...) and is genuinely painful to install into a
project venv on Windows. Containerising it means:

* the heavy install is one-shot and cached in the image layer,
* the host venv stays small,
* CI can build the image once and reuse the layer,
* the generated CSVs are checked into `verification/data/` so every
  unit test reads cached numbers without ever touching WISDEM.

## Layout

    Dockerfile                          -- python:3.11-slim + wisdem==4.2.1
    run.py                              -- entrypoint, reads /in/<name>.yaml
    compose.yml                         -- bind-mounted convenience runner
    inputs/
      beaupoil_2026.yaml                -- Beaupoil RAWES rotor (default config)
      _make_beaupoil_2026_yaml.py       -- regenerate from rotors/beaupoil_2026/
      nrel_phase_vi.yaml                -- NREL Phase VI Sequence S
      _make_nrel_phase_vi_yaml.py       -- regenerate from rotors/nrel_phase_vi/
    ../data/
      beaupoil_2026.csv                 -- CCBlade output (committed)
      nrel_phase_vi_seqS.csv            -- CCBlade output (committed once generated)

## Usage

Build once:

```bash
cd verification/ccblade_docker
docker compose build      # or: docker build -t aero-ccblade:4.2.1 .
```

Run the default Beaupoil sweep (4 blades, 25 operating points
covering V_wind 5..16 m/s x Omega 100..300 rpm at pitch 0 deg):

```bash
docker compose run --rm ccblade
# -> writes verification/data/beaupoil_2026.csv
```

Run a different config -- drop a `<name>.yaml` into `inputs/` and:

```bash
docker compose run --rm ccblade --config /in/<name>.yaml
```

Without compose, equivalent invocations (Git Bash on Windows needs
MSYS_NO_PATHCONV=1 so the `/in` path is not rewritten to a Windows
path):

```bash
docker build -t aero-ccblade:4.2.1 .
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "$(pwd)/inputs":/in:ro \
  -v "$(pwd)/../data":/out \
  aero-ccblade:4.2.1 --config /in/beaupoil_2026.yaml
```

## Output CSV schema

One row per operating point in the config:

| Column      | Units | Meaning |
|-------------|-------|---------|
| U_wind_ms   | m/s   | Free-stream wind speed (input) |
| Omega_rpm   | rpm   | Rotor speed (input) |
| pitch_deg   | deg   | Blade pitch (input, positive = pitch-to-feather) |
| T_N         | N     | Rotor thrust (axial force on shaft) |
| Q_Nm        | N*m   | Aerodynamic shaft torque |
| P_W         | W     | Aerodynamic power |
| CP          | --    | Power coefficient |
| CT          | --    | Thrust coefficient |
| CQ          | --    | Torque coefficient |

CCBlade normalises CT/CP/CQ by `rho * A * U_wind^N` per turbine
convention (not `rho * A * (Omega*R)^N`); be careful when comparing
to helicopter-convention coefficients.

## Config YAML schema

See the header of `run.py` for the full spec. Minimum example:

```yaml
name: my_rotor
rotor:
  n_blades: 2
  R_hub_m: 1.257
  R_tip_m: 5.029
air:
  rho: 1.225
blade_stations:
  - {r_m: 1.257, chord_m: 0.737, twist_deg: 22.00}
  - ...
airfoil_polar:
  Re: 1.0e6
  points:
    - {alpha_deg: -13.5, cl: -0.768, cd: 0.048}
    - ...
operating_points:
  - {U_wind_ms: 7.0, Omega_rpm: 72.0, pitch_deg: 3.0}
  - ...
```

The XFOIL header rows in the source polar CSV are skipped by the
helper script; only `Alpha,Cl,Cd` columns are used.

## Regenerating an input YAML

If the chord/twist table or polar source is updated:

```bash
.venv/Scripts/python verification/ccblade_docker/inputs/_make_beaupoil_2026_yaml.py
.venv/Scripts/python verification/ccblade_docker/inputs/_make_nrel_phase_vi_yaml.py
```

Both run on the host (no Docker), read the relevant CSVs / rotor.yaml
under `rotors/`, and emit a fresh YAML alongside themselves.

## When the CSV is regenerated

Re-run the container any time:

* the WISDEM/CCBlade version pin changes,
* the rotor fixture in `inputs/nrel_phase_vi.yaml` is updated,
* a new operating point is added to the sweep.

Commit the resulting `verification/data/nrel_phase_vi_seqS.csv` so
unit tests stay deterministic.
