# dynbem

**Rotor aerodynamics library -- from fast BEM inflow models to a
full vortex-particle free-wake solver, in one code path.**

`dynbem` is a rotor-aerodynamics library covering three model tiers, all
sharing the same API, coordinate conventions, and blade geometry format:

- **Level 1** -- quasi-static BEM (`"bem"`): single-call blade-element
  momentum, numerically valid from hover through windmill-brake state,
  autorotation, and forward flight via the Ning windmill Brent solver.
- **Level 2** -- dynamic-inflow BEM (`"pitt_peters"`, `"oye"`): same BEM
  kernel, augmented with a dynamic wake-state ODE that captures inflow
  transients. Pitt-Peters uses the Peters L-matrix with Glauert wake-skew
  coupling; Oye uses per-annulus filtered momentum inflow (OpenFAST DBEMT
  formulation).
- **Level 3** -- free-wake VPM (`"vpm"`): a vortex-particle method that
  represents the rotor wake explicitly as a time-evolving cloud of
  regularized vorticity. No wake-geometry assumptions; wake skew, VRS
  onset, and autorotation torque sign emerge from the physics rather than
  empirical patches. Cost ~10 ms/step (Barnes-Hut accelerated).

All three tiers are numerically valid across the **full operating envelope** --
helicopter hover, axial climb, axial descent, vortex-ring state (VRS),
windmill-brake state (WBS), autorotation, and wind-turbine power extraction
-- without switching equations or sign conventions between regimes.

The math core is a pure-Rust crate ([`dynbem_rs/`](dynbem_rs/), no pyo3 /
numpy / file IO) wrapped by a thin PyO3 + maturin binding crate
([`dynbem/`](dynbem/)) which is the publishable Python package.

All models share a tabulated polar interpolator and plug into the same
`AeroModel` trait (Rust).
The repo also includes a flight-envelope sweep driver
(`envelope/compute_map.py`), a cyclic-trim solver
([`dynbem_rs/src/trim.rs`](dynbem_rs/src/trim.rs)), and a point-mass +
cyclic-pitch attitude simulator.
For empirical validation against published rotor data (Castles-Gray
TN-2474 vertical descent, Caradonna-Tung TM-81232 hover CT and
spanwise CL, Harrington TN-2318 full-scale hover, Wheatley & Hood
TR 515 forward-flight autorotation), see
[EMPIRICAL_VALIDATION.md](docs/EMPIRICAL_VALIDATION.md).

Coordinates are NED throughout; rotor rotation is CCW-from-above
(American helicopter convention).

### RAWES applications

`dynbem` is particularly well-suited for **rotor-as-wind-energy-systems
(RAWES)** -- kites, autogiros, and other free-rotors extracting power from
wind. The library's coverage of autorotation, windmill-brake state,
oblique descent, and wake geometry (especially VPM Level 3) makes it
ideal for modeling the full flight envelope of energy-harvesting rotors.
See [VPM_DESIGN.md](docs/VPM_DESIGN.md) Section 11 for the roadmap on RAWES
fidelity improvements.

## Install

**Quick start** (all platforms):

```
./setup.sh           # POSIX shell or Windows git-bash/WSL
setup.cmd            # Windows cmd/PowerShell
```

Both scripts check prerequisites (uv, cargo), then run `uv sync --group dev`,
which creates the `.venv`, installs dependencies, and builds the Rust
extension via maturin. See [setup.sh](setup.sh) for details.

**Manual setup** (equivalent to what the scripts run):

The repo is a uv workspace + Cargo workspace. The Rust extension is built
automatically (via maturin) by `uv sync`:

```
uv sync                 # builds dynbem (Rust extension) editable
uv sync --group dev     # also installs pytest, maturin, build, twine
uv run pytest tests/    # run the test suite
```

The publishable Python package is [`dynbem/`](dynbem/) (Rust-backed via
PyO3). Inside a non-uv environment you can still `pip install -e dynbem/`
— maturin will pick up `dynbem/pyproject.toml` and compile the extension
against the sibling [`dynbem_rs/`](dynbem_rs/) crate. Requires a working
Rust toolchain (`rustup` stable).

## Usage

```python
import numpy as np
import dynbem

# Load rotor definition from YAML. Parsing happens in Rust (dynbem_rs)
# via PyO3 bindings; pure-Rust callers can use RotorDefinition::from_yaml_file(path).
defn   = dynbem.rotor_definition.load("rotors/castles_gray_6ft/rotor.yaml")
model  = dynbem.create_aero(defn, model="pitt_peters")  # or "oye", "bem", "vpm"
state  = model.initial_rotor_state()

omega = 125.7   # rad/s -- caller owns mechanical state
inputs = dynbem.RotorInputs(
    collective_rad=0.14,
    tilt_lon=0.0, tilt_lat=0.0,                # swashplate (helicopter-standard signs)
    R_hub=np.eye(3),
    v_hub_world=np.zeros(3),
    wind_world=np.zeros(3),
    t=0.0,
    omega_rad_s=omega,                         # rotor speed passed in each call
)
# All models advance with step() -- the same API across all three tiers:
result, state = model.step(inputs, state, dt, integration_method="semi_implicit")
# integration_method: "semi_implicit" | "explicit" | "exponential"  (ignored for "bem")
# result.F_world, result.m_hub_world, result.M_spin, result.Q_spin

# Mechanical ODE lives in the caller:
from dynbem.mechanical import omega_derivative
omega += dt * omega_derivative(result.Q_spin, motor_torque_Nm, I_ode_kgm2)

# Level-3 VPM: same step() API; the free-wake particle cloud advances each call.
# dt typically T_rev/36-72; average over several revolutions for steady-state loads.
vpm   = dynbem.create_aero(defn, model="vpm")
state = vpm.initial_rotor_state()
result, state = vpm.step(inputs, state, dt)   # no integration_method needed
# result fields are the same: F_world, m_hub_world, Q_spin, M_spin

# compute_forces() is also available when you want to manage integration yourself:
# result, derivative = model.compute_forces(inputs, state)
```

For the full API reference — all classes, fields, keyword arguments, and
return types — see **[API.md](docs/API.md)**.

## Flight envelope sweep

```
run_map.cmd                                           # quick grid, saves to out\map.npz, plots to out\
run_map.cmd --full --save out\map.npz --plot out\     # full grid
uv run python -m envelope.compute_map --help
```

## Tests and validation

```
uv run pytest tests/ -q
```

**Validation report** (release build required -- VPM is 50-100x slower in debug):
```
cargo run --release -p validation_rs
# prints one CHECK line per data point; exits 0 if all pass
```

If `uv` is not on `PATH` in your shell, run pytest with the workspace
interpreter directly:

```
c:/repos/aero/.venv/Scripts/python.exe -m pytest tests/ -q
```

The `tests/` directory contains unit tests, validation scripts against
published rotor data, and end-to-end force-balance / frame-transform
checks. Whole-dataset validation sweeps against each paper live in
[`verification/`](verification/) and are imported by the matching
`tests/test_<paper>_<quantity>.py` in sampled mode -- one source of
BEM-driver logic, fast tests, and a full-sweep script you can re-run to
refresh aggregate bounds. For which papers and tables the models are
checked against, what the achieved variance is, and the physical
reasons for any residual bias, see
[EMPIRICAL_VALIDATION.md](docs/EMPIRICAL_VALIDATION.md).

## Performance baseline (Criterion)

The Rust crate includes a stable Criterion benchmark suite at
[`dynbem_rs/benches/model_kernels.rs`](dynbem_rs/benches/model_kernels.rs).
Use this as the baseline before and after any performance refactor.

### Bench groups

- `solve_bem_element`: single-element BEM hot kernel cost
- `sweep_scalar/prescribed/...`: scalar psi-loop sweep at fixed `(n_psi, n_elements)`
- `models_compute_forces/{bem,pitt_peters,oye}`: model-level `compute_forces` cost

### Recommended commands

```
cargo bench -p dynbem_rs --bench model_kernels --no-run
cargo bench -p dynbem_rs --bench model_kernels -- "sweep_scalar|models_compute_forces" --sample-size 20 --measurement-time 3
```

For change validation, run the same command before/after your patch and
compare medians in each benchmark group.

### Quick throughput snapshot

For a quick local sanity check, the Python-facing Rust hot-path harness
[`dynbem/benchmarks/bench_rust_only.py`](dynbem/benchmarks/bench_rust_only.py)
reports minimum call time across repeated trials. On a recent local run
with 36 azimuth stations and 15 radial elements, approximate throughput was:

| Model | Hover | Forward / cyclic |
|---|---:|---:|
| `bem` | about 50k iter/sec | about 1.3k iter/sec |
| `pitt_peters` | about 3M iter/sec | about 200k iter/sec |
| `oye` | about 200k iter/sec | about 185k-200k iter/sec |
| `vpm` | ~100 steps/sec | ~100 steps/sec |

VPM throughput is roughly 10 ms/step (Barnes-Hut accelerated, ~800-3000 particles).
A typical steady-state run is 6-10 revolutions at 48 steps/rev, so about 3-5 seconds.

Treat these as rough scale numbers, not a release guarantee; use Criterion
before/after comparisons for performance-sensitive changes.

## External profiling harness

A standalone profiling binary is included at
[`dynbem_rs/bin/profile_kernels.rs`](dynbem_rs/bin/profile_kernels.rs).
It builds with normal package builds (see `[[bin]]` in
[`dynbem_rs/Cargo.toml`](dynbem_rs/Cargo.toml)).

```
cargo build --release -p dynbem_rs
./target/release/profile_kernels.exe oye
./target/release/profile_kernels.exe pitt_peters
./target/release/profile_kernels.exe solve_bem_element
```

Release profiles in the workspace keep debuginfo enabled so external
profilers can resolve symbols.

## AI tooling notes

- `CLAUDE.md` is the Claude Code instruction file and remains useful if
  Claude-based agents are in your workflow.
- GitHub tooling (Copilot coding agent / GitHub CLI agent workflows)
  does not use `CLAUDE.md` as its default instruction file.
- For GitHub-side defaults, use `AGENTS.md` (repo-level coding-agent
  instructions), and optionally `.github/copilot-instructions.md` for
  Copilot-specific repository guidance.

Recommendation: keep `CLAUDE.md` if you use Claude tools, but add and
maintain an `AGENTS.md` so GitHub agent flows pick up the same policy.

---

## Design documentation

| Document | Content |
|---|---|
| [BEM_COMMON.md](docs/BEM_COMMON.md) | Coordinate system, kinematics, force kernel, QS BEM, VRS, servo-flap, output assembly, model comparison |
| [PITT_PETERS_DESIGN.md](docs/PITT_PETERS_DESIGN.md) | Pitt-Peters 3-state dynamic inflow (formal math + implementation notes) |
| [OYE_DESIGN.md](docs/OYE_DESIGN.md) | Oye 2-stage annular dynamic inflow |
| [VPM_DESIGN.md](docs/VPM_DESIGN.md) | VPM free-wake solver |

## Research sources

Extracted tables and figures from primary literature live under
`Research/`. Each paper subfolder uses the convention
`page_NN_<description>.md` so extractions trace back to their source
page image.

- **CaradonnaTung/** --- NASA TM-81232 (1981). 2-blade NACA 0012 hover CT data at theta_c = 5/8/12 deg.
- **Buhl_NREL_TP500_36834/** --- NREL TP-500-36834 (2005). Windmill Brake State correction extending Glauert.
- **Castles_TN2474/** --- NACA TN-2474 (Castles & Gray, 1951). Induced velocity in hover/descent.
- **Harrington_TN2318/** --- NACA TN-2318 (Harrington, 1951). Hover CT vs CP polars for two full-scale rotors.

## License

MIT — see `LICENSE`.
