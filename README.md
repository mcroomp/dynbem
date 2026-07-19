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
`AeroModel` trait (Rust). The repo also includes a flight-envelope sweep
driver (`envelope/compute_map.py`), a cyclic-trim solver, and a
point-mass attitude simulator with cyclic-pitch control.

Coordinates are NED throughout; rotor rotation is CCW-from-above
(American helicopter convention). For empirical validation against
published rotor data see [EMPIRICAL_VALIDATION.md](docs/EMPIRICAL_VALIDATION.md).

## Energy-harvesting and free-rotor applications

The full-envelope coverage -- autorotation, windmill-brake state, oblique
descent, and vortex-ring state -- combined with VPM Level 3's explicit
wake geometry makes `dynbem` well suited for kites, autogiros, and other
free-rotors extracting power from wind. No regime switching or empirical
patches are needed across the power-extraction to hover-flight transition.
See [VPM_DESIGN.md](docs/VPM_DESIGN.md) Section 11 for the VPM roadmap.

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

# YAML loading is pure Python (yaml.safe_load); the resulting RotorDefinition
# holds lean Rust wrappers for the math fields.
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
    omega_rad_s=omega,                         # rotor speed passed in each call
    rho_kg_m3=1.225,
)
# All models advance with step() -- the same API across all three tiers:
result, state = model.step(inputs, state, dt, integration_method="semi_implicit")
# integration_method: "semi_implicit" | "explicit" | "exponential"  (ignored for quasi static and VPM)
# result.F_world, result.m_hub_world, result.M_spin, result.Q_spin

# Mechanical ODE lives in the caller. Coulomb bearing friction is always
# a parameter (default 0.0). Prefer semi_implicit_step_omega for
# fixed-timestep production use (unconditionally stable in the aero term).
from dynbem.mechanical import omega_derivative
omega += dt * omega_derivative(omega, result.Q_spin, motor_torque_Nm, I_ode_kgm2)

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

## Tests and validation

```
uv run pytest tests/ -q
```

**Validation report** (release build required -- VPM is 50-100x slower in debug):
```
cargo run --release -p validation_rs
# prints one CHECK line per data point; exits 0 if all pass
```

If `uv` is not on `PATH`, run pytest with the workspace interpreter directly:

```
.venv/Scripts/python.exe -m pytest tests/ -q    # Windows
.venv/bin/python -m pytest tests/ -q            # POSIX
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

## Performance

Approximate throughput (36 azimuth stations, 15 radial elements):

| Model | Hover | Forward / cyclic |
|---|---:|---:|
| `bem` | about 50k iter/sec | about 1.3k iter/sec |
| `pitt_peters` | about 3M iter/sec | about 200k iter/sec |
| `oye` | about 200k iter/sec | about 185k-200k iter/sec |
| `vpm` | ~100 steps/sec | ~100 steps/sec |

VPM is ~10 ms/step (Barnes-Hut accelerated, ~800-3000 particles). A typical
steady-state run is 6-10 revolutions at 48 steps/rev -- roughly 3-5 seconds.

These are rough scale numbers. For before/after comparisons on
performance-sensitive changes, use the Criterion suite at
[`dynbem_rs/benches/model_kernels.rs`](dynbem_rs/benches/model_kernels.rs).
An external profiling binary is also available at
[`dynbem_rs/bin/profile_kernels.rs`](dynbem_rs/bin/profile_kernels.rs)
(built with `cargo build --release -p dynbem_rs`; debuginfo is kept
enabled in the release profile for symbol resolution).

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

## AI tooling notes

- `CLAUDE.md` is the Claude Code instruction file for Claude-based agents.
- GitHub tooling (Copilot, GitHub CLI agent) uses `AGENTS.md` as its
  default instruction file instead.
- Both files carry the same conventions; keep both if you use both tools.

## License

MIT — see `LICENSE`.
