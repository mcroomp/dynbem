# dynbem

**Dynamic blade-element momentum rotor aerodynamics — helicopter and
wind-turbine modes in one code path.**

`dynbem` is a rotor-aerodynamics library built around a multi-element
blade-element-momentum (BEM) solver coupled to dynamic-inflow models. It is
designed to be numerically valid across the **full operating envelope** —
helicopter hover, axial climb, axial descent, vortex-ring state (VRS),
windmill-brake state (WBS), autorotation, and wind-turbine power extraction
— without switching equations or sign conventions between regimes.

The math core is a pure-Rust crate ([`dynbem_rs/`](dynbem_rs/), no pyo3 /
numpy / file IO) wrapped by a thin PyO3 + maturin binding crate
([`dynbem/`](dynbem/)) which is the publishable Python package.

Two dynamic-inflow models are provided:

- **Pitt-Peters** (three-state global ν₀/ν_s/ν_c) — with the Peters
  L-matrix, Glauert wake-skew via the off-diagonal coupling, and the
  Leishman empirical VRS polynomial baked into the uniform-inflow state.
- **Øye 2-stage annular** — per-annulus filtered momentum inflow (the
  OpenFAST DBEMT formulation), independent across radii and numerically
  stable at high advance ratios where Pitt-Peters becomes stiff.

Both models share a tabulated polar interpolator and a common BEM ψ-loop
kernel ([`dynbem_rs/src/bem_common.rs`](dynbem_rs/src/bem_common.rs)) and
plug into the same `AeroModel` trait (Rust).
The repo also includes a flight-envelope sweep driver
(`envelope/compute_map.py`), a cyclic-trim solver
([`dynbem_rs/src/trim.rs`](dynbem_rs/src/trim.rs)), and a point-mass +
cyclic-pitch attitude simulator.
For empirical validation against published rotor data (Castles-Gray
TN-2474 vertical descent, Caradonna-Tung TM-81232 hover CT and
spanwise CL, Harrington TN-2318 full-scale hover, Wheatley & Hood
TR 515 forward-flight autorotation), see
[EMPIRICAL_VALIDATION.md](EMPIRICAL_VALIDATION.md).

Coordinates are NED throughout; rotor rotation is CCW-from-above
(American helicopter convention).

## Install

**Quick start** (all platforms):

```
./setup.sh           # POSIX shell or Windows git-bash/WSL
setup.cmd            # Windows cmd/PowerShell thin wrapper (calls bash internally)
```

Both scripts check prerequisites (Python 3.10+, cargo, C compiler) upfront,
create a virtual environment, install dependencies, and build the Rust extension
via maturin. See [setup.sh](setup.sh) for details on supported platforms.

**Manual setup** (if you prefer not to use the setup scripts):

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
model  = dynbem.create_aero(defn, model="pitt_peters")  # or "oye", "bem"
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
result, derivative = model.compute_forces(inputs, state)
# result.F_world, result.M_orbital, result.M_spin, result.Q_spin
# derivative carries d/dt of the dynamic-inflow states (lambda_0/c/s or W/W_int)
# Mechanical ODE lives in the caller:
#   from dynbem.mechanical import omega_derivative
#   omega += dt * omega_derivative(result.Q_spin, motor_torque_Nm, I_ode_kgm2)
```

For the full API reference — all classes, fields, keyword arguments, and
return types — see **[API.md](API.md)**.

## Flight envelope sweep

```
run_map.cmd                                           # quick grid, saves to out\map.npz, plots to out\
run_map.cmd --full --save out\map.npz --plot out\     # full grid
uv run python -m envelope.compute_map --help
```

## Tests

```
uv run pytest tests/ -q
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
[EMPIRICAL_VALIDATION.md](EMPIRICAL_VALIDATION.md).

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

## Coordinate system — NED

This project uses **NED (North-East-Down)** throughout, without exception:

- X = North, Y = East, Z = Down
- Gravity acts in the **+Z** direction
- Rotor thrust (upward lift) is **negative Z** in world frame: `F_world[2] < 0`
- Wind blowing upward (driving a flying turbine) is **negative Z** in world frame
- `R_hub` rotates from hub frame → NED world frame

### Reading literature — coordinate trap

Most helicopter and wind-turbine literature uses one of:
- **SAE / helicopter**: X forward, Y right, Z down (body frame, not world NED)
- **Wind-turbine (IEC 61400)**: X downwind, Y lateral, Z up (**ENU-like**)
- **Aeronautics (NED)**: X North, Y East, Z Down

When adapting equations or sign conventions from papers, always check
which frame the authors use. Windmill-brake-state and axial-induction
literature (Glauert, Buhl) often defines positive inflow **upward**
(opposing thrust), which is **negative Z** here. Flip signs accordingly.

### Inflow sign convention (NED)

For a rotor disk lying in the XY-plane (hub pointing down):

- `lambda` (inflow ratio) is positive when flow passes through the disk
  from above (downward, +Z direction), i.e. in **normal rotor mode**
  (helicopter hover).
- In **windmill / autorotation mode** the wind drives flow upward (−Z),
  so `lambda` is **negative** when the rotor is in energy-harvesting mode.
- Collective pitch `theta_0 > 0` pitches blade leading edge up
  (toward −Z thrust).

---

## BEM solver design — critical notes

### Hover-safe inflow iteration

The standard wind-turbine BEM uses the induction factor `a = v_i / V_inf`,
which **collapses to zero in hover** (`V_inf = 0`). This code instead
iterates on the **total inflow ratio** `λ_r = v_a / (Ω·R)`, where `v_a`
is the total axial velocity at the disk (external freestream + induced).

The combined momentum-BEM equation at each annulus is:

    k·(λ_r² + x²) = λ_r·(λ_r − λ_c)

where `k = σ_r·cn / (8·F)`, `x = r/R`, and `λ_c = v_climb / (Ω·R)`.

This quadratic is solved per iteration step; `v_climb = 0` in hover is
handled naturally (gives the standard hover solution
`λ_r = x·sqrt(k/(1−k))`).

### v_climb sign convention (internal BEM)

`v_climb = dot(v_rel_world, hub_axis_ned)` (no negation):

- `v_climb > 0`: air flows **downward** through disk (helicopter climb /
  normal inflow)
- `v_climb = 0`: hover
- `v_climb < 0`: air flows **upward** through disk (autorotation /
  flying wind turbine)

### Root selection in the momentum-BEM quadratic

The quadratic has two roots. Selection is by operating mode:

- Helicopter / hover (`λ_c ≥ 0`): take the **positive** root (`λ_r > 0`)
- Turbine / autorotation (`λ_c < 0`): take the **negative** root
  (`λ_r < 0`)

### Autorotation torque sign

In autorotation (upward wind, `λ_c < 0`):
- `λ_r < 0` → `φ < 0` → `ct = cl·sin(φ) − cd·cos(φ) < 0` → `Q_total < 0`
- `d_omega = (−Q_total + Q_motor) / I` → positive angular acceleration ✓

In powered/hover mode (`λ_c ≥ 0`):
- `Q_total > 0` (aerodynamic drag on rotor) → `d_omega < 0` without
  motor torque ✓

### Force direction

`F_world = −T_total · hub_axis_ned`

`T_total` is always positive for a rotor generating lift (cn > 0 in
both modes). With `hub_axis_ned = [0, 0, 1]` for a level rotor:
`F_world[2] = −T_total < 0` (upward). ✓

---

## Pitt-Peters design notes (`dynbem_rs/src/pitt_peters.rs`)

### State interpretation

`λ_0` (and `λ_c`, `λ_s`) is the **induced** inflow ratio `v_i / (ΩR)`,
not the total inflow. The total axial flow seen by each blade element
is:

    λ_total = λ_0 + λ_climb    where  λ_climb = v_climb / (ΩR) < 0 in descent

This must be computed inside the blade element loop — **do not pass
only `λ_0`**. Without the freestream term the blade never sees
net-upward flow in WBS, so CQ never goes negative and autorotation is
suppressed entirely.

### VRS polynomial sign convention

The Leishman (2000) polynomial uses descent-positive
λ₂ = V_descent / V_h:

    λ₁/V_h = 1 + 1.125·λ₂ − 1.372·λ₂² + 1.718·λ₂³ − 0.655·λ₂⁴

This is NOT the form with coefficients
(−1.125, −1.372, −1.718, −0.655), which applies when the argument is
V_climb/V_h (negative for descent). The two forms are equivalent;
this code uses descent-positive throughout.

### V_T floor

`V_T = |v_climb + v_0|` → 0 in the middle of VRS (upward freestream ≈
downward induced). A floor of `1e-2 · max(ΩR, 1)` prevents
`τ_0 → ∞` and division by zero. This is physically reasonable:
`τ_0 → large` in VRS is correct (slow, unsteady response), and the
exact floor value doesn't matter for stability.

### Why CT still rises in deep VRS

At λ₂ ≈ 1.5, the Leishman polynomial gives
`λ_0_ss ≈ 2 · V_h/ΩR`. Combined with `λ_climb ≈ −1.5 · V_h/ΩR`,
the net blade inflow `λ_total ≈ 0.5 · V_h/ΩR` is less than hover, so
AoA increases and CT rises. The real VRS has recirculating wakes that
further restrict net throughflow; the 1-D polynomial captures the mean
induced velocity but not the 3-D blockage. This is a known limitation
of all momentum-based VRS models.

---

## Øye design notes (`dynbem_rs/src/oye.rs`)

### State interpretation

`W[i]` and `W_int[i]` are induced inflow ratios `v_i / (Ω·R)` **per
annulus**, not global harmonics. The total axial flow at annulus `i`
seen by the blade is `λ_total[i] = λ_climb + W[i]` (compare with
Pitt-Peters' `λ_total = λ_climb + λ_0 + x·(λ_c·cos ψ + λ_s·sin ψ)`).

`W` is what the blade actually reads in the ψ-loop. `W_int` is the
intermediate filter stage between the quasi-steady target `W_qs[i]`
and `W`. Both arrays have length `n_elements`.

### Quasi-steady target

`W_qs[i]` is solved per annulus from Glauert momentum balance using
the rotor-mean `µ_T = V_T / Ω·R`:

    W_qs[i] = dCT/dx[i] / (4·x[i]·µ_T)
    where  V_T = √(v_edge² + (v_climb + v_0_mean)²)

This linear (in `W_qs`) form is what Pitt-Peters effectively uses in
its aggregate `λ_0_ss = T / (2ρA·V_T·ΩR)`. The pure axial-momentum
form `4·x·λ_r·W = dCT/dx` is unstable in forward flight (small λ_r in
descent makes W blow up) and was rejected during development.

### Why no L matrix

Annulus-local: each `W[i]` evolves independently, driven only by
`W_qs[i]` from its own annulus. Cross-annulus coupling happens only
through the rotor-mean `µ_T` in the τ formulas and `V_h` in the VRS
override. There's no analogue of Pitt-Peters' `−L_off·C_T` term that
feeds total thrust into the cyclic harmonics, so no BEM-driven
feedback loop and no associated stiffness — at the cost of not
modelling cyclic inflow harmonics at all.

### Time constants

`τ₁` is rotor-mean (depends on `a_avg`, not per-annulus); `τ₂(r)`
varies with radius. With `dt = 5 ms` and a 1 m rotor at `V_∞ ~ 10 m/s`,
`τ₁ ~ 0.1 s` and `τ₂ ~ 0.04 s` — both well above the envelope's outer
`dt`, so the semi-implicit Euler in `envelope/point_mass.py` is gentle
damping at most.

### Cyclic input

Cyclic pitch flows through the same `cyclic_coeffs` → `θ(ψ) =
collective + θ_1c·cos ψ + θ_1s·sin ψ` path as Pitt-Peters; the ψ-loop
produces correct hub moments. What's *missing* compared to
Pitt-Peters: the cyclic-driven hub moment doesn't develop a
counter-acting inflow harmonic (no `λ_c`/`λ_s` states), so the
steady-state moment is over-predicted vs Pitt-Peters at hover.
Cyclic *control* (sign and order-of-magnitude) is right; cyclic
*inflow damping* is absent.

---

## Research sources

Extracted tables and figures from primary literature live under
`Research/`. Each paper subfolder uses the convention
`page_NN_<description>.md` so extractions trace back to their source
page image.

- **CaradonnaTung/** — NASA TM-81232 (1981). 2-blade NACA 0012 hover
  CT data at θc = 5°/8°/12°. Primary BEM validation source. No CP /
  torque data.
- **Buhl_NREL_TP500_36834/** — NREL TP-500-36834 (2005). Windmill
  Brake State correction extending Glauert. Used for the WBS
  quadratic.
- **Castles_TN2474/** — NACA TN-2474 (Castles & Gray, 1951). Induced
  velocity in hover/descent — experimental basis for the Leishman VRS
  polynomial.
- **Harrington_TN2318/** — NACA TN-2318 (Harrington, 1951). Hover
  CT vs CP polars for two full-scale rotors. Candidate dataset for
  CP-CT polar validation.

## License

MIT — see `LICENSE`.
