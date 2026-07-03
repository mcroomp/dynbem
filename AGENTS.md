# AGENTS.md — AI assistant instructions

This is the primary instruction file for agent-based tooling in this repository.

## Instruction file policy

- Treat `AGENTS.md` as the single source of truth for agent instructions.
- When instruction text needs to be added, removed, or updated, edit
  `AGENTS.md` only.
- Do not add instruction content to `CLAUDE.md`; keep it as a redirect.

Human-facing docs (install, usage, coordinate conventions, design notes,
implementation roadmap, research sources) live in [README.md](README.md).
**Read it first** — most of what you'd want to know about this codebase is
there, and this file does not repeat it.

Empirical validation — which papers/tables back each model, the
achieved variance vs published data, and why any residual bias exists
— lives in [EMPIRICAL_VALIDATION.md](docs/EMPIRICAL_VALIDATION.md). Read
that before changing anything in the BEM / Pitt-Peters / Øye signs or
coefficients.

This file holds only directives that are specifically for you (the AI
assistant) and would be noise in the README.

## Rotor rotation direction — CCW from above (American convention)

**This project uses the American helicopter convention: the rotor spins
counter-clockwise (CCW) when viewed from ABOVE the rotor disk** (Bell /
Sikorsky / Boeing standard). Robinson, Eurocopter/Airbus, and most
European/Russian designs spin CW from above — **do not** use those as a
sign reference.

### What "from above" means in NED

The observer is physically above the rotor and looks *down* at it. In
NED, "above" is the −Z direction and the observer's line of sight is
along **+Z**. With +X=North up on the page and +Y=East to the right,
CCW-from-above traces N → W → S → E (i.e. +X → −Y → −X → +Y).

### Azimuth ψ and tangential direction (hub frame)

- ψ is the rotor azimuth, measured from the +X (hub forward / North-ish)
  axis, **increasing in the direction of blade motion** — so ψ increases
  CCW when viewed from above.
- Blade radial unit vector at ψ (hub frame):
  `r_hat = [cos(ψ), −sin(ψ), 0]`
- Blade tangential unit vector (direction the tip is moving) at ψ
  (hub frame):
  `t_hat = [−sin(ψ), −cos(ψ), 0]`
- `dψ/dt = ω > 0` in normal powered flight.

### v_t_extra (tangential apparent wind)

The tangential apparent wind at a blade element, from first principles:

    v_t = (v_blade − v_air) · t_hat = ω·r − v_inplane · t_hat

so `v_t_extra = −v_inplane · t_hat`. With the CCW t_hat above, in hub
frame components:

    v_t_extra = +v_in_hub_x · sin(ψ) + v_in_hub_y · cos(ψ)

### Things this convention determines

- **Advancing side** is on the **right** (+Y / East) of the hub in
  forward flight along +X. CCW from above ⇒ at the +Y position the
  blade is moving toward +X (into the relative wind), so +Y is the
  advancing side.
- **Lateral cyclic λ_s sign** and **cyclic phase lag** direction
  (≈90° phase lag in the direction of rotation, i.e. CCW from above).
- **Coriolis / gyroscopic moment signs** on the airframe — flipping
  rotation direction flips these.
- **Tail rotor anti-torque direction** (if/when one is modelled): main
  rotor CCW from above ⇒ reaction torque on fuselage is CW from above
  ⇒ tail rotor pushes the tail to the right (American convention).

If you ever need to change the rotation direction, do it in exactly
one place — the definition of `t_hat` — and re-derive the matching
`v_t_extra`, hub-moment, and Pitt-Peters L-matrix signs. Do not flip
the sign of ω.

## Cyclic pitch convention

Inputs `RotorInputs.tilt_lon`, `RotorInputs.tilt_lat` are **swashplate
tilt angles** (rad). The mapping to blade pitch lives in
`dynbem_rs/src/cyclic.rs::cyclic_coeffs()` (Rust core; exposed to
Python as `dynbem.cyclic_coeffs`) and goes through the rotor's
`ControlProperties.swashplate_pitch_gain_rad` (gain) and
`swashplate_phase_deg` (phase φ):

    θ(ψ) = collective + θ_1c·cos(ψ) + θ_1s·sin(ψ)

with

    θ_1c = gain · (−tilt_lon·cos φ − tilt_lat·sin φ)
    θ_1s = gain · (−tilt_lon·sin φ + tilt_lat·cos φ)

Sign convention is **helicopter-standard**:

- `tilt_lon > 0`  →  **nose-down** disk (forward stick)
- `tilt_lat > 0`  →  **roll right**

This assumes no flap dynamics — blade pitch directly sets local thrust,
no 90° precession. The mapping was derived for our ψ=0-at-+X (nose),
CCW-from-above convention: `tilt_lon > 0` peaks pitch at ψ=π (tail),
giving more thrust at the back → nose-down moment via the hub moment
integral below.

When the rotor uses `PitchActuation::ServoFlap`, both swashplate
collective and cyclic are interpreted as servo-flap commands and the
feathering response replaces the direct swashplate pitch path. The
feathering dynamics introduce an intrinsic 90 deg lag at 1/rev and the
solver does NOT compensate this internally. If axis-preserving control
behavior is needed, apply phase correction in the controller or via
swashplate phase configuration.

`control = None` defaults: gain = 1, φ = 0 → `tilt_lon, tilt_lat` are
direct blade-pitch amplitudes with helicopter-standard signs.

## Hub-frame aero moments

See [PITT_PETERS_DESIGN.md](docs/PITT_PETERS_DESIGN.md) for the full coefficient
definitions, the cross-product derivation, and the BladeAD sign differences
(our C_L_hub = -BladeAD C_Mx, our C_M_hub = +BladeAD C_My).

## Pitt-Peters inflow model

Implementation: `dynbem_rs/src/pitt_peters.rs`.
See [PITT_PETERS_DESIGN.md](docs/PITT_PETERS_DESIGN.md) for the full design:
L-matrix sign translation from Peters' Nikolsky lecture, state interpretation,
mass-flow parameter choice, wind-axis rotation, and VRS notes.
**Read it before touching any Pitt-Peters signs or coefficients.**

## Shared BEM infrastructure (`dynbem_rs/src/bem_common.rs`)

`QuasiStaticBEM`, `PittPetersModel`, and `OyeBEMModel` (in
`dynbem_rs/src/quasi_static_bem.rs`, `pitt_peters.rs`, `oye.rs`) all
delegate to the helpers in
`dynbem_rs/src/bem_common.rs`:

- `PolarTable` — contiguous-array polar tabulation.
- `RadialGrid` — one-time radial geometry caching (r_mid, x_mid, chord,
  twist per station).
- `vrs_lambda1` (in `dynbem_rs/src/common.rs`) — Leishman VRS polynomial.
- `kinematics()` → `Kinematics` — once-per-call hub-frame setup
  (`omega_r`, `hub_axis`, `v_climb`, `v_inplane`, `v_edge`,
  `v_inplane_hub`, `mu`). Identical across all three models.
- `vrs_regime()` → `VrsRegime` — `(v_h, lam2, in_vrs)` from
  `(T, v_climb, ρ, A)`. Shared by Pitt-Peters and Øye.
- `assemble_result()` — builds `AeroResult` (F_world, M_orbital, Q_spin,
  M_spin) from `(T, Q, Mx_hub, My_hub)` + hub axes. Used by all three.
- `element_force()` — `#[inline(always)]` per-element BEM integrand
  returning `(dT, dQ)` given `(v_a, v_t, col_psi, twist, …, polar)`.
  Used by every radial inner loop.
- `PsiKernel` trait + `run_psi_loop()` — the single ψ × r kernel used
  by Pitt-Peters and Øye. Each model implements `PsiKernel` for its own
  `lam_local(i, cos ψ, sin ψ)` formula and (optionally) the
  `on_element` per-element callback. Monomorphized over `K: PsiKernel`
  with `#[inline(always)]` on the trait methods, so codegen is
  identical to a hand-rolled loop — there is no `dyn` here, the trait
  is used as a *static interface*.

Reach for these helpers when adding a new model rather than duplicating
the math. The earlier guidance in this section ("don't unify the
ψ-loop kernels") was written before the helpers existed and reflected
a worry about closure/indirection overhead that doesn't apply to monomorphized
Rust generics — empirical timing (see `dynbem/benchmarks/bench_rust_only.py`)
confirms zero perf cost from the trait abstraction.

## Oye dynamic inflow model

Implementation: `dynbem_rs/src/oye.rs`.
See [OYE_DESIGN.md](docs/OYE_DESIGN.md) for the full design: per-annulus filter
state interpretation, W_qs momentum target, why the axial form was rejected,
and what Oye cannot model.
**Read it before touching any Oye filter parameters or signs.**

## Kaman servo-flap modeling (Beaupoil rotor)

Blade pitch actuation is modelled by the `PitchActuation` enum (in
`dynbem_rs/src/rotor_definition.rs`), with two variants:

- `DirectMechanical` (default): the swashplate sets blade pitch directly.
- `ServoFlap(ServoFlapActuation)`: a trailing-edge servo-flap drives a
  passive feathering DOF; the feathering response replaces the direct
  swashplate pitch path.

Servo-flap forcing is active through the servo-flap actuation path:

- Geometry / parameters: `ServoFlapActuation` (mechanical: inertia,
  damper, AC offset) holding a `ServoFlapGeometry` (flap C_M_delta and
  span limits) in `dynbem_rs/src/rotor_definition.rs`
- Dynamics + solve: `dynbem_rs/src/servoflap.rs`
- Call sites: `dynbem_rs/src/pitt_peters.rs`, `dynbem_rs/src/oye.rs`,
  `dynbem_rs/src/quasi_static_bem.rs`

Model scope (current):

- Quasi-static 1/rev harmonic feathering solve
- Mechanical pitch-bearing damping
- Optional aerodynamic spring from AC offset
- Servo-flap aerodynamic pitching-moment forcing
- Servo mode path split: in `PitchActuation::ServoFlap`, both collective
  and cyclic are interpreted as flap commands and direct
  swashplate-to-blade pitch is disabled

Known limitations:

- No direct sectional `dCL/d_delta * delta_f` lift increment yet
- No servo actuator lag model yet
- DC flap path uses quasi-static approximation; for AC-on-axis configurations
  (`k_aero ~= 0`) it falls back to collective pass-through to retain authority
- `mu` cross-term in flap forcing uses scalar advance ratio only

The old `kaman_flap:` block under `control:` remains metadata-only and
is not consumed by the Rust aero solvers. Use the `pitch_actuation:`
YAML block (with a `servoflap:` sub-block) to enable active servo-flap
feathering dynamics.

## Quasi-static blade flapping (hub moment reduction)

`FlapProperties` (in `dynbem_rs/src/rotor_definition.rs`) models
out-of-plane blade flexibility as an equivalent spring-hinge. The blade
absorbs most aerodynamic pitching/rolling moment via deflection; only a
fraction reaches the airframe (hub).

Parameters:

- `I_blade_flap_kgm2` -- blade flap inertia about virtual hinge [kg*m^2]
- `omega_nr_rad_s` -- non-rotating flap natural frequency [rad/s]
  (K_beta = I_b * omega_NR^2). 0.0 = freely hinged (no spring).

Physics:

    nu_beta^2 = 1 + (omega_NR / Omega)^2
    hub_moment_factor = (nu_beta^2 - 1) / nu_beta^2

- Freely hinged (omega_NR=0): factor=0, no moment transfer.
- Rigid blade (omega_NR >> Omega): factor->1, full moment transfer.
- Typical hingeless rotor: nu_beta ~ 1.05-1.15, factor ~ 0.05-0.15.

Implementation:

- `apply_flap_reduction()` in `dynbem_rs/src/bem_common.rs` scales
  `mx_hub, my_hub` by the factor before `assemble_result`.
- Applied in all three models (Pitt-Peters, Oye, quasi-static BEM).
- The inflow ODE (Pitt-Peters) still uses full aerodynamic moments
  (the wake responds to disk loading, not what the airframe sees).
- Thrust and torque are unchanged by flapping.

YAML schema:

    flap:
      I_blade_flap_kgm2: 0.012
      omega_nr_rad_s: 8.0

Absent `flap:` section = rigid blade (no moment reduction), preserving
backward compatibility.

Tests: `tests/test_flap_hinge.py`.

## Windmill solver: non-axial v_t_extra extension

The quasi-static BEM (`dynbem_rs/src/quasi_static_bem.rs`) uses a Ning
2014 Brent's-method windmill solver (Section 9.1 of docs/BEM_COMMON.md) for energy-
extracting elements. This solver was originally derived for pure axial
(wind-turbine) flow. Our extension makes it work inside the azimuth-
resolved psi-loop with in-plane wind -- something OpenFAST's AeroDyn
BEMT does not do (it only runs the windmill path in the axial-flow code).

### The fix

Each blade element at azimuth psi has tangential velocity
`v_t = Omega*r + v_t_extra(psi)` where `v_t_extra` comes from forward
flight (Section 6 of docs/BEM_COMMON.md). The windmill Brent residual must include
this term:

    g(phi) = sin(phi) * (1+a') * lam_tilde + cos(phi) * (1-a) = 0

with `lam_tilde = (Omega*r + v_t_extra) / |v_climb|` instead of the
axial-only `lam_r = Omega*r / |v_climb|`.

Without this correction, when `|v_t_extra|` is large (high advance ratio
+ descent), the solver finds spurious roots near phi=0 where a->1. The
old workaround was an `allow_windmill: v_edge < |v_climb|` threshold that
disabled the windmill solver entirely in oblique descent -- effective but
overly conservative, preventing correct windmill operation at those
azimuths where the bracket does exist.

### Bracket-existence as the natural filter

The corrected solver checks whether the residual changes sign across
phi in (-pi/2, -epsilon). If both endpoints have the same sign, no
windmill root exists for that element/azimuth and the helicopter
quadratic is used instead. This is physics-based: the bracket closes
naturally on the retreating side where v_t_extra opposes Omega*r, and
opens on the advancing side where the element is genuinely windmilling.

### Why this matters

Autorotating kites, autogyros in forward flight, and any rotor
descending at oblique angles need the windmill solver active in the
psi-loop. Without the v_t_extra correction these operating points
produce thrust sign flips at certain rotor speeds (the helicopter
quadratic picks the wrong root branch). With it, the crossover from
helicopter to windmill mode is smooth and azimuthally resolved.

Commits: `acd3183`, `9249e2b`. Tests: `tests/test_rawes_ic_aero_sign.py`.

## Do not revert work without explicit instructions

If a test fails, a build breaks, or run_map blows up, **do not respond
by deleting or reverting the code that produced the failure unless the
user has told you to**. The first move is to understand *why* — read
the code, instrument, reason about it. Reverting silently throws away
work the user has chosen to keep, and the failure is usually fixable
in place.

Examples that are NOT a license to revert: a single test failure, a
"this used to work" report, a regression you introduced, your own
prior edits looking wrong in hindsight. Examples that ARE: the user
says "revert it", "drop that change", "go back to X". When uncertain,
ask before reverting.

This rule has bitten before — see [memory feedback-no-silent-reverts].

## Workflow

- **Temporary output files**: when a script, command, or tool needs to
  write scratch output (plots, CSV dumps, debug logs, pytest redirects,
  etc.) always write it under `tmp/` in the repo root. That directory is
  git-ignored. Never write temp files directly into the repo root or any
  tracked directory.

- **Python**: this repo is a uv workspace + Cargo workspace. `uv sync`
  from the repo root builds the Rust extension (via maturin) and
  installs `dynbem` editable; `uv sync --group dev` also pulls pytest +
  maturin + build + twine. Run anything Python through `uv run` (e.g.
  `uv run pytest tests/ -q`, `uv run python -m envelope.compute_map`).
  Don't create a `.venv\` by hand or `pip install` globally -- uv owns
  the environment.
- **Rust**: the math core is `dynbem_rs/` (pure Rust, no pyo3 / numpy /
  file IO). The PyO3 + maturin glue is `dynbem/`. `cargo test
  --workspace` runs the Rust unit tests; the authoritative regression
  suite is `uv run pytest tests/ -q`, which exercises the full
  Rust-backed Python API.
  - **VPM tests must be run in release mode.** The VPM rotor marches
    hundreds of steps per test case and is ~50-100x slower in debug.
    Always use `cargo test --release -p dynbem_rs -- vpm` (or
    `cargo test --release --workspace`) when running or adding tests
    that exercise `VpmRotor`. Running VPM tests in debug is not wrong
    but takes 10+ minutes and is not practical as part of a normal
    edit-test loop.
- **Shell**: always use the Bash tool. Do not switch to the PowerShell
  tool -- its quoting and Unicode handling have bitten this project's
  output (em-dashes render as the replacement glyph, `Select-Object`
  piping breaks on array args, etc.). If a one-liner is awkward in
  bash, write a short script under the appropriate dir and run it
  through bash instead.
- **CRITICAL -- ASCII only in new Python / CSV / Markdown content.** No
  Greek letters, no em-dashes, no degree signs, no smart quotes, no
  subscripts/superscripts, no plus-minus or less-equal glyphs. Use plain
  ASCII transliterations: `theta`, `Omega`, `psi`, `lambda`, `sigma`,
  `mu`, `deg`, `<=`, `+/-`, `--`, `"..."`. The Windows console codepage
  mangles non-ASCII output (em-dash renders as a replacement character),
  `extract_tables.py` transliterates everything for the CSV mirror
  anyway, and grep / sed / diff are noticeably less reliable on
  mixed-encoding text. Applies to source code, string literals, print
  output, comments, docstrings, CSV cells, and Markdown bodies. Existing
  non-ASCII content in legacy docstrings and Research/ table titles can
  stay until it is edited for another reason; do not introduce new
  instances.
  - **Exception -- design documents may use GitHub Flavored Markdown
    math.** For standalone design / documentation Markdown (e.g.
    `*_DESIGN.md`) that is meant to be read rendered on GitHub, GFM math
    markup (`$...$` inline and `$$...$$` or ` ```math ` blocks) is
    allowed and encouraged where it makes the equations more readable.
    This carve-out is only for the math notation inside such docs -- the
    ASCII rule still applies to source code, string literals, print
    output, comments, docstrings, and CSV cells. Keep the surrounding
    prose ASCII (no em-dashes or smart quotes); the exception covers the
    LaTeX math spans, not the whole file.
  - **GFM math traps to avoid in design docs.**
    - `\operatorname{...}` is not supported by GitHub's math renderer.
      Use `\mathrm{...}` instead (e.g. `\mathrm{atan2}`).
    - Multiple `_` subscripts in a single inline `$...$` span confuse
      GitHub's Markdown parser, which eats the underscores as italic
      markers before the math renderer runs. Any inline expression with
      two or more `_` subscripts (e.g.
      `$\mathbf{u}_\text{ind} = \mathbf{u}_\text{far}$`) should be
      promoted to a display block (`$$...$$` or ` ```math `).
      A single subscript in inline math is fine; the problem only
      triggers when there are two or more `_` in the same `$...$` span.
    - `\text{...}` and `\mathrm{...}` are both supported. Prefer
      `\text{...}` for word-labels inside equations and `\mathrm{...}`
      for function names (atan2, sin, etc.).
    - Do NOT put underscore-containing identifiers (e.g. `tilt_lon`)
      inside `\text{...}` or `\mathrm{...}`. GitHub's validator rejects
      `_` inside text-mode commands (`'_' allowed only in math mode`).
      Instead, end the `$...$` span before the identifier and write it
      as an inline code span: `$\theta_{1c} = -$\`tilt_lon\``.
- **Coordinate frame**: NED everywhere. See README "Coordinate system" —
  the "coordinate trap" section especially matters when you adapt
  equations from a paper, because most rotor literature uses a different
  frame and the sign flips are easy to miss.
- **Sign conventions**: before changing any inflow / thrust / torque
  sign, re-read the README "BEM solver design" and "Pitt-Peters design
  notes" sections. The signs are load-bearing and were tuned to make
  hover, climb, descent, VRS, and autorotation all work in one code
  path.
- **Validation tests pair with `verification/` scripts.** When a test
  in `tests/` checks the model against a published dataset, do not
  duplicate the BEM-call + comparison loop inside the test file.
  Instead: factor the loop into `verification/<paper>_<quantity>.py`
  as an importable function that accepts a `sample` argument, then
  have both the script's `main()` and the unit test call it. The unit
  test runs with a small `sample` (fast, fits in the pytest budget)
  and asserts on the returned aggregate; the verification script with
  no sample is the authoritative whole-dataset sweep used to
  re-baseline bounds. This keeps the per-test BEM-driver logic in
  exactly one place and prevents the spot tests and the survey from
  drifting apart.

## When extending the aero models

New aero models live in the Rust core (`dynbem_rs/`) and are exposed
to Python via pyo3 wrappers in `dynbem/src/wrappers.rs`. The full
recipe is in [`dynbem_rs/CLAUDE.md`](dynbem_rs/CLAUDE.md) ("Adding a
new aero model"); the short version:

- Implement the `AeroModel` trait in `dynbem_rs/src/aero_model.rs` for
  the new struct (`fn compute_forces(&self, inputs, state) -> (AeroResult, RotorState)`).
  Don't break existing call sites.
- Reuse `dynbem_rs/src/bem_common.rs` (`PolarTable`, `RadialGrid`) and
  `dynbem_rs/src/common.rs` (`vrs_lambda1`, the numerical floors).
  Hot-path kinematics and result assembly stay inline (see the
  "Shared BEM infrastructure" section above).
- Add `inflow_taus(inputs, state) -> Vec<f64>` (via `RotorStateExt`)
  returning the time constant for each state component (`f64::INFINITY`
  for quasi-static states). The envelope integrator's
  semi-implicit damping needs this.
- Define the model-specific state struct in that model's aero module and add
  the matching `RotorStateExt` impl there (`quasi_static_bem.rs`,
  `pitt_peters.rs`, or `oye.rs`). The `RotorStateExt` trait itself is
  declared in `dynbem_rs/src/aero_model.rs`.
  `RotorStateExt` serializes **inflow states only** via
  `get_inflow()` / `set_inflow(Vec<f64>)`.  There are no mechanical
  fields in any state struct — `omega_rad_s` is passed by the caller
  through `RotorInputs` every call, and the caller advances omega via
  `dynbem.mechanical.omega_derivative` externally.
- Add a `PyFoo` newtype in `dynbem/src/wrappers.rs` and an `AeroAny`
  variant in `dynbem/src/trim_py.rs`. Wire the new model into
  `create_aero` in `dynbem/python/dynbem/factory.py` with a stable
  string name.
- Validation data lives under `Research/`. Add a
  `tests/test_<model>.py` and, if appropriate, a `val_step*.py`
  script that compares against a specific paper's data.
- Don't store derived results inside `Research/` — that directory is
  for source-paper extractions only.
- **`Research/extract_tables.py`** converts every markdown table under
  `Research/` into an ASCII CSV under `Research/csv/<Paper>/…`,
  mirroring the source folder structure (Greek letters and subscripts
  transliterated to plain names — `mu`, `alpha`, `theta0.75R`,
  `DeltaCQ`, etc.). The script is idempotent: run it whenever a
  table extraction is updated to keep `Research/csv/` in sync. Tests
  can import the CSVs directly (e.g. with `numpy.genfromtxt` or
  `csv.DictReader`) instead of re-parsing the markdown. If a paper
  introduces a Unicode character not yet handled, add it to
  `_ASCII_MAP` in the script and re-run — the run prints a warning for
  any character it falls back to dropping.
- **The markdown table is the single source of truth.** If you find an
  error in a generated CSV (bad cell, wrong column name, scrambled
  row), **fix it in the source `.md` file under `Research/<Paper>/`**
  and then re-run `python Research/extract_tables.py` to regenerate
  the CSVs. Never edit a file under `Research/csv/` directly — those
  files are derived output and will be overwritten on the next run.
  When in doubt about a value, re-check the high-resolution `.png` of
  the source page that the `.md` cites (the file the extraction was
  originally made from); correct the `.md`, regenerate, and verify
  the consistency-check assertions in any related test still pass.

## Subfolder CLAUDE.md files

- `Research/CLAUDE.md` — extraction conventions for paper sources and
  the `extract_tables.py` MD→CSV converter described above.
- `Research/CaradonnaTung/CLAUDE.md` — Caradonna-Tung page index, CT
  tables, validation notes.
- `Research/Peters_Nikolsky_2008/CLAUDE.md` — canonical Pitt-Peters
  formulation (L matrix, M matrix, V mass-flow, forcing sign convention)
  from David Peters' Nikolsky lecture.
- `docs/PITT_PETERS_DESIGN.md` — implementation design: sign translation,
  state interpretation, mass-flow choice, wind-axis rotation, VRS notes.
  **Read before touching Pitt-Peters.**
- `docs/OYE_DESIGN.md` — implementation design: per-annulus states, W_qs
  momentum target, what Oye cannot model. **Read before touching Oye.**
- `dynbem/CLAUDE.md` — public `dynbem` Python package (PyO3 glue + Python
  compat shim). Drop-in replacement for the legacy pure-Python dynbem.
- `dynbem_rs/CLAUDE.md` — pure-Rust math core (no pyo3 / numpy / file IO).
  Module map, hot-path conventions, numerical floors.

Defer to those when working inside the respective directories.
