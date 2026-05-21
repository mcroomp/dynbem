# CLAUDE.md — AI assistant instructions

Human-facing docs (install, usage, coordinate conventions, design notes,
implementation roadmap, research sources) live in [README.md](README.md).
**Read it first** — most of what you'd want to know about this codebase is
there, and this file does not repeat it.

Empirical validation — which papers/tables back each model, the
achieved variance vs published data, and why any residual bias exists
— lives in [EMPIRICAL_VALIDATION.md](EMPIRICAL_VALIDATION.md). Read
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
`dynbem/cyclic.py:cyclic_coeffs()` and goes through the rotor's
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

If a future model adds full flap-dynamics ODE, set φ ≈ +90° (rad
internally, deg in the YAML) so the user's `tilt_lon`/`tilt_lat`
command a *disk* tilt rather than a thrust asymmetry directly.

`control = None` defaults: gain = 1, φ = 0 → `tilt_lon, tilt_lat` are
direct blade-pitch amplitudes with helicopter-standard signs.

## Hub-frame aero moments

In the ψ-loop, each blade element contributes per-azimuth thrust `dT`
in the hub-axis (−Z hub) direction. With `r_pos = r·r_hat(ψ)` and
`F = dT · (−ẑ_hub)`:

    dM_hub = r_pos × F = r · dT · [sin(ψ), cos(ψ), 0]

i.e. `Mx_hub = Σ r·dT·sin(ψ)`, `My_hub = Σ r·dT·cos(ψ)` (averaged over
ψ in the model). These are then rotated to world via `R_hub` and
returned in `AeroResult.M_orbital`.

Coefficient form used in Pitt-Peters:

    C_T     = T_total / (ρ·A·(ΩR)²)
    C_L_hub = Mx_hub  / (ρ·A·(ΩR)²·R)     # rolling moment coefficient
    C_M_hub = My_hub  / (ρ·A·(ΩR)²·R)     # pitching moment coefficient

These differ in sign from BladeAD: BladeAD uses ψ=0 at the tail and
`dMy = −r·cos(ψ)·dT`, so when porting formulas from BladeAD:

- our λ_c = − BladeAD λ_c
- our λ_s = − BladeAD λ_s
- our C_M_hub = + BladeAD C_My
- our C_L_hub = − BladeAD C_Mx

Signs of `M_x → roll-right`, `M_y → pitch-up` follow the standard
NED body-frame right-hand rule (q = ω_y > 0 ⇒ nose pitches up).

## Pitt-Peters inflow ODE — as implemented

**Canonical reference**: David Peters' own Nikolsky Lecture, JAHS 54(1):011001
(2009), saved at [Research/Peters_Nikolsky_2008/](Research/Peters_Nikolsky_2008/)
with sign-translation notes. Eqs 7–11 of that paper define the model.

Peters' L matrix (his Eq 10, with X = tan(χ/2), state ordering (ν_0, ν_s, ν_c)):

    [L] = | 1/2          0          −15π·X/64 |
          | 0            2(1+X²)     0         |
          | 15π·X/64     0           2(1−X²)   |

with forcing `{C_T, −C_L, −C_M}` and his ψ=0 at the tail.

After translating to our ψ=0-at-+X convention (our λ_c = −ν_c, our λ_s = −ν_s)
and using `µ_T = √(µ² + λ_total²)` as the mass-flow scaling, the steady-state
targets are:

    µ_T   = √(µ² + λ_total²)            # mass-flow non-dim
    χ     = atan2(µ_inplane, |λ_total|) # wake skew angle
    L_off = (15π/64) · tan(χ/2)
    L_cc  = 4·cos(χ) / (1 + cos χ)      # = 2(1−X²), Peters Eq 10
    L_ss  = 4 / (1 + cos χ)             # = 2(1+X²), Peters Eq 10

    λ_0_ss = C_T/(2·µ_T)              +  L_off · C_M_hub / µ_T
    λ_c_ss = (−L_off · C_T            +  L_cc  · C_M_hub) / µ_T
    λ_s_ss = (                            L_ss  · C_L_hub) / µ_T

Time constants (Peters Eq 9 apparent mass `M = diag(8/(3π), 16/(45π), 16/(45π))`):
`τ_0 = 8R/(3π·V_T)`, `τ_cs = 16R/(45π·V_T)`.

The `−L_off · C_T / µ_T` term in `λ_c_ss` is the Pitt-Peters
cross-coupling — it produces Glauert wake-skew naturally from thrust
forcing. The closed-form Glauert tilt has been removed; do not re-add
it (would double-count).

The cross-coupling term `+L_off · C_M_hub / µ_T` in `λ_0_ss` shifts
uniform inflow in response to cyclic pitching moment — a higher-order
effect, small in practice but the formulation is symmetric.

VRS region (`v_climb < 0` and `0 < V_descent/V_h < 2`) still overrides
`λ_0_ss` with the Leishman empirical polynomial — momentum theory
doesn't apply in a recirculating wake, so the cross-coupling is also
skipped in that regime.

### Mass-flow parameter: µ_T vs Peters' V

We use `µ_T = √(µ² + λ_total²)` (classical Glauert). Peters uses
`V = (µ² + (λ+ν)(λ+2ν)) / √(µ² + (λ+ν)²)` (his Eq 8). They agree in
high-speed forward flight but differ by 2× in hover. Our µ_T reproduces
classical Glauert hover `λ_0 = √(C_T/2)`; Peters' V gives `√(C_T)/2`
(factor √2 different — possibly a C_T normalization convention in his
paper). The L-matrix STRUCTURE matches Peters exactly; only the scalar
scaling differs. Swapping to Peters' V would need validation against
hover data — defer until needed.

### Shared BEM infrastructure (`dynbem/_bem_common.py`)

Both `PittPetersModel(JIT)` and `OyeBEMModel` import from
`dynbem/_bem_common.py`:

- `vrs_lambda1` — Leishman VRS empirical polynomial
- `_interp_polar` — JIT polar lookup (`@njit`, called from inside both
  models' ψ-loop kernels)
- `build_polar_arrays` — one-time tabulation of any `AirfoilPolar` onto
  contiguous numba arrays
- `radial_grid` — one-time radial geometry caching

Hot-path kinematics (`Omega_R`, `hub_axis`, `v_climb`, `v_edge`,
`v_inplane_hub`) and `AeroResult` assembly are **deliberately inline**
in each model's `compute_forces` to avoid Python function-call
overhead in the per-step loop. Each is ~10 lines and identical
across models — duplicate it deliberately rather than abstract.

If you find yourself wanting to share more across models: the ψ-loop
kernels are the largest duplicated structure, but they have
model-specific `lam_local(r, ψ)` formulas and numba's `@njit`
doesn't compose closures cleanly. Don't try to unify them.

### Wind-axis rotation — NOT applied (limitation)

The L matrix is diagonal-plus-off-diagonal in **wind axes**, but the
current code treats `(C_M_hub, C_L_hub)` as if already in wind axes
(i.e. assumes in-plane wind along hub −X).  Exact for axial flight and
pure-longitudinal forward flight; approximate for oblique flight
`µ_y ≠ 0`.

An earlier implementation rotated forcing/inflow by
`β = atan2(v_in_hub_y, −v_in_hub_x)` and was rotationally covariant.
It was reverted because it destabilised the tethered-rotor envelope
(`envelope.compute_map`) at descent + edgewise wind operating points
via a nonlinear feedback `λ_c → BEM(lam_local) → C_L_hub → λ_s_ss`.
Implicit Euler on the λ states alone (now applied in
`envelope/point_mass.py`) is necessary but not sufficient — the BEM
loop's λ_c sensitivity also needs damping before the rotation can be
re-introduced.

## Øye 2-stage annular dynamic inflow (`dynbem/oye.py`)

`OyeBEMModel` is the **annulus-local** alternative to Pitt-Peters
implemented for the same project, with a deliberately different
state structure.

Per radial annulus `i`:

    τ₁ · dW_int[i]/dt + W_int[i] = W_qs[i] + k · τ₁ · dW_qs[i]/dt
    τ₂(r) · dW[i]/dt + W[i]     = W_int[i]

`W` is what the blade reads in the ψ-loop; `W_int` is the
intermediate filter stage between the momentum target `W_qs` and `W`.
`k = 0.6` (empirical, OpenFAST default). Treats `dW_qs/dt = 0`
across each outer step — DBEMT_Mod=1 equivalent.

`W_qs` per annulus from Glauert momentum (linear form):

    W_qs[i] = dCT/dx[i] / (4·x[i]·µ_T)

with rotor-mean `µ_T = √(µ² + (λ_climb + v_0_mean)²) / Ω_R`. The
pure axial-momentum form `4·x·λ_r·W = dCT/dx` was tried first and
was unstable in forward flight — see the docstring in
`dynbem.oye._solve_W_qs`.

### Why this exists alongside Pitt-Peters

Pitt-Peters couples C_T, C_M_hub, C_L_hub globally into all three
inflow harmonics via the L matrix → BEM-driven feedback that's stiff
at high advance ratios + descent. Øye's per-annulus filters are
independent → no global feedback → numerically stable in regimes
that need adaptive time-stepping with Pitt-Peters. This is exactly
the trade-off OpenFAST's DBEMT made.

### What Øye CAN'T do

- **No cyclic inflow harmonics**: there's no λ_c/λ_s state, so the
  inflow doesn't develop an asymmetric tilt in response to cyclic
  pitching/rolling moments. Cyclic *control* still works (hub moments
  respond correctly to `tilt_lon`/`tilt_lat`), but the cyclic
  *inflow feedback* that reduces steady-state moment in Pitt-Peters
  is absent. `tests/test_cyclic.py::test_cyclic_inflow_reduces_hub_moment`
  doesn't apply.
- **No wake-skew off-diagonal**: the `-L_off·C_T` term that produces
  Glauert wake skew from thrust forcing in Pitt-Peters has no
  analogue here. Wake skew has to come from the BEM ψ-loop's
  asymmetric loading alone.

### Sign conventions (same as Pitt-Peters)

- `W > 0` for hover / helicopter (induced flow downward through disk)
- `W > 0` in autorotation too (induction *slows* the upward freestream,
  but in the same NED-+Z direction it would push in helicopter mode).
  `λ_total[i] = λ_climb + W[i]` matches Pitt-Peters'
  `λ_climb + λ_0`.

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

- **Python**: use the venv at `.venv\`. Activate with
  `.venv\Scripts\activate` (Windows), or invoke directly via
  `.venv\Scripts\python` / `.venv\Scripts\pytest`. Don't install packages
  globally or create a new venv.
- **Shell**: always use the Bash tool with `.venv/Scripts/python` (forward
  slashes work fine on this Windows checkout). Do not switch to the
  PowerShell tool -- its quoting and Unicode handling have bitten this
  project's output (em-dashes render as the replacement glyph,
  `Select-Object` piping breaks on array args, etc.). If a one-liner is
  awkward in bash, write a short script under the appropriate dir and
  run it through bash instead.
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

- New levels (e.g. Peters-He, polar-grid BEM) plug in behind the
  `AeroBase` interface in `dynbem/__init__.py`. Don't break existing
  call sites — keep
  `compute_forces(inputs, state) -> (AeroResult, RotorState)`.
- Reuse `dynbem/_bem_common.py` for polar tabulation, radial-grid
  setup, the VRS polynomial, and the JIT polar interpolator.
  Hot-path kinematics and result assembly stay inline (see the
  "Shared BEM infrastructure" section above).
- Add `inflow_taus(inputs, state) -> np.ndarray` returning the time
  constant for each state component (`np.inf` for mechanical / quasi-
  static states). The envelope integrator's semi-implicit damping
  needs this; the default in `AeroBase` returns all-infinity, which
  is wrong for any dynamic-inflow model.
- Add a new `RotorState` subclass in `dynbem/rotor_state.py` with
  `to_array` / `from_array`. Convention: mechanical states `ω, ψ` are
  the **last two** entries — `arr[-2] = omega_rad_s`,
  `arr[-1] = spin_angle_rad`. The envelope's clipping and recovery
  code relies on this.
- Wire the new model into `create_aero` in `dynbem/__init__.py` with a
  stable string name. Add docs to the factory docstring.
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
- `Research/Peters_Nikolsky_2008/CLAUDE.md` — **canonical Pitt-Peters
  formulation** (L matrix, M matrix, V mass-flow, forcing sign
  convention) from David Peters' Nikolsky lecture. Read this before
  touching any Pitt-Peters signs or coefficients.
- `dynbem/CLAUDE.md` — public `dynbem` Python package (PyO3 glue + Python
  compat shim). Drop-in replacement for the legacy pure-Python dynbem.
- `dynbem_rs/CLAUDE.md` — pure-Rust math core (no pyo3 / numpy / file IO).
  Module map, hot-path conventions, numerical floors.
- `dynbem_old/` — read-only reference copy of the legacy pure-Python
  implementation. **Do not import or depend on it in new code.** It exists
  so the original algorithms remain readable side-by-side with the Rust port.

Defer to those when working inside the respective directories.
