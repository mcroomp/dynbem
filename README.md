# dynbem

**Dynamic blade-element momentum rotor aerodynamics — helicopter and
wind-turbine modes in one code path.**

`dynbem` is a Python rotor-aerodynamics library built around a multi-element
blade-element-momentum (BEM) solver coupled to dynamic-inflow models. It is
designed to be numerically valid across the **full operating envelope** —
helicopter hover, axial climb, axial descent, vortex-ring state (VRS),
windmill-brake state (WBS), autorotation, and wind-turbine power extraction
— without switching equations or sign conventions between regimes.

Two dynamic-inflow models are provided:

- **Pitt-Peters** (three-state global ν₀/ν_s/ν_c) — Numba-JIT-compiled,
  with the Peters L-matrix, Glauert wake-skew via the off-diagonal
  coupling, and the Leishman empirical VRS polynomial baked into the
  uniform-inflow state.
- **Øye 2-stage annular** — per-annulus filtered momentum inflow (the
  OpenFAST DBEMT formulation), independent across radii and numerically
  stable at high advance ratios where Pitt-Peters becomes stiff.

Both models share a tabulated polar interpolator and a common BEM ψ-loop
kernel, and they plug into the same `AeroBase` interface. The package also
includes a flight-envelope sweep driver (`envelope/compute_map.py`), a
cyclic-trim solver, and a point-mass + cyclic-pitch attitude simulator.
For empirical validation against published rotor data (Castles-Gray
TN-2474 vertical descent, Caradonna-Tung TM-81232 hover CT and
spanwise CL, Harrington TN-2318 full-scale hover, Wheatley & Hood
TR 515 forward-flight autorotation), see
[EMPIRICAL_VALIDATION.md](EMPIRICAL_VALIDATION.md).

Coordinates are NED throughout; rotor rotation is CCW-from-above
(American helicopter convention).

## Install

The package is set up with a standard `pyproject.toml` and can be installed
into any environment:

```
pip install -e .
```

For the bundled `.venv` on Windows, run `setup.cmd` from the repo root —
it creates `.venv\`, upgrades pip, and installs `requirements.txt`.
Activate with `.venv\Scripts\activate`, or invoke directly via
`.venv\Scripts\python` / `.venv\Scripts\pytest`.

## Usage

```python
import numpy as np
import dynbem

defn   = dynbem.rotor_definition.load("rotors/castles_gray_6ft/rotor.yaml")
model  = dynbem.create_aero(defn, model="pitt_peters_jit")  # or "pitt_peters", "oye", "bem"
state  = model.initial_rotor_state()

inputs = dynbem.RotorInputs(
    collective_rad=0.14,
    tilt_lon=0.0, tilt_lat=0.0,                # swashplate (helicopter-standard signs)
    R_hub=np.eye(3),
    v_hub_world=np.zeros(3),
    wind_world=np.zeros(3),
    t=0.0,
)
result, derivative = model.compute_forces(inputs, state)
# result.F_world, result.M_orbital, result.M_spin, result.Q_spin
# derivative is a RotorState with d/dt of every state field
```

`dynbem/__init__.py` lists the full public surface — models (`BEMModel`,
`PittPetersModel`, `PittPetersModelJIT`, the `create_aero` factory),
inputs/outputs (`RotorInputs`, `AeroResult`), state types
(`QuasiStaticRotorState`, `PittPetersRotorState`), polars (`AirfoilPolar`,
`LinearPolar`), rotor-definition types (`RotorDefinition`,
`BladeGeometry`, `AirfoilProperties`, `ControlProperties`, etc.), and the
`vrs_lambda1` helper. Cyclic mapping (`tilt_lon`/`tilt_lat` →
blade-pitch coefficients) lives in [dynbem/cyclic.py](dynbem/cyclic.py).

## Flight envelope sweep

```
run_map.cmd                  # quick grid, saves to out\map.npz, plots to out\
run_map.cmd --full           # full grid
python -m envelope.compute_map --help
```

## Tests

```
.venv\Scripts\pytest
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

## Implementation roadmap

The model is built in phases from simple to state-of-the-art, each a
drop-in upgrade behind the same `AeroBase` interface.

### Level 1 — Multi-element quasi-static BEM ✅ DONE

- Multi-element BEM loop (radial quadrature over `n_elements` annuli)
- Hover-safe inflow iteration on `λ_r` (not wind-turbine induction
  factor `a`)
- Per-element Prandtl **tip + hub** loss `F = F_tip · F_hub` (both
  factors exported from `dynbem.bem`)
- Glauert / Buhl Windmill Brake State correction (quadratic root
  selection)
- Forward-flight ψ-loop: per-azimuth blade pitch (cyclic), tangential
  wind projection (advancing/retreating), and in-plane hub moment
  accumulation (M_orbital)
- `dω/dt = (Q_aero + Q_motor) / I_ode` — rotor speed integrated as ODE
  state
- `dψ/dt = ω` — spin angle integrated as ODE state
- Returns `QuasiStaticRotorState` derivative
- **Validation**: see [EMPIRICAL_VALIDATION.md](EMPIRICAL_VALIDATION.md).

### Level 2 — Pitt-Peters 3-state dynamic inflow ✅ DONE

- `PittPetersModel` (numpy) and `PittPetersModelJIT` (Numba-compiled, same
  physics) in `dynbem/pitt_peters.py` / `dynbem/pitt_peters_jit.py`
- Prescribed-inflow blade element loop with per-element Prandtl tip + hub
  loss; blade sees `λ_total = λ_0 + v_climb/ΩR` (induced state +
  freestream), so WBS and autorotation work correctly
- Pitt-Peters ODE in matrix form (Peters 2009 Eq 7, hub axes):
  `[M] dλ/dt + V·[L]⁻¹ λ = forcing` with
  `M = diag(8/(3π), 16/(45π), 16/(45π))` → `τ_0 = 8R/(3πV_T)`,
  `τ_cs = 16R/(45πV_T)`.
- Steady-state targets follow the canonical L matrix (Peters Eq 10) with
  X = tan(χ/2), translated to our ψ=0-at-+X convention:
  ```
  λ_0_ss = C_T/(2·µ_T)        + (15π·X/64) · C_M_hub / µ_T
  λ_c_ss = −(15π·X/64) · C_T  + 4·cos(χ)/(1+cos χ) · C_M_hub) / µ_T
  λ_s_ss =                      4/(1+cos χ)         · C_L_hub  / µ_T
  ```
  The `−(15π·X/64)·C_T` cross-coupling in λ_c_ss is the Pitt-Peters term
  that produces Glauert wake-skew naturally from thrust forcing — no
  closed-form Glauert tilt needed.
- Cyclic input (`tilt_lon`, `tilt_lat`) wired through both models:
  blade pitch `θ(ψ) = collective + θ_1c·cos(ψ) + θ_1s·sin(ψ)` with
  helicopter-standard signs (`tilt_lon > 0` → nose-down,
  `tilt_lat > 0` → roll right). See `dynbem/cyclic.py` and CLAUDE.md.
- In-plane hub moments returned via `AeroResult.M_orbital`
  (`Mx_hub, My_hub` accumulated in the ψ-loop) — needed for cyclic to
  produce vehicle attitude response in the outer loop.
- **VRS empirical correction** (Leishman 2000, fit to Castles-Gray
  data): in 0 < λ₂ < 2,
  `λ_0_ss` comes from the polynomial
  `λ₁/V_h = 1 + 1.125λ₂ − 1.372λ₂² + 1.718λ₂³ − 0.655λ₂⁴`
  rather than momentum theory, preventing the Level-1 CT blow-up in VRS.
  Cross-coupling is also skipped in the VRS regime.
- **Canonical reference**: Peters, D.A. (2009), "How Dynamic Inflow
  Survives in the Competitive World of Rotorcraft Aerodynamics: The
  Alexander Nikolsky Honorary Lecture," *JAHS* 54(1):011001. PDF and
  extraction notes in `Research/Peters_Nikolsky_2008/`.
- **Validation**: see [EMPIRICAL_VALIDATION.md](EMPIRICAL_VALIDATION.md).
- **Known limitations**:
  - VRS CT still rises to ~2× nominal in deep VRS (λ₂ ≈ 1.5–2) at
    fixed θ; real rotor stays near nominal (paper: θ barely adjusts).
    The Leishman polynomial shifts the operating point but doesn't
    fully suppress it.
  - Autorotation torque crossing at V/ΩR ≈ 0.14 vs paper's 0.083.
  - Mass-flow scaling uses `µ_T = √(µ²+λ²)` (classical Glauert) rather
    than Peters' Eq 8 `V = (µ²+(λ+ν)(λ+2ν))/√(µ²+(λ+ν)²)`. They agree
    in high-speed forward flight but differ by 2× in hover. Switching
    would need validation against a hover dataset.
  - Wind-axis rotation of the L-matrix is NOT applied; oblique flight
    `µ_y ≠ 0` is approximate.  Exact for axial and pure-longitudinal
    flight.  A previous implementation was reverted because it
    destabilised the tethered-rotor envelope — see CLAUDE.md.

### Level 2 alt — Øye 2-stage annular dynamic inflow ✅ DONE

- `OyeBEMModel` in `dynbem/oye.py` (Numba-compiled ψ-loop)
- **Annulus-local** inflow: each radial annulus has its own pair of
  first-order lag filters `(W_int, W)` chasing the quasi-steady
  momentum target `W_qs`.  No global L-matrix; no λ_c/λ_s harmonic
  states.
- Two time constants per annulus (Øye 1990, OpenFAST AD Theory §6.3.4):
  ```
  τ₁ = 1.1 / (1 − 1.3·min(a, 0.5)) · R / V_∞
  τ₂(r) = (0.39 − 0.26·(r/R)²) · τ₁
  τ₁·dW_int/dt + W_int = W_qs + k·τ₁·dW_qs/dt
  τ₂·dW/dt + W = W_int
  ```
  with empirical coupling `k = 0.6`.  DBEMT_Mod=1 equivalent
  (`dW_qs/dt = 0` across each outer step — exact for envelope sweeps).
- W_qs per annulus from Glauert momentum balance using rotor-mean
  `µ_T = V_T / Ω·R`:  `W_qs[i] = dCT/dx[i] / (4·x[i]·µ_T)`.
- Same VRS override (Leishman polynomial) as Pitt-Peters for
  `0 < V_descent/V_h < 2` — applied uniformly across annuli.
- Same cyclic-pitch wiring (`tilt_lon` / `tilt_lat` → per-ψ blade
  pitch) and same in-plane hub moments returned via `M_orbital`.
- **Why this alongside Pitt-Peters**: Pitt-Peters' L-matrix couples
  thrust + hub moments back into all three inflow harmonics
  globally, which produces a stiff BEM-driven feedback at high
  advance ratios and in descent + edgewise wind.  Øye's annulus-local
  filters are independent → no feedback loop → numerically stable in
  the same regimes that needed adaptive time-stepping with
  Pitt-Peters.  OpenFAST's DBEMT uses the same Øye-style formulation
  for this reason.
- **Trade-off**: no harmonic inflow states means the inflow doesn't
  develop a `λ_c`-like tilt in response to cyclic pitching moments,
  so `tests/test_cyclic.py::test_cyclic_inflow_reduces_hub_moment`
  (which checks PP's specific feedback mechanism) doesn't apply.
  Cyclic *control* still works (hub moments respond correctly to
  swashplate inputs), but cyclic *inflow feedback* is absent.
- **Validation**: see [EMPIRICAL_VALIDATION.md](EMPIRICAL_VALIDATION.md).
- **References**:
  - Øye, S. (1990).  A simple vortex model.  IEA Symposium.
  - Snel, H. & Schepers, J.G. (1995).  Joint investigation of dynamic
    inflow effects.  ECN.
  - OpenFAST AeroDyn Theory v3.5, §6.3.4 (DBEMT).

### Level 3 — Peters-He finite-state dynamic inflow (state of the art)

- 9-state (or higher-order) Peters-He inflow model
- Requires new `PetersHeRotorState` dataclass
- Captures higher harmonics of the inflow distribution
- Best accuracy for maneuvering flight and aeroelastic coupling
- **Validation**: see [EMPIRICAL_VALIDATION.md](EMPIRICAL_VALIDATION.md).

### Forward flight (applies to all levels) — implemented

- Oblique inflow: advance ratio `µ = V_edge / (Ω·R)` ≠ 0
- Blade azimuth-dependent velocity in the BEM loop (`n_psi=36` stations
  by default, triggered when `µ > 0.01`, cyclic input is nonzero, or
  cyclic inflow state is nonzero)
- In-plane hub moments `Mx_hub`, `My_hub` returned in `AeroResult.M_orbital`
- Pitt-Peters L matrix off-diagonal `−L_off·C_T` produces Glauert
  wake-skew naturally from thrust forcing (exact for axial and
  pure-longitudinal flight; approximate for oblique `µ_y ≠ 0`)

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

## Pitt-Peters design notes (`dynbem/pitt_peters.py`)

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

## Øye design notes (`dynbem/oye.py`)

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
