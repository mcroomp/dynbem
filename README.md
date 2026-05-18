# aero

Rotor aerodynamics for a flying-wind-turbine simulator. Multi-element BEM
with hover-safe inflow iteration, Pitt-Peters dynamic inflow with VRS
correction, and a flight-envelope sweep driver.

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
import aero

defn = aero.RotorDefinition(...)        # geometry, inertia, controls
model = aero.PittPetersModel(defn)      # or aero.BEMModel(defn)
state = model.initial_rotor_state()
result, dstate = model.compute_forces(inputs, state)
```

See `aero/__init__.py` for the full public surface (`AeroBase`,
`BEMModel`, `PittPetersModel`, `RotorInputs`, `RotorDefinition`, polars,
state types).

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

The `tests/` directory contains both unit tests (`test_*.py`) and
validation scripts against published rotor data (`val_step*.py`,
`validate_table_i.py`).

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
- Per-element Prandtl tip loss
- Glauert / Buhl Windmill Brake State correction (quadratic root
  selection)
- `dω/dt = (Q_aero + Q_motor) / I_ode` — rotor speed integrated as ODE
  state
- `dψ/dt = ω` — spin angle integrated as ODE state
- Returns `QuasiStaticRotorState` derivative
- **Validation target**: Caradonna-Tung rotor (NASA TM-81232, 1981) —
  2-blade NACA 0012, CT vs collective at 5°/8°/12°
  (see `Research/CaradonnaTung/` for CT tables and test notes)

### Level 2 — Pitt-Peters 3-state dynamic inflow ✅ DONE

- `PittPetersModel` in `aero/pitt_peters.py`
- Prescribed-inflow blade element loop: blade sees
  `λ_total = λ_0 + v_climb/ΩR` (induced state + freestream), so WBS
  and autorotation work correctly
- Pitt-Peters ODE: `dλ_0/dt = (λ_0_ss − λ_0) / τ_0`  with
  `τ_0 = 8R/(3π V_T)`
- Cyclic states `λ_c, λ_s` decay to zero for axial flight
  (forward-flight azimuth integration deferred to Level 3)
- **VRS empirical correction** (Leishman 2000, fit to Castles-Gray
  data): in 0 < λ₂ < 2,
  `λ_0_ss` comes from the polynomial
  `λ₁/V_h = 1 + 1.125λ₂ − 1.372λ₂² + 1.718λ₂³ − 0.655λ₂⁴`
  rather than momentum theory, preventing the Level-1 CT blow-up in VRS.
- Apparent-mass time constants from Peters & HaQuang (1988):
  `τ_0 = 8R/(3π V_T)`, `τ_cs = 16R/(45π V_T)`
- **Validation**: `tests/test_pitt_peters.py` — hover CT vs Level-1
  BEM, VRS no-blow-up, WBS autorotation sign, first-order inflow lag
- **Known limitations**:
  - VRS CT still rises to ~2× nominal in deep VRS (λ₂ ≈ 1.5–2) at
    fixed θ; real rotor stays near nominal (paper: θ barely adjusts).
    The Leishman polynomial shifts the operating point but doesn't
    fully suppress it.
  - Autorotation torque crossing at V/ΩR ≈ 0.14 vs paper's 0.083.
  - Forward flight (µ ≠ 0): λ_c/λ_s targets are zero; azimuth
    integration needed.

### Level 3 — Peters-He finite-state dynamic inflow (state of the art)

- 9-state (or higher-order) Peters-He inflow model
- Requires new `PetersHeRotorState` dataclass
- Captures higher harmonics of the inflow distribution
- Best accuracy for maneuvering flight and aeroelastic coupling
- **Validation target**: Caradonna-Tung unsteady / forward-flight data

### Forward flight (applies to all levels)

- Oblique inflow: advance ratio `µ = V_edge / (Ω·R)` ≠ 0
- Blade azimuth-dependent velocity in the BEM loop
- Required for the flying turbine in translating flight

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

## Pitt-Peters design notes (`aero/pitt_peters.py`)

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
