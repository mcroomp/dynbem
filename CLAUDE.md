# CLAUDE.md — AI assistant instructions

Human-facing docs (install, usage, coordinate conventions, design notes,
implementation roadmap, research sources) live in [README.md](README.md).
**Read it first** — most of what you'd want to know about this codebase is
there, and this file does not repeat it.

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

### v_t_extra and Glauert cyclic targets

The tangential apparent wind at a blade element is, from first
principles:

    v_t = (v_blade − v_air) · t_hat = ω·r − v_inplane · t_hat

so `v_t_extra = −v_inplane · t_hat`. With the CCW t_hat above, in hub
frame components:

    v_t_extra = +v_in_hub_x · sin(ψ) + v_in_hub_y · cos(ψ)

Glauert cyclic steady-state targets under this convention (max inflow
at the back of disk, i.e. at ψ such that `r_hat(ψ) = −v_hub/|v_hub|`):

    λ_c_ss = +mu_x · tan(χ/2)
    λ_s_ss = −mu_y · tan(χ/2)

where `mu_x, mu_y = v_inplane_hub / Ω_R` (note: these are
*wind-relative* advance ratios, **opposite sign** to the vehicle's
forward-speed advance ratio used in most textbooks — that's why the
λ_c_ss sign here looks inverted vs Leishman/Johnson). The λ_c / λ_s
asymmetry (one `+`, one `−`) is real: it follows from
`r_hat_y = −sin(ψ)` having a negation that `r_hat_x = +cos(ψ)` does
not.

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
`v_t_extra` and Glauert signs. Do not flip the sign of ω.

## Workflow

- **Python**: use the venv at `.venv\`. Activate with
  `.venv\Scripts\activate` (Windows), or invoke directly via
  `.venv\Scripts\python` / `.venv\Scripts\pytest`. Don't install packages
  globally or create a new venv.
- **Coordinate frame**: NED everywhere. See README "Coordinate system" —
  the "coordinate trap" section especially matters when you adapt
  equations from a paper, because most rotor literature uses a different
  frame and the sign flips are easy to miss.
- **Sign conventions**: before changing any inflow / thrust / torque
  sign, re-read the README "BEM solver design" and "Pitt-Peters design
  notes" sections. The signs are load-bearing and were tuned to make
  hover, climb, descent, VRS, and autorotation all work in one code
  path.

## When extending the aero models

- New levels (e.g. Peters-He) plug in behind the `AeroBase` interface
  in `aero/__init__.py`. Don't break existing call sites — keep
  `compute_forces(inputs, state) -> (AeroResult, RotorState)`.
- Validation data lives under `Research/`. When adding a new model
  level, add a `tests/test_<model>.py` and, if appropriate, a
  `val_step*.py` script that compares against a specific paper's data.
- Don't store derived results inside `Research/` — that directory is
  for source-paper extractions only.

## Subfolder CLAUDE.md files

- `Research/CLAUDE.md` — extraction conventions for paper sources.
- `Research/CaradonnaTung/CLAUDE.md` — Caradonna-Tung page index, CT
  tables, validation notes.

Defer to those when working inside the respective directories.
