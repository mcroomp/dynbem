# dynbem -- a rotor aero model for the RAWES work

An overview of a blade-element rotor code built to support
Christof's airborne wind energy project. I focus on the **quasi-static
BEM**, the part I understand best and use most.

> Before I start: I'm a software/engineering developer, not an
> aerodynamicist. I built this with a lot of AI assistance and by
> leaning on published rotor papers for the physics and the validation
> numbers. So read this as "here is what the tool does and where
> I've checked it," not as me claiming domain expertise. Corrections
> welcome.

---

## Slide 1 -- What this is, and why I built it

- The goal was a practical rotor aero model to support **Christof's
  RAWES work** -- a tethered rotary kite steered with cyclic pitch.
- The system-level AWE models I read (e.g. De Schutter 2018) represent
  the rotor with a **single actuator-disk induction factor** `a`. That
  is elegant for optimal control, but I wanted something that resolves
  the blade a bit more, so we could look at loads and at off-axis
  (tilted) operation.
- `dynbem` is the result: a small **blade-element** rotor code, written
  in Rust with a Python API, quick enough to evaluate many times inside
  a trim or sweep.
- I'm not presenting a new method -- it's a careful implementation of
  textbook/paper methods. The value I'm hoping to offer is a working,
  checked tool, and an honest account of its limits.

---

## Slide 2 -- The quasi-static BEM

This is the model I understand best, so it's the one I'll concentrate on.

- It's **blade-element momentum theory**: slice the blade along its span,
  work out the local airflow and forces on each slice, and step around
  the rotor azimuth so a tilted rotor's advancing/retreating difference
  shows up.
- **"Quasi-static"** means there's no memory between calls -- each annulus
  solves its own steady inflow inside the calculation. So one call is just
  a function of the operating point in, forces out. (As I understand it,
  that's the same spirit as the single `a` in the actuator-disk model,
  just solved locally per slice instead of once for the whole disk.)
- It handles **both** a rotor pushing air (helicopter-like) and a rotor
  freewheeling on the wind (autorotation / windmill-like), in one code
  path. For the windmill side I followed Ning's 2014 method. In this
  application the two modes map onto distinct phases of flight rather
  than a single "turbine" duty (see next slide).
- It also includes the standard **tip/hub loss** corrections, real
  **airfoil data tables**, and cyclic pitch input.
- There's also a **quasi-static blade-flapping correction** (an
  equivalent spring-hinge that bleeds off most of the hub moment) and
  support for **Kaman-style servo-flaps** (a 1/rev feathering solve) --
  another feature I haven't seen in mainstream rotor BEM codes. Both are
  quasi-static approximations
  with known gaps -- e.g. the servo-flap path doesn't yet add a direct
  sectional lift increment or an actuator-lag model -- so I'd flag them
  as "implemented and plausible" rather than "validated."

---

## Slide 3 -- How the flight phases map onto the model

You all know the pumping cycle, so just the part that touches the aero
model:

- **Helicopter (thrusting) mode is only takeoff and landing** -- an
  external motor spins the rotor up and releases it, and the same case
  handles the final flare. The rotor only actively pushes air there.
- **The working flight is autorotation** (windmill mode); the harvest is
  the winch reel-out/reel-in delta, not a rotor torque path. Having both
  modes in one code path means no model switch between the phases.
- **The key reason I think BEM matters here is non-axial inflow.** A
  steered kite flies with its rotor axis well off the wind, so the disk
  sees strong in-plane flow and a large advancing/retreating asymmetry.
  My understanding is that solvers like **OpenFAST (AeroDyn BEMT) only
  run the windmill/energy-extraction branch in axial flow** -- the
  azimuth-resolved windmill solve here was specifically extended to work
  with in-plane wind inside the `psi`-loop, which is the regime this
  application sits in. I'd really like to know if I've mischaracterised
  that.
- Because it carries **no state**, the aero call can be used directly in
  a steady trim solve, and the **cyclic pitch** used to steer maps onto
  the model's existing pitch inputs.

---

## Slide 4 -- How it relates to the actuator-disk model you may know

Same basic idea (steady momentum balance, no time history), at a finer
resolution. Corrections to this comparison welcome.

| Aspect              | Actuator disk (AWE system model) | This quasi-static BEM           |
|---------------------|----------------------------------|---------------------------------|
| Induction           | Single number `a` for the rotor  | Solved per spanwise slice       |
| Along the blade     | Not resolved                     | Full span + tip/hub loss        |
| Around the rotor    | Not resolved (axial)             | Stepped through azimuth         |
| Energy extraction   | Momentum only, axial             | Windmill solver + in-plane wind |
| Airfoil             | Linear lift, simple drag         | Measured polar tables, stall    |
| Hub moments         | Not modeled                      | Computed from the blade forces  |

My honest summary: this is a **more detailed version of the same
closure**, not a different theory. It adds resolution beyond a single
`a`, while staying compatible with how the AWE models are set up.
The one place it may go beyond a mainstream BEM code like OpenFAST is the
**windmill solve with in-plane wind** (the non-axial row above) -- again,
a claim I'd like checked.

---

## Slide 5 -- Adding unsteady inflow: Pitt-Peters

When the rotor state changes quickly (gusts, fast cyclic, manoeuvres),
the steady inflow assumption breaks down. Pitt-Peters is the first of two
optional **dynamic-inflow** models that relax it.

- It carries **three inflow states** -- one uniform plus two cyclic
  (a fore-aft and a side-to-side tilt of the inflow across the disk) --
  and lets them **lag** the loads through first-order ODEs rather than
  snapping to equilibrium.
- A small coupling matrix ties thrust and the two hub moments into those
  inflow states, so it reproduces **wake skew** (the inflow leaning back
  as the rotor tilts) and the way cyclic inflow **eases the hub moment**
  a steered rotor would otherwise feel.
- I followed Peters' own formulation for the matrix and time constants;
  the signs were the fiddly part and are documented in the repo.
- Honest caveat: it can get **numerically stiff** at high advance ratio
  combined with descent -- exactly some of the tilted, autorotating
  points this application cares about -- which is part of why a second
  model exists.

---

## Slide 6 -- A steadier alternative: Oye

Oye is the second dynamic-inflow option, chosen specifically for the
regimes where Pitt-Peters becomes stiff.

- Instead of disk-wide harmonics, it puts a **two-stage time filter on
  each annulus** independently: the momentum target feeds an intermediate
  state, which feeds the inflow the blade actually sees.
- Because the annuli don't talk to each other, there's **no global
  feedback loop**, so it stays numerically well-behaved in descent and
  high-advance-ratio cases that strain Pitt-Peters. This is the same
  trade-off OpenFAST's DBEMT makes (I used its default filter constant).
- The cost: it has **no cyclic-inflow harmonics and no wake-skew term**,
  so it won't reproduce the inflow-driven hub-moment relief that
  Pitt-Peters does. Cyclic *control* still works; the cyclic inflow
  *feedback* doesn't.
- So the rough rule I'm working with: **Pitt-Peters for richer physics,
  Oye for stability** -- and quasi-static when I just need steady loads.
  I'd welcome views on which is the right default here.

---

## Slide 7 -- Where I've checked it, and what I'm unsure about

- **Checked against published rotor data**: Caradonna-Tung, NREL Phase VI,
  Wheatley autorotation, Castles-Gray, and cross-compared with the open
  codes CCBlade / XROTOR / OpenFAST. Details in `verification/` and
  `EMPIRICAL_VALIDATION.md`. I weight these checks heavily because the
  reference numbers are not mine.
- **Using it**: `create_aero(rotor, model="quasi_static")`, then
  `compute_forces(inputs, state)` returns forces, torque, and hub moments.
- **If unsteady effects matter**, the same geometry can run under the two
  dynamic-inflow models above (Pitt-Peters, Oye) -- I leaned on the
  literature for those and would treat them more cautiously.
- **What I'd most like your input on**: whether the quasi-static
  assumptions hold for a **steadily autorotating, tilted** rotor under a
  reel-out/reel-in cycle, where the induction model is weakest, and which
  validation cases I should add for this application.

---

## The library

The Python package is published and open source (MIT):

- **Install**: `pip install dynbem`
- **PyPI**: https://pypi.org/project/dynbem/
- **Source / issues**: https://github.com/mcroomp/dynbem

**Reference**: `dynbem` -- Rust-backed BEM / Pitt-Peters / Oye
dynamic-inflow rotor models, v0.4.0. MIT License.
https://github.com/mcroomp/dynbem

---

*Built with AI assistance, on top of published methods, to support
Christof's project. Code in `dynbem_rs/` (Rust) and `dynbem/` (Python);
background in README.md and EMPIRICAL_VALIDATION.md.*
