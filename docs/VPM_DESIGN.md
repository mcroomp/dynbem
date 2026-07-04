# Vortex Particle Method (VPM) -- design doc

Status: **experimental / in progress**. The wake engine, the forward-flight
rotor coupling, and the blade flap DOF are implemented. The model is
exposed as a first-class `"vpm"` option in `create_aero()` and is validated
against measured Wheatley TR-515 forward-flight autorotation data and
classical rotor theory (Section 8).

This document describes what is built, the math it implements, the
conventions it follows, its measured performance, and the roadmap. It is
the design reference for anyone extending the VPM path. For the BEM /
Pitt-Peters / Oye models and the shared sign/frame conventions, see
[../AGENTS.md](../AGENTS.md) and [CLAUDE.md](../dynbem_rs/CLAUDE.md).

Math is written in GitHub Flavored Markdown LaTeX (`$...$` inline,
`$$...$$` and ` ```math ` blocks), per the design-document carve-out in
[../AGENTS.md](../AGENTS.md). The surrounding prose stays ASCII (no
em-dashes or smart quotes).

---

## 1. Motivation and scope

The BEM family (quasi-static BEM, Pitt-Peters, Oye) models the rotor inflow
with momentum theory plus empirical regime patches: the Leishman VRS
polynomial, the Ning windmill Brent solver, and the Pitt-Peters wake-skew
L-matrix. These are cheap (~0.1 ms/call) but the wake geometry is never
represented -- it is baked into the inflow model.

The VPM represents the wake explicitly as a cloud of regularized vortex
particles and computes the induced velocity field directly. The motivation
is that several effects the BEM models handle with separate empirical
patches are, in the VPM formulation, consequences of the same wake
convection rather than special cases:

- **Vortex ring state** would arise from wake recirculation rather than
  from an empirical polynomial.
- **Descent / windmill / autorotation** need no separate solver branch; the
  wake convects in whatever direction the induced velocities dictate.
- **Forward-flight wake skew** is geometric (the wake convects downstream)
  rather than an L-matrix off-diagonal.

These are properties of the method, not claims about the present
implementation -- Section 8 is explicit that the descent / VRS regime is
not validated here. The cost is several orders of magnitude more compute
per operating point (Section 7), so the VPM path is aimed at wake-fidelity
questions (VRS, blade-vortex interaction, rotor-rotor / rotor-body
interference), not at replacing the inflow models in trim sweeps.

What VPM still needs from a sectional model: it is inviscid and only
produces the induced velocity field. Sectional $C_l$ / $C_d$ still come from
the same airfoil polar tables as BEM, so stall, reversed flow, and
compressibility remain airfoil-table concerns.

---

## 2. Mathematical formulation

Standard VPM following Winckelmans & Leonard (1993) and the modern
rotorcraft formulation (Alvarez & Ning 2020 / FLOWVPM).

### 2.1 Vorticity discretization

The wake vorticity field is a sum of regularized particles:

$$
\boldsymbol{\omega}(\mathbf{x}) = \sum_p \boldsymbol{\alpha}_p \,
\zeta_{\sigma_p}(\mathbf{x} - \mathbf{x}_p)
$$

- $\boldsymbol{\alpha}_p$ -- particle vector strength, the integral of
  vorticity it carries. Units $\mathrm{m^3/s}$ (circulation times length).
- $\zeta_{\sigma_p}$ -- radially symmetric regularization kernel, core
  size $\sigma_p$.

### 2.2 Regularized Biot-Savart law

$$
\mathbf{u}(\mathbf{x}) = \frac{1}{4\pi} \sum_p g(\rho_p)\,
\frac{\boldsymbol{\alpha}_p \times (\mathbf{x} - \mathbf{x}_p)}{r_p^{\,3}}
$$

with $r_p = |\mathbf{x} - \mathbf{x}_p|$, $\rho_p = r_p / \sigma_p$, and the
low-order algebraic smoothing function

$$
g(\rho) = \frac{\rho^3 \left(\rho^2 + \tfrac{5}{2}\right)}{(\rho^2 + 1)^{5/2}}
$$

(paired with $\zeta(\rho) = \tfrac{15}{8\pi}(\rho^2 + 1)^{-7/2}$). As
$\rho \to \infty$, $g \to 1$ and the singular Biot-Savart law is recovered;
as $\rho \to 0$, $g \sim \tfrac{5}{2}\rho^3$ so the velocity stays finite.

### 2.3 Singularity-free rearrangement (what is evaluated)

The implementation never divides by $r$. Folding $\rho^3 = r^3/\sigma^3$ into
$g/r^3$ gives a kernel that depends only on $\rho^2$ and $\sigma$:

$$
K(\rho, \sigma) = \frac{g(\rho)}{r^3}
= \frac{\rho^2 + \tfrac{5}{2}}{\sigma^3 (\rho^2 + 1)^{5/2}}
$$

$$
\mathbf{u}(\mathbf{x}) = \frac{1}{4\pi} \sum_p K(\rho_p, \sigma_p)\,
\big[\boldsymbol{\alpha}_p \times (\mathbf{x} - \mathbf{x}_p)\big]
$$

Consequences that make this vectorization-friendly:

- $K$ is finite for all $r$, with $K(0, \sigma) = \tfrac{5}{2}/\sigma^3$.
- The self term ($\mathbf{x} = \mathbf{x}_p$) contributes exactly zero
  because $\boldsymbol{\alpha}_p \times \mathbf{0} = \mathbf{0}$. No masks,
  no branches, no division by $r$.
- $(\rho^2 + 1)^{5/2}$ is evaluated as $b^2\sqrt{b}$ with $b = \rho^2 + 1$
  -- one square root, the rest multiplies.

The source particle's own $\sigma$ is used (source-core convention).

### 2.4 Sign convention

A particle strength along $+z$, $\boldsymbol{\alpha} = (0, 0, A)$ with
$A > 0$, induces velocity in $+y$ at a point on the $+x$ axis --
counter-clockwise circulation about $+z$, the right-hand rule for vorticity
along $+z$. This is the sign relation the far-field unit test enforces
(Section 8.1).

The rotor coupling uses the project-wide NED frame and CCW-from-above
convention (see AGENTS.md): hub axis $+Z$, with

$$
\hat{\mathbf{r}}(\psi) = [\cos\psi,\ -\sin\psi,\ 0], \qquad
\hat{\mathbf{t}}(\psi) = [-\sin\psi,\ -\cos\psi,\ 0].
$$

---

## 3. Data structures

The particle field is stored **structure-of-arrays** (SoA), single
precision. Logically it is seven parallel arrays of length `N`, indexed by
particle `p`:

    position   px[p], py[p], pz[p]
    strength   ax[p], ay[p], az[p]     # vector strength alpha_p
    core       sigma[p]

Every per-particle quantity is a pure function of index `p` alone, so all
seven arrays can be read, written, and reduced independently -- the layout
is chosen precisely so that any loop over `p` is a data-parallel map (or
map-reduce) with no cross-particle dependence in the inner kernel.

The field exposes a `particles()` iterator that yields `(pos, strength, sigma)`
tuples for each particle. This is the read path used by the `vpm_viz`
visualiser and any downstream analysis.

Rationale:

- **SoA, not array-of-structs**: a data-parallel kernel over `p` touches
  one component array contiguously, so a width-`W` vector lane consumes `W`
  adjacent particles per step. AoS interleaves components and would force
  gathers, defeating vectorization.
- **f32**: the evaluation is compute-bound and the wake spans meters, so
  single precision is adequate here and roughly doubles vector throughput.
  Many VPM codes run in f64; this choice is local to the engine and
  revisitable. The rest of the codebase is f64.

A scalar f64 reduction of the same sum is retained as the readable
specification and the reference the vectorized path is checked against --
not the hot path.

---

## 4. Numerical methods (implemented)

### 4.1 Velocity evaluation -- direct O(N^2), data-parallel

Given `M` target points and `N` source particles, evaluate the regularized
Biot-Savart sum (Section 2.3) at every target. Self-evaluation (wake on
itself) is the special case where the targets are the source positions.

As a logical operation this is a **map over targets of a reduction over
sources**. The induced velocity at target $t$ is

$$
\mathbf{u}_t = \frac{1}{4\pi} \sum_s K(\rho_{ts}, \sigma_s)\,
\big[\boldsymbol{\alpha}_s \times (\mathbf{x}_t - \mathbf{x}_s)\big],
\qquad \rho_{ts} = \frac{|\mathbf{x}_t - \mathbf{x}_s|}{\sigma_s}
$$

- The outer map over $t$ is independent: targets never interact, so it maps
  onto both levels of CPU data-parallelism -- vector lanes within a core
  (x86 AVX, ARM64 NEON/SVE) and separate cores across the loop -- with no
  synchronization.
- The inner reduction over $s$ is a sum, hence associative -- it may be
  reassociated freely into lane-wise partial sums that are combined at the
  end (this is what lets $W$ sources be processed per vector step).
- The kernel $K$ is branch-free and finite for all separations (Section
  2.3), and the self term contributes exactly zero, so the reduction needs
  no masking or special-casing of $t = s$.
- The source count is padded up to a multiple of the vector width with
  inert particles (zero strength, unit core) so the tail of the reduction
  needs no scalar remainder path; inert particles add exactly zero.

The only floating-point subtlety is that the lane-wise reassociation gives
results that differ from a strict left-to-right scalar sum at the rounding
level; the scalar f64 reduction (Section 3) is the reference this is
checked against.

### 4.2 Time integration -- RK2 (midpoint)

One convection step of the free wake, freestream $\mathbf{U}_\infty$, step
$dt$, with $\mathbf{u}_\text{ind}(\cdot)$ the evaluation of Section 4.1:

```math
\begin{aligned}
\mathbf{u}_1 &= \mathbf{u}_\text{ind}(\mathbf{x}) \\
\mathbf{x}_\text{mid} &= \mathbf{x} + \tfrac{1}{2}\,dt\,(\mathbf{u}_1 + \mathbf{U}_\infty) \\
\mathbf{u}_2 &= \mathbf{u}_\text{ind}(\mathbf{x}_\text{mid}) \\
\mathbf{x} &\leftarrow \mathbf{x} + dt\,(\mathbf{u}_2 + \mathbf{U}_\infty)
\end{aligned}
```

Both velocity evaluations are the data-parallel map of Section 4.1; the
position updates are an elementwise (per-particle) parallel map. Strengths
are held constant (convection only -- stretching and diffusion are not
modelled yet, Section 6).

---

## 5. Rotor coupling I -- axial free-wake precursor (example)

The axial precursor is a minimal, standard lifting-line free-wake coupling
used to sanity check VPM thrust/torque against BEM. It is the pedagogical
precursor to the module in Section 5.5; it is retained because its
trailing-only, axisymmetric structure isolates the trailed-wake
contribution.

### 5.1 Scope and simplifications (all standard for a first coupling)

- **Axial flow only** -- no cyclic, azimuthally symmetric loading.
- **Trailing vorticity only** -- at steady state $d\Gamma/dt \to 0$, so the
  shed (temporal) vorticity vanishes and the trailed (radial-gradient)
  vorticity is the whole wake. This is the classic free-tip-vortex wake
  (Landgrebe / Bagai-Leishman) discretized as particles.
- **Lifting-line loads** via Kutta-Joukowski: $\Gamma = \tfrac{1}{2} U c\,C_l$.
- **Near-wake self-induction only** -- no bound-vortex self term (standard
  for a straight lifting line).

Shared infrastructure: the radial grid (r_mid, chord, twist per station),
the airfoil polar table (Cl/Cd interpolation), and the VPM velocity
evaluation and RK2 advection routines (Sections 4.1-4.2).

### 5.2 Per-step algorithm

```mermaid
flowchart TD
    A[probe wake-induced velocity at blade-0 stations] --> B[blade-element loads]
    B --> C[Gamma_i = 0.5 U c Cl, under-relax]
    C --> D[shed trailing particles from all blades]
    D --> E[RK2 advect free wake]
    E --> F[FIFO-truncate oldest wake beyond max_particles cap]
    F --> G[accumulate thrust/torque over final rev]
    G --> A
```

Blade-element step (station $i$, azimuth $\psi_0$), with
$U = \sqrt{U_a^2 + U_t^2}$:

```math
\begin{aligned}
\mathbf{v}_\text{air} &= \mathbf{U}_\infty + \mathbf{u}_\text{ind},
   \quad \mathbf{U}_\infty = [0,\,0,\,v_\text{climb}] \\
\mathbf{v}_\text{blade} &= \omega r\,\hat{\mathbf{t}}(\psi_0) \\
\mathbf{U}_\text{rel} &= \mathbf{v}_\text{air} - \mathbf{v}_\text{blade} \\
U_a &= \mathbf{U}_\text{rel}\cdot\hat{\mathbf{z}}
   \quad (\text{axial, } +Z) \\
U_t &= -\,\mathbf{U}_\text{rel}\cdot\hat{\mathbf{t}}
   \quad (\text{tangential, along } -\hat{\mathbf{t}}) \\
\phi &= \mathrm{atan2}(U_a,\, U_t) \\
\alpha &= (\theta_\text{coll} + \theta_{\mathrm{tw},i}) - \phi \\
(C_l, C_d) &= \text{polar}(\alpha) \\
dL &= \tfrac{1}{2}\rho U^2 c\,C_l\,dr, \qquad
dD = \tfrac{1}{2}\rho U^2 c\,C_d\,dr \\
dT &= dL\cos\phi - dD\sin\phi \\
dQ &= (dL\sin\phi + dD\cos\phi)\,r \\
\Gamma_i &= \tfrac{1}{2} U c\,C_l
\end{aligned}
```

Trailing shed (edge $j = 0 \ldots n$, both blades; $\Gamma = 0$ outside the
blade):

```math
\begin{aligned}
\Gamma_{\mathrm{trail},j} &= \Gamma_{j-1} - \Gamma_j \\
\mathbf{s}_j &= \big(-\omega\, r_\text{edge}\,\hat{\mathbf{t}}
   + v_\text{climb}\,\hat{\mathbf{z}}\big)\,dt \\
\boldsymbol{\alpha}_j &= \Gamma_{\mathrm{trail},j}\,\mathbf{s}_j
   \quad \text{at } r_\text{edge}\,\hat{\mathbf{r}}(\psi)
\end{aligned}
```

The tip edge carries `Gamma_tip` (the tip vortex, strongest); the
root edge carries `-Gamma_root`.

### 5.3 Convergence and truncation

- Bound circulation is under-relaxed:
  $\Gamma \leftarrow \Gamma + \text{relax}\,(\Gamma_\text{new} - \Gamma)$.
- Wake is FIFO-truncated to
  $N_\text{wake\,rev} \cdot N_\text{steps/rev} \cdot N_\text{shed/step}$
  particles (drop oldest).
- The rotor is marched 8 revolutions; thrust/torque are averaged over the
  final revolution.

### 5.4 Default parameters

| Parameter | Value |
|---|---|
| stations | 20 |
| steps per rev | 24 (15 deg) |
| wake age retained | 4 revs |
| core size sigma | 0.18 m |
| relaxation | 0.35 |
| total revs | 8 |
| implicit lifting-line | on |
| tip clustering (cosine) | on |
| local core sizing | on |

Steady-state particle count $N \approx 2 \cdot 21 \cdot 24 \cdot 4 = 4032$.

---

## 5.5 Rotor coupling II -- unsteady forward-flight free-wake (library)

The forward-flight coupling is an unsteady lifting-line free-wake coupling
(in the spirit of Leishman 2006, ch. 10; Bagai & Leishman 1995) that
supports cyclic pitch and arbitrary in-plane (edgewise / crosswind) inflow.
It removes the two structural restrictions of Section 5 -- azimuthal
symmetry and steady bound circulation -- while reusing the same radial-grid
/ polar-table geometry and the VPM engine.

### 5.5.1 State and time-march

The rotor is marched at fixed azimuthal step $d\psi = 2\pi / N_\psi$,
$dt = d\psi / \omega$. Each blade $b$ sits at
$\psi_b(n) = n\,d\psi + b\,\tfrac{2\pi}{N_b}$. The carried state (`VpmRotorState`) has
three fields:

- **`wake`** -- the shared free-wake particle field (`Option<ParticleField>`).
- **`gamma`** -- the previous-step relaxed bound circulation
  $\Gamma_i^{\,n-1}$ per blade, needed for the shed term.
- **`psi`** -- the current blade-0 azimuth (rad), accumulated across
  `step_one` calls so each call sheds particles from the correct rotor
  position rather than always restarting at $\psi = 0$.
  `march` always starts at $\psi = 0$ and ignores this field.

The solution sought is the **periodic** limit cycle, not a fixed point: loads
are phase-averaged over the final revolution once the wake has settled
($N_\text{settle\,rev}$ revolutions of spin-up).

### 5.5.2 Sectional aerodynamics -- implicit lifting-line (per blade, per step)

The bound circulation is the solution of a **nonlinear (Prandtl) lifting-line
matching problem**, solved implicitly each step. The induced velocity a
section sees is split into two parts:

- a **far field** $\mathbf{u}_\text{far}$ -- the shed particle wake from all
  previous steps, evaluated once per step (Section 4.1) and **frozen**
  during the solve;
- a **near field** $\mathbf{u}_\text{near}(\Gamma)$ -- the trailing vortex
  sheet just behind the blade, linear in the current-step circulation
  $\Gamma$ (Section 5.5.3).

With $U = \sqrt{U_a^2 + U_t^2}$ and

$$\mathbf{u}_\text{ind} = \mathbf{u}_\text{far} + \mathbf{u}_\text{near}(\Gamma):$$

```math
\begin{aligned}
\mathbf{v}_\text{air} &= \mathbf{U}_\infty + \mathbf{u}_\text{far} + \mathbf{u}_\text{near}(\Gamma)
   \quad (\mathbf{U}_\infty = \text{3D hub-frame freestream}) \\
\mathbf{v}_\text{blade} &= \omega r\,\hat{\mathbf{t}}(\psi_b) \\
\mathbf{U}_\text{rel} &= \mathbf{v}_\text{air} - \mathbf{v}_\text{blade} \\
U_a &= \mathbf{U}_\text{rel}\cdot\hat{\mathbf{z}}
   \quad (\text{through-disk, } +Z) \\
U_t &= -\,\mathbf{U}_\text{rel}\cdot\hat{\mathbf{t}}
   \quad (\text{in-plane, along } -\hat{\mathbf{t}}) \\
\phi &= \mathrm{atan2}(U_a,\, U_t)
   \quad (\text{inflow angle}) \\
\theta &= \theta_\text{coll} + \theta_{\mathrm{tw},i}
   + \theta_{1c}\cos\psi_b + \theta_{1s}\sin\psi_b
   \quad (\text{cyclic, Section 5.5.4}) \\
\alpha &= \theta - \phi \\
(C_l, C_d) &= \text{polar}(\alpha) \\
\Gamma_i &= \tfrac{1}{2} U c\,C_l
   \quad (\text{Kutta-Joukowski})
\end{aligned}
```

Because `u_near` depends on $\Gamma$ and $\Gamma$ depends on
$\alpha$ (hence on `u_near`), the last two lines are a coupled
fixed point. It is solved by relaxed iteration, seeded from the previous
step's circulation:

$$\Gamma_i \leftarrow \Gamma_i + \text{relax}\,\left(\tfrac12 U c\,C_l - \Gamma_i\right),$$

repeated until $\Gamma$ stops changing (a few tens of iterations at most).
The far field is not re-probed inside the loop, so each iteration is only
the small near-field matrix-vector product plus a polar lookup -- negligible
against the $O(N^2)$ particle evaluation.

Lift comes from the **polar**, so (unlike a thin-airfoil vortex-lattice) the
bound vortex carries no self term; only the *trailing* wake enters
`u_near`. This is the classical viscous lifting-line, in the
Bagai-Leishman free-wake / FLOWVPM tradition. Setting
`nonlinear_lifting_line = false` drops `u_near` and reduces
the solve to a single relaxed Kutta-Joukowski pass (the original pointwise
scheme), retained as a baseline.

Why it matters: the tip unloading ("tip loss") is the trailing tip vortex
driving down the local angle of attack. With only the lagged, core-smoothed
particle wake the section never feels its own just-shed trailing vorticity
within the step, so the load does not taper -- the source of the hover
thrust bias in Section 8. The implicit near wake supplies exactly that
within-step downwash, so the roll-off is computed rather than imposed. In
practice this moves the coarse-grid hover bias by a few percent on its own;
cashing in the full taper also needs tip-clustered (cosine) spanwise spacing
plus a core size scaled to the local segment length (`tip_clustering` and
`local_core`). These two compose: a local core is a no-op on a uniform grid
(every segment is the same width), but once the stations cluster at the tip
the shorter tip segments carry proportionally smaller cores that stop
over-smoothing the tip vortex. With both enabled the measured hover bias
falls to under one percent (Section 8).

"Quasi-steady" here means the *section* uses the static polar at the
instantaneous $\alpha$ -- there is no Theodorsen / Leishman-Beddoes
attached-flow lag and no dynamic-stall model. The only unsteadiness carried
is the **wake memory**: shed vorticity plus free convection. At the reduced
frequencies of a 1/rev cyclic input ($k = \omega c / (2 U_\text{tip})$ is
small) the wake response should be the dominant unsteady effect; this is
not adequate for high-$k$ phenomena (blade-vortex interaction pressure
signatures, deep dynamic stall).

### 5.5.3 Wake emission -- trailed and shed vorticity

The near-wake vortex sheet behind a lifting line splits into two
components whose strengths follow from the spanwise and temporal gradients
of bound circulation. Both are emitted each step, one particle per feature:

**Trailed** (streamwise, from $\partial\Gamma/\partial r$). At each radial
edge $j = 0 \ldots n$ ($\Gamma = 0$ outside the blade) the shed-off
circulation is $\Gamma_{j-1} - \Gamma_j$. It is laid down as a filament
aligned with the local relative wind, length $|\mathbf{U}_\text{rel}|\,dt$:

```math
\begin{aligned}
\mathbf{s}_j &= \mathbf{U}_{\text{rel},\,\text{edge}(j)}\,dt
   \quad (\text{edge value = mean of neighbours}) \\
\boldsymbol{\alpha}_j &= (\Gamma_{j-1} - \Gamma_j)\,\mathbf{s}_j
   \quad \text{emitted at } r_\text{edge}(j)\,\hat{\mathbf{r}}(\psi_b)
\end{aligned}
```

The tip edge carries `Gamma_tip` (the tip vortex), the root edge
`-Gamma_root`.

**Shed** (spanwise, from $\partial\Gamma/\partial t$). Kelvin's theorem
applied to a material contour enclosing the blade section and its near wake
requires that any change in bound circulation over $dt$ be balanced by an
equal and opposite spanwise vortex deposited into the wake. Per station:

$$
\boldsymbol{\alpha}_{\mathrm{shed},i} = -\left(\Gamma_i^{\,n} - \Gamma_i^{\,n-1}\right) dr\,\hat{\mathbf{r}}(\psi_b)
$$

i.e. a spanwise segment (parallel to the bound vortex) of length $dr$
carrying circulation $-(\Gamma_i^{\,n} - \Gamma_i^{\,n-1})$, emitted at the
station. This term vanishes identically at the axisymmetric steady state of
Section 5, which is why the axial precursor can omit it; it is non-zero and
load-bearing whenever bound circulation varies with azimuth -- every
cyclic or crosswind case.

All emitted particles use a single core size `sigma` and convect with the
free wake (Section 4.2), so the wake skew in forward flight is geometric --
no prescribed skew angle.

**Near-field influence and the one-step handoff.** The same trailed segments
$\mathbf{s}_j$ are what the implicit solve of Section 5.5.2 uses for the
near field `u_near`: each edge $j$ carries a straight vortex
filament from $r_\text{edge}(j)\,\hat{\mathbf{r}}$ along $\mathbf{s}_j$ with
circulation $\Gamma_{j-1} - \Gamma_j$, and its regularized Biot-Savart
velocity at the station control points gives the influence coefficients
$\mathbf{B}_{ij}$ (a finite-segment kernel with the same core $\sigma$). The
near-field velocity is then $\mathbf{u}_{\text{near},i} = \sum_j \mathbf{B}_{ij}\,(\Gamma_{j-1}-\Gamma_j)$,
linear in $\Gamma$, so the leg geometry is fixed once per step (from the
frozen far-field relative wind) and the solve iterates only the strengths.

The legs are exactly one convection step long, so they represent *this*
step's trailed vorticity and nothing else. At the end of the step the same
vorticity is emitted as the trailed particle above and convects with the
far wake; next step the leg is regenerated fresh. Each step's trailed sheet
is therefore counted once -- as a near-field filament during its own step,
then as a far-field particle thereafter -- with no double counting at the
handoff.

### 5.5.4 Cyclic pitch

The harmonics $(\theta_{1c}, \theta_{1s})$ come from the same
swashplate-tilt-to-blade-pitch harmonic mapping the BEM models use (with the
rotor's control gain / phase). With the repo default gain 1 / phase 0 this
reduces to $\theta_{1c} = -$`tilt_lon`,
$\theta_{1s} = +$`tilt_lat`, giving the helicopter-standard signs
(`tilt_lon > 0` nose-down, `tilt_lat > 0` roll-right)
verified by the moment-sign tests below.

### 5.5.5 Hub load integration

Summed over blades and stations, then phase-averaged over the final
revolution (AGENTS.md hub-frame convention, thrust along $-Z$):

```math
\begin{aligned}
dT &= dL\cos\phi - dD\sin\phi \\
T &= \sum dT \\
Q &= \sum r\,(dL\sin\phi + dD\cos\phi) \\
M_x &= \sum r\,dT\,\sin\psi \quad (\text{roll-right positive}) \\
M_y &= \sum r\,dT\,\cos\psi \quad (\text{pitch-up positive})
\end{aligned}
```

Blade flapping and servo-flap feathering are modeled as optional per-blade
time-domain DOFs (Section 5.6). The rigid-blade path remains the default when
those properties are not supplied.

### 5.5.6 Simplifications and validity envelope

- **Rigid blade by default; optional flap + feather DOFs** -- with no
  `FlapProperties` and direct-mechanical pitch the blade is rigid (pitch
  sets local loading directly, no precession), consistent with the BEM
  convention. Per-blade flap and feathering DOFs can be enabled (Section
  5.6); they are off unless the rotor supplies the matching properties.
- **Quasi-steady polar** -- static `Cl/Cd`, no unsteady-airfoil model
  (Section 5.5.2).
- **One-step near wake** -- the implicit lifting-line trailing sheet extends
  one convection step before handing off to particles (Section 5.5.3); a
  longer near-wake zone would sharpen the tip taper but needs a near/far
  aging split.
- **No reversed-flow / full-360 polar** -- the sectional table is not
  extended into the reversed-flow region, so validity is capped at
  moderate advance ratio (retreating-side reverse-flow patch stays small).
- **Tip-clustered spacing, locally scaled core** -- spanwise stations use
  cosine (tip/root) clustering and each shed/trailed feature carries a core
  scaled to its local segment width (`tip_clustering`, `local_core`, both on
  by default; Section 5.5.2). There is still no core spreading or vortex
  stretching in the free wake -- inherited from the engine (Section 6).

Default resolution mirrors Section 5.4; a coarse preset (12 steps/rev, 2
wake revs, 3 settle revs) is used by the acceptance tests.

## 5.6 Rotor coupling III -- per-blade rigid-blade DOFs (flap + feather)

Because the forward-flight coupling is time-marched per blade, blade
structural DOFs can be integrated in lockstep with the wake instead of being
reduced to a static factor (as the BEM path does for flap). Each blade
carries up to four extra states -- flap ($\beta$, $\dot\beta$) and feather
($\theta_f$, $\dot\theta_f$) -- advanced once per step. Both are off
unless the rotor supplies the matching properties, and the rigid path is
unchanged when they are.

**Flap DOF** (`FlapProperties`, enabled by `config.flap_dynamics`). Out-of-
plane rigid flapping about an equivalent hinge, $\beta > 0$ = flap up (tip
toward $-Z$). Two couplings into the wake, both of which only a free-wake
method can represent:

- *Flap-rate AoA damping*: $\dot\beta$ adds $+r\dot\beta$ to the section
  axial velocity, so flapping up lowers the local angle of attack. The
  aerodynamic flap damping thus emerges from the loads -- there is no
  analytical $\gamma/8$ term.
- *Out-of-plane wake geometry*: the blade sits at $z = -r\beta$, so the
  bound line and every shed/trailed particle is deposited at that height;
  the coned/tilted wake falls out naturally.

The integrated ODE is purely structural/inertial, forced by the aero flap
moment $M_\beta$ (the span sum of $r\,dF_z$, already computed in the loads
loop):

$$I_\beta\,\ddot\beta + I_\beta(\Omega^2 + \omega_{NR}^2)\,\beta = M_\beta,$$

advanced with symplectic Euler. $\omega_{NR}$ is the non-rotating flap
frequency (0 = freely hinged, $\nu_\beta = 1$).

**Feather DOF** (`PitchActuation::ServoFlap`, Kaman path). A passive
feathering rotation driven by a trailing-edge servo-flap. In servo mode the
swashplate collective/cyclic are reinterpreted as flap deflection commands
`delta_f`; the flap's pitching moment `M_servo` drives feathering,
and `theta_f` **replaces** the direct swashplate-to-pitch path in the
section angle of attack. This is the **feathering + damper** architecture:
the blade rides a pitch bearing and feathers freely, restrained by the
mechanical bearing damper `C_theta` (the only dissipation) and the
**aerodynamic spring** `k_aero` from any AC offset:

$$I_\theta\,\ddot\theta_f + C_\theta\,\dot\theta_f + k_\text{aero}\,\theta_f
= M_\text{servo} + M_\text{camber},$$

integrated semi-implicitly (implicit on the damper) for unconditional
stability regardless of damper strength. The aerodynamic spring is

$$k_\text{aero} = \tfrac{1}{2}\rho\,\omega^2 C_{L\alpha}\,\cdot ac\_offset
\cdot \int c\,r^2\,dr,$$

a nose-down restoring torque from the extra lift acting at the AC a distance
`ac_offset` aft of the feathering axis. `M_servo` is accumulated over
the flap span from the true local dynamic pressure, and `M_camber`
(from `blade_Cm_AC`) sets the DC trim. With `ac_offset = 0` (Kaman ideal, axis
at AC) the damper alone sets the cyclic phase lag and DC trim comes from
`M_camber`. There is no artificial control-stiffness spring -- all
constants are physical and measurable. (The alternative torsional-twist
servo-flap architecture is not modelled yet.)

Both DOFs share the per-blade state vector and compose (a flapping,
feathering blade tilts its wake *and* changes its pitch). The flap harmonics
this produces are compared against classical theory in Section 8.4.

---

## 6. Not yet implemented (next steps, priority order)

1. **Exploit the independent target loop across cores** -- the map over
   targets (Section 4.1) has no cross-target dependence, so on top of the
   per-lane vectorization already in place it can be spread across CPU
   cores. The same independence is what the vector ISAs exploit within a
   core (x86 AVX, ARM64 NEON/SVE), so the two levels compose. Expect close
   to linear speedup with core count until the evaluation becomes
   memory-bound.
2. **Vortex stretching** $d\boldsymbol{\alpha}/dt = (\boldsymbol{\alpha}\cdot\nabla)\,\mathbf{u}$ -- needed for 3D
   circulation conservation once the wake distorts. The algebraic kernel has
   an analytic velocity gradient; add it as a second accumulator.
3. **Viscous diffusion** (core spreading or PSE) so long-lived hover/VRS
   wakes do not stay artificially coherent.
4. **Reversed-flow full-360 polar** -- extend the sectional polar to the
   retreating/reversed-flow region for high advance ratio.
5. **Promote the coupling to a first-class inflow model** with the wake as
   carried, serialized state, registered in the model factory alongside the
   BEM / Pitt-Peters / Oye models and exposed through the Python bindings.

Done since first draft: shed (temporal) vorticity and forward-flight cyclic
+ crosswind coupling (Section 5.5); implicit (Prandtl) lifting-line with the
trailing near-wake solved within each step (Section 5.5.2); tip-clustered
(cosine) spanwise spacing and locally scaled core sizing, which together
close the residual hover thrust bias (Section 8); a monopole Barnes-Hut
O(N log N) evaluator (`induced_at_points_bh` / `advect_rk2_bh`) for when N
outgrows the direct path -- see Section 7.1; `ParticleField::particles()`
iterator for downstream wake access; `VpmRotor::step_one(fc, state, dt)`
for per-frame animation stepping with correct `psi` accumulation across
calls (the `psi` field in `VpmRotorState` was added to fix a bug where
all new particles were shed at azimuth 0 when calling `step_one` repeatedly);
`vpm_viz/` standalone egui/eframe visualiser showing the 30-deg crosswind
wake animated in real time (top-view XY + side-view XZ, viridis colormap
by log10 vorticity magnitude). Per-blade rigid-blade DOFs -- flap
($\beta$, $\dot\beta$) and servo-flap feathering ($\theta_f$,
$\dot\theta_f$) -- integrated in the azimuth march (Section 5.6), with the
flap harmonics checked against classical theory (Section 8.4).

---

## 7. Performance (measured, release)

Per-step cost at a non-axial operating point (forward 12 m/s + 5 m/s edgewise
+ descent 2 m/s). The BEM-family models return a converged answer per
`compute_forces` call; one VPM step advances the free wake by a single azimuth
increment (an O(N^2) Biot-Savart probe + RK2 advect over the whole cloud).

**Algorithm comparison** (ms/step). The VPM row is an optimal operating point:
a reasonable wake size (N=5,000) evaluated with the Barnes-Hut tree, parallel.

| Model | ms/step | x Oye |
|---|---:|---:|
| Oye (2-stage annular filter) | 0.105 | 1x |
| Pitt-Peters (3-state L-matrix) | 0.107 | 1x |
| Quasi-static BEM | 10.95 | 104x |
| VPM (Barnes-Hut, parallel, N=5,000) | 54 | 514x |

Oye and Pitt-Peters are algebraic inflow relations (near-free). The BEM is
~100x slower from its per-station Brent root-find. The VPM is not an inflow
model but a wake-resolving tool: its cost is dominated by the Biot-Savart
evaluation and scales with the particle count N.

**VPM cost vs particle count** (ms/step, matched N on the same settled wake,
`rotor_profile vpm-direct,vpm-bh --long`):

| N | direct seq | direct par | BH seq | BH par |
|---:|---:|---:|---:|---:|
| 2,000 | 21 | 11 | 39 | 15 |
| 5,000 | 198 | 59 | 113 | 53 |
| 10,000 | 778 | 227 | 197 | 108 |
| 16,000 | 2,633 | 428 | 708 | 139 |
| 32,000 | 8,750 | 2,265 | 1,058 | 337 |

The direct sum is O(N^2) (cost quadruples per doubling of N); the Barnes-Hut
tree (`VpmRotorConfig::barnes_hut`, off by default) lumps distant particle
clusters into single equivalent vortices and grows O(N log N), so it overtakes
the direct sum as the wake grows -- comparable at N=5k, ~3x faster at N=16k,
~7x at N=32k. Rayon parallelism (default) gives ~2-4x over sequential for the
direct sum and ~2x for the tree (`--seq`/`--par` isolate it). Absolute ms are
machine-dependent; the ratios are the point.

### 7.1 Barnes-Hut tree (O(N log N))

A box of particles seen from far away induces nearly the same velocity as a
single vortex carrying its net circulation -- the regularized kernel's far
field is the singular Biot-Savart law ($K \to 1/r^3$ for $r \gg \sigma$). The
tree (`induced_at_points_bh`, `advect_rk2_bh`) sums nearby vorticity
particle-by-particle and replaces each distant box with one lumped vortex,
evaluated with the *same* eight-wide SIMD kernel. A box of width $s$ at distance
$d$ is lumped when $s/d < \theta$ (`bh_theta`, default 0.5; `theta = 0` recovers
the exact sum). It is flattened into a pre-order node array with escape pointers
(stackless walk, no recursion) and particles are reordered into leaf-contiguous
eight-padded blocks, so traversal is light scalar work and the arithmetic stays
vectorized. Engaged only when `barnes_hut` is set and the wake reaches
`bh_min_particles`.

---

## 8. Validation

The VPM is validated by the `validation_rs` suite: each check builds a VPM
rotor and asserts its loads against either closed-form theory (Section 8.1) or
measured wind-tunnel data (Section 8.2). Run it with:

```
cargo run --release -p validation_rs [filter]
```

which writes `tmp/theory_report.txt` -- one `CHECK` line per data point
(`PASS | FAIL | INFO`). Underneath the coupling, the Biot-Savart engine has
its own `dynbem_rs` unit tests (`cargo test`): kernel far field vs
$A/(4\pi d^2)$, zero self-induction, SIMD vs f64 reference ($<10^{-3}$),
vortex-ring self-propagation (within 30% of the Kelvin speed), and Barnes-Hut
vs direct ($<5\%$ of peak at $\theta = 0.5$).

### 8.1 Theoretical validation (VPM vs closed-form theory)

| Validation (`validation_rs`) | Reference | What it checks | Result |
|---|---|---|---|
| `blade_element_hover` | Combined BEMT, hover (Leishman ch. 3) | Hover thrust coefficient vs closed form | within ~15-25% |
| `climb_momentum` | Axial-climb momentum theory | `C_T` ~ 2 lam_i (lam_i + lam_c); loads fall with climb | PASS (monotone) |
| `glauert_forward_inflow` | Glauert forward-flight inflow | Disk inflow + wake-skew angle | skew <1.2%; inflow <26% |
| `wake_skew` | Wake-skew geometry | Skew grows with mu; covariant under X/Y rotation | PASS |
| `prandtl_tip_loss` | Directional | Tip-loss flag reduces global loads | PASS |
| `autorotation` | Directional | Negative-torque branch reached in descent + edgewise | PASS |
| `flapping_harmonics` | Bramwell / Seddon flap theory | Coning `a0`, longitudinal flap `a1` vs closed form | `a0`, `a1` within ~14% |
| `cyclic_sign` | Directional (AGENTS.md) | Collective monotone; cyclic tilt signs | PASS |
| `flap_directional` | Directional | Flap coning in hover; hub-moment relief under cyclic | PASS |
| `servo_flap` | Directional | Kaman servo-flap feathering (zero / collective / cyclic) | PASS |
| `cyclic_phase_servo` | Directional | Direct-mech pitching `My` vs servo-flap rolling `Mx` | PASS |

Hover thrust tracks BEMT to within ~15-25%; coning `a0` and longitudinal
flapping `a1` match the closed forms to ~14% (the disk tilts back by the right
amount). Lateral flapping `b1` is currently under-predicted -- see the TODO.

### 8.2 Empirical validation (VPM vs measured data)

| Validation (`validation_rs`) | Dataset | What it checks | Result |
|---|---|---|---|
| `hover_castles_gray` | Castles-Gray NACA TN-2474 Table V | Hover thrust vs measured | within ~8-13% |
| `measured_companions` | CG 1600 rpm, CG descent, Wheatley TR-515 | Each theory module anchored to a measured point | PASS |
| `vpm_forward_flight_empirical` | Wheatley & Hood NACA TR-515 Tables III/IV | Forward-flight autorotation lift sweep | <= ~8% to mu = 0.32 |

`validation_rs` also validates the QS / Pitt-Peters / Oye BEM models against
the same Castles-Gray data (`hover_ct_empirical`, `hover_cq_empirical`,
`descent_cq_empirical`) -- those check the BEM family, not the VPM, and are
outside the scope of this doc.

### 8.3 TODO

- **Lateral flapping `b1`** -- VPM under-predicts vs uniform-inflow theory and
  the measured Wheatley phase; needs digitized data.
- **Descent / VRS regime** -- the trailing-only, truncated wake under-resolves
  the recirculating VRS wake; descent loads are not validated.
- **Reversed-flow region** at high mu -- no full-360 polar yet.
- **Trimmed forward-flight loads** -- no trim closure.
- **Absolute hub moments** -- quantitative anchor still missing (trend-only).
- **BVI-sensitive cases** -- not attempted.

---

## 9. Integration points with the codebase

- Frame / sign conventions: NED, CCW-from-above -- see AGENTS.md. The
  coupling's r_hat / t_hat match the shared BEM kinematics.
- Geometry / polars: `VpmRotor<P: Polar>` is generic over the polar type (same
  pattern as `PittPetersModel<P>`, `OyeBEMModel<P>`, `QuasiStaticBEM<P>`).
  The Python-side `VpmRotor()` factory dispatches to `_VpmRotorLinear` or
  `_VpmRotorTabulated` based on the polar instance, matching the BEM factory
  pattern.
- When promoted to a first-class model (Section 6, item 6): the wake becomes
  carried, serialized inflow state, and the model is registered in the model
  factory alongside the others and exposed through the Python bindings.

---

## 10. References

- Winckelmans, G. and Leonard, A. (1993). "Contributions to Vortex Particle
  Methods for the Computation of Three-Dimensional Incompressible Unsteady
  Flows." J. Comput. Phys. 109(2), pp. 247-273.
  DOI: [10.1006/jcph.1993.1216](https://doi.org/10.1006/jcph.1993.1216)
  -- the algebraic regularization kernel.
- Cottet, G.-H. and Koumoutsakos, P. (2000). "Vortex Methods: Theory and
  Practice." Cambridge University Press. ISBN 978-0-521-62 178-3.
  -- general VPM theory, stretching, diffusion.
- Alvarez, E. J. and Ning, A. (2020). "High-Fidelity Modeling of Multirotor
  Aerodynamic Interactions for Aircraft Design." AIAA J. 58(10), pp. 4385-4400.
  DOI: [10.2514/1.J059178](https://doi.org/10.2514/1.J059178)
  -- reformulated VPM / FLOWVPM, modern rotorcraft usage.
- Leishman, J. G. (2006). "Principles of Helicopter Aerodynamics," 2nd ed.
  Cambridge University Press. ISBN 978-0-521-85860-1.
  -- lifting line, Kutta-Joukowski, VRS, tip vortices.
- Bagai, A. and Leishman, J. G. (1995). "Rotor Free-Wake Modeling using a
  Pseudo-Implicit Technique -- Including Comparisons with Experimental Data."
  J. American Helicopter Society 40(3), pp. 29-41.
  DOI: [10.4050/jahs.40.29](https://doi.org/10.4050/jahs.40.29)
  -- rotor free-wake conventions.
- Landgrebe, A. J. (1972). "The Wake Geometry of a Hovering Helicopter Rotor
  and its Influence on Rotor Performance." J. American Helicopter Society
  17(4), pp. 3-15.
  DOI: [10.4050/jahs.17.4.3](https://doi.org/10.4050/jahs.17.4.3)
  -- prescribed / trailing wake geometry.

---

## 11. Gaps for RAWES application

RAWES (`c:\repos\windpower`) is a tethered 4-blade autorotating rotor kite. The rotor disk
is tilted ~67 deg from vertical at the nominal operating point (disk tilt from wind ξ ≈ 29 deg;
tether elevation β ≈ 8 deg); it operates in continuous autorotation with no engine, and the
pumping cycle alternates between a reel-out (climbing, windmill, power-generating) phase and a
reel-in (descending, reduced-power) phase. The production sim uses `dynbem quasi_static` BEM at
400 Hz; the VPM would be an optional higher-fidelity path for specific operating-point analysis,
not a real-time drop-in. What needs to change before VPM produces useful results for RAWES:

### 11.1 Autorotation / windmill torque balance

The RAWES rotor runs in **autorotation**: net shaft torque is zero (or slightly positive to
reel-in drag). The current VPM marches at a fixed `omega_rad_s` from `FlightCondition`; it does
not close the torque balance. For a meaningful autorotating operating point the caller needs to
drive omega to the zero-torque equilibrium, either by running an outer Newton iteration over
`omega_rad_s` until `VpmRotorResult.torque == 0`, or by coupling the VPM to the omega_spin ODE
from `windpower/simulation/dynamics.py` through the `AeroModel::step` interface.

The windmill quadrant (negative `v_climb`, advancing-tip regime, where BEM switches to the
Ning/Brent windmill solver) is present in the quasi-static BEM path but has never been exercised
through the VPM's shedding loop. The VPM's lifting-line solve uses `Kutta-Joukowski + polar`
on every blade element without a separate windmill branch; the polar and circulation solve should
naturally handle windmilling if the inflow angle `phi` goes past 90 deg (rotor acts as turbine),
but this has not been validated. Section 8 explicitly flags the descent/windmill regime as
unvalidated.

### 11.2 Heavily oblique disk tilt (ξ ~ 29-70 deg from wind)

In hover or gentle forward flight the disk is roughly perpendicular to the gravity vector and
the in-plane wind is a small fraction of the tip speed. For RAWES the disk may be tilted up to
~67 deg from vertical (ξ ≈ 29 deg from wind, or steeper during reel-in). The VPM `v_hub` is
already given in hub frame with `v_hub[2]` as axial (through-disk) and `v_hub[0:2]` as
in-plane, so the frame decomposition is correct for any tilt. What is NOT validated:

- **High advance ratio**: RAWES tip speed ~52 m/s at 270 RPM; 10 m/s wind at ξ ≈ 29 deg gives
  in-plane component ~5 m/s → mu ≈ 0.10. This is modest. At reel-in tilt (ξ ≈ 55-70 deg) the
  in-plane component grows to ~8-9 m/s → mu ≈ 0.15-0.18, which is within the validated range
  of the crosswind tests but at the upper end.
- **Near-axial reversed flow on the retreating side** has not been tested at these advance
  ratios. The lifting-line solve uses a static polar with no reversed-flow extension (Section
  5.5.6); if retreating-blade AoA exceeds stall, the polar will clip but not diverge.
- **Tether attachment below the disk plane**: the tether hangs below the hub. At steep tilt the
  tether direction is nearly in-plane with the disk. This affects the force balance (the net
  aerodynamic force must balance tether tension + gravity component) but does not change the
  rotor aerodynamics internally; it is a coupling issue, not a VPM model issue.

### 11.3 Nonlinear SG6042 airfoil polar

The RAWES blade uses the **SG6042** low-Reynolds-number airfoil (Re ~ 127,000 at 50% span,
10 m/s wind). This is a cambered, low-Re profile with:
- Non-zero CL0 (zero-lift angle ≈ -4.92 deg)
- Nonlinear Cl vs alpha above ~8 deg
- Pronounced laminar separation bubble and hysteresis near stall
- Thin blade (10% thickness) -- different Cl_alpha slope from NACA profiles

The VPM currently uses a `LinearPolar` (constant Cl_alpha, symmetric stall). SG6042 at
operating Re requires a tabulated `XFoilPolar` or XFLR5 polar (Cl, Cd, Cm vs alpha at the
relevant Re range). `VpmRotor<P: Polar>` is now generic over the polar type (same as the
BEM-family models), so passing a `TabulatedPolar` works directly -- no internal re-sampling
step. The polar also needs to cover at least +-20 deg to handle the large AoA excursions in
autorotation.

### 11.4 Kaman servo-flap pitch actuation

The RAWES blades use a **Kaman-type trailing-edge servo flap** for pitch control
(`PitchActuation::ServoFlap` in `dynbem_rs`). The VPM implements this path in
`dynbem_rs/src/vpm_rotor.rs` as a per-blade feathering DOF integrated each
step. In servo mode, swashplate collective/cyclic are reinterpreted as
flap-deflection commands and the solved feathering angle replaces direct
swashplate pitch in the section AoA. The servo-flap path:
- Introduces a first-order feathering lag (~90 deg phase at 1/rev)
- Attenuates high-frequency cyclic authority
- Requires the `ServoFlapActuation` / `ServoFlapGeometry` parameters from the RAWES YAML
  (`beaupoil_2026.yaml`)

The VPM's `march_window` passes `fc.collective_rad`, `fc.tilt_lon`, `fc.tilt_lat` through
`cyclic_coeffs` to form the harmonic commands. In direct-mechanical mode those
commands map directly to blade pitch; in servo mode they drive flap deflection
and feathering dynamics. Servo mode is active whenever the rotor carries a
`ServoFlapActuation` (`pitch_actuation: servoflap`). The feathering ODE is

    I_theta * theta'' + C_theta * theta' + k_aero * theta = M_servo

where `C_theta` is the bearing damper (the only dissipation) and `k_aero` is the
**aerodynamic spring** from the AC offset,
`k_aero = 0.5*rho*omega^2*cl_alpha*ac_offset*Int(c r^2 dr)` (a pitch-up makes more
lift at the AC, a distance `ac_offset` aft of the feathering axis, giving a
nose-down restoring torque). With `ac_offset = 0` (axis at AC, Kaman ideal) the
damper alone sets the ~90 deg cyclic phase lag and DC trim comes from the blade
camber moment `blade_Cm_AC`. There is no artificial control-stiffness spring; all
constants are physical and measurable. (The alternative torsional-twist servo-flap
architecture -- elastic blade twist against spar stiffness -- is not modelled yet.)

### 11.5 Variable omega coupling to the spin ODE

RAWES omega evolves in time as a state variable (`omega_spin` ODE in `dynamics.py`). The VPM
`FlightCondition.omega_rad_s` is a fixed parameter per call. For proper coupling the caller
must update `omega_rad_s` each step from the physics ODE. The `AeroModel::step` interface in
`dynbem_rs` already supports this (it rebuilds `FlightCondition` from `RotorInputs` each call).
VPM's `step()` advances the wake by exactly one step of duration `dt` per call (the caller
supplies `dt` and drives the loop), so `dpsi = omega * dt` tracks the current omega directly.
A large omega change mid-simulation simply changes the azimuth swept per call; this is
acceptable provided `dt` is small relative to the spin acceleration timescale (which it is:
at `dt = 1/400 s` and 270 RPM, `dpsi` = ~4 deg per call).

### 11.6 State serialization

The RAWES simulation uses `to_dict()` / `from_dict()` on all state objects for
checkpointing (e.g. `steady_state_starting.json`). `VpmRotorState` carries a
`ParticleField` (SoA of f32 arrays) and a `gamma` matrix; neither is currently
serializable. For the VPM to be a drop-in for the `dynbem` factory, its state needs
`get_inflow()` / `set_inflow()` that round-trips through a flat vector, or alternatively a
bespoke `to_dict()` that serializes the raw particle arrays to JSON/binary.

### 11.7 Non-inertial hub acceleration

The RAWES hub orbits the tether anchor at ~0.105 rad/s with radius ~15 m → centripetal
acceleration ≈ 0.16 m/s^2 (≈ 1.7% g). The VPM treats the hub as inertially fixed; blade
element inertial corrections for hub acceleration are not modelled and are below the accuracy
target at current orbit rates, but would matter if the orbit radius or rate increases.

### Priority order for RAWES integration

| Priority | Item |
|----------|------|
| 1 | Tabulated SG6042 polar (Section 11.3) -- affects every operating point |
| 2 | Zero-torque autorotation omega iteration (Section 11.1) -- needed for any valid OP |
| 3 | ServoFlap parameter identification and validation for RAWES (Section 11.4) |
| 4 | Validated windmill / reel-in regime (Section 11.1) -- only needed for reel-in analysis |
| 5 | State serialization (Section 11.6) -- needed for checkpointing/factory drop-in |
