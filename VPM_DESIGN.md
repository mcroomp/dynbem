# Vortex Particle Method (VPM) -- design doc

Status: **experimental / in progress**. The wake engine and its unit tests
are implemented. Two rotor couplings exist: a pedagogical axial precursor
and the forward-flight, cyclic- and crosswind-capable coupling (Section
5.5). Neither is yet exposed as a first-class inflow model with the wake
carried as serialized state (Section 6, item 6).

This document describes what is built, the math it implements, the
conventions it follows, its measured performance, and the roadmap. It is
the design reference for anyone extending the VPM path. For the BEM /
Pitt-Peters / Oye models and the shared sign/frame conventions, see
[../AGENTS.md](../AGENTS.md) and [CLAUDE.md](CLAUDE.md).

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
    E --> F[FIFO-truncate oldest wake beyond n_wake_rev]
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

The tip edge carries $\Gamma_\text{tip}$ (the tip vortex, strongest); the
root edge carries $-\Gamma_\text{root}$.

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
$\psi_b(n) = n\,d\psi + b\,\tfrac{2\pi}{N_b}$. The carried state is the free
wake (a single shared particle field) plus the previous-step relaxed bound
circulation $\Gamma_i^{\,n-1}$ per blade, needed for the shed term. The
solution sought is the **periodic** limit cycle, not a fixed point: loads
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

Because $\mathbf{u}_\text{near}$ depends on $\Gamma$ and $\Gamma$ depends on
$\alpha$ (hence on $\mathbf{u}_\text{near}$), the last two lines are a coupled
fixed point. It is solved by relaxed iteration, seeded from the previous
step's circulation:
$\Gamma_i \leftarrow \Gamma_i + \text{relax}\,(\tfrac12 U c\,C_l - \Gamma_i)$,
repeated until $\Gamma$ stops changing (a few tens of iterations at most).
The far field is not re-probed inside the loop, so each iteration is only
the small near-field matrix-vector product plus a polar lookup -- negligible
against the $O(N^2)$ particle evaluation.

Lift comes from the **polar**, so (unlike a thin-airfoil vortex-lattice) the
bound vortex carries no self term; only the *trailing* wake enters
$\mathbf{u}_\text{near}$. This is the classical viscous lifting-line, in the
Bagai-Leishman free-wake / FLOWVPM tradition. Setting
`nonlinear_lifting_line = false` drops $\mathbf{u}_\text{near}$ and reduces
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

The tip edge carries $\Gamma_\text{tip}$ (the tip vortex), the root edge
$-\Gamma_\text{root}$.

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
near field $\mathbf{u}_\text{near}$: each edge $j$ carries a straight vortex
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

Blade flapping is not modelled (rigid blade). Hub moments are therefore the
full aerodynamic moments; the flap-hinge frequency-ratio moment reduction
used by the BEM models could be layered on later if a flapping rotor is
needed.

### 5.5.6 Simplifications and validity envelope

- **Rigid blade, no flap dynamics** -- consistent with the BEM convention
  in this repo; blade pitch sets local loading directly, no 90 deg
  precession.
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
outgrows the direct path -- see Section 7.1.

---

## 7. Performance (measured, single-core, release)

Direct O(N^2) evaluation, vectorized eight sources per lane (map over
targets on a single core, no multi-core parallelism).

| N | us / velocity eval | ms / RK2 advect step |
|---:|---:|---:|
| 500 | 487 | 1.2 |
| 2,000 | 13,793 | 24.2 |
| 5,000 | 75,822 | 153.9 |
| 10,000 | 335,229 | 1,133.8 |
| 20,000 | 1,550,205 | 5,066.5 |

Scaling follows the expected O(N^2) (5k -> 10k is 4.0x N -> 4.4x time; the
excess over 4x is cache pressure).

For reference, the per-step cost of each model at a non-axial operating point
(forward 12 m/s + 5 m/s edgewise crosswind) on the same machine. In a
time-domain simulation you call a BEM-family model once per integration step
(it returns a converged answer, 30 elem x 72 psi), whereas one VPM step
advances the free wake by one azimuth increment (an O(N^2) Biot-Savart probe
+ RK2 advect over the whole cloud):

| Model | cost / step | vs Oye |
|---|---:|---:|
| Oye (annular 2-stage filter) | 0.153 ms | 1x |
| Pitt-Peters (3-state L-matrix) | 0.165 ms | ~1.1x |
| Quasi-static BEM | 19.9 ms | ~130x |
| VPM free-wake, direct (run-avg march, N up to ~7,900) | ~0.47 s | ~3.1 x 10^3 |
| VPM free-wake, Barnes-Hut theta=0.5 (full-cloud step, N~7,900) | ~0.10 s | ~660x |

All four are measured single-core release-mode wall-clock. The BEM-family
rows are the per-call cost (min of many calls). The direct VPM row is the
run-average per step: a full `simulate()` at the default resolution (24
steps/rev, 4 wake revs, 6 settle revs, 20 stations) marches 168 steps in
79.7 s (min of 3), so ~0.47 s/step. That average understates the
steady-state step, because the wake grows from empty to the ~7,900-particle
cap over the first ~4 revs, so early steps are cheap and the full-cloud steps
at the end are the expensive ones.

The Barnes-Hut row is that expensive end-of-march step, at the ~7,900 cap: an
RK2 advect is two velocity evals, and at N=8,000 the tree (`theta = 0.5`)
runs each in ~50.5 ms vs the direct path's ~406 ms (Section 7.1), so the
full-cloud step falls from ~0.81 s to ~0.10 s -- an ~8x cut on exactly the
steps that dominate the march. The tree engages only above `bh_min_particles`
and is off by default, so the direct run-average is still the default-config
figure; a tree-enabled march sees the same empty-to-full growth curve, so its
run-average would fall by a similar factor once the cloud is large.

Oye and Pitt-Peters are within ~10% of each other; both evaluate a fixed set
of algebraic inflow relations per (element, azimuth). The quasi-static BEM is
~130x slower because each (element, azimuth) station runs an iterative
root-finder (Brent) on the momentum / blade-element balance, and non-axial
wind makes it worse, not better: the swept azimuth exposes more sections to
the turbulent-wake / reversed-flow branches where the solver has to iterate,
so the per-call cost rises above its axial figure rather than falling. The
VPM is ~3 orders of magnitude beyond even the BEM per step (and ~5 orders per
converged operating point): it is not an inflow model but a wake-fidelity
tool that resolves the actual wake geometry, paying the O(N^2) particle cost
at every one of the ~170 marched steps.

**Per operating point** the gap is larger still: the BEM-family models are
already converged after one step, while the VPM needs the full ~170-step
march to reach a periodic wake -- Pitt-Peters ~0.17 ms vs the VPM at 79.7 s,
roughly 5 x 10^5 x. These are worst-case figures: a
single core, per-lane vectorization only, and the direct O(N^2) path. The
target loop is independent, so spreading it across cores should cut the
constant by something approaching the core count until memory bandwidth
limits it, and the Barnes-Hut tree (Section 7.1) changes the N scaling
itself -- ~8x at N=8,000 and widening with N.

Simulation guidance (bench rotor, omega=120, one rev = 52 ms): a
Tier-A config (~1,200 particles, 24 steps/rev, 3 revs) costs ~7 ms/step,
about 3x slower than real time on one core; spreading the target loop
across cores would be expected to bring it near real time. A Tier-B config
(~4,300 particles, 36 steps/rev, 4 revs) is ~115 ms/step, an offline
setting.

### 7.1 Barnes-Hut tree evaluator (O(N log N))

**Why the tree is faster (the physical picture).** In the direct evaluation
every particle induces a velocity on every other particle. With $N$ particles
that is $N \times N$ pairwise Biot-Savart evaluations: double the size of the
wake and the work quadruples. Most of that work is spent on particles that
are far apart, and there the detail is wasted -- a compact clump of vortex
particles seen from far away induces very nearly the same velocity as a
*single* vortex carrying the clump's total strength, because the regularized
kernel's far field is just the ordinary Biot-Savart law of the net
circulation ($K \to 1/r^3$ once $r \gg \sigma$). The tree exploits this: it
sorts the wake into a nested hierarchy of boxes, and when it evaluates the
velocity at a point it sums the *nearby* vorticity particle-by-particle (where
the detail matters) but replaces each *distant* box with one equivalent
lumped vortex (where it does not). Each evaluation point only ever has to
touch a handful of near particles plus a modest number of far lumps -- a count
that grows like $\log N$ rather than $N$ -- so the total cost grows like
$N \log N$ instead of $N^2$. The bigger the wake, the more this pays off (see
the table below: ~5x faster at 4,000 particles, ~13x at 16,000).

How "far enough to lump" is decided is the one knob: a box is treated as a
single lumped vortex when it looks small from the evaluation point, i.e. when
its width divided by its distance is below the opening angle `theta`. Smaller
`theta` lumps less aggressively (more boxes summed in full) -- more accurate,
slower; larger `theta` lumps more -- faster, coarser. `theta ~ 0.5` keeps the
approximation error to a few percent of the peak velocity.

The direct evaluator above is O(N^2). For large wakes there is a
monopole Barnes-Hut path (`induced_at_points_bh`, `induced_velocities_bh`,
`advect_rk2_bh`) that groups distant sources into a single super-particle
and evaluates the far field with the *same* SIMD kernel as a real particle
-- the regularized kernel's far field is the singular Biot-Savart law
($K \to 1/r^3$ for $r \gg \sigma$), so a far cell is just a source carrying
the lumped strength $\sum \boldsymbol{\alpha}$ placed at the
$|\boldsymbol{\alpha}|$-weighted cell centre with a representative core size.
The tree is flattened into a compact pre-order node array with escape
pointers (no child arrays, no recursion in the walk), and the particles are
reordered at build time into leaf-contiguous, eight-padded blocks. Per target
the stackless walk runs the SIMD kernel directly over each accepted near
leaf's packed block -- there is no per-target near-field gather -- while the
accepted far cells are batched into one small padded buffer and kernelled in
a single tail pass. Traversal is O(N log N) of light scalar work; the heavy
arithmetic stays vectorized.

Precisely, a cell of width $s$ seen at distance $d$ is accepted as far when
$s / d < \theta$, so `theta = 0` lumps nothing and recovers the exact direct
sum. The opening angle is exposed through `VpmRotorConfig` (`barnes_hut`,
`bh_theta`, `bh_min_particles`) and is **off by default** -- the direct sum
is the reference path; the tree engages only when `barnes_hut` is set and the
wake reaches `bh_min_particles` (below that the tree-build overhead is not
worth it). A dipole term could be added later for more accuracy at a given
`theta`.

Measured direct O(N^2) vs tree (`theta = 0.5`) on the profiler's wake-like
cloud (`bh_profile`, single core, release, min-noise 12 s runs). Both paths
go through the same measurement harness (`bh_profile <N> 0 <s>` for direct,
`bh_profile <N> 0.5 <s>` for the tree), so the speedup is an apples-to-apples
back-to-back ratio on the same machine state:

| N | direct O(N^2) [ms] | tree theta=0.5 [ms] | speedup |
|---:|---:|---:|---:|
| 4,000 | 102.7 | 22.0 | 4.7x |
| 8,000 | 406.2 | 50.5 | 8.0x |
| 16,000 | 1,566.9 | 117.9 | 13.3x |

The direct column is clean O(N^2) (each doubling of N is ~4x time); the tree
column grows ~2.3x per doubling (close to O(N log N)), so the crossover keeps
widening with N -- the tree is worth more the larger the wake. The absolute
ms figures are harness/machine specific (this cloud, this core), but the
relative speedup is robust because direct and tree are measured in the same
session.

---

## 8. Validation status

### 8.1 Unit tests (all passing)

- **Far field** matches $A / (4\pi d^2)$ within 1%, direction $+y$ for $+z$
  vorticity.
- **Self-induction** of a lone particle is zero.
- **Vectorized vs reference**: the vectorized path matches the f64 scalar
  reference on a random cloud to $< 10^{-3}$ relative.
- **Vortex-ring self-propagation**: a thin ring self-propagates along its
  axis within 30% of the Kelvin thin-core speed
  $U = \dfrac{\Gamma}{4\pi R}\left[\ln\dfrac{8R}{a} - \dfrac{1}{4}\right]$,
  transverse velocity negligible.
- **Advection**: one RK2 step convects a ring downstream.
- **Barnes-Hut vs direct**: the monopole tree (`theta = 0.5`) matches the
  direct sum to $< 5\%$ of peak velocity on a random cloud, converges to the
  direct result as `theta -> 0`, and its RK2 advect step moves a thin ring
  downstream to within 5% of the direct integrator.

### 8.2 Thrust/torque vs BEM (axial precursor)

Rotor R=1 m, 2 blades, collective 8 deg, omega=120 rad/s.

| Regime | BEM thrust [N] | VPM thrust [N] | delta |
|---|---:|---:|---:|
| hover | 215.6 | 238.2 | +10.5% |
| climb +3 | 183.2 | 199.9 | +9.1% |
| climb +6 | 143.4 | 159.2 | +11.0% |
| descent -3 | 533.4 | 230.5 | -56.8% |
| descent -6 | 407.9 | 375.6 | -7.9% |

Read:

- **Hover and climb agree to ~10%**, with the correct trend (thrust falls
  as climb rate rises). The consistent +10% bias has a plausible
  explanation: BEM has Prandtl tip-loss on (this axial precursor does not),
  plus coarse azimuth resolution, a hand-picked `sigma`, and no within-step
  self-induction -- all pushing the same way. Read as a consistency check
  this is about what one would expect from a correctly wired coupling; it is
  not an accuracy claim, and the ~10% is not decomposed into those
  contributions. The forward-flight module (Section 5.5) adds the implicit
  lifting-line self-induction the precursor lacks: on a comparable hover case
  (24 steps/rev) it lowers the analogous bias from about +12% (pointwise) to
  about +9% (uniform lifting-line). Adding tip-clustered (cosine) spanwise
  spacing plus a locally scaled core closes the rest -- the two levers
  compose (a local core does nothing on a uniform grid, but once the spacing
  clusters at the tip it stops over-smoothing the tip vortex), bringing the
  measured hover bias to under +1% (216.5 N vs a 214.9 N BEM anchor).
- **Descent diverges**, in the regime where BEM switches to its empirical
  VRS/windmill model (descent -3 has descent/v_h ~ 0.55, inside the 0-2 VRS
  band; BEM thrust spikes and its torque goes negative = windmilling). VPM
  gives a smooth, finite result with torque approaching zero (which is what
  autorotation would look like) at descent -6 -- but see the caveat.

**Caveat:** the VPM descent numbers are not validated. The trailing-only
wake, 4-rev truncation, and coarse resolution under-resolve the
recirculating VRS wake. Descent here means "the two methods diverge in the
regime where BEM is empirical", not "VPM is right there". Hover and climb
are the comparisons worth trusting; descent is not.

### 8.3 Forward-flight coupling (acceptance tests, all passing)

These are the acceptance tests for the Section 5.5 module. They are run at
the coarse resolution preset (Section 5.5.6), so they assert **directional /
qualitative** physics (signs, monotonicity, boundedness) rather than tight
quantitative agreement -- the quantitative anchor is the axial reduction.

- **Axial reduction**: with zero cyclic and zero in-plane wind the shed
  term vanishes and the model must collapse onto the axial case -- hover
  thrust within 30% of the quasi-static BEM, and hub moments
  $|M_x|, |M_y| < 0.10\,T R$ (axisymmetry). This is the one quantitative
  check and the regression guard on the coupling wiring.
- **Collective**: $dT/d(\text{collective}) > 0$.
- **Longitudinal cyclic**: `tilt_lon > 0` $\Rightarrow M_y < 0$ and
  $|M_y|$ grows from the hover baseline (correct cyclic phasing / hub-moment
  sign).
- **Lateral cyclic**: `tilt_lat > 0` $\Rightarrow M_x > 0$.
- **Crosswind**: 8 m/s edgewise inflow -- thrust stays finite and positive,
  torque finite, the advancing/retreating asymmetry develops a non-zero hub
  moment, and the wake centroid convects downstream ($x_\text{centroid} > 0$).
  This exercises the shed term and the free wake-skew, and guards against the
  numerical blow-up that a wrong-signed shed contribution would produce.

**Not yet covered:** quantitative forward-flight loads against measured
data (e.g. a trimmed rotor at a published advance ratio), reversed-flow
regime, and any BVI-sensitive case. The forward-flight module is validated
for *trend correctness and stability*, not yet for absolute accuracy.

---

## 9. Integration points with the codebase

- Frame / sign conventions: NED, CCW-from-above -- see AGENTS.md. The
  coupling's r_hat / t_hat match the shared BEM kinematics.
- Geometry / polars: the shared radial grid, polar table, and airfoil-polar
  interface used by the BEM models.
- When promoted to a first-class model (Section 6, item 6): the wake becomes
  carried, serialized inflow state, and the model is registered in the model
  factory alongside the others and exposed through the Python bindings.

---

## 10. References

- Winckelmans, G. and Leonard, A. (1993). "Contributions to Vortex Particle
  Methods for the Computation of Three-Dimensional Incompressible Unsteady
  Flows." J. Comput. Phys. 109(2). -- the algebraic regularization kernel.
- Cottet, G.-H. and Koumoutsakos, P. (2000). "Vortex Methods: Theory and
  Practice." Cambridge. -- general VPM theory, stretching, diffusion.
- Alvarez, E. J. and Ning, A. (2020). "High-Fidelity Modeling of Multirotor
  Aerodynamic Interactions for Aircraft Design." AIAA J. -- reformulated
  VPM / FLOWVPM, modern rotorcraft usage.
- Leishman, J. G. (2006). "Principles of Helicopter Aerodynamics," 2nd ed.
  Cambridge. -- lifting line, Kutta-Joukowski, VRS, tip vortices.
- Bagai, A. and Leishman, J. G. (1995). "Rotor Free-Wake Modeling using a
  Pseudo-Implicit Technique." J. Aircraft. -- rotor free-wake conventions.
- Landgrebe, A. J. (1972). "The Wake Geometry of a Hovering Helicopter Rotor
  and its Influence on Rotor Performance." JAHS. -- prescribed / trailing
  wake geometry.
