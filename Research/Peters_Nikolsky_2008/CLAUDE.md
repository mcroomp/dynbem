# Peters Nikolsky Lecture -- JAHS 54(1):011001 (2009)

Peters, David A. (2009). "How Dynamic Inflow Survives in the Competitive
World of Rotorcraft Aerodynamics: The Alexander Nikolsky Honorary Lecture,"
*Journal of the American Helicopter Society* 54(1):011001.
DOI: 10.4050/JAHS.54.011001.

Source PDF: `Peters_Nikolsky_Lecture_JAHS_2009.pdf` (this folder).
Plain-text extraction: `Peters_Nikolsky_Lecture_JAHS_2009.txt` -- but the
EQUATIONS in that .txt are mangled by the PDF text extractor. Read the PDF
images for formulas; the math below was transcribed directly from rendered
page images (article pages 011001-5 and 011001-6).

This is the **canonical reference** for the classical Pitt-Peters 3-state
dynamic-inflow model, by the original developer himself. It narrates the
whole lineage (Amer -> Sissingh -> Curtiss -> Carpenter/Fridovitch ->
Ormiston/Peters -> Pitt -> He -> Prasad -> Morillo/Yu/Hsieh -> Makinen)
and gives the Pitt model in closed form (Eqs 7-11).

--------------------------------------------------------------------------

## 1. Conventions used in the paper

- **psi origin:** standard helicopter -- `psi = 0 over the tail`,
  increasing in the direction of rotation. For a US (CCW-from-above)
  rotor this puts psi = pi/2 on the **right** (advancing side for +X
  flight). This is OPPOSITE to our code convention (psi = 0 at +X/nose);
  see section 7 for the sign translation.
- **Inflow harmonic decomposition (Eq 11):**
  `upsilon = nu_0 + nu_s * r * sin(psi) + nu_c * r * cos(psi)`
  with `upsilon` the total induced velocity at a disk point (normalized
  on tip speed) and `r` the nondimensional radius (r/R).
- **State ordering:** `(nu_0, nu_s, nu_c)` -- uniform, side-to-side
  (sin), fore-to-aft (cos).
- **Forcing:** `C_T` thrust, `C_L` roll moment, `C_M` pitch moment
  (aerodynamic loading contributions only).
- **Nondimensional time:** `tau = Omega*t`; `()* = d()/d_tau`.
- **mu** advance ratio, **lambda** climb ratio (axial freestream
  non-dim), **nu** average induced flow ratio (= nu_0 at the disk),
  **chi** wake skew angle, **X = tan(chi/2)** the wake-skew parameter
  (X in [0, 1] as chi goes 0 deg in hover to 90 deg in edgewise flow).

--------------------------------------------------------------------------

## 2. The genus: general dynamic-wake form (Eq 1)

Any dynamic-wake / dynamic-inflow model has the finite-state ODE form

    [M] d{v_n}/dt + [C]{v_n} = {F_m}                                   (1)

- `v_n` are flow-field states,
- `[M]` apparent-mass matrix,
- `[C]` matrix of influence coefficients,
- `{F_m}` depends on the blade loading.

For the 3-state Pitt model, `[C] = V*[L]^-1` (see Eq 7).

## 3. Lock number and Curtiss equivalent lift-curve slope (Eqs 2-3)

    gamma = rho * a * c * R^4 / I_y                                    (2)

(Lock number: ratio of aerodynamic to inertial forces; `a` = lift-curve
slope, `c` = chord, `R` = radius, `I_y` = blade flap inertia.)

Curtiss/Shupe quasi-steady result -- inflow reduces the effective
lift-curve slope (and hence the effective Lock number):

    a*/a = gamma*/gamma = [1 + sigma*a/(8V)]^-1                        (3)

where `sigma` is solidity and `V` the mass-flow parameter. This is the
"lift-deficiency from inflow" that agrees with the Loewy function at
integer-multiple frequencies. (Not used directly in our BEM code, but it
is WHY dynamic inflow matters -- it knocks ~50% off blade aerodynamic
effectiveness, doubling roll/pitch damping.)

--------------------------------------------------------------------------

## 4. Hover apparent-mass equations (Eqs 4-6)

Derived by Ormiston/Peters using the apparent mass + apparent inertia of
an impermeable disk (Carpenter-Fridovitch) for the time delays:

    [8/(3*pi)]   d_nu_0/d_tau + [2V]   nu_0 =  C_T                     (4)
    [16/(45*pi)] d_nu_s/d_tau + [V/2]  nu_s = -C_L                     (5)
    [16/(45*pi)] d_nu_c/d_tau + [V/2]  nu_c = -C_M                     (6)

These are the hover special case (X = 0, chi = 0) of Eqs 7-10:
with X=0, `[L] = diag(1/2, 2, 2)`, so `V*[L]^-1 = diag(2V, V/2, V/2)`,
matching the stiffness terms above exactly. Internal consistency check
that the L-matrix reduces to the impermeable-disk hover model.

--------------------------------------------------------------------------

## 5. The Pitt model (Eqs 7-11) -- THE canonical formulation

### Eq 7 -- the ODE

    [M] {nu_0*; nu_s*; nu_c*} + V*[L]^-1 {nu_0; nu_s; nu_c}
        = {C_T; -C_L; -C_M}                                           (7)

Note the **negative signs** on C_L and C_M in the forcing.

### Eq 8 -- mass-flow parameter V (wake contraction)

    V = ( mu^2 + (lambda + nu)(lambda + 2*nu) )
        / sqrt( mu^2 + (lambda + nu)^2 )                              (8)

written in terms of advance ratio mu, climb ratio lambda, and average
induced flow nu. NOTE: `lambda` (climb) and `nu` (induced) are DISTINCT
variables; in hover lambda = 0 but nu != 0. Limits:
- Hover (mu=0, lambda=0): `V = (nu)(2nu)/sqrt(nu^2) = 2*nu`.
- High speed (mu >> lambda+nu): `V -> mu`.

### Eq 9 -- apparent-mass matrix (impermeable disk)

    [M] = diag( 8/(3*pi), 16/(45*pi), 16/(45*pi) )                    (9)

All diagonal, all positive. With `tau = Omega*t`, these give the
dimensional time constants (M entry divided by the corresponding
diagonal of `V*[L]^-1`; see section 6).

### Eq 10 -- influence-coefficient matrix L

With `X = tan(chi/2)`:

    [L] = |  1/2          0            -(15*pi/64)*X |
          |  0            2(1 + X^2)    0            |
          |  (15*pi/64)*X 0            2(1 - X^2)    |                 (10)

- The (0,2)/(2,0) off-diagonal is **anti-symmetric** (opposite signs).
  This couples `C_M -> nu_0` AND `C_T -> nu_c` -- the latter is the
  Pitt "thrust-to-tilt" cross-coupling that yields Glauert wake skew
  naturally.
- Limiting forms:
  - Axial (X=0, chi=0): `L = diag(1/2, 2, 2)`. No cross-coupling.
  - Edgewise (X=1, chi=90 deg): `L[1,1]=4`, `L[2,2]=0`. The fore-aft
    harmonic loses its restoring term.

### Eq 11 -- inflow distribution

    upsilon(r, psi) = nu_0 + nu_s*r*sin(psi) + nu_c*r*cos(psi)        (11)

Uniform + linear-in-r, first harmonic in psi.

--------------------------------------------------------------------------

## 6. Steady state, inverse-L closed form, and time constants

At steady state (`d/d_tau = 0`), Eq 7 gives
`{nu} = (1/V) [L] {C_T; -C_L; -C_M}` because `(V*L^-1)^-1 = L/V`. So the
**steady-state targets** are just `[L]` applied to the forcing, over V:

    nu_0_ss = ( (1/2)*C_T          - (15*pi/64)*X*(-C_M) ) / V
    nu_s_ss = ( 2(1+X^2)*(-C_L) )                          / V
    nu_c_ss = ( (15*pi/64)*X*C_T  + 2(1-X^2)*(-C_M) )      / V

This is exactly the structure our code implements (after the sign
translation in section 7).

Time constants follow from `[M]` divided by the diagonal of `V*[L]^-1`.
In dimensional form (with R for length and an effective mass-flow speed
`V_T = V*Omega*R` in m/s):

    tau_0  = 8*R / (3*pi*V_T)        (uniform / collective mode)
    tau_cs = 16*R / (45*pi*V_T)      (both cyclic modes)

These match what `pitt_peters.rs` uses (`v_mf` plays the role of `V_T`).

### Eq 12 -- harmonic time-constant scaling (informational)

    Dynamic wake:    T = 0.75 / (m + 3/2)
    Loewy at r=3/4:  T = 0.75 / m

For the lowest mode of the m-th harmonic. m=0 gives T=0.50, m=1 gives
T=0.30 (dynamic wake). Confirms the apparent-mass time constants against
Loewy theory; not used numerically in our 3-state code.

--------------------------------------------------------------------------

## 7. Mapping to our code (sign / axis translation)

Our convention: **psi = 0 at +X (hub-frame nose), CCW from above** (see
top-level `AGENTS.md` "Rotor rotation direction"). Peters uses psi = 0 at
the tail. A psi-rotation by pi gives:

- our `lambda_c = -` Peters `nu_c`
- our `lambda_s = -` Peters `nu_s`
- our `C_M_hub  = +` Peters `C_M`  (pitch)
- our `C_L_hub  = +` Peters `C_L`  (roll)

After translation, with `X = tan(chi/2)`,
`l_off = (15*pi/64)*X`, `l_cc = 2(1-X^2) = 4*cos(chi)/(1+cos chi)`,
`l_ss = 2(1+X^2) = 4/(1+cos chi)`, the steady-state targets become:

    lambda_0_ss = C_T/(2V)        + l_off * C_M_hub / V
    lambda_c_ss = (-l_off * C_T   + l_cc  * C_M_hub) / V
    lambda_s_ss = (                  l_ss  * C_L_hub) / V

which is what `compute_forces` evaluates (in wind axes -- see flag D).
The `-l_off*C_T` term is the wake-skew cross-coupling; do NOT also add a
closed-form Glauert tilt (it would double-count).

--------------------------------------------------------------------------

## 8. Mass-flow: Peters V (Eq 8) vs our Glauert v_mf -- RECONCILED

Our code uses the **classical Glauert** mass-flow

    v_mf = sqrt(mu^2 + lambda_total^2)      (lambda_total = climb + nu)

instead of Peters' Eq 8 `V`. They agree at high speed (both -> mu) but
differ at low speed. The correct hover comparison, from the paper's own
algebra:

- **Peters V (Eq 8), hover:** `V = 2*nu`. Steady Eq 4: `2V*nu_0 = C_T`
  => `4*nu_0^2 = C_T` => `nu_0 = sqrt(C_T/4) = (1/2)*sqrt(C_T)`.
- **Glauert v_mf, hover:** `v_mf = nu_0`. Steady: `C_T/(2*v_mf) = nu_0`
  => `2*nu_0^2 = C_T` => `nu_0 = sqrt(C_T/2)`.
- **Ratio Glauert/Peters = sqrt(2)** (~1.414), and `V = 2*v_mf` in hover.

So the mass-flow *parameter* differs by a factor 2 in hover; the
resulting uniform inflow differs by sqrt(2). The L-matrix STRUCTURE
(diagonal + anti-symmetric off-diagonal) is identical either way.

Note: in Peters' Eq 8, `lambda` is the CLIMB ratio (= 0 in hover) and
`nu` is the induced flow -- they are DISTINCT variables. Conflating them
(setting `lambda = nu`) would give a spurious `V = 3*nu` / sqrt(3) hover
ratio; the correct hover limit is `V = 2*nu`, ratio sqrt(2).

### Why we keep Glauert anyway

- It reproduces the classical Glauert hover inflow `nu_0 = sqrt(C_T/2)`.
- Li (UC Davis MS, 2020), a finite-state Peters-He implementation, also
  falls back to the Glauert form `lambda_0 = C_T/(2*sqrt(mu^2+lambda^2))`
  for the uniform-inflow trim (his Eq 14), corroborating the choice.
- The 46-point empirical regression suite (Castles-Gray hover + descent)
  is calibrated to the Glauert form.

Swapping to Peters' V would rescale every hover time constant and inflow
magnitude and require re-validation; defer until there is hover inflow
data that specifically discriminates the two.

--------------------------------------------------------------------------

## 9. What this paper validates in our implementation

- Off-diagonal cross-coupling magnitude `15*pi/64 * tan(chi/2)`  [OK]
- Diagonal cyclic gains `2(1+X^2)` (nu_s) and `2(1-X^2)` (nu_c)   [OK]
- Apparent-mass `M = diag(8/(3pi), 16/(45pi), 16/(45pi))`, giving
  `tau_0 = 8R/(3*pi*V_T)`, `tau_cs = 16R/(45*pi*V_T)`             [OK]
- Forcing sign convention `{C_T, -C_L, -C_M}` (we translate via the
  psi rotation in section 7)                                       [OK]
- Inflow distribution uniform + linear-in-r, first harmonic in psi [OK]

## 10. Refinements mentioned (NOT in our 3-state model)

These are later extensions; our code is the 3-state Pitt level only.

- **Peters-He generalized finite-state** (Refs 20-21): all harmonics +
  radial shape functions; contains Pitt as the 3-state truncation. Our
  Oye annular model is a different (DBEMT-style) route to higher fidelity.
- **Wake curvature** (Prasad/Zhao, Refs 24-26): adds pitch-rate terms to
  `[L]` to fix off-axis (pitch->roll) coupling sign. Not modelled.
- **Off-rotor flow field** (Morillo/Yu/Hsieh, Refs 27-29): velocity
  potentials as states -> all 3 flow components on/off disk; ground
  effect. Not modelled.
- **Swirl velocity** (Makinen, Ref 30): mass-matrix correction for swirl
  kinetic energy, needed for high-inflow propeller power optimization.
  Not modelled.
