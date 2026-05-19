# Peters Nikolsky Lecture — JAHS 54(1):011001 (2009)

Peters, David A. (2009). "How Dynamic Inflow Survives in the Competitive
World of Rotorcraft Aerodynamics: The Alexander Nikolsky Honorary Lecture,"
*Journal of the American Helicopter Society* 54(1):011001.
DOI: 10.4050/JAHS.54.011001.

Source PDF: `Peters_Nikolsky_Lecture_JAHS_2009.pdf` (this folder).

This is the **canonical reference** for the classical Pitt-Peters 3-state
dynamic-inflow model, by the original developer himself. Cited by Pitt &
Peters (1981, refs 87/103 in this paper) and Peters & HaQuang (1988, which
is a practical-applications follow-up to the same model).

## Sign / convention summary

Peters' conventions in this paper:
- **ψ origin:** Standard helicopter — `ψ = 0 over the tail`, increasing
  in the direction of rotation. For a US (CCW-from-above) rotor this puts
  ψ = π/2 on the **right** (advancing for +X flight).
- **Inflow harmonic decomposition** (Eq 11):
  `υ = ν_0 + ν_s · r · sin(ψ) + ν_c · r · cos(ψ)`
- **State ordering:** `(ν_0, ν_s, ν_c)` — uniform, sin-coefficient, cos-coefficient.
- **C_T = thrust, C_L = roll moment, C_M = pitch moment** (aerodynamic
  loading contributions only; nomenclature p.1 of the lecture).

## The canonical equations

### Eq 7 — Pitt-Peters ODE

    [M]{ν̇_0; ν̇_s; ν̇_c} + V · [L]⁻¹ · {ν_0; ν_s; ν_c} = {C_T; −C_L; −C_M}

Note the **negative signs** on C_L and C_M in the forcing.

### Eq 8 — Mass-flow parameter

    V = (μ² + (λ+ν)(λ+2ν)) / √(μ² + (λ+ν)²)

where μ = advance ratio, λ = climb ratio (axial freestream non-dim),
ν = average induced flow ratio (= ν_0 at the disk).

Limits:
- Hover (μ=0, λ=0): V = 2·ν_0
- High-speed (μ ≫ λ+ν): V → μ
- Climb (λ ≫ μ): V → λ + 2ν (approx, depends on relative magnitudes)

This V is *not* the same as `µ_T = √(μ² + (λ+ν)²)` — they differ by
`ν·(λ+ν)/µ_T`, significant only at low speed.

### Eq 9 — Apparent-mass matrix

    [M] = | 8/(3π)    0          0         |
          | 0         16/(45π)   0         |
          | 0         0          16/(45π)  |

All diagonal entries positive. These give the time constants
`τ_0 = 8R/(3π·V_T)` and `τ_cs = 16R/(45π·V_T)` (with R for dimensionalization,
V_T as effective mass-flow speed in m/s).

### Eq 10 — Influence-coefficient matrix L

With X = tan(χ/2) (wake-skew parameter):

    [L] = | 1/2              0              −15π·X/64 |
          | 0                2·(1+X²)        0         |
          | 15π·X/64         0               2·(1−X²)  |

The off-diagonal (0,2) ↔ (2,0) is **anti-symmetric** (opposite signs).
This couples C_M forcing → ν_0 response **and** C_T forcing → ν_c response
— the latter is the famous Pitt-Peters "thrust-to-tilt" cross-coupling
that produces Glauert wake-skew naturally.

In limiting forms:
- Axial (X=0, χ=0): L = diag(1/2, 2, 2). No cross-coupling. All cyclic
  decouples; ν_0 = C_T/(2V) (matches momentum theory with Peters' V).
- Pure edgewise (X=1, χ=π/2): L[1,1] = 4, L[2,2] = 0. Longitudinal
  inflow harmonic vanishes (no cyclic damping in pure edgewise).

### Eq 11 — Inflow distribution

    υ(r, ψ) = ν_0 + ν_s · r · sin(ψ) + ν_c · r · cos(ψ)

where r is non-dim radius (r/R), ψ is azimuth in Peters' convention.

## How this maps to our code

Our convention differs from Peters' on ψ origin: **we use ψ = 0 at +X
(hub-frame nose), CCW from above**. See `[CLAUDE.md](../../CLAUDE.md)`
"Rotor rotation direction".

ψ-rotation by π gives:
- our λ_c = − Peters' ν_c
- our λ_s = − Peters' ν_s
- our C_L (roll moment in our hub frame) = + Peters' C_L (same physical sense)
- our C_M (pitch moment) = + Peters' C_M

After translation, our steady-state targets become:

    λ_0_ss = C_T/(2·V) + (15π·X/64) · C_M_hub / V
    λ_c_ss = (−15π·X/64) · C_T / V  +  4·cos(χ)/(1+cos χ) · C_M_hub / V
    λ_s_ss =                          +  4/(1+cos χ)       · C_L_hub / V

(Using `2(1−X²) = 4·cos(χ)/(1+cos χ)` and `2(1+X²) = 4/(1+cos χ)`.)

## Open question on V vs µ_T

Our code currently uses `µ_T = √(µ² + λ_total²)` (classical Glauert
mass-flow) rather than Peters' V (Eq 8). They agree in high-speed
forward flight but differ in hover by a factor of 2 (V = 2·µ_T in hover).

Trade-off:
- µ_T matches **classical Glauert** hover inflow `ν_0 = √(C_T/2)`.
- Peters' V gives `ν_0 = √(C_T)/2` in hover (smaller by √2).

The √2 difference suggests there's a normalization convention I haven't
fully pinned down — Peters' C_T might be 2× the standard
`T/(ρ·A·(ΩR)²)`, or his ν might refer to a different point in the wake.
Our code is internally consistent with `µ_T` and reproduces classical
Glauert; replacing µ_T → V would require verifying against an
independent hover dataset (e.g., Castles-Gray) to see which gives the
better fit. **Defer this swap until validation data is in hand.**

The L-matrix *structure* (diagonal + anti-symmetric off-diagonal) is
matched exactly.

## Things this paper validates in our implementation

- Off-diagonal cross-coupling magnitude `15π/64 · tan(χ/2)` ✓
- Diagonal cyclic gains: `4·cos(χ)/(1+cos χ)` for λ_c, `4/(1+cos χ)` for λ_s ✓
- Apparent-mass matrix `M = diag(8/(3π), 16/(45π), 16/(45π))` — matches
  our `τ_0 = 8R/(3π·V_T)`, `τ_cs = 16R/(45π·V_T)` ✓
- The forcing vector `{C_T, −C_L, −C_M}` sign convention ✓ (we handle
  sign translation via ψ rotation)
- Inflow distribution as uniform + linear in r, harmonic in ψ ✓

## Things this paper flags as potential refinements

- Use V (Eq 8) instead of µ_T for mass-flow scaling — would give
  different hover cyclic magnitudes.
- The "C" coupling matrix in Eq 1 (`[M]{dν/dt} + [C]{ν} = {F}` general
  form) — for the 3-state Pitt-Peters this is `V·[L]⁻¹`, but generalized
  finite-state inflow (Peters-He, ref 18+) extends to higher harmonics.
  Future Level-3 work.
