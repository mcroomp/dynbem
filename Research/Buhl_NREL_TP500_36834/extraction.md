# Buhl (2005) — NREL/TP-500-36834

**Title:** A New Empirical Relationship between Thrust Coefficient and Induction Factor for the Turbulent Windmill State
**Author:** Marshall L. Buhl, Jr.
**Report:** NREL/TP-500-36834, August 2005
**PDF:** 36834.pdf (12 pages, 315 KB)
**Source:** https://docs.nrel.gov/docs/fy05osti/36834.pdf

---

## What this paper is (and is not)

This is a **purely mathematical derivation** — it contains **no new experimental data**.
The experimental data referenced is Lock et al. (1926) autogyro/airscrew tests, digitized
from Figure 3.16 of Burton et al. *Wind Energy Handbook* (2001). That data has very large
scatter. Buhl explicitly notes: *"it does not account for the spread in the data, which is a
problem that needs further investigation."*

---

## Background: the problem

Classical momentum theory gives:

    CT = 4a(1 − a)                                          (1)

This is valid for induction factor a < ~0.4. For a > 0.5 it violates the derivation assumptions
(the "turbulent windmill state"). Glauert (1926) fit experimental data with an empirical parabola:

    CT = 0.6427 · (0.889 − (a − 0.143)² / 0.0203)          (2)  [Glauert/Eggleston form]

When you include tip/hub losses F in the classical equation:

    CT = 4Fa(1 − a)                                         (3)

a **numerical gap** opens between this curve and Glauert's Eq (2) for F < 1. BEM iterative
solvers can get stuck in this gap (no valid root).

---

## Buhl's derivation

Fit a quadratic `CT = b0 + b1·a + b2·a²` with three constraints:
1. Same value as Eq (3) at a = 0.4:  `CT = 0.96F`
2. Same slope as Eq (3) at a = 0.4:  `dCT/da = 0.8F`
3. CT = 2.0 at a = 1.0

Solving gives:

    b2 = 50/9 − 4F
    b1 = 4F − 40/9
    b0 = 8/9

**Buhl equation** (apply for a > 0.4):

    CT = 8/9 + (4F − 40/9)·a + (50/9 − 4F)·a²             (18)

This curve:
- Touches `CT = 4Fa(1−a)` tangentially at a = 0.4 (no gap)
- Passes through CT = 2.0 at a = 1.0
- Reduces to the F=1 case: `CT = 8/9 − 8a/9 + 50a²/9 − 4a²`... wait:
  - F=1: `CT = 8/9 + (4 − 40/9)a + (50/9 − 4)a²`
         `= 8/9 − 4a/9 + 14a²/9`

---

## Transition point

The switch from Eq (3) to Eq (18) occurs at **a = 0.4**:

    if a <= 0.4:  CT = 4Fa(1 − a)
    if a >  0.4:  CT = 8/9 + (4F − 40/9)·a + (50/9 − 4F)·a²

In terms of the rotor inflow ratio λ_r used in this project:
`a = 1 − λ_r/λ_c` (wind turbine convention, where λ_c < 0 for upward wind in NED),
so a = 0.4 corresponds to a moderately loaded turbine condition.

---

## Validation data

Lock, C.N.H.; Bateman, H.; Townsend, H.C.H. (1926). *An Extension of the Vortex Theory of
Airscrews with Applications to Airscrews of Small Pitch, Including Experimental Results.*
ARC R&M No. 1014. London: HMSO.

- Autogyro/airscrew tests, 1920s
- Large scatter in CT vs. a — no single clean curve fits
- Digitized in Burton et al. (2001) Fig 3.16; that is the version Buhl used

---

## References in the paper

- Lock et al. (1926) — original experimental data
- Glauert (1926) — ARC R&M No. 1026 — original empirical fit
- Burton, Sharpe, Jenkins, Bossanyi (2001) *Wind Energy Handbook*, Wiley, pp. 66–68
- Eggleston & Stoddard (1987) *Wind Turbine Engineering Design*, pp. 30–35, 58
- Manwell, McGowan, Rogers (2002) *Wind Energy Explained*, pp. 120–121
- Wilson (1994) in Spera ed. *Wind Turbine Technology*, pp. 231–232

---

## Relevance to this project

The Buhl correction is the standard WBS fix in BEM codes (including AeroDyn/OpenFAST).
It is already implemented in this codebase's BEM solver.

**Key limitation for flying-turbine validation:** There is no modern, clean experimental
dataset covering deep windmill brake state (a > 0.5, λ₂ > 2). The correction is backed
only by 1920s data with large scatter. Any BEM validation in deep WBS must acknowledge
this uncertainty in the reference itself.
