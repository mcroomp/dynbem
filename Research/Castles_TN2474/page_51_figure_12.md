# Figure 12 — λ₁ vs λ₂, 6-ft Constant-Chord Untwisted Blades
**Source: page_51.png (paper p.51), NACA TN-2474**

Graph: normalised induced velocity λ₁ = v_i/V_h vs normalised descent rate λ₂ = V_c/V_h.
X-axis runs right-to-left: λ₂ = 0 (hover) on the right, λ₂ = 2.8 on the left.
Two datasets: circles = λ₁ from thrust measurements; triangles = λ₁ from torque measurements.
Two reference curves: "Simple Momentum Theory" hyperbola; "Ideal Autorotation Line" (λ₁ = 1).

Confidence: MODERATE — values read from printed graph.

---

## Momentum theory formula (WBS regime, λ₂ > 2)

λ₁ = λ₂/2 − √(λ₂²/4 − 1)

Analytical values and measured scatter from Figure 12:

| λ₂  | λ₁ theory | λ₁ data centre | data scatter range |
|-----|-----------|---------------|--------------------|
| 2.0 | 1.000     | ~1.00         | 0.85–1.15          |
| 2.5 | 0.500     | ~0.55         | 0.50–0.70          |
| 3.0 | 0.382     | ~0.40         | 0.35–0.48          |
| 4.0 | 0.268     | ~0.28         | 0.25–0.32          |

WBS data (λ₂ > 2.0) tracks momentum theory to within ~20%; scatter shrinks at higher λ₂.

---

## VRS regime (λ₂ ≈ 0.5–2.0) — DO NOT use for BEM validation

The peak of the data cloud occurs around λ₂ ≈ 1.2–1.5, λ₁ ≈ 2.6–3.0.
Momentum theory gives λ₁ ≈ 1.4–1.9 in this range — massive under-prediction.
This is the Vortex Ring State where momentum theory breaks down entirely.
BEM (which is based on momentum theory) is not valid here.

---

## Analytical verification of momentum theory values

λ₁ = λ₂/2 − √(λ₂²/4 − 1):
- λ₂ = 2.0: 1.0 − √(0) = **1.000** ✓
- λ₂ = 2.5: 1.25 − √(0.5625) = 1.25 − 0.750 = **0.500** ✓
- λ₂ = 3.0: 1.5 − √(1.25) = 1.5 − 1.118 = **0.382** ✓
- λ₂ = 4.0: 2.0 − √(3.0) = 2.0 − 1.732 = **0.268** ✓

---

## Cross-validation: WBS momentum theory — VERIFIED ✓

λ₁ = λ₂/2 − √(λ₂²/4 − 1):
- λ₂ = 2.0: 1.0 − √(0)    = **1.000** ✓
- λ₂ = 2.5: 1.25 − √(0.5625) = 1.25 − 0.750 = **0.500** ✓
- λ₂ = 3.0: 1.5 − √(1.25)  = 1.5 − 1.118  = **0.382** ✓
- λ₂ = 4.0: 2.0 − √(3.0)   = 2.0 − 1.732  = **0.268** ✓

Measured data centres from Figure 12 are within ~20% of these values in WBS.

---

## Cross-validation: NED sign convention mapping — VERIFIED ✓

Paper convention: λ₂ > 0 = descent = upward flow through disk (air opposes thrust).
NED convention: v_climb < 0 = upward flow (same physical direction).
Mapping: **v_climb = −λ₂ · V_h**

This is consistent with the CLAUDE.md sign convention: `v_climb < 0` in autorotation /
flying-turbine mode, `lambda` (inflow ratio) negative in energy-harvesting mode.

WBS entry condition: λ₂ > 2.0 → **v_climb < −2·V_h**.
At CT = 0.004, 1000 rpm (V_h ≈ 4.29 m/s): WBS requires v_climb < −8.6 m/s.

---

## NED translation for BEM test setup

Paper convention: λ₂ > 0 = descent (upward flow through disk).
NED convention: v_climb < 0 = upward flow.
Mapping: **v_climb = −λ₂ · V_h**

For CT = 0.004 at 1200 rpm (V_h ≈ 4.29 m/s):

| λ₂  | NED v_climb (m/s) | V/ΩR (at 1200 rpm) |
|-----|------------------|-------------------|
| 2.0 | −8.58            | −0.090            |
| 2.5 | −10.73           | −0.112            |
| 3.0 | −12.87           | −0.134            |
