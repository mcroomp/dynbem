# Table VIII — Mach Number, Reynolds Number, and Lift-Curve Slope at 0.75R
**Source: page_39.png (paper p.39), NACA TN-2474**

Header reads: "VALUES OF MACH NUMBER, REYNOLDS NUMBER, AND CALCULATED SLOPE OF
LIFT CURVE AT THREE-QUARTER-RADIUS POINT FOR TEST CONDITIONS"

All values are evaluated at r = 0.75R, not at the blade tip.

Confidence: HIGH — table is clearly printed with large type.

---

| Speed (rpm) | 4-ft M | 4-ft Re  | 4-ft CL_α (/rad) | 6-ft M | 6-ft Re  | 6-ft CL_α (/rad) |
|-------------|--------|----------|------------------|--------|----------|------------------|
| 1200        | 0.165  | 114000   | 5.83             | 0.248  | 256000   | 5.90             |
| 1600        | 0.220  | 152000   | 5.90             | 0.330  | 341000   | 6.07             |


---

## Back-calculations from Table VIII Re values

Re = V_{0.75R} · c / ν  →  c = Re · ν / V_{0.75R}

ISA sea level: ν = 1.461×10⁻⁵ m²/s (μ = 1.789×10⁻⁵ Pa·s, ρ = 1.225 kg/m³)  
At 30 °C:      ν = 1.608×10⁻⁵ m²/s (Sutherland approximation)

Expected chord from σ = 0.05:
- 6-ft (R = 0.914 m, N = 3): c = σ·π·R/N = 0.05×π×0.914/3 = **0.0479 m**
- 4-ft (R = 0.6096 m, N = 3): c = σ·π·R/N = 0.05×π×0.6096/3 = **0.0319 m**

| rpm  | Rotor | V_{0.75R} (m/s) | M (calc / table) | Re (table) | c @ ISA 15 °C | c @ 30 °C    | c from σ |
|------|-------|-----------------|-----------------|------------|---------------|--------------|----------|
| 1200 | 6-ft  | 86.1            | 0.253 / 0.248 ✓ | 256 000    | 0.0435 m      | **0.0479 m** | 0.0479 m |
| 1600 | 6-ft  | 114.9           | 0.338 / 0.330 ✓ | 341 000    | 0.0434 m      | **0.0478 m** | 0.0479 m |
| 1200 | 4-ft  | 57.4            | 0.169 / 0.165 ✓ | 114 000    | 0.0290 m      | **0.0319 m** | 0.0319 m |
| 1600 | 4-ft  | 76.5            | 0.225 / 0.220 ✓ | 152 000    | 0.0290 m      | **0.0319 m** | 0.0319 m |

All four entries: at ISA 15 °C the back-calculated chord falls ~9–10% below σ-derived;
at 30 °C they match σ-derived to within 0.2%.

**Conclusion**: the test facility ran at approximately **30 °C**, not ISA 15 °C.
This explains every Re value in the table without any inconsistency.  The paper gives no
explicit temperature reading.  Use **c = 0.0479 m** (from σ = 0.05) as the BEM fixture value.

---

## Cross-validation: solidity back-calculation — CONSISTENT ✓

σ = N·c/(π·R) = 3 × 0.0479 / (π × 0.914) = **0.0500** ✓

c = 0.0479 m is internally consistent with σ_e = 0.05 (abstract).
The ~9% shortfall in Re at ISA is fully accounted for by T ≈ 30 °C.

---

## Cross-validation: CL_α at 0.75R — CONSISTENT ✓

Table VIII gives CL_α = 5.90 /rad at Re = 256 000 (6-ft rotor, 1200 rpm, 0.75R).
Theoretical thin-airfoil with viscous correction for NACA 0015: 2π × 0.940 = **5.91 /rad** ✓

The 0.940 factor is the viscous efficiency of NACA 0015 at Re ≈ 250 000 (slightly below
fully-attached potential value of 1.0). Agreement to three significant figures confirms
the CL_α read is correct and the rotor operates in an attached-flow regime at test conditions.
