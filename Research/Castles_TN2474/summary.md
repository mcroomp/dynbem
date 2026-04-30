# NACA TN-2474 Summary
## "Empirical Study of the Induced-Velocity Distribution Function for a Model Helicopter Rotor in Vertical Flight Including Autorotation"
### Castles, W. Jr. & Gray, R.B. — NACA, September 1951

**File naming**: page_NN.png = paper page NN.  page_cover1.png / page_cover2.png = the two
unnumbered cover pages at the front of the PDF.

---

## 1. Purpose and scope

Provides an empirical curve of normalised induced velocity (λ₁ = v_i/V_h) against
normalised climb/descent rate (λ₂ = V_c/V_h) spanning hover through autorotation for
three model-rotor blade configurations at two diameters.  The data cover hover (λ₂ = 0),
normal climb, Vortex Ring State (VRS, λ₂ ≈ 0.5–2.0), and Windmill Brake State (WBS, λ₂ > 2.0).

**Motivation**: Prior wind-tunnel data (Lock et al. 1926, ref. 7) used to construct
Glauert's empirical 1/f vs 1/F curve were in serious disagreement with flight-test data
(Stewart 1948, ref. 2).  Castles & Gray identified three sources of error in the earlier
tests and corrected all three:

1. Closed test section with energy ratio >> 1 (present tests: open jet, ER ≈ 0.7 via screen)
2. Equivalent free-stream velocity taken from a static-pressure wall tap (present: tunnel
   propeller speed calibration)
3. No correction for dynamic blade twist (present: dynamic twist calculated and subtracted)

Results are presented in both λ₁ vs λ₂ and 1/f vs 1/F forms, with flow-visualisation
photographs at six λ₂ values from hover through WBS.

**Coordinate note for this project**: the paper defines λ₂ as positive in *descent*
(V_c > 0 = air flowing *upward* through disk, opposing rotor thrust).  In this project's
NED convention that corresponds to v_climb < 0.  Level-1 BEM is only valid in
WBS (λ₂ > 2.0, i.e. v_climb < −2·V_h).

---

## 2. Rotor configurations tested

| Label | Diameter | N blades | Chord (derived) | Taper | Twist | σ (nominal) |
|-------|----------|----------|-----------------|-------|-------|-------------|
| 6-ft constant-chord | 1.829 m | 3 | 0.0479 m | none | 0°  | 0.050 |
| 6-ft 3/1 tapered    | 1.829 m | 3 | varies   | 3:1  | 0°  | ~0.05 |
| 6-ft 12° twist      | 1.829 m | 3 | 0.0479 m | none | 12° | ~0.05 |
| 4-ft constant-chord | 1.219 m | 3 | 0.0319 m | none | 0°  | 0.050 |

All blades: NACA 0015.  Effective solidity σ_e = 0.05 (stated in abstract, paper §"Model rotor blades").
Chord back-calculation: c = σ·π·R/N
- 6-ft: 0.050·π·0.914/3 = **0.0479 m** ✓
- 4-ft: 0.050·π·0.6096/3 = **0.0319 m** ✓

**Blade construction** (paper §"Model rotor blades"): solid alloy-steel leading edge back
to ~quarter-chord; constant-chord blades had hollow magnesium trailing-edge section;
twisted and tapered blades had solid laminated mahogany trailing edge.  Blade airfoil/CG/EA
all coincident at quarter-chord (pitch-change axis).

**Hub and root geometry** (paper §"Model rotor blades" + Figure 3 direct read):
All four rotors use the same test-stand hub, so absolute root dimensions are identical:
- Hub fairing radius: 3.25" = 0.083 m
- Blade shank start: 4.85" = 0.123 m (tapered, no airfoil contour)
- True aerodynamic contour starts: **6.09"** = **0.155 m** — paper text; Figure 3 confirms
- % R: 16.9% for 6-ft rotor; 25.4% for 4-ft rotor
- **BEM root_cutout_m = 0.155 m for both rotors** — see Section 7

Blade geometry drawing: Figure 3 (page_42.png) → [page_42_figure_3.md](page_42_figure_3.md)

---

## 3. Normalisation convention

- V_h = √(T / (2ρA))  — hover induced velocity [m/s]
- λ₁ = V_i / V_h  (induced velocity, ≥ 0)
- λ₂ = V / V_h    (**positive = descent** in this paper, V = descent velocity)
- ΔC_Q = C_Q − C_Q(zero thrust, zero descent) — torque increment above profile-drag baseline

Conversion to Glauert's coordinates (paper eqs. 12–13):
- 1/f = λ₂²
- 1/F = (λ₁ − λ₂)²

Simple momentum theory (paper eq. 15):

| Regime | λ₂ | Formula |
|--------|----|---------|
| Climb (NED v_climb > 0) | < 0 | λ₁ = λ₂/2 + √(λ₂²/4 + 1) |
| Hover | 0 | λ₁ = 1 |
| WBS | > 2 | λ₁ = λ₂/2 − √(λ₂²/4 − 1) |

"Ideal autorotation line": λ₁ = 1 (constant; energy balance V_i·T = 0 net power).

---

## 4. Key findings (from abstract and concluding remarks)

**General observation 1**: Present λ₁ values are *lower* than Glauert's curve at hover and
small descent rates, and *higher* at larger descent rates.  Good agreement with full-scale
flight tests at the hover and autorotation ends; peak VRS values exceed full-scale
(attributed to inability of full-scale aircraft to maintain steady VRS).

**General observation 2 — 3/1 taper**: Slightly reduces λ₁ at hover/small descent.
Increases ideal autorotation rate of descent by ~3%.

**General observation 3 — 12° linear twist**: Slightly reduces λ₁ at hover.  Increases
ideal autorotation descent rate by ~10%.  Raises peak λ₁ by ~24% and shifts the peak to
~17% higher λ₂.  Significantly larger force fluctuations in the VRS region.

**General observation 4**: No significant differences in λ₁ vs λ₂ curves due to variations
in C_T, rotor speed, or rotor diameter — within test range.  *This justifies using one
rotor and one CT value for Level-1 BEM validation.*

**No fluctuations in autorotation**: Forces and moments were steady for all rotors in the
autorotation range (λ₂ > ~2).

**Dynamic twist corrections** (paper §"Reduction of Data"):
- 6-ft constant-chord at 1600 rpm: θ₀.₇₅R = 0.963 × θ_root
- 6-ft constant-chord at 1200 rpm: θ₀.₇₅R = 0.978 × θ_root
- 4-ft constant-chord at 1600 rpm: θ₀.₇₅R = 0.960 × θ_root
- 4-ft constant-chord at 1200 rpm: θ₀.₇₅R = 0.965 × θ_root
- (Tapered and twisted blades have different corrections — see paper text)

Table VIII lift-curve slopes (0.75R, two-dimensional corrected for Re and Mach):
- 6-ft at 1200 rpm: a = 5.90 /rad (Re = 256 000, M = 0.248)
- 6-ft at 1600 rpm: a = 6.07 /rad (Re = 341 000, M = 0.330)

---

## 4b. Page-by-page contents

| File | Content |
|------|---------|
| page_cover1.png | Title page (NACA TN-2474) |
| page_cover2.png | Verso / blank |
| page_01.png | Abstract: 4 rotors, σ_e ≈ 0.05, NACA 0015, empirical λ₁ vs λ₂ relation |
| page_02–12.png | Introduction, symbols list, test apparatus, procedure |
| page_13–29.png | Body text: BEM background, Glauert's curve, test methodology, results discussion |
| page_30–31.png | Appendix equations (A22–A27), radial-distribution derivation |
| page_31.png | References (7 entries: Glauert 1926, Stewart 1948, Gessow 1945, …) |
| **page_32.png** | **Table I** — 6-ft constant-chord untwisted: V/ΩR, θ₀.₇₅R, ΔCQ, λ₂, λ₁(thrust), λ₁(torque) |
| **page_33.png** | **Table II** — 6-ft 3/1 tapered blades: same columns |
| **page_34.png** | **Table III** — 6-ft 12° linear twist: same columns |
| **page_35.png** | **Table IV** — 4-ft constant-chord untwisted: same columns |
| **page_36.png** | **Table V** — Hover data, 6-ft constant-chord untwisted: θ₀.₇₅R, CT, ΔCQ |
| page_37.png | Table VI — Hover data, 6-ft 3/1 tapered |
| page_38.png | Table VII — Hover data, 4-ft constant-chord untwisted |
| page_39.png | **Table VIII** — Mach, Reynolds, CL_α at test conditions |
| page_40.png | Figure 1 — Photograph: test stand with 6-ft rotor installed |
| page_41.png | Figure 2 — Schematic cross-section of test stand |
| **page_42.png** | **Figure 3** — Principal blade dimensions (all four blade types to scale); see [page_42_figure_3.md](page_42_figure_3.md) |
| page_43.png | Figure 4 — θ₀.₇₅R vs V/ΩR, 6-ft constant-chord (CT/σ = 0.04 / 0.08 / 0.10) |
| page_44.png | Figure 5 — θ₀.₇₅R vs V/ΩR, 6-ft 3/1 tapered |
| page_45.png | Figure 6 — θ₀.₇₅R vs V/ΩR, 6-ft 12° twist |
| page_46.png | Figure 7 — θ₀.₇₅R vs V/ΩR, 4-ft constant-chord |
| **page_47.png** | **Figure 8** — ΔCQ vs V/ΩR, 6-ft constant-chord (CT/σ = 0.04 / 0.08 / 0.10); autorotation crossings visible |
| page_48.png | Figure 9 — ΔCQ vs V/ΩR, 6-ft 3/1 tapered |
| page_49.png | Figure 10 — ΔCQ vs V/ΩR, 6-ft 12° twist (hump at VRS clearly visible) |
| page_50.png | Figure 11 — ΔCQ vs V/ΩR, 4-ft constant-chord |
| **page_51.png** | **Figure 12** — KEY: λ₁ vs λ₂, 6-ft constant-chord; momentum theory hyperbola, ideal autorotation line, measured scatter |
| page_52.png | Figure 13 — λ₁ vs λ₂, 6-ft 3/1 tapered |
| page_53.png | Figure 14 — λ₁ vs λ₂, 6-ft 12° twist |
| page_54.png | Figure 15 — λ₁ vs λ₂, 4-ft constant-chord |
| **page_55.png** | **Figure 16** — All data overlaid: Castles-Gray data, Glauert empirical curve (Ref 1), full-scale thrust (Ref 2), autorotation data (Ref 3) |
| page_56.png | Figure 17 — Same data as 1/f vs 1/F coordinates |
| page_57.png | Figure 18 — Flow visualisation (smoke/tufts) at λ₂ = 0 (hover) |
| page_58.png | Figure 19 — Flow at λ₂ ≈ 0.3: wake begins to curl back |
| page_59.png | Figure 20 — Flow at λ₂ ≈ 1.0: vortex ring state; large recirculation bubble |
| page_60.png | Figure 21 — Flow at λ₂ ≈ 1.35: VRS maximum; chaotic vortex core |
| page_61.png | Figure 22 — Flow at λ₂ ≈ 1.7: turbulence subsiding |
| page_62.png | Figure 23 — Flow at λ₂ ≈ 2.0: flow stabilising into WBS |
| page_63.png | Figure 24 — Photograph of tufts at λ₂ = 0.3 |
| **page_64.png** | **Figure 25** — CT vs θ₀.₇₅R (hover comparison), 6-ft constant-chord; data + theory line |
| page_65.png | Figure 26 — CT vs θ₀.₇₅R (hover), 6-ft 3/1 tapered |
| page_66.png | Figure 27 — CT vs θ₀.₇₅R (hover), 4-ft constant-chord |
| page_67.png | Figure 28 — ΔCQ vs CT (hover performance polar), 6-ft constant-chord |
| page_68.png | Figure 29 — ΔCQ vs CT (hover), 6-ft 3/1 tapered |
| page_69.png | Figure 30 — ΔCQ vs CT (hover), 4-ft constant-chord |
| page_70.png | Figure 31 — Actuator-disk flow pattern schematic |
| page_71.png | Figure 32 — Vortex-ring-state transition schematic |
| page_72.png | Figure 33 — Hypothetical wind-tunnel schematic |

---

## 5. Digitised data tables

Full table files (with confidence annotations and cross-checks):
- [page_36_table_v.md](page_36_table_v.md) — Hover CT, 6-ft constant-chord (HIGH confidence)
- [page_37_table_vi.md](page_37_table_vi.md) — Hover CT, 6-ft 3/1 tapered (MODERATE)
- [page_38_table_vii.md](page_38_table_vii.md) — Hover CT, 4-ft constant-chord (MODERATE)
- [page_39_table_viii.md](page_39_table_viii.md) — Mach/Re/CL_α at 0.75R (HIGH)
- [page_32_table_i.md](page_32_table_i.md) — Descent data, 6-ft constant-chord (LOW)

### Table V — Hover CT, 6-ft constant-chord untwisted (HIGH)
→ [page_36_table_v.md](page_36_table_v.md)

### Table VIII — Mach/Re/CL_α at 0.75R (HIGH)
→ [page_39_table_viii.md](page_39_table_viii.md)

### Table I — Descent data, 6-ft constant-chord (LOW — dense table)
→ [page_32_table_i.md](page_32_table_i.md) — WBS λ₁/λ₂ better read from Figure 12.

### Table II — Descent data, 6-ft 3/1 tapered (LOW — dense table)
→ [page_33_table_ii.md](page_33_table_ii.md) — λ₁/λ₂ better read from Figure 13 (page_52.png).

### Table III — Descent data, 6-ft 12° twist (LOW — dense table)
→ [page_34_table_iii.md](page_34_table_iii.md) — λ₁/λ₂ better read from Figure 14 (page_53.png).

### Table IV — Descent data, 4-ft constant-chord (LOW — dense table)
→ [page_35_table_iv.md](page_35_table_iv.md) — λ₁/λ₂ better read from Figure 15 (page_54.png).

### Figure 8 — ΔCQ vs V/ΩR autorotation crossings (MODERATE — graph read)
→ [page_47_figure_8.md](page_47_figure_8.md)

### Figure 12 — λ₁ vs λ₂, WBS regime (MODERATE — graph read)
→ [page_51_figure_12.md](page_51_figure_12.md)

### Figure 3 — Principal blade dimensions (HIGH — direct read from extracted PNGs)
→ [page_42_figure_3.md](page_42_figure_3.md) — chord confirmed 1.89" = 0.04801 m; second Section AA dimension = quarter-chord / pitch axis location (25% chord)

---

## 6. Corrections and cross-validations

Cross-validations and corrections are recorded in each data file:

- [page_32_table_i.md](page_32_table_i.md) — table numbering correction (Table I = 6-ft, not 4-ft)
- [page_36_table_v.md](page_36_table_v.md) — Table V vs Figure 25 cross-check (CONSISTENT ✓)
- [page_39_table_viii.md](page_39_table_viii.md) — solidity back-calc and CL_α check (CONSISTENT ✓)
- [page_47_figure_8.md](page_47_figure_8.md) — autorotation crossing vs V_h calculation (CONSISTENT ✓)
- [page_51_figure_12.md](page_51_figure_12.md) — WBS momentum theory verification and NED mapping (VERIFIED ✓)

### Test facility temperature — ~30 °C (inferred)

The paper gives no explicit air temperature.  Tests were conducted at **Georgia Institute of
Technology, Atlanta, GA** (paper closing, May 1950).  Back-calculating chord from Table VIII
Re values using ISA 15 °C (ν = 1.461×10⁻⁵ m²/s) gives c = 0.0434 m — 9% below the
σ-derived c = 0.0479 m.  At 30 °C (ν = 1.608×10⁻⁵ m²/s) the back-calculated chord matches
σ-derived to within 0.2% across all four table entries (both rotors, both speeds).

**Practical impact on BEM validation**: none.  CT and CQ are density-normalised; ρ cancels
in the coefficient comparison.  The finding confirms c = 0.0479 m is correct and that
Re = 256 000 is the true test condition (not a misprint).

---

## 7. BEM validation implications

### Valid test cases for Level-1 BEM

1. **Hover CT vs collective** (Table V, page_36.png): CT ≈ 0.004 at θ ≈ 8.5°.
   BEM expected to over-predict by ~30–45% (same inviscid bias as Caradonna-Tung).
2. **WBS torque sign**: λ₂ > 2 → Q_total < 0 (rotor harvests power).
3. **Autorotation crossing**: Q_total = 0 at V/ΩR ≈ 0.08 (CT/σ ≈ 0.08).
   BEM should reach Q = 0 within ±25% of this V/ΩR value.

### Invalid regime for Level-1 BEM

VRS (0.5 < λ₂ < 2.0): momentum theory fails; BEM results are non-physical. Do not compare.

### Root cutout discrepancy — affects CT and CQ levels

Paper text: aerodynamic blade starts at r = 6.09 inches = **0.155 m** (16.9% R).
Current `rotor.yaml` uses root_cutout_m = 0.10 m (~10.9% R) estimated from Figure 3.

The larger root cutout reduces disk area and blade planform area, lowering computed CT and CQ.
Impact at first order: disk area A = π(R² − r_root²).
- With 0.10 m: A = π(0.914² − 0.10²) = 2.594 m²
- With 0.155 m: A = π(0.914² − 0.155²) = 2.549 m²  (−1.7%)

This is a small correction to CT/CQ levels but should be fixed for accuracy.
**Recommended**: update `rotor.yaml` root_cutout_m to 0.155 (from paper text), retire the
figure-3 estimate.

---

## 8. Rotor fixture for WBS tests

```python
RotorDefinition(
    blade=BladeGeometry(
        n_blades=3,
        radius_m=0.914,           # R = 3 ft
        root_cutout_m=0.155,      # 6.09" from axis per paper §"Model rotor blades" (16.9% R)
        chord_m=0.0479,           # σ = N·c/(π·R) = 0.050 ✓
        twist_deg=0.0,
        n_elements=30,
    ),
    airfoil=AirfoilProperties(
        Re_design=256_000,
        CL0=0.0,
        CL_alpha_per_rad=5.90,    # Table VIII (page_39.png), 6-ft at 1200 rpm, 0.75R
        CD0=0.012,                # NACA 0015 at Re ≈ 256k
        alpha_stall_deg=12.0,
        tip_loss=True,
    ),
    autorotation=AutorotationProperties(I_ode_kgm2=1.0),
    name="Castles-Gray-6ft",
)
# Test RPM: 1000 → Ω = 104.7 rad/s, ΩR = 95.7 m/s, M_tip = 0.281
# CT/σ = 0.08 target: θ ≈ 8.5° (Table V Run 15 interpolation)
# WBS threshold: V_c > 2·V_h ≈ 8.6 m/s → V/ΩR > 0.090
```
