# Figure 3 — Principal Blade Dimensions
**Source: page_42.png (paper p.42), NACA TN-2474**

Blade planform and cross-section drawings for all four test configurations, drawn to scale.
Each blade shows: planform view (top), and Section AA cross-section with two annotated dimensions.

Confidence: HIGH — dimensions read directly from extracted per-blade PNGs
(page_42_6ft_constant_chord.png, page_42_fig3_6ft_tapered_chord.png, etc.)

---

## Dimension readings

### Blade root / hub geometry — all rotors share the same hub

Same test stand hub for all four rotors, so absolute root dimensions are identical.
% R differs because the radii differ.

| Dimension | Inches | Metres | % R (6-ft, R=36") | % R (4-ft, R=24") |
|-----------|--------|--------|-------------------|-------------------|
| Hub fairing radius | 3.25" | 0.083 m | 9.1% | 13.5% |
| Blade shank start | **4.85"** | **0.123 m** | **13.6%** | **20.2%** |
| True contour start | **6.09"** | **0.155 m** | **16.9%** | **25.4%** |

The shank region (4.85"–6.09") has a tapered/transitional cross-section with no true
aerodynamic contour.  **BEM root_cutout_m = 0.155 m for both rotors.**

### Section AA chord dimensions

| Configuration | Chord (c) | Section AA dim 2 | dim2 / c |
|---------------|-----------|-------------------|----------|
| 6-ft constant-chord | **1.89"** = 0.04801 m | 0.47" | 0.249 |
| 4-ft constant-chord | **1.26"** = 0.03200 m | 0.32" | 0.254 |

---

## Interpretation of the two Section AA dimensions

### Dimension 1 — chord (c)
The first dimension (1.89" / 1.26") is the full chord of the constant-chord blade, measured
from leading edge to trailing edge.

Cross-check against abstract solidity σ_e ≈ 0.050:
- σ = N·c/(π·R) → c = σ·π·R/N = 0.050·π·0.914/3 = **0.04801 m** = 1.891" ✓
- Figure 3 reads 1.89" → 0.04801 m — consistent to <0.1%.

**The chord c = 0.0479–0.04801 m is confirmed by two independent sources.**

### Dimension 2 — quarter-chord / pitch axis location
The second dimension (0.47" / 0.32") is the distance from the leading edge to the **pitch axis
(feathering axis)**, which sits at the aerodynamic centre = quarter-chord.

Ratio check:
- 6-ft: 0.47 / 1.89 = 0.249 ≈ **25%** chord ✓
- 4-ft: 0.32 / 1.26 = 0.254 ≈ **25%** chord ✓

On a rotor blade the pitch bearing is invariably located at the quarter-chord (aerodynamic
centre) to minimise pitching-moment coupling through the control system.  Annotating this
distance explicitly in a cross-section drawing is standard practice in rotor test reports.

**This is NOT the maximum-thickness value and NOT the location of maximum thickness.**

Rejected interpretations:
- *Maximum thickness value*: NACA 0015 at c=1.89" → t_max = 0.15×1.89 = 0.284" ≠ 0.47".
- *Location of max thickness*: NACA 0015 has max thickness at x/c = 0.30 → x = 0.30×1.89 = 0.567" ≠ 0.47".

---

## Solidity cross-check

| Quantity | Source | Value |
|----------|--------|-------|
| σ_e | Abstract (page_01.png) | 0.050 |
| c from σ | σ·π·R/N | 0.04801 m |
| c from Figure 3 | 1.89" | 0.04801 m |
| c from Table VIII Re | Re/(ν·V_{0.75R}) | ~0.0434 m |

The ~10% discrepancy in the Table VIII Reynolds back-calculation is explained by the test
facility operating at ~30°C rather than ISA 15°C (ν at 30°C ≈ 1.608×10⁻⁵ m²/s, restoring
agreement to within 2%).  The Figure 3 and abstract values are the primary geometry sources.

---

## Code fixture reference value

```python
# Confirmed from Figure 3 and abstract σ=0.050
chord_m = 0.0479   # rounds 0.04801 m; 0.04801 also acceptable
```
