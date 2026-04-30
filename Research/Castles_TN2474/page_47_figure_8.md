# Figure 8 — Variation of Torque Coefficient, 6-ft Constant-Chord Untwisted Blades
**Source: page_47.png (paper p.47), NACA TN-2474**

Graph: ΔCQ vs V/ΩR for CT/σ = 0.04 (circles), 0.08 (squares), 0.10 (triangles).
ΔCQ = CQ_flight − CQ_hover.  ΔCQ < 0 = rotor in autorotation / energy-harvesting mode.

Confidence: MODERATE — values read from printed graph; ±0.005 in V/ΩR axis.

---

## Autorotation crossings (ΔCQ = 0)

| CT/σ | CT (approx) | V/ΩR where ΔCQ = 0 |
|------|-------------|---------------------|
| 0.04 | 0.0020      | ~0.070              |
| 0.08 | 0.0040      | ~0.083              |
| 0.10 | 0.0050      | ~0.090              |

Cross-check: at CT = 0.004, 1200 rpm (ΩR = 95.7 m/s):
- V_h = √(T/(2ρA)) = √(118/(2×1.225×2.625)) ≈ 4.29 m/s
- λ₂ = 2 corresponds to V_c = 8.58 m/s → V/ΩR = 0.090
- Autorotation should occur at λ₂ slightly < 2 (VRS/WBS boundary), consistent with 0.083 ✓

## Hover ΔCQ values (V/ΩR = 0), read from graph

| CT/σ | CT (approx) | ΔCQ at hover |
|------|-------------|-------------|
| 0.04 | 0.0020      | ~0.000008   |
| 0.08 | 0.0040      | ~0.000026   |
| 0.10 | 0.0050      | ~0.000040   |

These are the absolute CQ values in hover (ΔCQ = 0 at V=0 would be the hover baseline,
so these are the differences accumulated over the range — actually ΔCQ at V/ΩR=0 is the
baseline CQ itself for the constant-CT curves).

## Cross-validation: autorotation crossing vs V_h — CONSISTENT ✓

Independent check for CT/σ = 0.08 (CT = 0.004) at 1200 rpm (ΩR = 95.7 m/s):

```
T   = CT × ρ × A × (ΩR)² = 0.004 × 1.225 × π × 0.914² × 95.7² ≈ 118 N
V_h = √(T / (2ρA)) = √(118 / (2 × 1.225 × π × 0.914²)) ≈ 4.29 m/s
λ₂ = 2 → V_c = 2 × 4.29 = 8.58 m/s → V/ΩR = 8.58 / 95.7 = 0.090
```

Figure 8 reads autorotation crossing at V/ΩR ≈ 0.083 for CT/σ = 0.08.
8% difference — well within graph-reading precision for a printed figure.
Autorotation occurs near the VRS/WBS boundary (λ₂ ≈ 2), consistent with 0.083 < 0.090.

---

## Maximum ΔCQ (deep descent, V/ΩR → 0.10+)

| CT/σ | Max negative ΔCQ (approx) |
|------|--------------------------|
| 0.04 | ~−0.000008               |
| 0.08 | ~−0.000006               |
| 0.10 | ~−0.000004               |

Note: The negative ΔCQ saturates as V/ΩR → 0.10–0.12; data points are sparse in WBS.
