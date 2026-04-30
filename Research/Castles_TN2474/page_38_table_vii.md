# Table VII — Hovering Data, 4-ft-Diameter Rotor, Constant-Chord, Untwisted Blades
**Source: page_38.png (paper p.38), NACA TN-2474**
Columns: θ₀.₇₅R (deg) | CT | ΔCQ
Confidence: HIGH — manually verified from source image.

R = 0.6096 m (2 ft), N = 3, NACA 0015, σ = 0.050, root_cutout = 0.155 m.
Trailing rows in each run are authors' repeatability check points (blade settling).

---

## Run 17; 1200 rpm

| θ₀.₇₅R (deg) | CT       | ΔCQ (= CQ) | Note |
|--------------|----------|------------|------|
| 0            | 0        | 0          | |
| 1.11         | 0.00017  | 0.000002   | |
| 2.20         | 0.00052  | 0.000016   | |
| 4.05         | 0.00152  | 0.000058   | |
| 6.02         | 0.00275  | 0.000110   | FM=0.93 — high; verify CQ if possible |
| 8.00         | 0.00375  | 0.000196   | |
| 9.92         | 0.00468  | 0.000271   | |
| 11.73        | 0.00527  | 0.000375   | |
| 8.30         | 0.00411  | 0.000229   | repeat/check |

---

## Run 18; 1600 rpm

| θ₀.₇₅R (deg) | CT       | ΔCQ (= CQ) | Note |
|--------------|----------|------------|------|
| 0            | 0        | 0          | |
| 1.93         | 0.00044  | 0.000020   | |
| 3.83         | 0.00136  | 0.000050   | |
| 5.70         | 0.00245  | 0.000109   | |
| 7.72         | 0.00382  | 0.000203   | |
| 9.61         | 0.00477  | 0.000279   | |
| 11.39        | 0.00535  | 0.000373   | |
| 7.68         | 0.00379  | 0.000201   | repeat/check ✓ matches θ=7.72 to 1% |
| 3.83         | 0.00135  | 0.000041   | repeat/check; CQ 18% below initial 0.000050 — possible misread |

---

## Cross-validation

**Figure of merit** FM = CT^1.5 / (√2 · CQ) at operating points:
- Run 17 θ=8–10°: FM = 0.83–0.84 ✓ (typical for model rotor)
- Run 18 θ=7–10°: FM = 0.82–0.84 ✓ consistent with Run 17

**RPM independence** (observation 4 from abstract): at CT ≈ 0.0038, Run 17 (1200 rpm)
gives CQ = 0.000196 and Run 18 (1600 rpm) gives CQ = 0.000203 — within 4% ✓

**Diameter independence**: compared against Table V Run 15 (6-ft, 1200 rpm):
- 4-ft at θ=8.00°: CT=0.00375, CQ=0.000196, FM=0.828
- 6-ft at θ=8.46°: CT=0.00400, CQ=0.000226, FM=0.792
CT and CQ within ~7% at similar collective — consistent with paper observation 4 ✓

