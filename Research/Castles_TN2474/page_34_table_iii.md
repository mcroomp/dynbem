# Table III — Summary of Data on 6-Foot-Diameter Rotor with Constant-Chord, Twisted Blades (12° Linear Washout)
**Source: page_34.png (paper p.34), NACA TN-2474**
Columns: V/ΩR | θ₀.₇₅R (deg) | ΔCQ | λ₂ | λ₁(thrust) | λ₁(torque)

Confidence: LOW — table is very dense (small font, multiple run blocks on one page).
Values below are best-effort reads; verify all against page_34.png before use.

ΔCQ = CQ_flight − CQ_hover (negative = energy-harvesting mode).
λ₂ = V_c / V_h (positive = descent in paper convention → NED v_climb < 0).
λ₁ = v_i / V_h.

Blade twist: 12° linear washout from hub to tip.  θ₀.₇₅R is the pitch at 75% radius
accounting for dynamic twist correction.

---

## Run blocks identified on page_34

Run numbers visible: Run 65, Run 66, Run 67, Run 68 (and others) at 1200 rpm and 1600 rpm.

### Run 65 — CT/σ ≈ 0.04, 1200 rpm (approximate reads)

| V/ΩR  | θ₀.₇₅R | ΔCQ       | λ₂   | λ₁(thr) | λ₁(trq) |
|-------|--------|-----------|------|---------|---------|
| 0     | ~5.3   | 0         | 0    | ~1.07   | ~1.03   |
| ~.021 | ~5.5   | ~+.000002 | ~.13 | ~1.07   | ~1.05   |
| ~.044 | ~5.8   | ~+.000006 | ~.27 | ~1.08   | ~1.08   |
| ~.066 | ~6.2   | ~+.000013 | ~.43 | ~1.11   | ~1.14   |
| ~.082 | ~6.6   | ~+.000016 | ~.57 | ~1.24   | ~1.29   |
| ~.092 | ~7.0   | ~+.000005 | ~.79 | ~1.63   | ~——     |
| ~.099 | ~7.3   | ~−.000012 | ~1.10| ~2.22   | ~——     |
| ~.103 | ~7.7   | ~−.000023 | ~1.39| ~2.67   | ~——     |
| ~.107 | ~8.0   | ~−.000031 | ~1.74| ~2.76   | ~——     |
| ~.110 | ~——    | ~−.000028 | ~2.11| ~1.70   | ~——     |

### Run 66 — CT/σ ≈ 0.08, 1200 rpm (approximate reads)

| V/ΩR  | θ₀.₇₅R | ΔCQ       | λ₂   | λ₁(thr) | λ₁(trq) |
|-------|--------|-----------|------|---------|---------|
| 0     | ~8.3   | 0         | 0    | ~1.06   | ~1.05   |
| ~.027 | ~8.6   | ~+.000009 | ~.13 | ~1.06   | ~1.06   |
| ~.052 | ~9.0   | ~+.000022 | ~.26 | ~1.07   | ~1.09   |
| ~.074 | ~9.5   | ~+.000041 | ~.40 | ~1.10   | ~1.14   |
| ~.087 | ~10.0  | ~+.000059 | ~.55 | ~1.22   | ~1.31   |
| ~.095 | ~10.5  | ~+.000048 | ~.73 | ~1.55   | ~1.69   |
| ~.100 | ~11.0  | ~+.000008 | ~1.00| ~2.09   | ~——     |
| ~.105 | ~11.5  | ~−.000038 | ~1.32| ~2.75   | ~——     |
| ~.107 | ~12.0  | ~−.000057 | ~1.61| ~3.10   | ~——     |
| ~.109 | ~——    | ~−.000058 | ~2.00| ~2.41   | ~——     |
| ~.112 | ~——    | ~−.000043 | ~2.36| ~1.60   | ~——     |

### Run 67 — CT/σ ≈ 0.04, 1600 rpm (approximate reads)

| V/ΩR  | θ₀.₇₅R | ΔCQ       | λ₂   | λ₁(thr) | λ₁(trq) |
|-------|--------|-----------|------|---------|---------|
| 0     | ~5.1   | 0         | 0    | ~1.07   | ~1.04   |
| ~.028 | ~5.3   | ~+.000003 | ~.17 | ~1.07   | ~1.06   |
| ~.055 | ~5.6   | ~+.000009 | ~.35 | ~1.08   | ~1.11   |
| ~.077 | ~5.9   | ~+.000017 | ~.51 | ~1.15   | ~1.21   |
| ~.091 | ~6.3   | ~+.000012 | ~.71 | ~1.43   | ~——     |
| ~.100 | ~6.6   | ~−.000011 | ~1.07| ~2.32   | ~——     |
| ~.106 | ~7.0   | ~−.000027 | ~1.55| ~2.82   | ~——     |
| ~.110 | ~——    | ~−.000026 | ~2.09| ~1.65   | ~——     |

### Run 68 — CT/σ ≈ 0.08, 1600 rpm (approximate reads)

| V/ΩR  | θ₀.₇₅R | ΔCQ       | λ₂   | λ₁(thr) | λ₁(trq) |
|-------|--------|-----------|------|---------|---------|
| 0     | ~7.9   | 0         | 0    | ~1.06   | ~1.05   |
| ~.033 | ~8.2   | ~+.000014 | ~.15 | ~1.06   | ~1.07   |
| ~.061 | ~8.7   | ~+.000036 | ~.30 | ~1.08   | ~1.11   |
| ~.083 | ~9.2   | ~+.000063 | ~.47 | ~1.12   | ~1.19   |
| ~.096 | ~9.8   | ~+.000073 | ~.65 | ~1.33   | ~1.46   |
| ~.102 | ~10.4  | ~+.000028 | ~.89 | ~1.85   | ~——     |
| ~.106 | ~11.0  | ~−.000024 | ~1.21| ~2.68   | ~——     |
| ~.109 | ~11.6  | ~−.000070 | ~1.58| ~3.15   | ~——     |
| ~.111 | ~——    | ~−.000073 | ~1.99| ~2.55   | ~——     |
| ~.113 | ~——    | ~−.000050 | ~2.39| ~1.68   | ~——     |

⚠️ All values above are uncertain — read from a dense small-print table.
The combined λ₁/λ₂ data for this configuration is better read from Figure 14 (page_53.png).

---

## Notes

- VRS peak λ₁ values (~3.1) are noticeably higher than for the constant-chord untwisted
  rotor (~2.6–2.8), consistent with the paper's conclusion that 12° twist increases peak
  induced velocity by ~24%.
- The autorotation crossing (ΔCQ = 0) occurs at higher V/ΩR than for the untwisted rotor,
  consistent with the paper's statement that twist increases autorotation descent rate ~10%.
- A chordwise bending fatigue failure occurred on one of the twisted blades during testing
  at CT ≈ 0.004 at 1600 rpm; hovering check run was not obtained for this configuration.
