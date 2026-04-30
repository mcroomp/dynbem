# Table II — Summary of Data on 6-Foot-Diameter Rotor Having Untwisted Blades with 3/1 Taper
**Source: page_33.png (paper p.33), NACA TN-2474**
Columns: V/ΩR | θ₀.₇₅R (deg) | ΔCQ | λ₂ | λ₁(thrust) | λ₁(torque)

Confidence: LOW — table is very dense (small font, multiple run blocks on one page).
Values below are best-effort reads; verify all against page_33.png before use.

ΔCQ = CQ_flight − CQ_hover (negative = rotor in autorotation / energy-harvesting mode).
λ₂ = V_c / V_h (positive = descent in paper convention → NED v_climb < 0).
λ₁ = v_i / V_h (induced velocity normalised by hover induced velocity).

---

## Run blocks identified on page_33

Multiple run blocks are present, grouped by CT and RPM.  Run numbers visible:
Run 60, Run 61, Run 70, Run 71 (and others) at 1200 rpm and 1600 rpm.

### Run 60 — CT/σ ≈ 0.04, 1200 rpm (approximate reads)

| V/ΩR  | θ₀.₇₅R | ΔCQ       | λ₂   | λ₁(thr) | λ₁(trq) |
|-------|--------|-----------|------|---------|---------|
| 0     | ~4.9   | 0         | 0    | ~1.07   | ~1.09   |
| ~.027 | ~5.1   | ~+.000001 | ~.17 | ~1.06   | ~1.09   |
| ~.050 | ~5.3   | ~+.000003 | ~.31 | ~1.07   | ~1.10   |
| ~.072 | ~5.6   | ~+.000007 | ~.47 | ~1.11   | ~1.14   |
| ~.088 | ~5.9   | ~−.000003 | ~.61 | ~1.22   | ~1.25   |
| ~.101 | ~6.2   | ~−.000016 | ~1.03| ~1.82   | ~——     |
| ~.109 | ~6.5   | ~−.000024 | ~1.28| ~2.10   | ~——     |

### Run 61 — CT/σ ≈ 0.08, 1200 rpm (approximate reads)

| V/ΩR  | θ₀.₇₅R | ΔCQ       | λ₂   | λ₁(thr) | λ₁(trq) |
|-------|--------|-----------|------|---------|---------|
| 0     | ~7.5   | 0         | 0    | ~1.05   | ~1.08   |
| ~.030 | ~7.7   | ~+.000005 | ~.14 | ~1.05   | ~1.09   |
| ~.057 | ~8.0   | ~+.000013 | ~.27 | ~1.07   | ~1.11   |
| ~.079 | ~8.4   | ~+.000025 | ~.42 | ~1.10   | ~1.14   |
| ~.090 | ~8.8   | ~+.000038 | ~.56 | ~1.22   | ~1.30   |
| ~.098 | ~9.1   | ~+.000003 | ~.77 | ~1.55   | ~1.68   |
| ~.104 | ~9.6   | ~−.000023 | ~1.07| ~2.18   | ~——     |
| ~.109 | ~10.0  | ~−.000041 | ~1.36| ~2.55   | ~——     |
| ~.113 | ~10.4  | ~−.000051 | ~1.74| ~2.60   | ~——     |
| ~.115 | ~——    | ~−.000041 | ~2.01| ~1.60   | ~——     |

### Run 70 — CT/σ ≈ 0.04, 1600 rpm (approximate reads)

| V/ΩR  | θ₀.₇₅R | ΔCQ       | λ₂   | λ₁(thr) | λ₁(trq) |
|-------|--------|-----------|------|---------|---------|
| 0     | ~4.7   | 0         | 0    | ~1.07   | ~1.10   |
| ~.035 | ~4.9   | ~+.000003 | ~.22 | ~1.07   | ~1.10   |
| ~.067 | ~5.1   | ~+.000008 | ~.43 | ~1.10   | ~1.14   |
| ~.086 | ~5.5   | ~+.000015 | ~.59 | ~1.21   | ~1.28   |
| ~.098 | ~5.8   | ~−.000009 | ~.89 | ~1.70   | ~——     |
| ~.107 | ~6.2   | ~−.000024 | ~1.31| ~2.22   | ~——     |
| ~.112 | ~6.5   | ~−.000031 | ~1.73| ~2.15   | ~——     |
| ~.113 | ~6.7   | ~−.000029 | ~2.00| ~1.47   | ~——     |

### Run 71 — CT/σ ≈ 0.08, 1600 rpm (approximate reads)

| V/ΩR  | θ₀.₇₅R | ΔCQ       | λ₂   | λ₁(thr) | λ₁(trq) |
|-------|--------|-----------|------|---------|---------|
| 0     | ~7.2   | 0         | 0    | ~1.05   | ~1.08   |
| ~.040 | ~7.5   | ~+.000012 | ~.19 | ~1.06   | ~1.09   |
| ~.073 | ~7.9   | ~+.000033 | ~.38 | ~1.09   | ~1.13   |
| ~.090 | ~8.3   | ~+.000053 | ~.55 | ~1.20   | ~1.28   |
| ~.100 | ~8.7   | ~+.000037 | ~.74 | ~1.52   | ~1.64   |
| ~.107 | ~9.2   | ~−.000010 | ~1.05| ~2.13   | ~——     |
| ~.111 | ~9.7   | ~−.000040 | ~1.41| ~2.60   | ~——     |
| ~.113 | ~10.0  | ~−.000048 | ~1.76| ~2.47   | ~——     |
| ~.115 | ~——    | ~−.000037 | ~2.05| ~1.58   | ~——     |

⚠️ All values above are uncertain — read from a dense small-print table.
The combined λ₁/λ₂ data for this configuration is better read from Figure 13 (page_52.png).

---

## Notes

- Taper 3/1 means root chord is 3× tip chord.
- The extended blade root chord c₀ is larger than the physical root chord.
- θ₀.₇₅R for this configuration incorporates the dynamic twist correction:
  θ₀.₇₅R = 0.936(θ_root − 7.79) at 1600 rpm
  θ₀.₇₅R = 0.964(θ_root − 7.79) at 1200 rpm
  (see paper page_11.png / page_12.png)
