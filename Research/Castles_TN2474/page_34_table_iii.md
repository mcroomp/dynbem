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

### Run 65 (CT/σ ≈ 0.04, 1200 rpm)

| V/ΩR   | θ₀.₇₅R | ΔCQ        | λ₂    | λ₁(T) | λ₁(Q) |
|--------|--------|------------|-------|-------|-------|
| 0      | 5.33   | 0.000000   | 0     | 1.07  | 1.03  |
| .0211  | 5.51   | 0.000002   | 0.13  | 1.07  | 1.05  |
| .0440  | 5.80   | 0.000006   | 0.27  | 1.08  | 1.08  |
| .0663  | 6.21   | 0.000013   | 0.43  | 1.11  | 1.14  |
| .0821  | 6.60   | 0.000016   | 0.57  | 1.24  | 1.29  |
| .0923  | 7.01   | 0.000005   | 0.79  | 1.63  | —     |
| .0991  | 7.33   | -0.000012  | 1.10  | 2.22  | —     |
| .1032  | 7.71   | -0.000023  | 1.39  | 2.67  | —     |
| .1070  | 8.00   | -0.000031  | 1.74  | 2.76  | —     |
| .1101  | 7.81   | -0.000028  | 2.11  | 1.70  | —     |

---

## Run 66 (CT/σ ≈ 0.08, 1200 rpm)

| V/ΩR   | θ₀.₇₅R | ΔCQ        | λ₂    | λ₁(T) | λ₁(Q) |
|--------|--------|------------|-------|-------|-------|
| 0      | 8.31   | 0.000000   | 0     | 1.06  | 1.05  |
| .0270  | 8.60   | 0.000009   | 0.13  | 1.06  | 1.06  |
| .0521  | 9.01   | 0.000022   | 0.26  | 1.07  | 1.09  |
| .0740  | 9.50   | 0.000041   | 0.40  | 1.10  | 1.14  |
| .0870  | 10.00  | 0.000059   | 0.55  | 1.22  | 1.31  |
| .0952  | 10.51  | 0.000048   | 0.73  | 1.55  | 1.69  |
| .1001  | 11.00  | 0.000008   | 1.00  | 2.09  | —     |
| .1050  | 11.50  | -0.000038  | 1.32  | 2.75  | —     |
| .1074  | 12.01  | -0.000057  | 1.61  | 3.10  | —     |
| .1092  | 11.80  | -0.000058  | 2.00  | 2.41  | —     |
| .1121  | 11.31  | -0.000043  | 2.36  | 1.60  | —     |

---

## Run 67 (CT/σ ≈ 0.04, 1600 rpm)

| V/ΩR   | θ₀.₇₅R | ΔCQ        | λ₂    | λ₁(T) | λ₁(Q) |
|--------|--------|------------|-------|-------|-------|
| 0      | 5.11   | 0.000000   | 0     | 1.07  | 1.04  |
| .0281  | 5.30   | 0.000003   | 0.17  | 1.07  | 1.06  |
| .0550  | 5.61   | 0.000009   | 0.35  | 1.08  | 1.11  |
| .0772  | 5.93   | 0.000017   | 0.51  | 1.15  | 1.21  |
| .0911  | 6.30   | 0.000012   | 0.71  | 1.43  | —     |
| .1000  | 6.60   | -0.000011  | 1.07  | 2.32  | —     |
| .1061  | 7.01   | -0.000027  | 1.55  | 2.82  | —     |
| .1100  | 6.81   | -0.000026  | 2.09  | 1.65  | —     |

---

## Run 68 (CT/σ ≈ 0.08, 1600 rpm)

| V/ΩR   | θ₀.₇₅R | ΔCQ        | λ₂    | λ₁(T) | λ₁(Q) |
|--------|--------|------------|-------|-------|-------|
| 0      | 7.90   | 0.000000   | 0     | 1.06  | 1.05  |
| .0330  | 8.21   | 0.000014   | 0.15  | 1.06  | 1.07  |
| .0612  | 8.70   | 0.000036   | 0.30  | 1.08  | 1.11  |
| .0831  | 9.20   | 0.000063   | 0.47  | 1.12  | 1.19  |
| .0960  | 9.81   | 0.000073   | 0.65  | 1.33  | 1.46  |
| .1023  | 10.40  | 0.000028   | 0.89  | 1.85  | —     |
| .1061  | 11.01  | -0.000024  | 1.21  | 2.68  | —     |
| .1090  | 11.60  | -0.000070  | 1.58  | 3.15  | —     |
| .1111  | 11.50  | -0.000073  | 1.99  | 2.55  | —     |
| .1130  | 11.00  | -0.000050  | 2.39  | 1.68  | —     |

---

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
