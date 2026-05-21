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

## Run 60 (CT/σ ≈ 0.04, 1200 rpm)

| V/ΩR   | θ₀.₇₅R | ΔCQ        | λ₂    | λ₁(T) | λ₁(Q) |
|--------|--------|------------|-------|-------|-------|
| 0      | 4.91   | 0.000000   | 0     | 1.07  | 1.09  |
| .0272  | 5.10   | 0.000001   | 0.17  | 1.06  | 1.09  |
| .0503  | 5.32   | 0.000003   | 0.31  | 1.07  | 1.10  |
| .0720  | 5.61   | 0.000007   | 0.47  | 1.11  | 1.14  |
| .0881  | 5.94   | -0.000003  | 0.61  | 1.22  | 1.25  |
| .1010  | 6.21   | -0.000016  | 1.03  | 1.82  | —     |
| .1091  | 6.50   | -0.000024  | 1.28  | 2.10  | —     |

---

## Run 61 (CT/σ ≈ 0.08, 1200 rpm)

| V/ΩR   | θ₀.₇₅R | ΔCQ        | λ₂    | λ₁(T) | λ₁(Q) |
|--------|--------|------------|-------|-------|-------|
| 0      | 7.53   | 0.000000   | 0     | 1.05  | 1.08  |
| .0301  | 7.72   | 0.000005   | 0.14  | 1.05  | 1.09  |
| .0570  | 8.03   | 0.000013   | 0.27  | 1.07  | 1.11  |
| .0791  | 8.41   | 0.000025   | 0.42  | 1.10  | 1.14  |
| .0903  | 8.80   | 0.000038   | 0.56  | 1.22  | 1.30  |
| .0981  | 9.13   | 0.000003   | 0.77  | 1.55  | 1.68  |
| .1042  | 9.60   | -0.000023  | 1.07  | 2.18  | —     |
| .1090  | 10.01  | -0.000041  | 1.36  | 2.55  | —     |
| .1133  | 10.40  | -0.000051  | 1.74  | 2.60  | —     |
| .1150  | 10.11  | -0.000041  | 2.01  | 1.60  | —     |

---

## Run 70 (CT/σ ≈ 0.04, 1600 rpm)

| V/ΩR   | θ₀.₇₅R | ΔCQ        | λ₂    | λ₁(T) | λ₁(Q) |
|--------|--------|------------|-------|-------|-------|
| 0      | 4.73   | 0.000000   | 0     | 1.07  | 1.10  |
| .0350  | 4.92   | 0.000003   | 0.22  | 1.07  | 1.10  |
| .0672  | 5.13   | 0.000008   | 0.43  | 1.10  | 1.14  |
| .0860  | 5.50   | 0.000015   | 0.59  | 1.21  | 1.28  |
| .0983  | 5.81   | -0.000009  | 0.89  | 1.70  | —     |
| .1071  | 6.20   | -0.000024  | 1.31  | 2.22  | —     |
| .1120  | 6.51   | -0.000031  | 1.73  | 2.15  | —     |
| .1134  | 6.70   | -0.000029  | 2.00  | 1.47  | —     |

---

## Run 71 (CT/σ ≈ 0.08, 1600 rpm)

| V/ΩR   | θ₀.₇₅R | ΔCQ        | λ₂    | λ₁(T) | λ₁(Q) |
|--------|--------|------------|-------|-------|-------|
| 0      | 7.20   | 0.000000   | 0     | 1.05  | 1.08  |
| .0403  | 7.51   | 0.000012   | 0.19  | 1.06  | 1.09  |
| .0731  | 7.90   | 0.000033   | 0.38  | 1.09  | 1.13  |
| .0900  | 8.33   | 0.000053   | 0.55  | 1.20  | 1.28  |
| .1002  | 8.70   | 0.000037   | 0.74  | 1.52  | 1.64  |
| .1070  | 9.21   | -0.000010  | 1.05  | 2.13  | —     |
| .1111  | 9.70   | -0.000040  | 1.41  | 2.60  | —     |
| .1134  | 10.00  | -0.000048  | 1.76  | 2.47  | —     |
| .1150  | 9.81   | -0.000037  | 2.05  | 1.58  | —     |

---

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
