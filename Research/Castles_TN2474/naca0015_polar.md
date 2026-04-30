# NACA 0015 Airfoil Polar Data
**Source: airfoiltools.com XFOIL predictions**
**Files: xf-naca0015-il-200000-n5.csv (NCrit=5), xf-naca0015-il-200000.csv (NCrit=9)**

---

## Which polar to use

The Castles-Gray tests were conducted in a closed-return wind tunnel with a turbulence
energy ratio reduced to ~0.7 by an 18×18-mesh screen.  NCrit=5 represents typical
wind-tunnel turbulence intensity; NCrit=9 is free-air.  **Use NCrit=5 for this dataset.**

Re=200,000 is the closest available to the test condition of Re≈256,000 at 0.75R
(1200 rpm, 6-ft rotor, ISA 15°C — see page_39_table_viii.md for details).

---

## Key values extracted from NCrit=5, Re=200,000

| Parameter | Value | Notes |
|-----------|-------|-------|
| CD0 (α=0°) | **0.01046** | Minimum drag; symmetric airfoil so α=0 = minimum |
| CD at α=2° | 0.01103 | Lightly loaded hover |
| CD at α=4° | 0.01273 | Moderate hover loading |
| CD at α=6° | 0.01538 | High hover loading |
| CD at α=8° | 0.01865 | Near stall |
| CL_alpha (0–2°) | **6.05 /rad** | = 0.1056/deg; consistent with Table VIII 5.90/rad |
| CL_alpha (0–4°) | 6.30 /rad | slight non-linearity above 3° |
| CL_max | 1.1583 | at α ≈ 15.5° |
| Alpha stall | **~15°** | gradual; CL starts dropping above 15.5° |
| Transition (α=0°) | 69% chord | both surfaces; NCrit=5 |

### NCrit=5 vs NCrit=9 comparison at α=0°

| NCrit | CD0    | Transition (top/bot) |
|-------|--------|----------------------|
| 5     | 0.01046 | 69% / 69% |
| 9     | 0.01071 | 81% / 81% |

Counter-intuitive result: NCrit=9 gives slightly higher CD0 despite later natural transition.
At Re=200k this is likely a laminar separation bubble near the trailing edge, which XFOIL
models as an extra drag source when NCrit=9 triggers late turbulent transition.

---

## Implications for BEM fixture

### CD0
Current fixture value: `CD0=0.012`
XFOIL NCrit=5 value: `CD0=0.01046` (12% lower)

The difference is most significant at low collective (α ≈ 0–2°) where profile drag
dominates the torque budget.  At low collective, σ·CD0/8 ≈ 7.5×10⁻⁵ using 0.012, vs
≈ 6.5×10⁻⁵ using 0.01046.  This directly caused the CQ overestimate at θ=4.91°.

**Recommended: update CD0 to 0.0105** (rounds XFOIL value; within Re uncertainty).

### CL_alpha
XFOIL gives 6.05/rad; Table VIII (actual measured) gives 5.90/rad.
**Keep 5.90/rad** — measured test data takes precedence over XFOIL prediction.

### Alpha stall
XFOIL (NCrit=5, Re=200k) gives CL_max at ~15.5°.
Current fixture: `alpha_stall_deg=12.0` — conservative but fine for rotor BEM since
blade sections rarely reach 12° angle of attack at the thrust levels tested.

---

## CD as a function of alpha (NCrit=5, Re=200,000) — selected points

| α (deg) | CL     | CD      | CL/CD |
|---------|--------|---------|-------|
| 0       | 0.0000 | 0.01046 | —     |
| 1       | 0.1049 | 0.01060 | 9.9   |
| 2       | 0.2112 | 0.01103 | 19.1  |
| 3       | 0.3220 | 0.01177 | 27.3  |
| 4       | 0.4399 | 0.01273 | 34.6  |
| 5       | 0.5678 | 0.01391 | 40.8  |
| 6       | 0.7028 | 0.01538 | 45.7  |
| 7       | 0.8287 | 0.01706 | 48.6  |
| 8       | 0.9075 | 0.01905 | 47.6  |
| 9       | 0.9508 | 0.02054 | 46.3  |
| 10      | 1.0070 | 0.02279 | 44.2  |
| 12      | 1.0899 | 0.02970 | 36.7  |
| 14      | 1.1454 | 0.04265 | 26.9  |
| 16      | 1.1557 | 0.06457 | 17.9  |

Note: CD is not constant — it rises significantly above 5°.  The BEM model uses a constant
CD0, which underestimates drag at high angles of attack (overestimates thrust at high
collective) and overestimates drag at low angles (overestimates torque at low collective).
A full polar lookup would improve accuracy across the CT range.
