# NACA TR 515 — Wheatley & Hood (1935)

**Citation.** Wheatley, J.B. & Hood, M.J. (1935). *Full-Scale Wind-Tunnel
Tests of a PCA-2 Autogiro Rotor.* NACA Technical Report 515.

**Source PDF.** `NACA-Report-515.pdf` (13 MB, from NTRS citation 19930091587).

**Companion paper.** NACA TR 487 (Wheatley 1934, in
`Research/Wheatley_NACA487/`) develops the analytical autogyro BEM theory;
this paper is the experimental dataset that validates it. The modern
re-analysis by Harris (in `Research/Harris_CR-2008-215370/`) digitizes
this exact data in rotor-axis coordinates.

## Rotor geometry (Pitcairn PCA-2)

- 4 blades, untwisted, NACA airfoils (22-3/4 in symmetric outer + 14-25/32 in
  cambered inner — see Harris §11.5 for the blade reconstruction)
- Diameter = 45 ft (R = 6.86 m)
- Tested at nominal Ω = 98.6, 118.7, 137.6, 147.9 rpm
  (ΩR ≈ 232, 280, 324, 348 ft/s)
- Pitch settings (collective): 0.8°, 1.9°, 2.7°
- Tunnel speeds 32–173 ft/s, giving μ ≈ 0.1–0.7
- No swashplate, no feathering — the rotor simply autorotates.

## Page index

| File | Content |
|---|---|
| `cover.png` | Cover (Report No. 515, 1935) |
| `symbols_front.png` | NACA Aeronautic Symbols front-matter |
| `page_01.png` … `page_08.png` | Body of paper (test procedure, instrumentation, methodology, conclusions) |
| `page_09.png` | **Tables I and II** (data; Tables I=pitch 1.9° exposed, II=pitch 0.8° faired) |
| `page_10.png` | **Tables III and IV** (data; III=pitch 1.9° faired, IV=pitch 2.7° faired) |
| `nomenclature.png` | Coordinate-axis diagram and force/moment definitions (bound at end of report) |
| `page_09_table_i.png`, `page_09_table_ii.png` | Cropped table images for transcription |
| `page_10_table_iii.png`, `page_10_table_iiv.png` | Same for page 10 (filename `iiv.png` is `iv.png`) |

## Data tables — extraction status

| Table | Pitch | Protuberances | Rows | Confidence | File |
|---|---|---|---|---|---|
| I   | 1.9° | Exposed | 60 (of 87 raw) | MODERATE — PDF text extraction + consistency filter | `page_09_table_i.md` |
| II  | 0.8° | Faired  | 73 (partial)   | MODERATE — same approach as Table I | `page_09_table_ii.md` |
| III | 1.9° | Faired  | 79             | **HIGH — manually transcribed from cropped PNG; 0/79 rows fail L/D consistency check** | `page_10_table_iii.md` |
| IV  | 2.7° | Faired  | 74             | **HIGH — manually transcribed; 0/74 rows fail L/D consistency check** | `page_10_table_iv.md` |

Tables I and II were partially recoverable from the PDF's OCR text layer.
Tables III and IV had heavy OCR damage on page 10 (the μ column was
entirely missing from the text layer), so they were transcribed by hand
from the high-resolution cropped PNGs and every row was internally
validated against L/D = CL/CD within 10%.

## Coordinate system in this paper

Coefficients are in **airplane axes**:
- CL = lift / (q · π R²)
- CD = drag / (q · π R²)
- L/D — lift-drag ratio
- C_l, C_m — rolling / pitching moment coefficients about hub
- C_y — lateral force coefficient

where q = ½ρV² (tunnel dynamic pressure) and π R² is rotor disk area.

This is **not** the rotor-axis (CT, CH, CY) form that BEM codes use.
To convert, rotate by the shaft angle of attack α (column 2). See
Harris CR-2008-215370 §5.3 for the explicit transformation and his
digitized rotor-axis figures (Figs 5-17 to 5-20).

## Validation use

For `dynbem` BEM validation against autorotation in forward flight:
- Tables III and IV (manually transcribed, HIGH confidence) are the
  primary reference.
- Convert tabulated (CL, CD, α) to (CT, CH) via the shaft rotation, then
  compare against `BEMModel.compute_forces` at the same operating point.
- Alternatively, validate directly against the rotor-axis curves Harris
  digitized in his Figs 5-17 to 5-20 — easier because the transformation
  is already done.
