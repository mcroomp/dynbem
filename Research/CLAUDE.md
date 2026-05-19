# Research — Source materials

## Extraction convention

All tables and figures extracted from a paper are saved as markdown files in the paper's subfolder,
named with the page image as a prefix: `page_NN_<description>.md`
(e.g. `page_32_table_i.md`, `page_42_figure_3.md`).
This keeps extractions traceable back to their source page image.

## CaradonnaTung/
NASA TM-81232 (Caradonna & Tung, 1981) — hover blade-pressure and wake-geometry study.
2-blade NACA 0012, R=1.143 m, σ=0.1062. CT data at θc=5°/8°/12°, Ω=1250–2500 rpm.
**No CP/torque data.** Primary BEM validation source for Level 1 and Level 3.
See [CaradonnaTung/CLAUDE.md](CaradonnaTung/CLAUDE.md) for page index, CT tables, and test notes.

## Buhl_NREL_TP500_36834/
NREL TP-500-36834 (Buhl, 2005) — Windmill Brake State correction extending Glauert's
empirical correction for a > 0.4. Used for the WBS quadratic in the BEM solver.

## Castles_TN2474/
NACA TN-2474 (Castles & Gray, 1951) — induced velocity measurements in helicopter hover/descent.
Contains the experimental data underlying the Leishman VRS polynomial used in Pitt-Peters.

## Harrington/ and Harrington_TN2318/
NACA TN-2318 (Harrington, 1951) — hover CT vs CP polars for two full-scale rotors.
Candidate dataset for CP-CT polar validation (this paper has torque data; Caradonna-Tung does not).

## Peters_Nikolsky_2008/
Peters, D.A. (2009). "How Dynamic Inflow Survives in the Competitive World of
Rotorcraft Aerodynamics: The Alexander Nikolsky Honorary Lecture," *JAHS*
54(1):011001. The **canonical reference for the Pitt-Peters 3-state dynamic
inflow model** — Eqs 7–11 give the ODE, V, M, L matrices, and inflow
distribution verbatim. See [Peters_Nikolsky_2008/CLAUDE.md] for the
sign-convention translation between Peters' conventions (ψ=0 at tail) and
ours (ψ=0 at hub +X).
