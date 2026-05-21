# NREL Phase VI — Unsteady Aerodynamics Experiment

The PCA-2 of wind turbines: a 2-blade, 10.058 m diameter, S809-airfoil
horizontal-axis wind turbine tested in the NASA Ames 80 ft x 120 ft
wind tunnel by NREL (May 2000). Validated against IEA Annex XX
blind-comparison.

## Source documents

| File | Reference | What's in it |
|---|---|---|
| `TP-500-29955/source.pdf` | Hand, Simms, Fingersh, Jager, Cotrell, Schreck, Larwood (Dec 2001). NREL/TP-500-29955. | The canonical test-configurations / instrumentation document, 310 pages. Includes rotor geometry, instrumentation, sequence definitions, file-by-file data inventory. |
| `TP-500-43508/source.pdf` | IEA Wind Annex XX Final Report (Dec 2008). NREL/TP-500-43508. | The blind-comparison report, 91 pages. Tabulated measured power, torque, and blade-element Cn at multiple radial stations for Sequence S and others. **The most directly usable artifact for BEM validation.** |
| `SR-440-6918/source.pdf` | Somers (1997). NREL/SR-440-6918, *Design and Experimental Results for the S809 Airfoil*. | Wind-tunnel polar paper for the S809 airfoil used on the Phase VI rotor, 103 pages. Coordinates, pressure distributions, lift/drag/moment polars at Re = 1.0e6 - 3.0e6. |

Each subdirectory holds the corresponding source PDF (`source.pdf`),
the extracted page PNGs (named by **printed** page number), and any
data MDs derived from this paper.

## Data archive status

The **raw time-series measurements are not freely downloadable**. As of
the most recent public NREL forum statement (2017), Lee Jay Fingersh
at NREL is the contact for the dataset. No GitHub mirror, no OSTI
release, no Zenodo upload. For our validation work we extract
operating-point data from the figures and tables in TP-500-43508 and
TP-500-29955 - the same workflow as Castles-Gray, Caradonna-Tung,
Wheatley.

## Page index — printed page numbers (PNG filename = printed page)

### TP-500-29955 (offset: pdf page = printed page + 10)

| Printed page | Content |
|---|---|
| 12 (`TP-500-29955/page_012.png`) | Blade geometry description: S809 airfoil, tapered + twisted, designed by Giguere & Selig (1999) |
| 21-22 | Sequence S, T, U descriptions: upwind, 3 deg tip pitch (S), 2 deg (T), 4 deg (U); 72 RPM; 5-25 m/s wind |
| 28 | 5-hole probe description (Sequence S has probes removed) |
| 63-66 | Appendix A — rotor geometry. **Table A-1 (printed 64) = chord/twist by radial station.** This is the key fixture-building table. |
| 202-203 | Table C-18: Sequence S wind-speed x yaw-angle test matrix |

### TP-500-43508 (offset: pdf page = printed page + 6)

| Printed page | Content |
|---|---|
| 4-5 (`TP-500-43508/page_004.png` / `005.png`) | Executive summary, list of participating organizations |
| 51 | **Fig. 4(a): measured LSS torque vs wind speed (Sequence S, axial flow).** Compares free-wake, BEM-with-Prandtl, BEM-without-Prandtl, and pressure-measured / strain-gauge measurements. |
| 52-55 | Power, root flap moment, root edge moment comparisons across the same wind speeds |
| 56-57 | Blade-element Cn, Cy at four radial stations vs wind speed |
| 77-80 | Sequence H / Sequence S detailed blade-element comparisons (multiple yaw angles) |

### SR-440-6918 (offset: pdf page = printed page + 5)

| Printed page | Content |
|---|---|
| `cover.png` | Title page (no printed number) |
| 5 (`SR-440-6918/page_005.png`) | Body intro |
| 10 | Principal results summary |
| 15 | **Table 2: S809 airfoil coordinates (x/c, z/c) for upper and lower surfaces** — needed to construct the airfoil shape |
| 16 | Table 3: pressure-tap orifice locations |
| 18 | Table 4: roughness grit specification for boundary-layer trip tests |
| 45 | Figure 10: measured lift coefficient vs angle of attack |
| 65 | Figure 13/14: drag polars at multiple Reynolds numbers |

## Validation roadmap

Goal: validate `dynbem.bem` against Sequence S measured power/torque
across the 5-25 m/s wind-speed range at 0 deg yaw, fixed pitch 3 deg.

1. **Build the rotor fixture** (`rotors/nrel_phase_vi/rotor.yaml`)
   - 2 blades, R = 5.029 m, 72 RPM (Ω ≈ 7.54 rad/s, ΩR ≈ 37.9 m/s)
   - Chord and twist by radial station from Table A-1 (TP-500-29955 p64)
   - S809 airfoil polar — most practical source is NREL's
     machine-readable `s809_clean.dat` from the OpenFAST distribution,
     or transcribed from Somers Figures 10/13.

2. **Extract Sequence S measured data tables** (TP-500-43508 pp.51-55)
   into MD then CSV via `extract_tables.py`. Need at least:
   - LSS torque vs wind speed at 0 deg yaw (Fig. 4a)
   - Power vs wind speed
   - Root flap moment vs wind speed

3. **Write `tests/test_nrel_phase_vi.py`** mirroring the Wheatley
   structure: at each tabulated wind speed, run the BEM in wind-turbine
   mode (`wind_world` along +X, shaft at 0 deg tilt, blade pitch +3 deg
   relative to chord line — sign opposite to helicopter convention!),
   assert torque/power within a tolerance band.

## Convention reminders

- Wind turbine has shaft **horizontal** (perpendicular to vertical),
  not vertical. R_hub will reflect that.
- In NED with wind blowing in +X: turbine extracts power, so Q_aero
  is positive in the rotation direction. Sign of CT is reversed from
  helicopter convention (thrust on rotor is downwind, not upward).
- The S809 is a **fixed-pitch stall-regulated** airfoil; above ~10 m/s
  the rotor goes into stall regulation. The BEM must handle the
  high-CD post-stall regime - this is exactly where Buhl's WBS quadratic
  is relevant, even though Phase VI is technically not in WBS.
