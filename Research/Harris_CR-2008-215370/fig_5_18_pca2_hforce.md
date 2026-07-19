# Fig. 5-18 -- PCA-2 Rotor H-Force Coefficient vs Advance Ratio

**Source:** Harris, F.D. (2008). *Introduction to Autogyros, Helicopters,
and Other V/STOL Aircraft* -- NASA/CR-2008-215370, Fig. 5-18 (p.88).
Re-analysis of the Wheatley & Hood PCA-2 wind-tunnel data (NACA TR 515),
converted to rotor-axis coefficients.

**Rotor:** Pitcairn PCA-2 autogyro, 4 blades, D = 45 ft (R = 6.86 m),
untwisted, autorotating (no shaft power). Same rotor as
`Research/Wheatley_Hood_NACA515/`.

## Coefficient definition (confirm against Harris symbol page)

Harris uses the American (Army/NACA) coefficient convention, NOT the
1/2-rho dynamic-pressure form:

    CH = H / (rho * pi * R^2 * (Omega*R)^2)

where H is the rotor in-plane (H-) force along the rotor-axis drag
direction, rho air density, R tip radius, Omega rotor speed. Fig. 5-18
plots the BARE CH (Fig. 5-25 by contrast plots CH/sigma). The exact
reference area/velocity should be verified against Harris's front-matter
symbol list (image-only in the PDF) before using these numbers
dimensionally. To recover a force: H = CH * rho * pi * R^2 * (Omega*R)^2.

Nominal test RPMs (from the paper): 98.6, 118.7, 137.6, 147.9, giving
Omega*R ~ 232, 280, 324, 348 ft/s (~70.7, 85.3, 98.8, 106.1 m/s).

## Confidence: LOW (visual digitization)

These points were read by eye from the rendered scatter plot, NOT from a
tabulated source. Estimated reading tolerance: mu +/- 0.01,
CH +/- 0.00003. Use as a trend/order-of-magnitude validation target. The
higher-fidelity primary data for this same rotor is the Wheatley & Hood
tables (`Research/Wheatley_Hood_NACA515/page_10_table_iii.md`, `_iv.md`),
which give CD/C_y in airplane axes at HIGH confidence; Fig. 5-18 is
Harris's rotor-axis re-projection of that data across all four RPM runs.

The clear physical trend (the point of the figure): CH rises roughly
linearly with advance ratio, ~0.0002 at mu ~ 0.14 to ~0.001 at mu ~ 0.7,
and is essentially independent of RPM (all four series collapse onto one
curve).

## Data (grouped by nominal RPM series)

Column key -- mu: advance ratio | CH: rotor-axis H-force coefficient |
rpm_series: nominal RPM marker series in the figure.

| mu    | CH       | rpm_series |
|-------|----------|------------|
| 0.13  | 0.00020  | 98.6       |
| 0.17  | 0.00024  | 98.6       |
| 0.205 | 0.00027  | 98.6       |
| 0.25  | 0.00030  | 98.6       |
| 0.27  | 0.00037  | 98.6       |
| 0.31  | 0.00039  | 98.6       |
| 0.345 | 0.00041  | 98.6       |
| 0.38  | 0.00049  | 98.6       |
| 0.41  | 0.00065  | 98.6       |
| 0.44  | 0.00057  | 98.6       |
| 0.46  | 0.00062  | 98.6       |
| 0.49  | 0.00066  | 98.6       |
| 0.53  | 0.00072  | 98.6       |
| 0.59  | 0.00078  | 98.6       |
| 0.635 | 0.00086  | 98.6       |
| 0.70  | 0.00097  | 98.6       |
| 0.14  | 0.00020  | 118.7      |
| 0.18  | 0.00021  | 118.7      |
| 0.22  | 0.00030  | 118.7      |
| 0.245 | 0.00040  | 118.7      |
| 0.285 | 0.00043  | 118.7      |
| 0.34  | 0.00047  | 118.7      |
| 0.375 | 0.00055  | 118.7      |
| 0.42  | 0.00062  | 118.7      |
| 0.50  | 0.00066  | 118.7      |
| 0.595 | 0.00077  | 118.7      |
| 0.605 | 0.00085  | 118.7      |
| 0.135 | 0.00007  | 137.6      |
| 0.205 | 0.00016  | 137.6      |
| 0.21  | 0.00028  | 137.6      |
| 0.24  | 0.00036  | 137.6      |
| 0.275 | 0.00042  | 137.6      |
| 0.30  | 0.00046  | 137.6      |
| 0.335 | 0.00050  | 137.6      |
| 0.37  | 0.00053  | 137.6      |
| 0.41  | 0.00060  | 137.6      |
| 0.45  | 0.00062  | 137.6      |
| 0.49  | 0.00070  | 137.6      |
| 0.52  | 0.00075  | 137.6      |
| 0.14  | 0.00018  | 147.9      |
| 0.175 | 0.00009  | 147.9      |
| 0.225 | 0.00031  | 147.9      |
| 0.28  | 0.00046  | 147.9      |
| 0.32  | 0.00054  | 147.9      |
| 0.35  | 0.00060  | 147.9      |
| 0.415 | 0.00060  | 147.9      |
| 0.135 | 0.00019  | transition |
| 0.72  | 0.00107  | transition |
