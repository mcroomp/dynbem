# Master Field Mapping Table

Comprehensive single-reference table of all Rust struct fields and their Python property names.

---

| **Struct Name** | **Rust Field** | **Python Property** | **Type** | **Match?** | **Notes** |
|-----------------|----------------|-------------------|---------|-----------|----------|
| **RotorInputs** | collective_rad | .collective_rad | f64 | ✓ | r/w |
| | tilt_lon | .tilt_lon | f64 | ✓ | r/w |
| | tilt_lat | .tilt_lat | f64 | ✓ | r/w |
| | R_hub | .R_hub | numpy (3,3) | ✓ | r/o |
| | v_hub_world | .v_hub_world | numpy (3,) | ✓ | r/o |
| | wind_world | .wind_world | numpy (3,) | ✓ | r/o |
| | rho_kg_m3 | .rho_kg_m3 | f64 | ✓ | r/w |
| | omega_rad_s | .omega_rad_s | f64 | ✓ | r/w |
| | *(Python only)* | .t | f64 | - | r/w; timestamp |
| **AeroResult** | F_world | .F_world | numpy (3,) | ✓ | r/o |
| | M_hub_world | .M_orbital | numpy (3,) | **✗** | **RENAMED** |
| | Q_spin | .Q_spin | f64 | ✓ | r/o |
| | M_spin | .M_spin | numpy (3,) | ✓ | r/o |
| **BladeGeometry** | n_blades | .n_blades | usize | ✓ | r/o |
| | radius_m | .radius_m | f64 | ✓ | r/o |
| | root_cutout_m | .root_cutout_m | f64 | ✓ | r/o |
| | chord_m | .chord_m | f64 | ✓ | r/o |
| | twist_deg | .twist_deg | f64 | ✓ | r/o |
| | n_elements | .n_elements | usize | ✓ | r/o |
| | tip_loss | .tip_loss | bool | ✓ | r/o |
| | r_stations_m | *(internal)* | Vec<f64> | ✓ | not exposed |
| | chord_stations_m | *(internal)* | Vec<f64> | ✓ | not exposed |
| | twist_stations_deg | *(internal)* | Vec<f64> | ✓ | not exposed |
| | *(computed)* | .span_m() | f64 | ✓ | method |
| | *(computed)* | .r_cp_m() | f64 | ✓ | method |
| | *(computed)* | .disk_area_m2() | f64 | ✓ | method |
| | *(computed)* | .solidity() | f64 | ✓ | method |
| | *(computed)* | .has_radial_stations() | bool | ✓ | method |
| | *(computed)* | .chord_at(r) | f64 | ✓ | method |
| | *(computed)* | .twist_at(r) | f64 | ✓ | method |
| **LinearPolarParameters** | CL0 | .CL0 | f64 | ✓ | r/o |
| | CL_alpha_per_rad | .CL_alpha_per_rad | f64 | ✓ | r/o |
| | CD0 | .CD0 | f64 | ✓ | r/o |
| | alpha_stall_deg | .alpha_stall_deg | f64 | ✓ | r/o; degrees |
| **LinearPolar** | CL0 | .CL0 | f64 | ✓ | r/o |
| | CL_alpha_per_rad | .CL_alpha_per_rad | f64 | ✓ | r/o |
| | CD0 | .CD0 | f64 | ✓ | r/o |
| | alpha_stall_rad | .alpha_stall_rad | f64 | ✓ | r/o; **radians** |
| | *(method)* | .cl_cd(alpha) | (f64, f64) | ✓ | returns (CL, CD) |
| | *(method)* | .cl_cd_arr(alpha) | (array, array) | ✓ | vectorized |
| **ControlProperties** | swashplate_pitch_gain_rad | .swashplate_pitch_gain_rad | f64 | ✓ | r/o |
| | swashplate_phase_deg | .swashplate_phase_deg | Option<f64> | ✓ | r/o |
| **ServoFlapGeometry** | C_M_delta_per_rad | .C_M_delta_per_rad | f64 | ✓ | r/o |
| | r_inner_m | .r_inner_m | f64 | ✓ | r/o |
| | r_outer_m | .r_outer_m | f64 | ✓ | r/o |
| **ServoFlapActuation** | I_theta_kgm2 | .I_theta_kgm2 | f64 | ✓ | r/o |
| | damper_Nms_per_rad | .damper_Nms_per_rad | f64 | ✓ | r/o |
| | ac_offset_m | .ac_offset_m | f64 | ✓ | r/o |
| | blade_Cm_AC | .blade_Cm_AC | f64 | ✓ | r/o |
| | flap | .flap | ServoFlapGeometry | ✓ | r/o; nested |
| **FlapProperties** | I_blade_flap_kgm2 | .I_blade_flap_kgm2 | f64 | ✓ | r/o |
| | omega_nr_rad_s | .omega_nr_rad_s | f64 | ✓ | r/o |
| | *(method)* | .hub_moment_factor(ω) | f64 | ✓ | computed |
| **RotorDefinition** | blade | .blade | BladeGeometry | ✓ | r/o; nested |
| | airfoil | .airfoil | LinearPolarParameters | ✓ | r/o; nested |
| | control | .control | Option<ControlProperties> | ✓ | r/o |
| | pitch_actuation | .servoflap | Option<ServoFlapActuation> | **✗** | **ENUM EXTRACT** |
| | flap | .flap | Option<FlapProperties> | ✓ | r/o |
| | name | .name | &str | ✓ | r/o |
| | description | .description | &str | ✓ | r/o |
| **QuasiStaticRotorState** | *(empty)* | .to_array() | numpy () | ✓ | empty array |
| | *(empty)* | .from_array(arr) | Self | ✓ | validates len==0 |
| **PittPetersRotorState** | lambda_0 | .lambda_0 | f64 | ✓ | r/w |
| | lambda_c | .lambda_c | f64 | ✓ | r/w |
| | lambda_s | .lambda_s | f64 | ✓ | r/w |
| | *(method)* | .to_array() | numpy (3,) | ✓ | [λ₀, λc, λs] |
| | *(method)* | .from_array(arr) | Self | ✓ | validates len==3 |
| | *(setters)* | .set_lambda_0(x) | None | ✓ | via setter |
| | *(setters)* | .set_lambda_c(x) | None | ✓ | via setter |
| | *(setters)* | .set_lambda_s(x) | None | ✓ | via setter |
| **OyeRotorState** | n_elements | *(internal)* | usize | ✓ | not exposed |
| | W_int | .W_int | numpy (n,) | ✓ | r/o; getter only |
| | W | .W | numpy (n,) | ✓ | r/o; getter only |
| | *(method)* | .to_array() | numpy (2n,) | ✓ | [W_int, W] flat |
| | *(method)* | .from_array(arr) | Self | ✓ | validates len==2n |
| | *(static)* | .zeros(n) | Self | ✓ | factory |
| **VpmRotorState** | *(TBD)* | *(TBD)* | ? | ? | see vpm/mod.rs |
| **Model Classes** | defn | .defn | RotorDefinition | ✓ | r/o |
| (all 6 variants) | n_psi_elements | .n_psi_elements | usize | ✓ | r/o |
| | *(method)* | .initial_rotor_state() | State | ✓ | factory |
| | *(method)* | .compute_forces(i, s) | (AeroResult, State) | ✓ | no advance |
| | *(method)* | .step(i, s, dt, method) | (AeroResult, State) | ✓ | integrates state |
| | *(method)* | .inflow_taus(i, s) | numpy (n_dof,) | ✓ | time constants |
| **OyeBEMModel** | coupling_k | .coupling_k | f64 | ✓ | r/o; Oye specific |

---

## Legend

| Symbol | Meaning |
|--------|---------|
| ✓ | Field names match exactly |
| **✗** | Field names differ (mismatch) |
| r/w | Read-write (getter + setter) |
| r/o | Read-only (getter only) |
| - | N/A (Python-only field) |
| *(method)* | Computed/helper method, not a direct field |
| *(internal)* | Stored in Rust but not exposed as Python property |
| *(Python only)* | Not in Rust struct; added by Python wrapper |
| **RENAMED** | Name changed intentionally for API consistency |
| **ENUM EXTRACT** | Extracted from enum variant match |

---

## Structs by Category

### Control/Input
- `RotorInputs` — 8 fields + 1 Python-only
- `ControlProperties` — 2 fields
- `LinearPolarParameters` — 4 fields
- `LinearPolar` — 4 fields + 2 methods

### Geometry
- `BladeGeometry` — 10 fields + 7 methods
- `ServoFlapGeometry` — 3 fields
- `FlapProperties` — 2 fields + 1 method

### Actuation
- `ServoFlapActuation` — 5 fields (nested flap)
- `RotorDefinition` — 7 fields (all nested or enums)

### Output
- `AeroResult` — 4 fields (⚠️ 1 renamed)

### State
- `QuasiStaticRotorState` — 0 fields (quasi-static)
- `PittPetersRotorState` — 3 fields + 2 methods + 3 setters
- `OyeRotorState` — 3 fields + 2 methods (W_int, W vectors)
- `VpmRotorState` — *(not detailed)*

### Models (6 variants)
- `_QuasiStaticBEM{Linear,Tabulated}` — same interface
- `_PittPetersModel{Linear,Tabulated}` — same interface
- `_OyeBEMModel{Linear,Tabulated}` — same interface + coupling_k
- `_VpmRotor{Linear,Tabulated}` — same interface

---

## Critical Implementation Notes

1. **AeroResult.M_orbital ← Rust M_hub_world**
   - Users MUST use `.M_orbital` in Python code
   - Rust side still uses `M_hub_world` internally
   - Intentional rename for clarity (moment is in orbital/hub frame)

2. **RotorDefinition.servoflap ← Rust pitch_actuation enum**
   - Property extracts `ServoFlapActuation` from `PitchActuation::ServoFlap(...)`
   - Returns `None` for `DirectMechanical` mode
   - No way to access the enum itself from Python

3. **LinearPolar vs LinearPolarParameters units**
   - `LinearPolarParameters.alpha_stall_deg` — degrees
   - `LinearPolar.alpha_stall_rad` — radians
   - Constructor auto-converts: `to_radians()` called in `LinearPolar::new()`

4. **RotorInputs.t storage**
   - Stored in `PyRotorInputs.t`, not in Rust `inner.RotorInputs`
   - Separate field in Python wrapper for time tracking
   - Not part of BEM/aerodynamic calculation

5. **OyeRotorState properties are read-only**
   - `.W_int` and `.W` are getters that return copies (Vec→numpy)
   - Modification requires `.set_inflow()` or state replacement
   - Per-annulus states; total length = 2 × n_elements

---

## Version Information

- **dynbem_rs version:** 0.5.0
- **dynbem Python version:** 0.5.0
- **PyO3 binding:** maturin-managed
- **Mapping validated:** v0.5.0 release (VPM theory validation + docs overhaul)
