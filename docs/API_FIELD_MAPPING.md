# dynbem Rust-to-Python API Field Mapping

Comprehensive analysis of all public structs, field names, and how they're exposed through Python bindings.

---

## Core I/O Types (aero_io.rs)

### Vec3

| Rust Field | Python Exposure | Match? | Notes |
|------------|-----------------|--------|-------|
| `[0]` (x) | numpy array access `[0]` | Y | Nested in Mat3/Vec3 arrays |
| `[1]` (y) | numpy array access `[1]` | Y | |
| `[2]` (z) | numpy array access `[2]` | Y | |

**Export:** Via `vec3_to_py()` helper converts to numpy 1D arrays

---

### Mat3

| Rust Field | Python Exposure | Match? | Notes |
|------------|-----------------|--------|-------|
| `[[0][0]]` (3x3) | numpy 2D array `shape=(3,3)` | Y | Converted via `mat3_to_py()` |

**Export:** Via `mat3_to_py()` helper converts to numpy 2D arrays

---

### RotorInputs (Most Used)

| Rust Field Name | Python Property Name | Match? | Type | Notes |
|-----------------|----------------------|--------|------|-------|
| `collective_rad` | `.collective_rad` | Y | f64 | Scalar, r/w |
| `tilt_lon` | `.tilt_lon` | Y | f64 | Scalar, r/w |
| `tilt_lat` | `.tilt_lat` | Y | f64 | Scalar, r/w |
| `R_hub` | `.R_hub` | Y | numpy (3,3) | Matrix, read-only |
| `v_hub_world` | `.v_hub_world` | Y | numpy (3,) | Vector, read-only |
| `wind_world` | `.wind_world` | Y | numpy (3,) | Vector, read-only |
| `rho_kg_m3` | `.rho_kg_m3` | Y | f64 | Scalar, r/w |
| `omega_rad_s` | `.omega_rad_s` | Y | f64 | Scalar, r/w |
| **NEW (Python only)** | `.t` | - | f64 | Timestamp (not in Rust RotorInputs) |

**PyClass Name:** `_dynbem.RotorInputs`  
**Python Constructor:**
```python
RotorInputs(collective_rad, tilt_lon, tilt_lat, R_hub, 
            v_hub_world, wind_world, omega_rad_s, t, rho_kg_m3)
```

**Key Mismatch:** `t` (timestamp) is stored in `PyRotorInputs.t`, not in Rust `inner`. Used for external time tracking.

---

### AeroResult

| Rust Field Name | Python Property Name | Match? | Type | Notes |
|-----------------|----------------------|--------|------|-------|
| `F_world` | `.F_world` | Y | numpy (3,) | World force vector |
| `M_hub_world` | `.M_orbital` | **N** | numpy (3,) | **RENAMED!** Rust uses M_hub_world |
| `Q_spin` | `.Q_spin` | Y | f64 | Scalar torque |
| `M_spin` | `.M_spin` | Y | numpy (3,) | Spin moment vector |

**PyClass Name:** `_dynbem.AeroResult`  
**Important:** The field `M_hub_world` is aliased as `.M_orbital` in Python for backward compat.

---

## Rotor Definition Types (rotor_definition.rs)

### BladeGeometry

| Rust Field | Python Property | Match? | Type | Notes |
|------------|-----------------|--------|------|-------|
| `n_blades` | `.n_blades` | Y | usize | Read-only |
| `radius_m` | `.radius_m` | Y | f64 | Read-only |
| `root_cutout_m` | `.root_cutout_m` | Y | f64 | Read-only |
| `chord_m` | `.chord_m` | Y | f64 | Read-only |
| `twist_deg` | `.twist_deg` | Y | f64 | Read-only |
| `n_elements` | `.n_elements` | Y | usize | Read-only |
| `tip_loss` | `.tip_loss` | Y | bool | Read-only |
| `r_stations_m` | *(internal)* | Y | Vec<f64> | Not exposed as property |
| `chord_stations_m` | *(internal)* | Y | Vec<f64> | Not exposed as property |
| `twist_stations_deg` | *(internal)* | Y | Vec<f64> | Not exposed as property |
| **Methods:** | `.span_m()` | Y | f64 | Computed: `radius_m - root_cutout_m` |
| | `.r_cp_m()` | Y | f64 | Computed: pressure center |
| | `.disk_area_m2()` | Y | f64 | Computed |
| | `.solidity()` | Y | f64 | Computed |
| | `.has_radial_stations()` | Y | bool | Method |
| | `.chord_at(r)` | Y | f64 | Interpolate |
| | `.twist_at(r)` | Y | f64 | Interpolate |

**PyClass Name:** `_dynbem.BladeGeometry`  
**Python Class Wrapper:** Also in `dynbem.rotor_definition.BladeGeometry` (subclass with YAML loading)

---

### LinearPolarParameters

| Rust Field | Python Property | Match? | Type | Notes |
|------------|-----------------|--------|------|-------|
| `CL0` | `.CL0` | Y | f64 | Lift coefficient at α=0 |
| `CL_alpha_per_rad` | `.CL_alpha_per_rad` | Y | f64 | dCL/dα [rad^-1] |
| `CD0` | `.CD0` | Y | f64 | Drag coefficient at α=0 |
| `alpha_stall_deg` | `.alpha_stall_deg` | Y | f64 | Stall angle [deg] |

**PyClass Name:** `_dynbem.LinearPolarParameters`  
**Note:** Rust stores `alpha_stall_deg` in degrees; Python also uses degrees (no conversion on property access, unlike LinearPolar which uses radians)

---

### LinearPolar (Derived from LinearPolarParameters)

| Rust Field | Python Property | Match? | Type | Notes |
|------------|-----------------|--------|------|-------|
| `CL0` | `.CL0` | Y | f64 | |
| `CL_alpha_per_rad` | `.CL_alpha_per_rad` | Y | f64 | |
| `CD0` | `.CD0` | Y | f64 | |
| `alpha_stall_rad` | `.alpha_stall_rad` | Y | f64 | **Units: RADIANS** |

**PyClass Name:** `_dynbem.LinearPolar`  
**Key Difference:** `alpha_stall_rad` uses radians (constructor argument is in radians)

---

### ControlProperties

| Rust Field | Python Property | Match? | Type | Notes |
|------------|-----------------|--------|------|-------|
| `swashplate_pitch_gain_rad` | `.swashplate_pitch_gain_rad` | Y | f64 | Gain [rad] |
| `swashplate_phase_deg` | `.swashplate_phase_deg` | Y | Option<f64> | Phase angle [deg] or None |

**PyClass Name:** `_dynbem.ControlProperties`

---

### ServoFlapGeometry

| Rust Field | Python Property | Match? | Type | Notes |
|------------|-----------------|--------|------|-------|
| `C_M_delta_per_rad` | `.C_M_delta_per_rad` | Y | f64 | Pitching moment coeff [rad^-1] |
| `r_inner_m` | `.r_inner_m` | Y | f64 | Flap inner radius [m] |
| `r_outer_m` | `.r_outer_m` | Y | f64 | Flap outer radius [m] |

**PyClass Name:** `_dynbem.ServoFlapGeometry`

---

### ServoFlapActuation

| Rust Field | Python Property | Match? | Type | Notes |
|------------|-----------------|--------|------|-------|
| `I_theta_kgm2` | `.I_theta_kgm2` | Y | f64 | Pitch inertia [kg⋅m²] |
| `damper_Nms_per_rad` | `.damper_Nms_per_rad` | Y | f64 | Damping [N⋅m⋅s/rad] |
| `ac_offset_m` | `.ac_offset_m` | Y | f64 | AC distance aft of feathering axis [m] (aero spring) |
| `blade_Cm_AC` | `.blade_Cm_AC` | Y | f64 | Blade zero-lift pitching moment coeff [-] (DC trim) |
| `flap` | `.flap` | Y | ServoFlapGeometry | Nested struct |

**PyClass Name:** `_dynbem.ServoFlapActuation`

---

### FlapProperties

| Rust Field | Python Property | Match? | Type | Notes |
|------------|-----------------|--------|------|-------|
| `I_blade_flap_kgm2` | `.I_blade_flap_kgm2` | Y | f64 | Flap inertia [kg⋅m²] |
| `omega_nr_rad_s` | `.omega_nr_rad_s` | Y | f64 | Natural frequency [rad/s] |
| **Method:** | `.hub_moment_factor(omega_rad_s)` | Y | f64 | Compute reduction factor |

**PyClass Name:** `_dynbem.FlapProperties`

---

### RotorDefinition

| Rust Field | Python Property | Match? | Type | Notes |
|------------|-----------------|--------|------|-------|
| `blade` | `.blade` | Y | BladeGeometry | Nested |
| `airfoil` | `.airfoil` | Y | LinearPolarParameters | Nested |
| `control` | `.control` | Y | Option<ControlProperties> | Nested or None |
| `pitch_actuation` | `.servoflap` | **N** | Option<ServoFlapActuation> | **Renamed!** Returns None or ServoFlapActuation |
| `flap` | `.flap` | Y | Option<FlapProperties> | Nested or None |
| `name` | `.name` | Y | &str | Read-only |
| `description` | `.description` | Y | &str | Read-only |

**PyClass Name:** `_dynbem.RotorDefinition`

**Field Rename:** Rust field `pitch_actuation: PitchActuation` (enum) is exposed as `.servoflap: Option<ServoFlapActuation>` in Python. The `.servoflap` getter returns `Some(ServoFlapActuation)` if mode is `ServoFlap(...)`, else `None`.

**Python Wrapper:** Also in `dynbem.rotor_definition.RotorDefinition` (subclass with YAML loading, metadata fields)

---

## Rotor State Types

### QuasiStaticRotorState

| Rust Field | Python Property | Match? | Notes |
|------------|-----------------|--------|-------|
| *(no fields)* | `.to_array()` | Y | Returns empty numpy array |
| | `.from_array(arr)` | Y | Validates arr.len() == 0 |

**PyClass Name:** `_dynbem.QuasiStaticRotorState`  
**Behavior:** Zero inflow DOFs (quasi-static). Serializes to/from empty array.

---

### PittPetersRotorState

| Rust Field | Python Property | Match? | Type | Notes |
|------------|-----------------|--------|------|-------|
| `lambda_0` | `.lambda_0` (r/w) | Y | f64 | Uniform inflow harmonic |
| `lambda_c` | `.lambda_c` (r/w) | Y | f64 | Cyclic cosine harmonic |
| `lambda_s` | `.lambda_s` (r/w) | Y | f64 | Cyclic sine harmonic |
| **Methods:** | `.to_array()` | Y | numpy (3,) | `[lambda_0, lambda_c, lambda_s]` |
| | `.from_array(arr)` | Y | Self | Validates arr.len() == 3 |

**PyClass Name:** `_dynbem.PittPetersRotorState`

---

### OyeRotorState

| Rust Field | Python Property | Match? | Type | Notes |
|------------|-----------------|--------|------|-------|
| `n_elements` | *(internal)* | Y | usize | Not exposed |
| `W_int` | `.W_int` (get) | Y | numpy (n,) | Intermediate filter state [per annulus] |
| `W` | `.W` (get) | Y | numpy (n,) | Current inflow state [per annulus] |
| **Methods:** | `.to_array()` | Y | numpy (2n,) | `[W_int[0..n], W[0..n]]` |
| | `.from_array(arr)` | Y | Self | Validates arr.len() == 2*n_elements |
| | `.zeros(n_elements)` | Y | cls method | Factory |

**PyClass Name:** `_dynbem.OyeRotorState`

**Note:** Properties are read-only; state is modified via `.set_inflow()` or replacement.

---

### VpmRotorState

| Rust Field | Notes |
|------------|-------|
| *(TBD - VPM-specific)* | Not fully detailed in this analysis; consult `dynbem_rs/src/vpm/mod.rs` |

**PyClass Name:** `_dynbem.VpmRotorState`

---

## Rotor Model Classes (Factory-Generated)

All models expose the following interface:

```python
model.defn -> RotorDefinition
model.n_psi_elements -> usize
model.initial_rotor_state() -> (QuasiStaticRotorState | PittPetersRotorState | OyeRotorState)
model.compute_forces(inputs, state) -> (AeroResult, state)
model.step(inputs, state, dt, integration_method=None) -> (AeroResult, state)
model.inflow_taus(inputs, state) -> numpy (n_dof,)
```

### QuasiStaticBEM (Quasi-Static BEM Level 1)

**PyClass Names:**
- `_dynbem._QuasiStaticBEMLinear`
- `_dynbem._QuasiStaticBEMTabulated`

**Factory:** `dynbem.QuasiStaticBEM(defn, polar=None, n_psi_elements=36)`

---

### PittPetersModel (Pitt-Peters Level 2)

**PyClass Names:**
- `_dynbem._PittPetersModelLinear`
- `_dynbem._PittPetersModelTabulated`

**Factory:** `dynbem.PittPetersModel(defn, polar=None, n_psi_elements=36)`

---

### OyeBEMModel (Oye Level 2)

**PyClass Names:**
- `_dynbem._OyeBEMModelLinear`
- `_dynbem._OyeBEMModelTabulated`

**Factory:** `dynbem.OyeBEMModel(defn, polar=None, n_psi_elements=36, coupling_k=0.6)`

**Additional Field:**
| Rust Field | Python Property | Match? |
|------------|-----------------|--------|
| `coupling_k` | `.coupling_k` | Y |

---

### VpmRotor (Free-Wake VPM Level 3)

**PyClass Names:**
- `_dynbem._VpmRotorLinear`
- `_dynbem._VpmRotorTabulated`

**Factory:** `dynbem.VpmRotor(defn, polar=None, **config)`

**Config Parameters:** `max_particles`, `sigma`, `relax`, `nonlinear_lifting_line`, `tip_clustering`, `local_core`, `barnes_hut`, `bh_theta`, `bh_min_particles`

---

## Python-Side Wrappers (dynbem/python/dynbem/)

### rotor_definition.py

Additional Python-only classes wrapping the Rust types with metadata fields:

| Class | Rust Backing | Extra Fields | Notes |
|-------|--------------|--------------|-------|
| `LinearPolarParameters` | `_RustLinearPolarParameters` | `polar_csv`, `Re_design` | YAML metadata |
| `BladeGeometry` | `_RustBladeGeometry` | *(none)* | Compat subclass |
| `ControlProperties` | `_RustControlProperties` | *(none)* | Direct re-export |
| `ServoFlapGeometry` | `_RustServoFlapGeometry` | *(none)* | Direct re-export |
| `ServoFlapActuation` | `_RustServoFlapActuation` | *(none)* | Direct re-export |
| `FlapProperties` | `_RustFlapProperties` | *(none)* | Direct re-export |
| `RotorDefinition` | `_RustRotorDefinition` | `inertia`, `kaman_flap`, `autorotation`, `validation_issues` | YAML metadata |

**Key:** Python classes hold a `_rust` attribute containing the lean Rust wrapper.

### rotor_state.py

Virtual base class registering concrete state types:

```python
RotorState.register(QuasiStaticRotorState)
RotorState.register(PittPetersRotorState)
RotorState.register(OyeRotorState)
```

---

## Summary of Field Name Mismatches

| Type | Rust Field | Python Property | Impact | Severity |
|------|------------|-----------------|--------|----------|
| `AeroResult` | `M_hub_world` | `M_orbital` | User-facing alias for moment | **MEDIUM** |
| `RotorDefinition` | `pitch_actuation: PitchActuation` | `servoflap: Option<ServoFlapActuation>` | Extracted from enum variant | **LOW** (intentional) |
| `RotorInputs` | *(N/A)* | `t: f64` | Python-only timestamp field | **LOW** (enhancement) |

---

## Key Constants & Enums

### IntegrationMethod (aero_model.rs)

Rust enum mapped to Python string parser:

```python
# Rust
ExplicitEuler
SemiImplicitEuler
ExponentialRelaxation

# Python accepts
"explicit", "explicit_euler"
"semi_implicit", "semi-implicit", "implicit"
"exponential", "exponential_relaxation", "exp"
```

**Default:** `"semi_implicit"`

---

## Usage Patterns from README & Examples

### Creating Inputs

```python
from dynbem import RotorInputs, Mat3, Vec3
import numpy as np

inputs = RotorInputs(
    collective_rad=0.1,
    tilt_lon=0.05,
    tilt_lat=0.0,
    R_hub=np.eye(3),  # 3x3 matrix
    v_hub_world=np.array([10.0, 0.0, -2.0]),
    wind_world=np.array([5.0, 1.0, 0.0]),
    omega_rad_s=30.0,
    t=0.0,
    rho_kg_m3=1.225,
)
```

### Reading Results

```python
aero_result, new_state = model.compute_forces(inputs, state)

F_world = aero_result.F_world  # numpy (3,)
M_moment = aero_result.M_orbital  # **ALIAS** for M_hub_world
Q_torque = aero_result.Q_spin  # scalar
M_spin = aero_result.M_spin  # numpy (3,)
```

### Pitt-Peters State Access

```python
state = model.initial_rotor_state()
print(state.lambda_0, state.lambda_c, state.lambda_s)

# Modify individually
state.set_lambda_0(0.15)

# Or serialize/deserialize
arr = state.to_array()  # [lambda_0, lambda_c, lambda_s]
new_state = PittPetersRotorState.from_array(arr)
```

---

## References

- **Rust Core:** `dynbem_rs/src/`
- **Python Bindings:** `dynbem/src/wrappers.rs`
- **Python Wrappers:** `dynbem/python/dynbem/`
- **Tests:** `tests/test_*.py`
- **AGENTS.md:** Sign conventions, NED frame
