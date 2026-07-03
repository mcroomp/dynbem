# Quick Reference: Struct Field Names

## Most-Used Types

### RotorInputs (Control inputs)
```
collective_rad       : f64 [rad]      <- r/w ->  .collective_rad
tilt_lon             : f64 [rad]      <- r/w ->  .tilt_lon
tilt_lat             : f64 [rad]      <- r/w ->  .tilt_lat
R_hub                : Mat3 (3x3)     <- r/o ->  .R_hub (numpy 3x3)
v_hub_world          : Vec3           <- r/o ->  .v_hub_world (numpy 3,)
wind_world           : Vec3           <- r/o ->  .wind_world (numpy 3,)
rho_kg_m3            : f64 [kg/m³]    <- r/w ->  .rho_kg_m3
omega_rad_s          : f64 [rad/s]    <- r/w ->  .omega_rad_s
(Python only)  t     : f64 [s]        <- r/w ->  .t
```

### AeroResult (Aerodynamic output)
```
F_world              : Vec3           -> .F_world (numpy 3,)
M_hub_world          : Vec3           -> .M_orbital ⚠️ RENAMED (numpy 3,)
Q_spin               : f64 [N⋅m]      -> .Q_spin
M_spin               : Vec3           -> .M_spin (numpy 3,)
```

### PittPetersRotorState (Inflow state)
```
lambda_0             : f64            <- r/w ->  .lambda_0
lambda_c             : f64            <- r/w ->  .lambda_c
lambda_s             : f64            <- r/w ->  .lambda_s
```

### OyeRotorState (Per-annulus inflow)
```
W_int: Vec<f64>      : Array[n]       <- r/o ->  .W_int (numpy n,)
W: Vec<f64>          : Array[n]       <- r/o ->  .W (numpy n,)
```

### BladeGeometry (Blade specs)
```
n_blades             : usize          -> .n_blades
radius_m             : f64 [m]        -> .radius_m
root_cutout_m        : f64 [m]        -> .root_cutout_m
chord_m              : f64 [m]        -> .chord_m
twist_deg            : f64 [deg]      -> .twist_deg
n_elements           : usize          -> .n_elements
tip_loss             : bool           -> .tip_loss
r_stations_m         : Vec<f64>       (internal, used for interpolation)
chord_stations_m     : Vec<f64>       (internal, used for interpolation)
twist_stations_deg   : Vec<f64>       (internal, used for interpolation)
```

### LinearPolarParameters (Airfoil coefficients)
```
CL0                  : f64            -> .CL0
CL_alpha_per_rad     : f64 [rad⁻¹]    -> .CL_alpha_per_rad
CD0                  : f64            -> .CD0
alpha_stall_deg      : f64 [deg]      -> .alpha_stall_deg
```

### ControlProperties (Swashplate)
```
swashplate_pitch_gain_rad  : f64 [rad]     -> .swashplate_pitch_gain_rad
swashplate_phase_deg       : Option<f64>   -> .swashplate_phase_deg
```

### FlapProperties (Blade flapping)
```
I_blade_flap_kgm2    : f64 [kg⋅m²]    -> .I_blade_flap_kgm2
omega_nr_rad_s       : f64 [rad/s]    -> .omega_nr_rad_s
```

### ServoFlapActuation (Servo-flap)
```
I_theta_kgm2         : f64 [kg⋅m²]    -> .I_theta_kgm2
damper_Nms_per_rad   : f64 [N⋅m⋅s/rad] -> .damper_Nms_per_rad
ac_offset_m          : f64 [m]        -> .ac_offset_m
control_stiffness_Nm_per_rad : f64    -> .control_stiffness_Nm_per_rad
flap                 : ServoFlapGeometry -> .flap
```

### ServoFlapGeometry
```
C_M_delta_per_rad    : f64 [rad⁻¹]    -> .C_M_delta_per_rad
r_inner_m            : f64 [m]        -> .r_inner_m
r_outer_m            : f64 [m]        -> .r_outer_m
```

### RotorDefinition
```
blade                : BladeGeometry  -> .blade
airfoil              : LinearPolarParameters -> .airfoil
control              : Option<ControlProperties> -> .control
pitch_actuation: PitchActuation -> .servoflap ⚠️ ENUM EXTRACT
flap                 : Option<FlapProperties> -> .flap
name                 : String         -> .name (read-only)
description          : String         -> .description (read-only)
```

---

## Abbreviations

| r/w | read/write (mutable property/setter) |
| r/o | read-only (getter only) |
| ⚠️  | name mismatch or special handling |

---

## PyO3 Class Names

All Rust types exposed as `_dynbem.<ClassName>` (module `dynbem._dynbem`):

- `LinearPolar`
- `TabulatedPolar`
- `BladeGeometry`
- `LinearPolarParameters`
- `ControlProperties`
- `ServoFlapGeometry`
- `ServoFlapActuation`
- `FlapProperties`
- `RotorDefinition`
- `RotorInputs`
- `AeroResult`
- `QuasiStaticRotorState`
- `PittPetersRotorState`
- `OyeRotorState`
- `VpmRotorState`
- `_QuasiStaticBEMLinear`
- `_QuasiStaticBEMTabulated`
- `_PittPetersModelLinear`
- `_PittPetersModelTabulated`
- `_OyeBEMModelLinear`
- `_OyeBEMModelTabulated`
- `_VpmRotorLinear`
- `_VpmRotorTabulated`

---

## Critical Mismatches

1. **AeroResult.M_hub_world → M_orbital**
   - Rust internal name is `M_hub_world`
   - Python property is named `M_orbital`
   - This is intentional for API consistency
   - Code must use `.M_orbital` when reading Python results

2. **RotorDefinition.pitch_actuation → servoflap**
   - Rust has enum `PitchActuation { DirectMechanical, ServoFlap(...) }`
   - Python property `.servoflap` extracts the `ServoFlapActuation`
   - Returns `ServoFlapActuation | None` (not the enum itself)

3. **RotorInputs.t (Python-only)**
   - Not stored in Rust `RotorInputs` struct
   - Stored separately in `PyRotorInputs.t`
   - Used for time tracking in Python integrators

4. **LinearPolarParameters.alpha_stall_deg vs LinearPolar.alpha_stall_rad**
   - `LinearPolarParameters`: degrees (.alpha_stall_deg)
   - `LinearPolar`: radians (.alpha_stall_rad)
   - Conversion happens in `LinearPolar::new()` or via factory

---

## Array Serialization

State types support `.to_array()` and `.from_array()` for checkpoint/restore:

| State Type | Array Shape | Contents |
|------------|-------------|----------|
| `QuasiStaticRotorState` | `()` (empty) | No states |
| `PittPetersRotorState` | `(3,)` | `[lambda_0, lambda_c, lambda_s]` |
| `OyeRotorState` | `(2*n,)` | `[W_int[0..n], W[0..n]]` |

---

## Common Constructor Signatures

```python
# Core geometry
BladeGeometry(
    n_blades, radius_m, root_cutout_m, chord_m,
    twist_deg=0.0, n_elements=20,
    r_stations_m=None, chord_stations_m=None, twist_stations_deg=None,
    tip_loss=True
)

# Airfoil
LinearPolarParameters(CL0, CL_alpha_per_rad, CD0, alpha_stall_deg)

# Control
ControlProperties(swashplate_pitch_gain_rad, swashplate_phase_deg=None)

# Flapping
FlapProperties(I_blade_flap_kgm2, omega_nr_rad_s=0.0)

# Servo-flap
ServoFlapActuation(
    I_theta_kgm2, damper_Nms_per_rad, flap,
    ac_offset_m=0.0, control_stiffness_Nm_per_rad=0.0
)

# Rotor definition
RotorDefinition(
    blade, airfoil, control, name, description,
    servoflap=None, flap=None
)

# Inputs
RotorInputs(
    collective_rad, tilt_lon, tilt_lat,
    R_hub, v_hub_world, wind_world,
    omega_rad_s, t, rho_kg_m3
)

# States
QuasiStaticRotorState()
PittPetersRotorState(lambda_0, lambda_c, lambda_s)
OyeRotorState(W_int, W)  # numpy arrays

# Models
QuasiStaticBEM(defn, polar=None, n_psi_elements=36)
PittPetersModel(defn, polar=None, n_psi_elements=36)
OyeBEMModel(defn, polar=None, n_psi_elements=36, coupling_k=0.6)
VpmRotor(defn, polar=None, **config)
```

---

## Integration Pattern

```python
# Typical time-stepping loop
state = model.initial_rotor_state()
for t in time_points:
    inputs.t = t
    aero_result, state = model.step(inputs, state, dt, "semi_implicit")
    
    # Read outputs
    thrust = aero_result.F_world
    moment = aero_result.M_orbital  # ← Use M_orbital, not M_hub_world
    torque = aero_result.Q_spin
```
