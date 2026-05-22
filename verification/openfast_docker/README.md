# OpenFAST AeroDyn cross-check (Docker)

Drives OpenFAST's `aerodyn_driver` standalone in a container so dynbem
can be compared against the NREL AeroDyn reference implementation
(specifically DBEMT_Mod=1, which is the same Oye 2-stage filter algorithm
implemented by `dynbem.OyeBEMModel`). Mirrors the
[ccblade_docker](../ccblade_docker/) layout.

## Status (initial spike)

**Infrastructure**: working. The Dockerfile, input-file generators, and
output parser all run end-to-end on Caradonna-Tung hover input. See
[Replicating the Caradonna-Tung run](#replicating-the-caradonna-tung-run)
below.

**Hover cross-check**: AeroDyn's BEM solver iterates on the wind-turbine
induction factor `a = v_axial / V_inf`. At `V_inf = 0` (helicopter hover)
that ratio is ill-defined, so AeroDyn effectively does no induced-inflow
correction and the blade sees the geometric pitch directly. Thrust is
over-predicted by 2.5x to 4x relative to the paper.

This is a known limitation of AeroDyn (it is designed for wind turbines,
which always have non-zero `V_inf`). It is *not* a bug in our driver
scaffolding. dynbem iterates on the total inflow ratio
`lambda_r = v_a/(Omega*R)` instead, which is well-defined at `V_inf=0` --
that is precisely why dynbem exists as a separate codebase.

**What does work** with this setup:

- **Axial climb / descent** at `V_inf` large enough for AeroDyn's
  induction factor to make sense (roughly `V_inf >= 2 * v_h`, where
  `v_h = sqrt(CT/2) * Omega*R` is the hover induced velocity).
- **Forward flight** at non-zero advance ratio.

Both are useful comparison regimes for `dynbem.OyeBEMModel` and
`dynbem.PittPetersModel`, which themselves handle hover via their own
solver path.

## What this does not (yet) cover

- Pure hover validation: needs a different solver convention than
  AeroDyn provides. Best handled by direct CT comparison against the
  Caradonna-Tung paper (which `tests/test_caradonna_tung*.py` already
  does for dynbem).
- Vortex-Ring State: AeroDyn does not model VRS at all.
- Forward-flight cyclic + hub-moment comparison: AeroDyn supports it
  but the driver `.dvr` setup is more elaborate; deferred.

## Files

    verification/openfast_docker/
      Dockerfile               FROM octue/openfast:3.5.3 (modern OpenFAST
                               build with /openfast/install/bin/aerodyn_driver),
                               adds python3 + pyyaml + numpy.
      compose.yml              docker compose service definition.
      run.py                   Reads YAML; writes a .dvr + AeroDyn .dat + blade
                               .dat + airfoil .dat per case under a workdir;
                               invokes aerodyn_driver; parses the .out
                               text output and averages the last revolution
                               into a per-rotor CSV under verification/data/.
      inputs/
        caradonna_tung_hover.yaml   Caradonna-Tung 2-blade NACA 0012 rotor,
                                    3 pitch settings at 1250 rpm. Useful
                                    for confirming infrastructure; thrust
                                    values do not match the paper because
                                    of the hover-BEM issue described above.
      README.md                This file.

## Replicating the Caradonna-Tung run

From the repo root:

    cd verification/openfast_docker
    docker compose build
    MSYS_NO_PATHCONV=1 docker compose run --rm openfast \
        --config /in/caradonna_tung_hover.yaml

(`MSYS_NO_PATHCONV=1` is needed on Windows + Git Bash to stop the shell
from rewriting `/in` to a Windows path. On WSL or Linux it can be
omitted.)

Output goes to `verification/data/caradonna_tung_hover_openfast.csv`
with one row per pitch case. Columns include the AeroDyn channels
(`RtAeroFxh`, `RtAeroMxh`, `RtAeroPwr`, ...) plus the input
`(U_wind_ms, Omega_rpm, pitch_deg)`.

## YAML schema

```yaml
name: <free-form, used to name the output CSV>

rotor:
  n_blades: int
  R_hub_m: float            # blade inboard cutout (m, from rotor centre)
  R_tip_m: float
  precone_deg: float        # optional, defaults 0

air:
  rho: float                # optional, defaults 1.225 kg/m^3
  mu:  float                # optional, defaults 1.81206e-5 Pa*s
  a:   float                # speed of sound; optional, defaults 340 m/s

blade_stations:             # sorted by r_m, in m from rotor centre
  - {r_m: float, chord_m: float, twist_deg: float}
  - ...

airfoil_polar:              # single Re, single table; all blade stations
                            # share this polar (refine to per-station
                            # polars later if needed)
  Re: float
  points:
    - {alpha_deg: float, cl: float, cd: float}
    - ...

wake:
  WakeMod:    int           # 1=BEMT (steady), 2=DBEMT (dynamic), 3=OLAF
  DBEMT_Mod:  int           # 1=const tau1 (matches dynbem.oye Mod=1),
                            # 2=time-dependent tau1, 3=const tau1 continuous
  tau1_const: float         # seconds; only used when DBEMT_Mod in {1,3}
  SkewMod:    int           # 1=uncoupled (axial), 2=Pitt-Peters skew
  TipLoss:    bool          # Prandtl tip loss
  HubLoss:    bool          # Prandtl hub loss

operating_points:
  - {U_wind_ms: float, Omega_rpm: float, pitch_deg: float, settle_s: float}
  - ...                     # settle_s is how long to run before averaging
                            # the last rotor revolution into the CSV row.
```

## Next steps

1. Add an axial climb YAML (e.g. `caradonna_tung_axial_climb.yaml` with
   `U_wind_ms ~= 10`) where AeroDyn's BEM does converge cleanly, and
   compare against `dynbem.OyeBEMModel` driven at the same operating
   point.
2. Add a forward-flight YAML for Pitt-Peters skew-correction
   comparison (AeroDyn `SkewMod=2`) -- this is the closest direct
   analogue of `dynbem.PittPetersModel`'s wake-skew off-diagonal term.
3. Decide whether to also wire DBEMT_Mod=2 (time-dependent tau1) which
   has no direct dynbem analogue but is the OpenFAST default.

## Caveats

- The container image (`octue/openfast:3.5.3`) is built single-precision;
  DBEMT_Init refuses to start with `dt < sqrt(eps) ~= 3.5e-4 s`. The
  driver enforces `dt >= 1 ms`. If you want to compare against a
  double-precision build, change the base image.
- Channel coverage from the standalone AeroDyn driver is narrower than
  full OpenFAST; in particular, blade-element output (per-station Cn,
  Ct, induced inflow) requires turning on the optional "all-blade-node"
  output section in the AeroDyn primary file. Not wired here yet -- only
  rotor-aggregate channels are pulled into the CSV.
- AeroDyn always reads three `ADBlFile` entries (one per max-supported
  blade) regardless of `NumBlades`; entries past `NumBl` are ignored.
  The driver always emits three.
