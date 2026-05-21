# CLAUDE.md - dynbem_rs (pure-Rust core)

Pure-Rust BEM / Pitt-Peters / Oye dynamic-inflow rotor models. **No
pyo3, no numpy, no file IO.** This crate is the math core; the Python
bindings live in the sibling [`dynbem/`](../dynbem/) glue crate.

**Parent doc:** [../CLAUDE.md](../CLAUDE.md) - sign conventions, NED
frame, rotor rotation direction, hub-frame moments, Pitt-Peters L-matrix,
Oye filter -- all load-bearing physics.

## Hard rules

1. **No `pyo3` / `numpy` / `pyyaml` / file IO imports.** If you need a
   Python-facing helper, file IO, or pickling, you're in the wrong
   crate -- add it in [`../dynbem/`](../dynbem/) instead.
2. **Public API stability matters.** The glue crate depends on every
   public field, struct, and function here. Renaming or moving things
   requires a matching edit there. Prefer additive changes.
3. **Sign conventions are NOT documented in this crate's source.**
   They live in [../CLAUDE.md](../CLAUDE.md). Refer to it; don't
   duplicate.

## Module map

    src/
    +-- lib.rs                public module declarations
    +-- aero_io.rs            Vec3, Mat3, RotorInputs, AeroResult
    +-- aero_model.rs         AeroModel trait + RotorStateExt
    +-- bem.rs                BEMModel + pub solve_bem_element + windmill solver
    +-- bem_common.rs         RadialGrid, PolarTable (shared by all models)
    +-- common.rs             numerical floors (EPS_*), VRS polynomial
    +-- cyclic.rs             swashplate -> theta_1c, theta_1s mapping
    +-- oye.rs                OyeBEMModel (annular 2-stage filter)
    +-- pitt_peters.rs        PittPetersModel (3-state L-matrix ODE)
    +-- polar.rs              LinearPolar, TabulatedPolar, PolarKind
    +-- rotor_definition.rs   Blade / Airfoil / Control / Inertia / etc.
    +-- rotor_state.rs        QuasiStatic / PittPeters / Oye state structs
    +-- trim.rs               solve_trim_cyclic<M>, relax_inflow<M> (generic)

## Adding a new aero model

1. Add `src/foo.rs` with the model struct and `impl AeroModel for FooModel`.
2. Add `FooRotorState` to `rotor_state.rs` and the `RotorStateExt` impl
   in `aero_model.rs`. **Mechanical states (omega, spin_angle) MUST be
   the last two entries in `to_vec()`** -- the trim integrator's
   omega-clamp relies on this.
3. Add `pub mod foo;` to `lib.rs`.
4. Add a wrapper newtype in [../dynbem/src/wrappers.rs](../dynbem/src/wrappers.rs)
   (mark it `subclass = true` if Python should auto-build the polar)
   plus the `AeroAny` variant in [../dynbem/src/trim_py.rs](../dynbem/src/trim_py.rs).

## Hot-path conventions

- Per-azimuth kinematics (`hub_axis`, `v_climb`, `v_edge`,
  `v_inplane_hub`) are **deliberately inline** in each model's
  `compute_forces`, not abstracted to a helper. The polar interp loop
  and the radial integration loop are the only big duplications;
  numba/LLVM don't compose closures cleanly across them.
- Plain `for` loops over `&[f64]`, no SIMD intrinsics. LLVM autovectorizes
  the if-converted bodies. See user-memory `feedback-rust-autovectorize-first`.
- Vec3/Mat3 are `Copy` newtypes around plain f64 arrays. Operators
  (`R * v`, `a - b`, `s * v`) lower to the same scalar FMA chains as
  hand-rolled index arithmetic.

## Numerical floors

All in `common.rs`:
- `EPS_DENOM` (1e-9) - generic division/ratio safety
- `EPS_OMEGA_R` (1e-6) - rotor-not-spinning threshold
- `MIN_LOSS_FACTOR` (1e-4) - Prandtl tip+hub loss floor
- `V_T_HOVER_FLOOR_FRAC` (1e-2) - mass-flow floor at hover/zero-thrust
- `VRS_DESCENT_THRESHOLD` (1e-3) - VRS detection guard against hover
  chattering
- `MU_T_FLOOR` (0.05) - Pitt-Peters L-matrix denominator floor

Empirical, tuned to make hover / climb / descent / VRS / autorotation
all stable in one code path. Don't change them without running the
full `tests/` suite (`pytest tests/ -q`).
