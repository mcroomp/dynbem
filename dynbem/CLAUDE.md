# CLAUDE.md - dynbem (PyO3 glue + Python shim)

This is the **public `dynbem` Python package**: PyO3 bindings wrapping
the pure-Rust [`dynbem_rs`](../dynbem_rs/) core, plus a small Python
compat shim that gives the package a drop-in-compatible surface with
the legacy pure-Python `dynbem` (now archived as [`dynbem_old/`](../dynbem_old/)
for reference only).

**Parent doc:** [../CLAUDE.md](../CLAUDE.md) - sign conventions, NED
frame, rotor rotation direction, all the load-bearing physics that this
package binds.

## Layout

    dynbem/
    +-- Cargo.toml             pyo3 cdylib named "_dynbem", depends on ../dynbem_rs
    +-- pyproject.toml         maturin -> builds python/dynbem/_dynbem.pyd
    +-- src/
    |   +-- lib.rs             #[pymodule] registration, free pyfunctions,
    |   |                      solve_bem_element pyfunction + BEMElementResult
    |   +-- conv.rs            numpy <-> Vec3/Mat3 marshalling (only numpy site)
    |   +-- wrappers.rs        PyFoo(pub Foo) newtypes for ~15 core types
    |   +-- trim_py.rs         AeroAny dispatch, TrimResult, trim/relax pyfunctions
    +-- python/dynbem/
    |   +-- __init__.py        top-level API + Python subclasses (see below)
    |   +-- bem.py             compat shim re-exporting dynbem.bem.*
    |   +-- cyclic.py          compat shim
    |   +-- factory.py         create_aero(), build_polar(), load_tabulated_polar()
    |   +-- oye.py             compat shim
    |   +-- pitt_peters.py     compat shim
    |   +-- polar.py           compat shim + AirfoilPolar tuple alias
    |   +-- rotor_definition.py  compat shim + YAML load() + ValidationIssue
    |   +-- rotor_state.py     compat shim + RotorState ABC
    |   +-- trim.py            solve_trim_cyclic / relax_inflow (thin re-exports)
    +-- tests/                 (intentionally empty -- repo-wide tests live in /tests)
    +-- benchmarks/
        +-- bench_rust_only.py

## Drop-in compatibility with legacy dynbem

The Rust pyclass surface is small and explicit; the legacy Python
package had richer ergonomics. The Python shim in `python/dynbem/`
fills the gap so existing `import dynbem` callers keep working:

- **Dotted submodule access** (`dynbem.bem.BEMModel`, `dynbem.polar.LinearPolar`,
  ...): each shim file re-exports the relevant names from the parent
  package and adds a few helpers.
- **YAML rotor loader** `dynbem.rotor_definition.load(path)` /
  `loads(text, base_dir=None)`: thin wrappers around the
  `_dynbem.load_rotor_yaml` / `loads_rotor_yaml` pyfunctions, which
  delegate to `dynbem_rs::rotor_yaml` (serde_yaml). Rust and Python
  share a single parser + schema, so pure-Rust callers can load the
  same `rotor.yaml` files via `RotorDefinition::from_yaml_file(path)`.
- **`.validate()` methods** on BladeGeometry / AirfoilProperties /
  RotorDefinition: monkey-patched onto the Rust pyclasses in
  `rotor_definition.py`. Return `list[ValidationIssue]`. If pyo3 ever
  freezes pyclass type-dicts, convert these to free functions and
  update the tests' `b.validate()` -> `validate_blade(b)`.
- **Virtual ABCs** `dynbem.AeroBase` and `dynbem.RotorState`: declared
  in `__init__.py` with `ABC.register(...)` so `isinstance(model, AeroBase)`
  works. The Rust classes don't actually inherit from these; the ABC is
  a marker only.
- **`solve_trim_cyclic` / `relax_inflow`**: thin re-exports of the Rust
  pyfunctions; both require a pre-built `RotorInputs` as the third argument.
- **Auto polar inference**: `BEMModel(defn=...)` without a polar argument
  builds one from `defn.airfoil` automatically (including `polar_csv`
  CSV loading). Implemented as Python subclasses of the Rust pyclasses
  in `__init__.py`. This requires `subclass=True` on the three Rust
  model pyclasses in `wrappers.rs`.

## When you add a new core API

If the new core feature requires Python ergonomics that don't naturally
fit in a pyclass (e.g. file IO, polymorphic kwargs, validation
warnings, abstract base class registration):

1. Expose a minimal Rust pyfunction or pyclass in `src/`.
2. Wrap with the convenience layer in `python/dynbem/`.
3. Add the dotted-submodule re-export if relevant.
4. **Do NOT add pyo3 / numpy to dynbem_rs (the pure-Rust core).** File
   IO is now allowed but only in a dedicated module (currently just
   `dynbem_rs::rotor_yaml`); the math modules stay free of `std::fs`
   and `serde` so they remain embeddable.

## Numpy boundary discipline

All numpy use is in `conv.rs` (3 helpers) and a handful of `to_array` /
`from_array` / array-returning getters in `wrappers.rs`. If you find
yourself reaching for `PyReadonlyArray*` anywhere else in Rust, route
it through `conv.rs` instead.

## When tests fail

The repo's tests live in `e:/repos/aero/tests/` (228 tests covering
both production code paths and verification scripts) and are the
authoritative regression suite. Run with:

    .venv\Scripts\python -m pytest tests/ -q

A `dynbem`-side change that breaks them is a real regression. The
`dynbem_old/` package is read-only reference material; do not depend on
it in new code or tests.
