# dynbem (PyO3 glue + Python package)

This is the **public `dynbem` Python package**: PyO3 bindings wrapping
the pure-Rust [`dynbem_rs`](../dynbem_rs/) core, plus a Python layer that
provides YAML loading, compat shims, and ergonomic helpers.

For sign conventions, NED frame, rotor rotation direction, hub-frame moments,
and Pitt-Peters L-matrix see [`../AGENTS.md`](../AGENTS.md).

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
        +-- __init__.py        top-level API + Python subclasses (see below)
        +-- bem.py             compat shim re-exporting dynbem.bem.*
        +-- cyclic.py          compat shim
        +-- factory.py         create_aero(), build_polar(), load_tabulated_polar()
        +-- mechanical.py      omega_derivative(), euler_step_omega() -- caller-owned spin ODE
        +-- oye.py             compat shim
        +-- pitt_peters.py     compat shim
        +-- polar.py           compat shim + AirfoilPolar tuple alias
        +-- rotor_definition.py  YAML load() / loads() in pure Python (yaml.safe_load)
        +-- rotor_state.py     compat shim + RotorState ABC
        +-- trim.py            solve_trim_cyclic / relax_inflow (thin re-exports)

## Python layer responsibilities

**YAML loading** is pure Python (`yaml.safe_load` in `rotor_definition.py`).
The resulting Python `RotorDefinition` holds a lean Rust `_dynbem.RotorDefinition`
wrapper (math fields only) as its `._rust` attribute. Call
`dynbem.rotor_definition.load(path)` or `loads(text, base_dir=None)`.

**Dotted submodule access** (`dynbem.bem.BEMModel`, `dynbem.polar.LinearPolar`,
etc.): each shim file re-exports the relevant names from the parent package
and adds a few helpers so existing `import dynbem` callers keep working.

**Auto polar inference**: `BEMModel(defn=...)` without a polar argument
builds one from `defn.airfoil` automatically (including `polar_csv` CSV
loading). Implemented as Python subclasses of the Rust pyclasses in
`__init__.py`. Requires `subclass=True` on the three Rust model pyclasses
in `wrappers.rs`.

**Mechanical ODE**: `mechanical.py` provides `omega_derivative(Q_aero,
motor_torque_Nm, I_ode_kgm2)` and `euler_step_omega`. The aero models are
pure aerodynamic -- they take `omega_rad_s` via `RotorInputs` and return
loads only. The caller owns the spin ODE.

**Virtual ABCs** `dynbem.AeroBase` and `dynbem.RotorState` are declared in
`__init__.py` with `ABC.register(...)` so `isinstance(model, AeroBase)`
works. The Rust classes do not actually inherit from these; the ABC is a
marker only.

**`solve_trim_cyclic` / `relax_inflow`**: thin re-exports of the Rust
pyfunctions; both require a pre-built `RotorInputs` as the third argument.

## Hard rules

1. **All numpy use is confined to `conv.rs`** (3 helpers) and a handful of
   `to_array` / `from_array` / array-returning getters in `wrappers.rs`.
   Do not reach for `PyReadonlyArray*` elsewhere in Rust -- route it through
   `conv.rs`.
2. **Do NOT add pyo3 / numpy to dynbem_rs.** The math core stays free of
   `std::fs` and `serde` (except the YAML fields in `rotor_definition.rs`
   which are populated by the Python loader, not by serde).
3. **Expose a minimal Rust pyclass/pyfunction first, then wrap in Python.**
   If ergonomics don't fit naturally in a pyclass (file IO, polymorphic
   kwargs, validation warnings, ABC registration), add the convenience in
   `python/dynbem/`.

## When you add a new core API

1. Expose a minimal Rust pyfunction or pyclass in `src/`.
2. Wrap with the convenience layer in `python/dynbem/`.
3. Add the dotted-submodule re-export if relevant.
4. Update `create_aero()` in `factory.py` with a stable model-name string.

## Tests

The repo-wide test suite lives in `../tests/` and covers both production
code paths and validation scripts. Run with:

    uv run pytest tests/ -q

A `dynbem`-side change that breaks them is a real regression.
