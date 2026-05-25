// Python-facing trim solver: AeroAny dispatch enum + the two pyfunctions.
// The actual trim math is in dynbem_rs::trim, generic over AeroModel.

use crate::wrappers::{
    PyOyeBEMModelLinear, PyOyeBEMModelTabulated, PyOyeRotorState, PyPittPetersModelLinear,
    PyPittPetersModelTabulated, PyPittPetersRotorState, PyQuasiStaticBEMLinear,
    PyQuasiStaticBEMTabulated, PyQuasiStaticRotorState, PyRotorInputs,
};
use dynbem_rs::polar::{LinearPolar, TabulatedPolar};
use dynbem_rs::trim::{relax_inflow, solve_trim_cyclic};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

enum AeroAny {
    QuasiStaticBEMLinear(dynbem_rs::quasi_static_bem::QuasiStaticBEM<LinearPolar>),
    QuasiStaticBEMTabulated(dynbem_rs::quasi_static_bem::QuasiStaticBEM<TabulatedPolar>),
    PittPetersLinear(dynbem_rs::pitt_peters::PittPetersModel<LinearPolar>),
    PittPetersTabulated(dynbem_rs::pitt_peters::PittPetersModel<TabulatedPolar>),
    OyeLinear(dynbem_rs::oye::OyeBEMModel<LinearPolar>),
    OyeTabulated(dynbem_rs::oye::OyeBEMModel<TabulatedPolar>),
}

impl AeroAny {
    fn from_py(obj: &Bound<'_, PyAny>) -> PyResult<Self> {
        if let Ok(m) = obj.extract::<PyQuasiStaticBEMLinear>() {
            return Ok(AeroAny::QuasiStaticBEMLinear(*m.0));
        }
        if let Ok(m) = obj.extract::<PyQuasiStaticBEMTabulated>() {
            return Ok(AeroAny::QuasiStaticBEMTabulated(*m.0));
        }
        if let Ok(m) = obj.extract::<PyPittPetersModelLinear>() {
            return Ok(AeroAny::PittPetersLinear(*m.0));
        }
        if let Ok(m) = obj.extract::<PyPittPetersModelTabulated>() {
            return Ok(AeroAny::PittPetersTabulated(*m.0));
        }
        if let Ok(m) = obj.extract::<PyOyeBEMModelLinear>() {
            return Ok(AeroAny::OyeLinear(*m.0));
        }
        if let Ok(m) = obj.extract::<PyOyeBEMModelTabulated>() {
            return Ok(AeroAny::OyeTabulated(*m.0));
        }
        Err(PyValueError::new_err(
            "aero must be a QuasiStaticBEM, PittPetersModel, or OyeBEMModel instance",
        ))
    }
}

enum FinalState {
    QuasiStatic(dynbem_rs::rotor_state::QuasiStaticRotorState),
    PittPeters(dynbem_rs::rotor_state::PittPetersRotorState),
    Oye(dynbem_rs::rotor_state::OyeRotorState),
}

#[pyclass(name = "TrimResult", module = "dynbem._dynbem")]
#[allow(non_snake_case)]
pub struct PyTrimResult {
    #[pyo3(get)]
    pub tilt_lon: f64,
    #[pyo3(get)]
    pub tilt_lat: f64,
    #[pyo3(get)]
    #[allow(non_snake_case)]
    pub Mx_residual: f64,
    #[pyo3(get)]
    #[allow(non_snake_case)]
    pub My_residual: f64,
    #[pyo3(get)]
    pub iterations: usize,
    #[pyo3(get)]
    pub converged: bool,
    final_state: FinalState,
}

#[pymethods]
impl PyTrimResult {
    #[getter]
    fn final_state<'py>(&self, py: Python<'py>) -> PyObject {
        match &self.final_state {
            FinalState::QuasiStatic(s) => PyQuasiStaticRotorState(s.clone()).into_py(py),
            FinalState::PittPeters(s) => PyPittPetersRotorState(s.clone()).into_py(py),
            FinalState::Oye(s) => PyOyeRotorState(s.clone()).into_py(py),
        }
    }
}

#[pyfunction]
#[pyo3(signature = (
    aero, state, base_inputs, *,
    target_moment = (0.0, 0.0),
    tilt_lon_init = 0.0, tilt_lat_init = 0.0,
    tilt_min = -0.2617993877991494,  // -15 deg
    tilt_max =  0.2617993877991494,  // +15 deg
    tolerance_Nm = 0.02, max_iterations = 50,
    probe_rad = 0.008726646259971648,  // 0.5 deg
    dt_relax = 0.005, n_inflow_relax = 100,
    n_settle = 0,
))]
#[allow(clippy::too_many_arguments)]
#[allow(non_snake_case)]
pub fn solve_trim_cyclic_py(
    aero: &Bound<'_, PyAny>,
    state: &Bound<'_, PyAny>,
    base_inputs: &PyRotorInputs,
    target_moment: (f64, f64),
    tilt_lon_init: f64,
    tilt_lat_init: f64,
    tilt_min: f64,
    tilt_max: f64,
    tolerance_Nm: f64,
    max_iterations: usize,
    probe_rad: f64,
    dt_relax: f64,
    n_inflow_relax: usize,
    n_settle: usize,
) -> PyResult<PyTrimResult> {
    let model = AeroAny::from_py(aero)?;
    macro_rules! do_trim {
        ($m:expr, $s_ty:ty, $fs_variant:ident) => {{
            let s = state.extract::<$s_ty>()?.0;
            let out = solve_trim_cyclic(
                &$m,
                s,
                &base_inputs.0,
                target_moment.0,
                target_moment.1,
                tilt_lon_init,
                tilt_lat_init,
                tilt_min,
                tilt_max,
                tolerance_Nm,
                max_iterations,
                probe_rad,
                dt_relax,
                n_inflow_relax,
                n_settle,
            );
            Ok(PyTrimResult {
                tilt_lon: out.tilt_lon,
                tilt_lat: out.tilt_lat,
                Mx_residual: out.mx_residual,
                My_residual: out.my_residual,
                iterations: out.iterations,
                converged: out.converged,
                final_state: FinalState::$fs_variant(out.final_state),
            })
        }};
    }
    match model {
        AeroAny::QuasiStaticBEMLinear(m) => do_trim!(m, PyQuasiStaticRotorState, QuasiStatic),
        AeroAny::QuasiStaticBEMTabulated(m) => do_trim!(m, PyQuasiStaticRotorState, QuasiStatic),
        AeroAny::PittPetersLinear(m) => do_trim!(m, PyPittPetersRotorState, PittPeters),
        AeroAny::PittPetersTabulated(m) => do_trim!(m, PyPittPetersRotorState, PittPeters),
        AeroAny::OyeLinear(m) => do_trim!(m, PyOyeRotorState, Oye),
        AeroAny::OyeTabulated(m) => do_trim!(m, PyOyeRotorState, Oye),
    }
}

#[pyfunction]
#[pyo3(signature = (aero, state, inputs, n_steps = 200, dt = 0.005))]
pub fn relax_inflow_py(
    py: Python<'_>,
    aero: &Bound<'_, PyAny>,
    state: &Bound<'_, PyAny>,
    inputs: &PyRotorInputs,
    n_steps: usize,
    dt: f64,
) -> PyResult<PyObject> {
    let model = AeroAny::from_py(aero)?;
    macro_rules! do_relax {
        ($m:expr, $s_ty:ty, $out_wrapper:expr) => {{
            let s = state.extract::<$s_ty>()?.0;
            let out = relax_inflow(&$m, s, &inputs.0, n_steps, dt);
            Ok($out_wrapper(out).into_py(py))
        }};
    }
    match model {
        AeroAny::QuasiStaticBEMLinear(m) => {
            do_relax!(m, PyQuasiStaticRotorState, PyQuasiStaticRotorState)
        }
        AeroAny::QuasiStaticBEMTabulated(m) => {
            do_relax!(m, PyQuasiStaticRotorState, PyQuasiStaticRotorState)
        }
        AeroAny::PittPetersLinear(m) => {
            do_relax!(m, PyPittPetersRotorState, PyPittPetersRotorState)
        }
        AeroAny::PittPetersTabulated(m) => {
            do_relax!(m, PyPittPetersRotorState, PyPittPetersRotorState)
        }
        AeroAny::OyeLinear(m) => do_relax!(m, PyOyeRotorState, PyOyeRotorState),
        AeroAny::OyeTabulated(m) => do_relax!(m, PyOyeRotorState, PyOyeRotorState),
    }
}
