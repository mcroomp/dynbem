// PyO3 wrapper newtypes around dynbem-core types. Each one is a tuple-struct
// holding the core value; #[pymethods] forwards to the core API and handles
// numpy marshalling at the boundary.

use crate::conv::{mat3_to_py, read_mat3, read_vec3, vec3_to_py};
use dynbem_rs as core_;
use dynbem_rs::polar::Polar as _;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::PyTypeInfo;

// ===========================================================================
// Polars
// ===========================================================================

#[pyclass(name = "LinearPolar", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyLinearPolar(pub core_::polar::LinearPolar);

#[pymethods]
impl PyLinearPolar {
    #[new]
    #[pyo3(signature = (CL0, CL_alpha_per_rad, CD0, alpha_stall_rad))]
    #[allow(non_snake_case)]
    fn new(CL0: f64, CL_alpha_per_rad: f64, CD0: f64, alpha_stall_rad: f64) -> Self {
        PyLinearPolar(core_::polar::LinearPolar::new(
            CL0,
            CL_alpha_per_rad,
            CD0,
            alpha_stall_rad,
        ))
    }

    #[getter]
    #[allow(non_snake_case)]
    fn CL0(&self) -> f64 {
        self.0.CL0
    }
    #[getter]
    #[allow(non_snake_case)]
    fn CL_alpha_per_rad(&self) -> f64 {
        self.0.CL_alpha_per_rad
    }
    #[getter]
    #[allow(non_snake_case)]
    fn CD0(&self) -> f64 {
        self.0.CD0
    }
    #[getter]
    fn alpha_stall_rad(&self) -> f64 {
        self.0.alpha_stall_rad
    }

    fn cl_cd(&self, alpha_rad: f64) -> (f64, f64) {
        self.0.cl_cd(alpha_rad)
    }

    fn cl_cd_arr<'py>(
        &self,
        py: Python<'py>,
        alpha_rad: PyReadonlyArray1<'py, f64>,
    ) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>)> {
        let a = alpha_rad.as_slice()?;
        let n = a.len();
        let mut cl = vec![0.0f64; n];
        let mut cd = vec![0.0f64; n];
        self.0.cl_cd_into(a, &mut cl, &mut cd);
        Ok((cl.into_pyarray_bound(py), cd.into_pyarray_bound(py)))
    }

    fn __repr__(&self) -> String {
        format!(
            "LinearPolar(CL0={}, CL_alpha_per_rad={}, CD0={}, alpha_stall_rad={})",
            self.0.CL0, self.0.CL_alpha_per_rad, self.0.CD0, self.0.alpha_stall_rad,
        )
    }

    /// Build a LinearPolar from an AirfoilProperties-like object (either the
    /// lean _dynbem.AirfoilProperties or the Python AirfoilProperties wrapper
    /// that holds a ._rust attribute).
    #[staticmethod]
    fn from_properties(airfoil: &Bound<'_, PyAny>) -> PyResult<Self> {
        // Try the lean Rust class first (direct extraction).
        let (cl0, cl_alpha, cd0, stall_deg) = if let Ok(a) = airfoil.extract::<PyAirfoilProperties>() {
            (a.0.CL0, a.0.CL_alpha_per_rad, a.0.CD0, a.0.alpha_stall_deg)
        } else {
            // Fall back to Python attribute access (Python AirfoilProperties wrapper).
            let cl0: f64 = airfoil.getattr("CL0")?.extract()?;
            let cl_alpha: f64 = airfoil.getattr("CL_alpha_per_rad")?.extract()?;
            let cd0: f64 = airfoil.getattr("CD0")?.extract()?;
            let stall: f64 = airfoil.getattr("alpha_stall_deg")?.extract()?;
            (cl0, cl_alpha, cd0, stall)
        };
        Ok(PyLinearPolar(core_::polar::LinearPolar::new(
            cl0, cl_alpha, cd0, stall_deg.to_radians(),
        )))
    }
}

#[pyclass(name = "TabulatedPolar", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyTabulatedPolar(pub core_::polar::TabulatedPolar);

#[pymethods]
impl PyTabulatedPolar {
    #[new]
    fn new<'py>(
        alpha_rad: PyReadonlyArray1<'py, f64>,
        cl: PyReadonlyArray1<'py, f64>,
        cd: PyReadonlyArray1<'py, f64>,
    ) -> PyResult<Self> {
        let a = alpha_rad.as_slice()?.to_vec();
        let l = cl.as_slice()?.to_vec();
        let d = cd.as_slice()?.to_vec();
        core_::polar::TabulatedPolar::new(a, l, d)
            .map(PyTabulatedPolar)
            .map_err(PyValueError::new_err)
    }

    fn cl_cd(&self, alpha_rad: f64) -> (f64, f64) {
        self.0.cl_cd(alpha_rad)
    }

    fn cl_cd_arr<'py>(
        &self,
        py: Python<'py>,
        alpha_rad: PyReadonlyArray1<'py, f64>,
    ) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>)> {
        let a = alpha_rad.as_slice()?;
        let n = a.len();
        let mut cl = vec![0.0f64; n];
        let mut cd = vec![0.0f64; n];
        self.0.cl_cd_into(a, &mut cl, &mut cd);
        Ok((cl.into_pyarray_bound(py), cd.into_pyarray_bound(py)))
    }

    fn __repr__(&self) -> String {
        format!("TabulatedPolar(n={})", self.0.alpha.len())
    }
}

pub fn extract_polar(obj: &Bound<'_, PyAny>) -> PyResult<core_::polar::PolarKind> {
    if let Ok(p) = obj.extract::<PyLinearPolar>() {
        return Ok(core_::polar::PolarKind::Linear(p.0));
    }
    if let Ok(p) = obj.extract::<PyTabulatedPolar>() {
        return Ok(core_::polar::PolarKind::Tabulated(p.0));
    }
    Err(PyValueError::new_err(
        "Expected LinearPolar or TabulatedPolar",
    ))
}

/// Build a default LinearPolar from a RotorDefinition's airfoil. Used when
/// the model constructor is called without an explicit polar argument
/// (mirrors the legacy Python dynbem behaviour where polar was optional).
pub fn default_polar_from_defn(
    defn: &core_::rotor_definition::RotorDefinition,
) -> core_::polar::PolarKind {
    let a = &defn.airfoil;
    core_::polar::PolarKind::Linear(core_::polar::LinearPolar::new(
        a.CL0,
        a.CL_alpha_per_rad,
        a.CD0,
        a.alpha_stall_deg.to_radians(),
    ))
}

/// Resolve the optional `polar` constructor argument: if provided, extract
/// it; otherwise build a LinearPolar from the rotor's airfoil properties.
pub fn resolve_polar(
    polar: Option<&Bound<'_, PyAny>>,
    defn: &core_::rotor_definition::RotorDefinition,
) -> PyResult<core_::polar::PolarKind> {
    match polar {
        Some(obj) => extract_polar(obj),
        None => Ok(default_polar_from_defn(defn)),
    }
}

// ===========================================================================
// Rotor definition pieces
// ===========================================================================

#[pyclass(name = "BladeGeometry", module = "dynbem._dynbem", subclass)]
#[derive(Clone, Debug)]
pub struct PyBladeGeometry(pub core_::rotor_definition::BladeGeometry);

#[pymethods]
impl PyBladeGeometry {
    #[new]
    #[pyo3(signature = (
        n_blades, radius_m, root_cutout_m, chord_m,
        twist_deg, n_elements,
        r_stations_m, chord_stations_m, twist_stations_deg,
    ))]
    fn new(
        n_blades: usize,
        radius_m: f64,
        root_cutout_m: f64,
        chord_m: f64,
        twist_deg: f64,
        n_elements: usize,
        r_stations_m: Vec<f64>,
        chord_stations_m: Vec<f64>,
        twist_stations_deg: Vec<f64>,
    ) -> Self {
        PyBladeGeometry(core_::rotor_definition::BladeGeometry {
            n_blades,
            radius_m,
            root_cutout_m,
            chord_m,
            twist_deg,
            n_elements,
            r_stations_m,
            chord_stations_m,
            twist_stations_deg,
        })
    }

    #[getter]
    fn n_blades(&self) -> usize {
        self.0.n_blades
    }
    #[getter]
    fn radius_m(&self) -> f64 {
        self.0.radius_m
    }
    #[getter]
    fn root_cutout_m(&self) -> f64 {
        self.0.root_cutout_m
    }
    #[getter]
    fn chord_m(&self) -> f64 {
        self.0.chord_m
    }
    #[getter]
    fn twist_deg(&self) -> f64 {
        self.0.twist_deg
    }
    #[getter]
    fn n_elements(&self) -> usize {
        self.0.n_elements
    }

    #[getter]
    fn span_m(&self) -> f64 {
        self.0.span_m()
    }
    #[getter]
    fn r_cp_m(&self) -> f64 {
        self.0.r_cp_m()
    }
    #[getter]
    fn disk_area_m2(&self) -> f64 {
        self.0.disk_area_m2()
    }
    #[getter]
    fn solidity(&self) -> f64 {
        self.0.solidity()
    }
    #[getter]
    fn has_radial_stations(&self) -> bool {
        self.0.has_radial_stations()
    }

    fn chord_at(&self, r: f64) -> f64 {
        self.0.chord_at(r)
    }
    fn twist_at(&self, r: f64) -> f64 {
        self.0.twist_at(r)
    }

    fn __reduce__(&self, py: Python<'_>) -> PyResult<(PyObject, PyObject)> {
        let cls: PyObject = Self::type_object_bound(py).into_any().unbind();
        let args = (
            self.0.n_blades,
            self.0.radius_m,
            self.0.root_cutout_m,
            self.0.chord_m,
            self.0.twist_deg,
            self.0.n_elements,
            self.0.r_stations_m.clone(),
            self.0.chord_stations_m.clone(),
            self.0.twist_stations_deg.clone(),
        )
            .into_py(py);
        Ok((cls, args))
    }
}

#[pyclass(name = "AirfoilProperties", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyAirfoilProperties(pub core_::rotor_definition::AirfoilProperties);

#[pymethods]
impl PyAirfoilProperties {
    #[new]
    #[pyo3(signature = (CL0, CL_alpha_per_rad, CD0, alpha_stall_deg, tip_loss))]
    #[allow(non_snake_case)]
    fn new(
        CL0: f64,
        CL_alpha_per_rad: f64,
        CD0: f64,
        alpha_stall_deg: f64,
        tip_loss: bool,
    ) -> Self {
        PyAirfoilProperties(core_::rotor_definition::AirfoilProperties {
            CL0,
            CL_alpha_per_rad,
            CD0,
            alpha_stall_deg,
            tip_loss,
        })
    }

}

#[pyclass(name = "ControlProperties", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyControlProperties(pub core_::rotor_definition::ControlProperties);

#[pymethods]
impl PyControlProperties {
    #[new]
    #[pyo3(signature = (swashplate_pitch_gain_rad, swashplate_phase_deg))]
    fn new(swashplate_pitch_gain_rad: f64, swashplate_phase_deg: Option<f64>) -> Self {
        PyControlProperties(core_::rotor_definition::ControlProperties {
            swashplate_pitch_gain_rad,
            swashplate_phase_deg,
        })
    }

}

#[pyclass(name = "RotorDefinition", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyRotorDefinition(pub core_::rotor_definition::RotorDefinition);

#[pymethods]
impl PyRotorDefinition {
    #[new]
    #[pyo3(signature = (blade, airfoil, control, name, description))]
    fn new(
        blade: PyBladeGeometry,
        airfoil: PyAirfoilProperties,
        control: Option<PyControlProperties>,
        name: String,
        description: String,
    ) -> Self {
        PyRotorDefinition(core_::rotor_definition::RotorDefinition {
            blade: blade.0,
            airfoil: airfoil.0,
            control: control.map(|c| c.0),
            name,
            description,
        })
    }

}

// ===========================================================================
// Rotor states
// ===========================================================================

#[pyclass(name = "QuasiStaticRotorState", module = "dynbem._dynbem")]
#[derive(Clone, Debug, Default)]
pub struct PyQuasiStaticRotorState(pub core_::rotor_state::QuasiStaticRotorState);

#[pymethods]
impl PyQuasiStaticRotorState {
    #[new]
    fn new() -> Self {
        PyQuasiStaticRotorState(core_::rotor_state::QuasiStaticRotorState)
    }

    fn to_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        Vec::<f64>::new().into_pyarray_bound(py)
    }

    fn from_array(&self, arr: PyReadonlyArray1<'_, f64>) -> PyResult<Self> {
        let a = arr.as_slice()?;
        if !a.is_empty() {
            return Err(PyValueError::new_err(format!(
                "QuasiStaticRotorState expects 0 states, got {}",
                a.len(),
            )));
        }
        Ok(PyQuasiStaticRotorState(core_::rotor_state::QuasiStaticRotorState))
    }
}

#[pyclass(name = "PittPetersRotorState", module = "dynbem._dynbem")]
#[derive(Clone, Debug, Default)]
pub struct PyPittPetersRotorState(pub core_::rotor_state::PittPetersRotorState);

#[pymethods]
impl PyPittPetersRotorState {
    #[new]
    #[pyo3(signature = (lambda_0, lambda_c, lambda_s))]
    fn new(lambda_0: f64, lambda_c: f64, lambda_s: f64) -> Self {
        PyPittPetersRotorState(core_::rotor_state::PittPetersRotorState {
            lambda_0,
            lambda_c,
            lambda_s,
        })
    }

    #[getter]
    fn lambda_0(&self) -> f64 {
        self.0.lambda_0
    }
    #[setter]
    fn set_lambda_0(&mut self, v: f64) {
        self.0.lambda_0 = v;
    }
    #[getter]
    fn lambda_c(&self) -> f64 {
        self.0.lambda_c
    }
    #[setter]
    fn set_lambda_c(&mut self, v: f64) {
        self.0.lambda_c = v;
    }
    #[getter]
    fn lambda_s(&self) -> f64 {
        self.0.lambda_s
    }
    #[setter]
    fn set_lambda_s(&mut self, v: f64) {
        self.0.lambda_s = v;
    }

    fn to_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vec![self.0.lambda_0, self.0.lambda_c, self.0.lambda_s].into_pyarray_bound(py)
    }

    fn from_array(&self, arr: PyReadonlyArray1<'_, f64>) -> PyResult<Self> {
        let a = arr.as_slice()?;
        if a.len() != 3 {
            return Err(PyValueError::new_err(format!(
                "PittPetersRotorState expects 3 states, got {}",
                a.len(),
            )));
        }
        Ok(PyPittPetersRotorState(
            core_::rotor_state::PittPetersRotorState {
                lambda_0: a[0],
                lambda_c: a[1],
                lambda_s: a[2],
            },
        ))
    }
}

#[pyclass(name = "OyeRotorState", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyOyeRotorState(pub core_::rotor_state::OyeRotorState);

#[pymethods]
impl PyOyeRotorState {
    #[new]
    #[pyo3(signature = (W_int, W))]
    #[allow(non_snake_case)]
    fn new<'py>(
        W_int: PyReadonlyArray1<'py, f64>,
        W: PyReadonlyArray1<'py, f64>,
    ) -> PyResult<Self> {
        let wi = W_int.as_slice()?.to_vec();
        let w = W.as_slice()?.to_vec();
        if wi.len() != w.len() {
            return Err(PyValueError::new_err("W_int and W must have equal length"));
        }
        let n = wi.len();
        let mut s = core_::rotor_state::OyeRotorState::zeros(n);
        s.W_int[..n].copy_from_slice(&wi);
        s.W[..n].copy_from_slice(&w);
        Ok(PyOyeRotorState(s))
    }

    #[staticmethod]
    fn zeros(n_elements: usize) -> Self {
        PyOyeRotorState(core_::rotor_state::OyeRotorState::zeros(n_elements))
    }

    #[getter]
    #[allow(non_snake_case)]
    fn W_int<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.0.w_int_slice().to_vec().into_pyarray_bound(py)
    }
    #[getter]
    #[allow(non_snake_case)]
    fn W<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.0.w_slice().to_vec().into_pyarray_bound(py)
    }

    fn to_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let n = self.0.n_elements;
        let mut v = Vec::with_capacity(2 * n);
        v.extend_from_slice(self.0.w_int_slice());
        v.extend_from_slice(self.0.w_slice());
        v.into_pyarray_bound(py)
    }

    fn from_array(&self, arr: PyReadonlyArray1<'_, f64>) -> PyResult<Self> {
        let a = arr.as_slice()?;
        let n_total = a.len();
        if n_total % 2 != 0 {
            return Err(PyValueError::new_err(format!(
                "OyeRotorState array length {} invalid; expected 2*n_elements",
                n_total,
            )));
        }
        let n = n_total / 2;
        let mut s = core_::rotor_state::OyeRotorState::zeros(n);
        s.W_int[..n].copy_from_slice(&a[..n]);
        s.W[..n].copy_from_slice(&a[n..2 * n]);
        Ok(PyOyeRotorState(s))
    }
}

// ===========================================================================
// RotorInputs / AeroResult
// ===========================================================================

#[pyclass(name = "RotorInputs", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyRotorInputs(pub core_::aero_io::RotorInputs);

#[pymethods]
impl PyRotorInputs {
    #[new]
    #[pyo3(signature = (
        collective_rad, tilt_lon, tilt_lat,
        R_hub, v_hub_world, wind_world,
        omega_rad_s, t, rho_kg_m3,
    ))]
    #[allow(non_snake_case)]
    fn new<'py>(
        collective_rad: f64,
        tilt_lon: f64,
        tilt_lat: f64,
        R_hub: PyReadonlyArray2<'py, f64>,
        v_hub_world: PyReadonlyArray1<'py, f64>,
        wind_world: PyReadonlyArray1<'py, f64>,
        omega_rad_s: f64,
        t: f64,
        rho_kg_m3: f64,
    ) -> PyResult<Self> {
        Ok(PyRotorInputs(core_::aero_io::RotorInputs {
            collective_rad,
            tilt_lon,
            tilt_lat,
            R_hub: read_mat3(R_hub, "R_hub")?,
            v_hub_world: read_vec3(v_hub_world, "v_hub_world")?,
            wind_world: read_vec3(wind_world, "wind_world")?,
            t,
            rho_kg_m3,
            omega_rad_s,
        }))
    }

    #[getter]
    fn collective_rad(&self) -> f64 {
        self.0.collective_rad
    }
    #[setter]
    fn set_collective_rad(&mut self, v: f64) {
        self.0.collective_rad = v;
    }
    #[getter]
    fn tilt_lon(&self) -> f64 {
        self.0.tilt_lon
    }
    #[setter]
    fn set_tilt_lon(&mut self, v: f64) {
        self.0.tilt_lon = v;
    }
    #[getter]
    fn tilt_lat(&self) -> f64 {
        self.0.tilt_lat
    }
    #[setter]
    fn set_tilt_lat(&mut self, v: f64) {
        self.0.tilt_lat = v;
    }
    #[getter]
    fn t(&self) -> f64 {
        self.0.t
    }
    #[setter]
    fn set_t(&mut self, v: f64) {
        self.0.t = v;
    }
    #[getter]
    fn rho_kg_m3(&self) -> f64 {
        self.0.rho_kg_m3
    }
    #[setter]
    fn set_rho_kg_m3(&mut self, v: f64) {
        self.0.rho_kg_m3 = v;
    }
    #[getter]
    fn omega_rad_s(&self) -> f64 {
        self.0.omega_rad_s
    }
    #[setter]
    fn set_omega_rad_s(&mut self, v: f64) {
        self.0.omega_rad_s = v;
    }

    #[getter]
    #[allow(non_snake_case)]
    fn R_hub<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f64>> {
        mat3_to_py(py, &self.0.R_hub)
    }
    #[getter]
    fn v_hub_world<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vec3_to_py(py, &self.0.v_hub_world)
    }
    #[getter]
    fn wind_world<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vec3_to_py(py, &self.0.wind_world)
    }
}

#[pyclass(name = "AeroResult", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyAeroResult(pub core_::aero_io::AeroResult);

#[pymethods]
impl PyAeroResult {
    #[getter]
    #[allow(non_snake_case)]
    fn F_world<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vec3_to_py(py, &self.0.F_world)
    }
    #[getter]
    #[allow(non_snake_case)]
    fn M_orbital<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vec3_to_py(py, &self.0.M_orbital)
    }
    #[getter]
    #[allow(non_snake_case)]
    fn Q_spin(&self) -> f64 {
        self.0.Q_spin
    }
    #[getter]
    #[allow(non_snake_case)]
    fn M_spin<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vec3_to_py(py, &self.0.M_spin)
    }

    fn __repr__(&self) -> String {
        format!(
            "AeroResult(F_world={:?}, Q_spin={})",
            self.0.F_world, self.0.Q_spin,
        )
    }
}

// ===========================================================================
// Aero models (BEM, Pitt-Peters, Oye)
// ===========================================================================

use dynbem_rs::aero_model::AeroModel as _;

#[pyclass(name = "QuasiStaticBEM", module = "dynbem._dynbem", subclass)]
#[derive(Clone)]
pub struct PyQuasiStaticBEM(pub Box<core_::quasi_static_bem::QuasiStaticBEM>);

#[pymethods]
impl PyQuasiStaticBEM {
    #[new]
    #[pyo3(signature = (defn, polar, n_psi_elements))]
    fn new(
        defn: PyRotorDefinition,
        polar: Option<&Bound<'_, PyAny>>,
        n_psi_elements: usize,
    ) -> PyResult<Self> {
        let polar = resolve_polar(polar, &defn.0)?;
        Ok(PyQuasiStaticBEM(Box::new(
            core_::quasi_static_bem::QuasiStaticBEM::build(defn.0, n_psi_elements, polar),
        )))
    }

    fn initial_rotor_state(&self) -> PyQuasiStaticRotorState {
        PyQuasiStaticRotorState(self.0.initial_state())
    }

    fn compute_forces(
        &self,
        inputs: &PyRotorInputs,
        state: &PyQuasiStaticRotorState,
    ) -> (PyAeroResult, PyQuasiStaticRotorState) {
        let (r, s) = self.0.compute_forces(&inputs.0, &state.0);
        (PyAeroResult(r), PyQuasiStaticRotorState(s))
    }

    fn inflow_taus<'py>(
        &self,
        py: Python<'py>,
        inputs: &PyRotorInputs,
        state: &PyQuasiStaticRotorState,
    ) -> Bound<'py, PyArray1<f64>> {
        self.0
            .inflow_taus(&inputs.0, &state.0)
            .into_pyarray_bound(py)
    }

    #[getter]
    fn defn(&self) -> PyRotorDefinition {
        PyRotorDefinition(self.0.defn.clone())
    }
    #[getter]
    fn n_psi_elements(&self) -> usize {
        self.0.n_psi_elements
    }
}

#[pyclass(name = "PittPetersModel", module = "dynbem._dynbem", subclass)]
#[derive(Clone)]
pub struct PyPittPetersModel(pub Box<core_::pitt_peters::PittPetersModel>);

#[pymethods]
impl PyPittPetersModel {
    #[new]
    #[pyo3(signature = (defn, polar, n_psi_elements))]
    fn new(
        defn: PyRotorDefinition,
        polar: Option<&Bound<'_, PyAny>>,
        n_psi_elements: usize,
    ) -> PyResult<Self> {
        let polar = resolve_polar(polar, &defn.0)?;
        Ok(PyPittPetersModel(Box::new(
            core_::pitt_peters::PittPetersModel::build(defn.0, n_psi_elements, polar),
        )))
    }

    fn initial_rotor_state(&self) -> PyPittPetersRotorState {
        PyPittPetersRotorState(self.0.initial_state())
    }

    fn compute_forces(
        &self,
        inputs: &PyRotorInputs,
        state: &PyPittPetersRotorState,
    ) -> (PyAeroResult, PyPittPetersRotorState) {
        let (r, s) = self.0.compute_forces(&inputs.0, &state.0);
        (PyAeroResult(r), PyPittPetersRotorState(s))
    }

    fn inflow_taus<'py>(
        &self,
        py: Python<'py>,
        inputs: &PyRotorInputs,
        state: &PyPittPetersRotorState,
    ) -> Bound<'py, PyArray1<f64>> {
        self.0
            .inflow_taus(&inputs.0, &state.0)
            .into_pyarray_bound(py)
    }

    #[getter]
    fn defn(&self) -> PyRotorDefinition {
        PyRotorDefinition(self.0.defn.clone())
    }
    #[getter]
    fn n_psi_elements(&self) -> usize {
        self.0.n_psi_elements
    }
}

#[pyclass(name = "OyeBEMModel", module = "dynbem._dynbem", subclass)]
#[derive(Clone)]
pub struct PyOyeBEMModel(pub Box<core_::oye::OyeBEMModel>);

#[pymethods]
impl PyOyeBEMModel {
    #[new]
    #[pyo3(signature = (defn, polar, n_psi_elements, coupling_k))]
    fn new(
        defn: PyRotorDefinition,
        polar: Option<&Bound<'_, PyAny>>,
        n_psi_elements: usize,
        coupling_k: f64,
    ) -> PyResult<Self> {
        let polar = resolve_polar(polar, &defn.0)?;
        let grid = core_::bem_common::RadialGrid::from_blade(&defn.0.blade);
        Ok(PyOyeBEMModel(Box::new(core_::oye::OyeBEMModel {
            defn: defn.0,
            n_psi_elements,
            coupling_k,
            polar,
            grid,
        })))
    }

    fn initial_rotor_state(&self) -> PyOyeRotorState {
        PyOyeRotorState(self.0.initial_state())
    }

    fn compute_forces(
        &self,
        inputs: &PyRotorInputs,
        state: &PyOyeRotorState,
    ) -> (PyAeroResult, PyOyeRotorState) {
        let (r, s) = self.0.compute_forces(&inputs.0, &state.0);
        (PyAeroResult(r), PyOyeRotorState(s))
    }

    fn inflow_taus<'py>(
        &self,
        py: Python<'py>,
        inputs: &PyRotorInputs,
        state: &PyOyeRotorState,
    ) -> Bound<'py, PyArray1<f64>> {
        self.0
            .inflow_taus(&inputs.0, &state.0)
            .into_pyarray_bound(py)
    }

    #[getter]
    fn defn(&self) -> PyRotorDefinition {
        PyRotorDefinition(self.0.defn.clone())
    }
    #[getter]
    fn n_psi_elements(&self) -> usize {
        self.0.n_psi_elements
    }
    #[getter]
    fn coupling_k(&self) -> f64 {
        self.0.coupling_k
    }
}
