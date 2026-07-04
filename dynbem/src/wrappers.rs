// PyO3 wrapper newtypes around dynbem-core types. Each one is a tuple-struct
// holding the core value; #[pymethods] forwards to the core API and handles
// numpy marshalling at the boundary.

use crate::conv::{mat3_to_py, read_mat3, read_vec3, vec3_to_py};
use dynbem_rs as core_;
use dynbem_rs::polar::Polar as _;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::PyTypeInfo;

fn parse_integration_method(
    method: Option<&str>,
) -> PyResult<core_::aero_model::IntegrationMethod> {
    match method.unwrap_or("semi_implicit").to_ascii_lowercase().as_str() {
        "semi_implicit" | "semi-implicit" | "semiimplicit" | "implicit" => {
            Ok(core_::aero_model::IntegrationMethod::SemiImplicitEuler)
        }
        "explicit" | "explicit_euler" | "explicit-euler" => {
            Ok(core_::aero_model::IntegrationMethod::ExplicitEuler)
        }
        "exponential" | "exponential_relaxation" | "exponential-relaxation" | "exp" => {
            Ok(core_::aero_model::IntegrationMethod::ExponentialRelaxation)
        }
        other => Err(PyValueError::new_err(format!(
            "Unknown integration_method '{other}'. Use 'semi_implicit', 'explicit', or 'exponential'."
        ))),
    }
}

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

    /// Build a LinearPolar from a LinearPolarParameters-like object (either the
    /// lean _dynbem.LinearPolarParameters or the Python LinearPolarParameters wrapper
    /// that holds a ._rust attribute).
    #[staticmethod]
    fn from_properties(airfoil: &Bound<'_, PyAny>) -> PyResult<Self> {
        // Try the lean Rust class first (direct extraction).
        let (cl0, cl_alpha, cd0, stall_deg) =
            if let Ok(a) = airfoil.extract::<PyLinearPolarParameters>() {
                (a.0.CL0, a.0.CL_alpha_per_rad, a.0.CD0, a.0.alpha_stall_deg)
            } else {
                // Fall back to Python attribute access (Python LinearPolarParameters wrapper).
                let cl0: f64 = airfoil.getattr("CL0")?.extract()?;
                let cl_alpha: f64 = airfoil.getattr("CL_alpha_per_rad")?.extract()?;
                let cd0: f64 = airfoil.getattr("CD0")?.extract()?;
                let stall: f64 = airfoil.getattr("alpha_stall_deg")?.extract()?;
                (cl0, cl_alpha, cd0, stall)
            };
        Ok(PyLinearPolar(core_::polar::LinearPolar::new(
            cl0,
            cl_alpha,
            cd0,
            stall_deg.to_radians(),
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

pub enum ResolvedPolar {
    Linear(core_::polar::LinearPolar),
    Tabulated(core_::polar::TabulatedPolar),
}

pub fn extract_polar(obj: &Bound<'_, PyAny>) -> PyResult<ResolvedPolar> {
    if let Ok(p) = obj.extract::<PyLinearPolar>() {
        return Ok(ResolvedPolar::Linear(p.0));
    }
    if let Ok(p) = obj.extract::<PyTabulatedPolar>() {
        return Ok(ResolvedPolar::Tabulated(p.0));
    }
    Err(PyValueError::new_err(
        "Expected LinearPolar or TabulatedPolar",
    ))
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
        tip_loss = true,
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
        tip_loss: bool,
    ) -> Self {
        PyBladeGeometry(core_::rotor_definition::BladeGeometry {
            n_blades,
            radius_m,
            root_cutout_m,
            chord_m,
            twist_deg,
            n_elements,
            tip_loss,
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
    fn tip_loss(&self) -> bool {
        self.0.tip_loss
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
            self.0.tip_loss,
        )
            .into_py(py);
        Ok((cls, args))
    }
}

#[pyclass(name = "LinearPolarParameters", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyLinearPolarParameters(pub core_::rotor_definition::LinearPolarParameters);

#[pymethods]
impl PyLinearPolarParameters {
    #[new]
    #[pyo3(signature = (CL0, CL_alpha_per_rad, CD0, alpha_stall_deg))]
    #[allow(non_snake_case)]
    fn new(CL0: f64, CL_alpha_per_rad: f64, CD0: f64, alpha_stall_deg: f64) -> Self {
        PyLinearPolarParameters(core_::rotor_definition::LinearPolarParameters {
            CL0,
            CL_alpha_per_rad,
            CD0,
            alpha_stall_deg,
        })
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
    fn alpha_stall_deg(&self) -> f64 {
        self.0.alpha_stall_deg
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

    #[getter]
    fn swashplate_pitch_gain_rad(&self) -> f64 {
        self.0.swashplate_pitch_gain_rad
    }
    #[getter]
    fn swashplate_phase_deg(&self) -> Option<f64> {
        self.0.swashplate_phase_deg
    }
}

#[pyclass(name = "ServoFlapGeometry", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyServoFlapGeometry(pub core_::rotor_definition::ServoFlapGeometry);

#[pymethods]
impl PyServoFlapGeometry {
    #[new]
    #[pyo3(signature = (C_M_delta_per_rad, r_inner_m, r_outer_m))]
    #[allow(non_snake_case)]
    fn new(C_M_delta_per_rad: f64, r_inner_m: f64, r_outer_m: f64) -> Self {
        PyServoFlapGeometry(core_::rotor_definition::ServoFlapGeometry {
            C_M_delta_per_rad,
            r_inner_m,
            r_outer_m,
        })
    }

    #[getter]
    #[allow(non_snake_case)]
    fn C_M_delta_per_rad(&self) -> f64 {
        self.0.C_M_delta_per_rad
    }
    #[getter]
    fn r_inner_m(&self) -> f64 {
        self.0.r_inner_m
    }
    #[getter]
    fn r_outer_m(&self) -> f64 {
        self.0.r_outer_m
    }
}

#[pyclass(name = "ServoFlapActuation", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyServoFlapActuation(pub core_::rotor_definition::ServoFlapActuation);

#[pymethods]
impl PyServoFlapActuation {
    #[new]
    #[pyo3(signature = (I_theta_kgm2, damper_Nms_per_rad, flap, ac_offset_m=0.0, blade_Cm_AC=0.0))]
    #[allow(non_snake_case)]
    fn new(
        I_theta_kgm2: f64,
        damper_Nms_per_rad: f64,
        flap: PyServoFlapGeometry,
        ac_offset_m: f64,
        blade_Cm_AC: f64,
    ) -> Self {
        PyServoFlapActuation(core_::rotor_definition::ServoFlapActuation {
            I_theta_kgm2,
            damper_Nms_per_rad,
            ac_offset_m,
            blade_Cm_AC,
            flap: flap.0,
        })
    }

    #[getter]
    #[allow(non_snake_case)]
    fn I_theta_kgm2(&self) -> f64 {
        self.0.I_theta_kgm2
    }
    #[getter]
    #[allow(non_snake_case)]
    fn damper_Nms_per_rad(&self) -> f64 {
        self.0.damper_Nms_per_rad
    }
    #[getter]
    fn ac_offset_m(&self) -> f64 {
        self.0.ac_offset_m
    }
    #[getter]
    #[allow(non_snake_case)]
    fn blade_Cm_AC(&self) -> f64 {
        self.0.blade_Cm_AC
    }
    #[getter]
    fn flap(&self) -> PyServoFlapGeometry {
        PyServoFlapGeometry(self.0.flap.clone())
    }
}

#[pyclass(name = "FlapProperties", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyFlapProperties(pub core_::rotor_definition::FlapProperties);

#[pymethods]
impl PyFlapProperties {
    #[new]
    #[pyo3(signature = (I_blade_flap_kgm2, omega_nr_rad_s))]
    #[allow(non_snake_case)]
    fn new(I_blade_flap_kgm2: f64, omega_nr_rad_s: f64) -> Self {
        PyFlapProperties(core_::rotor_definition::FlapProperties {
            I_blade_flap_kgm2,
            omega_nr_rad_s,
        })
    }

    #[getter]
    #[allow(non_snake_case)]
    fn I_blade_flap_kgm2(&self) -> f64 {
        self.0.I_blade_flap_kgm2
    }
    #[getter]
    fn omega_nr_rad_s(&self) -> f64 {
        self.0.omega_nr_rad_s
    }

    /// Compute the hub moment reduction factor at a given rotor speed.
    fn hub_moment_factor(&self, omega_rad_s: f64) -> f64 {
        self.0.hub_moment_factor(omega_rad_s)
    }
}

#[pyclass(name = "RotorDefinition", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyRotorDefinition(pub core_::rotor_definition::RotorDefinition);

#[pymethods]
impl PyRotorDefinition {
    #[new]
    #[pyo3(signature = (blade, airfoil, control, name, description, servoflap=None, flap=None))]
    fn new(
        blade: PyBladeGeometry,
        airfoil: PyLinearPolarParameters,
        control: Option<PyControlProperties>,
        name: String,
        description: String,
        servoflap: Option<PyServoFlapActuation>,
        flap: Option<PyFlapProperties>,
    ) -> Self {
        let pitch_actuation = match servoflap {
            Some(s) => core_::rotor_definition::PitchActuation::ServoFlap(s.0),
            None => core_::rotor_definition::PitchActuation::DirectMechanical,
        };
        PyRotorDefinition(core_::rotor_definition::RotorDefinition {
            blade: blade.0,
            airfoil: airfoil.0,
            control: control.map(|c| c.0),
            pitch_actuation,
            flap: flap.map(|f| f.0),
            name,
            description,
        })
    }

    #[getter]
    fn blade(&self) -> PyBladeGeometry {
        PyBladeGeometry(self.0.blade.clone())
    }
    #[getter]
    fn airfoil(&self) -> PyLinearPolarParameters {
        PyLinearPolarParameters(self.0.airfoil.clone())
    }
    #[getter]
    fn name(&self) -> &str {
        &self.0.name
    }
    #[getter]
    fn description(&self) -> &str {
        &self.0.description
    }
    #[getter]
    fn control(&self) -> Option<PyControlProperties> {
        self.0.control.clone().map(PyControlProperties)
    }
    /// ServoFlapActuation when in servo-flap mode, else None (direct mechanical).
    #[getter]
    fn servoflap(&self) -> Option<PyServoFlapActuation> {
        match &self.0.pitch_actuation {
            core_::rotor_definition::PitchActuation::ServoFlap(act) => {
                Some(PyServoFlapActuation(act.clone()))
            }
            core_::rotor_definition::PitchActuation::DirectMechanical => None,
        }
    }
    /// FlapProperties for quasi-static blade flapping, or None (rigid blade).
    #[getter]
    fn flap(&self) -> Option<PyFlapProperties> {
        self.0.flap.clone().map(PyFlapProperties)
    }
}

// ===========================================================================
// Rotor states
// ===========================================================================

#[pyclass(name = "QuasiStaticRotorState", module = "dynbem._dynbem")]
#[derive(Clone, Debug, Default)]
pub struct PyQuasiStaticRotorState(pub core_::quasi_static_bem::QuasiStaticRotorState);

#[pymethods]
impl PyQuasiStaticRotorState {
    #[new]
    fn new() -> Self {
        PyQuasiStaticRotorState(core_::quasi_static_bem::QuasiStaticRotorState)
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
        Ok(PyQuasiStaticRotorState(
            core_::quasi_static_bem::QuasiStaticRotorState,
        ))
    }
}

#[pyclass(name = "PittPetersRotorState", module = "dynbem._dynbem")]
#[derive(Clone, Debug, Default)]
pub struct PyPittPetersRotorState(pub core_::pitt_peters::PittPetersRotorState);

#[pymethods]
impl PyPittPetersRotorState {
    #[new]
    #[pyo3(signature = (lambda_0, lambda_c, lambda_s))]
    fn new(lambda_0: f64, lambda_c: f64, lambda_s: f64) -> Self {
        PyPittPetersRotorState(core_::pitt_peters::PittPetersRotorState {
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
            core_::pitt_peters::PittPetersRotorState {
                lambda_0: a[0],
                lambda_c: a[1],
                lambda_s: a[2],
            },
        ))
    }
}

#[pyclass(name = "OyeRotorState", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyOyeRotorState(pub core_::oye::OyeRotorState);

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
        let mut s = core_::oye::OyeRotorState::zeros(n);
        s.W_int[..n].copy_from_slice(&wi);
        s.W[..n].copy_from_slice(&w);
        Ok(PyOyeRotorState(s))
    }

    #[staticmethod]
    fn zeros(n_elements: usize) -> Self {
        PyOyeRotorState(core_::oye::OyeRotorState::zeros(n_elements))
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
        let mut s = core_::oye::OyeRotorState::zeros(n);
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
pub struct PyRotorInputs {
    pub inner: core_::aero_io::RotorInputs,
}

#[pymethods]
impl PyRotorInputs {
    #[new]
    #[pyo3(signature = (
        collective_rad, tilt_lon, tilt_lat,
        R_hub, v_hub_world, wind_world,
        omega_rad_s, rho_kg_m3,
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
        rho_kg_m3: f64,
    ) -> PyResult<Self> {
        Ok(PyRotorInputs {
            inner: core_::aero_io::RotorInputs {
                collective_rad,
                tilt_lon,
                tilt_lat,
                R_hub: read_mat3(R_hub, "R_hub")?,
                v_hub_world: read_vec3(v_hub_world, "v_hub_world")?,
                wind_world: read_vec3(wind_world, "wind_world")?,
                rho_kg_m3,
                omega_rad_s,
            },
        })
    }

    #[getter]
    fn collective_rad(&self) -> f64 {
        self.inner.collective_rad
    }
    #[setter]
    fn set_collective_rad(&mut self, v: f64) {
        self.inner.collective_rad = v;
    }
    #[getter]
    fn tilt_lon(&self) -> f64 {
        self.inner.tilt_lon
    }
    #[setter]
    fn set_tilt_lon(&mut self, v: f64) {
        self.inner.tilt_lon = v;
    }
    #[getter]
    fn tilt_lat(&self) -> f64 {
        self.inner.tilt_lat
    }
    #[setter]
    fn set_tilt_lat(&mut self, v: f64) {
        self.inner.tilt_lat = v;
    }
    #[getter]
    fn rho_kg_m3(&self) -> f64 {
        self.inner.rho_kg_m3
    }
    #[setter]
    fn set_rho_kg_m3(&mut self, v: f64) {
        self.inner.rho_kg_m3 = v;
    }
    #[getter]
    fn omega_rad_s(&self) -> f64 {
        self.inner.omega_rad_s
    }
    #[setter]
    fn set_omega_rad_s(&mut self, v: f64) {
        self.inner.omega_rad_s = v;
    }

    #[getter]
    #[allow(non_snake_case)]
    fn R_hub<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f64>> {
        mat3_to_py(py, &self.inner.R_hub)
    }
    #[getter]
    fn v_hub_world<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vec3_to_py(py, &self.inner.v_hub_world)
    }
    #[getter]
    fn wind_world<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vec3_to_py(py, &self.inner.wind_world)
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
    fn m_hub_world<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vec3_to_py(py, &self.0.M_hub_world)
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
// Aero models (BEM, Pitt-Peters, Oye) — two variants each (Linear / Tabulated)
// ===========================================================================

use dynbem_rs::aero_model::AeroModel as _;

// ---------------------------------------------------------------------------
// QuasiStaticBEM
// ---------------------------------------------------------------------------

#[pyclass(name = "_QuasiStaticBEMLinear", module = "dynbem._dynbem", subclass)]
#[derive(Clone)]
pub struct PyQuasiStaticBEMLinear(
    pub Box<core_::quasi_static_bem::QuasiStaticBEM<core_::polar::LinearPolar>>,
);

#[pymethods]
impl PyQuasiStaticBEMLinear {
    #[new]
    #[pyo3(signature = (defn, polar, n_psi_elements))]
    fn new(defn: PyRotorDefinition, polar: PyLinearPolar, n_psi_elements: usize) -> Self {
        PyQuasiStaticBEMLinear(Box::new(core_::quasi_static_bem::QuasiStaticBEM::build(
            defn.0,
            n_psi_elements,
            polar.0,
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
        let (r, s) = self.0.compute_forces(&inputs.inner, &state.0);
        (PyAeroResult(r), PyQuasiStaticRotorState(s))
    }
    #[pyo3(signature = (inputs, state, dt, integration_method=None))]
    fn step(
        &self,
        inputs: &PyRotorInputs,
        state: &PyQuasiStaticRotorState,
        dt: f64,
        integration_method: Option<&str>,
    ) -> PyResult<(PyAeroResult, PyQuasiStaticRotorState)> {
        let method = parse_integration_method(integration_method)?;
        let (r, s) = self.0.step(&inputs.inner, &state.0, dt, method);
        Ok((PyAeroResult(r), PyQuasiStaticRotorState(s)))
    }
    fn inflow_taus<'py>(
        &self,
        py: Python<'py>,
        inputs: &PyRotorInputs,
        state: &PyQuasiStaticRotorState,
    ) -> Bound<'py, PyArray1<f64>> {
        self.0
            .inflow_taus(&inputs.inner, &state.0)
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

#[pyclass(name = "_QuasiStaticBEMTabulated", module = "dynbem._dynbem", subclass)]
#[derive(Clone)]
pub struct PyQuasiStaticBEMTabulated(
    pub Box<core_::quasi_static_bem::QuasiStaticBEM<core_::polar::TabulatedPolar>>,
);

#[pymethods]
impl PyQuasiStaticBEMTabulated {
    #[new]
    #[pyo3(signature = (defn, polar, n_psi_elements))]
    fn new(defn: PyRotorDefinition, polar: PyTabulatedPolar, n_psi_elements: usize) -> Self {
        PyQuasiStaticBEMTabulated(Box::new(core_::quasi_static_bem::QuasiStaticBEM::build(
            defn.0,
            n_psi_elements,
            polar.0,
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
        let (r, s) = self.0.compute_forces(&inputs.inner, &state.0);
        (PyAeroResult(r), PyQuasiStaticRotorState(s))
    }
    #[pyo3(signature = (inputs, state, dt, integration_method=None))]
    fn step(
        &self,
        inputs: &PyRotorInputs,
        state: &PyQuasiStaticRotorState,
        dt: f64,
        integration_method: Option<&str>,
    ) -> PyResult<(PyAeroResult, PyQuasiStaticRotorState)> {
        let method = parse_integration_method(integration_method)?;
        let (r, s) = self.0.step(&inputs.inner, &state.0, dt, method);
        Ok((PyAeroResult(r), PyQuasiStaticRotorState(s)))
    }
    fn inflow_taus<'py>(
        &self,
        py: Python<'py>,
        inputs: &PyRotorInputs,
        state: &PyQuasiStaticRotorState,
    ) -> Bound<'py, PyArray1<f64>> {
        self.0
            .inflow_taus(&inputs.inner, &state.0)
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

// -----
// PittPetersModel
// ---------------------------------------------------------------------------

#[pyclass(name = "_PittPetersModelLinear", module = "dynbem._dynbem", subclass)]
#[derive(Clone)]
pub struct PyPittPetersModelLinear(
    pub Box<core_::pitt_peters::PittPetersModel<core_::polar::LinearPolar>>,
);

#[pymethods]
impl PyPittPetersModelLinear {
    #[new]
    #[pyo3(signature = (defn, polar, n_psi_elements))]
    fn new(defn: PyRotorDefinition, polar: PyLinearPolar, n_psi_elements: usize) -> Self {
        PyPittPetersModelLinear(Box::new(core_::pitt_peters::PittPetersModel::build(
            defn.0,
            n_psi_elements,
            polar.0,
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
        let (r, s) = self.0.compute_forces(&inputs.inner, &state.0);
        (PyAeroResult(r), PyPittPetersRotorState(s))
    }
    #[pyo3(signature = (inputs, state, dt, integration_method=None))]
    fn step(
        &self,
        inputs: &PyRotorInputs,
        state: &PyPittPetersRotorState,
        dt: f64,
        integration_method: Option<&str>,
    ) -> PyResult<(PyAeroResult, PyPittPetersRotorState)> {
        let method = parse_integration_method(integration_method)?;
        let (r, s) = self.0.step(&inputs.inner, &state.0, dt, method);
        Ok((PyAeroResult(r), PyPittPetersRotorState(s)))
    }
    fn inflow_taus<'py>(
        &self,
        py: Python<'py>,
        inputs: &PyRotorInputs,
        state: &PyPittPetersRotorState,
    ) -> Bound<'py, PyArray1<f64>> {
        self.0
            .inflow_taus(&inputs.inner, &state.0)
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

#[pyclass(
    name = "_PittPetersModelTabulated",
    module = "dynbem._dynbem",
    subclass
)]
#[derive(Clone)]
pub struct PyPittPetersModelTabulated(
    pub Box<core_::pitt_peters::PittPetersModel<core_::polar::TabulatedPolar>>,
);

#[pymethods]
impl PyPittPetersModelTabulated {
    #[new]
    #[pyo3(signature = (defn, polar, n_psi_elements))]
    fn new(defn: PyRotorDefinition, polar: PyTabulatedPolar, n_psi_elements: usize) -> Self {
        PyPittPetersModelTabulated(Box::new(core_::pitt_peters::PittPetersModel::build(
            defn.0,
            n_psi_elements,
            polar.0,
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
        let (r, s) = self.0.compute_forces(&inputs.inner, &state.0);
        (PyAeroResult(r), PyPittPetersRotorState(s))
    }
    #[pyo3(signature = (inputs, state, dt, integration_method=None))]
    fn step(
        &self,
        inputs: &PyRotorInputs,
        state: &PyPittPetersRotorState,
        dt: f64,
        integration_method: Option<&str>,
    ) -> PyResult<(PyAeroResult, PyPittPetersRotorState)> {
        let method = parse_integration_method(integration_method)?;
        let (r, s) = self.0.step(&inputs.inner, &state.0, dt, method);
        Ok((PyAeroResult(r), PyPittPetersRotorState(s)))
    }
    fn inflow_taus<'py>(
        &self,
        py: Python<'py>,
        inputs: &PyRotorInputs,
        state: &PyPittPetersRotorState,
    ) -> Bound<'py, PyArray1<f64>> {
        self.0
            .inflow_taus(&inputs.inner, &state.0)
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

// -----
// OyeBEMModel
// ---------------------------------------------------------------------------

#[pyclass(name = "_OyeBEMModelLinear", module = "dynbem._dynbem", subclass)]
#[derive(Clone)]
pub struct PyOyeBEMModelLinear(pub Box<core_::oye::OyeBEMModel<core_::polar::LinearPolar>>);

#[pymethods]
impl PyOyeBEMModelLinear {
    #[new]
    #[pyo3(signature = (defn, polar, n_psi_elements, coupling_k))]
    fn new(
        defn: PyRotorDefinition,
        polar: PyLinearPolar,
        n_psi_elements: usize,
        coupling_k: f64,
    ) -> Self {
        PyOyeBEMModelLinear(Box::new(core_::oye::OyeBEMModel::build_with_k(
            defn.0,
            n_psi_elements,
            polar.0,
            coupling_k,
        )))
    }
    fn initial_rotor_state(&self) -> PyOyeRotorState {
        PyOyeRotorState(self.0.initial_state())
    }
    fn compute_forces(
        &self,
        inputs: &PyRotorInputs,
        state: &PyOyeRotorState,
    ) -> (PyAeroResult, PyOyeRotorState) {
        let (r, s) = self.0.compute_forces(&inputs.inner, &state.0);
        (PyAeroResult(r), PyOyeRotorState(s))
    }
    #[pyo3(signature = (inputs, state, dt, integration_method=None))]
    fn step(
        &self,
        inputs: &PyRotorInputs,
        state: &PyOyeRotorState,
        dt: f64,
        integration_method: Option<&str>,
    ) -> PyResult<(PyAeroResult, PyOyeRotorState)> {
        let method = parse_integration_method(integration_method)?;
        let (r, s) = self.0.step(&inputs.inner, &state.0, dt, method);
        Ok((PyAeroResult(r), PyOyeRotorState(s)))
    }
    fn inflow_taus<'py>(
        &self,
        py: Python<'py>,
        inputs: &PyRotorInputs,
        state: &PyOyeRotorState,
    ) -> Bound<'py, PyArray1<f64>> {
        self.0
            .inflow_taus(&inputs.inner, &state.0)
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

#[pyclass(name = "_OyeBEMModelTabulated", module = "dynbem._dynbem", subclass)]
#[derive(Clone)]
pub struct PyOyeBEMModelTabulated(pub Box<core_::oye::OyeBEMModel<core_::polar::TabulatedPolar>>);

#[pymethods]
impl PyOyeBEMModelTabulated {
    #[new]
    #[pyo3(signature = (defn, polar, n_psi_elements, coupling_k))]
    fn new(
        defn: PyRotorDefinition,
        polar: PyTabulatedPolar,
        n_psi_elements: usize,
        coupling_k: f64,
    ) -> Self {
        PyOyeBEMModelTabulated(Box::new(core_::oye::OyeBEMModel::build_with_k(
            defn.0,
            n_psi_elements,
            polar.0,
            coupling_k,
        )))
    }
    fn initial_rotor_state(&self) -> PyOyeRotorState {
        PyOyeRotorState(self.0.initial_state())
    }
    fn compute_forces(
        &self,
        inputs: &PyRotorInputs,
        state: &PyOyeRotorState,
    ) -> (PyAeroResult, PyOyeRotorState) {
        let (r, s) = self.0.compute_forces(&inputs.inner, &state.0);
        (PyAeroResult(r), PyOyeRotorState(s))
    }
    #[pyo3(signature = (inputs, state, dt, integration_method=None))]
    fn step(
        &self,
        inputs: &PyRotorInputs,
        state: &PyOyeRotorState,
        dt: f64,
        integration_method: Option<&str>,
    ) -> PyResult<(PyAeroResult, PyOyeRotorState)> {
        let method = parse_integration_method(integration_method)?;
        let (r, s) = self.0.step(&inputs.inner, &state.0, dt, method);
        Ok((PyAeroResult(r), PyOyeRotorState(s)))
    }
    fn inflow_taus<'py>(
        &self,
        py: Python<'py>,
        inputs: &PyRotorInputs,
        state: &PyOyeRotorState,
    ) -> Bound<'py, PyArray1<f64>> {
        self.0
            .inflow_taus(&inputs.inner, &state.0)
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

// ---------------------------------------------------------------------------
// VpmRotor (forward-flight free-wake VPM)
//
// Unlike the BEM-family models, VpmRotor is a time-marching free-wake solver:
// it has no valid single-shot instantaneous evaluation, so `compute_forces`
// raises instead of running. Advance it with `step(inputs, state, dt)` and let
// the caller own the settling loop. The whole wake is carried on the state
// object (VpmRotorState), NOT in the inflow vector -- `to_array` is empty.
//
// Like the BEM-family models, VpmRotor is generic over the polar type
// (VpmRotor<P: Polar>), so there are two Python classes -- _VpmRotorLinear and
// _VpmRotorTabulated -- selected by the polar passed in. The Python-side
// VpmRotor() factory dispatches on the polar instance.
// ---------------------------------------------------------------------------

#[pyclass(name = "VpmRotorState", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyVpmRotorState(pub core_::vpm_rotor::VpmRotorState);

#[pymethods]
impl PyVpmRotorState {
    #[new]
    fn new() -> Self {
        PyVpmRotorState(core_::vpm_rotor::VpmRotorState::default())
    }

    #[staticmethod]
    #[pyo3(signature = (n_elements=0))]
    fn zeros(n_elements: usize) -> Self {
        // VPM state has no scalar inflow DOFs sized by n; the wake starts empty.
        let _ = n_elements;
        PyVpmRotorState(core_::vpm_rotor::VpmRotorState::default())
    }

    /// Number of wake particles currently carried (diagnostic).
    #[getter]
    fn n_particles(&self) -> usize {
        self.0.wake.as_ref().map_or(0, |w| w.len())
    }

    /// Inflow serialization. VPM carries no scalar inflow DOFs (the free wake
    /// lives on this object), so the inflow vector is always empty.
    fn to_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        Vec::<f64>::new().into_pyarray_bound(py)
    }

    /// Round-trip the inflow vector (which must be empty for VPM). The wake is
    /// preserved -- pass the state object itself back into `step` to continue
    /// the wake; array serialization only ever covers scalar inflow DOFs.
    fn from_array(&self, arr: PyReadonlyArray1<'_, f64>) -> PyResult<Self> {
        let a = arr.as_slice()?;
        if !a.is_empty() {
            return Err(PyValueError::new_err(
                "VpmRotorState.from_array expects an empty array (no scalar inflow DOFs)",
            ));
        }
        Ok(self.clone())
    }
}

// Build a VpmRotorConfig from the individual pyo3 signature parameters. Shared
// by the Linear and Tabulated wrapper constructors so the two classes stay in
// lock-step.
#[allow(clippy::too_many_arguments)]
fn build_vpm_config(
    max_particles: usize,
    sigma: f32,
    relax: f64,
    nonlinear_lifting_line: bool,
    tip_clustering: bool,
    local_core: bool,
    barnes_hut: bool,
    bh_theta: f32,
    bh_min_particles: usize,
) -> core_::vpm_rotor::VpmRotorConfig {
    core_::vpm_rotor::VpmRotorConfig {
        max_particles,
        sigma,
        relax,
        nonlinear_lifting_line,
        tip_clustering,
        local_core,
        barnes_hut,
        bh_theta,
        bh_min_particles,
        flap_dynamics: true,
        use_rayon: true,
        use_scalar_nan_check: false,
    }
}

#[pyclass(name = "_VpmRotorLinear", module = "dynbem._dynbem", subclass)]
#[derive(Clone)]
pub struct PyVpmRotorLinear(pub Box<core_::vpm_rotor::VpmRotor<core_::polar::LinearPolar>>);

#[pymethods]
impl PyVpmRotorLinear {
    #[new]
    #[pyo3(signature = (
        defn, polar,
        max_particles=4800, sigma=0.18, relax=0.35,
        nonlinear_lifting_line=true, tip_clustering=true, local_core=true,
        barnes_hut=false, bh_theta=0.5, bh_min_particles=2048,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        defn: PyRotorDefinition,
        polar: PyLinearPolar,
        max_particles: usize,
        sigma: f32,
        relax: f64,
        nonlinear_lifting_line: bool,
        tip_clustering: bool,
        local_core: bool,
        barnes_hut: bool,
        bh_theta: f32,
        bh_min_particles: usize,
    ) -> Self {
        let config = build_vpm_config(
            max_particles,
            sigma,
            relax,
            nonlinear_lifting_line,
            tip_clustering,
            local_core,
            barnes_hut,
            bh_theta,
            bh_min_particles,
        );
        let ctrl = defn.0.control_gains();
        PyVpmRotorLinear(Box::new(core_::vpm_rotor::VpmRotor::new(
            &defn.0, polar.0, ctrl, config,
        )))
    }

    fn initial_rotor_state(&self) -> PyVpmRotorState {
        PyVpmRotorState(self.0.initial_state())
    }

    fn compute_forces(
        &self,
        _inputs: &PyRotorInputs,
        _state: &PyVpmRotorState,
    ) -> PyResult<(PyAeroResult, PyVpmRotorState)> {
        Err(PyRuntimeError::new_err(
            "VpmRotor has no single-shot compute_forces: a free-wake VPM rotor \
             must be advanced in time. Call step(inputs, state, dt) repeatedly \
             (the caller owns the settling loop).",
        ))
    }

    #[pyo3(signature = (inputs, state, dt, integration_method=None))]
    fn step(
        &self,
        inputs: &PyRotorInputs,
        state: &PyVpmRotorState,
        dt: f64,
        integration_method: Option<&str>,
    ) -> PyResult<(PyAeroResult, PyVpmRotorState)> {
        let method = parse_integration_method(integration_method)?;
        let (r, s) = self.0.step(&inputs.inner, &state.0, dt, method);
        Ok((PyAeroResult(r), PyVpmRotorState(s)))
    }
}

#[pyclass(name = "_VpmRotorTabulated", module = "dynbem._dynbem", subclass)]
#[derive(Clone)]
pub struct PyVpmRotorTabulated(pub Box<core_::vpm_rotor::VpmRotor<core_::polar::TabulatedPolar>>);

#[pymethods]
impl PyVpmRotorTabulated {
    #[new]
    #[pyo3(signature = (
        defn, polar,
        max_particles=4800, sigma=0.18, relax=0.35,
        nonlinear_lifting_line=true, tip_clustering=true, local_core=true,
        barnes_hut=false, bh_theta=0.5, bh_min_particles=2048,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        defn: PyRotorDefinition,
        polar: PyTabulatedPolar,
        max_particles: usize,
        sigma: f32,
        relax: f64,
        nonlinear_lifting_line: bool,
        tip_clustering: bool,
        local_core: bool,
        barnes_hut: bool,
        bh_theta: f32,
        bh_min_particles: usize,
    ) -> Self {
        let config = build_vpm_config(
            max_particles,
            sigma,
            relax,
            nonlinear_lifting_line,
            tip_clustering,
            local_core,
            barnes_hut,
            bh_theta,
            bh_min_particles,
        );
        let ctrl = defn.0.control_gains();
        PyVpmRotorTabulated(Box::new(core_::vpm_rotor::VpmRotor::new(
            &defn.0, polar.0, ctrl, config,
        )))
    }

    fn initial_rotor_state(&self) -> PyVpmRotorState {
        PyVpmRotorState(self.0.initial_state())
    }

    fn compute_forces(
        &self,
        _inputs: &PyRotorInputs,
        _state: &PyVpmRotorState,
    ) -> PyResult<(PyAeroResult, PyVpmRotorState)> {
        Err(PyRuntimeError::new_err(
            "VpmRotor has no single-shot compute_forces: a free-wake VPM rotor \
             must be advanced in time. Call step(inputs, state, dt) repeatedly \
             (the caller owns the settling loop).",
        ))
    }

    #[pyo3(signature = (inputs, state, dt, integration_method=None))]
    fn step(
        &self,
        inputs: &PyRotorInputs,
        state: &PyVpmRotorState,
        dt: f64,
        integration_method: Option<&str>,
    ) -> PyResult<(PyAeroResult, PyVpmRotorState)> {
        let method = parse_integration_method(integration_method)?;
        let (r, s) = self.0.step(&inputs.inner, &state.0, dt, method);
        Ok((PyAeroResult(r), PyVpmRotorState(s)))
    }
}
