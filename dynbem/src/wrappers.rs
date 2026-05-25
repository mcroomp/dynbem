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

    /// Build a LinearPolar from an AirfoilProperties (CL0, CL_alpha_per_rad,
    /// CD0, alpha_stall_deg). Mirrors the legacy Python dynbem helper.
    #[staticmethod]
    fn from_properties(airfoil: PyAirfoilProperties) -> Self {
        let a = &airfoil.0;
        PyLinearPolar(core_::polar::LinearPolar::new(
            a.CL0,
            a.CL_alpha_per_rad,
            a.CD0,
            a.alpha_stall_deg.to_radians(),
        ))
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

#[pyclass(name = "KamanFlap", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyKamanFlap(pub core_::rotor_definition::KamanFlap);

#[pymethods]
impl PyKamanFlap {
    #[new]
    #[pyo3(signature = (
        chord_fraction = None, span_start_m = None, span_end_m = None,
        tau = None, CM_gamma_per_rad = None, swashplate_load_fraction = None, notes = String::new()
    ))]
    #[allow(non_snake_case)]
    fn new(
        chord_fraction: Option<f64>,
        span_start_m: Option<f64>,
        span_end_m: Option<f64>,
        tau: Option<f64>,
        CM_gamma_per_rad: Option<f64>,
        swashplate_load_fraction: Option<f64>,
        notes: String,
    ) -> Self {
        PyKamanFlap(core_::rotor_definition::KamanFlap {
            chord_fraction,
            span_start_m,
            span_end_m,
            tau,
            CM_gamma_per_rad,
            swashplate_load_fraction,
            notes,
        })
    }

    #[getter]
    fn chord_fraction(&self) -> Option<f64> {
        self.0.chord_fraction
    }
    #[getter]
    fn span_start_m(&self) -> Option<f64> {
        self.0.span_start_m
    }
    #[getter]
    fn span_end_m(&self) -> Option<f64> {
        self.0.span_end_m
    }
    #[getter]
    fn tau(&self) -> Option<f64> {
        self.0.tau
    }
    #[getter]
    #[allow(non_snake_case)]
    fn CM_gamma_per_rad(&self) -> Option<f64> {
        self.0.CM_gamma_per_rad
    }
    #[getter]
    fn swashplate_load_fraction(&self) -> Option<f64> {
        self.0.swashplate_load_fraction
    }
    #[getter]
    fn notes(&self) -> String {
        self.0.notes.clone()
    }

    fn __reduce__(&self, py: Python<'_>) -> PyResult<(PyObject, PyObject)> {
        let cls: PyObject = Self::type_object_bound(py).into_any().unbind();
        let args = (
            self.0.chord_fraction,
            self.0.span_start_m,
            self.0.span_end_m,
            self.0.tau,
            self.0.CM_gamma_per_rad,
            self.0.swashplate_load_fraction,
            self.0.notes.clone(),
        )
            .into_py(py);
        Ok((cls, args))
    }
}

#[pyclass(name = "BladeGeometry", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyBladeGeometry(pub core_::rotor_definition::BladeGeometry);

#[pymethods]
impl PyBladeGeometry {
    #[new]
    #[pyo3(signature = (
        n_blades, radius_m, root_cutout_m, chord_m,
        twist_deg = 0.0, n_elements = 10,
        r_stations_m = Vec::<f64>::new(),
        chord_stations_m = Vec::<f64>::new(),
        twist_stations_deg = Vec::<f64>::new(),
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
    #[pyo3(signature = (
        Re_design, CL0, CL_alpha_per_rad, CD0, alpha_stall_deg,
        tip_loss = true, name = String::new(), source = String::new(),
        polar_csv = None, CD_structural = 0.0, Re_operating = None,
    ))]
    #[allow(non_snake_case)]
    fn new(
        Re_design: i64,
        CL0: f64,
        CL_alpha_per_rad: f64,
        CD0: f64,
        alpha_stall_deg: f64,
        tip_loss: bool,
        name: String,
        source: String,
        polar_csv: Option<String>,
        CD_structural: f64,
        Re_operating: Option<i64>,
    ) -> Self {
        PyAirfoilProperties(core_::rotor_definition::AirfoilProperties {
            Re_design,
            CL0,
            CL_alpha_per_rad,
            CD0,
            alpha_stall_deg,
            tip_loss,
            name,
            source,
            polar_csv,
            CD_structural,
            Re_operating,
        })
    }

    #[getter]
    #[allow(non_snake_case)]
    fn Re_design(&self) -> i64 {
        self.0.Re_design
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
    #[getter]
    fn tip_loss(&self) -> bool {
        self.0.tip_loss
    }
    #[getter]
    fn name(&self) -> String {
        self.0.name.clone()
    }
    #[getter]
    fn source(&self) -> String {
        self.0.source.clone()
    }
    #[getter]
    fn polar_csv(&self) -> Option<String> {
        self.0.polar_csv.clone()
    }
    #[getter]
    #[allow(non_snake_case)]
    fn CD_structural(&self) -> f64 {
        self.0.CD_structural
    }
    #[getter]
    #[allow(non_snake_case)]
    fn Re_operating(&self) -> Option<i64> {
        self.0.Re_operating
    }

    fn __reduce__(&self, py: Python<'_>) -> PyResult<(PyObject, PyObject)> {
        let cls: PyObject = Self::type_object_bound(py).into_any().unbind();
        let args = (
            self.0.Re_design,
            self.0.CL0,
            self.0.CL_alpha_per_rad,
            self.0.CD0,
            self.0.alpha_stall_deg,
            self.0.tip_loss,
            self.0.name.clone(),
            self.0.source.clone(),
            self.0.polar_csv.clone(),
            self.0.CD_structural,
            self.0.Re_operating,
        )
            .into_py(py);
        Ok((cls, args))
    }
}

#[pyclass(name = "InertiaProperties", module = "dynbem._dynbem")]
#[derive(Clone, Debug, Default)]
pub struct PyInertiaProperties(pub core_::rotor_definition::InertiaProperties);

#[pymethods]
impl PyInertiaProperties {
    #[new]
    #[pyo3(signature = (
        mass_kg = None, I_body_kgm2 = Vec::<f64>::new(),
        I_spin_kgm2 = None, blade_mass_kg = None,
        stationary_assembly_mass_kg = None,
        spinning_hub_shell_mass_kg = None,
        I_blade_flap_kgm2 = None,
    ))]
    #[allow(non_snake_case)]
    fn new(
        mass_kg: Option<f64>,
        I_body_kgm2: Vec<f64>,
        I_spin_kgm2: Option<f64>,
        blade_mass_kg: Option<f64>,
        stationary_assembly_mass_kg: Option<f64>,
        spinning_hub_shell_mass_kg: Option<f64>,
        I_blade_flap_kgm2: Option<f64>,
    ) -> Self {
        PyInertiaProperties(core_::rotor_definition::InertiaProperties {
            mass_kg,
            I_body_kgm2,
            I_spin_kgm2,
            blade_mass_kg,
            stationary_assembly_mass_kg,
            spinning_hub_shell_mass_kg,
            I_blade_flap_kgm2,
        })
    }

    #[getter]
    fn mass_kg(&self) -> Option<f64> {
        self.0.mass_kg
    }
    #[getter]
    #[allow(non_snake_case)]
    fn I_spin_kgm2(&self) -> Option<f64> {
        self.0.I_spin_kgm2
    }
    #[getter]
    fn blade_mass_kg(&self) -> Option<f64> {
        self.0.blade_mass_kg
    }
    #[getter]
    fn stationary_assembly_mass_kg(&self) -> Option<f64> {
        self.0.stationary_assembly_mass_kg
    }
    #[getter]
    fn spinning_hub_shell_mass_kg(&self) -> Option<f64> {
        self.0.spinning_hub_shell_mass_kg
    }
    #[getter]
    #[allow(non_snake_case)]
    fn I_blade_flap_kgm2(&self) -> Option<f64> {
        self.0.I_blade_flap_kgm2
    }
    #[getter]
    #[allow(non_snake_case)]
    fn I_body_kgm2(&self) -> Vec<f64> {
        self.0.I_body_kgm2.clone()
    }

    fn __reduce__(&self, py: Python<'_>) -> PyResult<(PyObject, PyObject)> {
        let cls: PyObject = Self::type_object_bound(py).into_any().unbind();
        let args = (
            self.0.mass_kg,
            self.0.I_body_kgm2.clone(),
            self.0.I_spin_kgm2,
            self.0.blade_mass_kg,
            self.0.stationary_assembly_mass_kg,
            self.0.spinning_hub_shell_mass_kg,
            self.0.I_blade_flap_kgm2,
        )
            .into_py(py);
        Ok((cls, args))
    }
}

#[pyclass(name = "ControlProperties", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyControlProperties(pub core_::rotor_definition::ControlProperties);

#[pymethods]
impl PyControlProperties {
    #[new]
    #[pyo3(signature = (
        swashplate_pitch_gain_rad,
        axle_attachment_length_m = None, K_cyc = None,
        swashplate_phase_deg = None, servo_slew_rate_deg_s = None,
        servo_travel_deg = None, kaman_flap = None,
    ))]
    #[allow(non_snake_case)]
    fn new(
        swashplate_pitch_gain_rad: f64,
        axle_attachment_length_m: Option<f64>,
        K_cyc: Option<f64>,
        swashplate_phase_deg: Option<f64>,
        servo_slew_rate_deg_s: Option<f64>,
        servo_travel_deg: Option<f64>,
        kaman_flap: Option<PyKamanFlap>,
    ) -> Self {
        PyControlProperties(core_::rotor_definition::ControlProperties {
            swashplate_pitch_gain_rad,
            axle_attachment_length_m,
            K_cyc,
            swashplate_phase_deg,
            servo_slew_rate_deg_s,
            servo_travel_deg,
            kaman_flap: kaman_flap.map(|k| k.0),
        })
    }

    #[getter]
    fn swashplate_pitch_gain_rad(&self) -> f64 {
        self.0.swashplate_pitch_gain_rad
    }
    #[getter]
    fn axle_attachment_length_m(&self) -> Option<f64> {
        self.0.axle_attachment_length_m
    }
    #[getter]
    #[allow(non_snake_case)]
    fn K_cyc(&self) -> Option<f64> {
        self.0.K_cyc
    }
    #[getter]
    fn swashplate_phase_deg(&self) -> Option<f64> {
        self.0.swashplate_phase_deg
    }
    #[getter]
    fn servo_slew_rate_deg_s(&self) -> Option<f64> {
        self.0.servo_slew_rate_deg_s
    }
    #[getter]
    fn servo_travel_deg(&self) -> Option<f64> {
        self.0.servo_travel_deg
    }
    #[getter]
    fn kaman_flap(&self) -> Option<PyKamanFlap> {
        self.0.kaman_flap.clone().map(PyKamanFlap)
    }

    fn __reduce__(&self, py: Python<'_>) -> PyResult<(PyObject, PyObject)> {
        let cls: PyObject = Self::type_object_bound(py).into_any().unbind();
        let args = (
            self.0.swashplate_pitch_gain_rad,
            self.0.axle_attachment_length_m,
            self.0.K_cyc,
            self.0.swashplate_phase_deg,
            self.0.servo_slew_rate_deg_s,
            self.0.servo_travel_deg,
            self.0.kaman_flap.clone().map(PyKamanFlap),
        )
            .into_py(py);
        Ok((cls, args))
    }
}

#[pyclass(name = "AutorotationProperties", module = "dynbem._dynbem")]
#[derive(Clone, Debug, Default)]
pub struct PyAutorotationProperties(pub core_::rotor_definition::AutorotationProperties);

#[pymethods]
impl PyAutorotationProperties {
    #[new]
    #[pyo3(signature = (I_ode_kgm2 = None, omega_min_rad_s = None, omega_eq_rad_s = None))]
    #[allow(non_snake_case)]
    fn new(
        I_ode_kgm2: Option<f64>,
        omega_min_rad_s: Option<f64>,
        omega_eq_rad_s: Option<f64>,
    ) -> Self {
        PyAutorotationProperties(core_::rotor_definition::AutorotationProperties {
            I_ode_kgm2,
            omega_min_rad_s,
            omega_eq_rad_s,
        })
    }

    #[getter]
    #[allow(non_snake_case)]
    fn I_ode_kgm2(&self) -> Option<f64> {
        self.0.I_ode_kgm2
    }
    #[getter]
    fn omega_min_rad_s(&self) -> Option<f64> {
        self.0.omega_min_rad_s
    }
    #[getter]
    fn omega_eq_rad_s(&self) -> Option<f64> {
        self.0.omega_eq_rad_s
    }

    fn __reduce__(&self, py: Python<'_>) -> PyResult<(PyObject, PyObject)> {
        let cls: PyObject = Self::type_object_bound(py).into_any().unbind();
        let args = (
            self.0.I_ode_kgm2,
            self.0.omega_min_rad_s,
            self.0.omega_eq_rad_s,
        )
            .into_py(py);
        Ok((cls, args))
    }
}

#[pyclass(name = "RotorDefinition", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
pub struct PyRotorDefinition(pub core_::rotor_definition::RotorDefinition);

#[pymethods]
impl PyRotorDefinition {
    #[new]
    #[pyo3(signature = (
        blade, airfoil,
        control = None,
        inertia = PyInertiaProperties::default(),
        autorotation = PyAutorotationProperties::default(),
        name = String::new(),
        description = String::new(),
    ))]
    fn new(
        blade: PyBladeGeometry,
        airfoil: PyAirfoilProperties,
        control: Option<PyControlProperties>,
        inertia: PyInertiaProperties,
        autorotation: PyAutorotationProperties,
        name: String,
        description: String,
    ) -> Self {
        PyRotorDefinition(core_::rotor_definition::RotorDefinition {
            blade: blade.0,
            airfoil: airfoil.0,
            control: control.map(|c| c.0),
            inertia: inertia.0,
            autorotation: autorotation.0,
            name,
            description,
        })
    }

    #[getter]
    fn blade(&self) -> PyBladeGeometry {
        PyBladeGeometry(self.0.blade.clone())
    }
    #[getter]
    fn airfoil(&self) -> PyAirfoilProperties {
        PyAirfoilProperties(self.0.airfoil.clone())
    }
    #[getter]
    fn control(&self) -> Option<PyControlProperties> {
        self.0.control.clone().map(PyControlProperties)
    }
    #[getter]
    fn inertia(&self) -> PyInertiaProperties {
        PyInertiaProperties(self.0.inertia.clone())
    }
    #[getter]
    fn autorotation(&self) -> PyAutorotationProperties {
        PyAutorotationProperties(self.0.autorotation.clone())
    }
    #[getter]
    fn name(&self) -> String {
        self.0.name.clone()
    }
    #[getter]
    fn description(&self) -> String {
        self.0.description.clone()
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

    fn __reduce__(&self, py: Python<'_>) -> PyResult<(PyObject, PyObject)> {
        let cls: PyObject = Self::type_object_bound(py).into_any().unbind();
        let args = (
            PyBladeGeometry(self.0.blade.clone()),
            PyAirfoilProperties(self.0.airfoil.clone()),
            self.0.control.clone().map(PyControlProperties),
            PyInertiaProperties(self.0.inertia.clone()),
            PyAutorotationProperties(self.0.autorotation.clone()),
            self.0.name.clone(),
            self.0.description.clone(),
        )
            .into_py(py);
        Ok((cls, args))
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
    #[pyo3(signature = (omega_rad_s = 0.0, spin_angle_rad = 0.0))]
    fn new(omega_rad_s: f64, spin_angle_rad: f64) -> Self {
        PyQuasiStaticRotorState(core_::rotor_state::QuasiStaticRotorState {
            omega_rad_s,
            spin_angle_rad,
        })
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
    fn spin_angle_rad(&self) -> f64 {
        self.0.spin_angle_rad
    }
    #[setter]
    fn set_spin_angle_rad(&mut self, v: f64) {
        self.0.spin_angle_rad = v;
    }

    fn to_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vec![self.0.omega_rad_s, self.0.spin_angle_rad].into_pyarray_bound(py)
    }

    fn from_array(&self, arr: PyReadonlyArray1<'_, f64>) -> PyResult<Self> {
        let a = arr.as_slice()?;
        if a.len() != 2 {
            return Err(PyValueError::new_err(format!(
                "QuasiStaticRotorState expects 2 states, got {}",
                a.len(),
            )));
        }
        Ok(PyQuasiStaticRotorState(
            core_::rotor_state::QuasiStaticRotorState {
                omega_rad_s: a[0],
                spin_angle_rad: a[1],
            },
        ))
    }
}

#[pyclass(name = "PittPetersRotorState", module = "dynbem._dynbem")]
#[derive(Clone, Debug, Default)]
pub struct PyPittPetersRotorState(pub core_::rotor_state::PittPetersRotorState);

#[pymethods]
impl PyPittPetersRotorState {
    #[new]
    #[pyo3(signature = (
        lambda_0 = 0.0, lambda_c = 0.0, lambda_s = 0.0,
        omega_rad_s = 0.0, spin_angle_rad = 0.0,
    ))]
    fn new(
        lambda_0: f64,
        lambda_c: f64,
        lambda_s: f64,
        omega_rad_s: f64,
        spin_angle_rad: f64,
    ) -> Self {
        PyPittPetersRotorState(core_::rotor_state::PittPetersRotorState {
            lambda_0,
            lambda_c,
            lambda_s,
            omega_rad_s,
            spin_angle_rad,
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
    #[getter]
    fn omega_rad_s(&self) -> f64 {
        self.0.omega_rad_s
    }
    #[setter]
    fn set_omega_rad_s(&mut self, v: f64) {
        self.0.omega_rad_s = v;
    }
    #[getter]
    fn spin_angle_rad(&self) -> f64 {
        self.0.spin_angle_rad
    }
    #[setter]
    fn set_spin_angle_rad(&mut self, v: f64) {
        self.0.spin_angle_rad = v;
    }

    fn to_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vec![
            self.0.lambda_0,
            self.0.lambda_c,
            self.0.lambda_s,
            self.0.omega_rad_s,
            self.0.spin_angle_rad,
        ]
        .into_pyarray_bound(py)
    }

    fn from_array(&self, arr: PyReadonlyArray1<'_, f64>) -> PyResult<Self> {
        let a = arr.as_slice()?;
        if a.len() != 5 {
            return Err(PyValueError::new_err(format!(
                "PittPetersRotorState expects 5 states, got {}",
                a.len(),
            )));
        }
        Ok(PyPittPetersRotorState(
            core_::rotor_state::PittPetersRotorState {
                lambda_0: a[0],
                lambda_c: a[1],
                lambda_s: a[2],
                omega_rad_s: a[3],
                spin_angle_rad: a[4],
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
    #[pyo3(signature = (W_int, W, omega_rad_s = 0.0, spin_angle_rad = 0.0))]
    #[allow(non_snake_case)]
    fn new<'py>(
        W_int: PyReadonlyArray1<'py, f64>,
        W: PyReadonlyArray1<'py, f64>,
        omega_rad_s: f64,
        spin_angle_rad: f64,
    ) -> PyResult<Self> {
        let wi = W_int.as_slice()?.to_vec();
        let w = W.as_slice()?.to_vec();
        if wi.len() != w.len() {
            return Err(PyValueError::new_err("W_int and W must have equal length"));
        }
        let n = wi.len();
        let mut s = core_::rotor_state::OyeRotorState::zeros(n, omega_rad_s);
        s.spin_angle_rad = spin_angle_rad;
        s.W_int[..n].copy_from_slice(&wi);
        s.W[..n].copy_from_slice(&w);
        Ok(PyOyeRotorState(s))
    }

    #[staticmethod]
    #[pyo3(signature = (n_elements, omega_rad_s = 0.0))]
    fn zeros(n_elements: usize, omega_rad_s: f64) -> Self {
        PyOyeRotorState(core_::rotor_state::OyeRotorState::zeros(
            n_elements,
            omega_rad_s,
        ))
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
    #[getter]
    fn omega_rad_s(&self) -> f64 {
        self.0.omega_rad_s
    }
    #[setter]
    fn set_omega_rad_s(&mut self, v: f64) {
        self.0.omega_rad_s = v;
    }
    #[getter]
    fn spin_angle_rad(&self) -> f64 {
        self.0.spin_angle_rad
    }
    #[setter]
    fn set_spin_angle_rad(&mut self, v: f64) {
        self.0.spin_angle_rad = v;
    }

    fn to_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let n = self.0.n_elements;
        let mut v = Vec::with_capacity(2 * n + 2);
        v.extend_from_slice(self.0.w_int_slice());
        v.extend_from_slice(self.0.w_slice());
        v.push(self.0.omega_rad_s);
        v.push(self.0.spin_angle_rad);
        v.into_pyarray_bound(py)
    }

    fn from_array(&self, arr: PyReadonlyArray1<'_, f64>) -> PyResult<Self> {
        let a = arr.as_slice()?;
        let n_total = a.len();
        if n_total < 4 || (n_total - 2) % 2 != 0 {
            return Err(PyValueError::new_err(format!(
                "OyeRotorState array length {} invalid; expected 2*n_elements + 2",
                n_total,
            )));
        }
        let n = (n_total - 2) / 2;
        let mut s = core_::rotor_state::OyeRotorState::zeros(n, a[n_total - 2]);
        s.spin_angle_rad = a[n_total - 1];
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
        t = 0.0, rho_kg_m3 = 1.225, motor_torque_Nm = 0.0,
    ))]
    #[allow(non_snake_case)]
    fn new<'py>(
        collective_rad: f64,
        tilt_lon: f64,
        tilt_lat: f64,
        R_hub: PyReadonlyArray2<'py, f64>,
        v_hub_world: PyReadonlyArray1<'py, f64>,
        wind_world: PyReadonlyArray1<'py, f64>,
        t: f64,
        rho_kg_m3: f64,
        motor_torque_Nm: f64,
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
            motor_torque_Nm,
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
    #[allow(non_snake_case)]
    fn motor_torque_Nm(&self) -> f64 {
        self.0.motor_torque_Nm
    }
    #[setter]
    #[allow(non_snake_case)]
    fn set_motor_torque_Nm(&mut self, v: f64) {
        self.0.motor_torque_Nm = v;
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
    #[pyo3(signature = (defn, polar = None, n_psi_elements = 36))]
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
    #[pyo3(signature = (defn, polar = None, n_psi_elements = 36))]
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
    #[pyo3(signature = (defn, polar = None, n_psi_elements = 36, coupling_k = core_::oye::OYE_K))]
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
