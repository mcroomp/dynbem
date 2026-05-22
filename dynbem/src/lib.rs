// dynbem-py: PyO3 glue crate. Wraps the pure-Rust dynbem-core API in
// #[pyclass] newtypes and registers the _dynbem_rs Python module.
// See crates/dynbem-py/CLAUDE.md.

#![allow(clippy::too_many_arguments)]

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

mod conv;
mod trim_py;
mod wrappers;

use trim_py::{relax_inflow_py, solve_trim_cyclic_py, PyTrimResult};
use wrappers::*;

#[pyfunction]
fn vrs_lambda1(lambda2: f64) -> f64 {
    dynbem_rs::common::vrs_lambda1(lambda2)
}

#[pyfunction]
#[pyo3(signature = (tilt_lon, tilt_lat, control = None))]
fn cyclic_coeffs(tilt_lon: f64, tilt_lat: f64, control: Option<PyControlProperties>) -> (f64, f64) {
    let gains = match control {
        None => dynbem_rs::cyclic::ControlGains::default(),
        Some(c) => {
            let phase = c.0.swashplate_phase_deg.unwrap_or(0.0).to_radians();
            dynbem_rs::cyclic::ControlGains {
                gain: c.0.swashplate_pitch_gain_rad,
                phase_rad: phase,
            }
        }
    };
    dynbem_rs::cyclic::cyclic_coeffs(tilt_lon, tilt_lat, gains)
}

#[pyfunction]
fn prandtl_tip_loss(n_blades: usize, x: f64, phi_rad: f64) -> f64 {
    dynbem_rs::quasi_static_bem::prandtl_tip_loss(n_blades, x, phi_rad)
}

#[pyfunction]
fn prandtl_hub_loss(n_blades: usize, x: f64, x_hub: f64, phi_rad: f64) -> f64 {
    dynbem_rs::quasi_static_bem::prandtl_hub_loss(n_blades, x, x_hub, phi_rad)
}

// ---------------------------------------------------------------------------
// solve_bem_element: per-annulus BEM solver, exposed for diagnostics +
// the spanwise-CL verification scripts. Mirrors the legacy Python
// dynbem.bem.solve_bem_element NamedTuple API.
// ---------------------------------------------------------------------------

#[pyclass(name = "BEMElementResult", module = "dynbem._dynbem")]
#[derive(Clone, Debug)]
#[allow(non_snake_case)]
pub struct PyBEMElementResult {
    #[pyo3(get)]
    pub lambda_r: f64,
    #[pyo3(get)]
    pub a_prime: f64,
    #[pyo3(get)]
    pub dT: f64,
    #[pyo3(get)]
    pub dQ: f64,
    #[pyo3(get)]
    pub momentum_residual: f64,
}

#[pyfunction]
#[pyo3(signature = (
    r, dr, chord, twist_rad, collective_rad, omega, v_climb, rho,
    n_blades, radius_m, polar, use_tip_loss,
    v_t_extra = 0.0, root_cutout_m = 0.0,
))]
fn solve_bem_element(
    r: f64,
    dr: f64,
    chord: f64,
    twist_rad: f64,
    collective_rad: f64,
    omega: f64,
    v_climb: f64,
    rho: f64,
    n_blades: usize,
    radius_m: f64,
    polar: &Bound<'_, PyAny>,
    use_tip_loss: bool,
    v_t_extra: f64,
    root_cutout_m: f64,
) -> PyResult<PyBEMElementResult> {
    let polar = wrappers::extract_polar(polar)?;
    let res = dynbem_rs::quasi_static_bem::solve_bem_element(
        r,
        dr,
        chord,
        twist_rad,
        collective_rad,
        omega,
        v_climb,
        rho,
        n_blades,
        radius_m,
        &polar,
        use_tip_loss,
        v_t_extra,
        root_cutout_m,
    );
    Ok(PyBEMElementResult {
        lambda_r: res.lambda_r,
        a_prime: res.a_prime,
        dT: res.d_t,
        dQ: res.d_q,
        momentum_residual: res.momentum_residual,
    })
}

#[pymodule]
fn _dynbem(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(vrs_lambda1, m)?)?;
    m.add_function(wrap_pyfunction!(cyclic_coeffs, m)?)?;
    m.add_function(wrap_pyfunction!(prandtl_tip_loss, m)?)?;
    m.add_function(wrap_pyfunction!(prandtl_hub_loss, m)?)?;
    m.add_function(wrap_pyfunction!(solve_bem_element, m)?)?;
    m.add_class::<PyBEMElementResult>()?;
    m.add_class::<PyLinearPolar>()?;
    m.add_class::<PyTabulatedPolar>()?;
    m.add_class::<PyKamanFlap>()?;
    m.add_class::<PyBladeGeometry>()?;
    m.add_class::<PyAirfoilProperties>()?;
    m.add_class::<PyInertiaProperties>()?;
    m.add_class::<PyControlProperties>()?;
    m.add_class::<PyAutorotationProperties>()?;
    m.add_class::<PyRotorDefinition>()?;
    m.add_class::<PyQuasiStaticRotorState>()?;
    m.add_class::<PyPittPetersRotorState>()?;
    m.add_class::<PyOyeRotorState>()?;
    m.add_class::<PyRotorInputs>()?;
    m.add_class::<PyAeroResult>()?;
    m.add_class::<PyQuasiStaticBEM>()?;
    m.add_class::<PyPittPetersModel>()?;
    m.add_class::<PyOyeBEMModel>()?;
    m.add_class::<PyTrimResult>()?;
    m.add_function(wrap_pyfunction!(solve_trim_cyclic_py, m)?)?;
    m.add_function(wrap_pyfunction!(relax_inflow_py, m)?)?;
    Ok(())
}
