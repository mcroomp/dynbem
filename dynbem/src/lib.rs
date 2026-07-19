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

#[pymodule]
fn _dynbem(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(vrs_lambda1, m)?)?;
    m.add_function(wrap_pyfunction!(cyclic_coeffs, m)?)?;
    m.add_function(wrap_pyfunction!(prandtl_tip_loss, m)?)?;
    m.add_function(wrap_pyfunction!(prandtl_hub_loss, m)?)?;
    m.add_class::<PyLinearPolar>()?;
    m.add_class::<PyTabulatedPolar>()?;
    m.add_class::<PyBladeGeometry>()?;
    m.add_class::<PyLinearPolarParameters>()?;
    m.add_class::<PyControlProperties>()?;
    m.add_class::<PyServoFlapGeometry>()?;
    m.add_class::<PyServoFlapActuation>()?;
    m.add_class::<PyFlapProperties>()?;
    m.add_class::<PyRotorDefinition>()?;
    m.add_class::<PyQuasiStaticRotorState>()?;
    m.add_class::<PyPittPetersRotorState>()?;
    m.add_class::<PyOyeRotorState>()?;
    m.add_class::<PyVpmRotorState>()?;
    m.add_class::<PyRotorInputs>()?;
    m.add_class::<PyAeroResult>()?;
    m.add_class::<PyQuasiStaticBEMLinear>()?;
    m.add_class::<PyQuasiStaticBEMTabulated>()?;
    m.add_class::<PyPittPetersModelLinear>()?;
    m.add_class::<PyPittPetersModelTabulated>()?;
    m.add_class::<PyOyeBEMModelLinear>()?;
    m.add_class::<PyOyeBEMModelTabulated>()?;
    m.add_class::<PyVpmRotorLinear>()?;
    m.add_class::<PyVpmRotorTabulated>()?;
    m.add_class::<PyTrimResult>()?;
    m.add_function(wrap_pyfunction!(solve_trim_cyclic_py, m)?)?;
    m.add_function(wrap_pyfunction!(relax_inflow_py, m)?)?;
    Ok(())
}
