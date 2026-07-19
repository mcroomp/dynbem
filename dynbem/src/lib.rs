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

/// Rigid-body rotor spin ODE derivative. See dynbem_rs::mechanical for the
/// full physics (Coulomb bearing friction is always a parameter; pass
/// bearing_friction_nm=0.0 for a frictionless bearing).
#[pyfunction]
#[pyo3(signature = (omega, q_aero, motor_torque_nm, i_ode_kgm2, bearing_friction_nm = 0.0))]
fn omega_derivative(
    omega: f64,
    q_aero: f64,
    motor_torque_nm: f64,
    i_ode_kgm2: f64,
    bearing_friction_nm: f64,
) -> f64 {
    dynbem_rs::mechanical::omega_derivative(
        omega,
        q_aero,
        motor_torque_nm,
        i_ode_kgm2,
        bearing_friction_nm,
    )
}

/// Semi-implicit (locally-frozen relaxation) step for the rigid-body spin
/// ODE. This is the single canonical, recommended integrator -- see
/// dynbem_rs::mechanical for the full derivation and stability properties.
#[pyfunction]
#[pyo3(signature = (omega, q_aero, motor_torque_nm, i_ode_kgm2, dt, bearing_friction_nm = 0.0))]
fn step_omega(
    omega: f64,
    q_aero: f64,
    motor_torque_nm: f64,
    i_ode_kgm2: f64,
    dt: f64,
    bearing_friction_nm: f64,
) -> f64 {
    dynbem_rs::mechanical::step_omega(
        omega,
        q_aero,
        motor_torque_nm,
        i_ode_kgm2,
        dt,
        bearing_friction_nm,
    )
}

#[pymodule]
fn _dynbem(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(cyclic_coeffs, m)?)?;
    m.add_function(wrap_pyfunction!(omega_derivative, m)?)?;
    m.add_function(wrap_pyfunction!(step_omega, m)?)?;
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
