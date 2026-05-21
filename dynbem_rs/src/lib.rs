// dynbem-core - Pure-Rust BEM / Pitt-Peters / Oye dynamic-inflow models.
// No pyo3 / numpy here. The Python bindings live in the dynbem-py glue
// crate. See ../../../CLAUDE.md for sign conventions, NED frame, the
// CCW-from-above rotor convention.

#![allow(clippy::too_many_arguments)]

pub mod aero_io;
pub mod aero_model;
pub mod bem;
pub mod bem_common;
pub mod common;
pub mod cyclic;
pub mod oye;
pub mod pitt_peters;
pub mod polar;
pub mod rotor_definition;
pub mod rotor_state;
pub mod trim;
