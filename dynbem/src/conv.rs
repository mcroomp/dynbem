// Numpy <-> core type marshalling helpers. The only place in the glue
// crate where numpy is touched directly.

use dynbem_rs::aero_io::{Mat3, Vec3};
use numpy::{
    IntoPyArray, PyArray1, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2,
    PyUntypedArrayMethods,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

pub fn read_vec3(arr: PyReadonlyArray1<'_, f64>, name: &str) -> PyResult<Vec3> {
    let s = arr.as_slice()?;
    if s.len() != 3 {
        return Err(PyValueError::new_err(format!(
            "{} must be a length-3 array",
            name
        )));
    }
    Ok(Vec3([s[0], s[1], s[2]]))
}

pub fn read_mat3(arr: PyReadonlyArray2<'_, f64>, name: &str) -> PyResult<Mat3> {
    let dims = arr.shape();
    if dims != [3, 3] {
        return Err(PyValueError::new_err(format!(
            "{} must be a 3x3 array",
            name
        )));
    }
    let v = arr.as_array();
    let mut m = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] = v[[i, j]];
        }
    }
    Ok(Mat3(m))
}

pub fn vec3_to_py<'py>(py: Python<'py>, v: &Vec3) -> Bound<'py, PyArray1<f64>> {
    vec![v.0[0], v.0[1], v.0[2]].into_pyarray_bound(py)
}

pub fn mat3_to_py<'py>(py: Python<'py>, m: &Mat3) -> Bound<'py, numpy::PyArray2<f64>> {
    let flat: Vec<f64> = m.0.iter().flatten().copied().collect();
    numpy::PyArray1::from_vec_bound(py, flat)
        .reshape([3, 3])
        .unwrap()
}
