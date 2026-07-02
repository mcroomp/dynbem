// VPM vs standard rotor theory -- closed-form agreement tests.
//
// Each module compares the free-wake VPM against a classical closed form in
// the regime where that closed form is valid, and documents which BEM element
// (MODEL.md section) the theory underpins -- so this suite is the record that
// the VPM reproduces the same standard theory the BEM models are built on.
//
// RUN IN RELEASE (the VPM is ~50-100x slower in debug):
//   cargo test --release --test theory
//
// Coverage (filled in over phases):
//   hover_castles_gray  -- hover C_T/C_Q vs measured Castles-Gray TN-2474
//   (blade_element_hover, bemt_combined, glauert_forward_inflow,
//    prandtl_tip_loss, climb_momentum, wake_skew, flapping_harmonics,
//    autorotation -- added in later phases)
#![allow(dead_code)]

mod common;

mod flapping_harmonics;
mod hover_castles_gray;
