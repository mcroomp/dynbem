//! Empirical validation suite for rotor aerodynamics models.
//!
//! Runs theory validation checks against published datasets and closed-form solutions.

pub mod report;
pub mod helpers;
pub mod checks;

pub use report::{Report, Status, Row};

/// Run all 16 theory validation checks and return the report.
pub fn run_theory_validation() -> Report {
    let mut report = Report::new();
    checks::run_all_checks(&mut report);
    report
}

