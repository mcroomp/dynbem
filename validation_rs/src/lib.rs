//! Empirical validation suite for rotor aerodynamics models.
//!
//! Runs theory validation checks against published datasets and closed-form solutions.

pub mod report;
pub mod helpers;
pub mod checks;

pub use report::{Report, Status, Row};

pub use checks::run_filtered_checks;

/// Run all 17 theory validation checks and return the report.
pub fn run_theory_validation() -> Report {
    let mut report = Report::new();
    checks::run_all_checks(&mut report);
    report
}

/// Run only the checks whose module name contains `filter`.
pub fn run_theory_validation_filtered(filter: &str) -> Report {
    let mut report = Report::new();
    checks::run_filtered_checks(&mut report, Some(filter));
    report
}

