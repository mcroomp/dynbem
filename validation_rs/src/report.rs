//! Structured reporting framework for validation results.

use std::io::Write;
use std::path::PathBuf;

#[derive(Debug)]
pub struct Row {
    pub module: &'static str,
    pub case: String,
    pub qty: &'static str,
    pub vpm: f64,
    pub reference: f64,
    pub err_pct: f64, // signed % error vs reference (NaN = no reference)
    pub tol_pct: f64, // tolerance threshold (NaN = info-only)
    pub status: Status,
    pub note: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Status {
    Pass,
    Fail,
    Info,
}

impl std::fmt::Display for Status {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Status::Pass => write!(f, "PASS"),
            Status::Fail => write!(f, "FAIL"),
            Status::Info => write!(f, "INFO"),
        }
    }
}

pub struct Report {
    pub rows: Vec<Row>,
    pub current_module: &'static str,
}

impl Report {
    pub fn new() -> Self {
        Self {
            rows: Vec::new(),
            current_module: "",
        }
    }

    pub fn begin_module(&mut self, name: &'static str, desc: &str) {
        self.current_module = name;
        println!();
        println!("=== MODULE {}  ({})", name, desc);
    }

    /// Record a checked quantity with a tolerance.
    pub fn check(
        &mut self,
        case: impl Into<String>,
        qty: &'static str,
        vpm: f64,
        reference: f64,
        tol_pct: f64,
    ) {
        let err_pct = if reference.abs() > 1e-15 {
            (vpm - reference) / reference * 100.0
        } else {
            f64::NAN
        };
        let pass = err_pct.abs() < tol_pct || err_pct.is_nan();
        let status = if pass { Status::Pass } else { Status::Fail };
        self.emit(
            case.into(),
            qty,
            vpm,
            reference,
            err_pct,
            tol_pct,
            status,
            None,
        );
    }

    /// Record an info-only quantity (no pass/fail assertion).
    pub fn info(&mut self, case: impl Into<String>, qty: &'static str, vpm: f64, reference: f64) {
        let err_pct = if reference.abs() > 1e-15 {
            (vpm - reference) / reference * 100.0
        } else {
            f64::NAN
        };
        self.emit(
            case.into(),
            qty,
            vpm,
            reference,
            err_pct,
            f64::NAN,
            Status::Info,
            None,
        );
    }

    /// Record a directional / boolean check.
    pub fn assert_bool(
        &mut self,
        case: impl Into<String>,
        qty: &'static str,
        vpm: f64,
        reference: f64,
        pass: bool,
        note: impl Into<String>,
    ) {
        let err_pct = if reference.abs() > 1e-15 {
            (vpm - reference) / reference * 100.0
        } else {
            f64::NAN
        };
        let status = if pass { Status::Pass } else { Status::Fail };
        self.emit(
            case.into(),
            qty,
            vpm,
            reference,
            err_pct,
            f64::NAN,
            status,
            Some(note.into()),
        );
    }

    fn emit(
        &mut self,
        case: String,
        qty: &'static str,
        vpm: f64,
        reference: f64,
        err_pct: f64,
        tol_pct: f64,
        status: Status,
        note: Option<String>,
    ) {
        let tol_str = if tol_pct.is_nan() {
            "NA".to_string()
        } else {
            format!("{:.0}%", tol_pct)
        };
        let err_str = if err_pct.is_nan() {
            "NA".to_string()
        } else {
            format!("{:+.1}%", err_pct)
        };
        let ref_str = if reference.is_nan() {
            "NA".to_string()
        } else {
            format!("{:.6}", reference)
        };
        let note_str = note.as_deref().unwrap_or("");

        println!(
            "  CHECK  module={}  case={:?}  qty={}  vpm={:.6}  ref={}  err={}  tol={}  {}  {}",
            self.current_module, case, qty, vpm, ref_str, err_str, tol_str, status, note_str
        );

        self.rows.push(Row {
            module: self.current_module,
            case,
            qty,
            vpm,
            reference,
            err_pct,
            tol_pct,
            status,
            note,
        });
    }

    pub fn summary(&self) -> (usize, usize, usize) {
        let total = self
            .rows
            .iter()
            .filter(|r| r.status != Status::Info)
            .count();
        let pass = self
            .rows
            .iter()
            .filter(|r| r.status == Status::Pass)
            .count();
        let fail = self
            .rows
            .iter()
            .filter(|r| r.status == Status::Fail)
            .count();
        (total, pass, fail)
    }

    pub fn write_file(&self, path: &PathBuf) {
        let mut f = std::fs::File::create(path).expect("create theory_report.txt");
        let (total, pass, fail) = self.summary();
        writeln!(f, "THEORY_REPORT  dynbem_rs  generated={}", chrono_now()).unwrap();
        writeln!(f, "SUMMARY  total={}  pass={}  fail={}", total, pass, fail).unwrap();
        writeln!(f).unwrap();
        for row in &self.rows {
            let tol = if row.tol_pct.is_nan() {
                "NA".to_string()
            } else {
                format!("{:.0}%", row.tol_pct)
            };
            let err = if row.err_pct.is_nan() {
                "NA".to_string()
            } else {
                format!("{:+.1}%", row.err_pct)
            };
            let note = row.note.as_deref().unwrap_or("");
            writeln!(
                f,
                "CHECK  module={}  case={:?}  qty={}  vpm={:.6}  ref={:.6}  err={}  tol={}  {}  {}",
                row.module, row.case, row.qty, row.vpm, row.reference, err, tol, row.status, note
            )
            .unwrap();
        }
    }
}

fn chrono_now() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let s = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs();
    format!("unix={}", s)
}
