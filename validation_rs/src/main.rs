// VPM theory-validation report binary.
//
// Runs all VPM-vs-standard-theory checks and writes a structured report to
// stdout and tmp/theory_report.txt.  The text format is optimised for AI
// analysis: every data point is on one key=value line; PASS/FAIL/INFO tokens
// are grep-friendly; section headers delimit the report.
//
// Build and run (RELEASE is mandatory -- VPM is ~50-100x slower in debug):
//   cargo run --release -p validation_rs

use validation_rs::{run_theory_validation, run_theory_validation_filtered};
use std::time::Instant;

fn main() {
    let t0 = Instant::now();

    // Optional filter: first positional argument is a substring of the module
    // name to run. Pass no argument to run all checks.
    //   cargo run --release -p validation_rs -- cyclic_phase_servo
    let filter: Option<String> = std::env::args().nth(1);

    println!("THEORY_REPORT  dynbem_rs  run_start=unix_seconds_{}  (RELEASE mode required)",
             std::time::SystemTime::now()
                 .duration_since(std::time::UNIX_EPOCH)
                 .unwrap()
                 .as_secs());
    if let Some(ref f) = filter {
        println!("Filter: running only checks matching '{f}'");
    }
    println!("Each CHECK line: module case qty vpm ref err tol PASS|FAIL|INFO");
    println!("Modules run sequentially; VPM is the ONLY model under test.");

    let report = match filter.as_deref() {
        Some(f) => run_theory_validation_filtered(f),
        None    => run_theory_validation(),
    };

    let elapsed = t0.elapsed();
    let (total, pass, fail) = report.summary();
    println!();
    println!("=== SUMMARY  total={}  pass={}  fail={}  elapsed={:.1}s", total, pass, fail, elapsed.as_secs_f64());

    // Write to tmp/
    let out_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent().unwrap()
        .join("tmp");
    std::fs::create_dir_all(&out_dir).ok();
    let out_path = out_dir.join("theory_report.txt");
    report.write_file(&out_path);
    println!("Report written -> {}", out_path.display());

    if fail > 0 {
        eprintln!("FAILED: {} check(s) did not meet tolerance", fail);
        std::process::exit(1);
    }
}
