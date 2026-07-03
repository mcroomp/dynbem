// Barnes-Hut VPM profiling harness for external profilers (Intel VTune, perf).
//
// The release profile already builds with full debug symbols
// (`debug = true`, `strip = "none"` in the workspace Cargo.toml), so:
//
//   Build: cargo build --release -p dynbem_rs --features parallel
//   Run:   ./target/release/bh_profile.exe [N] [theta] [seconds] [seq|par]
//
// Then point VTune at target/release/bh_profile.exe (Hotspots analysis).
// The interesting symbols to look for in the breakdown are:
//   - dynbem_rs::vpm::build_node        (octree construction)
//   - dynbem_rs::vpm::FlatTree::flatten (pre-order flatten + leaf packing)
//   - dynbem_rs::vpm::FlatTree::evaluate(stackless tree walk)
//   - dynbem_rs::vpm::accumulate_chunks (the SIMD Biot-Savart kernel)
// Time in build_node/flatten/evaluate = tree overhead; time in
// accumulate_chunks = actual Biot-Savart arithmetic.
//
// Passing theta <= 0 runs the direct O(N^2) path (`induced_velocities`)
// instead of the tree, for an apples-to-apples baseline on the same cloud.
//
// The optional 4th argument selects the execution path when the binary was
// built with --features parallel:
//   par  (default) -- Rayon parallel outer target loop
//   seq            -- single-threaded path (for comparison)
//
// Defaults: N = 8000 particles, theta = 0.5, run for 10 seconds of wall clock.

#[cfg(feature = "parallel")]
use dynbem_rs::vpm_rotor::{induced_velocities_bh_seq, induced_velocities_seq};
use dynbem_rs::vpm_rotor::{induced_velocities, induced_velocities_bh, ParticleField};
use std::env;
use std::time::{Duration, Instant};

/// Build a deterministic rotor-wake-like cloud: `n` particles laid on a few
/// trailing helices with a little jitter, so the spatial distribution (and
/// hence the octree shape) resembles a real free wake rather than a uniform
/// cube.
fn wake_cloud(n: usize) -> ParticleField {
    let mut state = 0x1234_5678u32;
    let mut rng = || {
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        (state >> 8) as f32 / (1u32 << 24) as f32 - 0.5
    };

    let n_helix = 3usize; // trailing tip vortices
    let radius = 1.0f32;
    let pitch = 0.15f32; // axial advance per revolution (fraction of R)
    let revs = 6.0f32;

    let mut f = ParticleField::with_capacity(n);
    for i in 0..n {
        let h = i % n_helix;
        let t = (i / n_helix) as f32 / (n / n_helix).max(1) as f32; // 0..1 along wake
        let phase = 2.0 * std::f32::consts::PI * h as f32 / n_helix as f32;
        let ang = phase + 2.0 * std::f32::consts::PI * revs * t;
        // Contract the wake slightly downstream, add jitter.
        let r = radius * (1.0 - 0.15 * t) + 0.02 * rng();
        let x = r * ang.cos() + 0.01 * rng();
        let y = r * ang.sin() + 0.01 * rng();
        let z = -pitch * revs * t + 0.01 * rng();

        // Strength tangent to the helix, decaying downstream.
        let mag = 0.02 * (1.0 - 0.5 * t);
        let ax = -mag * ang.sin();
        let ay = mag * ang.cos();
        let az = 0.1 * mag;

        f.push([x, y, z], [ax, ay, az], 0.05 + 0.02 * t);
    }
    f
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(8000);
    let theta: f32 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(0.5);
    let seconds: f64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(10.0);
    let force_seq: bool = args.get(4).map(|s| s == "seq").unwrap_or(false);

    let field = wake_cloud(n);
    let direct = theta <= 0.0;

    #[cfg(feature = "parallel")]
    let par_label = if force_seq { "sequential (forced)" } else { "parallel (Rayon)" };
    #[cfg(not(feature = "parallel"))]
    let par_label = "sequential (no parallel feature)";

    eprintln!(
        "bh_profile: N = {}, {}, {}, target = {} s",
        field.len(),
        if direct {
            "direct O(N^2)".to_string()
        } else {
            format!("theta = {}", theta)
        },
        par_label,
        seconds
    );

    let eval = |f: &ParticleField| {
        #[cfg(feature = "parallel")]
        if force_seq {
            return if direct {
                induced_velocities_seq(f)
            } else {
                induced_velocities_bh_seq(f, theta)
            };
        }
        if direct {
            induced_velocities(f)
        } else {
            induced_velocities_bh(f, theta)
        }
    };

    // Warm up once so allocation / first-touch isn't in the sampled window.
    let mut checksum = 0.0f64;
    let warm = eval(&field);
    checksum += warm[0][0] as f64;

    let budget = Duration::from_secs_f64(seconds);
    let start = Instant::now();
    let mut iters = 0u64;
    while start.elapsed() < budget {
        let u = eval(&field);
        // Fold a value from the result so the call can't be optimized away.
        checksum += u[iters as usize % u.len()][0] as f64;
        iters += 1;
    }
    let elapsed = start.elapsed().as_secs_f64();

    eprintln!(
        "bh_profile: {} iterations in {:.2} s = {:.2} ms per eval (checksum {:.6})",
        iters,
        elapsed,
        1e3 * elapsed / iters as f64,
        checksum
    );
}
