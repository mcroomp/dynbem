// Ballpark comparison: quasi-static BEM vs a minimal axial free-wake VPM.
//
// This is NOT a validation study -- it is a sanity check that the VPM
// rotor coupling produces thrust/torque in the same ballpark as the BEM
// family across the axial regimes (hover, climb, descent), where momentum
// theory is cleanest and the two should agree to within modelling error.
//
// Scope / simplifications (all standard for a first free-wake coupling):
//   * Axial flow only (no cyclic, azimuthally symmetric loading).
//   * Trailing vorticity only. At steady state dGamma/dt -> 0, so the shed
//     (temporal) vorticity vanishes and the trailed (radial-gradient)
//     vorticity is the whole wake -- exactly the classic free-tip-vortex
//     wake (Landgrebe / Bagai-Leishman) discretized as particles.
//   * Lifting-line loads via Kutta-Joukowski: Gamma = 1/2 U c Cl.
//   * Blade self-induction from the near wake only (no bound-vortex self
//     term, standard for a straight lifting line).
//
// Build/run: cargo run --release -p dynbem_rs --example vpm_vs_bem_axial

use dynbem_rs::aero_io::{Mat3, RotorInputs, Vec3};
use dynbem_rs::aero_model::AeroModel;
use dynbem_rs::bem_common::{PolarTable, RadialGrid};
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::quasi_static_bem::{QuasiStaticBEM, QuasiStaticRotorState};
use dynbem_rs::rotor_definition::{
    BladeGeometry, LinearPolarParameters, PitchActuation, RotorDefinition,
};
use dynbem_rs::vpm::{advect_rk2, induced_at_points, ParticleField};
use std::f64::consts::PI;

// ---- shared rotor definition (same as the profiling harness rotor) --------

const N_BLADES: usize = 2;
const R_TIP: f64 = 1.0;
const R_ROOT: f64 = 0.2;
const CHORD: f64 = 0.06;
const TWIST_DEG: f64 = 2.0;
const CL_ALPHA: f64 = 5.7;
const CD0: f64 = 0.01;
const ALPHA_STALL_DEG: f64 = 15.0;
const RHO: f64 = 1.225;
const OMEGA: f64 = 120.0;
const COLLECTIVE_DEG: f64 = 8.0;
const N_STATIONS: usize = 20;

fn rotor_definition() -> RotorDefinition {
    RotorDefinition {
        blade: BladeGeometry {
            n_blades: N_BLADES,
            radius_m: R_TIP,
            root_cutout_m: R_ROOT,
            chord_m: CHORD,
            twist_deg: TWIST_DEG,
            n_elements: N_STATIONS,
            tip_loss: true,
            r_stations_m: Vec::new(),
            chord_stations_m: Vec::new(),
            twist_stations_deg: Vec::new(),
        },
        airfoil: LinearPolarParameters {
            CL0: 0.0,
            CL_alpha_per_rad: CL_ALPHA,
            CD0,
            alpha_stall_deg: ALPHA_STALL_DEG,
        },
        control: None,
        pitch_actuation: PitchActuation::DirectMechanical,
        flap: None,
        name: "cmp_rotor".to_string(),
        description: "vpm-vs-bem axial comparison".to_string(),
    }
}

fn linear_polar() -> LinearPolar {
    LinearPolar::new(0.0, CL_ALPHA, CD0, ALPHA_STALL_DEG.to_radians())
}

// ---- BEM reference --------------------------------------------------------

/// Returns (thrust_up_N, torque_Nm) from the quasi-static BEM at a given
/// climb rate (m/s, positive = ascending).
fn bem_thrust_torque(v_climb: f64) -> (f64, f64) {
    let defn = rotor_definition();
    let bem = QuasiStaticBEM::build(defn, 72, linear_polar());
    let state = QuasiStaticRotorState;
    let inputs = RotorInputs {
        collective_rad: COLLECTIVE_DEG.to_radians(),
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        R_hub: Mat3::eye(),
        // Ascending rotor moves in -Z (up in NED); v_rel = wind - v_hub gives
        // through-disk axial component +v_climb, matching the VPM sign below.
        v_hub_world: Vec3::new(0.0, 0.0, -v_climb),
        wind_world: Vec3::zero(),
        rho_kg_m3: RHO,
        omega_rad_s: OMEGA,
    };
    let (res, _) = bem.compute_forces(&inputs, &state);
    // Thrust "up" is the -Z component of the world force.
    (-res.F_world[2], res.Q_spin)
}

// ---- Minimal axial free-wake VPM -----------------------------------------

struct AxialVpm {
    grid: RadialGrid,
    polar: PolarTable,
    sigma: f32,
    n_steps_per_rev: usize,
    n_wake_rev: usize,
    relax: f64,
}

/// Hub-frame unit vectors at azimuth psi (CCW-from-above convention).
#[inline]
fn r_hat(psi: f64) -> [f64; 3] {
    [psi.cos(), -psi.sin(), 0.0]
}
#[inline]
fn t_hat(psi: f64) -> [f64; 3] {
    [-psi.sin(), -psi.cos(), 0.0]
}

impl AxialVpm {
    fn new() -> Self {
        let defn = rotor_definition();
        Self {
            grid: RadialGrid::from_blade(&defn.blade),
            polar: PolarTable::from_polar(&linear_polar()),
            sigma: 0.18,
            n_steps_per_rev: 24,
            n_wake_rev: 4,
            relax: 0.35,
        }
    }

    /// March the free wake to steady state and return (thrust_up_N,
    /// torque_Nm) averaged over the final revolution.
    fn thrust_torque(&self, v_climb: f64) -> (f64, f64) {
        let n = self.grid.n_elements;
        let dpsi = 2.0 * PI / self.n_steps_per_rev as f64;
        let dt = dpsi / OMEGA;
        let collective = COLLECTIVE_DEG.to_radians();

        let mut wake = ParticleField::new();
        let mut gamma = vec![0.0f64; n]; // relaxed bound circulation, blade 0
        let shed_per_step = N_BLADES * (n + 1);
        let max_particles = self.n_wake_rev * self.n_steps_per_rev * shed_per_step;

        let total_steps = 8 * self.n_steps_per_rev;
        let avg_start = total_steps - self.n_steps_per_rev;
        let mut t_acc = 0.0;
        let mut q_acc = 0.0;
        let mut avg_count = 0usize;

        for step in 0..total_steps {
            let psi0 = step as f64 * dpsi;

            // --- probe wake-induced velocity at blade-0 station midpoints ---
            let (tx, ty, tz): (Vec<f32>, Vec<f32>, Vec<f32>) = (0..n)
                .map(|i| {
                    let r = self.grid.r_mid[i];
                    let rh = r_hat(psi0);
                    ((r * rh[0]) as f32, (r * rh[1]) as f32, 0.0f32)
                })
                .fold((vec![], vec![], vec![]), |mut acc, (x, y, z)| {
                    acc.0.push(x);
                    acc.1.push(y);
                    acc.2.push(z);
                    acc
                });
            let ind = induced_at_points(&wake, &tx, &ty, &tz);

            // --- blade-element loads on blade 0 -----------------------------
            let th = t_hat(psi0);
            let mut new_gamma = vec![0.0f64; n];
            let mut thrust = 0.0;
            let mut torque = 0.0;
            for i in 0..n {
                let r = self.grid.r_mid[i];
                let c = self.grid.chord[i];
                let twist = self.grid.twist_rad[i];

                // Total air velocity = freestream (through disk, +Z = v_climb)
                // plus wake-induced velocity.
                let u_ind = [ind[i][0] as f64, ind[i][1] as f64, ind[i][2] as f64];
                let v_air = [u_ind[0], u_ind[1], v_climb + u_ind[2]];
                // Blade velocity = omega r along +t_hat.
                let vb = OMEGA * r;
                // Relative wind U_rel = v_air - v_blade.
                let u_rel = [
                    v_air[0] - vb * th[0],
                    v_air[1] - vb * th[1],
                    v_air[2] - vb * th[2],
                ];
                // Axial and tangential components of the relative wind.
                let u_a = u_rel[2]; // along +Z
                let u_t = -(u_rel[0] * th[0] + u_rel[1] * th[1]); // along -t_hat
                let u_mag = (u_a * u_a + u_t * u_t).sqrt().max(1e-6);
                let phi = u_a.atan2(u_t);
                let alpha = (collective + twist) - phi;
                let (cl, cd) = self.polar.interp(alpha);

                let q_dyn = 0.5 * RHO * u_mag * u_mag;
                let dl = q_dyn * c * cl * self.grid.dr; // lift per element
                let dd = q_dyn * c * cd * self.grid.dr; // drag per element
                thrust += dl * phi.cos() - dd * phi.sin();
                torque += (dl * phi.sin() + dd * phi.cos()) * r;

                new_gamma[i] = 0.5 * u_mag * c * cl; // Kutta-Joukowski
            }
            thrust *= N_BLADES as f64;
            torque *= N_BLADES as f64;

            // Under-relax the bound circulation for convergence.
            for i in 0..n {
                gamma[i] += self.relax * (new_gamma[i] - gamma[i]);
            }

            // --- shed trailing vorticity from every blade -------------------
            // Edge j trailing circulation = Gamma_{j-1} - Gamma_j, with
            // Gamma outside the blade = 0. Segment vector = U_rel * dt
            // (filament aligned with the local relative wind).
            for b in 0..N_BLADES {
                let psi_b = psi0 + b as f64 * 2.0 * PI / N_BLADES as f64;
                let rh = r_hat(psi_b);
                let thb = t_hat(psi_b);
                for j in 0..=n {
                    let g_in = if j == 0 { 0.0 } else { gamma[j - 1] };
                    let g_out = if j == n { 0.0 } else { gamma[j] };
                    let g_trail = g_in - g_out;
                    // Edge radius.
                    let r_edge = R_ROOT + j as f64 * self.grid.dr;
                    // Relative wind at the edge (axial + tangential only).
                    let vb = OMEGA * r_edge;
                    let seg = [
                        (-vb * thb[0]) * dt,
                        (-vb * thb[1]) * dt,
                        (v_climb) * dt,
                    ];
                    let alpha = [
                        (g_trail * seg[0]) as f32,
                        (g_trail * seg[1]) as f32,
                        (g_trail * seg[2]) as f32,
                    ];
                    let pos = [
                        (r_edge * rh[0]) as f32,
                        (r_edge * rh[1]) as f32,
                        0.0f32,
                    ];
                    wake.push(pos, alpha, self.sigma);
                }
            }

            // --- convect the free wake --------------------------------------
            advect_rk2(&mut wake, [0.0, 0.0, v_climb as f32], dt as f32);

            // --- truncate the oldest wake (FIFO) ----------------------------
            if wake.len() > max_particles {
                let excess = wake.len() - max_particles;
                drain_front(&mut wake, excess);
            }

            if step >= avg_start {
                t_acc += thrust;
                q_acc += torque;
                avg_count += 1;
            }
        }

        (t_acc / avg_count as f64, q_acc / avg_count as f64)
    }
}

/// Remove the oldest `k` particles (front of the SoA arrays).
fn drain_front(f: &mut ParticleField, k: usize) {
    f.drain_front(k);
}

fn main() {
    println!(
        "Axial comparison: rotor R={} m, {} blades, collective {} deg, omega {} rad/s\n",
        R_TIP, N_BLADES, COLLECTIVE_DEG, OMEGA
    );
    println!(
        "{:>10}  {:>18}  {:>18}  {:>18}  {:>18}",
        "regime", "BEM thrust [N]", "VPM thrust [N]", "BEM torque [Nm]", "VPM torque [Nm]"
    );

    let vpm = AxialVpm::new();
    let cases = [
        ("hover", 0.0f64),
        ("climb +3", 3.0),
        ("climb +6", 6.0),
        ("descent -3", -3.0),
        ("descent -6", -6.0),
    ];
    for (name, vc) in cases {
        let (tb, qb) = bem_thrust_torque(vc);
        let (tv, qv) = vpm.thrust_torque(vc);
        println!(
            "{:>10}  {:>18.2}  {:>18.2}  {:>18.3}  {:>18.3}",
            name, tb, tv, qb, qv
        );
    }
}
