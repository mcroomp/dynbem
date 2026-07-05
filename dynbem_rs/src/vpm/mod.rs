// Forward-flight free-wake VPM rotor coupling.
//
// Extends the axial coupling to cyclic pitch and crosswind / edgewise flow,
// following standard unsteady lifting-line free-wake practice (Leishman,
// "Principles of Helicopter Aerodynamics"; Bagai & Leishman 1995).
// See VPM_DESIGN.md and validation_rs for comprehensive theory validation tests.
//
// What this adds over the axial coupling:
//   * Per-blade loading and induced-velocity probing (no azimuthal-symmetry
//     shortcut) -- every blade at its own azimuth.
//   * Cyclic pitch via `cyclic::cyclic_coeffs` (repo swashplate convention).
//   * Full 3D freestream (through-disk + in-plane), so the advancing/
//     retreating asymmetry and the skewed wake fall out naturally.
//   * Shed (temporal) vorticity  -dGamma/dt  in addition to trailed
//     (radial-gradient) vorticity. This is the load-bearing addition: once
//     bound circulation varies with azimuth (any cyclic or crosswind case)
//     Kelvin's theorem requires the shed term, or the wake is wrong.
//   * Hub pitch/roll moments  Mx = sum r dT sin(psi),  My = sum r dT cos(psi)
//     (AGENTS.md hub-frame convention), so cyclic/crosswind are observable.
//
// Frame: NED, CCW-from-above (AGENTS.md). Hub axis +Z, thrust is "up" (-Z).
// Everything here is evaluated in the hub frame (caller supplies the
// freestream already rotated into hub axes).

use crate::aero_io::{AeroResult, RotorInputs};
use crate::aero_model::{AeroModel, IntegrationMethod, RotorStateExt};
use crate::bem_common::{assemble_result, kinematics, Kinematics};
use crate::cyclic::ControlGains;
use crate::polar::Polar;
use crate::rotor_definition::{
    FlapProperties, PitchActuation, RotorDefinition, ServoFlapActuation,
};

// Wake-engine submodules. `common` is the shared foundation (particle field +
// Biot-Savart velocity kernel + Barnes-Hut + the classic convection step);
// `reformulated` is the rVPM engine; `merge` and `aging` are optional per-step
// wake operations; `march` holds the per-sub-step time-advance loop. This
// module -- the rotor coupling -- is the public face of `vpm` and selects
// between the engines via `WakeEngine`.
mod aging;
mod common;
mod march;
mod merge;
mod reformulated;

// `induced_at_points` backs the public `induced_velocities_at_points` helper
// below; the rest of the engine API is used from `march`.
use common::induced_at_points;
use std::f64::consts::PI;

// Engine primitives re-exported through `dynbem_rs::vpm::…` (used by the
// standalone profiler binaries and any external kernel benchmark). Bringing
// them in with `pub use` also puts them in scope for this module's own use.
pub use common::{
    advect_rk2, induced_velocities, induced_velocities_bh_seq, induced_velocities_seq,
    ParticleField,
};

/// Which free-wake evolution engine the rotor drives the wake with. Chosen by
/// the caller through [`VpmRotorConfig::wake_engine`]; `vpm` (this module) is
/// the common
/// coupling (shedding, lifting line, DOFs) and dispatches the per-step wake
/// convection to one of these.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum WakeEngine {
    /// Classic convection-only VPM ([`crate::vpm::advect_rk2`]): strengths and
    /// core sizes are frozen; supports Barnes-Hut and the merge/aging knobs.
    #[default]
    ClassicVpm,
    /// Reformulated VPM ([`crate::vpm_r::advect_rvpm`]): particle strength and
    /// core size evolve by vortex stretching (Alvarez & Ning). Direct O(N^2)
    /// evaluation; the merge/aging/Barnes-Hut knobs do not apply.
    ReformulatedVpm,
}

/// Solver resolution / wake settings.
#[derive(Clone, Copy, Debug)]
pub struct VpmRotorConfig {
    /// Maximum number of wake particles retained. Older particles are
    /// discarded (front-truncation) once this cap is reached. Size this to
    /// cover the desired wake age: a rough guide is
    /// `n_blades * (2*n_elements + 1) * steps_per_rev * wake_revs`.
    pub max_particles: usize,
    /// Particle core size (m).
    pub sigma: f32,
    /// Under-relaxation factor on bound circulation.
    pub relax: f64,
    /// Solve the bound circulation with an implicit (Prandtl / Weissinger)
    /// lifting-line, including the within-step trailing near-wake downwash.
    /// When false, falls back to the original pointwise Kutta-Joukowski pass
    /// (no within-step self-induction).
    pub nonlinear_lifting_line: bool,
    /// Cluster the spanwise stations toward the tip (and root) with cosine
    /// spacing, to resolve the bound-circulation roll-off. When false, the
    /// stations are uniformly spaced.
    pub tip_clustering: bool,
    /// Scale each feature's particle / filament core size with the local
    /// radial spacing (smaller cores where stations are clustered), instead
    /// of using the single global `sigma` everywhere.
    pub local_core: bool,
    /// Use the Barnes-Hut O(N log N) wake evaluator instead of the direct
    /// O(N^2) sum (approximate; controlled by `bh_theta`). Only engaged once
    /// the wake reaches `bh_min_particles`.
    pub barnes_hut: bool,
    /// Barnes-Hut opening angle (smaller = more accurate + slower).
    pub bh_theta: f32,
    /// Particle count below which the direct sum is used even with
    /// `barnes_hut` on (the tree overhead is not worth it for small wakes).
    pub bh_min_particles: usize,
    /// Model per-blade rigid flap dynamics (beta, beta_dot) when the rotor
    /// definition supplies `FlapProperties`. When false (or no FlapProperties)
    /// blades are rigid in flap and the wake stays in the disk plane.
    pub flap_dynamics: bool,
    /// Use Rayon parallelism for wake induction and convection. When false,
    /// uses sequential implementations for debugging or single-threaded execution.
    pub use_rayon: bool,
    /// Debug: route every wake induction call through the scalar f64 reference
    /// path that asserts each per-pair contribution is finite. Panics with
    /// source/target indices and values at the first NaN. Very slow -- O(N^2)
    /// non-vectorized. Default false.
    pub use_scalar_nan_check: bool,
    /// Population control: merge small, coherent, far-field wake cells into
    /// single equivalent particles (tree-collapse; see `vpm::merge_particles`).
    /// Off by default -- when off the wake is only FIFO-truncated at
    /// `max_particles`, preserving the original behaviour.
    pub merge_wake: bool,
    /// Run a merge pass every `merge_every` sub-steps (0 disables merging even
    /// when `merge_wake` is true).
    pub merge_every: usize,
    /// Skip merging until the wake has at least this many particles.
    pub merge_min_particles: usize,
    /// Merge threshold: collapse a cell when `2*half <= merge_kappa*sigma_rep`.
    /// Keep consistent with `bh_theta` (~0.5-1.0).
    pub merge_kappa: f32,
    /// Merge coherence gate: collapse only when `|sum_a| >= merge_chi_min*wsum`
    /// (near 1.0 protects thin tip vortices).
    pub merge_chi_min: f32,
    /// Merge region gate: only collapse cells whose centre is at least this far
    /// (m) from the hub origin. 0 merges everywhere.
    pub merge_region_dist: f32,
    /// Viscous core spreading: grow each particle core by `sigma^2 += 2*nu*dt`
    /// per sub-step (m^2/s). 0 disables. Conserves circulation.
    pub core_spread_nu: f64,
    /// Wake strength fade: decay each particle strength to 1/e over this many
    /// revolutions (non-conservative, models wake breakdown). 0 disables.
    pub strength_decay_tau_rev: f64,
    /// Which wake-evolution engine to use (classic convection-only VPM or the
    /// reformulated VPM with strength/size evolution). Default `ClassicVpm`.
    pub wake_engine: WakeEngine,
    /// Reformulated VPM only: Pedrizzetti relaxation factor in [0,1]. Each
    /// sub-step realigns particle strengths with the local regularized
    /// vorticity (magnitude-conserving), the primary stabilizer that keeps
    /// inviscid rVPM from diverging. 0 disables. Ignored by `ClassicVpm`.
    pub rvpm_relax: f64,
    /// Reformulated VPM only: viscous / subfilter-scale (SFS) eddy viscosity
    /// (m^2/s) applied as core spreading (sigma^2 += 2*nu*dt) inside the rVPM
    /// step. Low-order SFS energy drain. 0 disables. Ignored by `ClassicVpm`.
    pub rvpm_nu: f64,
}

impl Default for VpmRotorConfig {
    fn default() -> Self {
        Self {
            max_particles: 4800,
            sigma: 0.18,
            relax: 0.35,
            nonlinear_lifting_line: true,
            tip_clustering: true,
            local_core: true,
            barnes_hut: false,
            bh_theta: 0.5,
            bh_min_particles: 2048,
            flap_dynamics: true,
            use_rayon: true,
            use_scalar_nan_check: false,
            merge_wake: false,
            merge_every: 0,
            merge_min_particles: 0,
            merge_kappa: 0.7,
            merge_chi_min: 0.9,
            merge_region_dist: 0.0,
            core_spread_nu: 0.0,
            strength_decay_tau_rev: 0.0,
            wake_engine: WakeEngine::ClassicVpm,
            rvpm_relax: 0.3,
            rvpm_nu: 0.0,
        }
    }
}

impl VpmRotorConfig {
    /// Small, fast preset for unit tests (coarse but sign-correct).
    pub fn fast_test() -> Self {
        Self {
            max_particles: 800,
            sigma: 0.2,
            relax: 0.4,
            nonlinear_lifting_line: true,
            tip_clustering: true,
            local_core: true,
            barnes_hut: false,
            bh_theta: 0.5,
            bh_min_particles: 2048,
            flap_dynamics: true,
            use_rayon: true,
            use_scalar_nan_check: false,
            merge_wake: false,
            merge_every: 0,
            merge_min_particles: 0,
            merge_kappa: 0.7,
            merge_chi_min: 0.9,
            merge_region_dist: 0.0,
            core_spread_nu: 0.0,
            strength_decay_tau_rev: 0.0,
            wake_engine: WakeEngine::ClassicVpm,
            rvpm_relax: 0.3,
            rvpm_nu: 0.0,
        }
    }
}

/// A single steady flight condition, in hub axes.
#[derive(Clone, Copy, Debug)]
pub struct FlightCondition {
    pub collective_rad: f64,
    pub tilt_lon: f64,
    pub tilt_lat: f64,
    /// Air velocity relative to the hub, hub frame. +Z is through-disk
    /// (axial); X/Y are in-plane (edgewise / crosswind).
    pub v_hub: [f64; 3],
    pub omega_rad_s: f64,
    pub rho: f64,
}

/// Cycle-averaged loads (mean over the final revolution).
#[derive(Clone, Copy, Debug)]
pub struct VpmRotorResult {
    /// Thrust "up" (-Z hub), N.
    pub thrust: f64,
    /// Shaft reaction torque (spin), N*m.
    pub torque: f64,
    /// Hub rolling moment (roll-right positive), N*m.
    pub mx_hub: f64,
    /// Hub pitching moment (pitch-up positive), N*m.
    pub my_hub: f64,
    /// Wake particle count at the end of the run.
    pub n_particles: usize,
    /// Mean wake position (for diagnostics, e.g. skew).
    pub wake_centroid: [f64; 3],
}

/// Forward-flight free-wake VPM rotor.
#[derive(Clone)]
pub struct VpmRotor<P: Polar> {
    polar: P,
    n_blades: usize,
    n_elements: usize,
    /// Radial edges (n+1), possibly cosine-clustered toward the tip.
    r_edge: Vec<f64>,
    /// Station centres (n) -- midpoints of adjacent edges.
    r_mid: Vec<f64>,
    /// Per-station radial width (n) = r_edge[i+1] - r_edge[i].
    dr: Vec<f64>,
    /// Per-station chord (n), sampled at r_mid.
    chord: Vec<f64>,
    /// Per-station twist (rad, n), sampled at r_mid.
    twist: Vec<f64>,
    /// Per-station shed-particle core size (n).
    sigma_mid: Vec<f32>,
    /// Per-edge trailed-particle / near-wake filament core size (n+1).
    sigma_edge: Vec<f32>,
    /// Optional rigid-flap properties (inertia + non-rotating flap frequency).
    /// `Some` enables per-blade flap dynamics when `config.flap_dynamics`.
    flap: Option<FlapProperties>,
    /// Optional servo-flap feathering actuation (Kaman path). `Some` enables
    /// the per-blade feathering DOF (feathering + damper architecture).
    feather: Option<ServoFlapActuation>,
    /// Blade section lift-curve slope [1/rad], from the rotor's airfoil.
    /// Used to size the aerodynamic feathering spring k_aero from ac_offset.
    cl_alpha: f64,
    /// Spanwise integral Sum(chord[i] * r_mid[i]^2 * dr[i]) [m^4], precomputed
    /// for the aerodynamic feathering spring: k_aero = 0.5*rho*omega^2*
    /// cl_alpha*ac_offset * feather_span_integral.
    feather_span_integral: f64,
    ctrl: ControlGains,
    config: VpmRotorConfig,
}

impl<P: Polar> VpmRotor<P> {
    pub fn new(
        defn: &RotorDefinition,
        polar: P,
        ctrl: ControlGains,
        config: VpmRotorConfig,
    ) -> Self {
        let blade = &defn.blade;
        let n = blade.n_elements;
        let r_root = blade.root_cutout_m;
        let r_tip = blade.radius_m;
        let span = r_tip - r_root;

        // Radial edges: cosine-clustered toward both ends when requested
        // (clusters at the tip and root, resolving the circulation roll-off),
        // otherwise uniform.
        let r_edge: Vec<f64> = (0..=n)
            .map(|k| {
                let f = if config.tip_clustering {
                    0.5 * (1.0 - (PI * k as f64 / n as f64).cos())
                } else {
                    k as f64 / n as f64
                };
                r_root + span * f
            })
            .collect();

        let r_mid: Vec<f64> = (0..n).map(|i| 0.5 * (r_edge[i] + r_edge[i + 1])).collect();
        let dr: Vec<f64> = (0..n).map(|i| r_edge[i + 1] - r_edge[i]).collect();
        let chord: Vec<f64> = r_mid.iter().map(|&r| blade.chord_at(r)).collect();
        let twist: Vec<f64> = r_mid
            .iter()
            .map(|&r| blade.twist_at(r).to_radians())
            .collect();

        // Local core sizing: scale the base sigma by the local radial spacing
        // relative to the uniform spacing, clamped to [0.5, 1.0] x base so the
        // wake stays well-overlapped where stations are coarse and the tip
        // vortex is not over-smoothed where they are fine.
        let dr_uniform = span / n as f64;
        let base = config.sigma;
        let core_for = |width: f64| -> f32 {
            if config.local_core {
                let s = (width / dr_uniform).clamp(0.5, 1.0) as f32;
                base * s
            } else {
                base
            }
        };
        let sigma_mid: Vec<f32> = dr.iter().map(|&w| core_for(w)).collect();
        let sigma_edge: Vec<f32> = (0..=n)
            .map(|k| {
                // Edge width = mean of adjacent station widths.
                let w = if k == 0 {
                    dr[0]
                } else if k == n {
                    dr[n - 1]
                } else {
                    0.5 * (dr[k - 1] + dr[k])
                };
                core_for(w)
            })
            .collect();

        // Spanwise integral for the aerodynamic feathering spring (see the
        // `feather_span_integral` field).
        let feather_span_integral: f64 =
            (0..n).map(|i| chord[i] * r_mid[i] * r_mid[i] * dr[i]).sum();

        Self {
            polar,
            n_blades: blade.n_blades,
            n_elements: n,
            r_edge,
            r_mid,
            dr,
            chord,
            twist,
            sigma_mid,
            sigma_edge,
            flap: defn.flap.clone(),
            feather: match &defn.pitch_actuation {
                PitchActuation::ServoFlap(act) => Some(act.clone()),
                PitchActuation::DirectMechanical => None,
            },
            cl_alpha: defn.airfoil.CL_alpha_per_rad,
            feather_span_integral,
            ctrl,
            config,
        }
    }

    /// Build the hub-frame flight condition (and the world<->hub kinematics)
    /// from generic `RotorInputs`. Shared by the `AeroModel` load evaluation
    /// (`compute_forces`) and the time-advance (`step`).
    fn flight_condition(&self, inputs: &RotorInputs) -> (FlightCondition, Kinematics) {
        let r_tip = *self
            .r_edge
            .last()
            .expect("r_edge always has n+1 >= 2 edges");
        let kin = kinematics(inputs, inputs.omega_rad_s, r_tip);
        let fc = FlightCondition {
            collective_rad: inputs.collective_rad,
            tilt_lon: inputs.tilt_lon,
            tilt_lat: inputs.tilt_lat,
            // Air velocity relative to the hub, hub axes: in-plane (X, Y) from
            // v_inplane_hub, axial (Z, through-disk) from v_climb.
            v_hub: [kin.v_inplane_hub[0], kin.v_inplane_hub[1], kin.v_climb],
            omega_rad_s: inputs.omega_rad_s,
            rho: inputs.rho_kg_m3,
        };
        (fc, kin)
    }

    /// March the free wake from empty to a periodic state and return
    /// cycle-averaged loads (cold start; no persisted wake).
    /// `dt` is the sub-step duration (seconds) -- use the same value as your
    /// controller loop (e.g. 1.0/400.0). `n_steps` is how many sub-steps to
    /// run in total; loads are averaged over the second half.
    pub fn simulate(&self, fc: &FlightCondition, dt: f64, n_steps: usize) -> VpmRotorResult {
        self.march(fc, None, dt, n_steps).0
    }

    /// Advance the free wake by exactly one sub-step of duration `dt`.
    /// Call this once per controller tick; read `out_state.wake` afterwards
    /// for the current particle cloud.
    pub fn step_one(
        &self,
        fc: &FlightCondition,
        state: &VpmRotorState,
        dt: f64,
    ) -> (VpmRotorResult, VpmRotorState) {
        self.march_window(fc, Some(state), state.psi, dt, 1, 1)
    }

    /// Settle the free wake for `n_steps` sub-steps of duration `dt`, then
    /// return loads averaged over the second half. When `warm` carries a
    /// persisted wake the march continues from it. Typical usage: call once
    /// at startup with enough steps to develop a periodic wake (~several
    /// hundred at 400 Hz for a 30 rad/s rotor), then drive with `step()`.
    pub fn march(
        &self,
        fc: &FlightCondition,
        warm: Option<&VpmRotorState>,
        dt: f64,
        n_steps: usize,
    ) -> (VpmRotorResult, VpmRotorState) {
        let avg_window = (n_steps / 2).max(1);
        self.march_window(fc, warm, 0.0, dt, n_steps, avg_window)
    }
}

/// Persistent state for the VPM rotor, carried between `compute_forces` calls.
///
/// The free-wake particle cloud and the relaxed bound circulation live here as
/// plain struct fields so the wake is *continued* across calls instead of being
/// rebuilt from empty each time. They are deliberately NOT exposed through the
/// `RotorStateExt` inflow vector: `get_inflow` returns an empty `Vec`, so the
/// envelope integrator sees zero inflow DOFs and never tries to damp, blend, or
/// serialize the wake (a particle cloud has no meaningful scalar-inflow tau).
#[derive(Clone, Debug, Default)]
pub struct VpmRotorState {
    /// Persisted free-wake particle cloud. `None` on the first call -> the
    /// march starts from an empty wake.
    pub wake: Option<ParticleField>,
    /// Persisted relaxed bound circulation per (blade, station). Seeds the next
    /// call's implicit circulation solve and its dGamma/dt shed vorticity.
    pub gamma: Option<Vec<Vec<f64>>>,
    /// Current blade-0 azimuth (rad). Accumulated by `step_one` so that each
    /// call sheds particles from the correct rotor position rather than always
    /// restarting at psi=0. `march` ignores this field and always starts at 0.
    pub psi: f64,
    /// Per-blade flap angle (rad), positive = flap up (tip toward -Z, the
    /// thrust side). `None` when flap dynamics are inactive.
    pub beta: Option<Vec<f64>>,
    /// Per-blade flap rate (rad/s). `None` when flap dynamics are inactive.
    pub beta_dot: Option<Vec<f64>>,
    /// Per-blade feathering angle (rad) about the pitch axis, driven by the
    /// servo-flap moment. `None` when the feathering DOF is inactive
    /// (direct-mechanical pitch; no servo-flap actuation).
    pub theta_f: Option<Vec<f64>>,
    /// Per-blade feathering rate (rad/s). `None` when inactive.
    pub theta_f_dot: Option<Vec<f64>>,
}

impl RotorStateExt for VpmRotorState {
    // The wake is persisted via the struct fields above, not the inflow vector,
    // so the inflow DOF count is zero.
    fn get_inflow(&self) -> Vec<f64> {
        Vec::new()
    }
    fn set_inflow(&mut self, arr: Vec<f64>) {
        debug_assert!(arr.is_empty());
    }
}

impl<P: Polar> AeroModel for VpmRotor<P> {
    type State = VpmRotorState;

    fn initial_state(&self) -> Self::State {
        VpmRotorState::default()
    }

    // inflow_taus: the trait default (all-infinity) is correct -- there are no
    // scalar inflow DOFs (the persistent wake lives in the struct fields).

    fn compute_forces(
        &self,
        _inputs: &RotorInputs,
        _state: &VpmRotorState,
    ) -> (AeroResult, VpmRotorState) {
        // VPM is a time-marching free-wake model: loads cannot be produced from
        // a single instantaneous evaluation the way the BEM-family models can.
        // The wake must be advanced in time, and the settling loop is the
        // caller's responsibility -- drive it with `step(dt)` repeatedly (N
        // sub-steps to settle) and read the loads from each `AeroResult`.
        //
        // This is the "assert false on the instantaneous path" contract: there
        // is no correct single-shot compute_forces for a free wake, so calling
        // it is a programming error. (`panic!` rather than `assert!(false, ..)`
        // only to avoid clippy::assertions_on_constants.)
        panic!(
            "VpmRotor::compute_forces is unsupported: a free-wake VPM rotor has \
             no single-shot instantaneous evaluation. Advance the wake with \
             AeroModel::step(dt) N times (the caller owns the settling loop); \
             for a self-contained steady solve use VpmRotor::simulate / march."
        );
    }

    /// Advance the free wake by `dt` seconds and return the loads averaged over
    /// the advanced window.
    ///
    /// This is the only valid time-advance for a VPM rotor: `compute_forces`
    /// panics because a free wake has no instantaneous single-shot evaluation.
    /// The caller settles the wake by calling `step` repeatedly (N sub-steps).
    ///
    /// VPM carries the whole free wake as its state, not scalar inflow DOFs, so
    /// the generic derivative-integration path in the default `AeroModel::step`
    /// (`lam_{n+1} = lam_n + dt*dlam`) is meaningless here and is deliberately
    /// replaced: the wake itself is marched directly. `method` is ignored.
    /// One sub-step of duration `dt` is taken per call -- align `dt` with your
    /// controller loop (e.g. 1/400 s at 400 Hz).
    fn step(
        &self,
        inputs: &RotorInputs,
        state: &VpmRotorState,
        dt: f64,
        _method: IntegrationMethod,
    ) -> (AeroResult, VpmRotorState) {
        // The derivative-integration path must never run for the VPM model:
        // there are no scalar inflow states to integrate. Enforce the invariant
        // that makes this override the only correct time-advance for the wake.
        assert!(
            state.get_inflow().is_empty(),
            "VPM step: derivative-integration path is invalid -- VpmRotorState \
             carries the free wake, not scalar inflow DOFs"
        );
        let (fc, kin) = self.flight_condition(inputs);
        // One sub-step per call; dt IS the sub-step duration.
        let (res, out_state) = self.march_window(&fc, Some(state), state.psi, dt, 1, 1);
        let result = assemble_result(
            res.thrust,
            res.torque,
            res.mx_hub,
            res.my_hub,
            kin.hub_axis,
            &inputs.R_hub,
        );
        (result, out_state)
    }
}

// ============================================================================
// Public API for validation / analysis
// ============================================================================

/// Compute induced velocities at arbitrary points in space from a wake
/// particle field. Used for analysis and validation (e.g., disk inflow sampling).
pub fn induced_velocities_at_points(
    wake: &ParticleField,
    tx: &[f32],
    ty: &[f32],
    tz: &[f32],
) -> Vec<[f32; 3]> {
    induced_at_points(wake, tx, ty, tz)
}
