// Rigid-body rotor spin ODE:
//
//     I * d(omega)/dt = motor_torque - Q_aero(omega) - Q_friction(omega)
//
// The aero models (QuasiStaticBEM, PittPetersModel, OyeBEMModel, VpmRotor)
// only ever return Q_aero (AeroResult::Q_spin) for the omega they were
// handed -- integrating omega forward in time is caller-owned (see the
// module doc comment on AeroModel::step in aero_model.rs). This module is
// the SINGLE canonical place that ODE is evaluated and stepped: every
// caller in this repo (production code, benches, tests) must call into
// these functions rather than re-deriving `(motor_torque - Q_spin) / I`
// locally -- that duplication is exactly how the sign/physics bugs this
// module fixed crept in before.
//
// Exactly two functions:
//   - `omega_derivative` -- the canonical ODE right-hand side. Anything
//     that needs d(omega)/dt (diagnostics, custom integrators, tests)
//     calls this instead of reimplementing the formula.
//   - `step_omega` -- the single canonical integrator (semi-implicit,
//     unconditionally stable in the aero term). Use this for all
//     production and test time-stepping; there is no separate
//     plain-forward-Euler stepper to choose between.
//
// Mirrors the Python API in dynbem.mechanical (omega_derivative,
// step_omega).
//
// Bearing friction is modelled as simple Coulomb (dry) friction: a
// constant-magnitude torque that always opposes the current rotation
// direction. It is a first-class parameter of every function here (not a
// separate optional add-on) -- callers that don't want it simply pass
// `bearing_friction_nm = 0.0`.

/// Coulomb (dry) bearing friction torque: constant magnitude, opposing the
/// current rotation direction. Zero at exactly `omega == 0` -- this ODE
/// does not model static breakaway torque (a rotor already at rest with no
/// aero/motor torque driving it stays at rest; see the zero-crossing clamp
/// in the stepper functions below for how a decelerating rotor comes to
/// rest in finite time despite friction being explicitly integrated).
#[inline]
fn friction_torque(omega: f64, bearing_friction_nm: f64) -> f64 {
    if omega > 0.0 {
        bearing_friction_nm
    } else if omega < 0.0 {
        -bearing_friction_nm
    } else {
        0.0
    }
}

/// Clamp a single integration step so that Coulomb friction (or any other
/// purely-decelerating term integrated explicitly) cannot flip the sign of
/// omega within one step. A constant-magnitude decelerating torque will
/// always eventually overshoot past exactly zero in a finite-step
/// integrator; physically the rotor simply comes to rest sometime during
/// that step and stays there (static friction), it does not reverse
/// direction. This clamp is the standard fix for that.
#[inline]
fn clamp_zero_crossing(old_omega: f64, new_omega: f64) -> f64 {
    if old_omega > 0.0 && new_omega < 0.0 {
        0.0
    } else if old_omega < 0.0 && new_omega > 0.0 {
        0.0
    } else {
        new_omega
    }
}

/// d(omega)/dt for the rotor rigid-body spin ODE.
///
/// `omega` is the current rotor speed [rad/s] (needed only to get the sign
/// of the friction torque right). `q_aero` is the aerodynamic reaction
/// torque (`AeroResult::Q_spin`; positive opposes rotation -- use it
/// directly, it is already the shaft reaction). `motor_torque_nm` is the
/// applied shaft torque (positive drives rotation in the direction of
/// omega). `i_ode_kgm2` is the rotor's polar moment of inertia about the
/// spin axis. `bearing_friction_nm` is the Coulomb friction torque
/// magnitude [N.m] (pass `0.0` for a frictionless bearing).
#[inline]
pub fn omega_derivative(
    omega: f64,
    q_aero: f64,
    motor_torque_nm: f64,
    i_ode_kgm2: f64,
    bearing_friction_nm: f64,
) -> f64 {
    let q_friction = friction_torque(omega, bearing_friction_nm);
    (motor_torque_nm - q_aero - q_friction) / i_ode_kgm2
}

/// Semi-implicit (locally-frozen relaxation) step for the spin ODE.
///
/// **This is the single canonical integrator for this ODE.** All
/// production code and tests should step omega through this function --
/// there is no separate plain-forward-Euler stepper to choose between.
///
/// Freezes the local aerodynamic damping coefficient over the step as
/// `tau = I * omega / Q_aero` (the same "freeze the coefficient/target over
/// dt" assumption used by [`crate::aero_model::IntegrationMethod::SemiImplicitEuler`]
/// for inflow states) and treats it implicitly, while motor torque and
/// Coulomb friction (which doesn't scale with omega, so it can't be
/// linearized the same way) are treated explicitly:
///
/// ```text
/// explicit_term = omega + dt * (motor_torque - Q_friction) / I
/// omega_new     = explicit_term / (1 + dt / tau)
/// ```
///
/// This is unconditionally stable in the aero term for any `dt > 0`: when
/// `Q_aero` opposes the current spin direction (the normal no-wind /
/// powered-flight case), the aero contribution alone can never overshoot
/// past zero or oscillate. The explicit friction/motor contribution is
/// subject to a zero-crossing clamp (a rotor decelerating under friction
/// alone comes to rest, it doesn't reverse). Reduces to plain explicit
/// Euler as `dt -> 0` or when there's no meaningful local aero damping to
/// linearize around (`omega == 0`, `Q_aero == 0`, or -- defensively --
/// `Q_aero` and `omega` have opposite signs, e.g. autorotation driving
/// torque, where there is no local *damping* to treat implicitly).
///
/// Windmill / autorotation mode (`Q_aero` and `omega` opposite signs, i.e.
/// the aerodynamic torque is adding power to the shaft rather than
/// resisting it) always falls into that last explicit-Euler branch: `tau`
/// would come out negative there, and it isn't safe to treat a *driving*
/// term implicitly with this frozen-coefficient trick (it would divide by
/// a shrinking-toward-zero-or-negative `1 + dt/tau` and blow up). Plain
/// explicit Euler is only conditionally stable (`dt` less than about twice
/// the true local spin-up timescale), same as any other explicit method,
/// but in a properly-coupled simulation `Q_aero` is recomputed from the
/// real aero model at the new `omega` every step -- so as `omega`
/// approaches the equilibrium where the driving and resisting torques
/// balance, the fallback's own inputs shrink and it settles smoothly. See
/// the `test_step_omega_windmill_*` tests below for multi-step
/// confirmation of this at a realistic fixed timestep.
#[inline]
pub fn step_omega(
    omega: f64,
    q_aero: f64,
    motor_torque_nm: f64,
    i_ode_kgm2: f64,
    dt: f64,
    bearing_friction_nm: f64,
) -> f64 {
    let q_friction = friction_torque(omega, bearing_friction_nm);
    let non_aero_term = omega + dt * (motor_torque_nm - q_friction) / i_ode_kgm2;

    let full_explicit_step = || {
        // Full explicit Euler step, including the aero term -- must not be
        // silently dropped just because there's no local damping
        // coefficient that can be treated implicitly this step.
        omega
            + dt * omega_derivative(
                omega,
                q_aero,
                motor_torque_nm,
                i_ode_kgm2,
                bearing_friction_nm,
            )
    };

    let raw = if omega.abs() < 1e-12 {
        // At omega == 0 there's no rotation to define a local damping
        // timescale from (`tau = I*omega/Q_aero` would be exactly zero, or
        // NaN if Q_aero is also zero) -- but Q_aero can still be nonzero
        // here (e.g. a windmill/autorotation torque acting on a rotor at
        // rest) and must not be dropped from the step.
        full_explicit_step()
    } else {
        let tau = i_ode_kgm2 * omega / q_aero;
        if !tau.is_finite() || tau <= 0.0 {
            // No positive local damping timescale to linearize around
            // (Q_aero == 0, or Q_aero/omega have opposite signs -- e.g. a
            // driving/autorotation torque). Fall back to the *full*
            // explicit Euler step (including the aero term) -- it must
            // not be silently dropped just because it can't be treated
            // implicitly this step.
            full_explicit_step()
        } else {
            non_aero_term / (1.0 + dt / tau)
        }
    };

    clamp_zero_crossing(omega, raw)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_omega_derivative_matches_formula_no_friction() {
        // d_omega = (motor_torque - Q_aero) / I
        assert_eq!(omega_derivative(5.0, 10.0, 0.0, 2.0, 0.0), -5.0);
        assert_eq!(omega_derivative(5.0, 0.0, 4.0, 2.0, 0.0), 2.0);
        assert_eq!(omega_derivative(5.0, 10.0, 10.0, 2.0, 0.0), 0.0);
    }

    #[test]
    fn test_omega_derivative_friction_opposes_positive_omega() {
        // omega > 0: friction subtracts (extra deceleration).
        let d = omega_derivative(5.0, 0.0, 0.0, 1.0, 3.0);
        assert_eq!(d, -3.0);
    }

    #[test]
    fn test_omega_derivative_friction_opposes_negative_omega() {
        // omega < 0: friction acts in +direction (still decelerating
        // toward zero, i.e. opposing the current rotation).
        let d = omega_derivative(-5.0, 0.0, 0.0, 1.0, 3.0);
        assert_eq!(d, 3.0);
    }

    #[test]
    fn test_omega_derivative_friction_zero_at_rest() {
        // omega == 0: no friction torque (no static breakaway modelled).
        let d = omega_derivative(0.0, 0.0, 7.0, 1.0, 3.0);
        assert_eq!(d, 7.0); // pure motor torque, no friction subtracted
    }

    #[test]
    fn test_step_omega_never_overshoots_past_zero() {
        // With Q_aero ~ k*omega^2, a huge dt in plain Euler would send
        // omega deeply negative. Semi-implicit must stay >= 0.
        let omega = 10.0;
        let k = 0.5; // Q_aero = k * omega^2 (quadratic drag model)
        let q_aero = k * omega * omega;
        let i_ode = 1.0;

        for &dt in &[0.01, 0.1, 1.0, 10.0, 1000.0] {
            let new_omega = step_omega(omega, q_aero, 0.0, i_ode, dt, 0.0);
            assert!(
                new_omega >= 0.0,
                "semi-implicit step went negative at dt={dt}: {new_omega}"
            );
            assert!(
                new_omega <= omega,
                "semi-implicit step increased omega at dt={dt}: {omega} -> {new_omega}"
            );
        }
    }

    #[test]
    fn test_step_omega_matches_exact_riccati_solution_for_frozen_k() {
        // For the frozen-coefficient ODE d(omega)/dt = -k*omega^2 (no motor
        // torque, no friction), the exact analytic solution over a step is
        // omega_new = omega / (1 + k*omega*dt). step_omega should
        // reproduce this exactly (it's the same algebraic form).
        let omega = 50.0;
        let k = 0.02;
        let q_aero = k * omega * omega;
        let i_ode = 1.0;
        let dt = 5.0;

        let exact = omega / (1.0 + k * omega * dt);
        let got = step_omega(omega, q_aero, 0.0, i_ode, dt, 0.0);
        assert!(
            (got - exact).abs() < 1e-9,
            "semi-implicit step {got} != exact frozen-k solution {exact}"
        );
    }

    #[test]
    fn test_step_omega_reduces_to_explicit_for_small_dt() {
        // For tiny dt, step_omega should match plain forward Euler on the
        // ODE right-hand side (omega_derivative is the single canonical
        // formula for that -- no separate stepper to compare against).
        let omega = 20.0;
        let q_aero = 3.0;
        let i_ode = 1.0;
        let dt = 1e-6;

        let explicit = omega + dt * omega_derivative(omega, q_aero, 0.0, i_ode, 0.0);
        let semi = step_omega(omega, q_aero, 0.0, i_ode, dt, 0.0);
        assert!(
            (explicit - semi).abs() < 1e-9,
            "step_omega ({semi}) should match explicit Euler ({explicit}) for tiny dt"
        );
    }

    #[test]
    fn test_step_omega_handles_zero_omega_and_zero_torque() {
        // Stopped rotor, no aero torque, no motor torque, no friction:
        // stays at rest.
        assert_eq!(step_omega(0.0, 0.0, 0.0, 1.0, 0.1, 0.0), 0.0);
    }

    #[test]
    fn test_step_omega_falls_back_for_driving_torque() {
        // Autorotation: Q_aero < 0 while omega > 0 (driving, not damping).
        // No positive local damping timescale exists; must fall back to
        // the full explicit-Euler derivative rather than blow up or flip
        // sign from a negative tau.
        let omega = 10.0;
        let q_aero = -2.0;
        let i_ode = 1.0;
        let dt = 0.1;

        let expected = omega + dt * omega_derivative(omega, q_aero, 0.0, i_ode, 0.0);
        let got = step_omega(omega, q_aero, 0.0, i_ode, dt, 0.0);
        assert_eq!(got, expected);
    }

    #[test]
    fn test_step_omega_applies_motor_torque_explicitly() {
        // Pure motor spin-up from rest, zero aero torque, zero friction:
        // reduces to explicit Euler on the forcing term.
        let omega = 0.0;
        let q_aero = 0.0;
        let motor_torque = 5.0;
        let i_ode = 2.0;
        let dt = 0.5;

        let expected = omega + dt * motor_torque / i_ode;
        let got = step_omega(omega, q_aero, motor_torque, i_ode, dt, 0.0);
        assert_eq!(got, expected);
    }

    #[test]
    fn test_step_omega_clamps_zero_crossing_from_friction() {
        // Near-zero omega, aero torque tiny, large friction + dt: must not
        // flip sign.
        let omega = 0.02;
        let q_aero = 1e-6;
        let bearing_friction_nm = 5.0;
        let i_ode = 1.0;
        let dt = 1.0;
        let new_omega = step_omega(omega, q_aero, 0.0, i_ode, dt, bearing_friction_nm);
        assert_eq!(new_omega, 0.0);
    }

    #[test]
    fn test_friction_allows_restart_from_rest_with_motor_torque() {
        // Once clamped to rest, a nonzero motor torque should still be
        // able to spin the rotor back up on the next step (friction is
        // exactly zero at omega == 0 -- no static breakaway threshold
        // modelled).
        let omega = 0.0;
        let motor_torque = 2.0;
        let bearing_friction_nm = 5.0;
        let i_ode = 1.0;
        let dt = 0.1;
        let new_omega = step_omega(omega, 0.0, motor_torque, i_ode, dt, bearing_friction_nm);
        assert!(
            new_omega > 0.0,
            "motor torque should restart rotor from rest"
        );
    }

    #[test]
    fn test_friction_adds_extra_deceleration_beyond_aero_alone() {
        // With friction, decay over a step should be at least as fast as
        // with zero friction (monotone comparison).
        let omega = 20.0;
        let q_aero = 4.0;
        let i_ode = 1.0;
        let dt = 0.1;

        let no_friction = step_omega(omega, q_aero, 0.0, i_ode, dt, 0.0);
        let with_friction = step_omega(omega, q_aero, 0.0, i_ode, dt, 2.0);
        assert!(
            with_friction < no_friction,
            "friction should add extra deceleration: {with_friction} should be < {no_friction}"
        );
    }

    #[test]
    fn test_omega_derivative_windmill_q_aero_adds_power() {
        // Q_aero < 0 while omega > 0 means the aerodynamic torque is
        // driving (adding power to) the shaft rather than resisting it --
        // the classic windmill/autorotation regime. domega/dt must be
        // positive (the rotor speeds up).
        let d = omega_derivative(10.0, -5.0, 0.0, 1.0, 0.0);
        assert!(
            d > 0.0,
            "windmill mode: negative Q_aero should accelerate positive omega, got d={d}"
        );
        assert_eq!(d, 5.0);
    }

    #[test]
    fn test_step_omega_windmill_multistep_converges_without_overshoot() {
        // Synthetic "windmill/autorotation" torque curve
        // Q_aero(omega) = k*(omega - omega_eq): negative (aero ADDS power,
        // spinning the rotor up) below the equilibrium speed, positive
        // (normal aerodynamic drag) above it -- the generic shape of a
        // real autorotation/windmill torque curve (it crosses zero at the
        // equilibrium tip speed) without needing a full BEM aero model. A
        // real caller recomputes Q_aero from the actual aero model at the
        // *new* omega every step -- this test does the same via the
        // closure below, so it exercises step_omega the way production
        // code actually drives it, not just a single frozen-Q_aero step.
        let k = 2.0; // N.m per rad/s of speed error
        let omega_eq = 50.0; // rad/s, equilibrium windmill speed
        let i_ode = 1.0;
        let dt = 0.0025; // 400 Hz -- the fixed timestep this sim uses
        let q_aero_of = |omega: f64| k * (omega - omega_eq);

        let mut omega = 5.0; // start well below equilibrium: aero is driving
        assert!(
            q_aero_of(omega) < 0.0,
            "test setup: must start in the windmill/driving regime"
        );

        let mut prev = omega;
        for _ in 0..4000 {
            omega = step_omega(omega, q_aero_of(omega), 0.0, i_ode, dt, 0.0);
            assert!(
                omega >= prev - 1e-9,
                "windmill spin-up must be monotonically increasing: {prev} -> {omega}"
            );
            assert!(
                omega <= omega_eq + 1e-6,
                "windmill spin-up must not overshoot equilibrium: {omega} > {omega_eq}"
            );
            prev = omega;
        }

        assert!(
            (omega - omega_eq).abs() < 1e-3,
            "windmill spin-up should converge to equilibrium, got {omega}"
        );
    }

    #[test]
    fn test_step_omega_windmill_spinup_from_rest() {
        // omega == 0 exactly hits the `omega.abs() < 1e-12` branch in
        // step_omega (no local damping coefficient to freeze at rest), so
        // the very first step is explicit -- confirm the rotor still
        // starts moving in the windmill regime and keeps converging.
        let k = 2.0;
        let omega_eq = 50.0;
        let i_ode = 1.0;
        let dt = 0.0025;
        let q_aero_of = |omega: f64| k * (omega - omega_eq);

        let mut omega = 0.0;
        let mut prev = omega;
        for _ in 0..4000 {
            omega = step_omega(omega, q_aero_of(omega), 0.0, i_ode, dt, 0.0);
            assert!(
                omega >= prev - 1e-9,
                "must not decelerate while windmilling from rest: {prev} -> {omega}"
            );
            prev = omega;
        }

        assert!(
            omega > 0.0,
            "windmill torque should spin the rotor up from rest"
        );
        assert!(
            (omega - omega_eq).abs() < 1e-3,
            "should converge near equilibrium, got {omega}"
        );
    }

    #[test]
    fn test_step_omega_windmill_equilibrium_shifts_with_friction() {
        // With Coulomb friction opposing rotation, the windmill can only
        // spin up to the point where the aero driving torque exactly
        // balances friction (Q_aero == -bearing_friction), which is BELOW
        // the frictionless zero-crossing equilibrium `omega_eq`.
        let k = 2.0;
        let omega_eq = 50.0;
        let bearing_friction_nm = 4.0;
        let i_ode = 1.0;
        let dt = 0.0025;
        let q_aero_of = |omega: f64| k * (omega - omega_eq);

        let mut omega = 5.0;
        let mut prev = omega;
        for _ in 0..8000 {
            omega = step_omega(omega, q_aero_of(omega), 0.0, i_ode, dt, bearing_friction_nm);
            assert!(
                omega >= prev - 1e-9,
                "windmill spin-up with friction must still be monotonic: {prev} -> {omega}"
            );
            prev = omega;
        }

        let expected_equilibrium = omega_eq - bearing_friction_nm / k;
        assert!(
            (omega - expected_equilibrium).abs() < 1e-2,
            "friction should shift the windmill equilibrium down to {expected_equilibrium}, got {omega}"
        );
        assert!(
            omega < omega_eq,
            "windmill equilibrium with friction must be below the frictionless equilibrium"
        );
    }

    #[test]
    fn test_step_omega_continuous_near_q_aero_zero_crossing() {
        // Directly probes continuity of step_omega as Q_aero crosses zero
        // at fixed omega. This is exactly the boundary between the
        // semi-implicit "damping" branch (tau valid) and the explicit
        // "driving/fallback" branch (tau invalid) -- a realistic worry
        // given this repo has previously had windmill/helicopter branch
        // discontinuity bugs in the BEM solver itself (see the
        // lambda_climb boundary fix). A tiny change in Q_aero must not
        // produce a large jump in the result.
        let omega = 35.0; // rad/s, realistic autorotation-range rotor speed
        let i_ode = 5.0;
        let dt = 0.0025; // 400 Hz
        let eps = 1e-4;

        let just_above = step_omega(omega, eps, 0.0, i_ode, dt, 0.0); // damping branch
        let at_zero = step_omega(omega, 0.0, 0.0, i_ode, dt, 0.0); // fallback branch
        let just_below = step_omega(omega, -eps, 0.0, i_ode, dt, 0.0); // fallback branch

        assert!(
            (just_above - at_zero).abs() < 1e-6,
            "branch switch at Q_aero=0 should be continuous: {just_above} vs {at_zero}"
        );
        assert!(
            (at_zero - just_below).abs() < 1e-6,
            "fallback branch should be continuous in Q_aero near zero: {at_zero} vs {just_below}"
        );
    }

    #[test]
    fn test_step_omega_no_chattering_from_sign_flipping_aero_torque_near_equilibrium() {
        // Realistic corner case: near the autorotative equilibrium, gusts
        // or BEM-solver noise can flip the sign of Q_aero step-to-step
        // even though omega itself barely moves. The branch switch
        // between semi-implicit and explicit-fallback must not cause
        // chattering or a blow-up.
        let i_ode = 5.0;
        let dt = 0.0025;
        let q_sequence = [0.05, -0.05, 0.05, -0.05, 0.02, -0.02, 0.0, 0.05, -0.05];

        let mut omega = 35.0;
        for &q in &q_sequence {
            let next = step_omega(omega, q, 0.0, i_ode, dt, 0.0);
            assert!(
                next.is_finite(),
                "chattering aero torque must not produce a non-finite result: {next}"
            );
            assert!(
                (next - omega).abs() < 0.01,
                "chattering aero torque should not cause a large single-step jump: {omega} -> {next}"
            );
            omega = next;
        }
    }

    #[test]
    fn test_step_omega_autorotation_entry_crosses_damping_to_driving_branch_smoothly() {
        // Realistic engine-failure-in-forward-flight scenario: the rotor
        // starts ABOVE the true (zero-torque) autorotative equilibrium
        // speed, so it is initially decelerated by ordinary aerodynamic
        // drag (Q_aero > 0, the semi-implicit "damping" branch). As it
        // slows, Q_aero crosses zero and goes negative (the windmill
        // torque now drives the rotor), landing in the explicit-fallback
        // "driving" branch, and it settles at the friction-shifted
        // equilibrium below the zero-torque point. Unlike the windmill
        // tests above (which start already in the driving branch), this
        // exercises the actual mid-simulation branch switch.
        let i_ode = 5.0;
        let c = 20.0; // N.m per rad/s of speed above/below the autorotative eq.
        let omega_auto = 35.0; // rad/s, true (zero-torque) autorotation speed
        let bearing_friction_nm = 2.0;
        let dt = 0.0025; // 400 Hz
        let q_aero_of = |omega: f64| c * (omega - omega_auto);

        let mut omega = 45.0; // powered-flight NR, above eq.: engine just failed
        assert!(
            q_aero_of(omega) > 0.0,
            "test setup: must start in the damping regime"
        );

        let mut prev = omega;
        let mut crossed_into_driving_branch = false;
        for _ in 0..4000 {
            let q = q_aero_of(omega);
            if q < 0.0 {
                crossed_into_driving_branch = true;
            }
            omega = step_omega(omega, q, 0.0, i_ode, dt, bearing_friction_nm);
            assert!(
                omega <= prev + 1e-9,
                "autorotation entry must decelerate monotonically: {prev} -> {omega}"
            );
            assert!(
                (prev - omega).abs() < 1.0,
                "must not take an unreasonably large single-step jump across the branch boundary: {prev} -> {omega}"
            );
            prev = omega;
        }

        assert!(
            crossed_into_driving_branch,
            "test setup: simulation should reach the windmill/driving regime below omega_auto"
        );

        let expected_equilibrium = omega_auto - bearing_friction_nm / c;
        assert!(
            (omega - expected_equilibrium).abs() < 1e-2,
            "should settle at the friction-shifted autorotative equilibrium {expected_equilibrium}, got {omega}"
        );
    }

    #[test]
    fn test_step_omega_settles_exactly_at_rest_and_stays_there() {
        // Realistic rotor-shutdown scenario: motor off, negligible
        // residual aero torque, only bearing friction decelerating the
        // rotor. Once it reaches exactly rest it must stay there
        // indefinitely -- no creeping negative from floating-point
        // residue, and no spurious restart from the omega==0 special
        // case.
        let i_ode = 5.0;
        let bearing_friction_nm = 3.0;
        let dt = 0.0025;
        let mut omega = 2.0; // rad/s, rotor coasting down after shutdown

        for _ in 0..2000 {
            omega = step_omega(omega, 0.0, 0.0, i_ode, dt, bearing_friction_nm);
        }

        assert_eq!(
            omega, 0.0,
            "rotor should come to rest and stay there, got {omega}"
        );
    }
}
