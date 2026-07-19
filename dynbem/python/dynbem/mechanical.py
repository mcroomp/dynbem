"""Mechanical (rigid-body) rotor ODE utilities.

The aero models (QuasiStaticBEM, PittPetersModel, OyeBEMModel) are now
pure aerodynamic: they take omega as part of RotorInputs and return only
inflow state derivatives. The caller owns the mechanical ODE.

These helpers implement the standard rigid-body spin-up equation, with
Coulomb (dry) bearing friction always a first-class parameter (pass
``bearing_friction_Nm=0.0`` for a frictionless bearing):

    I * d(omega)/dt = motor_torque - Q_aero - Q_friction(omega)
    d(spin_angle)/dt = omega

The math is implemented in Rust (``dynbem_rs::mechanical``); this module
is a thin wrapper that also carries the ``spin_angle`` bookkeeping.

This is the single canonical way to advance the mechanical ODE -- there
are exactly two functions: ``omega_derivative`` (the raw ODE right-hand
side, for diagnostics or custom integrators) and ``step_omega`` (the one
recommended integrator, semi-implicit and unconditionally stable in the
aerodynamic term: it never overshoots omega past zero from pure aero
drag). ``step_omega`` also applies a zero-crossing clamp so that Coulomb
friction (or any other purely-decelerating explicit term) cannot flip the
sign of omega within a single step -- the rotor comes to rest instead of
reversing direction.

Usage example::

    result, inflow_deriv = aero.compute_forces(inputs, state)
    omega, spin_angle = step_omega(
        omega, spin_angle, result.Q_spin, motor_torque_Nm, I_ode_kgm2, dt,
    )
    inputs.omega_rad_s = omega
"""

from ._dynbem import (
    omega_derivative as _omega_derivative_rust,
    step_omega as _step_omega_rust,
)

__all__ = ["omega_derivative", "step_omega"]


def omega_derivative(
    omega: float,
    Q_aero: float,
    motor_torque_Nm: float,
    I_ode_kgm2: float,
    bearing_friction_Nm: float = 0.0,
) -> float:
    """Return d(omega)/dt for the rotor rigid-body spin ODE.

    Parameters
    ----------
    omega:
        Current rotor speed [rad/s] (needed only to get the sign of the
        bearing friction torque right).
    Q_aero:
        Aerodynamic reaction torque on the rotor shaft [N.m].
        Positive Q_aero opposes rotation (drag convention): use
        ``AeroResult.Q_spin`` directly (it is already the shaft reaction).
    motor_torque_Nm:
        Applied shaft torque from motor or generator [N.m].
        Positive drives rotation in the direction of omega.
    I_ode_kgm2:
        Rotor polar moment of inertia about the spin axis [kg.m^2].
    bearing_friction_Nm:
        Coulomb (dry) bearing friction torque magnitude [N.m]. Always
        opposes the current rotation direction; exactly zero at
        ``omega == 0`` (no static breakaway threshold modelled). Pass
        ``0.0`` for a frictionless bearing.

    Returns
    -------
    float
        d(omega)/dt [rad/s^2].
    """
    return _omega_derivative_rust(omega, Q_aero, motor_torque_Nm, I_ode_kgm2, bearing_friction_Nm)


def step_omega(
    omega: float,
    spin_angle: float,
    Q_aero: float,
    motor_torque_Nm: float,
    I_ode_kgm2: float,
    dt: float,
    bearing_friction_Nm: float = 0.0,
) -> tuple:
    """Semi-implicit (locally-frozen relaxation) step for the spin ODE.

    This is the single canonical integrator for the mechanical ODE -- use
    it for all production and test time-stepping. Freezes the local
    aerodynamic damping coefficient over the step (the same "freeze the
    coefficient/target over dt" idea used for inflow states' semi-implicit
    integration) and treats it implicitly, so it is unconditionally stable
    in the aerodynamic term: pure aero drag can never send omega past zero
    or oscillate for any ``dt > 0``. Motor torque and Coulomb friction are
    treated explicitly and are subject to a zero-crossing clamp (the rotor
    comes to rest, it doesn't reverse direction). See
    ``dynbem_rs::mechanical::step_omega`` for the full derivation.

    Parameters
    ----------
    omega:
        Current rotor speed [rad/s].
    spin_angle:
        Current rotor azimuth (spin angle) [rad].
    Q_aero:
        Aerodynamic reaction torque [N.m] (from AeroResult.Q_spin).
    motor_torque_Nm:
        Applied shaft torque [N.m].
    I_ode_kgm2:
        Rotor polar moment of inertia [kg.m^2].
    dt:
        Timestep [s].
    bearing_friction_Nm:
        Coulomb (dry) bearing friction torque magnitude [N.m]. Pass
        ``0.0`` for a frictionless bearing.

    Returns
    -------
    (omega_new, spin_angle_new) : (float, float)
    """
    omega_new = _step_omega_rust(
        omega, Q_aero, motor_torque_Nm, I_ode_kgm2, dt, bearing_friction_Nm
    )
    spin_angle_new = spin_angle + dt * omega
    return omega_new, spin_angle_new
