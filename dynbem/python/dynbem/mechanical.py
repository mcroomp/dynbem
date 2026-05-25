"""Mechanical (rigid-body) rotor ODE utilities.

The aero models (QuasiStaticBEM, PittPetersModel, OyeBEMModel) are now
pure aerodynamic: they take omega as part of RotorInputs and return only
inflow state derivatives. The caller owns the mechanical ODE.

These helpers implement the standard rigid-body spin-up equation:

    I * d(omega)/dt = motor_torque - Q_aero
    d(spin_angle)/dt = omega

Usage example::

    result, inflow_deriv = aero.compute_forces(inputs, state)
    d_omega = omega_derivative(result.Q_spin, motor_torque_Nm, I_ode_kgm2)
    omega += dt * d_omega
    spin_angle += dt * omega
    inputs.omega_rad_s = omega
"""

__all__ = ["omega_derivative", "euler_step_omega"]


def omega_derivative(Q_aero: float, motor_torque_Nm: float, I_ode_kgm2: float) -> float:
    """Return d(omega)/dt for the rotor rigid-body spin ODE.

    Parameters
    ----------
    Q_aero:
        Aerodynamic reaction torque on the rotor shaft [N.m].
        Positive Q_aero opposes rotation (drag convention): use
        ``AeroResult.Q_spin`` directly (it is already the shaft reaction).
    motor_torque_Nm:
        Applied shaft torque from motor or generator [N.m].
        Positive drives rotation in the direction of omega.
    I_ode_kgm2:
        Rotor polar moment of inertia about the spin axis [kg.m^2].

    Returns
    -------
    float
        d(omega)/dt [rad/s^2].
    """
    return (motor_torque_Nm - Q_aero) / I_ode_kgm2


def euler_step_omega(
    omega: float,
    spin_angle: float,
    Q_aero: float,
    motor_torque_Nm: float,
    I_ode_kgm2: float,
    dt: float,
) -> tuple:
    """Forward-Euler step for the rigid-body spin ODE.

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

    Returns
    -------
    (omega_new, spin_angle_new) : (float, float)
    """
    d_omega = omega_derivative(Q_aero, motor_torque_Nm, I_ode_kgm2)
    omega_new = omega + dt * d_omega
    spin_angle_new = spin_angle + dt * omega
    return omega_new, spin_angle_new
