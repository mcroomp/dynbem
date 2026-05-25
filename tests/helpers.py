"""Shared test utility functions.

Production code has no defaults.  Tests use these helpers to build
standard objects without repeating boilerplate everywhere.
All numeric physics values here are explicit; callers can override any
parameter they care about.
"""
from __future__ import annotations

import math

import numpy as np

from dynbem import (
    LinearPolarParameters,
    BladeGeometry,
    OyeBEMModel,
    PittPetersModel,
    PittPetersRotorState,
    QuasiStaticBEM,
    RotorInputs,
    build_polar,
)
from dynbem.rotor_definition import ControlProperties, RotorDefinition


# ---------------------------------------------------------------------------
# BladeGeometry factory
# ---------------------------------------------------------------------------

def make_blade(
    n_blades,
    radius_m,
    root_cutout_m,
    chord_m,
    twist_deg=0.0,
    n_elements=20,
    r_stations_m=None,
    chord_stations_m=None,
    twist_stations_deg=None,
) -> BladeGeometry:
    """BladeGeometry with explicit required fields; sensible test defaults
    for twist_deg and n_elements."""
    return BladeGeometry(
        n_blades=n_blades,
        radius_m=radius_m,
        root_cutout_m=root_cutout_m,
        chord_m=chord_m,
        twist_deg=twist_deg,
        n_elements=n_elements,
        r_stations_m=r_stations_m,
        chord_stations_m=chord_stations_m,
        twist_stations_deg=twist_stations_deg,
    )


# ---------------------------------------------------------------------------
# LinearPolarParameters factory
# ---------------------------------------------------------------------------

def make_airfoil(
    CL0,
    CL_alpha_per_rad,
    CD0,
    alpha_stall_deg,
    Re_design=None,
) -> LinearPolarParameters:
    """LinearPolarParameters with explicit required fields."""
    return LinearPolarParameters(
        CL0=CL0,
        CL_alpha_per_rad=CL_alpha_per_rad,
        CD0=CD0,
        alpha_stall_deg=alpha_stall_deg,
        Re_design=Re_design,
    )


# ---------------------------------------------------------------------------
# ControlProperties factory
# ---------------------------------------------------------------------------

def make_control(swashplate_pitch_gain_rad, swashplate_phase_deg=None) -> ControlProperties:
    """ControlProperties with required gain; phase defaults to None (0 deg)."""
    return ControlProperties(
        swashplate_pitch_gain_rad=swashplate_pitch_gain_rad,
        swashplate_phase_deg=swashplate_phase_deg,
    )


# ---------------------------------------------------------------------------
# Standard Caradonna-Tung NACA-0012 airfoil parameters
# ---------------------------------------------------------------------------

def ct_airfoil(**overrides) -> LinearPolarParameters:
    """NACA 0012 airfoil (Caradonna-Tung parameters). Override any field."""
    p = dict(
        CL0=0.0,
        CL_alpha_per_rad=2 * math.pi,
        CD0=0.008,
        alpha_stall_deg=15.0,
        Re_design=1_000_000,
    )
    p.update(overrides)
    return make_airfoil(**p)


# ---------------------------------------------------------------------------
# Model builders — all required, no defaults
# ---------------------------------------------------------------------------

def make_bem(defn, n_psi_elements=36) -> QuasiStaticBEM:
    """QuasiStaticBEM with auto-built polar."""
    return QuasiStaticBEM(defn, build_polar(defn.airfoil), n_psi_elements)


def make_pitt_peters(defn, n_psi_elements=36) -> PittPetersModel:
    """PittPetersModel with auto-built polar."""
    return PittPetersModel(defn, build_polar(defn.airfoil), n_psi_elements)


def make_oye(defn, n_psi_elements=36, coupling_k=0.6) -> OyeBEMModel:
    """OyeBEMModel with auto-built polar."""
    return OyeBEMModel(defn, build_polar(defn.airfoil), n_psi_elements, coupling_k)


# ---------------------------------------------------------------------------
# RotorInputs builders
# ---------------------------------------------------------------------------

def hover_inputs(
    collective_deg: float,
    omega_rad_s: float,
    rho_kg_m3: float = 1.225,
    t: float = 0.0,
) -> RotorInputs:
    """Pure hover: no translation, no wind, hub aligned with world NED."""
    return RotorInputs(
        collective_rad=math.radians(collective_deg),
        tilt_lon=0.0,
        tilt_lat=0.0,
        R_hub=np.eye(3),
        v_hub_world=np.zeros(3),
        wind_world=np.zeros(3),
        omega_rad_s=omega_rad_s,
        t=t,
        rho_kg_m3=rho_kg_m3,
    )


def pp_state(lambda_0: float = 0.0, lambda_c: float = 0.0, lambda_s: float = 0.0) -> PittPetersRotorState:
    """PittPetersRotorState with explicit inflow; all three default to 0.0."""
    return PittPetersRotorState(lambda_0, lambda_c, lambda_s)
