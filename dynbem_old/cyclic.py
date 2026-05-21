"""Cyclic pitch / swashplate utilities — shared by all aero models.

Convention (see CLAUDE.md "Rotor rotation direction"):
  CCW from above, ψ=0 at +X (hub-frame nose).
  θ_cyclic(ψ) = θ_1c · cos(ψ) + θ_1s · sin(ψ)

Sign convention for the swashplate inputs (helicopter standard):
  tilt_lon > 0  ⇒  nose-down disk (forward stick)
  tilt_lat > 0  ⇒  roll right

Hardware:
  gain  = control.swashplate_pitch_gain_rad   (blade pitch per unit swashplate tilt)
  phi   = control.swashplate_phase_deg        (mechanical phase advance; rad in this module)

Mapping (derived assuming no flap dynamics — blade pitch directly sets local
thrust, so the cyclic→moment phase is 0°, not 90°):
  tilt_lon > 0 → peak pitch at ψ=π   → θ_cyclic ∝ -cos(ψ)
  tilt_lat > 0 → peak pitch at ψ=π/2 → θ_cyclic ∝ +sin(ψ)

After absorbing the swashplate phase rotation:
  θ_1c = gain · (-tilt_lon · cos φ - tilt_lat · sin φ)
  θ_1s = gain · (-tilt_lon · sin φ + tilt_lat · cos φ)

For a future model with full flap-ODE dynamics, set φ ≈ +90° (give or take Lock-
number-dependent fudge) so that the user's tilt_lon/tilt_lat command a *disk*
tilt rather than a thrust asymmetry.
"""
from __future__ import annotations

import math
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from .rotor_definition import ControlProperties


def cyclic_coeffs(
    tilt_lon: float,
    tilt_lat: float,
    control: Optional["ControlProperties"] = None,
) -> tuple[float, float]:
    """Return (θ_1c, θ_1s) such that θ_cyclic(ψ) = θ_1c·cos(ψ) + θ_1s·sin(ψ).

    Defaults gain=1, phase=0 when ``control`` is None — direct blade-pitch
    amplitudes, with helicopter-standard signs.
    """
    if control is None:
        gain = 1.0
        phi = 0.0
    else:
        gain = float(control.swashplate_pitch_gain_rad)
        phi_deg = control.swashplate_phase_deg
        phi = math.radians(float(phi_deg)) if phi_deg is not None else 0.0

    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)
    theta_1c = gain * (-tilt_lon * cos_phi - tilt_lat * sin_phi)
    theta_1s = gain * (-tilt_lon * sin_phi + tilt_lat * cos_phi)
    return theta_1c, theta_1s
