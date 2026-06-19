"""RAWES IC aero sign regression.

This test captures the exact input/output mismatch observed from the windpower
IC warmup diagnostic on 2026-06-19.

Scenario copied from ``E:/repos/windpower``:

* rotor: ``rotors/beaupoil_2026/rotor.yaml``
* NED/FRD convention: ``R_hub[:, 2]`` is body_z, down through the disk
  and toward the tether anchor.
* clean IC position: ``pos=[0.0, 90.41543463617201, -42.72059432582967]``
* wind: ``[0.0, 10.0, 0.0]`` NED, so +East is downwind.
* hub velocity: zero.
* collective: ``-0.18 rad``.
* cyclic: zero for the sign check.

Captured windpower diagnostic outputs for this exact attitude:

QuasiStaticBEM at omega=53.161687 rad/s returned::

    F_world = [0.0, -1029.260, +486.317]
    F_body  = [0.0, 0.0, +1138.368]

Pitt-Peters at omega=38.132161 rad/s returned::

    F_world = [0.0, +481.351, -227.434]
    F_body  = [0.0, 0.0, -532.377]

With this project's convention, body_z points hub-to-anchor/upwind/downward.
A thrust-like rotor force for this IC should be roughly ``-body_z``: downwind
and upward.  Therefore ``F dot body_z`` should be negative and ``F_downwind``
should be positive.  The quasi-static model currently flips that sign.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from dynbem import RotorInputs, create_aero
from dynbem.rotor_definition import load as load_rotor


_ROTOR_YAML = Path(__file__).parent.parent / "rotors" / "beaupoil_2026" / "rotor.yaml"
_RHO = 1.225
_COLLECTIVE_RAD = -0.18
_WIND_NED = np.array([0.0, 10.0, 0.0])
_VEL_NED = np.zeros(3)

# Columns are body axes expressed in NED.  This is the clean wind-aligned
# RAWES IC attitude from windpower's test_generate_ic.py.
_R_RAWES_IC = np.array(
    [
        [0.0, -1.0, -0.0],
        [0.42720594325829603, 0.0, -0.9041543463617201],
        [0.9041543463617201, 0.0, 0.4272059432582967],
    ],
    dtype=float,
)


def _inputs(*, omega_rad_s: float, tilt_lon: float = 0.0, tilt_lat: float = 0.0) -> RotorInputs:
    return RotorInputs(
        collective_rad=_COLLECTIVE_RAD,
        tilt_lon=tilt_lon,
        tilt_lat=tilt_lat,
        R_hub=_R_RAWES_IC,
        v_hub_world=_VEL_NED,
        wind_world=_WIND_NED,
        omega_rad_s=omega_rad_s,
        t=45.0,
        rho_kg_m3=_RHO,
    )


def _force_projections(force_world: np.ndarray) -> dict[str, float]:
    body_z = _R_RAWES_IC[:, 2]
    return {
        "downwind": float(force_world[1]),
        "up": float(-force_world[2]),
        "body_z": float(np.dot(force_world, body_z)),
    }


def test_quasi_static_rawes_ic_force_points_against_body_z() -> None:
    """Quasi-static BEM should not flip RAWES IC thrust along body_z."""
    rotor = load_rotor(str(_ROTOR_YAML))
    model = create_aero(rotor, model="quasi_static")

    result, _ = model.compute_forces(
        _inputs(omega_rad_s=53.161687),
        model.initial_rotor_state(),
    )
    proj = _force_projections(np.asarray(result.F_world, dtype=float))

    assert proj["body_z"] < 0.0, (
        "RAWES IC force sign is flipped for quasi_static BEM.\n"
        f"F_world={np.asarray(result.F_world, dtype=float).tolist()}\n"
        f"F_body={(_R_RAWES_IC.T @ np.asarray(result.F_world, dtype=float)).tolist()}\n"
        f"projections={proj}\n"
        "Expected thrust-like force roughly -body_z: positive downwind/up, "
        "negative F dot body_z."
    )
    assert proj["downwind"] > 0.0
    assert proj["up"] > 0.0


def test_pitt_peters_rawes_ic_force_has_expected_sign() -> None:
    """Dynamic inflow reference for the same RAWES IC sign convention."""
    rotor = load_rotor(str(_ROTOR_YAML))
    model = create_aero(rotor, model="pitt_peters")

    result, _ = model.compute_forces(
        _inputs(omega_rad_s=38.132161),
        model.initial_rotor_state(),
    )
    proj = _force_projections(np.asarray(result.F_world, dtype=float))

    assert proj["body_z"] < 0.0
    assert proj["downwind"] > 0.0
    assert proj["up"] > 0.0


@pytest.mark.parametrize(
    "tilt_lon,tilt_lat",
    [
        (0.0, 0.022616),  # quasi-static trim cyclic captured from windpower
    ],
)
def test_quasi_static_rawes_ic_sign_with_captured_trim_cyclic(
    tilt_lon: float,
    tilt_lat: float,
) -> None:
    """The sign should remain thrust-like with the captured trim cyclic too."""
    rotor = load_rotor(str(_ROTOR_YAML))
    model = create_aero(rotor, model="quasi_static")

    result, _ = model.compute_forces(
        _inputs(omega_rad_s=53.161687, tilt_lon=tilt_lon, tilt_lat=tilt_lat),
        model.initial_rotor_state(),
    )
    proj = _force_projections(np.asarray(result.F_world, dtype=float))

    assert proj["body_z"] < 0.0, (
        f"tilt_lon={tilt_lon:+.6f} tilt_lat={tilt_lat:+.6f} "
        f"F_world={np.asarray(result.F_world, dtype=float).tolist()} "
        f"F_body={(_R_RAWES_IC.T @ np.asarray(result.F_world, dtype=float)).tolist()} "
        f"projections={proj}"
    )