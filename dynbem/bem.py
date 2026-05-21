"""Level 1 BEM solver — multi-element, NED frame.

Coordinate system: NED (North-East-Down). Rotation: CCW from above
(American helicopter convention; see CLAUDE.md).
  - Rotor hub axis points in the +Z (down) direction for a level hover.
  - Thrust opposes gravity: F_world[2] < 0 (upward force = negative Z).
  - v_climb > 0: axial freestream flows downward through disk (helicopter climb).
  - v_climb < 0: axial freestream flows upward through disk (autorotation/turbine).
  - v_climb = 0: hover — rotor generates its own induced inflow.

Inflow iteration uses the total inflow ratio λ_r = v_a/(Ω·R) rather than the
wind-turbine induction-factor form (which degenerates at v_climb = 0 / hover).

Momentum-BEM derivation
-----------------------
Momentum annulus (with combined Prandtl tip + hub loss F = F_tip · F_hub):
    dCT/dx = 4·F·x·λ_r·(λ_r − λ_c)

Blade element (N blades, chord c, local solidity σ_r = N·c/(2π·r)):
    dCT/dx = σ_r·x·cn·(λ_r² + x²)

Cancel x, define k = σ_r·cn/(4·F):
    k·(λ_r² + x²) = λ_r·(λ_r − λ_c)   ← the iteration equation
    (k−1)·λ_r² + λ_c·λ_r + k·x² = 0   ← quadratic form

Forward flight / cyclic
-----------------------
When edgewise advance ratio mu > 0.01 OR cyclic input is nonzero, the
model integrates over n_psi azimuth stations. At each ψ:
  - Per-azimuth blade pitch θ(ψ) = collective + θ_1c·cos ψ + θ_1s·sin ψ
    (cyclic coefficients via dynbem.cyclic.cyclic_coeffs).
  - In-plane wind projected onto t_hat for v_t_extra.
  - In-plane hub moments Mx_hub, My_hub accumulated for AeroResult.M_orbital.

Integration contract
--------------------
compute_forces() returns (AeroResult, derivative) where derivative holds dstate/dt.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import NamedTuple

import numpy as np

from . import AeroBase, AeroResult, RotorInputs
from .cyclic import cyclic_coeffs
from .rotor_definition import RotorDefinition
from .rotor_state import QuasiStaticRotorState, RotorState
from .polar import AirfoilPolar, LinearPolar, TabulatedPolar


# ---------------------------------------------------------------------------
# Prandtl tip / hub losses (public — used in component tests)
# ---------------------------------------------------------------------------

def prandtl_tip_loss(n_blades: int, x: float, phi_rad: float) -> float:
    """Prandtl tip-loss factor F_tip at normalised radius x = r/R.

    Returns F ∈ (0, 1].  F → 1 far from tip; F → 0 at the tip for small phi.

    n_blades  number of blades
    x         r/R, normalised radius ∈ (0, 1]
    phi_rad   flow angle from rotor plane (rad); may be negative in turbine mode
    """
    if abs(phi_rad) < 1e-9 or x >= 1.0:
        return 1.0
    f = (n_blades / 2.0) * (1.0 - x) / (x * abs(math.sin(phi_rad)))
    return (2.0 / math.pi) * math.acos(min(1.0, math.exp(-f)))


def prandtl_hub_loss(n_blades: int, x: float, x_hub: float, phi_rad: float) -> float:
    """Prandtl hub-loss factor F_hub at normalised radius x = r/R.

    Accounts for root-vortex losses near the hub cutout, mirroring the
    tip-loss correction. Form follows Glauert/Prandtl with the hub
    radius substituted for the tip radius:

        f_hub = (N/2) · (x − x_hub) / (x_hub · |sin φ|)
        F_hub = (2/π) · arccos(exp(−f_hub))

    Returns F → 1 far from the hub, F → 0 at the hub cutout.

    n_blades  number of blades
    x         r/R, normalised radius (must satisfy x > x_hub)
    x_hub     root_cutout / R, normalised hub radius ∈ (0, 1)
    phi_rad   flow angle from rotor plane (rad); may be negative in turbine mode
    """
    if abs(phi_rad) < 1e-9 or x <= x_hub or x_hub <= 0.0:
        return 1.0
    f = (n_blades / 2.0) * (x - x_hub) / (x_hub * abs(math.sin(phi_rad)))
    return (2.0 / math.pi) * math.acos(min(1.0, math.exp(-f)))


# ---------------------------------------------------------------------------
# BEM element (public — exposed for unit testing)
# ---------------------------------------------------------------------------

_MAX_BEM_ITER = 60
_BEM_TOL = 1e-7


class BEMElementResult(NamedTuple):
    """Converged state and forces for one blade-element annulus."""
    lambda_r: float   # total axial inflow ratio v_a/(Ω·R)
    a_prime: float    # tangential induction factor
    dT: float         # thrust contribution (N) — positive = upward for level rotor
    dQ: float         # reaction torque contribution (N·m) — positive = drag on rotor
    momentum_residual: float  # |4F·λ_r·(λ_r−λ_c) − σ_r·cn·(λ_r²+x²)| at convergence


def _solve_bem_element_windmill(
    r: float,
    dr: float,
    chord: float,
    twist_rad: float,
    collective_rad: float,
    omega: float,
    v_climb: float,
    rho: float,
    n_blades: int,
    radius_m: float,
    polar: AirfoilPolar,
    use_tip_loss: bool,
    root_cutout_m: float,
) -> BEMElementResult | None:
    """Wind-turbine BEM iteration for axial upflow (v_climb < 0).

    Solves the standard windmill momentum balance
        a / (1 - a) = sigma_r * Cn / (4 F sin^2(phi))
    with helicopter-convention angles throughout (positive
    collective+twist = pitch-to-stall, phi < 0 in upflow).  Returns
    `None` when the iteration does not converge to a valid windmill
    state (0 < a < 0.5 and AoA below stall) -- in that case the caller
    should fall back to the helicopter momentum quadratic.

    "Valid windmill state" means the rotor extracts axial momentum
    from the wind (induction factor a > 0, axial flow at disk is in
    same direction as freestream but reduced in magnitude).  Excludes
    propeller mode (a < 0, rotor accelerates the flow) and deep WBS
    (a > 0.5, momentum theory breaks down -- Glauert/Buhl correction
    would be needed but is out of scope for this branch).
    """
    if v_climb >= -1e-9:
        return None  # axial windmill regime requires upflow

    U = -v_climb  # positive axial freestream speed
    Omega_R = omega * radius_m
    if Omega_R < 1e-6:
        return None
    x = r / radius_m
    x_hub = root_cutout_m / radius_m if radius_m > 0.0 else 0.0
    sigma_r = n_blades * chord / (2.0 * math.pi * r)
    theta = collective_rad + twist_rad  # helicopter convention

    a = 0.3       # typical windmill starting point
    a_prime = 0.0
    F_final = 1.0

    for _ in range(_MAX_BEM_ITER):
        # NED conventions: v_a points in +Z direction (downward).  Axial
        # wind is upflow, so v_a at disk is negative.  (1 - a) is the
        # induction-reduced axial speed; we apply the NED sign by
        # multiplying by -U.
        v_a = -(1.0 - a) * U
        v_t = (1.0 + a_prime) * omega * r
        if v_t < 1e-9:
            return None

        phi = math.atan2(v_a, v_t)  # negative for windmill (v_a < 0)
        alpha = theta - phi          # helicopter convention
        cl, cd = polar.cl_cd(alpha)

        # No pre-stall AoA guard: the windmill momentum equation depends
        # only on Cn (not on dCn/dalpha), so as long as Cn > 0 and the
        # iteration converges, the polar's post-stall plateau is fine.
        # Helicopter autorotation operates routinely above stall (high
        # AoA, high Cl, high Cd) and we want the windmill solver to own
        # those cases too.

        if use_tip_loss:
            F = (prandtl_tip_loss(n_blades, x, phi)
                 * prandtl_hub_loss(n_blades, x, x_hub, phi))
        else:
            F = 1.0
        F = max(F, 1e-4)
        F_final = F

        # Helicopter Cn/Ct convention.  At negative phi, cos(phi) =
        # cos(|phi|) and sin(phi) = -sin(|phi|), so the formula matches
        # the textbook turbine `Cn = Cl cos|phi| + Cd sin|phi|` exactly.
        cn = cl * math.cos(phi) - cd * math.sin(phi)
        ct = cl * math.sin(phi) + cd * math.cos(phi)

        sin2 = math.sin(phi) ** 2
        if cn <= 1e-9 or sin2 < 1e-12:
            return None

        # Classical windmill momentum for a <= 0.4 reduces algebraically
        # to:
        #   a / (1 - a) = sigma_r Cn / (4 F sin^2 phi)
        # (the (1-a)^2 factor on the blade-element side and the
        # 4F a (1-a) on the momentum side cancel one (1-a) cleanly).
        denom_a = 4.0 * F * sin2 / (sigma_r * cn) + 1.0
        a_new = 1.0 / denom_a

        if a_new > 0.4:
            # Buhl (NREL/TP-500-36834 Eq 18) replaces the classical
            # 4F a (1-a) momentum side with an empirical quadratic that
            # bridges 0.4 < a < ~1.  The blade-element side stays
            # C_T_BE = sigma_r Cn (1-a)^2 / sin^2(phi), so at
            # convergence:
            #   8/9 + (4F - 40/9) a + (50/9 - 4F) a^2 = K (1-a)^2
            # with K = sigma_r Cn / sin^2(phi).  Rearranged in a:
            #   (50/9 - 4F - K) a^2 + (4F - 40/9 + 2K) a + (8/9 - K) = 0
            # We use the current iterate's K (it depends on phi which
            # depends on a, so the whole thing is a fixed-point).
            #
            # The Glauert / Buhl pieces match at a = 0.4 in both value
            # and first derivative, so transitions between the
            # classical branch above and the Buhl branch here are
            # smooth.
            K = sigma_r * cn / sin2
            A = 50.0 / 9.0 - 4.0 * F - K
            B = 4.0 * F - 40.0 / 9.0 + 2.0 * K
            C = 8.0 / 9.0 - K
            disc = B * B - 4.0 * A * C
            if disc < 0.0 or abs(A) < 1e-12:
                return None
            sq = math.sqrt(disc)
            a_lower = (-B - sq) / (2.0 * A)
            a_upper = (-B + sq) / (2.0 * A)
            # Pick the physical root in [0.4, 1.0].
            candidates = [
                cand for cand in (a_lower, a_upper)
                if 0.4 <= cand <= 1.0
            ]
            if not candidates:
                return None
            # If both roots are in range, take the smaller (closer to
            # the classical branch).
            a_new = min(candidates)

        # Tangential induction: a' / (1+a') = sigma_r Ct / (4 F sin phi cos phi)
        sin_cos = math.sin(phi) * math.cos(phi)
        if abs(ct) > 1e-9 and abs(sin_cos) > 1e-9:
            denom_ap = 4.0 * F * sin_cos / (sigma_r * ct) - 1.0
            a_prime_new = 1.0 / denom_ap if abs(denom_ap) > 1e-8 else 0.0
        else:
            a_prime_new = 0.0

        # Buhl is calibrated up to a ~ 1; beyond that (flow reversal at
        # the disk) the helicopter quadratic's empirical fit is the
        # right tool.
        if not (0.0 <= a_new <= 1.0):
            return None
        a_prime_new = max(-0.5, min(0.5, a_prime_new))

        converged = (abs(a_new - a) < _BEM_TOL
                     and abs(a_prime_new - a_prime) < _BEM_TOL)
        a = 0.5 * a + 0.5 * a_new
        a_prime = 0.5 * a_prime + 0.5 * a_prime_new
        if converged:
            break
    else:
        return None  # did not converge

    # Final state and forces.
    v_a_f = -(1.0 - a) * U
    v_t_f = (1.0 + a_prime) * omega * r
    v_rel = math.sqrt(v_a_f ** 2 + v_t_f ** 2)
    phi_f = math.atan2(v_a_f, v_t_f)
    alpha_f = theta - phi_f
    cl_f, cd_f = polar.cl_cd(alpha_f)
    cn_f = cl_f * math.cos(phi_f) - cd_f * math.sin(phi_f)
    ct_f = cl_f * math.sin(phi_f) + cd_f * math.cos(phi_f)
    if cn_f <= 0.0:
        return None

    q_dyn = 0.5 * rho * v_rel ** 2 * chord * dr * n_blades
    # dT, dQ in helicopter convention -- exactly the same expressions
    # the helicopter branch uses.  Sign emerges naturally from cn_f /
    # ct_f at phi < 0:
    #   cn_f = cl cos(phi) - cd sin(phi)  -> positive in upflow
    #     (lift normal to rotor plane, in -NED-Z = same direction as
    #     reaction thrust the rotor feels from the wind).
    #   ct_f = cl sin(phi) + cd cos(phi)  -> negative at phi < 0 when
    #     lift dominates: lift's tangential component drives the rotor
    #     (extraction), so dQ < 0 in dynbem convention as expected.
    dT = q_dyn * cn_f
    dQ = q_dyn * ct_f * r

    # Map a -> lambda_r in dynbem NED convention.
    lambda_r = v_a_f / Omega_R
    # Track momentum residual for diagnostics (always ~0 at windmill
    # convergence by construction).
    _ = F_final  # tip-loss factor preserved only for diagnostic prints
    return BEMElementResult(lambda_r, a_prime, dT, dQ, 0.0)


def solve_bem_element(
    r: float,
    dr: float,
    chord: float,
    twist_rad: float,
    collective_rad: float,
    omega: float,
    v_climb: float,
    rho: float,
    n_blades: int,
    radius_m: float,
    polar: AirfoilPolar,
    use_tip_loss: bool,
    v_t_extra: float = 0.0,
    root_cutout_m: float = 0.0,
) -> BEMElementResult:
    """Solve one blade-element annulus; return converged state and forces.

    v_climb   External axial freestream in the hub-axis direction (NED):
                > 0  air from above (helicopter climb)
                = 0  hover
                < 0  air from below (autorotation / flying wind turbine)
    v_t_extra In-plane wind tangential contribution at this azimuth (m/s).
              Positive adds to blade tangential velocity (advancing side).
              Zero (default) gives the standard axial-flight result.

    dT > 0: thrust opposing inflow (upward for level rotor = −Z in NED).
    dQ > 0: reaction torque opposing rotor spin.
    """
    x = r / radius_m
    x_hub = root_cutout_m / radius_m if radius_m > 0.0 else 0.0
    Omega_R = omega * radius_m
    if Omega_R < 1e-6:
        return BEMElementResult(0.0, 0.0, 0.0, 0.0, 0.0)

    sigma_r = n_blades * chord / (2.0 * math.pi * r)
    theta = collective_rad + twist_rad
    lambda_c = v_climb / Omega_R

    if lambda_c >= 0.0:
        lambda_r = max(lambda_c + 0.03, 0.02)
    else:
        lambda_r = min(lambda_c * 0.85, -0.02)

    a_prime = 0.0
    F_final = 1.0

    for _ in range(_MAX_BEM_ITER):
        v_a = lambda_r * Omega_R
        v_t = omega * r * (1.0 + a_prime) + v_t_extra
        if v_t < 1e-9:
            break

        phi = math.atan2(v_a, v_t)
        alpha = theta - phi
        cl, cd = polar.cl_cd(alpha)

        if use_tip_loss:
            F = (
                prandtl_tip_loss(n_blades, x, phi)
                * prandtl_hub_loss(n_blades, x, x_hub, phi)
            )
        else:
            F = 1.0
        F = max(F, 1e-4)
        F_final = F

        cn = cl * math.cos(phi) - cd * math.sin(phi)
        ct = cl * math.sin(phi) + cd * math.cos(phi)

        # Momentum-BEM quadratic: k = σ_r·cn/(4·F)
        k = sigma_r * cn / (4.0 * F)

        if abs(k - 1.0) > 1e-6:
            disc = max(0.0, lambda_c**2 - 4.0 * (k - 1.0) * k * x**2)
            sq = math.sqrt(disc)
            denom = 2.0 * (k - 1.0)
            r1 = (-lambda_c + sq) / denom
            r2 = (-lambda_c - sq) / denom
            if lambda_c >= 0.0:
                lambda_r_new = r2 if r2 > 0.0 else r1
            else:
                lambda_r_new = r1 if r1 < 0.0 else r2
        else:
            if abs(lambda_c) > 1e-8:
                lambda_r_new = -k * x**2 / lambda_c
            else:
                lambda_r_new = x * math.sqrt(max(0.0, k))

        lambda_r_new = max(-2.0, min(2.0, lambda_r_new))

        sc = math.sin(phi) * math.cos(phi)
        if abs(sc) > 1e-8 and abs(ct) > 1e-10:
            ap_denom = 4.0 * F * sc / (sigma_r * ct) - 1.0
            a_prime_new = (1.0 / ap_denom) if abs(ap_denom) > 1e-8 else 0.0
            a_prime_new = max(-0.5, min(0.5, a_prime_new))
        else:
            a_prime_new = 0.0

        converged = (
            abs(lambda_r_new - lambda_r) < _BEM_TOL
            and abs(a_prime_new - a_prime) < _BEM_TOL
        )
        lambda_r = 0.5 * lambda_r + 0.5 * lambda_r_new
        a_prime = 0.5 * a_prime + 0.5 * a_prime_new
        if converged:
            break

    # Recompute final forces at converged state
    v_a_f = lambda_r * Omega_R
    v_t_f = omega * r * (1.0 + a_prime) + v_t_extra
    v_rel = math.sqrt(v_a_f**2 + v_t_f**2)
    phi_f = math.atan2(v_a_f, v_t_f)
    alpha_f = theta - phi_f
    cl_f, cd_f = polar.cl_cd(alpha_f)
    cn_f = cl_f * math.cos(phi_f) - cd_f * math.sin(phi_f)
    ct_f = cl_f * math.sin(phi_f) + cd_f * math.cos(phi_f)

    q_dyn = 0.5 * rho * v_rel**2 * chord * dr * n_blades
    dT = q_dyn * cn_f
    dQ = q_dyn * ct_f * r

    # Momentum-BEM balance residual — should be near zero at convergence.
    # 4·F·λ_r·(λ_r−λ_c) == σ_r·cn·(λ_r²+x²)
    momentum_side = 4.0 * F_final * lambda_r * (lambda_r - lambda_c)
    blade_side = sigma_r * cn_f * (lambda_r**2 + x**2)
    residual = abs(momentum_side - blade_side)

    return BEMElementResult(lambda_r, a_prime, dT, dQ, residual)


# ---------------------------------------------------------------------------
# BEM model
# ---------------------------------------------------------------------------

@dataclass
class BEMModel(AeroBase):
    """Multi-element BEM rotor model (Level 1).

    NED frame throughout.  Returns QuasiStaticRotorState derivatives.

    n_psi_elements controls azimuth discretisation for forward flight
    (mu > 0.01).  36 points (10 deg steps) is accurate to ~1% for mu < 0.4.
    Set higher for smoother torque curves; lower for faster envelope sweeps.
    """

    defn: RotorDefinition
    n_psi_elements: int = 36
    _polar: AirfoilPolar = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if self.defn.airfoil.polar_csv is not None:
            self._polar = TabulatedPolar.from_csv(self.defn.airfoil.polar_csv)
        else:
            self._polar = LinearPolar.from_properties(self.defn.airfoil)

    def initial_rotor_state(self) -> QuasiStaticRotorState:
        return QuasiStaticRotorState()

    def compute_forces(
        self,
        inputs: RotorInputs,
        state: RotorState,
    ) -> tuple[AeroResult, RotorState]:
        assert isinstance(state, QuasiStaticRotorState)

        blade = self.defn.blade
        airfoil = self.defn.airfoil
        omega = state.omega_rad_s

        hub_axis_ned = inputs.R_hub @ np.array([0.0, 0.0, 1.0])
        v_rel_world = inputs.wind_world - inputs.v_hub_world
        v_climb = float(np.dot(v_rel_world, hub_axis_ned))

        # In-plane (edgewise) wind component and advance ratio
        v_inplane = v_rel_world - v_climb * hub_axis_ned
        v_edge = float(np.linalg.norm(v_inplane))
        Omega_R = max(omega * blade.radius_m, 1e-6)
        mu = v_edge / Omega_R

        n = blade.n_elements
        r_root, r_tip = blade.root_cutout_m, blade.radius_m
        dr = (r_tip - r_root) / n
        r_stations = np.linspace(r_root + 0.5 * dr, r_tip - 0.5 * dr, n)
        # Per-station chord / twist: interpolated from blade.r_stations_m
        # if the YAML provides them, otherwise falls back to the scalar
        # blade.chord_m / blade.twist_deg.  See BladeGeometry.chord_at.
        chord_per_station = np.array([blade.chord_at(float(r)) for r in r_stations])
        twist_per_station_rad = np.array(
            [math.radians(blade.twist_at(float(r))) for r in r_stations])

        T_total = 0.0
        Q_total = 0.0
        # In-plane hub-frame moments from non-axisymmetric thrust
        # (forward flight and/or cyclic).  Per-element contribution:
        #   dM_hub = r·dT·[sin(ψ), cos(ψ), 0]
        Mx_hub = 0.0
        My_hub = 0.0

        # Cyclic pitch → θ_cyclic(ψ) = θ_1c·cos(ψ) + θ_1s·sin(ψ).
        theta_1c, theta_1s = cyclic_coeffs(
            inputs.tilt_lon, inputs.tilt_lat, self.defn.control
        )
        has_cyclic = abs(theta_1c) + abs(theta_1s) > 1e-12

        if (mu > 0.01 or has_cyclic) and omega > 1.0:
            # Forward flight and/or cyclic: average blade forces over n_psi
            # azimuth stations.  At each ψ, project the in-plane wind onto the
            # blade tangential direction to get v_t_extra, add cyclic to the
            # collective, then solve the full BEM element.  Retreating-side
            # elements where total tangential velocity reverses are skipped
            # (standard BEM approximation for mu < 0.5).
            n_psi = self.n_psi_elements
            for i_psi in range(n_psi):
                psi = 2.0 * math.pi * i_psi / n_psi
                sin_psi = math.sin(psi)
                cos_psi = math.cos(psi)
                # CCW-from-above rotation (American convention), ψ=0 at +X.
                # See CLAUDE.md "Rotor rotation direction".
                # Tangential direction (blade tip's motion) in hub frame:
                #   t_hat = [-sin(ψ), -cos(ψ), 0]
                t_hat_ned = inputs.R_hub @ np.array(
                    [-sin_psi, -cos_psi, 0.0]
                )
                # v_t = ω·r − v_inplane · t_hat, so v_t_extra = −v_inplane·t_hat.
                v_t_extra = -float(np.dot(v_inplane, t_hat_ned))

                col_psi = inputs.collective_rad + theta_1c * cos_psi + theta_1s * sin_psi

                for i_r, r in enumerate(r_stations):
                    # Skip reverse-flow region
                    if omega * float(r) + v_t_extra <= 0.0:
                        continue
                    elem = solve_bem_element(
                        r=float(r), dr=dr,
                        chord=float(chord_per_station[i_r]),
                        twist_rad=float(twist_per_station_rad[i_r]),
                        collective_rad=col_psi,
                        omega=omega, v_climb=v_climb,
                        rho=inputs.rho_kg_m3,
                        n_blades=blade.n_blades, radius_m=blade.radius_m,
                        polar=self._polar, use_tip_loss=airfoil.tip_loss,
                        v_t_extra=v_t_extra,
                        root_cutout_m=blade.root_cutout_m,
                    )
                    T_total += elem.dT
                    Q_total += elem.dQ
                    Mx_hub  += float(r) * elem.dT * sin_psi
                    My_hub  += float(r) * elem.dT * cos_psi

            # Average over azimuth.
            T_total /= n_psi
            Q_total /= n_psi
            Mx_hub  /= n_psi
            My_hub  /= n_psi

        else:
            # Axial flight (hover / climb / descent / windmill), no
            # cyclic: axisymmetric.  When v_climb < 0 (wind blows
            # axially through the disk) try the wind-turbine windmill
            # iteration first -- it has a valid root in the
            # energy-extraction regime that the helicopter momentum
            # quadratic does not.  Fall back to the helicopter quadratic
            # element-by-element when windmill fails (a outside [0, 1],
            # Cn <= 0, no-convergence).
            for i_r, r in enumerate(r_stations):
                elem = None
                if v_climb < -1e-9:
                    elem = _solve_bem_element_windmill(
                        r=float(r), dr=dr,
                        chord=float(chord_per_station[i_r]),
                        twist_rad=float(twist_per_station_rad[i_r]),
                        collective_rad=inputs.collective_rad,
                        omega=omega, v_climb=v_climb, rho=inputs.rho_kg_m3,
                        n_blades=blade.n_blades, radius_m=blade.radius_m,
                        polar=self._polar, use_tip_loss=airfoil.tip_loss,
                        root_cutout_m=blade.root_cutout_m,
                    )
                if elem is None:
                    elem = solve_bem_element(
                        r=float(r), dr=dr,
                        chord=float(chord_per_station[i_r]),
                        twist_rad=float(twist_per_station_rad[i_r]),
                        collective_rad=inputs.collective_rad,
                        omega=omega, v_climb=v_climb, rho=inputs.rho_kg_m3,
                        n_blades=blade.n_blades, radius_m=blade.radius_m,
                        polar=self._polar, use_tip_loss=airfoil.tip_loss,
                        root_cutout_m=blade.root_cutout_m,
                    )
                T_total += elem.dT
                Q_total += elem.dQ

        F_world = -T_total * hub_axis_ned
        M_orbital = inputs.R_hub @ np.array([Mx_hub, My_hub, 0.0])
        M_spin = Q_total * hub_axis_ned

        I_ode = (
            self.defn.autorotation.I_ode_kgm2
            if self.defn.autorotation.I_ode_kgm2 is not None
            else 1.0
        )
        d_omega = (-Q_total + inputs.motor_torque_Nm) / I_ode
        d_spin_angle = omega

        result = AeroResult(
            F_world=F_world,
            M_orbital=M_orbital,
            Q_spin=Q_total,
            M_spin=M_spin,
        )
        derivative = QuasiStaticRotorState(
            omega_rad_s=d_omega,
            spin_angle_rad=d_spin_angle,
        )
        return result, derivative

