# Oye 2-Stage Annular Dynamic Inflow

Implementation: [`dynbem_rs/src/oye.rs`](../dynbem_rs/src/oye.rs)

Formal math (state equations, time constants, ODE):
[BEM_COMMON.md section 11](BEM_COMMON.md#11-oye-2-stage-annular-dynamic-inflow)

References: Oye (1990), Snel & Schepers (1995), OpenFAST AeroDyn Theory v3.5
section 6.3.4 (DBEMT_Mod=1).

---

## Overview

Per-annulus alternative to Pitt-Peters. Each radial annulus carries two
first-order filter states (W_int[i], W[i]) that relax the local induction toward
a Glauert momentum target. Annuli are independent -- no global coupling.

**Why this exists alongside Pitt-Peters.** Pitt-Peters couples C_T, C_M_hub, and
C_L_hub globally through the L-matrix. That BEM-driven feedback is numerically
stiff at high advance ratios and in descent + edgewise wind. Oye's per-annulus
filters have no global coupling, so they stay numerically stable in exactly those
regimes. This is the same trade-off OpenFAST's DBEMT made.

**Advantages.** No global feedback -- stable in high-advance-ratio + descent
regimes that require adaptive time steps with Pitt-Peters. Resolves an arbitrary
radial inflow distribution (not just uniform + linear). Two-stage filter matches
measured wind-turbine inflow-lag data well.

**Disadvantages.** No azimuthal harmonic states (lambda_c/lambda_s) -- cyclic
inflow feedback is absent. State size grows with the radial grid (2*N_r). Still
needs the empirical VRS override in the vortex-ring state.

---

## State Interpretation

Per radial annulus i: two first-order filter states W_int[i] and W[i].
Total state dimension: 2*N_r.

W[i] is what the blade reads in the psi-loop. W_int[i] is the intermediate
filter stage between the momentum target W_qs[i] and W[i].

Total axial flow at annulus i:

$$\lambda_\text{local}(i) = \lambda_\text{climb} + W[i]$$

This mirrors Pitt-Peters' `lambda_total = lambda_climb + lambda_0`. **Both terms
must appear.** Without `lambda_climb` the blade never sees net-upward flow in WBS
and autorotation is suppressed.

**W sign convention** (same as Pitt-Peters lambda_0):
- W > 0 in hover / helicopter (induced flow downward through disk)
- W > 0 in autorotation too: induction slows the upward freestream in the same
  +Z (downward in NED) direction it would push in helicopter mode

---

## Quasi-Steady Target

W_qs[i] is solved per annulus from Glauert momentum (linear form):

$$W_{qs,i} = \frac{dC_T/dx|_i}{4\,x_i\,F_i\,\mu_T}$$

with rotor-mean $\mu_T = \sqrt{\mu^2 + (\lambda_\text{climb} + v_{0,\text{mean}})^2} / (\Omega R)$.

**Why the linear form, not the axial-momentum form.** The pure axial form
`4*x*lambda_r*W = dCT/dx` was tried first and was unstable in forward flight: at
descent + edgewise wind, small `lambda_r` in the denominator makes W blow up. The
linear form above is what Pitt-Peters uses in its aggregate
`lambda_0_ss = T / (2*rho*A*V_T*Omega*R)` and is numerically stable. See the
comment block above `solve_w_qs` in oye.rs.

---

## Two-Stage Filter ODE

Time constants:

$$\tau_1 = \frac{1.1}{1 - 1.3\min(a,\,0.5)} \cdot \frac{R}{V_\text{mf}}, \qquad \tau_2(r) = \bigl(0.39 - 0.26\,(r/R)^2\bigr)\,\tau_1$$

The coupling constant k=0.6 (OpenFAST default) enters through a dW_qs/dt term,
which is set to zero across each outer step (DBEMT_Mod=1 approximation):

$$\dot{W}_{\text{int},i} \approx \frac{W_{qs,i} - W_{\text{int},i}}{\tau_1}, \qquad \dot{W}_i = \frac{W_{\text{int},i} - W_i}{\tau_2(r_i)}$$

tau_1 is rotor-mean (depends on a_avg, not per-annulus); tau_2(r) varies with
radius. Both time constants are well above the envelope's outer dt, so the
semi-implicit Euler in `envelope/point_mass.py` is gentle damping at most.

---

## What Oye Cannot Do

**No cyclic inflow harmonics.** There is no lambda_c/lambda_s state, so the
inflow does not develop an asymmetric tilt in response to cyclic pitching or
rolling moments. Cyclic *control* still works (hub moments respond correctly to
tilt_lon/tilt_lat), but the cyclic *inflow feedback* that reduces steady-state
moment in Pitt-Peters is absent.
`tests/test_cyclic.py::test_cyclic_inflow_reduces_hub_moment` does not apply.

**No wake-skew off-diagonal.** The `-L_off*C_T` term that produces Glauert wake
skew from thrust forcing in Pitt-Peters has no analogue. Wake skew must come from
the BEM psi-loop's asymmetric loading alone.

**No radial coupling.** Each W[i] evolves independently. Cross-annulus coupling
enters only through the rotor-mean mu_T in the time-constant formulas and V_h in
the VRS override.

---

## Cyclic Input

Cyclic pitch flows through the same `cyclic_coeffs` path as Pitt-Peters; the
psi-loop produces correct hub moments. What is missing compared to Pitt-Peters:
the cyclic-driven hub moment does not develop a counter-acting inflow harmonic
(no lambda_c/lambda_s states), so the steady-state moment is over-predicted vs
Pitt-Peters at hover. Cyclic *control* (sign and order-of-magnitude) is right;
cyclic *inflow damping* is absent.
