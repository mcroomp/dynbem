# Pitt-Peters 3-State Dynamic Inflow

Implementation: [`dynbem_rs/src/pitt_peters.rs`](../dynbem_rs/src/pitt_peters.rs)

Formal math (state equations, L-matrix coefficients, ODE):
[BEM_COMMON.md section 10](BEM_COMMON.md#10-pitt-peters-3-state-dynamic-inflow)

Canonical reference for the original formulation:
[`Research/Peters_Nikolsky_2008/CLAUDE.md`](../Research/Peters_Nikolsky_2008/CLAUDE.md)
(L-matrix, M matrix, V mass-flow, forcing sign convention from David Peters,
"How Dynamic Inflow Survives in the Competitive World of Rotorcraft Aerodynamics",
JAHS 54(1):011001 (2009)). **Read that file before touching any signs or
coefficients.**

---

## Overview

Three global inflow harmonics -- uniform, longitudinal tilt, lateral tilt --
relax toward a momentum-theory steady state through Peters' apparent-mass time
constants. Thrust, rolling moment, and pitching moment couple to the three states
through the L-matrix.

**Advantages.** Captures dynamic-inflow lag and the cyclic inflow/hub-moment
feedback that quasi-static BEM misses. The L_off cross-term reproduces Glauert
wake skew naturally. Only three states -- cheap to integrate, standard for
rotor flight dynamics and trim.

**Disadvantages.** The global L-matrix coupling is numerically stiff at high
advance ratio and in descent + edgewise wind, demanding small or adaptive time
steps. Radial inflow shape is fixed (uniform + linear). Momentum theory breaks
down in VRS, requiring the empirical override.

---

## Hub-Frame Moment Conventions

In the psi-loop each blade element contributes thrust dT along -z_hub. With
r = r*r_hat(psi):

    dM_hub = r_pos x F = r * dT * [sin(psi), cos(psi), 0]

so Mx_hub = sum r*dT*sin(psi) (roll), My_hub = sum r*dT*cos(psi) (pitch).

Non-dimensional form:

$$C_T = \frac{T}{\rho A (\Omega R)^2}, \quad C_L = \frac{M_{x,\text{hub}}}{\rho A (\Omega R)^2 R}, \quad C_M = \frac{M_{y,\text{hub}}}{\rho A (\Omega R)^2 R}$$

$M_x > 0$ is roll-right, $M_y > 0$ is pitch-up (NED right-hand rule).

**BladeAD sign differences** (BladeAD uses psi=0 at the tail, dMy = -r*cos(psi)*dT):

- our lambda_c = -BladeAD lambda_c
- our lambda_s = -BladeAD lambda_s
- our C_M_hub = +BladeAD C_My
- our C_L_hub = -BladeAD C_Mx

---

## State Interpretation

`lambda_0` (and `lambda_c`, `lambda_s`) is the **induced** inflow ratio
`v_i / (Omega*R)`, not the total inflow. The total axial flow seen by each blade
element is:

$$\lambda_\text{total} = \lambda_0 + \lambda_\text{climb}$$

where `lambda_climb = v_climb / (Omega*R) < 0` in descent. **Both terms must
appear inside the BEM sweep.** Without the freestream term the blade never sees
net-upward flow in WBS, so CQ never goes negative and autorotation is suppressed.

---

## Peters' L-Matrix -- Sign Translation

Peters' own convention (Eq 10 of the Nikolsky lecture, psi=0 at the tail,
state ordering (nu_0, nu_s, nu_c), X = tan(chi/2)):

```text
[L] = | 1/2           0           -15*pi*X/64 |
      | 0             2*(1+X^2)    0           |
      | 15*pi*X/64    0            2*(1-X^2)   |
```

Forcing vector: `{C_T, -C_L, -C_M}`.

Our convention: psi=0 at +X (nose). Translation: our lambda_c = -nu_c,
our lambda_s = -nu_s. After translation, with

    L_off = (15*pi/64) * tan(chi/2)
    L_cc  = 4*cos(chi) / (1 + cos(chi))   # = 2*(1-X^2)
    L_ss  = 4 / (1 + cos(chi))            # = 2*(1+X^2)

the steady-state targets are:

$$\lambda_{0,ss} = \frac{C_T}{2\mu_T} + \frac{L_\text{off}\,C_M}{\mu_T}$$

$$\lambda_{c,ss} = \frac{-L_\text{off}\,C_T + L_{cc}\,C_M}{\mu_T}$$

$$\lambda_{s,ss} = \frac{L_{ss}\,C_L}{\mu_T}$$

where $\mu_T = \sqrt{\mu^2 + \lambda_\text{total}^2}$ (see below).

The `-L_off*C_T` term in `lambda_c_ss` is the wake-skew cross-coupling -- it
produces the Glauert inflow tilt naturally from thrust forcing. The closed-form
Glauert tilt has been removed; **do not re-add it** (would double-count).

VRS region (`v_climb < 0`, `0 < V_descent/V_h < 2`) overrides `lambda_0_ss`
with the Leishman empirical polynomial. The cross-coupling is also skipped in
that regime (momentum theory does not apply in a recirculating wake).

---

## ODE

Peters' apparent mass M = diag(8/(3*pi), 16/(45*pi), 16/(45*pi)) gives time
constants:

$$\tau_0 = \frac{8R}{3\pi V_T}, \qquad \tau_c = \tau_s = \frac{16R}{45\pi V_T}$$

State derivatives returned each call:

$$\dot\lambda_0 = \frac{\lambda_{0,ss} - \lambda_0}{\tau_0}, \quad \dot\lambda_c = \frac{\lambda_{c,ss} - \lambda_c}{\tau_c}, \quad \dot\lambda_s = \frac{\lambda_{s,ss} - \lambda_s}{\tau_s}$$

**V_T floor:** `1e-2 * max(Omega*R, 1)` prevents `tau_0 -> inf` in the middle of
VRS where upward freestream approximately cancels downward induced flow. The floor
value does not affect stability; `tau_0 -> large` in VRS is physically correct.

---

## Mass-Flow Parameter: mu_T vs Peters' V

We use `mu_T = sqrt(mu^2 + lambda_total^2)` (classical Glauert). Peters uses a
different `V` (his Eq 8) that agrees at high speed but differs by ~2x in hover:

- our mu_T: reproduces classical Glauert hover `lambda_0 = sqrt(C_T/2)` (correct)
- Peters' V: gives `sqrt(C_T)/2` (factor sqrt(2) different -- likely a C_T
  normalization convention in his paper)

The L-matrix structure matches Peters exactly; only the scalar scaling differs.
Swapping to Peters' V requires re-validation against hover data -- defer until
needed.

---

## Wind-Axis Rotation

The L-matrix is diagonal-plus-off-diagonal in **wind axes**. The code rotates
into wind axes before applying steady-state targets and rotates derivatives back,
making the model rotationally covariant for oblique flight (mu_y != 0).

Let `beta = atan2(v_y_hub, -v_x_hub)` (beta_wind; beta=0 means pure longitudinal
wind). Then:

- Aerodynamic forcing is rotated to wind axes:
  `C_M_wind = cos(beta)*C_M + sin(beta)*C_L`
  `C_L_wind = -sin(beta)*C_M + cos(beta)*C_L`
- Current (lambda_c, lambda_s) states are rotated the same way
- Relaxation `(ss - current)/tau` is evaluated entirely in wind axes
- Derivatives are rotated back to hub-frame coordinates:
  `d_lam_c = cos(beta)*d_lam_c_wind - sin(beta)*d_lam_s_wind`
  `d_lam_s = sin(beta)*d_lam_c_wind + cos(beta)*d_lam_s_wind`

**Stability history.** Wind-axis rotation was reverted once because it destabilised
the tethered-rotor envelope at descent + edgewise operating points via the nonlinear
feedback `lambda_c -> BEM(lam_local) -> C_L_hub -> lambda_s_ss`. The fix is the
semi-implicit Euler damping in `envelope/point_mass.py::_step_state_semi_implicit`
-- each state is damped by `(1 + dt/tau)^-1` using per-state time constants from
`inflow_taus()`. **Do not remove the semi-implicit damping and expect wind-axis
rotation to stay stable.**

Covered by `tests/test_pitt_peters.py::TestPittPetersObliqueFlow`.

---

## VRS Notes

**Polynomial sign.** The Leishman (2000) polynomial is descent-positive
(lambda_2 = V_descent/V_h):

$$\lambda_1/V_h = 1 + 1.125\lambda_2 - 1.372\lambda_2^2 + 1.718\lambda_2^3 - 0.655\lambda_2^4$$

This is NOT the form with coefficients (-1.125, -1.372, -1.718, -0.655), which
applies when the argument is V_climb/V_h (negative for descent). Both forms are
equivalent; this code uses descent-positive throughout.

**Why CT still rises in deep VRS.** At lambda_2 ~ 1.5 the polynomial gives
`lambda_0_ss ~ 2*V_h/(Omega*R)`. With `lambda_climb ~ -1.5*V_h/(Omega*R)`, net
blade inflow `lambda_total ~ 0.5*V_h/(Omega*R)` -- less than hover, so AoA
increases and CT rises. The real VRS has recirculating wakes that further restrict
throughflow; the 1-D polynomial captures mean induced velocity but not 3-D
blockage. This is a known limitation of all momentum-based VRS models.
