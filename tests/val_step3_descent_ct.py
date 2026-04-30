"""Step 3: descent CT — does BEM find the right collective to maintain CT_nominal?

For each Table I row with V > 0, we:
  (a) Ask: what does BEM predict CT if we hold theta FIXED at the hover value?
  (b) Ask: what theta does BEM need to match CT_nominal at this descent rate?
      (bisection inversion)
  (c) Compare that required theta_BEM to the paper's theta.

This tells us whether the BEM's theta-vs-descent curve matches the paper,
independent of absolute CT scale errors.
"""
import math, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent.parent))
from aero.bem import BEMModel
from aero import RotorInputs
import aero.rotor_definition as rotor_definition
from aero.rotor_state import QuasiStaticRotorState

defn = rotor_definition.load(str(Path(__file__).parent.parent / "rotors/castles_gray_6ft/rotor.yaml"))
model = BEMModel(defn=defn)
R = defn.blade.radius_m
RHO = 1.225

def run_bem(theta_deg, rpm, v_climb_ms):
    omega = rpm * math.pi / 30.0
    A = math.pi * R**2
    inp = RotorInputs(collective_rad=math.radians(theta_deg),
                      tilt_lon=0.0, tilt_lat=0.0, R_hub=np.eye(3),
                      v_hub_world=np.zeros(3),
                      wind_world=np.array([0.0, 0.0, v_climb_ms]), t=0.0)
    res, _ = model.compute_forces(inp, QuasiStaticRotorState(omega_rad_s=omega))
    T = -res.F_world[2]
    CT = T / (RHO * A * (omega*R)**2)
    CQ = res.Q_spin / (RHO * A * (omega*R)**2 * R)
    return CT, CQ

def find_theta_for_ct(CT_target, rpm, v_climb_ms, lo=-5.0, hi=20.0, tol=1e-4):
    """Bisect on theta to find CT = CT_target."""
    for _ in range(60):
        mid = 0.5*(lo + hi)
        ct, _ = run_bem(mid, rpm, v_climb_ms)
        if ct > CT_target:
            hi = mid
        else:
            lo = mid
        if hi - lo < tol:
            break
    return 0.5*(lo + hi)

# Use Run 5 (CT=0.004, 1200rpm) as representative — HIGH confidence, dense data
# theta_hover = 8.86
rows_run5 = [
    # (V/OR, theta_paper, lam2)
    (0.0000, 8.86, 0.00),
    (0.0176, 8.70, 0.40),
    (0.0388, 8.40, 0.89),
    (0.0552, 8.53, 1.26),
    (0.0670, 8.46, 1.53),
    (0.0734, 4.60, 1.67),
    (0.0790, 3.99, 1.81),
    (0.0856, 2.58, 1.96),
    (0.0949, 1.18, 2.17),
]
CT_nom, rpm = 0.004, 1200

print("Run 5  CT=0.004  1200rpm")
print(f"{'lam2':>6}  {'V/OR':>6}  {'theta_paper':>12}  {'theta_BEM':>10}  {'dtheta':>8}  {'CT@theta_p':>12}")
print("-"*68)
omega = rpm * math.pi / 30.0
for (v_or, theta_p, lam2) in rows_run5:
    v_climb = -(v_or * omega * R)
    theta_b = find_theta_for_ct(CT_nom, rpm, v_climb)
    ct_at_paper_theta, _ = run_bem(theta_p, rpm, v_climb)
    dtheta = theta_b - theta_p
    print(f"{lam2:>6.2f}  {v_or:>6.4f}  {theta_p:>12.2f}  {theta_b:>10.2f}  {dtheta:>+8.2f}  {ct_at_paper_theta:>12.5f}")

print()
# Now do the same for Run 9 (CT=0.002, 1600rpm) including WBS rows
rows_run9 = [
    (0.0000, 5.13, 0.00),
    (0.0136, 5.26, 0.43),
    (0.0303, 4.88, 0.98),
    (0.0463, 5.01, 1.50),
    (0.0650, -0.11, 2.10),
    (0.0720, -1.32, 2.32),
    (0.0763, -1.66, 2.46),
]
CT_nom9, rpm9 = 0.002, 1600
omega9 = rpm9 * math.pi / 30.0

print("Run 9  CT=0.002  1600rpm")
print(f"{'lam2':>6}  {'V/OR':>6}  {'theta_paper':>12}  {'theta_BEM':>10}  {'dtheta':>8}  {'CT@theta_p':>12}")
print("-"*68)
for (v_or, theta_p, lam2) in rows_run9:
    v_climb = -(v_or * omega9 * R)
    theta_b = find_theta_for_ct(CT_nom9, rpm9, v_climb)
    ct_at_paper_theta, _ = run_bem(theta_p, rpm9, v_climb)
    dtheta = theta_b - theta_p
    print(f"{lam2:>6.2f}  {v_or:>6.4f}  {theta_p:>12.2f}  {theta_b:>10.2f}  {dtheta:>+8.2f}  {ct_at_paper_theta:>12.5f}")
