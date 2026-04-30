"""Step 4: hold theta fixed at hover value, vary descent rate.
Shows how BEM predicts CT changes with descent at constant collective.
Paper expectation (VRS): CT barely changes (author barely adjusts theta).
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
    A_ = math.pi * R**2
    omega_ = rpm * math.pi / 30.0
    CT = T / (RHO * A_ * (omega_*R)**2)
    CQ = res.Q_spin / (RHO * A_ * (omega_*R)**2 * R)
    return CT, CQ

# Use Run 9: CT=0.002, 1600rpm, theta_hover=5.13
# Compare paper theta (adjusts to maintain CT) vs fixed theta=5.13
rpm, theta_hover, CT_nom = 1600, 5.13, 0.002
omega = rpm * math.pi / 30.0

rows = [
    (0.0000, 5.13, 0.00),
    (0.0136, 5.26, 0.43),
    (0.0303, 4.88, 0.98),
    (0.0463, 5.01, 1.50),
    (0.0650, -0.11, 2.10),
    (0.0720, -1.32, 2.32),
    (0.0763, -1.66, 2.46),
]

print("Run 9  CT=0.002  1600rpm  -- CT at fixed hover theta=5.13 vs paper theta")
print(f"{'lam2':>6}  {'paper_theta':>12}  {'paper_CT_nom':>13}  {'CT_fixed_theta':>15}  {'CT_ratio':>10}")
print("-"*72)
for (v_or, theta_p, lam2) in rows:
    v_climb = -(v_or * omega * R)
    ct_fixed, _ = run_bem(theta_hover, rpm, v_climb)
    ratio = ct_fixed / CT_nom
    print(f"{lam2:>6.2f}  {theta_p:>12.2f}  {CT_nom:>13.5f}  {ct_fixed:>15.5f}  {ratio:>10.2f}x")

print()
print("Interpretation:")
print("  paper_theta stays ~5° in VRS => real rotor CT barely changes with descent")
print("  CT_fixed_theta >> CT_nom in BEM => BEM greatly over-predicts descent thrust")
print("  Root cause: BEM has no VRS model; momentum theory breaks down for 0 < lam2 < ~2")
