"""Step 1: hover CT validation against Table I hover rows (V/OR = 0).
Compare BEM CT to paper's nominal CT at the given collective.
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

def ct(theta_deg, rpm):
    omega = rpm * math.pi / 30.0
    A = math.pi * R**2
    inp = RotorInputs(collective_rad=math.radians(theta_deg),
                      tilt_lon=0.0, tilt_lat=0.0, R_hub=np.eye(3),
                      v_hub_world=np.zeros(3), wind_world=np.zeros(3), t=0.0)
    res, _ = model.compute_forces(inp, QuasiStaticRotorState(omega_rad_s=omega))
    return -res.F_world[2] / (RHO * A * (omega*R)**2)

# Table I hover rows: (run, CT_nominal, rpm, theta_0.75R)
hover_rows = [
    ( 3, 0.004, 1200,  8.83),
    ( 4, 0.004, 1200,  8.66),
    ( 5, 0.004, 1200,  8.86),
    ( 6, 0.002, 1200,  5.33),
    ( 9, 0.002, 1600,  5.13),
    (14, 0.005, 1600, 10.06),
    (32, 0.002, 1200,  5.32),
    (34, 0.002, 1600,  5.33),
    (35, 0.004, 1600,  8.68),
    (36, 0.004, 1200,  8.49),
    (38, 0.005, 1600, 10.19),
]

print(f"{'Run':>4}  {'CT_nom':>8}  {'RPM':>4}  {'theta':>6}  {'CT_pred':>8}  {'err%':>7}")
print("-" * 52)
errs = []
for run, CT_nom, rpm, theta in hover_rows:
    CT_p = ct(theta, rpm)
    e = (CT_p - CT_nom) / CT_nom * 100
    errs.append(e)
    print(f"{run:>4}  {CT_nom:>8.5f}  {rpm:>4}  {theta:>6.2f}  {CT_p:>8.5f}  {e:>+7.1f}%")
print(f"\nMean err: {sum(errs)/len(errs):+.1f}%  RMSE: {math.sqrt(sum(e**2 for e in errs)/len(errs)):.1f}%")
