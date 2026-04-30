"""Step 2: hover delta-CQ validation against Table I hover rows.
DeltaCQ = CQ(theta) - CQ(theta=0), compared to paper's measured DeltaCQ.
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

def cq(theta_deg, rpm):
    omega = rpm * math.pi / 30.0
    A = math.pi * R**2
    inp = RotorInputs(collective_rad=math.radians(theta_deg),
                      tilt_lon=0.0, tilt_lat=0.0, R_hub=np.eye(3),
                      v_hub_world=np.zeros(3), wind_world=np.zeros(3), t=0.0)
    res, _ = model.compute_forces(inp, QuasiStaticRotorState(omega_rad_s=omega))
    return res.Q_spin / (RHO * A * (omega*R)**2 * R)

cq0_1200 = cq(0.0, 1200)
cq0_1600 = cq(0.0, 1600)

# (run, CT_nom, rpm, theta, delta_CQ_measured)
# Runs 6 and 38 hover have no measured CQ — skipped
hover_rows = [
    ( 3, 0.004, 1200,  8.83, 0.000282),
    ( 4, 0.004, 1200,  8.66, 0.000285),
    ( 5, 0.004, 1200,  8.86, 0.000292),
    ( 9, 0.002, 1600,  5.13, 0.000062),
    (14, 0.005, 1600, 10.06, 0.000414),
    (32, 0.002, 1200,  5.32, 0.000087),
    (34, 0.002, 1600,  5.33, 0.000087),
    (35, 0.004, 1600,  8.68, 0.000271),
    (36, 0.004, 1200,  8.49, 0.000242),
]

print(f"{'Run':>4}  {'CT':>5}  {'RPM':>4}  {'theta':>6}  {'dCQ_meas':>10}  {'dCQ_pred':>10}  {'err%':>7}")
print("-" * 66)
errs = []
for run, CT_nom, rpm, theta, dcq_meas in hover_rows:
    cq0 = cq0_1200 if rpm == 1200 else cq0_1600
    dcq_pred = cq(theta, rpm) - cq0
    e = (dcq_pred - dcq_meas) / dcq_meas * 100
    errs.append(e)
    print(f"{run:>4}  {CT_nom:.3f}  {rpm:>4}  {theta:>6.2f}  {dcq_meas:>10.6f}  {dcq_pred:>10.6f}  {e:>+7.1f}%")
print(f"\nMean err: {sum(errs)/len(errs):+.1f}%  RMSE: {math.sqrt(sum(e**2 for e in errs)/len(errs)):.1f}%")
