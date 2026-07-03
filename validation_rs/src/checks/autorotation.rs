use crate::helpers::*;
use crate::report::Report;
use dynbem_rs::vpm_rotor::FlightCondition;

pub fn check_autorotation(r: &mut Report) {
    r.begin_module(
        "autorotation",
        "Directional: VPM reaches negative-torque branch in descent+edgewise",
    );
    let defn = theory_rotor(12, 0.0);
    let rotor = make_rotor(&defn);
    let candidates = [
        (2.0_f64, 0.18_f64, -2.0_f64),
        (2.0_f64, 0.22_f64, -3.0_f64),
        (1.5_f64, 0.25_f64, -3.0_f64),
        (1.0_f64, 0.28_f64, -4.0_f64),
    ];
    let mut found = false;
    let mut best_cq = f64::INFINITY;
    for &(col, mu, vz) in &candidates {
        let v = mu * tip_speed();
        let fc = FlightCondition {
            collective_rad: col.to_radians(),
            tilt_lon: 0.0,
            tilt_lat: 0.0,
            v_hub: [v, 0.0, vz],
            omega_rad_s: OMEGA,
            rho: RHO,
        };
        let (res, _s) = settle(&rotor, &fc, 8);
        let cq = c_q(res.torque);
        let ct = c_t(res.thrust);
        best_cq = best_cq.min(cq);
        let case = format!("col={col:.1} mu={mu:.2} vz={vz:.1}");
        r.info(&case, "CQ", cq, f64::NAN);
        r.info(&case, "CT", ct, f64::NAN);
        if cq < -5e-6 {
            found = true;
            break;
        }
    }
    r.assert_bool(
        "candidates",
        "found_negative_torque",
        best_cq,
        0.0,
        found,
        &format!("no autorotation point found (best CQ={best_cq:+.5})"),
    );
}
