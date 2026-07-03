use crate::helpers::*;
use crate::report::Report;

pub fn check_hover_castles_gray(r: &mut Report) {
    r.begin_module("hover_castles_gray", "Hover vs measured Castles-Gray NACA TN-2474 Table V");
    struct M {
        theta: f64,
        rpm: f64,
        ct: f64,
        cq: f64,
    }
    let meas = [
        M { theta: 8.46, rpm: 1200.0, ct: 0.00400, cq: 0.000226 },
        M { theta: 10.29, rpm: 1200.0, ct: 0.00488, cq: 0.000342 },
    ];
    let defn = castles_gray_rotor(10);
    let rotor = make_rotor(&defn);
    for m in &meas {
        let omega = omega_from_rpm(m.rpm);
        let fc = hover_fc_omega(m.theta, omega);
        let (res, _s) = settle(&rotor, &fc, 10);
        let ct = ct_at(res.thrust, omega, R_TIP);
        let cq = cq_at(res.torque, omega, R_TIP);
        let fm = ct.powf(1.5) / (2f64.sqrt() * cq.max(1e-9));
        let fm_meas = m.ct.powf(1.5) / (2f64.sqrt() * m.cq);
        let case = format!("theta={:.2} rpm={:.0}", m.theta, m.rpm);
        r.check(&case, "CT", ct, m.ct, 25.0);
        r.info(&case, "CQ", cq, m.cq);
        r.info(&case, "FM", fm, fm_meas);
    }
}
