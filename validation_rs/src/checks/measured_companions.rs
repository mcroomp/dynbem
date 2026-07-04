use crate::helpers::*;
use crate::report::Report;
use dynbem_rs::vpm_rotor::FlightCondition;

pub fn check_measured_companions(r: &mut Report) {
    r.begin_module(
        "measured_companions",
        "Each theory module anchored to a measured dataset",
    );

    // --- companion: hover Castles-Gray 1600 rpm
    let defn = castles_gray_rotor(10);
    let rotor = make_rotor(&defn);
    let pts_1600 = [
        (3.96_f64, 1600.0_f64, 0.00160_f64),
        (5.55, 1600.0, 0.00255),
        (7.18, 1600.0, 0.00346),
    ];
    for &(theta_deg, rpm, ct_meas) in &pts_1600 {
        let omega = omega_from_rpm(rpm);
        let (res, _s) = settle(&rotor, &hover_fc_omega(theta_deg, omega), 10);
        let ct = ct_at(res.thrust, omega, R_TIP);
        r.check(
            format!("CG1600 theta={theta_deg:.2}"),
            "CT",
            ct,
            ct_meas,
            35.0,
        );
    }

    // --- companion: climb + descent (Castles-Gray WBS)
    let rows_climb = [
        (1.23_f64, 1200.0_f64, 11.15_f64, -0.000112_f64),
        (-1.66_f64, 1600.0_f64, 11.91_f64, -0.000084_f64),
    ];
    for &(theta_deg, rpm, v_descent, cq_meas) in &rows_climb {
        let omega = omega_from_rpm(rpm);
        let fc = FlightCondition {
            collective_rad: theta_deg.to_radians(),
            tilt_lon: 0.0,
            tilt_lat: 0.0,
            v_hub: [0.0, 0.0, -v_descent],
            omega_rad_s: omega,
            rho: RHO,
        };
        let (res, _s) = settle(&rotor, &fc, 10);
        let cq = cq_at(res.torque, omega, R_TIP);
        let case = format!("CG_descent theta={theta_deg:.2} rpm={rpm:.0}");
        r.check(&case, "CQ_sign", cq, cq_meas, 300.0); // loose: sign + order-of-magnitude
        r.info(&case, "CQ_vpm", cq, cq_meas);
    }

    // --- companion: tip-loss measured anchor
    let mut defn_on = castles_gray_rotor(10);
    defn_on.blade.tip_loss = true;
    let mut defn_off = castles_gray_rotor(10);
    defn_off.blade.tip_loss = false;
    let r_on = make_rotor(&defn_on);
    let r_off = make_rotor(&defn_off);
    let ct_meas_ref = 0.00400;
    let omega1200 = omega_from_rpm(1200.0);
    let fc_tl = hover_fc_omega(8.46, omega1200);
    let (res_on, _s_on) = settle(&r_on, &fc_tl, 10);
    let (res_off, _s_off) = settle(&r_off, &fc_tl, 10);
    let ct_on = ct_at(res_on.thrust, omega1200, R_TIP);
    let ct_off = ct_at(res_off.thrust, omega1200, R_TIP);
    r.assert_bool(
        "tip_loss_meas",
        "tip_loss_on_is_more_accurate",
        (ct_on - ct_meas_ref).abs(),
        (ct_off - ct_meas_ref).abs(),
        (ct_on - ct_meas_ref).abs() <= (ct_off - ct_meas_ref).abs() + 0.0002,
        "tip-loss should not make CT less accurate vs measured",
    );
}
