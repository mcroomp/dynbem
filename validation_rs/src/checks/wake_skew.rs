use crate::helpers::*;
use crate::report::Report;
use dynbem_rs::vpm::FlightCondition;

pub fn check_wake_skew(r: &mut Report) {
    r.begin_module(
        "wake_skew",
        "Wake skew grows with mu; covariant under X/Y rotation",
    );
    let defn = theory_rotor(12, 0.0);
    let rotor = make_rotor(&defn);
    let mut chis = Vec::new();
    for &mu in &[0.08_f64, 0.16, 0.24] {
        let fc = forward_fc(8.5, mu);
        let (res, state) = settle(&rotor, &fc, 8);
        let chi = wake_skew_angle(&state, &res);
        chis.push(chi);
        r.info(format!("mu={mu:.2}"), "chi_deg", chi.to_degrees(), f64::NAN);
    }
    // Monotone increase (loose 2-deg slack each step).
    r.assert_bool(
        "mu_sweep",
        "chi_increases_01",
        chis[1],
        chis[0],
        chis[1] > chis[0] - 2f64.to_radians(),
        "chi not growing with mu",
    );
    r.assert_bool(
        "mu_sweep",
        "chi_increases_12",
        chis[2],
        chis[1],
        chis[2] > chis[1] - 2f64.to_radians(),
        "chi not growing with mu",
    );
    // Covariance: X vs Y flight direction at mu=0.20.
    let mu = 0.20;
    let v = mu * tip_speed();
    let fc_x = FlightCondition {
        collective_rad: 8.5f64.to_radians(),
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        v_hub: [v, 0.0, 0.0],
        omega_rad_s: OMEGA,
        rho: RHO,
    };
    let fc_y = FlightCondition {
        collective_rad: 8.5f64.to_radians(),
        tilt_lon: 0.0,
        tilt_lat: 0.0,
        v_hub: [0.0, v, 0.0],
        omega_rad_s: OMEGA,
        rho: RHO,
    };
    let (res_x, st_x) = settle(&rotor, &fc_x, 8);
    let (res_y, st_y) = settle(&rotor, &fc_y, 8);
    let chi_x = wake_skew_angle(&st_x, &res_x);
    let chi_y = wake_skew_angle(&st_y, &res_y);
    let dchi = (chi_x - chi_y).abs().to_degrees();
    let hx = (res_x.wake_centroid[0].powi(2) + res_x.wake_centroid[1].powi(2)).sqrt();
    let hy = (res_y.wake_centroid[0].powi(2) + res_y.wake_centroid[1].powi(2)).sqrt();
    let dh = (hx - hy).abs() / hx.max(1e-9);
    r.check("covariance_mu=0.20", "chi_xy_diff_deg", dchi, 0.0, 8.0);
    r.check(
        "covariance_mu=0.20",
        "horiz_shift_diff_pct",
        dh * 100.0,
        0.0,
        30.0,
    );
}
