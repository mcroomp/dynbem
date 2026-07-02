// Hover -- VPM vs measured data (Castles & Gray, NACA TN-2474, Table V).
//
// This is the same measured 6-ft rotor the BEM hover tests use
// (hover_empirical.rs / hover_cq_empirical.rs). Grounding the VPM hover check
// in MEASUREMENT rather than closed-form momentum theory also adjudicates
// whether a disagreement is a VPM error or a theory limitation.
//
// Measured, 1200 rpm:
//   theta (deg) | C_T     | C_Q      | FM = C_T^1.5 / (sqrt(2) C_Q)
//   8.46        | 0.00400 | 0.000226 | 0.79
//   10.29       | 0.00488 | 0.000342 | 0.71
//
// C_T is the trustworthy hover quantity. C_Q (hence FM) is the known
// under-resolved quantity for a truncated free wake in the recirculating
// hover regime (VPM_DESIGN sec 8.2); it is reported and loosely bounded, not
// asserted tightly.

use crate::common::*;

struct Meas {
    theta: f64,
    rpm: f64,
    ct: f64,
    cq: f64,
}

const CG: &[Meas] = &[
    Meas {
        theta: 8.46,
        rpm: 1200.0,
        ct: 0.00400,
        cq: 0.000226,
    },
    Meas {
        theta: 10.29,
        rpm: 1200.0,
        ct: 0.00488,
        cq: 0.000342,
    },
];

#[test]
fn hover_thrust_vs_castles_gray() {
    let defn = castles_gray_rotor(10);
    let r = defn.blade.radius_m;
    let rotor = make_rotor(&defn);

    for m in CG {
        let omega = omega_from_rpm(m.rpm);
        let fc = hover_fc_omega(m.theta, omega);
        let (res, _s) = settle(&rotor, &fc, 10);

        let ct = ct_at(res.thrust, omega, r);
        let cq = cq_at(res.torque, omega, r);
        let fm = ct.powf(1.5) / (2f64.sqrt() * cq);
        let fm_meas = m.ct.powf(1.5) / (2f64.sqrt() * m.cq);
        let ct_err = (ct - m.ct).abs() / m.ct;
        let cq_err = (cq - m.cq).abs() / m.cq;

        eprintln!(
            "CG theta={:.2} rpm={:.0}: C_T={:.5} (meas {:.5}, {:.0}%)  \
             C_Q={:.6} (meas {:.6}, {:.0}%)  FM={:.2} (meas {:.2})",
            m.theta,
            m.rpm,
            ct,
            m.ct,
            ct_err * 100.0,
            cq,
            m.cq,
            cq_err * 100.0,
            fm,
            fm_meas
        );

        // Thrust is the hover quantity the VPM should get right.
        assert!(
            ct_err < 0.25,
            "C_T {ct:.5} vs measured {:.5} ({:.0}% off, expected < 25%)",
            m.ct,
            ct_err * 100.0
        );
        // Torque / power: the VPM under-resolves hover induced power somewhat,
        // but on the real rotor it stays within ~25% of measurement.
        assert!(
            cq_err < 0.30,
            "C_Q {cq:.6} vs measured {:.6} ({:.0}% off, expected < 30%)",
            m.cq,
            cq_err * 100.0
        );
        // Figure of merit must be physical (< 1) and in a sane hover band.
        assert!(
            fm > 0.5 && fm < 1.0,
            "figure of merit {fm:.2} outside the physical hover band [0.5, 1.0]"
        );
    }
}
