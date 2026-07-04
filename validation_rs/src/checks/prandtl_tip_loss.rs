use crate::helpers::*;
use crate::report::Report;

pub fn check_prandtl_tip_loss(r: &mut Report) {
    r.begin_module(
        "prandtl_tip_loss",
        "Tip-loss flag must reduce global loads (directional)",
    );
    let mut defn_on = theory_rotor(12, 0.0);
    defn_on.blade.tip_loss = true;
    let mut defn_off = theory_rotor(12, 0.0);
    defn_off.blade.tip_loss = false;
    let rotor_on = make_rotor(&defn_on);
    let rotor_off = make_rotor(&defn_off);
    let fc = hover_fc(9.0);
    let (res_on, _s_on) = settle(&rotor_on, &fc, 8);
    let (res_off, _s_off) = settle(&rotor_off, &fc, 8);
    let ct_on = c_t(res_on.thrust);
    let ct_off = c_t(res_off.thrust);
    let cq_on = c_q(res_on.torque);
    let cq_off = c_q(res_off.torque);
    r.assert_bool(
        "hover_col=9",
        "CT_tip_loss_not_higher",
        ct_on,
        ct_off,
        ct_on <= ct_off + 1e-6,
        "tip-loss must not increase CT",
    );
    r.assert_bool(
        "hover_col=9",
        "CQ_tip_loss_not_higher",
        cq_on,
        cq_off,
        cq_on <= cq_off + 1e-6,
        "tip-loss must not increase CQ",
    );
    r.info("hover_col=9", "CT_on", ct_on, ct_off);
    r.info("hover_col=9", "CQ_on", cq_on, cq_off);
}
