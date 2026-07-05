//! Validation check modules, one per file for clarity and organization.

mod autorotation;
mod blade_element_hover;
mod climb_momentum;
mod cyclic_phase_servo;
mod cyclic_sign;
mod descent_empirical;
mod flap_directional;
mod flapping_harmonics;
mod glauert_forward_inflow;
mod hover_castles_gray;
mod hover_cq_empirical;
mod hover_empirical;
mod measured_companions;
mod prandtl_tip_loss;
mod servo_flap;
mod vpm_aging_survey;
mod vpm_empirical;
mod wake_skew;

pub use autorotation::check_autorotation;
pub use blade_element_hover::check_blade_element_hover;
pub use climb_momentum::check_climb_momentum;
pub use cyclic_phase_servo::check_cyclic_phase_servo;
pub use cyclic_sign::check_cyclic_sign;
pub use descent_empirical::check_descent_empirical;
pub use flap_directional::check_flap_directional;
pub use flapping_harmonics::check_flapping_harmonics;
pub use glauert_forward_inflow::check_glauert_forward_inflow;
pub use hover_castles_gray::check_hover_castles_gray;
pub use hover_cq_empirical::check_hover_cq_empirical;
pub use hover_empirical::check_hover_empirical;
pub use measured_companions::check_measured_companions;
pub use prandtl_tip_loss::check_prandtl_tip_loss;
pub use servo_flap::check_servo_flap;
pub use vpm_aging_survey::check_vpm_aging_survey;
pub use vpm_empirical::check_vpm_empirical;
pub use wake_skew::check_wake_skew;

use crate::report::Report;

/// Run all 17 validation checks in order.
pub fn run_all_checks(report: &mut Report) {
    run_filtered_checks(report, None);
}

/// Run only the checks whose module name contains `filter` (case-sensitive).
/// Pass `None` to run all checks.
pub fn run_filtered_checks(report: &mut Report, filter: Option<&str>) {
    macro_rules! maybe {
        ($name:expr, $fn:expr) => {
            if filter.map_or(true, |f| $name.contains(f)) {
                $fn(report);
            }
        };
    }
    maybe!("blade_element_hover", check_blade_element_hover);
    maybe!("hover_castles_gray", check_hover_castles_gray);
    maybe!("climb_momentum", check_climb_momentum);
    maybe!("prandtl_tip_loss", check_prandtl_tip_loss);
    maybe!("glauert_forward_inflow", check_glauert_forward_inflow);
    maybe!("wake_skew", check_wake_skew);
    maybe!("flapping_harmonics", check_flapping_harmonics);
    maybe!("autorotation", check_autorotation);
    maybe!("measured_companions", check_measured_companions);
    maybe!("cyclic_sign", check_cyclic_sign);
    maybe!("servo_flap", check_servo_flap);
    maybe!("cyclic_phase_servo", check_cyclic_phase_servo);
    maybe!("flap_directional", check_flap_directional);
    maybe!("hover_empirical", check_hover_empirical);
    maybe!("hover_cq_empirical", check_hover_cq_empirical);
    maybe!("descent_empirical", check_descent_empirical);
    maybe!("vpm_empirical", check_vpm_empirical);

    // Long opt-in survey: runs the VPM aging config across ALL empirical
    // datasets. Excluded from run_all (filter=None) because it is slow; run it
    // explicitly with `-- vpm_aging_survey`.
    if filter == Some("vpm_aging_survey") {
        check_vpm_aging_survey(report);
    }
}
