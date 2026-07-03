//! Validation check modules, one per file for clarity and organization.

mod blade_element_hover;
mod hover_castles_gray;
mod climb_momentum;
mod prandtl_tip_loss;
mod glauert_forward_inflow;
mod wake_skew;
mod flapping_harmonics;
mod autorotation;
mod measured_companions;
mod cyclic_sign;
mod servo_flap;
mod flap_directional;
mod hover_empirical;
mod hover_cq_empirical;
mod descent_empirical;
mod vpm_empirical;

pub use blade_element_hover::check_blade_element_hover;
pub use hover_castles_gray::check_hover_castles_gray;
pub use climb_momentum::check_climb_momentum;
pub use prandtl_tip_loss::check_prandtl_tip_loss;
pub use glauert_forward_inflow::check_glauert_forward_inflow;
pub use wake_skew::check_wake_skew;
pub use flapping_harmonics::check_flapping_harmonics;
pub use autorotation::check_autorotation;
pub use measured_companions::check_measured_companions;
pub use cyclic_sign::check_cyclic_sign;
pub use servo_flap::check_servo_flap;
pub use flap_directional::check_flap_directional;
pub use hover_empirical::check_hover_empirical;
pub use hover_cq_empirical::check_hover_cq_empirical;
pub use descent_empirical::check_descent_empirical;
pub use vpm_empirical::check_vpm_empirical;

use crate::report::Report;

/// Run all 16 validation checks in order.
pub fn run_all_checks(report: &mut Report) {
    check_blade_element_hover(report);
    check_hover_castles_gray(report);
    check_climb_momentum(report);
    check_prandtl_tip_loss(report);
    check_glauert_forward_inflow(report);
    check_wake_skew(report);
    check_flapping_harmonics(report);
    check_autorotation(report);
    check_measured_companions(report);
    check_cyclic_sign(report);
    check_servo_flap(report);
    check_flap_directional(report);
    check_hover_empirical(report);
    check_hover_cq_empirical(report);
    check_descent_empirical(report);
    check_vpm_empirical(report);
}
