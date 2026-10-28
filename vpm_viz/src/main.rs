// VPM Rotor particle-cloud visualiser -- adjustable wind, interactive.
//
// Workflow:
//   1. Opens a two-panel egui window; the rotor starts at OMEGA_MIN (a
//      floor -- omega=0 gives mu=inf) with zero drive torque.
//   2. Each frame calls VpmRotor::step_one, extracts the new particle
//      cloud from the returned state, and redraws -- giving a continuous
//      forward-play animation of the skewed helical wake.
//   3. Wind, cyclic, collective, and drive torque are live keyboard
//      controls (see the on-screen legend / startup banner): Z/C wind
//      speed (uncapped), A/D wind angle (0=axial, 90=crosswind), arrow
//      keys cyclic, Q/E collective, [ / ] drive torque (uncapped, floor 0).
//
// Build & run:
//   cd vpm_viz
//   cargo run --release

use dynbem_rs::cyclic::ControlGains;
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::rotor_definition::{BladeGeometry, LinearPolarParameters, PitchActuation, RotorDefinition};
use dynbem_rs::vpm::{FlightCondition, VpmRotor, VpmRotorConfig, VpmRotorResult, VpmRotorState};
use eframe::egui::{self, Color32, FontId, Pos2, Rect, Sense, Stroke, Vec2};

// ---------------------------------------------------------------------------
// Rotor / flight parameters
// ---------------------------------------------------------------------------

const R_TIP:           f64 = 1.0;
const R_ROOT:          f64 = 0.2;
const CHORD:           f64 = 0.06;
const TWIST_DEG:       f64 = 2.0;
const CL_ALPHA:        f64 = 5.7;
const CD0:             f64 = 0.01;
const ALPHA_STALL_DEG: f64 = 15.0;
const N_BLADES:        usize = 2;
const N_STATIONS:      usize = 16;
const RHO:             f64 = 1.225;
const OMEGA_MIN:       f64 = 5.0;    // floor -- omega=0 gives mu=inf
const I_ROTOR:         f64 = 0.05;   // effective blade inertia (kg*m^2)
const DT:              f64 = 1.0 / 1200.0;
const COLLECTIVE_DEG:  f64 = 8.0;    // initial collective (deg), Q/E-controlled thereafter
// Azimuth (within the XY/disk plane) of the wind's horizontal component.
// The wind's mix of axial vs horizontal is set by wind_angle_deg (A/D),
// not by this constant -- see build_fc.
const WIND_AZIMUTH_DEG: f64 = 30.0;

// Manual keyboard controls -- rates are per second, applied while the key
// is held so a tap nudges and holding ramps smoothly.
const WIND_RATE_M_S2:        f64 = 40.0;  // wind speed change rate (m/s^2), no upper cap
const WIND_PX_PER_MS:        f32 = 5.0;   // wind arrow scale (px per m/s), clamped on screen
const ANGLE_MIN_DEG:         f64 = 0.0;   // 0 = fully axial (along shaft)
const ANGLE_MAX_DEG:         f64 = 90.0;  // 90 = fully in-plane (crosswind)
const ANGLE_RATE_DEG_S:      f64 = 60.0;  // wind angle change rate (deg/s)
const CYCLIC_MAX_DEG:        f64 = 10.0;  // max swashplate cyclic tilt (deg)
const CYCLIC_RATE_DEG_S:     f64 = 80.0;  // cyclic tilt change rate (deg/s)
const COLLECTIVE_MIN_DEG:    f64 = 0.0;
const COLLECTIVE_MAX_DEG:    f64 = 14.0;
const COLLECTIVE_RATE_DEG_S: f64 = 60.0;  // collective change rate (deg/s)
const DRIVE_TORQUE_RATE_NM_S: f64 = 20.0; // [ / ] drive torque change rate (N*m/s), no upper cap

// ---------------------------------------------------------------------------
// Viridis colormap (22-point approximation)
// ---------------------------------------------------------------------------

fn viridis(t: f32) -> Color32 {
    #[rustfmt::skip]
    const TABLE: &[[f32; 3]] = &[
        [0.267, 0.005, 0.329], [0.278, 0.051, 0.375], [0.280, 0.098, 0.420],
        [0.272, 0.145, 0.459], [0.254, 0.193, 0.496], [0.231, 0.240, 0.526],
        [0.204, 0.287, 0.551], [0.175, 0.333, 0.571], [0.147, 0.381, 0.587],
        [0.121, 0.428, 0.596], [0.100, 0.475, 0.601], [0.087, 0.522, 0.598],
        [0.094, 0.567, 0.587], [0.129, 0.613, 0.565], [0.197, 0.656, 0.534],
        [0.291, 0.697, 0.494], [0.404, 0.734, 0.443], [0.527, 0.767, 0.382],
        [0.654, 0.797, 0.308], [0.782, 0.822, 0.224], [0.906, 0.843, 0.131],
        [0.993, 0.906, 0.144],
    ];
    let t = t.clamp(0.0, 1.0);
    let idx = t * (TABLE.len() - 1) as f32;
    let lo = idx.floor() as usize;
    let hi = (lo + 1).min(TABLE.len() - 1);
    let f = idx - lo as f32;
    let r = (TABLE[lo][0] * (1.0 - f) + TABLE[hi][0] * f) * 255.0;
    let g = (TABLE[lo][1] * (1.0 - f) + TABLE[hi][1] * f) * 255.0;
    let b = (TABLE[lo][2] * (1.0 - f) + TABLE[hi][2] * f) * 255.0;
    Color32::from_rgb(r as u8, g as u8, b as u8)
}

// ---------------------------------------------------------------------------
// Simulation helpers
// ---------------------------------------------------------------------------

fn build_rotor() -> VpmRotor<LinearPolar> {
    let defn = RotorDefinition {
        blade: BladeGeometry {
            n_blades:          N_BLADES,
            radius_m:          R_TIP,
            root_cutout_m:     R_ROOT,
            chord_m:           CHORD,
            twist_deg:         TWIST_DEG,
            n_elements:        N_STATIONS,
            tip_loss:          true,
            r_stations_m:      Vec::new(),
            chord_stations_m:  Vec::new(),
            twist_stations_deg: Vec::new(),
        },
        airfoil: LinearPolarParameters {
            CL0:              0.0,
            CL_alpha_per_rad: CL_ALPHA,
            CD0,
            alpha_stall_deg:  ALPHA_STALL_DEG,
        },
        control:         None,
        pitch_actuation: PitchActuation::DirectMechanical,
        flap:            None,
        name:            "viz_rotor".to_string(),
        description:     "VPM crosswind visualisation rotor".to_string(),
    };
    let polar  = LinearPolar::new(0.0, CL_ALPHA, CD0, ALPHA_STALL_DEG.to_radians());
    let ctrl   = ControlGains::default();
    let config = VpmRotorConfig {
        max_particles:  3_000,
        sigma:           0.18,
        relax:           0.35,
        nonlinear_lifting_line: true,
        tip_clustering:  true,
        local_core:      true,
        barnes_hut:      false,
        bh_theta:        0.5,
        bh_min_particles: 2048,
        flap_dynamics:   true,
        ..VpmRotorConfig::default()
    };
    VpmRotor::new(&defn, polar, ctrl, config)
}

// Builds the flight condition. wind_angle_deg mixes between fully axial
// (0 deg, wind along the shaft/Z axis -- climb/descent-like flow through
// the disk) and fully in-plane crosswind (90 deg, at WIND_AZIMUTH_DEG
// within the disk plane).
fn build_fc(omega: f64, v_wind: f64, wind_angle_deg: f64, tilt_lon: f64, tilt_lat: f64, collective_deg: f64) -> FlightCondition {
    let az    = WIND_AZIMUTH_DEG.to_radians();
    let elev  = wind_angle_deg.to_radians();
    let horiz = v_wind * elev.sin();
    let axial = v_wind * elev.cos();
    FlightCondition {
        collective_rad: collective_deg.to_radians(),
        tilt_lon,
        tilt_lat,
        v_hub: [horiz * az.cos(), horiz * az.sin(), axial],
        omega_rad_s:    omega,
        rho:            RHO,
    }
}

// Extract (pos, log10_alpha_mag) pairs from the current state.
fn extract_particles(state: &VpmRotorState) -> (Vec<([f32; 3], f32)>, f32, f32) {
    let wake = match state.wake.as_ref() {
        Some(w) => w,
        None    => return (Vec::new(), -6.0, -3.0),
    };
    let mut pts = Vec::with_capacity(wake.len());
    let mut lo = f32::MAX;
    let mut hi = f32::MIN;
    for (pos, strength, _) in wake.particles() {
        let mag  = (strength[0]*strength[0] + strength[1]*strength[1] + strength[2]*strength[2]).sqrt();
        let logm = mag.max(1e-12_f32).log10();
        pts.push((pos, logm));
        if logm < lo { lo = logm; }
        if logm > hi { hi = logm; }
    }
    if pts.is_empty() { lo = -6.0; hi = -3.0; }
    (pts, lo, hi)
}

// ---------------------------------------------------------------------------
// Drawing helpers
// ---------------------------------------------------------------------------

fn draw_arrow(painter: &egui::Painter, origin: Pos2, tip: Pos2, stroke: Stroke) {
    painter.line_segment([origin, tip], stroke);
    let dir  = (tip - origin).normalized();
    let perp = Vec2::new(-dir.y, dir.x);
    let head = 10.0_f32;
    painter.line_segment([tip, tip - dir * head + perp * (head * 0.4)], stroke);
    painter.line_segment([tip, tip - dir * head - perp * (head * 0.4)], stroke);
}

// Scale (dx, dy) down (preserving direction) so an arrow drawn from the
// origin never exceeds the given per-direction margins. Keeps the force
// arrows fully on-screen regardless of instantaneous magnitude (thrust
// transients at cold-start otherwise draw past the panel edge).
fn clamp_arrow_extent(dx: f32, dy: f32, max_up: f32, max_down: f32, max_side: f32) -> (f32, f32) {
    let ratio_y = if dy.abs() > 1e-6 {
        (if dy < 0.0 { max_up } else { max_down }) / dy.abs()
    } else {
        f32::INFINITY
    };
    let ratio_x = if dx.abs() > 1e-6 { max_side / dx.abs() } else { f32::INFINITY };
    let ratio = ratio_y.min(ratio_x).min(1.0);
    (dx * ratio, dy * ratio)
}

fn draw_panel_bg(painter: &egui::Painter, rect: Rect, cx: f32, cy: f32, scale: f32) {
    painter.rect_filled(rect, 0.0, Color32::from_rgb(10, 10, 20));
    let gc = Color32::from_rgb(35, 35, 55);
    for i in -5i32..=5 {
        let x = cx + i as f32 * scale;
        let y = cy + i as f32 * scale;
        if x >= rect.left() && x <= rect.right() {
            painter.line_segment([Pos2::new(x, rect.top()), Pos2::new(x, rect.bottom())], Stroke::new(0.5, gc));
        }
        if y >= rect.top() && y <= rect.bottom() {
            painter.line_segment([Pos2::new(rect.left(), y), Pos2::new(rect.right(), y)], Stroke::new(0.5, gc));
        }
    }
    painter.rect_stroke(rect, 0.0, Stroke::new(1.0, Color32::from_rgb(80, 80, 100)));
}

// ---------------------------------------------------------------------------
// App
// ---------------------------------------------------------------------------

struct VpmVizApp {
    rotor:        VpmRotor<LinearPolar>,
    state:        VpmRotorState,
    last_result:  VpmRotorResult,
    step_count:   u64,
    omega:        f64,   // current rotor speed (rad/s), integrated each frame
    v_wind:       f64,   // current wind speed (m/s), manual Z/C control (uncapped)
    wind_angle_deg: f64, // wind elevation (deg): 0=axial, 90=crosswind, manual A/D control
    tilt_lon:     f64,   // current longitudinal swashplate tilt (rad), manual Up/Down control
    tilt_lat:     f64,   // current lateral swashplate tilt (rad), manual Left/Right control
    collective_deg: f64, // current collective pitch (deg), manual Q/E control
    drive_torque: f64,   // current motor drive torque (N*m), manual [/] control, starts at 0

    // cached display data, refreshed from state.wake after each step_one call
    particles:    Vec<([f32; 3], f32)>,
    log_min:      f32,
    log_max:      f32,

    // EMA-smoothed force values used ONLY for the arrow drawings below --
    // the raw per-step result carries real 1/rev (and 2/rev, N_BLADES=2)
    // ripple from the crosswind + cyclic input that would otherwise make
    // the arrows jump every frame. Physics/omega integration always uses
    // the raw self.last_result values, never these.
    thrust_disp:  f64,
    fx_disp:      f64,
    fy_disp:      f64,
}

impl VpmVizApp {
    fn new(rotor: VpmRotor<LinearPolar>) -> Self {
        let state = VpmRotorState::default();
        let last_result = VpmRotorResult {
            thrust: 0.0, torque: 0.0, mx_hub: 0.0, my_hub: 0.0,
            fx_hub: 0.0, fy_hub: 0.0,
            n_particles: 0, wake_centroid: [0.0; 3],
        };
        let (particles, log_min, log_max) = extract_particles(&state);
        Self {
            rotor, state, last_result, step_count: 0,
            omega: OMEGA_MIN, v_wind: 0.0, wind_angle_deg: ANGLE_MIN_DEG,
            tilt_lon: 0.0, tilt_lat: 0.0, collective_deg: COLLECTIVE_DEG,
            drive_torque: 0.0,
            particles, log_min, log_max,
            thrust_disp: 0.0, fx_disp: 0.0, fy_disp: 0.0,
        }
    }

    fn color_for(&self, logm: f32) -> Color32 {
        let range = self.log_max - self.log_min;
        let t = if range > 1e-6 { (logm - self.log_min) / range } else { 0.5 };
        viridis(t)
    }

    // Top-down view (XY).  screen_x = cx + scale*Y,  screen_y = cy - scale*X
    fn draw_top_view(&self, ui: &mut egui::Ui) {
        let (resp, painter) = ui.allocate_painter(Vec2::splat(480.0), Sense::hover());
        let rect  = resp.rect;
        let scale = 110.0_f32;
        let cx    = rect.center().x;
        let cy    = rect.center().y;

        draw_panel_bg(&painter, rect, cx, cy, scale);

        for (pos, logm) in &self.particles {
            let sx = cx + pos[1] * scale;
            let sy = cy - pos[0] * scale;
            if rect.contains(Pos2::new(sx, sy)) {
                painter.circle_filled(Pos2::new(sx, sy), 2.5, self.color_for(*logm));
            }
        }

        // Rotor disk outline (faint circle).
        painter.circle_stroke(Pos2::new(cx, cy), R_TIP as f32 * scale,
            Stroke::new(1.0, Color32::from_rgba_unmultiplied(180, 180, 220, 70)));
        // Rotating blade spans. r_hat(psi) = [cos(psi), -sin(psi), 0];
        // screen: sx = cx + Y*scale = cx - sin(psi)*r*scale
        //         sy = cy - X*scale = cy - cos(psi)*r*scale
        let psi = self.state.psi as f32;
        let r_root_px = R_ROOT as f32 * scale;
        let r_tip_px  = R_TIP  as f32 * scale;
        for b in 0..N_BLADES {
            let psi_b = psi + b as f32 * std::f32::consts::PI;
            let (sp, cp) = psi_b.sin_cos();
            let root = Pos2::new(cx - sp * r_root_px, cy - cp * r_root_px);
            let tip  = Pos2::new(cx - sp * r_tip_px,  cy - cp * r_tip_px);
            painter.line_segment([root, tip], Stroke::new(3.0, Color32::from_rgb(255, 255, 255)));
        }
        painter.circle_filled(Pos2::new(cx, cy), 4.0, Color32::WHITE);

        let ao  = Pos2::new(rect.left() + 40.0, rect.bottom() - 45.0);
        let az    = WIND_AZIMUTH_DEG.to_radians();
        let elev  = self.wind_angle_deg.to_radians();
        let horiz = self.v_wind * elev.sin();
        // Only the horizontal component is visible in this top-down (XY)
        // view -- the axial (Z) component collapses to a point here. Same
        // screen mapping as the particles (sx = cx + Y*scale, sy = cy - X*scale).
        let raw_dx = (horiz * az.sin()) as f32 * WIND_PX_PER_MS;
        let raw_dy = -(horiz * az.cos()) as f32 * WIND_PX_PER_MS;
        let raw_len = raw_dx.hypot(raw_dy);
        let max_len = 90.0_f32;
        let (dx, dy) = if raw_len > max_len && raw_len > 1e-6 {
            (raw_dx * max_len / raw_len, raw_dy * max_len / raw_len)
        } else {
            (raw_dx, raw_dy)
        };
        draw_arrow(&painter, ao, ao + Vec2::new(dx, dy),
            Stroke::new(2.0, Color32::from_rgb(64, 210, 255)));
        painter.text(ao + Vec2::new(dx, dy) + Vec2::new(6.0, 0.0),
            egui::Align2::LEFT_CENTER,
            &format!("wind {:.1} m/s  elev {:.0}deg", self.v_wind, self.wind_angle_deg),
            FontId::proportional(11.0), Color32::from_rgb(64, 210, 255));

        // In-plane hub force (H-force): the horizontal (X,Y) component of the
        // same net rotor force vector whose vertical (Z) component is drawn
        // as "thrust" in the side view. Same screen mapping as the particles
        // above (sx = cx + Y*scale, sy = cy - X*scale), scaled by N/px so the
        // arrow is visible without swamping the disk.
        let fx = self.fx_disp as f32;
        let fy = self.fy_disp as f32;
        let h_mag = fx.hypot(fy);
        if h_mag > 1e-3 {
            let n_per_px = 0.15_f32;
            let max_extent = rect.width().min(rect.height()) * 0.5 - 18.0;
            let (dx, dy) = clamp_arrow_extent(fy / n_per_px, -fx / n_per_px, max_extent, max_extent, max_extent);
            let tip = Pos2::new(cx + dx, cy + dy);
            draw_arrow(&painter, Pos2::new(cx, cy), tip,
                Stroke::new(2.0, Color32::from_rgb(255, 200, 60)));
            painter.text(tip + Vec2::new(6.0, 0.0), egui::Align2::LEFT_CENTER,
                &format!("H = {:.1} N", h_mag),
                FontId::proportional(12.0), Color32::from_rgb(255, 200, 60));
        }

        painter.text(rect.min + Vec2::new(8.0, 8.0), egui::Align2::LEFT_TOP,
            "Top view  (XY, looking -Z)", FontId::proportional(13.0), Color32::from_rgb(200, 200, 200));
        painter.text(Pos2::new(rect.right()-6.0, cy), egui::Align2::RIGHT_CENTER,
            "+Y", FontId::proportional(11.0), Color32::from_rgb(130, 130, 150));
        painter.text(Pos2::new(cx+4.0, rect.top()+24.0), egui::Align2::LEFT_CENTER,
            "+X", FontId::proportional(11.0), Color32::from_rgb(130, 130, 150));
    }

    // Side view (XZ).  screen_x = cx + scale*X,  screen_y = cy + scale*Z
    fn draw_side_view(&self, ui: &mut egui::Ui) {
        let (resp, painter) = ui.allocate_painter(Vec2::splat(480.0), Sense::hover());
        let rect  = resp.rect;
        let scale = 110.0_f32;
        let cx    = rect.center().x;
        let cy    = rect.top() + 110.0;

        draw_panel_bg(&painter, rect, cx, cy, scale);

        for (pos, logm) in &self.particles {
            let sx = cx + pos[0] * scale;
            let sy = cy + pos[2] * scale;
            if rect.contains(Pos2::new(sx, sy)) {
                painter.circle_filled(Pos2::new(sx, sy), 2.5, self.color_for(*logm));
            }
        }

        // Faint disk edge-on line.
        let r_px = R_TIP as f32 * scale;
        painter.line_segment(
            [Pos2::new(cx - r_px, cy), Pos2::new(cx + r_px, cy)],
            Stroke::new(1.0, Color32::from_rgba_unmultiplied(180, 180, 220, 60)));
        // Rotating blades projected onto XZ: Z=0, X = R*cos(psi_b).
        let psi = self.state.psi as f32;
        for b in 0..N_BLADES {
            let psi_b = psi + b as f32 * std::f32::consts::PI;
            let cp = psi_b.cos();
            let root_sx = cx + cp * R_ROOT as f32 * scale;
            let tip_sx  = cx + cp * R_TIP  as f32 * scale;
            painter.line_segment(
                [Pos2::new(root_sx, cy), Pos2::new(tip_sx, cy)],
                Stroke::new(3.0, Color32::from_rgb(255, 255, 255)));
        }
        painter.circle_filled(Pos2::new(cx, cy), 4.0, Color32::WHITE);

        // Wind-direction indicator: the real (X, Z) projection of the wind
        // vector in this XZ (looking +Y) view -- axial component (Z) plus
        // the horizontal component's X-projection (its Y-projection is
        // into/out of the page here, not visible). Scaled + clamped the
        // same way as the top-view wind arrow.
        let az    = WIND_AZIMUTH_DEG.to_radians();
        let elev  = self.wind_angle_deg.to_radians();
        let horiz = self.v_wind * elev.sin();
        let axial = self.v_wind * elev.cos();
        let wind_x = horiz * az.cos();
        let wo     = Pos2::new(rect.left() + 40.0, rect.bottom() - 45.0);
        let raw_dx = wind_x as f32 * WIND_PX_PER_MS;
        let raw_dy = axial as f32 * WIND_PX_PER_MS;
        let raw_len = raw_dx.hypot(raw_dy);
        let max_len = 90.0_f32;
        let (dx, dy) = if raw_len > max_len && raw_len > 1e-6 {
            (raw_dx * max_len / raw_len, raw_dy * max_len / raw_len)
        } else {
            (raw_dx, raw_dy)
        };
        draw_arrow(&painter, wo, wo + Vec2::new(dx, dy),
            Stroke::new(2.0, Color32::from_rgb(64, 210, 255)));
        painter.text(wo + Vec2::new(dx, dy) + Vec2::new(6.0, 0.0),
            egui::Align2::LEFT_CENTER,
            &format!("wind {:.1} m/s  elev {:.0}deg", self.v_wind, self.wind_angle_deg),
            FontId::proportional(11.0), Color32::from_rgb(64, 210, 255));

        // Net force arrow in the XZ plane, from the ACTUAL computed force
        // (not the swashplate tilt_lon command): X = fx_hub (in-plane
        // H-force component, the same quantity arrowed in the top view),
        // Z = -thrust (thrust is "up" = -Z per AGENTS.md). Screen XZ:
        // sx = cx + X*scale, sy = cy + Z*scale (+Z is down).
        let fx = self.fx_disp;
        let fz = -self.thrust_disp;
        let n_per_px = 2.5_f64;
        let max_up   = (cy - rect.top() - 18.0).max(12.0);
        let max_down = (rect.bottom() - cy - 18.0).max(12.0);
        let max_side = (rect.width() * 0.5 - 18.0).max(12.0);
        let (dx, dy) = clamp_arrow_extent(
            (fx / n_per_px) as f32, (fz / n_per_px) as f32, max_up, max_down, max_side);
        let tip = Pos2::new(cx + dx, cy + dy);
        draw_arrow(&painter, Pos2::new(cx, cy), tip,
            Stroke::new(2.0, Color32::from_rgb(255, 200, 60)));
        painter.text(tip + Vec2::new(6.0, 0.0), egui::Align2::LEFT_CENTER,
            &format!("T = {:.0} N", self.thrust_disp),
            FontId::proportional(12.0), Color32::from_rgb(255, 200, 60));

        painter.text(rect.min + Vec2::new(8.0, 8.0), egui::Align2::LEFT_TOP,
            "Side view  (XZ, looking +Y)", FontId::proportional(13.0), Color32::from_rgb(200, 200, 200));
        painter.text(Pos2::new(rect.right()-6.0, cy), egui::Align2::RIGHT_CENTER,
            "+X (fwd)", FontId::proportional(11.0), Color32::from_rgb(130, 130, 150));
        painter.text(Pos2::new(cx+4.0, rect.bottom()-8.0), egui::Align2::LEFT_BOTTOM,
            "+Z (down)", FontId::proportional(11.0), Color32::from_rgb(130, 130, 150));
    }

    fn draw_colorbar(&self, ui: &mut egui::Ui) {
        let bar_w = 340.0_f32;
        let bar_h = 18.0_f32;
        let (resp, painter) = ui.allocate_painter(Vec2::new(bar_w + 130.0, bar_h + 36.0), Sense::hover());
        let rect  = resp.rect;
        let bx0   = rect.left() + 60.0;
        let by0   = rect.top()  + 10.0;
        let bar_rect = Rect::from_min_size(Pos2::new(bx0, by0), Vec2::new(bar_w, bar_h));

        let n_seg = 80usize;
        for i in 0..n_seg {
            let t = i as f32 / (n_seg - 1) as f32;
            painter.rect_filled(
                Rect::from_min_size(Pos2::new(bx0 + i as f32 * bar_w / n_seg as f32, by0),
                    Vec2::new(bar_w / n_seg as f32 + 1.0, bar_h)),
                0.0, viridis(t));
        }
        painter.rect_stroke(bar_rect, 0.0, Stroke::new(1.0, Color32::from_rgb(120, 120, 120)));

        let lc = Color32::from_rgb(150, 150, 150);
        let ft = FontId::proportional(11.0);
        for &(frac, val) in &[(0.0_f32, self.log_min), (0.5, 0.5*(self.log_min+self.log_max)), (1.0, self.log_max)] {
            painter.text(Pos2::new(bx0 + frac * bar_w, by0 + bar_h + 4.0),
                egui::Align2::CENTER_TOP, &format!("{:.1}", val), ft.clone(), lc);
        }
        painter.text(Pos2::new(bx0 + bar_w + 8.0, by0 + bar_h * 0.5),
            egui::Align2::LEFT_CENTER, "log10(|alpha| m^3/s)", FontId::proportional(11.0), lc);
    }
}

impl eframe::App for VpmVizApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Manual controls -- checked every frame so a held key ramps smoothly.
        //   Z / C           -- wind speed down / up (no upper limit)
        //   A / D           -- wind angle: axial (0 deg) <-> crosswind (90 deg)
        //   Up / Down       -- cyclic longitudinal tilt (tilt_lon)
        //   Left / Right    -- cyclic lateral tilt (tilt_lat)
        //   Q / E           -- collective down / up
        //   [ / ]           -- drive torque down / up (starts at 0, no upper limit)
        let cyclic_rate_rad_s = CYCLIC_RATE_DEG_S.to_radians();
        ctx.input(|i| {
            if i.key_down(egui::Key::Z) { self.v_wind -= WIND_RATE_M_S2 * DT; }
            if i.key_down(egui::Key::C) { self.v_wind += WIND_RATE_M_S2 * DT; }
            if i.key_down(egui::Key::A) { self.wind_angle_deg -= ANGLE_RATE_DEG_S * DT; }
            if i.key_down(egui::Key::D) { self.wind_angle_deg += ANGLE_RATE_DEG_S * DT; }
            if i.key_down(egui::Key::ArrowUp)    { self.tilt_lon += cyclic_rate_rad_s * DT; }
            if i.key_down(egui::Key::ArrowDown)  { self.tilt_lon -= cyclic_rate_rad_s * DT; }
            if i.key_down(egui::Key::ArrowRight) { self.tilt_lat += cyclic_rate_rad_s * DT; }
            if i.key_down(egui::Key::ArrowLeft)  { self.tilt_lat -= cyclic_rate_rad_s * DT; }
            if i.key_down(egui::Key::Q) { self.collective_deg -= COLLECTIVE_RATE_DEG_S * DT; }
            if i.key_down(egui::Key::E) { self.collective_deg += COLLECTIVE_RATE_DEG_S * DT; }
            if i.key_down(egui::Key::OpenBracket)  { self.drive_torque -= DRIVE_TORQUE_RATE_NM_S * DT; }
            if i.key_down(egui::Key::CloseBracket) { self.drive_torque += DRIVE_TORQUE_RATE_NM_S * DT; }
        });
        self.v_wind = self.v_wind.max(0.0);
        self.wind_angle_deg = self.wind_angle_deg.clamp(ANGLE_MIN_DEG, ANGLE_MAX_DEG);
        let cyclic_max_rad = CYCLIC_MAX_DEG.to_radians();
        self.tilt_lon = self.tilt_lon.clamp(-cyclic_max_rad, cyclic_max_rad);
        self.tilt_lat = self.tilt_lat.clamp(-cyclic_max_rad, cyclic_max_rad);
        self.collective_deg = self.collective_deg.clamp(COLLECTIVE_MIN_DEG, COLLECTIVE_MAX_DEG);
        self.drive_torque = self.drive_torque.max(0.0);

        // Build flight condition from current dynamic omega/wind, advance one step.
        let fc = build_fc(self.omega, self.v_wind, self.wind_angle_deg, self.tilt_lon, self.tilt_lat, self.collective_deg);
        let (result, new_state) = self.rotor.step_one(&fc, &self.state, DT);
        self.state       = new_state;
        self.last_result = result;
        self.step_count += 1;

        // Smooth the force values for the arrow display only (see field
        // doc comment) -- omega integration below still uses the raw
        // self.last_result.torque.
        let disp_alpha = 0.06;
        self.thrust_disp += (self.last_result.thrust - self.thrust_disp) * disp_alpha;
        self.fx_disp     += (self.last_result.fx_hub  - self.fx_disp)     * disp_alpha;
        self.fy_disp     += (self.last_result.fy_hub  - self.fy_disp)     * disp_alpha;

        // Integrate omega: manual drive torque minus aerodynamic resistive
        // torque, via the single canonical mechanical ODE stepper (no
        // bearing friction modelled in this viz). Floor only -- no upper
        // clamp, so omega free-runs to whatever the torque balance gives.
        self.omega = dynbem_rs::mechanical::step_omega(
            self.omega,
            self.last_result.torque,
            self.drive_torque,
            I_ROTOR,
            DT,
            0.0,
        )
        .max(OMEGA_MIN);

        let (pts, lo, hi) = extract_particles(&self.state);
        self.particles = pts;
        // Gently track the colour range so it does not jump on transients.
        self.log_min = self.log_min * 0.95 + lo * 0.05;
        self.log_max = self.log_max * 0.95 + hi * 0.05;

        ctx.set_visuals(egui::Visuals::dark());

        egui::CentralPanel::default()
            .frame(egui::Frame::none().fill(Color32::from_rgb(8, 8, 16)))
            .show(ctx, |ui| {
                ui.vertical_centered(|ui| {
                    ui.add_space(8.0);
                    ui.label(egui::RichText::new(format!(
                        "VPM Rotor  --  (V = {:.1} m/s,  wind elev = {:.0} deg,  omega = {:.1} rad/s,  R = {:.1} m)",
                        self.v_wind, self.wind_angle_deg, self.omega, R_TIP))
                        .size(17.0).color(Color32::WHITE));
                    ui.label(egui::RichText::new(format!(
                        "step {}   T = {:.1} N   Q = {:.2} N*m   Mx = {:.2} N*m   My = {:.2} N*m   Qdrive = {:.2} N*m   {} particles",
                        self.step_count, self.last_result.thrust, self.last_result.torque,
                        self.last_result.mx_hub, self.last_result.my_hub,
                        self.drive_torque, self.particles.len()))
                        .size(12.5).color(Color32::from_rgb(175, 175, 175)));
                    ui.label(egui::RichText::new(format!(
                        "cyclic: lon = {:+.2} deg   lat = {:+.2} deg   collective = {:.1} deg   wind = {:.1} m/s",
                        self.tilt_lon.to_degrees(), self.tilt_lat.to_degrees(),
                        self.collective_deg, self.v_wind))
                        .size(11.5).color(Color32::from_rgb(130, 190, 255)));
                    ui.label(egui::RichText::new(
                        "controls:  Z/C = wind speed (uncapped)   A/D = wind angle (axial<->crosswind)   arrows = cyclic (up/down lon, left/right lat)   Q/E = collective   [ / ] = drive torque (uncapped)")
                        .size(11.5).color(Color32::from_rgb(150, 210, 150)));
                    ui.add_space(6.0);
                });

                ui.horizontal(|ui| {
                    ui.add_space(8.0);
                    self.draw_top_view(ui);
                    ui.add_space(12.0);
                    self.draw_side_view(ui);
                });

                ui.add_space(6.0);
                ui.vertical_centered(|ui| { self.draw_colorbar(ui); });
            });

        // Keep animating -- request the next frame immediately.
        ctx.request_repaint();
    }
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

fn main() {
    println!("==========================================================");
    println!("  VPM Rotor Visualiser  --  adjustable wind (axial <-> crosswind), cold start");
    println!("  omega: starts at floor {:.0} rad/s (no upper limit), drive torque starts at 0 N*m (I={:.2} kg*m^2)",
        OMEGA_MIN, I_ROTOR);
    println!("  controls:  Z/C = wind speed (uncapped)   A/D = wind angle ({:.0}-{:.0} deg, axial->crosswind)   arrows = cyclic (+/-{:.0} deg)   Q/E = collective ({:.0}-{:.0} deg)   [ / ] = drive torque (uncapped)",
        ANGLE_MIN_DEG, ANGLE_MAX_DEG, CYCLIC_MAX_DEG, COLLECTIVE_MIN_DEG, COLLECTIVE_MAX_DEG);
    println!("  Opening window...");

    let rotor = build_rotor();

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1020.0, 660.0])
            .with_title("VPM Rotor -- adjustable wind, cold start (animated)"),
        ..Default::default()
    };

    eframe::run_native(
        "VPM Rotor Visualisation",
        options,
        Box::new(|_cc| Ok(Box::new(VpmVizApp::new(rotor)))),
    )
    .expect("Failed to open window");
}
