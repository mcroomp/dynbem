// VPM Rotor particle-cloud visualiser -- 30-degree crosswind, animated.
//
// Workflow:
//   1. Settles the free wake synchronously before the window opens
//      (same march as before, prints progress to the terminal).
//   2. Opens a two-panel egui window.
//   3. Each frame calls VpmRotor::step_one, extracts the new particle
//      cloud from the returned state, and redraws -- giving a continuous
//      forward-play animation of the skewed helical wake.
//
// Build & run:
//   cd vpm_viz
//   cargo run --release

use dynbem_rs::cyclic::ControlGains;
use dynbem_rs::polar::LinearPolar;
use dynbem_rs::rotor_definition::{BladeGeometry, LinearPolarParameters, PitchActuation, RotorDefinition};
use dynbem_rs::vpm_rotor::{FlightCondition, VpmRotor, VpmRotorConfig, VpmRotorResult, VpmRotorState};
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
const OMEGA:           f64 = 120.0;
const COLLECTIVE_DEG:  f64 = 8.0;
const CROSSWIND_DEG:   f64 = 30.0;
const V_INPLANE:       f64 = 10.0;

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
        n_steps_per_rev: 18,
        n_wake_rev:      3,
        n_settle_rev:    4,
        sigma:           0.18,
        relax:           0.35,
        nonlinear_lifting_line: true,
        tip_clustering:  true,
        local_core:      true,
        barnes_hut:      false,
        bh_theta:        0.5,
        bh_min_particles: 2048,
    };
    VpmRotor::new(&defn, polar, ctrl, config)
}

fn make_fc() -> FlightCondition {
    let ang = CROSSWIND_DEG.to_radians();
    FlightCondition {
        collective_rad: COLLECTIVE_DEG.to_radians(),
        tilt_lon:       0.0,
        tilt_lat:       0.0,
        v_hub: [V_INPLANE * ang.cos(), V_INPLANE * ang.sin(), 0.0],
        omega_rad_s:    OMEGA,
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
    rotor:       VpmRotor<LinearPolar>,
    fc:          FlightCondition,
    state:       VpmRotorState,
    last_result: VpmRotorResult,
    step_count:  u64,

    // cached display data, refreshed from state.wake after each step_one call
    particles:   Vec<([f32; 3], f32)>,
    log_min:     f32,
    log_max:     f32,
}

impl VpmVizApp {
    fn new(
        rotor:       VpmRotor<LinearPolar>,
        fc:          FlightCondition,
        state:       VpmRotorState,
        last_result: VpmRotorResult,
    ) -> Self {
        let (particles, log_min, log_max) = extract_particles(&state);
        Self { rotor, fc, state, last_result, step_count: 0, particles, log_min, log_max }
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

        painter.circle_stroke(Pos2::new(cx, cy), R_TIP as f32 * scale,
            Stroke::new(1.5, Color32::from_rgba_unmultiplied(220, 220, 220, 110)));
        painter.line_segment(
            [Pos2::new(cx, cy), Pos2::new(cx, cy - R_TIP as f32 * scale)],
            Stroke::new(1.5, Color32::from_rgb(220, 220, 220)));
        painter.line_segment(
            [Pos2::new(cx, cy), Pos2::new(cx, cy + R_TIP as f32 * scale)],
            Stroke::new(1.5, Color32::from_rgb(220, 220, 220)));
        painter.circle_filled(Pos2::new(cx, cy), 4.0, Color32::WHITE);

        let ang = CROSSWIND_DEG as f32 * std::f32::consts::PI / 180.0;
        let ao  = Pos2::new(rect.left() + 40.0, rect.bottom() - 45.0);
        draw_arrow(&painter, ao, ao + Vec2::new(ang.sin() * 50.0, -ang.cos() * 50.0),
            Stroke::new(2.0, Color32::from_rgb(64, 210, 255)));
        painter.text(ao + Vec2::new(ang.sin() * 58.0, -ang.cos() * 58.0),
            egui::Align2::LEFT_CENTER, &format!("wind {:.0}deg", CROSSWIND_DEG),
            FontId::proportional(11.0), Color32::from_rgb(64, 210, 255));

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

        let r_px = R_TIP as f32 * scale;
        painter.line_segment(
            [Pos2::new(cx - r_px, cy), Pos2::new(cx + r_px, cy)],
            Stroke::new(2.5, Color32::from_rgba_unmultiplied(220, 220, 220, 180)));
        painter.circle_filled(Pos2::new(cx, cy), 4.0, Color32::WHITE);

        draw_arrow(&painter, Pos2::new(cx, cy), Pos2::new(cx, cy - 60.0),
            Stroke::new(2.0, Color32::from_rgb(255, 200, 60)));
        painter.text(Pos2::new(cx+6.0, cy-64.0), egui::Align2::LEFT_CENTER,
            &format!("T = {:.0} N", self.last_result.thrust),
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
        // Advance one sub-step and refresh the particle cache.
        let (result, new_state) = self.rotor.step_one(&self.fc, &self.state);
        self.state       = new_state;
        self.last_result = result;
        self.step_count += 1;

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
                        "VPM Rotor  --  {:.0}-deg Crosswind  (V = {:.0} m/s,  omega = {:.0} rad/s,  R = {:.1} m)",
                        CROSSWIND_DEG, V_INPLANE, OMEGA, R_TIP))
                        .size(17.0).color(Color32::WHITE));
                    ui.label(egui::RichText::new(format!(
                        "step {}   T = {:.1} N   Q = {:.2} N*m   Mx = {:.2} N*m   My = {:.2} N*m   {} particles",
                        self.step_count, self.last_result.thrust, self.last_result.torque,
                        self.last_result.mx_hub, self.last_result.my_hub,
                        self.particles.len()))
                        .size(12.5).color(Color32::from_rgb(175, 175, 175)));
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
    println!("  VPM Rotor Visualiser  --  {:.0}-deg Crosswind (animated)", CROSSWIND_DEG);
    println!("  Settling free wake (may take ~30-60 s in release)...");

    let rotor = build_rotor();
    let fc    = make_fc();
    let (result, state) = rotor.march(&fc, None);

    println!("  Settled: T = {:.1} N  Q = {:.2} N*m  {} particles",
        result.thrust, result.torque, result.n_particles);
    println!("  Opening window -- animation playing...");

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1020.0, 660.0])
            .with_title("VPM Rotor -- 30-deg Crosswind (animated)"),
        ..Default::default()
    };

    eframe::run_native(
        "VPM Rotor Visualisation",
        options,
        Box::new(|_cc| Ok(Box::new(VpmVizApp::new(rotor, fc, state, result)))),
    )
    .expect("Failed to open window");
}
