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
const OMEGA_MIN:       f64 = 5.0;    // start here -- omega=0 gives mu=inf
const OMEGA_MAX:       f64 = 120.0;  // target rotor speed (rad/s)
const I_ROTOR:         f64 = 0.05;   // effective blade inertia (kg*m^2)
const Q_DRIVE:         f64 = 8.0;    // constant drive torque for spin-up (N*m)
const V_WIND_MAX:      f64 = 10.0;   // final wind speed (m/s)
const WIND_RAMP_STEPS: u64 = 3000;   // steps to full wind (3000/400 = 7.5 s)
const CYCLIC_AMP:      f64 = 0.04;   // swashplate tilt amplitude (rad, ~2.3 deg)
const CYCLIC_PERIOD_S: f64 = 2.0;    // one full cyclic rotation every 2 seconds
const DT:              f64 = 1.0 / 400.0;
const COLLECTIVE_DEG:  f64 = 8.0;
const CROSSWIND_DEG:   f64 = 30.0;

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

fn build_fc(omega: f64, v_wind: f64, tilt_lon: f64, tilt_lat: f64) -> FlightCondition {
    let ang = CROSSWIND_DEG.to_radians();
    FlightCondition {
        collective_rad: COLLECTIVE_DEG.to_radians(),
        tilt_lon,
        tilt_lat,
        v_hub: [v_wind * ang.cos(), v_wind * ang.sin(), 0.0],
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
    v_wind:       f64,   // current wind speed (m/s), ramped 0->V_WIND_MAX
    cyclic_phase: f64,   // angle (rad), advances 2pi every CYCLIC_PERIOD_S
    tilt_lon:     f64,   // current longitudinal swashplate tilt (rad)
    tilt_lat:     f64,   // current lateral swashplate tilt (rad)

    // cached display data, refreshed from state.wake after each step_one call
    particles:    Vec<([f32; 3], f32)>,
    log_min:      f32,
    log_max:      f32,
}

impl VpmVizApp {
    fn new(rotor: VpmRotor<LinearPolar>) -> Self {
        let state = VpmRotorState::default();
        let last_result = VpmRotorResult {
            thrust: 0.0, torque: 0.0, mx_hub: 0.0, my_hub: 0.0,
            n_particles: 0, wake_centroid: [0.0; 3],
        };
        let (particles, log_min, log_max) = extract_particles(&state);
        Self {
            rotor, state, last_result, step_count: 0,
            omega: OMEGA_MIN, v_wind: 0.0,
            cyclic_phase: 0.0, tilt_lon: 0.0, tilt_lat: 0.0,
            particles, log_min, log_max,
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

        let ang = CROSSWIND_DEG as f32 * std::f32::consts::PI / 180.0;
        let ao  = Pos2::new(rect.left() + 40.0, rect.bottom() - 45.0);
        draw_arrow(&painter, ao, ao + Vec2::new(ang.sin() * 50.0, -ang.cos() * 50.0),
            Stroke::new(2.0, Color32::from_rgb(64, 210, 255)));
        painter.text(ao + Vec2::new(ang.sin() * 58.0, -ang.cos() * 58.0),
            egui::Align2::LEFT_CENTER,
            &format!("wind {:.1} m/s  {:.0}deg", self.v_wind, CROSSWIND_DEG),
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

        // Thrust arrow: tilted by tilt_lon in XZ (nose-down = forward = +X screen).
        // Hub-frame thrust direction: X = -sin(tilt_lon), Z = -cos(tilt_lon).
        // Screen XZ: sx = cx + X*scale, sy = cy + Z*scale (+Z is down).
        let tl = self.tilt_lon as f32;
        let arrow_len = 60.0_f32;
        let tip = Pos2::new(cx - tl.sin() * arrow_len, cy - tl.cos() * arrow_len);
        draw_arrow(&painter, Pos2::new(cx, cy), tip,
            Stroke::new(2.0, Color32::from_rgb(255, 200, 60)));
        painter.text(tip + Vec2::new(6.0, 0.0), egui::Align2::LEFT_CENTER,
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
        // Advance cyclic phase (one full rotation per CYCLIC_PERIOD_S).
        self.cyclic_phase += DT * 2.0 * std::f64::consts::PI / CYCLIC_PERIOD_S;
        self.tilt_lon = CYCLIC_AMP * self.cyclic_phase.cos();
        self.tilt_lat = CYCLIC_AMP * self.cyclic_phase.sin();

        // Build flight condition from current dynamic omega/wind, advance one step.
        let fc = build_fc(self.omega, self.v_wind, self.tilt_lon, self.tilt_lat);
        let (result, new_state) = self.rotor.step_one(&fc, &self.state, DT);
        self.state       = new_state;
        self.last_result = result;
        self.step_count += 1;

        // Integrate omega: constant drive torque minus aerodynamic resistive torque.
        let d_omega = (Q_DRIVE - self.last_result.torque) / I_ROTOR * DT;
        self.omega = (self.omega + d_omega).clamp(OMEGA_MIN, OMEGA_MAX);

        // Ramp wind from 0 to V_WIND_MAX over WIND_RAMP_STEPS steps.
        self.v_wind = (V_WIND_MAX * (self.step_count as f64 / WIND_RAMP_STEPS as f64))
            .min(V_WIND_MAX);

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
                        "VPM Rotor  --  {:.0}-deg Crosswind  (V = {:.1} m/s,  omega = {:.1} rad/s,  R = {:.1} m)",
                        CROSSWIND_DEG, self.v_wind, self.omega, R_TIP))
                        .size(17.0).color(Color32::WHITE));
                    ui.label(egui::RichText::new(format!(
                        "step {}   T = {:.1} N   Q = {:.2} N*m   Mx = {:.2} N*m   My = {:.2} N*m   {} particles",
                        self.step_count, self.last_result.thrust, self.last_result.torque,
                        self.last_result.mx_hub, self.last_result.my_hub,
                        self.particles.len()))
                        .size(12.5).color(Color32::from_rgb(175, 175, 175)));
                    ui.label(egui::RichText::new(format!(
                        "cyclic: lon = {:+.2} deg   lat = {:+.2} deg   phase = {:.0} deg",
                        self.tilt_lon.to_degrees(), self.tilt_lat.to_degrees(),
                        self.cyclic_phase.to_degrees() % 360.0))
                        .size(11.5).color(Color32::from_rgb(130, 190, 255)));
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
    println!("  VPM Rotor Visualiser  --  {:.0}-deg Crosswind, cold start", CROSSWIND_DEG);
    println!("  omega: {:.0} -> {:.0} rad/s (drive torque {:.0} N*m, I={:.2} kg*m^2)",
        OMEGA_MIN, OMEGA_MAX, Q_DRIVE, I_ROTOR);
    println!("  wind:  0 -> {:.0} m/s over {} steps ({:.1} s)",
        V_WIND_MAX, WIND_RAMP_STEPS, WIND_RAMP_STEPS as f64 * DT);
    println!("  Opening window...");

    let rotor = build_rotor();

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1020.0, 660.0])
            .with_title("VPM Rotor -- 30-deg Crosswind, cold start (animated)"),
        ..Default::default()
    };

    eframe::run_native(
        "VPM Rotor Visualisation",
        options,
        Box::new(|_cc| Ok(Box::new(VpmVizApp::new(rotor)))),
    )
    .expect("Failed to open window");
}
