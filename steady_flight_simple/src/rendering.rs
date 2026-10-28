use bevy::prelude::*;
use bevy::render::mesh::{Indices, PrimitiveTopology};
use bevy::render::render_asset::RenderAssetUsages;

use crate::simulation::{hub_rotation, SimRunner};

const WORLD_SCALE: f32 = 0.3;
const HUB_RADIUS: f32 = 0.2 * WORLD_SCALE;
const BLADE_THICKNESS: f32 = 0.04;
const DISK_THICKNESS: f32 = 0.02;
const BLADE_WIDTH: f32 = 0.06;
const TETHER_THICKNESS: f32 = 0.03;
const AXLE_OFFSET_M: f32 = 0.4;
const ROTOR_R_ROOT_M: f32 = 0.5;
const ROTOR_R_TIP_M: f32 = 2.5;
const BLADE_ATTACH_R_M: f32 = 0.0;
const BLADE_VISIBLE_SPAN_M: f32 = ROTOR_R_TIP_M - BLADE_ATTACH_R_M;
const ROTOR_N_BLADES: usize = 4;

// Wake particle cloud rendering. The whole cloud is a single mesh of
// camera-facing billboard quads: one opaque draw call, per-vertex colours, and
// each particle is a small screen-facing square (visible and cheap). Rebuilt on
// the CPU each frame -- ~5-6k quads is trivial for the GPU.
const MAX_PARTICLES: usize = 6000;
// Billboard half-size in Bevy world units (constant world size).
const PARTICLE_HALF_SIZE: f32 = 0.015;

#[derive(Component)]
pub(crate) struct RotorRoot;

#[derive(Component)]
pub(crate) struct Blade {
    phase_rad: f32,
}

#[derive(Component)]
pub(crate) struct TetherVisual;

#[derive(Component)]
struct DiskOuter;

#[derive(Component)]
struct DiskCutout;

#[derive(Component)]
pub(crate) struct FollowCamera;

// Camera framing: a side-on view. The eye sits off to one side (along the Bevy
// +Z / East axis) at roughly hub height so the rotor disk is seen edge-on and
// the tether/wake are viewed from the side.
const CAM_DISTANCE: f32 = ROTOR_R_TIP_M * WORLD_SCALE * 4.8;
const CAM_ELEVATION: f32 = ROTOR_R_TIP_M * WORLD_SCALE * 0.15;

#[derive(Component)]
pub(crate) struct WakeCloud;

fn ned_to_bevy(v: [f64; 3]) -> Vec3 {
    // NED -> Bevy (X forward, Y up, Z right-ish)
    Vec3::new(
        (v[0] as f32) * WORLD_SCALE,
        (-(v[2] as f32)) * WORLD_SCALE,
        (v[1] as f32) * WORLD_SCALE,
    )
}

fn ned_dir_to_bevy(v: [f64; 3]) -> Vec3 {
    // Direction-only NED -> Bevy mapping (no world scale).
    Vec3::new(v[0] as f32, -(v[2] as f32), v[1] as f32)
}

fn rotor_quat_in_bevy(roll: f64, pitch: f64) -> Quat {
    let r_hub = hub_rotation(roll, pitch);

    // R_hub columns are hub basis vectors expressed in NED.
    let hub_x_ned = [r_hub.0[0][0], r_hub.0[1][0], r_hub.0[2][0]];
    let hub_y_ned = [r_hub.0[0][1], r_hub.0[1][1], r_hub.0[2][1]];
    let hub_z_ned = [r_hub.0[0][2], r_hub.0[1][2], r_hub.0[2][2]];

    // Rotor mesh local frame convention:
    // local +X -> hub +X
    // local +Y -> hub +Z (spin axis)
    // local +Z -> hub -Y
    let local_x_world = ned_dir_to_bevy(hub_x_ned).normalize_or_zero();
    let local_y_world = ned_dir_to_bevy(hub_z_ned).normalize_or_zero();
    let local_z_world =
        ned_dir_to_bevy([-hub_y_ned[0], -hub_y_ned[1], -hub_y_ned[2]]).normalize_or_zero();

    let rot = bevy::math::Mat3::from_cols(local_x_world, local_y_world, local_z_world);
    Quat::from_mat3(&rot)
}

pub fn setup_scene(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    commands.spawn(PointLightBundle {
        point_light: PointLight {
            intensity: 80000.0,
            shadows_enabled: true,
            ..default()
        },
        transform: Transform::from_xyz(18.0, 20.0, 12.0),
        ..default()
    });

    commands.spawn((
        FollowCamera,
        Camera3dBundle {
            transform: Transform::from_xyz(0.0, CAM_ELEVATION, CAM_DISTANCE)
                .looking_at(Vec3::ZERO, Vec3::Y),
            ..default()
        },
    ));

    // Origin marker.
    commands.spawn(PbrBundle {
        mesh: meshes.add(Mesh::from(Cuboid::from_size(Vec3::splat(0.25)))),
        material: materials.add(StandardMaterial {
            base_color: Color::srgb(1.0, 0.8, 0.2),
            ..default()
        }),
        ..default()
    });

    // Ground plane for depth cue.
    commands.spawn(PbrBundle {
        mesh: meshes.add(Mesh::from(Cuboid::from_size(Vec3::new(120.0, 0.1, 120.0)))),
        material: materials.add(StandardMaterial {
            base_color: Color::srgb(0.08, 0.10, 0.11),
            perceptual_roughness: 0.95,
            ..default()
        }),
        transform: Transform::from_xyz(0.0, -0.05, 0.0),
        ..default()
    });

    let rotor_root = commands.spawn((RotorRoot, SpatialBundle::default())).id();

    commands.entity(rotor_root).with_children(|parent| {
        // Windpower-style annulus disk: blue translucent outer disk plus inner cutout.
        parent.spawn((
            DiskOuter,
            PbrBundle {
                mesh: meshes.add(Mesh::from(Cylinder {
                    radius: ROTOR_R_TIP_M * WORLD_SCALE,
                    half_height: DISK_THICKNESS * 0.5,
                })),
                material: materials.add(StandardMaterial {
                    base_color: Color::srgba(0.27, 0.53, 1.0, 0.20),
                    alpha_mode: AlphaMode::Blend,
                    ..default()
                }),
                ..default()
            },
        ));

        parent.spawn((
            DiskCutout,
            PbrBundle {
                mesh: meshes.add(Mesh::from(Cylinder {
                    radius: ROTOR_R_ROOT_M * WORLD_SCALE,
                    half_height: DISK_THICKNESS,
                })),
                material: materials.add(StandardMaterial {
                    base_color: Color::srgb(0.09, 0.11, 0.13),
                    ..default()
                }),
                ..default()
            },
        ));

        // Hub sphere (windpower-style: light-grey sphere, radius 0.2 m).
        parent.spawn(PbrBundle {
            mesh: meshes.add(Mesh::from(Sphere::new(HUB_RADIUS))),
            material: materials.add(StandardMaterial {
                base_color: Color::srgb(0.82, 0.82, 0.84),
                ..default()
            }),
            ..default()
        });

        // Rotor axis marker (green), similar to windpower body_z arrow cue.
        parent.spawn(PbrBundle {
            mesh: meshes.add(Mesh::from(Cuboid::from_size(Vec3::new(0.04, 1.2, 0.04)))),
            material: materials.add(StandardMaterial {
                base_color: Color::srgb(0.0, 1.0, 0.27),
                ..default()
            }),
            transform: Transform::from_xyz(0.0, 0.6, 0.0),
            ..default()
        });

        // Windpower-style per-blade coloring and root->tip geometry.
        let blade_colors = [
            Color::srgb(1.0, 0.44, 0.26),
            Color::srgb(0.30, 0.82, 0.88),
            Color::srgb(1.0, 0.84, 0.31),
            Color::srgb(0.81, 0.58, 0.85),
        ];
        for i in 0..ROTOR_N_BLADES {
            let phase = (i as f32) * std::f32::consts::FRAC_PI_2;
            parent.spawn((
                Blade { phase_rad: phase },
                PbrBundle {
                    mesh: meshes.add(Mesh::from(Cuboid::from_size(Vec3::new(
                        BLADE_VISIBLE_SPAN_M * WORLD_SCALE,
                        BLADE_THICKNESS,
                        BLADE_WIDTH,
                    )))),
                    material: materials.add(StandardMaterial {
                        base_color: blade_colors[i],
                        ..default()
                    }),
                    ..default()
                },
            ));
        }
    });

    commands.spawn((
        TetherVisual,
        PbrBundle {
            mesh: meshes.add(Mesh::from(Cuboid::from_size(Vec3::ONE))),
            material: materials.add(StandardMaterial {
                base_color: Color::srgb(0.88, 0.88, 0.90),
                ..default()
            }),
            ..default()
        },
    ));

    // The wake is a single billboard-quad mesh: one opaque draw call, per-vertex
    // colours. It starts empty and is rewritten each frame by `update_particles`.
    let cloud_mesh = meshes.add(empty_point_cloud());
    let cloud_material = materials.add(StandardMaterial {
        base_color: Color::WHITE,
        unlit: true,
        cull_mode: None,
        double_sided: true,
        ..default()
    });
    commands.spawn((
        WakeCloud,
        PbrBundle {
            mesh: cloud_mesh,
            material: cloud_material,
            ..default()
        },
    ));
}

/// An empty mesh with `TriangleList` topology kept in the main world so we can
/// rewrite its vertex attributes each frame with camera-facing billboard quads.
fn empty_point_cloud() -> Mesh {
    let mut mesh = Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::MAIN_WORLD | RenderAssetUsages::RENDER_WORLD,
    );
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, Vec::<[f32; 3]>::new());
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, Vec::<[f32; 3]>::new());
    mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, Vec::<[f32; 2]>::new());
    mesh.insert_attribute(Mesh::ATTRIBUTE_COLOR, Vec::<[f32; 4]>::new());
    mesh.insert_indices(Indices::U32(Vec::new()));
    mesh
}

/// Map a normalized strength `t` in [0, 1] to a colour: blue (weak) -> cyan ->
/// green -> yellow -> red (strong). A simple piecewise "jet"-style ramp.
fn strength_colormap(t: f32) -> [f32; 4] {
    let t = t.clamp(0.0, 1.0);
    let (r, g, b) = if t < 0.25 {
        let s = t / 0.25;
        (0.0, s, 1.0)
    } else if t < 0.5 {
        let s = (t - 0.25) / 0.25;
        (0.0, 1.0, 1.0 - s)
    } else if t < 0.75 {
        let s = (t - 0.5) / 0.25;
        (s, 1.0, 0.0)
    } else {
        let s = (t - 0.75) / 0.25;
        (1.0, 1.0 - s, 0.0)
    };
    [r, g, b, 1.0]
}

pub fn update_visuals(
    runner: Res<SimRunner>,
    mut q: ParamSet<(
        Query<&mut Transform, With<RotorRoot>>,
        Query<(&Blade, &mut Transform)>,
        Query<&mut Transform, With<TetherVisual>>,
    )>,
) {
    let sim = &runner.sim;
    let hub_pos = ned_to_bevy(sim.position_ned());
    let rotor_q = rotor_quat_in_bevy(sim.roll_rad(), sim.pitch_rad());
    let axis_world = rotor_q * Vec3::Y;

    if let Ok(mut rotor_tf) = q.p0().get_single_mut() {
        rotor_tf.translation = hub_pos;
        rotor_tf.rotation = rotor_q;
    }

    // Spin blades around local Y and root them at the hub center so the box
    // extends straight outward along its own (rotated) long axis.
    let blade_center_r = (BLADE_ATTACH_R_M + BLADE_VISIBLE_SPAN_M * 0.5) * WORLD_SCALE;
    for (blade, mut tf) in &mut q.p1() {
        let a = sim.rotor_phase_rad() as f32 + blade.phase_rad;
        let spin = Quat::from_rotation_y(a);
        tf.rotation = spin;
        tf.translation = (spin * Vec3::X) * blade_center_r;
    }

    if let Ok(mut tether_tf) = q.p2().get_single_mut() {
        let a = Vec3::ZERO;
        // Attach at axle-bottom-like offset instead of hub center.
        let b = hub_pos - axis_world.normalize_or_zero() * (AXLE_OFFSET_M * WORLD_SCALE);
        let d = b - a;
        let len = d.length().max(1e-6);
        let mid = a + d * 0.5;
        let dir = d / len;
        tether_tf.translation = mid;
        tether_tf.rotation = Quat::from_rotation_arc(Vec3::Y, dir);
        tether_tf.scale = Vec3::new(TETHER_THICKNESS, len, TETHER_THICKNESS);
    }
}

pub fn follow_camera(runner: Res<SimRunner>, mut cam_q: Query<&mut Transform, With<FollowCamera>>) {
    let sim = &runner.sim;
    let hub_pos = ned_to_bevy(sim.position_ned());

    if let Ok(mut cam_tf) = cam_q.get_single_mut() {
        // Keep a fixed 3/4 offset relative to the hub so it stays large on screen.
        // Side-on view: offset along the horizontal +Z axis at hub height.
        let eye = hub_pos + Vec3::new(0.0, CAM_ELEVATION, CAM_DISTANCE);
        cam_tf.translation = eye;
        cam_tf.look_at(hub_pos, Vec3::Y);
    }
}

/// Map a hub-frame position (metres) into the rotor's local Bevy frame, matching
/// the axis convention used by `rotor_quat_in_bevy`:
/// local +X <- hub +X, local +Y <- hub +Z, local +Z <- hub -Y.
fn hub_point_to_local_bevy(p: [f32; 3]) -> Vec3 {
    Vec3::new(p[0], p[2], -p[1]) * WORLD_SCALE
}

pub fn update_particles(
    runner: Res<SimRunner>,
    cloud_q: Query<&Handle<Mesh>, With<WakeCloud>>,
    cam_q: Query<&Transform, With<FollowCamera>>,
    mut meshes: ResMut<Assets<Mesh>>,
) {
    let sim = &runner.sim;
    let particles = sim.wake_particles_hub();
    let hub_pos = ned_to_bevy(sim.position_ned());
    let rotor_q = rotor_quat_in_bevy(sim.roll_rad(), sim.pitch_rad());

    let Ok(mesh_handle) = cloud_q.get_single() else {
        return;
    };
    let Some(mesh) = meshes.get_mut(mesh_handle) else {
        return;
    };

    // Camera-facing basis so each quad always faces the viewer.
    let (right, up) = match cam_q.get_single() {
        Ok(t) => (Vec3::from(t.right()), Vec3::from(t.up())),
        Err(_) => (Vec3::X, Vec3::Y),
    };
    let rx = right * PARTICLE_HALF_SIZE;
    let uy = up * PARTICLE_HALF_SIZE;
    let normal = {
        let n = right.cross(up).normalize_or_zero();
        [n.x, n.y, n.z]
    };

    let n = particles.len().min(MAX_PARTICLES);

    // Normalize strength by the current maximum so the colourmap always spans
    // the live range (with a small floor to avoid divide-by-zero early on).
    let max_mag = particles
        .iter()
        .map(|(_pos, mag)| *mag)
        .fold(0.0_f32, f32::max)
        .max(1e-6);

    let mut positions = Vec::with_capacity(n * 4);
    let mut normals = Vec::with_capacity(n * 4);
    let mut uvs = Vec::with_capacity(n * 4);
    let mut colors = Vec::with_capacity(n * 4);
    let mut indices = Vec::with_capacity(n * 6);

    for (i, &(pos, mag)) in particles.iter().take(n).enumerate() {
        let c = hub_pos + rotor_q * hub_point_to_local_bevy(pos);
        let col = strength_colormap(mag / max_mag);

        let corners = [c - rx - uy, c + rx - uy, c + rx + uy, c - rx + uy];
        for corner in corners {
            positions.push([corner.x, corner.y, corner.z]);
            normals.push(normal);
            colors.push(col);
        }
        uvs.push([0.0, 0.0]);
        uvs.push([1.0, 0.0]);
        uvs.push([1.0, 1.0]);
        uvs.push([0.0, 1.0]);

        let base = (i * 4) as u32;
        indices.extend_from_slice(&[base, base + 1, base + 2, base, base + 2, base + 3]);
    }

    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, positions);
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, normals);
    mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, uvs);
    mesh.insert_attribute(Mesh::ATTRIBUTE_COLOR, colors);
    mesh.insert_indices(Indices::U32(indices));
}
