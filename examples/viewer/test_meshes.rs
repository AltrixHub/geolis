//! Test mesh data for the viewer example.
//!
//! Replace or extend this module to visualise different meshes.

use std::sync::Arc;

use geolis::math::Point3;
use geolis::tessellation::{StrokeStyle, TessellateStroke, TriangleMesh};
use revion_core::{RawMesh2D, RawMesh2DId, RawMesh3D, RawMesh3DId, RawVertex2D, RawVertex3D};
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

/// Converts a Geolis `TriangleMesh` into a Revion `RawMesh2D`, consuming the mesh.
///
/// Positions are projected onto the XY plane (Z is dropped).
#[allow(clippy::cast_possible_truncation, clippy::needless_pass_by_value)]
fn into_raw_mesh_2d(mesh: TriangleMesh, color: Color) -> RawMesh2D {
    let vertices: Vec<RawVertex2D> = mesh
        .vertices
        .iter()
        .zip(mesh.uvs.iter())
        .map(|(pos, uv)| RawVertex2D::new([pos.x as f32, pos.y as f32], [uv.x as f32, uv.y as f32]))
        .collect();

    let indices: Vec<u32> = mesh
        .indices
        .iter()
        .flat_map(|tri| tri.iter().copied())
        .collect();

    RawMesh2D::new(vertices, indices, color)
}

/// Converts a Geolis `TriangleMesh` into a Revion `RawMesh3D`, consuming the mesh.
#[allow(clippy::cast_possible_truncation, clippy::needless_pass_by_value)]
fn into_raw_mesh_3d(mesh: TriangleMesh, color: Color) -> RawMesh3D {
    let vertices: Vec<RawVertex3D> = mesh
        .vertices
        .iter()
        .zip(mesh.normals.iter())
        .zip(mesh.uvs.iter())
        .map(|((pos, nrm), uv)| RawVertex3D {
            position: [pos.x as f32, pos.y as f32, pos.z as f32],
            normal: [nrm.x as f32, nrm.y as f32, nrm.z as f32],
            uv: [uv.x as f32, uv.y as f32],
        })
        .collect();

    let indices: Vec<u32> = mesh
        .indices
        .iter()
        .flat_map(|tri| tri.iter().copied())
        .collect();

    RawMesh3D::new(vertices, indices, color)
}

/// Register a set of sample meshes into `storage`.
#[allow(clippy::too_many_lines)]
pub fn register_test_meshes(storage: &MeshStorage) {
    // (points, width, closed, color)
    let strokes: &[(&[Point3], f64, bool, Color)] = &[
        // 1) Zigzag — alternating sharp (~30°) angles
        (
            &[
                Point3::new(-8.0, -2.0, 0.0),
                Point3::new(-7.0, 1.0, 0.0),
                Point3::new(-6.0, -2.0, 0.0),
                Point3::new(-5.0, 1.0, 0.0),
                Point3::new(-4.0, -2.0, 0.0),
            ],
            0.2,
            false,
            Color::rgb(255, 80, 80), // red
        ),
        // 2) Very acute hairpin (~20°) — tests miter clamping
        (
            &[
                Point3::new(-3.0, -2.0, 0.0),
                Point3::new(-2.0, 2.0, 0.0),
                Point3::new(-2.8, -1.8, 0.0),
            ],
            0.25,
            false,
            Color::rgb(255, 200, 50), // yellow
        ),
        // 3) Obtuse-angle chain (~135° turns) — gentle bends
        (
            &[
                Point3::new(-1.0, -2.0, 0.0),
                Point3::new(0.0, -1.0, 0.0),
                Point3::new(1.5, -1.2, 0.0),
                Point3::new(2.5, -0.2, 0.0),
                Point3::new(4.0, -0.5, 0.0),
            ],
            0.2,
            false,
            Color::rgb(100, 200, 255), // blue
        ),
        // 4) Right-angle staircase (90° turns)
        (
            &[
                Point3::new(-1.0, 0.5, 0.0),
                Point3::new(0.0, 0.5, 0.0),
                Point3::new(0.0, 1.5, 0.0),
                Point3::new(1.0, 1.5, 0.0),
                Point3::new(1.0, 2.5, 0.0),
                Point3::new(2.0, 2.5, 0.0),
            ],
            0.15,
            false,
            Color::rgb(100, 255, 150), // green
        ),
        // 5) Star shape (closed) — very sharp tips (~36°)
        (
            &[
                Point3::new(5.0, 2.0, 0.0),
                Point3::new(5.6, 0.2, 0.0),
                Point3::new(7.5, 0.2, 0.0),
                Point3::new(6.0, -0.8, 0.0),
                Point3::new(6.6, -2.5, 0.0),
                Point3::new(5.0, -1.5, 0.0),
                Point3::new(3.4, -2.5, 0.0),
                Point3::new(4.0, -0.8, 0.0),
                Point3::new(2.5, 0.2, 0.0),
                Point3::new(4.4, 0.2, 0.0),
            ],
            0.12,
            true,
            Color::rgb(255, 100, 255), // magenta
        ),
        // 6) U-turn (~170° reversal) — near-parallel segments
        (
            &[
                Point3::new(3.0, 1.0, 0.0),
                Point3::new(4.5, 2.5, 0.0),
                Point3::new(3.1, 1.2, 0.0),
            ],
            0.2,
            false,
            Color::rgb(255, 160, 50), // orange
        ),
        // 7) Wide stroke with mixed angles — tests thick miter
        (
            &[
                Point3::new(-8.0, 3.0, 0.0),
                Point3::new(-6.5, 4.5, 0.0),
                Point3::new(-5.0, 3.0, 0.0),
                Point3::new(-3.5, 4.5, 0.0),
                Point3::new(-2.0, 3.0, 0.0),
            ],
            0.5,
            false,
            Color::rgb(180, 130, 255), // purple
        ),
    ];

    for &(points, width, closed, color) in strokes {
        if let Ok(style) = StrokeStyle::new(width) {
            let op = TessellateStroke::new(points.to_vec(), style, closed);
            if let Ok(mesh) = op.execute() {
                storage.upsert_2d(RawMesh2DId::new(), Arc::new(into_raw_mesh_2d(mesh, color)));
            }
            let op = TessellateStroke::new(points.to_vec(), style, closed);
            if let Ok(mesh) = op.execute() {
                storage.upsert_3d(RawMesh3DId::new(), Arc::new(into_raw_mesh_3d(mesh, color)));
            }
        }
    }
}
