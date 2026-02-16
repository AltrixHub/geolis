pub mod basic_strokes;
pub mod offset_intersection;
pub mod polyline_offset;
pub mod stroke_joins;

use std::sync::Arc;

use geolis::math::Point3;
use geolis::tessellation::{StrokeStyle, TessellateStroke, TriangleMesh};
use revion_core::{RawMesh2D, RawMesh2DId, RawMesh3D, RawMesh3DId, RawVertex2D, RawVertex3D};
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

/// All available pattern names.
pub const PATTERNS: &[&str] = &[
    "stroke_joins",
    "basic_strokes",
    "polyline_offset",
    "offset_intersection",
];

/// Register meshes for the named pattern. Returns `true` if found.
pub fn register(storage: &MeshStorage, name: &str) -> bool {
    match name {
        "stroke_joins" => {
            stroke_joins::register(storage);
            true
        }
        "basic_strokes" => {
            basic_strokes::register(storage);
            true
        }
        "polyline_offset" => {
            polyline_offset::register(storage);
            true
        }
        "offset_intersection" => {
            offset_intersection::register(storage);
            true
        }
        _ => false,
    }
}

// ── Shared utilities ────────────────────────────────────────────────

/// Converts a Geolis `TriangleMesh` into a Revion `RawMesh2D`.
#[allow(clippy::cast_possible_truncation, clippy::needless_pass_by_value)]
pub fn into_raw_mesh_2d(mesh: TriangleMesh, color: Color) -> RawMesh2D {
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

/// Converts a Geolis `TriangleMesh` into a Revion `RawMesh3D`.
#[allow(clippy::cast_possible_truncation, clippy::needless_pass_by_value)]
pub fn into_raw_mesh_3d(mesh: TriangleMesh, color: Color) -> RawMesh3D {
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

/// Tessellate a stroke and register both 2D and 3D meshes.
pub fn register_stroke(
    storage: &MeshStorage,
    points: &[Point3],
    style: StrokeStyle,
    closed: bool,
    color: Color,
) {
    let op = TessellateStroke::new(points.to_vec(), style, closed);
    if let Ok(mesh) = op.execute() {
        storage.upsert_2d(RawMesh2DId::new(), Arc::new(into_raw_mesh_2d(mesh, color)));
    }
    let op = TessellateStroke::new(points.to_vec(), style, closed);
    if let Ok(mesh) = op.execute() {
        storage.upsert_3d(RawMesh3DId::new(), Arc::new(into_raw_mesh_3d(mesh, color)));
    }
}
