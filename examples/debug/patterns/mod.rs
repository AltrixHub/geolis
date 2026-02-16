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

/// Register a numeric label as a 7-segment display mesh at `(x, y)`.
///
/// `text` may contain digits `0`–`9`; other characters are skipped.
/// `size` controls the height of each digit character.
#[allow(clippy::cast_possible_truncation, clippy::many_single_char_names)]
pub fn register_label(
    storage: &MeshStorage,
    x: f64,
    y: f64,
    text: &str,
    size: f64,
    color: Color,
) {
    let digit_w = size * 0.6;
    let digit_h = size;
    let thickness = size * 0.12;
    let gap = size * 0.2;

    let mut verts_2d: Vec<RawVertex2D> = Vec::new();
    let mut verts_3d: Vec<RawVertex3D> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();

    let mut cursor_x = x;
    for ch in text.chars() {
        let segs = digit_segments(ch);
        if segs == 0 {
            cursor_x += digit_w + gap;
            continue;
        }
        for bit in 0..7u8 {
            if segs & (1 << bit) == 0 {
                continue;
            }
            let (rx, ry, rw, rh) = segment_rect(bit, cursor_x, y, digit_w, digit_h, thickness);
            let base = u32::try_from(verts_2d.len()).unwrap_or(0);

            let min = [rx as f32, ry as f32];
            let max = [(rx + rw) as f32, (ry + rh) as f32];

            verts_2d.push(RawVertex2D::new([min[0], min[1]], [0.0, 0.0]));
            verts_2d.push(RawVertex2D::new([max[0], min[1]], [0.0, 0.0]));
            verts_2d.push(RawVertex2D::new([max[0], max[1]], [0.0, 0.0]));
            verts_2d.push(RawVertex2D::new([min[0], max[1]], [0.0, 0.0]));

            let nrm = [0.0_f32, 0.0, 1.0];
            let uv = [0.0_f32, 0.0];
            verts_3d.push(RawVertex3D { position: [min[0], min[1], 0.0], normal: nrm, uv });
            verts_3d.push(RawVertex3D { position: [max[0], min[1], 0.0], normal: nrm, uv });
            verts_3d.push(RawVertex3D { position: [max[0], max[1], 0.0], normal: nrm, uv });
            verts_3d.push(RawVertex3D { position: [min[0], max[1], 0.0], normal: nrm, uv });

            indices.extend_from_slice(&[base, base + 1, base + 2, base, base + 2, base + 3]);
        }
        cursor_x += digit_w + gap;
    }

    if !verts_2d.is_empty() {
        let mesh_2d = RawMesh2D::new(verts_2d, indices.clone(), color);
        storage.upsert_2d(RawMesh2DId::new(), Arc::new(mesh_2d));
        let mesh_3d = RawMesh3D::new(verts_3d, indices, color);
        storage.upsert_3d(RawMesh3DId::new(), Arc::new(mesh_3d));
    }
}

/// 7-segment bitmask: bit0=a(top), bit1=b(top-right), bit2=c(bottom-right),
/// bit3=d(bottom), bit4=e(bottom-left), bit5=f(top-left), bit6=g(middle).
fn digit_segments(ch: char) -> u8 {
    match ch {
        '0' => 0b0011_1111,
        '1' => 0b0000_0110,
        '2' => 0b0101_1011,
        '3' => 0b0100_1111,
        '4' => 0b0110_0110,
        '5' => 0b0110_1101,
        '6' => 0b0111_1101,
        '7' => 0b0000_0111,
        '8' => 0b0111_1111,
        '9' => 0b0110_1111,
        _ => 0,
    }
}

/// Rectangle `(x, y, width, height)` for a 7-segment segment within a digit cell.
#[allow(clippy::many_single_char_names)]
fn segment_rect(
    seg: u8,
    x: f64,
    y: f64,
    width: f64,
    height: f64,
    thick: f64,
) -> (f64, f64, f64, f64) {
    let half = height * 0.5;
    match seg {
        0 => (x, y + height - thick, width, thick),   // a: top
        1 => (x + width - thick, y + half, thick, half), // b: top-right
        2 => (x + width - thick, y, thick, half),     // c: bottom-right
        3 => (x, y, width, thick),                     // d: bottom
        4 => (x, y, thick, half),                      // e: bottom-left
        5 => (x, y + half, thick, half),               // f: top-left
        6 => (x, y + half - thick * 0.5, width, thick), // g: middle
        _ => (0.0, 0.0, 0.0, 0.0),
    }
}
