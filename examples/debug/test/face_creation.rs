use geolis::math::{Point2, Point3, Vector3};
use geolis::tessellation::{StrokeStyle, TriangleMesh};
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_face, register_label, register_stroke};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const GRAY: Color = Color::rgb(180, 180, 180);
const GREEN: Color = Color::rgb(100, 200, 100);
const BLUE: Color = Color::rgb(100, 150, 255);
const UP: Vector3 = Vector3::new(0.0, 0.0, 1.0);

fn make_style() -> Option<StrokeStyle> {
    StrokeStyle::new(0.05).ok()
}

/// Ground truth: triangle → exactly 1 triangle, 3 vertices.
fn triangle_mesh(bx: f64, by: f64) -> TriangleMesh {
    TriangleMesh {
        vertices: vec![
            Point3::new(bx, by, 0.0),
            Point3::new(bx + 4.0, by, 0.0),
            Point3::new(bx + 2.0, by + 4.0, 0.0),
        ],
        normals: vec![UP; 3],
        uvs: vec![Point2::new(0.0, 0.0), Point2::new(4.0, 0.0), Point2::new(2.0, 4.0)],
        indices: vec![[0, 1, 2]],
    }
}

/// Ground truth: square → 2 triangles, 4 vertices.
fn square_mesh(bx: f64, by: f64) -> TriangleMesh {
    TriangleMesh {
        vertices: vec![
            Point3::new(bx, by, 0.0),
            Point3::new(bx + 4.0, by, 0.0),
            Point3::new(bx + 4.0, by + 4.0, 0.0),
            Point3::new(bx, by + 4.0, 0.0),
        ],
        normals: vec![UP; 4],
        uvs: vec![
            Point2::new(0.0, 0.0),
            Point2::new(4.0, 0.0),
            Point2::new(4.0, 4.0),
            Point2::new(0.0, 4.0),
        ],
        indices: vec![[0, 1, 2], [0, 2, 3]],
    }
}

pub fn register(storage: &MeshStorage) {
    let Some(style) = make_style() else {
        return;
    };

    // Case 1: Triangle
    let bx = -12.0;
    let by = 2.0;
    register_label(storage, bx - 2.0, by + 5.5, "1", LABEL_SIZE, LABEL_COLOR);
    let tri_pts = [
        Point3::new(bx, by, 0.0),
        Point3::new(bx + 4.0, by, 0.0),
        Point3::new(bx + 2.0, by + 4.0, 0.0),
    ];
    register_stroke(storage, &tri_pts, style, true, GRAY);
    register_face(storage, triangle_mesh(bx, by), GREEN);

    // Case 2: Square
    let bx = -4.0;
    let by = 2.0;
    register_label(storage, bx - 2.0, by + 5.5, "2", LABEL_SIZE, LABEL_COLOR);
    let sq_pts = [
        Point3::new(bx, by, 0.0),
        Point3::new(bx + 4.0, by, 0.0),
        Point3::new(bx + 4.0, by + 4.0, 0.0),
        Point3::new(bx, by + 4.0, 0.0),
    ];
    register_stroke(storage, &sq_pts, style, true, GRAY);
    register_face(storage, square_mesh(bx, by), BLUE);

    // Case 3: L-shape (concave) — 4 triangles, 6 vertices
    let bx = 4.0;
    let by = 2.0;
    register_label(storage, bx - 2.0, by + 5.5, "3", LABEL_SIZE, LABEL_COLOR);
    let l_pts = [
        Point3::new(bx, by, 0.0),
        Point3::new(bx + 4.0, by, 0.0),
        Point3::new(bx + 4.0, by + 2.0, 0.0),
        Point3::new(bx + 2.0, by + 2.0, 0.0),
        Point3::new(bx + 2.0, by + 4.0, 0.0),
        Point3::new(bx, by + 4.0, 0.0),
    ];
    register_stroke(storage, &l_pts, style, true, GRAY);
    let l_mesh = TriangleMesh {
        vertices: vec![
            Point3::new(bx, by, 0.0),
            Point3::new(bx + 4.0, by, 0.0),
            Point3::new(bx + 4.0, by + 2.0, 0.0),
            Point3::new(bx + 2.0, by + 2.0, 0.0),
            Point3::new(bx + 2.0, by + 4.0, 0.0),
            Point3::new(bx, by + 4.0, 0.0),
        ],
        normals: vec![UP; 6],
        uvs: vec![
            Point2::new(0.0, 0.0),
            Point2::new(4.0, 0.0),
            Point2::new(4.0, 2.0),
            Point2::new(2.0, 2.0),
            Point2::new(2.0, 4.0),
            Point2::new(0.0, 4.0),
        ],
        indices: vec![[0, 1, 2], [0, 2, 3], [0, 3, 5], [3, 4, 5]],
    };
    register_face(storage, l_mesh, GREEN);

    // Case 4: Square with hole — outer 6x6, inner 3x3 centered
    let bx = 14.0;
    let by = 2.0;
    register_label(storage, bx - 2.0, by + 7.5, "4", LABEL_SIZE, LABEL_COLOR);
    let outer_pts = [
        Point3::new(bx, by, 0.0),
        Point3::new(bx + 6.0, by, 0.0),
        Point3::new(bx + 6.0, by + 6.0, 0.0),
        Point3::new(bx, by + 6.0, 0.0),
    ];
    let inner_pts = [
        Point3::new(bx + 1.5, by + 1.5, 0.0),
        Point3::new(bx + 4.5, by + 1.5, 0.0),
        Point3::new(bx + 4.5, by + 4.5, 0.0),
        Point3::new(bx + 1.5, by + 4.5, 0.0),
    ];
    register_stroke(storage, &outer_pts, style, true, GRAY);
    register_stroke(storage, &inner_pts, style, true, GRAY);
    // For ground truth, we show the outline only (mesh is complex to hand-compute)
    // The algorithm output pattern (Step 3) will show the actual tessellation.
}
