//! Offset self-intersection visualization using `PolylineOffset2D`.
//!
//! Shows closed polygons (T-shape, cross) and open polylines (cross with
//! 180-degree reversals) with algorithmically computed offsets.
//! For hardcoded ground truth, see `offset_intersection_test`.
//!
//! ## Colors
//!
//! - **Gray** (thin): original shape
//! - **Green**: positive offset (inward for closed CCW polygons, left for open)
//! - **Blue**: negative offset (outward for closed CCW polygons, right for open)

use geolis::math::Point3;
use geolis::operations::offset::PolylineOffset2D;
use geolis::tessellation::StrokeStyle;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::register_stroke;

/// Stroke width for original polygons.
const BASE_WIDTH: f64 = 0.015;
/// Stroke width for offset results.
const OFFSET_WIDTH: f64 = 0.04;

/// Register offset test meshes.
pub fn register(storage: &MeshStorage) {
    // ── T-shape cases ────────────────────────────────────────────────
    let t = t_shape_points();

    register_closed_offsets(storage, &t, -12.0, 6.0, &[0.3, -0.3]);
    register_closed_offsets(storage, &t, 2.0, 6.0, &[0.6, -0.6]);
    register_closed_offsets(storage, &t, -12.0, -3.0, &[0.8, -0.8]);
    register_closed_offsets(storage, &t, 2.0, -3.0, &[1.5, -1.5]);

    // ── Closed cross-shape cases ──────────────────────────────────────
    let c = cross_shape_points();

    register_closed_offsets(storage, &c, -12.0, -16.0, &[0.5, -0.5]);
    register_closed_offsets(storage, &c, 2.0, -16.0, &[1.5, -1.5]);

    // ── Open cross with 180-degree reversals ──────────────────────────
    let oc = open_cross_points();

    register_open_offsets(storage, &oc, -12.0, -30.0, &[0.3, -0.3]);
    register_open_offsets(storage, &oc, 2.0, -30.0, &[0.5, -0.5]);
}

/// T-shape CCW vertices at origin.
fn t_shape_points() -> Vec<Point3> {
    vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(10.0, 0.0, 0.0),
        Point3::new(10.0, 1.0, 0.0),
        Point3::new(7.0, 1.0, 0.0),
        Point3::new(7.0, 6.0, 0.0),
        Point3::new(3.0, 6.0, 0.0),
        Point3::new(3.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ]
}

/// Closed cross shape CCW vertices at origin.
fn cross_shape_points() -> Vec<Point3> {
    vec![
        Point3::new(3.0, 0.0, 0.0),
        Point3::new(7.0, 0.0, 0.0),
        Point3::new(7.0, 3.0, 0.0),
        Point3::new(10.0, 3.0, 0.0),
        Point3::new(10.0, 5.0, 0.0),
        Point3::new(7.0, 5.0, 0.0),
        Point3::new(7.0, 10.0, 0.0),
        Point3::new(3.0, 10.0, 0.0),
        Point3::new(3.0, 5.0, 0.0),
        Point3::new(0.0, 5.0, 0.0),
        Point3::new(0.0, 3.0, 0.0),
        Point3::new(3.0, 3.0, 0.0),
    ]
}

/// Open cross / plus — 4 arms from center with 180-degree reversals.
///
/// Path: left arm → center → up arm → center → right arm → center → down arm.
fn open_cross_points() -> Vec<Point3> {
    vec![
        Point3::new(-1.5, 0.0, 0.0),
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 1.5, 0.0),
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.5, 0.0, 0.0),
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, -1.5, 0.0),
    ]
}

/// Translate all points by `(bx, by)`.
fn translate(points: &[Point3], bx: f64, by: f64) -> Vec<Point3> {
    points
        .iter()
        .map(|p| Point3::new(p.x + bx, p.y + by, 0.0))
        .collect()
}

/// Register original closed polygon + computed offsets for given distances.
fn register_closed_offsets(
    storage: &MeshStorage,
    base_points: &[Point3],
    bx: f64,
    by: f64,
    distances: &[f64],
) {
    register_offsets(storage, base_points, bx, by, distances, true);
}

/// Register original open polyline + computed offsets for given distances.
fn register_open_offsets(
    storage: &MeshStorage,
    base_points: &[Point3],
    bx: f64,
    by: f64,
    distances: &[f64],
) {
    register_offsets(storage, base_points, bx, by, distances, false);
}

/// Register original shape + computed offsets for given distances.
fn register_offsets(
    storage: &MeshStorage,
    base_points: &[Point3],
    bx: f64,
    by: f64,
    distances: &[f64],
    closed: bool,
) {
    let translated = translate(base_points, bx, by);

    let color_original = Color::rgb(180, 180, 180);
    let color_positive = Color::rgb(100, 220, 100);
    let color_negative = Color::rgb(80, 140, 255);

    // Original shape (gray, thin).
    if let Ok(style) = StrokeStyle::new(BASE_WIDTH) {
        register_stroke(storage, &translated, style, closed, color_original);
    }

    // Computed offsets.
    // Open polyline offsets now produce closed polygon outlines,
    // so always render the result as closed.
    for &dist in distances {
        let op = PolylineOffset2D::new(translated.clone(), dist, closed);
        if let Ok(offset_pts) = op.execute() {
            let color = if dist > 0.0 {
                color_positive
            } else {
                color_negative
            };
            // Offset results are always closed polygons (open → both-sides outline).
            let result_closed = true;
            if let Ok(style) = StrokeStyle::new(OFFSET_WIDTH) {
                register_stroke(storage, &offset_pts, style, result_closed, color);
            }
        }
    }
}
