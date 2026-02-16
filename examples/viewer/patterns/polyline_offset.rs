//! Polyline offset visualization — original polylines and their offsets.
//!
//! Left column: open polylines (including cross / T / X with 180-degree
//! reversals). Right column: closed polygons (concave cross / T / arrow outlines).

use geolis::math::Point3;
use geolis::operations::offset::PolylineOffset2D;
use geolis::tessellation::StrokeStyle;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::register_stroke;

/// Stroke width for offset result polylines.
const STROKE_WIDTH: f64 = 0.04;
/// Thinner stroke for original (base) polylines.
const BASE_STROKE_WIDTH: f64 = 0.015;

/// Register polyline offset demonstration meshes.
#[allow(clippy::too_many_lines)]
pub fn register(storage: &MeshStorage) {
    // ── Left column: Open polylines ─────────────────────────────────

    // Straight line
    register_offset_pair(
        storage,
        &[
            Point3::new(-12.0, 9.0, 0.0),
            Point3::new(-5.0, 9.0, 0.0),
        ],
        false,
        &[0.4, -0.4],
    );

    // L-shape (90 degree)
    register_offset_pair(
        storage,
        &[
            Point3::new(-12.0, 6.0, 0.0),
            Point3::new(-8.0, 6.0, 0.0),
            Point3::new(-8.0, 8.0, 0.0),
        ],
        false,
        &[0.4, -0.4],
    );

    // U-shape (two 90 degree turns)
    register_offset_pair(
        storage,
        &[
            Point3::new(-12.0, 5.0, 0.0),
            Point3::new(-12.0, 3.0, 0.0),
            Point3::new(-8.0, 3.0, 0.0),
            Point3::new(-8.0, 5.0, 0.0),
        ],
        false,
        &[0.4, -0.4],
    );

    // Zigzag (sharp turns)
    register_offset_pair(
        storage,
        &[
            Point3::new(-12.0, 0.0, 0.0),
            Point3::new(-10.0, 1.5, 0.0),
            Point3::new(-8.0, 0.0, 0.0),
            Point3::new(-6.0, 1.5, 0.0),
            Point3::new(-4.0, 0.0, 0.0),
        ],
        false,
        &[0.3, -0.3],
    );

    // Hairpin (near-180 degree turn)
    register_offset_pair(
        storage,
        &[
            Point3::new(-12.0, -1.5, 0.0),
            Point3::new(-7.0, -1.0, 0.0),
            Point3::new(-12.0, -0.5, 0.0),
        ],
        false,
        &[0.2, -0.2],
    );

    // Staircase (right angles)
    register_offset_pair(
        storage,
        &[
            Point3::new(-12.0, -4.0, 0.0),
            Point3::new(-10.0, -4.0, 0.0),
            Point3::new(-10.0, -3.0, 0.0),
            Point3::new(-8.0, -3.0, 0.0),
            Point3::new(-8.0, -2.0, 0.0),
            Point3::new(-6.0, -2.0, 0.0),
        ],
        false,
        &[0.3, -0.3],
    );

    // S-curve
    register_offset_pair(
        storage,
        &[
            Point3::new(-12.0, -6.0, 0.0),
            Point3::new(-10.0, -5.0, 0.0),
            Point3::new(-8.0, -6.0, 0.0),
            Point3::new(-6.0, -7.0, 0.0),
            Point3::new(-4.0, -6.0, 0.0),
        ],
        false,
        &[0.35, -0.35],
    );

    // Cross / plus — 4 arms from center (180-degree reversals)
    register_offset_pair(
        storage,
        &[
            Point3::new(-10.0, -9.0, 0.0),
            Point3::new(-8.5, -9.0, 0.0),
            Point3::new(-8.5, -7.5, 0.0),
            Point3::new(-8.5, -9.0, 0.0),
            Point3::new(-7.0, -9.0, 0.0),
            Point3::new(-8.5, -9.0, 0.0),
            Point3::new(-8.5, -10.5, 0.0),
        ],
        false,
        &[0.3, -0.3],
    );

    // T-junction — stem meets crossbar (reversal at junction)
    register_offset_pair(
        storage,
        &[
            Point3::new(-4.0, -8.0, 0.0),
            Point3::new(-2.0, -8.0, 0.0),
            Point3::new(-2.0, -9.5, 0.0),
            Point3::new(-2.0, -8.0, 0.0),
            Point3::new(0.0, -8.0, 0.0),
        ],
        false,
        &[0.25, -0.25],
    );

    // W-shape (multiple acute angles)
    register_offset_pair(
        storage,
        &[
            Point3::new(-12.0, -12.0, 0.0),
            Point3::new(-10.5, -10.0, 0.0),
            Point3::new(-9.0, -12.0, 0.0),
            Point3::new(-7.5, -10.0, 0.0),
            Point3::new(-6.0, -12.0, 0.0),
        ],
        false,
        &[0.3, -0.3],
    );

    // ── Right column: Closed polygons ───────────────────────────────

    // Closed square
    register_offset_pair(
        storage,
        &[
            Point3::new(1.0, 7.0, 0.0),
            Point3::new(5.0, 7.0, 0.0),
            Point3::new(5.0, 11.0, 0.0),
            Point3::new(1.0, 11.0, 0.0),
        ],
        true,
        &[0.5, -0.5],
    );

    // Closed triangle
    register_offset_pair(
        storage,
        &[
            Point3::new(7.0, 7.0, 0.0),
            Point3::new(12.0, 7.0, 0.0),
            Point3::new(9.5, 11.0, 0.0),
        ],
        true,
        &[0.5, -0.5],
    );

    // Closed cross / plus outline (12 vertices, 8 concave corners)
    register_offset_pair(
        storage,
        &[
            Point3::new(2.0, 3.5, 0.0),
            Point3::new(3.5, 3.5, 0.0),
            Point3::new(3.5, 2.0, 0.0),
            Point3::new(4.5, 2.0, 0.0),
            Point3::new(4.5, 3.5, 0.0),
            Point3::new(6.0, 3.5, 0.0),
            Point3::new(6.0, 4.5, 0.0),
            Point3::new(4.5, 4.5, 0.0),
            Point3::new(4.5, 6.0, 0.0),
            Point3::new(3.5, 6.0, 0.0),
            Point3::new(3.5, 4.5, 0.0),
            Point3::new(2.0, 4.5, 0.0),
        ],
        true,
        &[0.15, -0.15, 0.3, -0.3],
    );

    // Closed T-shape outline (8 vertices, concave)
    register_offset_pair(
        storage,
        &[
            Point3::new(8.0, 2.0, 0.0),
            Point3::new(13.0, 2.0, 0.0),
            Point3::new(13.0, 3.0, 0.0),
            Point3::new(11.5, 3.0, 0.0),
            Point3::new(11.5, 6.0, 0.0),
            Point3::new(9.5, 6.0, 0.0),
            Point3::new(9.5, 3.0, 0.0),
            Point3::new(8.0, 3.0, 0.0),
        ],
        true,
        &[0.2, -0.2],
    );

    // Closed L-shape (concave)
    register_offset_pair(
        storage,
        &[
            Point3::new(1.0, -2.0, 0.0),
            Point3::new(5.0, -2.0, 0.0),
            Point3::new(5.0, 1.5, 0.0),
            Point3::new(3.5, 1.5, 0.0),
            Point3::new(3.5, -0.5, 0.0),
            Point3::new(1.0, -0.5, 0.0),
        ],
        true,
        &[0.25, -0.25],
    );

    // Closed diamond
    register_offset_pair(
        storage,
        &[
            Point3::new(9.5, -2.0, 0.0),
            Point3::new(12.0, 0.0, 0.0),
            Point3::new(9.5, 2.0, 0.0),
            Point3::new(7.0, 0.0, 0.0),
        ],
        true,
        &[0.4, -0.4],
    );

    // Closed arrow / chevron (concave, thin at tips — use small offset)
    register_offset_pair(
        storage,
        &[
            Point3::new(1.0, -5.0, 0.0),
            Point3::new(5.0, -3.5, 0.0),
            Point3::new(1.0, -2.5, 0.0),
            Point3::new(2.5, -3.5, 0.0),
        ],
        true,
        &[0.1, -0.1],
    );

    // Closed narrow rectangle (wall section, multi-offset)
    register_offset_pair(
        storage,
        &[
            Point3::new(7.0, -4.0, 0.0),
            Point3::new(13.0, -4.0, 0.0),
            Point3::new(13.0, -3.5, 0.0),
            Point3::new(7.0, -3.5, 0.0),
        ],
        true,
        &[0.1, -0.1, 0.2, -0.2],
    );

    // Closed H-shape outline (concave, 12 vertices)
    register_offset_pair(
        storage,
        &[
            Point3::new(1.0, -9.0, 0.0),
            Point3::new(2.0, -9.0, 0.0),
            Point3::new(2.0, -7.5, 0.0),
            Point3::new(4.0, -7.5, 0.0),
            Point3::new(4.0, -9.0, 0.0),
            Point3::new(5.0, -9.0, 0.0),
            Point3::new(5.0, -5.5, 0.0),
            Point3::new(4.0, -5.5, 0.0),
            Point3::new(4.0, -7.0, 0.0),
            Point3::new(2.0, -7.0, 0.0),
            Point3::new(2.0, -5.5, 0.0),
            Point3::new(1.0, -5.5, 0.0),
        ],
        true,
        &[0.15, -0.15],
    );

    // Closed star (5-pointed, sharp angles)
    register_offset_pair(
        storage,
        &[
            Point3::new(9.5, -5.0, 0.0),
            Point3::new(10.1, -6.8, 0.0),
            Point3::new(12.0, -6.8, 0.0),
            Point3::new(10.5, -8.0, 0.0),
            Point3::new(11.1, -9.8, 0.0),
            Point3::new(9.5, -8.8, 0.0),
            Point3::new(7.9, -9.8, 0.0),
            Point3::new(8.5, -8.0, 0.0),
            Point3::new(7.0, -6.8, 0.0),
            Point3::new(8.9, -6.8, 0.0),
        ],
        true,
        &[0.15, -0.15],
    );
}

/// Register the original polyline and its offset variants.
fn register_offset_pair(
    storage: &MeshStorage,
    points: &[Point3],
    closed: bool,
    distances: &[f64],
) {
    let color_original = Color::rgb(180, 180, 180);
    let colors_positive = [Color::rgb(100, 200, 255), Color::rgb(60, 160, 220)];
    let colors_negative = [Color::rgb(255, 130, 100), Color::rgb(220, 100, 70)];

    // Original polyline (gray, thin).
    if let Ok(style) = StrokeStyle::new(BASE_STROKE_WIDTH) {
        register_stroke(storage, points, style, closed, color_original);
    }

    // Offset variants.
    let mut pos_idx = 0_usize;
    let mut neg_idx = 0_usize;
    for &dist in distances {
        let op = PolylineOffset2D::new(points.to_vec(), dist, closed);
        if let Ok(offset_pts) = op.execute() {
            let color = if dist > 0.0 {
                let c = colors_positive[pos_idx % colors_positive.len()];
                pos_idx += 1;
                c
            } else {
                let c = colors_negative[neg_idx % colors_negative.len()];
                neg_idx += 1;
                c
            };
            // Offset results are always closed polygons (open → both-sides outline).
            if let Ok(style) = StrokeStyle::new(STROKE_WIDTH) {
                register_stroke(storage, &offset_pts, style, true, color);
            }
        }
    }
}
