//! Offset self-intersection ground truth — hand-computed expected results only.
//!
//! Draws original polygons and their manually computed correct inward/outward offsets
//! WITHOUT using `PolylineOffset2D`. This establishes visual ground truth
//! before comparing with algorithm output.
//!
//! ## T-shape geometry
//!
//! ```text
//!       ┌────┐
//!       │stem│  width=4 (x: 3..7), height=5 (y: 1..6)
//!       │    │
//! ┌─────┴────┴─────┐
//! │      bar        │  width=10 (x: 0..10), height=1 (y: 0..1)
//! └────────────────-┘
//! ```
//!
//! CCW vertices: (0,0)→(10,0)→(10,1)→(7,1)→(7,6)→(3,6)→(3,1)→(0,1)
//!
//! ## Cross shape geometry
//!
//! ```text
//!          ┌────┐
//!          │    │  vertical arm: x:3..7, y:0..10 (width=4)
//!     ┌────┴────┴────┐
//!     │  horiz arm   │  horizontal arm: x:0..10, y:3..5 (height=2)
//!     └────┬────┬────┘
//!          │    │
//!          └────┘
//! ```
//!
//! CCW vertices: (3,0)→(7,0)→(7,3)→(10,3)→(10,5)→(7,5)→(7,10)→(3,10)→(3,5)→(0,5)→(0,3)→(3,3)
//!
//! ## Colors
//!
//! - **Gray**: original polygon
//! - **Green**: expected correct inward offset (hand-computed vertices)
//! - **Blue**: expected correct outward offset (hand-computed vertices)

use geolis::math::Point3;
use geolis::tessellation::StrokeStyle;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::register_stroke;

/// Stroke width for original polygons.
const BASE_STROKE_WIDTH: f64 = 0.03;
/// Stroke width for expected (ground truth) offset.
const EXPECTED_STROKE_WIDTH: f64 = 0.03;

/// Register offset ground-truth test meshes.
pub fn register(storage: &MeshStorage) {
    // ── T-shape Case 1: d=0.3 — bar and stem both survive ──────────
    register_t_ground_truth(
        storage,
        -12.0,
        6.0,
        0.3,
        // Inward: full T-shape, each edge inset by 0.3
        &[
            p(-12.0 + 0.3, 6.0 + 0.3),       // bar bottom-left
            p(-12.0 + 9.7, 6.0 + 0.3),        // bar bottom-right
            p(-12.0 + 9.7, 6.0 + 0.7),        // bar top-right
            p(-12.0 + 6.7, 6.0 + 0.7),        // junction right
            p(-12.0 + 6.7, 6.0 + 5.7),        // stem top-right
            p(-12.0 + 3.3, 6.0 + 5.7),        // stem top-left
            p(-12.0 + 3.3, 6.0 + 0.7),        // junction left
            p(-12.0 + 0.3, 6.0 + 0.7),        // bar top-left
        ],
    );

    // ── T-shape Case 2: d=0.6 — bar collapsed, stem survives ───────
    register_t_ground_truth(
        storage,
        2.0,
        6.0,
        0.6,
        // Inward: stem rectangle only (bar collapsed)
        &[
            p(2.0 + 3.6, 6.0 + 0.6),         // bottom-left
            p(2.0 + 6.4, 6.0 + 0.6),         // bottom-right
            p(2.0 + 6.4, 6.0 + 5.4),         // top-right
            p(2.0 + 3.6, 6.0 + 5.4),         // top-left
        ],
    );

    // ── T-shape Case 3: d=0.8 — bar collapsed more, stem still OK ──
    register_t_ground_truth(
        storage,
        -12.0,
        -3.0,
        0.8,
        // Inward: stem rectangle only
        &[
            p(-12.0 + 3.8, -3.0 + 0.8),      // bottom-left
            p(-12.0 + 6.2, -3.0 + 0.8),       // bottom-right
            p(-12.0 + 6.2, -3.0 + 5.2),       // top-right
            p(-12.0 + 3.8, -3.0 + 5.2),       // top-left
        ],
    );

    // ── T-shape Case 4: d=1.5 — bar very collapsed, stem narrow ────
    register_t_ground_truth(
        storage,
        2.0,
        -3.0,
        1.5,
        // Inward: narrow stem rectangle
        &[
            p(2.0 + 4.5, -3.0 + 1.5),        // bottom-left
            p(2.0 + 5.5, -3.0 + 1.5),         // bottom-right
            p(2.0 + 5.5, -3.0 + 4.5),         // top-right
            p(2.0 + 4.5, -3.0 + 4.5),         // top-left
        ],
    );

    // ── Cross Case 1: d=0.5 — all arms survive ─────────────────────
    register_cross_ground_truth(
        storage,
        -12.0,
        -16.0,
        0.5,
        // Inward: full 12-vertex cross, each edge inset by 0.5
        &[
            p(-12.0 + 3.5, -16.0 + 0.5),
            p(-12.0 + 6.5, -16.0 + 0.5),
            p(-12.0 + 6.5, -16.0 + 3.5),
            p(-12.0 + 9.5, -16.0 + 3.5),
            p(-12.0 + 9.5, -16.0 + 4.5),
            p(-12.0 + 6.5, -16.0 + 4.5),
            p(-12.0 + 6.5, -16.0 + 9.5),
            p(-12.0 + 3.5, -16.0 + 9.5),
            p(-12.0 + 3.5, -16.0 + 4.5),
            p(-12.0 + 0.5, -16.0 + 4.5),
            p(-12.0 + 0.5, -16.0 + 3.5),
            p(-12.0 + 3.5, -16.0 + 3.5),
        ],
    );

    // ── Cross Case 2: d=1.5 — horizontal arms collapsed ────────────
    register_cross_ground_truth(
        storage,
        2.0,
        -16.0,
        1.5,
        // Inward: vertical arm rectangle only (horiz arm height=2, collapsed)
        &[
            p(2.0 + 4.5, -16.0 + 1.5),
            p(2.0 + 5.5, -16.0 + 1.5),
            p(2.0 + 5.5, -16.0 + 8.5),
            p(2.0 + 4.5, -16.0 + 8.5),
        ],
    );

    // ── Open cross with 180-degree reversals ────────────────────────
    //
    // Open cross centered at (bx, by):
    //   left(-1.5, 0) → center(0,0) → up(0, 1.5) → center → right(1.5, 0) → center → down(0, -1.5)
    //
    // Both d>0 and d<0 produce the same result: a closed 12-vertex cross outline
    // at distance |d| from the original strokes. The offset of an open polyline
    // traces both sides and caps, forming a closed polygon.

    // ── Open cross Case 1: d=0.3 at (-12, -30) ─────────────────────
    register_open_cross_ground_truth(storage, -12.0, -30.0, 0.3);

    // ── Open cross Case 2: d=0.5 at (2, -30) ───────────────────────
    register_open_cross_ground_truth(storage, 2.0, -30.0, 0.5);
}

/// Shorthand to create a `Point3` at z=0.
fn p(x: f64, y: f64) -> Point3 {
    Point3::new(x, y, 0.0)
}

/// T-shape original points (CCW) at base position (bx, by).
fn t_shape(bx: f64, by: f64) -> [Point3; 8] {
    [
        p(bx, by),                // 0: bottom-left
        p(bx + 10.0, by),         // 1: bottom-right
        p(bx + 10.0, by + 1.0),   // 2: bar top-right
        p(bx + 7.0, by + 1.0),    // 3: junction right
        p(bx + 7.0, by + 6.0),    // 4: stem top-right
        p(bx + 3.0, by + 6.0),    // 5: stem top-left
        p(bx + 3.0, by + 1.0),    // 6: junction left
        p(bx, by + 1.0),          // 7: bar top-left
    ]
}

/// Outward offset of the T-shape (each edge moved outward by d).
fn t_shape_outward(bx: f64, by: f64, d: f64) -> [Point3; 8] {
    [
        p(bx - d, by - d),                    // 0: bottom-left
        p(bx + 10.0 + d, by - d),             // 1: bottom-right
        p(bx + 10.0 + d, by + 1.0 + d),       // 2: bar top-right
        p(bx + 7.0 + d, by + 1.0 + d),        // 3: junction right
        p(bx + 7.0 + d, by + 6.0 + d),        // 4: stem top-right
        p(bx + 3.0 - d, by + 6.0 + d),        // 5: stem top-left
        p(bx + 3.0 - d, by + 1.0 + d),        // 6: junction left
        p(bx - d, by + 1.0 + d),              // 7: bar top-left
    ]
}

/// Cross shape original points (CCW) at base position (bx, by).
///
/// ```text
///          ┌────┐
///          │    │  vertical arm: x:3..7, y:0..10
///     ┌────┴────┴────┐
///     │  horiz arm   │  horizontal arm: x:0..10, y:3..5
///     └────┬────┬────┘
///          │    │
///          └────┘
/// ```
fn cross_shape(bx: f64, by: f64) -> [Point3; 12] {
    [
        p(bx + 3.0, by),          // 0
        p(bx + 7.0, by),          // 1
        p(bx + 7.0, by + 3.0),    // 2
        p(bx + 10.0, by + 3.0),   // 3
        p(bx + 10.0, by + 5.0),   // 4
        p(bx + 7.0, by + 5.0),    // 5
        p(bx + 7.0, by + 10.0),   // 6
        p(bx + 3.0, by + 10.0),   // 7
        p(bx + 3.0, by + 5.0),    // 8
        p(bx, by + 5.0),          // 9
        p(bx, by + 3.0),          // 10
        p(bx + 3.0, by + 3.0),    // 11
    ]
}

/// Outward offset of the cross shape (each edge moved outward by d).
fn cross_shape_outward(bx: f64, by: f64, d: f64) -> [Point3; 12] {
    [
        p(bx + 3.0 - d, by - d),
        p(bx + 7.0 + d, by - d),
        p(bx + 7.0 + d, by + 3.0 - d),
        p(bx + 10.0 + d, by + 3.0 - d),
        p(bx + 10.0 + d, by + 5.0 + d),
        p(bx + 7.0 + d, by + 5.0 + d),
        p(bx + 7.0 + d, by + 10.0 + d),
        p(bx + 3.0 - d, by + 10.0 + d),
        p(bx + 3.0 - d, by + 5.0 + d),
        p(bx - d, by + 5.0 + d),
        p(bx - d, by + 3.0 - d),
        p(bx + 3.0 - d, by + 3.0 - d),
    ]
}

/// Register T-shape: original (gray) + inward expected (green) + outward expected (blue).
fn register_t_ground_truth(
    storage: &MeshStorage,
    bx: f64,
    by: f64,
    d: f64,
    inward: &[Point3],
) {
    let color_original = Color::rgb(180, 180, 180);
    let color_inward = Color::rgb(100, 220, 100);
    let color_outward = Color::rgb(80, 140, 255);

    let original = t_shape(bx, by);
    if let Ok(style) = StrokeStyle::new(BASE_STROKE_WIDTH) {
        register_stroke(storage, &original, style, true, color_original);
    }
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, inward, style, true, color_inward);
    }
    let outward = t_shape_outward(bx, by, d);
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, &outward, style, true, color_outward);
    }
}

/// Register cross shape: original (gray) + inward expected (green) + outward expected (blue).
fn register_cross_ground_truth(
    storage: &MeshStorage,
    bx: f64,
    by: f64,
    d: f64,
    inward: &[Point3],
) {
    let color_original = Color::rgb(180, 180, 180);
    let color_inward = Color::rgb(100, 220, 100);
    let color_outward = Color::rgb(80, 140, 255);

    let original = cross_shape(bx, by);
    if let Ok(style) = StrokeStyle::new(BASE_STROKE_WIDTH) {
        register_stroke(storage, &original, style, true, color_original);
    }
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, inward, style, true, color_inward);
    }
    let outward = cross_shape_outward(bx, by, d);
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, &outward, style, true, color_outward);
    }
}

/// Open cross original points at center (cx, cy) with arm length 1.5.
fn open_cross(cx: f64, cy: f64) -> [Point3; 7] {
    [
        p(cx - 1.5, cy),
        p(cx, cy),
        p(cx, cy + 1.5),
        p(cx, cy),
        p(cx + 1.5, cy),
        p(cx, cy),
        p(cx, cy - 1.5),
    ]
}

/// Closed 12-vertex cross outline at distance d from the open cross strokes.
///
/// Both d>0 and d<0 offsets produce this same closed polygon.
///
/// ```text
///     (-d,+1.5)───(+d,+1.5)
///         │           │
///     (-d,+d)─────(+d,+d)
///         │           │
/// (-1.5,-d)           (+1.5,+d)
/// (-1.5,+d)           (+1.5,-d)
///         │           │
///     (-d,-d)─────(+d,-d)
///         │           │
///     (-d,-1.5)───(+d,-1.5)
/// ```
fn open_cross_outline(cx: f64, cy: f64, d: f64) -> [Point3; 12] {
    [
        p(cx - 1.5, cy - d),   //  0: left arm bottom-left
        p(cx - 1.5, cy + d),   //  1: left arm top-left (cap)
        p(cx - d, cy + d),     //  2: center TL
        p(cx - d, cy + 1.5),   //  3: up arm top-left
        p(cx + d, cy + 1.5),   //  4: up arm top-right (cap)
        p(cx + d, cy + d),     //  5: center TR
        p(cx + 1.5, cy + d),   //  6: right arm top-right
        p(cx + 1.5, cy - d),   //  7: right arm bottom-right (cap)
        p(cx + d, cy - d),     //  8: center BR
        p(cx + d, cy - 1.5),   //  9: down arm bottom-right
        p(cx - d, cy - 1.5),   // 10: down arm bottom-left (cap)
        p(cx - d, cy - d),     // 11: center BL
    ]
}

/// Register open cross: original (gray) + d>0 (green) + d<0 (blue).
///
/// Both offsets produce the same closed cross outline, so only one color is visible.
fn register_open_cross_ground_truth(
    storage: &MeshStorage,
    cx: f64,
    cy: f64,
    d: f64,
) {
    let color_original = Color::rgb(180, 180, 180);
    let color_positive = Color::rgb(100, 220, 100);
    let color_negative = Color::rgb(80, 140, 255);

    let original = open_cross(cx, cy);
    if let Ok(style) = StrokeStyle::new(BASE_STROKE_WIDTH) {
        register_stroke(storage, &original, style, false, color_original);
    }

    // d>0 offset → closed cross outline
    let outline_pos = open_cross_outline(cx, cy, d);
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, &outline_pos, style, true, color_positive);
    }

    // d<0 offset → same closed cross outline
    let outline_neg = open_cross_outline(cx, cy, d);
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, &outline_neg, style, true, color_negative);
    }
}

