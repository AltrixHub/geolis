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
//!
//! ## Line-based shapes (cases 16+)
//!
//! Cases 16 onward use a different pattern: the base is **open line segments**
//! (not a closed polygon), and the offset creates a **closed polygon outline**
//! around those lines. Both d>0 and d<0 produce the same outline.
//!
//! - **Gray**: original base lines (open strokes)
//! - **Green/Blue**: expected closed outline at distance d

use std::f64::consts::SQRT_2;

use geolis::math::Point3;
use geolis::tessellation::StrokeStyle;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_label, register_stroke};

/// Stroke width for original polygons.
const BASE_STROKE_WIDTH: f64 = 0.03;
/// Stroke width for expected (ground truth) offset.
const EXPECTED_STROKE_WIDTH: f64 = 0.03;
/// Label size and color.
const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);

/// Register offset ground-truth test meshes.
pub fn register(storage: &MeshStorage) {
    register_t_cases(storage);
    register_cross_cases(storage);
    register_open_cross_cases(storage);
    register_inverted_t_cases(storage);
    register_comb_cases(storage);
    register_diamond_cases(storage);
    register_vnotch_cases(storage);
    register_x_cross_cases(storage);
    register_double_cross_cases(storage);
    register_fork_cases(storage);
}

fn register_t_cases(storage: &MeshStorage) {
    // Case 1: d=0.3 — bar and stem both survive
    register_label(storage, -14.0, 12.5, "1", LABEL_SIZE, LABEL_COLOR);
    register_t_ground_truth(
        storage, -12.0, 6.0, 0.3,
        &[
            p(-12.0 + 0.3, 6.0 + 0.3),
            p(-12.0 + 9.7, 6.0 + 0.3),
            p(-12.0 + 9.7, 6.0 + 0.7),
            p(-12.0 + 6.7, 6.0 + 0.7),
            p(-12.0 + 6.7, 6.0 + 5.7),
            p(-12.0 + 3.3, 6.0 + 5.7),
            p(-12.0 + 3.3, 6.0 + 0.7),
            p(-12.0 + 0.3, 6.0 + 0.7),
        ],
    );
    // Case 2: d=0.6 — bar collapsed, stem survives
    register_label(storage, 0.0, 12.5, "2", LABEL_SIZE, LABEL_COLOR);
    register_t_ground_truth(
        storage, 2.0, 6.0, 0.6,
        &[
            p(2.0 + 3.6, 6.0 + 0.6),
            p(2.0 + 6.4, 6.0 + 0.6),
            p(2.0 + 6.4, 6.0 + 5.4),
            p(2.0 + 3.6, 6.0 + 5.4),
        ],
    );
    // Case 3: d=0.8 — bar collapsed more, stem still OK
    register_label(storage, -14.0, 3.5, "3", LABEL_SIZE, LABEL_COLOR);
    register_t_ground_truth(
        storage, -12.0, -3.0, 0.8,
        &[
            p(-12.0 + 3.8, -3.0 + 0.8),
            p(-12.0 + 6.2, -3.0 + 0.8),
            p(-12.0 + 6.2, -3.0 + 5.2),
            p(-12.0 + 3.8, -3.0 + 5.2),
        ],
    );
    // Case 4: d=1.5 — bar very collapsed, stem narrow
    register_label(storage, 0.0, 3.5, "4", LABEL_SIZE, LABEL_COLOR);
    register_t_ground_truth(
        storage, 2.0, -3.0, 1.5,
        &[
            p(2.0 + 4.5, -3.0 + 1.5),
            p(2.0 + 5.5, -3.0 + 1.5),
            p(2.0 + 5.5, -3.0 + 4.5),
            p(2.0 + 4.5, -3.0 + 4.5),
        ],
    );
}

fn register_cross_cases(storage: &MeshStorage) {
    // Case 5: d=0.5 — all arms survive
    register_label(storage, -14.0, -5.5, "5", LABEL_SIZE, LABEL_COLOR);
    register_cross_ground_truth(
        storage, -12.0, -16.0, 0.5,
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
    // Case 6: d=1.5 — horizontal arms collapsed
    register_label(storage, 0.0, -5.5, "6", LABEL_SIZE, LABEL_COLOR);
    register_cross_ground_truth(
        storage, 2.0, -16.0, 1.5,
        &[
            p(2.0 + 4.5, -16.0 + 1.5),
            p(2.0 + 5.5, -16.0 + 1.5),
            p(2.0 + 5.5, -16.0 + 8.5),
            p(2.0 + 4.5, -16.0 + 8.5),
        ],
    );
}

fn register_open_cross_cases(storage: &MeshStorage) {
    register_label(storage, -14.0, -28.0, "7", LABEL_SIZE, LABEL_COLOR);
    register_open_cross_ground_truth(storage, -12.0, -30.0, 0.3);
    register_label(storage, 0.0, -28.0, "8", LABEL_SIZE, LABEL_COLOR);
    register_open_cross_ground_truth(storage, 2.0, -30.0, 0.5);
}

fn register_inverted_t_cases(storage: &MeshStorage) {
    // Case 9: d=0.3 — bar and stem both survive (8 vertices)
    register_label(storage, -14.0, -37.5, "9", LABEL_SIZE, LABEL_COLOR);
    register_inverted_t_ground_truth(
        storage, -12.0, -44.0, 0.3,
        &[
            p(-12.0 + 3.3, -44.0 + 0.3),
            p(-12.0 + 6.7, -44.0 + 0.3),
            p(-12.0 + 6.7, -44.0 + 5.3),
            p(-12.0 + 9.7, -44.0 + 5.3),
            p(-12.0 + 9.7, -44.0 + 5.7),
            p(-12.0 + 0.3, -44.0 + 5.7),
            p(-12.0 + 0.3, -44.0 + 5.3),
            p(-12.0 + 3.3, -44.0 + 5.3),
        ],
    );
    // Case 10: d=0.6 — bar collapsed, stem rectangle only
    register_label(storage, -1.5, -37.5, "10", LABEL_SIZE, LABEL_COLOR);
    register_inverted_t_ground_truth(
        storage, 2.0, -44.0, 0.6,
        &[
            p(2.0 + 3.6, -44.0 + 0.6),
            p(2.0 + 6.4, -44.0 + 0.6),
            p(2.0 + 6.4, -44.0 + 5.4),
            p(2.0 + 3.6, -44.0 + 5.4),
        ],
    );
}

fn register_comb_cases(storage: &MeshStorage) {
    // Case 11: d=0.3 — all features survive (12 vertices)
    register_label(storage, -14.0, -50.5, "11", LABEL_SIZE, LABEL_COLOR);
    register_comb_ground_truth(
        storage, -12.0, -58.0, 0.3,
        &[
            p(-12.0 + 0.3, -58.0 + 0.3),
            p(-12.0 + 9.7, -58.0 + 0.3),
            p(-12.0 + 9.7, -58.0 + 0.7),
            p(-12.0 + 8.7, -58.0 + 0.7),
            p(-12.0 + 8.7, -58.0 + 6.7),
            p(-12.0 + 7.3, -58.0 + 6.7),
            p(-12.0 + 7.3, -58.0 + 0.7),
            p(-12.0 + 2.7, -58.0 + 0.7),
            p(-12.0 + 2.7, -58.0 + 3.7),
            p(-12.0 + 1.3, -58.0 + 3.7),
            p(-12.0 + 1.3, -58.0 + 0.7),
            p(-12.0 + 0.3, -58.0 + 0.7),
        ],
    );
    // Case 12: d=0.6 — bar collapses, both teeth survive as separate rectangles
    register_label(storage, -1.5, -50.5, "12", LABEL_SIZE, LABEL_COLOR);
    register_comb_ground_truth(
        storage, 2.0, -58.0, 0.6,
        &[
            p(2.0 + 7.6, -58.0 + 0.6),
            p(2.0 + 8.4, -58.0 + 0.6),
            p(2.0 + 8.4, -58.0 + 6.4),
            p(2.0 + 7.6, -58.0 + 6.4),
        ],
    );
    // Short tooth also survives: column x:1..3, y:0..4, inward by 0.6
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(
            storage,
            &[
                p(2.0 + 1.6, -58.0 + 0.6),
                p(2.0 + 2.4, -58.0 + 0.6),
                p(2.0 + 2.4, -58.0 + 3.4),
                p(2.0 + 1.6, -58.0 + 3.4),
            ],
            style,
            true,
            Color::rgb(100, 220, 100),
        );
    }
}

fn register_diamond_cases(storage: &MeshStorage) {
    // Case 13: d=1.0 — diamond shrinks/expands uniformly (no self-intersection)
    register_label(storage, -14.0, -64.5, "13", LABEL_SIZE, LABEL_COLOR);
    register_diamond_ground_truth(storage, -12.0, -75.0, 1.0);
}

fn register_vnotch_cases(storage: &MeshStorage) {
    // Case 14: d=0.5 — notch survives (7 vertices, diagonal inward offset)
    register_label(storage, 0.0, -70.5, "14", LABEL_SIZE, LABEL_COLOR);
    register_vnotch_ground_truth(
        storage, 2.0, -75.0, 0.5,
        &[
            p(2.0 + 0.5, -75.0 + 0.5),
            p(2.0 + 7.5, -75.0 + 0.5),
            p(2.0 + 7.5, -75.0 + 3.5),
            p(2.0 + 6.0 + 0.5 * (SQRT_2 - 1.0), -75.0 + 3.5),
            p(2.0 + 4.0, -75.0 + 2.0 - 0.5 * SQRT_2),
            p(2.0 + 2.0 + 0.5 * (1.0 - SQRT_2), -75.0 + 3.5),
            p(2.0 + 0.5, -75.0 + 3.5),
        ],
    );
    // Case 15: d=1.0 — notch collapses into bottom → two trapezoids survive
    register_label(storage, 14.0, -70.5, "15", LABEL_SIZE, LABEL_COLOR);
    register_vnotch_ground_truth(
        storage, 16.0, -75.0, 1.0,
        &[
            p(16.0 + 3.0 + SQRT_2, -75.0 + 1.0),
            p(16.0 + 7.0, -75.0 + 1.0),
            p(16.0 + 7.0, -75.0 + 3.0),
            p(16.0 + 5.0 + SQRT_2, -75.0 + 3.0),
        ],
    );
    // Left trapezoid also survives: bounded by x=1, y=1, y=3, left notch wall y=-x+6-√2
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(
            storage,
            &[
                p(16.0 + 1.0, -75.0 + 1.0),
                p(16.0 + 5.0 - SQRT_2, -75.0 + 1.0),
                p(16.0 + 3.0 - SQRT_2, -75.0 + 3.0),
                p(16.0 + 1.0, -75.0 + 3.0),
            ],
            style,
            true,
            Color::rgb(100, 220, 100),
        );
    }
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

// ── Inverted T-shape ──────────────────────────────────────────────

/// Inverted T-shape: stem at bottom (x:3..7, y:0..5), bar at top (x:0..10, y:5..6).
///
/// ```text
///  ┌────────────────-┐
///  │      bar        │  x: 0..10, y: 5..6
///  └─────┬────┬─────-┘
///        │stem│  x: 3..7, y: 0..5
///        │    │
///        └────┘
/// ```
///
/// CCW vertices: (3,0)→(7,0)→(7,5)→(10,5)→(10,6)→(0,6)→(0,5)→(3,5)
fn inverted_t_shape(bx: f64, by: f64) -> [Point3; 8] {
    [
        p(bx + 3.0, by),
        p(bx + 7.0, by),
        p(bx + 7.0, by + 5.0),
        p(bx + 10.0, by + 5.0),
        p(bx + 10.0, by + 6.0),
        p(bx, by + 6.0),
        p(bx, by + 5.0),
        p(bx + 3.0, by + 5.0),
    ]
}

/// Outward offset of the inverted T-shape.
fn inverted_t_outward(bx: f64, by: f64, d: f64) -> [Point3; 8] {
    [
        p(bx + 3.0 - d, by - d),
        p(bx + 7.0 + d, by - d),
        p(bx + 7.0 + d, by + 5.0 - d),
        p(bx + 10.0 + d, by + 5.0 - d),
        p(bx + 10.0 + d, by + 6.0 + d),
        p(bx - d, by + 6.0 + d),
        p(bx - d, by + 5.0 - d),
        p(bx + 3.0 - d, by + 5.0 - d),
    ]
}

/// Register inverted T-shape ground truth.
fn register_inverted_t_ground_truth(
    storage: &MeshStorage,
    bx: f64,
    by: f64,
    d: f64,
    inward: &[Point3],
) {
    let color_original = Color::rgb(180, 180, 180);
    let color_inward = Color::rgb(100, 220, 100);
    let color_outward = Color::rgb(80, 140, 255);

    let original = inverted_t_shape(bx, by);
    if let Ok(style) = StrokeStyle::new(BASE_STROKE_WIDTH) {
        register_stroke(storage, &original, style, true, color_original);
    }
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, inward, style, true, color_inward);
    }
    let outward = inverted_t_outward(bx, by, d);
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, &outward, style, true, color_outward);
    }
}

// ── Comb shape ────────────────────────────────────────────────────

/// Comb: bar x:0..10 y:0..1, left tooth x:1..3 y:1..4, right tooth x:7..9 y:1..7.
///
/// ```text
///            ┌──┐
///  ┌──┐     │  │  right tooth (tall): x:7..9, y:1..7
///  │  │     │  │
///  └──┴─────┴──┘
///     bar: x:0..10, y:0..1
/// ```
///
/// CCW vertices (12):
///   (0,0)→(10,0)→(10,1)→(9,1)→(9,7)→(7,7)→(7,1)→(3,1)→(3,4)→(1,4)→(1,1)→(0,1)
fn comb_shape(bx: f64, by: f64) -> [Point3; 12] {
    [
        p(bx, by),
        p(bx + 10.0, by),
        p(bx + 10.0, by + 1.0),
        p(bx + 9.0, by + 1.0),
        p(bx + 9.0, by + 7.0),
        p(bx + 7.0, by + 7.0),
        p(bx + 7.0, by + 1.0),
        p(bx + 3.0, by + 1.0),
        p(bx + 3.0, by + 4.0),
        p(bx + 1.0, by + 4.0),
        p(bx + 1.0, by + 1.0),
        p(bx, by + 1.0),
    ]
}

/// Outward offset of the comb shape.
fn comb_outward(bx: f64, by: f64, d: f64) -> [Point3; 12] {
    [
        p(bx - d, by - d),
        p(bx + 10.0 + d, by - d),
        p(bx + 10.0 + d, by + 1.0 + d),
        p(bx + 9.0 + d, by + 1.0 + d),
        p(bx + 9.0 + d, by + 7.0 + d),
        p(bx + 7.0 - d, by + 7.0 + d),
        p(bx + 7.0 - d, by + 1.0 + d),
        p(bx + 3.0 + d, by + 1.0 + d),
        p(bx + 3.0 + d, by + 4.0 + d),
        p(bx + 1.0 - d, by + 4.0 + d),
        p(bx + 1.0 - d, by + 1.0 + d),
        p(bx - d, by + 1.0 + d),
    ]
}

/// Register comb shape ground truth.
fn register_comb_ground_truth(
    storage: &MeshStorage,
    bx: f64,
    by: f64,
    d: f64,
    inward: &[Point3],
) {
    let color_original = Color::rgb(180, 180, 180);
    let color_inward = Color::rgb(100, 220, 100);
    let color_outward = Color::rgb(80, 140, 255);

    let original = comb_shape(bx, by);
    if let Ok(style) = StrokeStyle::new(BASE_STROKE_WIDTH) {
        register_stroke(storage, &original, style, true, color_original);
    }
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, inward, style, true, color_inward);
    }
    let outward = comb_outward(bx, by, d);
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, &outward, style, true, color_outward);
    }
}

// ── Diamond (rotated square) ──────────────────────────────────────

/// Diamond: vertices at cardinal directions, side length 5√2.
///
/// ```text
///        (5,10)
///       /      \
///    (0,5)    (10,5)
///       \      /
///        (5,0)
/// ```
///
/// CCW vertices: (5,0)→(10,5)→(5,10)→(0,5)
fn diamond_shape(bx: f64, by: f64) -> [Point3; 4] {
    [
        p(bx + 5.0, by),
        p(bx + 10.0, by + 5.0),
        p(bx + 5.0, by + 10.0),
        p(bx, by + 5.0),
    ]
}

/// Inward offset of the diamond. Each vertex moves inward by d√2.
fn diamond_inward(bx: f64, by: f64, d: f64) -> [Point3; 4] {
    let r = d * SQRT_2;
    [
        p(bx + 5.0, by + r),
        p(bx + 10.0 - r, by + 5.0),
        p(bx + 5.0, by + 10.0 - r),
        p(bx + r, by + 5.0),
    ]
}

/// Outward offset of the diamond.
fn diamond_outward(bx: f64, by: f64, d: f64) -> [Point3; 4] {
    let r = d * SQRT_2;
    [
        p(bx + 5.0, by - r),
        p(bx + 10.0 + r, by + 5.0),
        p(bx + 5.0, by + 10.0 + r),
        p(bx - r, by + 5.0),
    ]
}

/// Register diamond ground truth (inward + outward computed from d).
fn register_diamond_ground_truth(storage: &MeshStorage, bx: f64, by: f64, d: f64) {
    let color_original = Color::rgb(180, 180, 180);
    let color_inward = Color::rgb(100, 220, 100);
    let color_outward = Color::rgb(80, 140, 255);

    let original = diamond_shape(bx, by);
    if let Ok(style) = StrokeStyle::new(BASE_STROKE_WIDTH) {
        register_stroke(storage, &original, style, true, color_original);
    }
    let inward = diamond_inward(bx, by, d);
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, &inward, style, true, color_inward);
    }
    let outward = diamond_outward(bx, by, d);
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, &outward, style, true, color_outward);
    }
}

// ── V-notch rectangle ────────────────────────────────────────────

/// Rectangle 8×4 with a 45° V-notch cut into the top.
///
/// ```text
/// ┌───┐   ┌───┐
/// │   └─V─┘   │  notch: (2,4) → (4,2) → (6,4), 45° walls
/// │           │
/// └───────────┘  rectangle: x:0..8, y:0..4
/// ```
///
/// CCW vertices (7):
///   (0,0)→(8,0)→(8,4)→(6,4)→(4,2)→(2,4)→(0,4)
fn vnotch_shape(bx: f64, by: f64) -> [Point3; 7] {
    [
        p(bx, by),
        p(bx + 8.0, by),
        p(bx + 8.0, by + 4.0),
        p(bx + 6.0, by + 4.0),
        p(bx + 4.0, by + 2.0),
        p(bx + 2.0, by + 4.0),
        p(bx, by + 4.0),
    ]
}

/// Outward offset of the V-notch rectangle.
///
/// Notch walls move outward (left-up / right-up), tip deepens by d√2.
fn vnotch_outward(bx: f64, by: f64, d: f64) -> [Point3; 7] {
    [
        p(bx - d, by - d),
        p(bx + 8.0 + d, by - d),
        p(bx + 8.0 + d, by + 4.0 + d),
        p(bx + 6.0 + d * (1.0 - SQRT_2), by + 4.0 + d),
        p(bx + 4.0, by + 2.0 + d * SQRT_2),
        p(bx + 2.0 + d * (SQRT_2 - 1.0), by + 4.0 + d),
        p(bx - d, by + 4.0 + d),
    ]
}

/// Register V-notch rectangle ground truth.
fn register_vnotch_ground_truth(
    storage: &MeshStorage,
    bx: f64,
    by: f64,
    d: f64,
    inward: &[Point3],
) {
    let color_original = Color::rgb(180, 180, 180);
    let color_inward = Color::rgb(100, 220, 100);
    let color_outward = Color::rgb(80, 140, 255);

    let original = vnotch_shape(bx, by);
    if let Ok(style) = StrokeStyle::new(BASE_STROKE_WIDTH) {
        register_stroke(storage, &original, style, true, color_original);
    }
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, inward, style, true, color_inward);
    }
    let outward = vnotch_outward(bx, by, d);
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, &outward, style, true, color_outward);
    }
}

// ── X-cross (two diagonal crossing lines) ───────────────────────────

fn register_x_cross_cases(storage: &MeshStorage) {
    // Case 16: d=0.5 — two diagonal crossing lines
    register_label(storage, -14.0, -77.0, "16", LABEL_SIZE, LABEL_COLOR);
    register_x_cross_ground_truth(storage, -8.0, -82.0, 3.0, 0.5);
}

/// Register X-cross: two diagonal crossing lines and their closed outline.
///
/// Base: two lines crossing at 90° on a 45° rotation:
///   - SW-NE diagonal: `(cx-a, cy-a)` to `(cx+a, cy+a)`
///   - NW-SE diagonal: `(cx-a, cy+a)` to `(cx+a, cy-a)`
///
/// Outline: 12-vertex polygon (union of two rotated rectangles at distance d).
fn register_x_cross_ground_truth(
    storage: &MeshStorage,
    cx: f64,
    cy: f64,
    a: f64,
    d: f64,
) {
    let color_original = Color::rgb(180, 180, 180);
    let color_positive = Color::rgb(100, 220, 100);
    let color_negative = Color::rgb(80, 140, 255);

    // Gray: two crossing diagonal lines
    let line_sw_ne = [p(cx - a, cy - a), p(cx + a, cy + a)];
    let line_nw_se = [p(cx - a, cy + a), p(cx + a, cy - a)];
    if let Ok(style) = StrokeStyle::new(BASE_STROKE_WIDTH) {
        register_stroke(storage, &line_sw_ne, style, false, color_original);
        register_stroke(storage, &line_nw_se, style, false, color_original);
    }

    // Green/Blue: closed outline (both d>0 and d<0 produce the same shape)
    let outline = x_cross_outline(cx, cy, a, d);
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, &outline, style, true, color_positive);
    }
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, &outline, style, true, color_negative);
    }
}

/// 12-vertex closed outline of two diagonal crossing lines at distance d.
///
/// The outline is the union of two rotated rectangles (one per diagonal arm),
/// each of width 2d centered on the arm line. The center diamond has vertices
/// at `(±D, 0)` and `(0, ±D)` where `D = d√2`.
///
/// Vertex layout (CCW):
/// ```text
///               NW cap
///            v10─────v9
///            /         \
///   SW    v0  v11   v8  v7    NE
///   cap   │     ·       │    cap
///         v1  v2    v5  v6
///            \         /
///            v3─────v4
///               SE cap
/// ```
fn x_cross_outline(cx: f64, cy: f64, a: f64, d: f64) -> [Point3; 12] {
    let h = d * SQRT_2 / 2.0; // d / √2
    let d2 = d * SQRT_2; // 2h
    [
        p(cx - a - h, cy - a + h), //  0: SW cap top-left
        p(cx - a + h, cy - a - h), //  1: SW cap bottom-right
        p(cx, cy - d2),            //  2: center bottom
        p(cx + a - h, cy - a - h), //  3: SE cap bottom-left
        p(cx + a + h, cy - a + h), //  4: SE cap top-right
        p(cx + d2, cy),            //  5: center right
        p(cx + a + h, cy + a - h), //  6: NE cap bottom-right
        p(cx + a - h, cy + a + h), //  7: NE cap top-left
        p(cx, cy + d2),            //  8: center top
        p(cx - a + h, cy + a + h), //  9: NW cap top-right
        p(cx - a - h, cy + a - h), // 10: NW cap bottom-left
        p(cx - d2, cy),            // 11: center left
    ]
}

// ── Double-cross (井 — 4 crossing lines) ─────────────────────────────

fn register_double_cross_cases(storage: &MeshStorage) {
    // Case 17: d=0.3 — all features survive
    register_label(storage, -4.0, -88.5, "17", LABEL_SIZE, LABEL_COLOR);
    register_double_cross_ground_truth(storage, -2.0, -100.0, 0.3);
    // Case 18: d=0.8 — features thinner but all survive
    register_label(storage, 12.0, -88.5, "18", LABEL_SIZE, LABEL_COLOR);
    register_double_cross_ground_truth(storage, 14.0, -100.0, 0.8);
}

/// Register 井 shape: 4 crossing lines and their closed outline.
///
/// Base: 4 line segments forming a hash pattern:
///   - Left vertical:    `x=3`, `y: 0..10`
///   - Right vertical:   `x=7`, `y: 0..10`
///   - Bottom horizontal: `y=3`, `x: 0..10`
///   - Top horizontal:    `y=7`, `x: 0..10`
///
/// Outline: 28-vertex polygon (outer boundary of union of 4 rectangles at
/// distance d). Note: the center region `(3+d, 3+d)` to `(7-d, 7-d)` is a
/// hole not reachable from outside, but represented as interior of the single
/// closed polygon.
fn register_double_cross_ground_truth(
    storage: &MeshStorage,
    bx: f64,
    by: f64,
    d: f64,
) {
    let color_original = Color::rgb(180, 180, 180);
    let color_positive = Color::rgb(100, 220, 100);
    let color_negative = Color::rgb(80, 140, 255);

    // Gray: 4 crossing lines
    let lv = [p(bx + 3.0, by), p(bx + 3.0, by + 10.0)];
    let rv = [p(bx + 7.0, by), p(bx + 7.0, by + 10.0)];
    let bh = [p(bx, by + 3.0), p(bx + 10.0, by + 3.0)];
    let uh = [p(bx, by + 7.0), p(bx + 10.0, by + 7.0)];
    if let Ok(style) = StrokeStyle::new(BASE_STROKE_WIDTH) {
        register_stroke(storage, &lv, style, false, color_original);
        register_stroke(storage, &rv, style, false, color_original);
        register_stroke(storage, &bh, style, false, color_original);
        register_stroke(storage, &uh, style, false, color_original);
    }

    // Green/Blue: closed outline (both d>0 and d<0 produce the same shape)
    let outline = double_cross_outline(bx, by, d);
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, &outline, style, true, color_positive);
    }
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, &outline, style, true, color_negative);
    }

    // Inner hole boundary — the center rectangle not reachable from outside.
    let inner = double_cross_inner_hole(bx, by, d);
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, &inner, style, true, color_positive);
    }
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, &inner, style, true, color_negative);
    }
}

/// 28-vertex closed outline of 4 crossing lines (井) at distance d.
///
/// Each line's buffer is a rectangle of width 2d. The union of 4 such
/// rectangles forms a 28-sided polygon (CCW) with 12 concave corners.
///
/// Lines: LV at `x=3`, RV at `x=7`, BH at `y=3`, UH at `y=7`.
/// All within a 10×10 region `(bx, by)` to `(bx+10, by+10)`.
fn double_cross_outline(bx: f64, by: f64, d: f64) -> [Point3; 28] {
    [
        p(bx + 3.0 - d, by),           //  0: LV bottom-left
        p(bx + 3.0 + d, by),           //  1: LV bottom-right
        p(bx + 3.0 + d, by + 3.0 - d), //  2: LV right → BH bottom
        p(bx + 7.0 - d, by + 3.0 - d), //  3: BH bottom → RV left
        p(bx + 7.0 - d, by),           //  4: RV bottom-left
        p(bx + 7.0 + d, by),           //  5: RV bottom-right
        p(bx + 7.0 + d, by + 3.0 - d), //  6: RV right → BH bottom
        p(bx + 10.0, by + 3.0 - d),    //  7: BH right cap bottom
        p(bx + 10.0, by + 3.0 + d),    //  8: BH right cap top
        p(bx + 7.0 + d, by + 3.0 + d), //  9: BH top → RV right
        p(bx + 7.0 + d, by + 7.0 - d), // 10: RV right → UH bottom
        p(bx + 10.0, by + 7.0 - d),    // 11: UH right cap bottom
        p(bx + 10.0, by + 7.0 + d),    // 12: UH right cap top
        p(bx + 7.0 + d, by + 7.0 + d), // 13: UH top → RV right
        p(bx + 7.0 + d, by + 10.0),    // 14: RV top-right
        p(bx + 7.0 - d, by + 10.0),    // 15: RV top-left
        p(bx + 7.0 - d, by + 7.0 + d), // 16: RV left → UH top
        p(bx + 3.0 + d, by + 7.0 + d), // 17: UH top → LV right
        p(bx + 3.0 + d, by + 10.0),    // 18: LV top-right
        p(bx + 3.0 - d, by + 10.0),    // 19: LV top-left
        p(bx + 3.0 - d, by + 7.0 + d), // 20: LV left → UH top
        p(bx, by + 7.0 + d),           // 21: UH left cap top
        p(bx, by + 7.0 - d),           // 22: UH left cap bottom
        p(bx + 3.0 - d, by + 7.0 - d), // 23: UH bottom → LV left
        p(bx + 3.0 - d, by + 3.0 + d), // 24: LV left → BH top
        p(bx, by + 3.0 + d),           // 25: BH left cap top
        p(bx, by + 3.0 - d),           // 26: BH left cap bottom
        p(bx + 3.0 - d, by + 3.0 - d), // 27: BH bottom → LV left
    ]
}

/// Inner hole of the 井 outline: the center rectangle between the 4 bars.
///
/// For offset d, the hole spans `(3+d, 3+d)` to `(7-d, 7-d)` relative to base.
/// This hole exists when `d < 2.0` (i.e., the vertical bars at x=3,7 and
/// horizontal bars at y=3,7 don't overlap in the center).
fn double_cross_inner_hole(bx: f64, by: f64, d: f64) -> [Point3; 4] {
    [
        p(bx + 3.0 + d, by + 3.0 + d), // bottom-left
        p(bx + 7.0 - d, by + 3.0 + d), // bottom-right
        p(bx + 7.0 - d, by + 7.0 - d), // top-right
        p(bx + 3.0 + d, by + 7.0 - d), // top-left
    ]
}

// ── Fork (又-shape / Y-fork from lines) ─────────────────────────────

fn register_fork_cases(storage: &MeshStorage) {
    // Case 19: d=0.5 — Y-shaped fork from 3 crossing lines
    register_label(storage, -14.0, -103.5, "19", LABEL_SIZE, LABEL_COLOR);
    register_fork_ground_truth(storage, -9.0, -114.0, 0.5);
}

/// Register fork: Y-shaped base lines and their closed outline.
///
/// Base: 3 line segments meeting at junction `(bx+5, by+4)`:
///   - Stem:         `(bx+5, by)` to `(bx+5, by+4)` (vertical)
///   - Left branch:  `(bx+5, by+4)` to `(bx, by+9)` (45° up-left)
///   - Right branch: `(bx+5, by+4)` to `(bx+10, by+9)` (45° up-right)
///
/// Outline: 9-vertex polygon (union of stem rectangle + two rotated branch
/// rectangles at distance d).
fn register_fork_ground_truth(
    storage: &MeshStorage,
    bx: f64,
    by: f64,
    d: f64,
) {
    let color_original = Color::rgb(180, 180, 180);
    let color_positive = Color::rgb(100, 220, 100);
    let color_negative = Color::rgb(80, 140, 255);

    // Gray: Y-shaped base lines (stem + two branches)
    let stem = [p(bx + 5.0, by), p(bx + 5.0, by + 4.0)];
    let left_branch = [p(bx + 5.0, by + 4.0), p(bx, by + 9.0)];
    let right_branch = [p(bx + 5.0, by + 4.0), p(bx + 10.0, by + 9.0)];
    if let Ok(style) = StrokeStyle::new(BASE_STROKE_WIDTH) {
        register_stroke(storage, &stem, style, false, color_original);
        register_stroke(storage, &left_branch, style, false, color_original);
        register_stroke(storage, &right_branch, style, false, color_original);
    }

    // Green/Blue: closed outline (both d>0 and d<0 produce the same shape)
    let outline = fork_outline(bx, by, d);
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, &outline, style, true, color_positive);
    }
    if let Ok(style) = StrokeStyle::new(EXPECTED_STROKE_WIDTH) {
        register_stroke(storage, &outline, style, true, color_negative);
    }
}

/// 9-vertex closed outline of a Y-fork at distance d.
///
/// The stem buffer is a vertical rectangle `(5±d, 0..4)`. Each branch buffer
/// is a rotated rectangle at 45°. The union forms a 9-sided polygon (CCW).
///
/// Junction at `(5, 4)`. Branch directions: left `(-1,1)/√2`, right `(1,1)/√2`.
///
/// Key points:
///   - `jy = 4 + d(1 - √2)`: y where stem side meets branch side (below junction)
///   - `(5, 4 + d√2)`: crotch where inner branch sides intersect (above junction)
///   - `h = d√2/2`: perpendicular offset projected onto axes for 45° lines
fn fork_outline(bx: f64, by: f64, d: f64) -> [Point3; 9] {
    let h = d * SQRT_2 / 2.0;
    let jy = 4.0 + d * (1.0 - SQRT_2);
    [
        p(bx + 5.0 - d, by),            //  0: stem bottom-left
        p(bx + 5.0 + d, by),            //  1: stem bottom-right
        p(bx + 5.0 + d, by + jy),       //  2: stem right → right branch right
        p(bx + 10.0 + h, by + 9.0 - h), //  3: right branch tip right cap
        p(bx + 10.0 - h, by + 9.0 + h), //  4: right branch tip left cap
        p(bx + 5.0, by + 4.0 + d * SQRT_2), // 5: crotch (inner sides intersect)
        p(bx + h, by + 9.0 + h),        //  6: left branch tip right cap
        p(bx - h, by + 9.0 - h),        //  7: left branch tip left cap
        p(bx + 5.0 - d, by + jy),       //  8: stem left → left branch left
    ]
}
