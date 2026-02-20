//! `WallOutline2D` algorithm output — wall outline from centerline networks.

use geolis::geometry::pline::{Pline, PlineVertex};
use geolis::math::Point3;
use geolis::operations::offset::WallOutline2D;
use geolis::tessellation::StrokeStyle;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_label, register_stroke};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const GRAY: Color = Color::rgb(160, 160, 160);
const GREEN: Color = Color::rgb(100, 220, 100);
const BLUE: Color = Color::rgb(100, 150, 255);

// ── Centerline definitions ──────────────────────────────────────────

fn single_line() -> Vec<(f64, f64)> {
    vec![(0.0, 0.0), (5.0, 0.0)]
}

fn l_shape() -> Vec<(f64, f64)> {
    vec![(0.0, 0.0), (5.0, 0.0), (5.0, 5.0)]
}

fn t_shape() -> Vec<(f64, f64)> {
    vec![(0.0, 3.0), (10.0, 3.0), (5.0, 3.0), (5.0, 0.0)]
}

fn cross_shape() -> Vec<(f64, f64)> {
    vec![(5.0, 0.0), (5.0, 10.0), (5.0, 5.0), (0.0, 5.0), (10.0, 5.0)]
}

fn double_cross() -> Vec<(f64, f64)> {
    vec![
        (3.0, 0.0), (3.0, 10.0), (3.0, 7.0), (0.0, 7.0), (10.0, 7.0),
        (7.0, 7.0), (7.0, 10.0), (7.0, 0.0), (7.0, 3.0), (10.0, 3.0), (0.0, 3.0),
    ]
}

fn y_fork() -> Vec<(f64, f64)> {
    let sin60 = std::f64::consts::FRAC_PI_3.sin();
    let cos60 = std::f64::consts::FRAC_PI_3.cos();
    let (cx, cy, arm) = (5.0, 5.0, 5.0);
    vec![
        (cx, cy + arm),
        (cx, cy),
        (cx - arm * sin60, cy - arm * cos60),
        (cx, cy),
        (cx + arm * sin60, cy - arm * cos60),
    ]
}

fn h_shape() -> Vec<(f64, f64)> {
    vec![
        (3.0, 0.0), (3.0, 10.0),
        (3.0, 5.0), (7.0, 5.0),
        (7.0, 0.0), (7.0, 10.0),
    ]
}

fn e_shape() -> Vec<(f64, f64)> {
    vec![
        (0.0, 3.0), (12.0, 3.0),
        (2.0, 3.0), (2.0, 8.0),
        (2.0, 3.0), (6.0, 3.0), (6.0, 8.0),
        (6.0, 3.0), (10.0, 3.0), (10.0, 8.0),
    ]
}

fn angled_cross() -> Vec<(f64, f64)> {
    let sin60 = std::f64::consts::FRAC_PI_3.sin();
    let (cx, cy, arm) = (5.0, 5.0, 5.0);
    vec![
        (0.0, 5.0), (10.0, 5.0),
        (cx, cy), (cx - arm * 0.5, cy - arm * sin60),
        (cx, cy), (cx + arm * 0.5, cy + arm * sin60),
    ]
}

fn square_room() -> Vec<(f64, f64)> {
    vec![(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
}

fn rectangle_room() -> Vec<(f64, f64)> {
    vec![(0.0, 0.0), (12.0, 0.0), (12.0, 8.0), (0.0, 8.0)]
}

fn l_room() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0), (5.0, 0.0), (5.0, 3.0),
        (3.0, 3.0), (3.0, 5.0), (0.0, 5.0),
    ]
}

/// Closed square room with corridor extending from bottom at x=5.
fn room_with_corridor() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0), (5.0, 0.0), (5.0, -5.0), (5.0, 0.0),
        (10.0, 0.0), (10.0, 10.0), (0.0, 10.0),
    ]
}

/// Closed square room divided by horizontal partition at y=5.
fn room_with_partition() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0), (10.0, 0.0), (10.0, 5.0), (0.0, 5.0),
        (10.0, 5.0), (10.0, 10.0), (0.0, 10.0),
    ]
}

/// Closed room with a diagonal wall from (-5,0) to (15,10) as a single line.
fn room_with_diagonal_wall() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0), (10.0, 0.0), (10.0, 7.5),
        (15.0, 10.0), (-5.0, 0.0),
        (0.0, 2.5), (0.0, 10.0), (10.0, 10.0),
        (10.0, 7.5), (0.0, 2.5),
    ]
}

/// Closed room with a wall penetrating through both sides at y=5.
fn room_with_penetrating_wall() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0), (10.0, 0.0), (10.0, 5.0), (13.0, 5.0),
        (10.0, 5.0), (0.0, 5.0), (10.0, 5.0),
        (10.0, 10.0), (0.0, 10.0), (0.0, 5.0),
        (-3.0, 5.0), (0.0, 5.0),
    ]
}

fn t_very_short_arm() -> Vec<(f64, f64)> {
    vec![(0.0, 3.0), (8.0, 3.0), (4.0, 3.0), (4.0, 3.5)]
}

fn t_arm_2d() -> Vec<(f64, f64)> {
    vec![(0.0, 3.0), (8.0, 3.0), (4.0, 3.0), (4.0, 5.0)]
}

fn cross_short() -> Vec<(f64, f64)> {
    vec![(4.0, 2.0), (4.0, 6.0), (4.0, 4.0), (2.0, 4.0), (6.0, 4.0)]
}

fn l_large_d() -> Vec<(f64, f64)> {
    vec![(0.0, 0.0), (2.0, 0.0), (2.0, 4.0)]
}

fn t_arm_eq_d() -> Vec<(f64, f64)> {
    vec![(0.0, 3.0), (8.0, 3.0), (4.0, 3.0), (4.0, 4.0)]
}

fn l_shape_45() -> Vec<(f64, f64)> {
    vec![(0.0, 0.0), (5.0, 0.0), (8.0, 3.0)]
}

fn t_diagonal_branch() -> Vec<(f64, f64)> {
    vec![(0.0, 3.0), (10.0, 3.0), (5.0, 3.0), (7.0, 5.0)]
}

fn y_mixed_junction() -> Vec<(f64, f64)> {
    vec![(0.0, 0.0), (5.0, 0.0), (5.0, 5.0), (5.0, 0.0), (8.0, 3.0)]
}

fn room_with_corner_stub() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0), (-3.0, -3.0), (0.0, 0.0),
        (8.0, 0.0), (8.0, 8.0), (0.0, 8.0),
    ]
}

fn room_with_corner_diagonal() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0), (10.0, 0.0), (10.0, 8.0), (8.0, 8.0),
        (11.0, 11.0), (-3.0, -3.0),
        (0.0, 0.0), (0.0, 8.0), (8.0, 8.0),
    ]
}

fn room_with_near_corner_stub() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.5), (-3.0, -2.5), (0.0, 0.5),
        (0.0, 8.0), (8.0, 8.0), (8.0, 0.0), (0.0, 0.0),
    ]
}

fn room_with_near_corner_diagonal() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0), (10.0, 0.0), (10.0, 8.0), (7.5, 8.0),
        (10.5, 11.0), (-3.0, -2.5),
        (0.0, 0.5), (0.0, 8.0), (7.5, 8.0),
        (0.0, 0.5),
    ]
}

// ── Drawing helper ──────────────────────────────────────────────────

fn draw_case(
    storage: &MeshStorage,
    pts: &[(f64, f64)],
    half_w: f64,
    closed: bool,
    bx: f64,
    by: f64,
) {
    // Centerline in gray.
    let center: Vec<Point3> = pts
        .iter()
        .map(|&(x, y)| Point3::new(x + bx, y + by, 0.0))
        .collect();
    if let Ok(s) = StrokeStyle::new(0.05) {
        register_stroke(storage, &center, s, closed, GRAY);
    }

    // Algorithm output.
    let pline = Pline {
        vertices: pts.iter().map(|&(x, y)| PlineVertex::line(x, y)).collect(),
        closed,
    };
    let wall = WallOutline2D::new(vec![pline], half_w);
    if let Ok(outlines) = wall.execute() {
        for (i, ol) in outlines.iter().enumerate() {
            let p: Vec<Point3> = ol
                .vertices
                .iter()
                .map(|v| Point3::new(v.x + bx, v.y + by, 0.0))
                .collect();
            let color = if i == 0 { GREEN } else { BLUE };
            if let Ok(s) = StrokeStyle::new(0.08) {
                register_stroke(storage, &p, s, ol.closed, color);
            }
        }
    }
}

// ── Registration ────────────────────────────────────────────────────

/// Register `wall_offset` pattern meshes.
#[allow(clippy::too_many_lines, clippy::type_complexity)]
pub fn register(storage: &MeshStorage) {
    // (centerline, half_width, closed, base_x, base_y, label_x, label_y)
    let cases: Vec<(Vec<(f64, f64)>, f64, bool, f64, f64, f64, f64)> = vec![
        (single_line(), 0.3, false, 0.0, 0.0, -1.5, 1.5),
        (l_shape(), 0.3, false, 10.0, 0.0, 8.5, 6.0),
        (t_shape(), 0.3, false, 22.0, 0.0, 20.5, 4.5),
        (cross_shape(), 0.3, false, 36.0, 0.0, 34.5, 11.5),
        (double_cross(), 0.3, false, 0.0, -16.0, -1.5, -4.5),
        (double_cross(), 0.8, false, 16.0, -16.0, 14.5, -4.5),
        (y_fork(), 0.5, false, 32.0, -16.0, 30.5, -4.5),
        (h_shape(), 0.3, false, 0.0, -32.0, -1.5, -20.5),
        (e_shape(), 0.3, false, 16.0, -32.0, 14.5, -22.5),
        (angled_cross(), 0.3, false, 32.0, -32.0, 30.5, -20.5),
        (square_room(), 0.3, true, 0.0, -48.0, -1.5, -36.5),
        (rectangle_room(), 0.3, true, 16.0, -48.0, 14.5, -36.5),
        (l_room(), 0.3, true, 32.0, -48.0, 30.5, -41.5),
        (room_with_corridor(), 0.3, true, 0.0, -64.0, -1.5, -52.5),
        (room_with_partition(), 0.3, true, 16.0, -64.0, 14.5, -52.5),
        (room_with_penetrating_wall(), 0.3, true, 32.0, -64.0, 30.5, -52.5),
        (room_with_diagonal_wall(), 0.3, true, 48.0, -64.0, 46.5, -52.5),
        (t_very_short_arm(), 1.0, false, 0.0, -80.0, -1.5, -74.5),
        (t_arm_2d(), 1.0, false, 16.0, -80.0, 14.5, -73.5),
        (cross_short(), 2.0, false, 32.0, -80.0, 30.5, -72.5),
        (l_large_d(), 2.5, false, 48.0, -80.0, 46.5, -74.5),
        (t_arm_eq_d(), 1.0, false, 64.0, -80.0, 62.5, -74.5),
        (l_shape_45(), 0.3, false, 0.0, -96.0, -1.5, -92.0),
        (t_diagonal_branch(), 0.3, false, 16.0, -96.0, 14.5, -90.0),
        (y_mixed_junction(), 0.3, false, 32.0, -96.0, 30.5, -90.0),
        (room_with_corner_stub(), 0.3, true, 48.0, -96.0, 46.5, -85.0),
        (room_with_corner_diagonal(), 0.3, true, 68.0, -96.0, 66.5, -82.0),
        (room_with_near_corner_stub(), 0.3, true, 0.0, -112.0, -1.5, -100.5),
        (room_with_near_corner_diagonal(), 0.3, true, 20.0, -112.0, 18.5, -98.0),
    ];

    for (i, (pts, hw, closed, bx, by, lx, ly)) in cases.iter().enumerate() {
        register_label(
            storage,
            *lx,
            *ly,
            &format!("{}", i + 1),
            LABEL_SIZE,
            LABEL_COLOR,
        );
        draw_case(storage, pts, *hw, *closed, *bx, *by);
    }
}
