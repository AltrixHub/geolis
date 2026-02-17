//! Ground truth for `wall_offset` — hand-computed expected outlines.

use geolis::math::Point3;
use geolis::tessellation::StrokeStyle;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_label, register_stroke};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const GRAY: Color = Color::rgb(160, 160, 160);
const GREEN: Color = Color::rgb(100, 220, 100);
const BLUE: Color = Color::rgb(100, 150, 255);

// ── Centerline definitions (same as patterns/wall_offset.rs) ────────

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

fn room_with_corridor() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0), (5.0, 0.0), (5.0, -5.0), (5.0, 0.0),
        (10.0, 0.0), (10.0, 10.0), (0.0, 10.0),
    ]
}

fn room_with_partition() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0), (10.0, 0.0), (10.0, 5.0), (0.0, 5.0),
        (10.0, 5.0), (10.0, 10.0), (0.0, 10.0),
    ]
}

fn room_with_diagonal_wall() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0), (10.0, 0.0), (10.0, 7.5),
        (15.0, 10.0), (-5.0, 0.0),
        (0.0, 2.5), (0.0, 10.0), (10.0, 10.0),
        (10.0, 7.5), (0.0, 2.5),
    ]
}

fn room_with_penetrating_wall() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0), (10.0, 0.0), (10.0, 5.0), (13.0, 5.0),
        (10.0, 5.0), (0.0, 5.0), (10.0, 5.0),
        (10.0, 10.0), (0.0, 10.0), (0.0, 5.0),
        (-3.0, 5.0), (0.0, 5.0),
    ]
}

// ── Expected outlines (hand-computed) ───────────────────────────────

/// Case 1: Single line d=0.3 → 4-vertex rectangle.
fn case_1_outer() -> Vec<(f64, f64)> {
    vec![(0.0, -0.3), (5.0, -0.3), (5.0, 0.3), (0.0, 0.3)]
}

/// Case 2: L-shape d=0.3 → 6-vertex L outline.
fn case_2_outer() -> Vec<(f64, f64)> {
    vec![
        (0.0, -0.3), (5.3, -0.3), (5.3, 5.0),
        (4.7, 5.0), (4.7, 0.3), (0.0, 0.3),
    ]
}

/// Case 3: T-shape d=0.3 → 9-vertex T outline.
fn case_3_outer() -> Vec<(f64, f64)> {
    vec![
        (4.7, 0.0), (5.3, 0.0), (5.3, 2.7), (10.0, 2.7), (10.0, 3.3),
        (5.0, 3.3), (0.0, 3.3), (0.0, 2.7), (4.7, 2.7),
    ]
}

/// Case 4: Cross d=0.3 → 12-vertex cross outline.
fn case_4_outer() -> Vec<(f64, f64)> {
    vec![
        (4.7, 0.0), (5.3, 0.0), (5.3, 4.7), (10.0, 4.7), (10.0, 5.3),
        (5.3, 5.3), (5.3, 10.0), (4.7, 10.0), (4.7, 5.3), (0.0, 5.3),
        (0.0, 4.7), (4.7, 4.7),
    ]
}

/// Case 5: 井-shape d=0.3 → 28-vertex outer boundary.
fn case_5_outer() -> Vec<(f64, f64)> {
    vec![
        (2.7, 0.0), (3.3, 0.0), (3.3, 2.7), (6.7, 2.7), (6.7, 0.0),
        (7.3, 0.0), (7.3, 2.7), (10.0, 2.7), (10.0, 3.3), (7.3, 3.3),
        (7.3, 6.7), (10.0, 6.7), (10.0, 7.3), (7.3, 7.3), (7.3, 10.0),
        (6.7, 10.0), (6.7, 7.3), (3.3, 7.3), (3.3, 10.0), (2.7, 10.0),
        (2.7, 7.3), (0.0, 7.3), (0.0, 6.7), (2.7, 6.7), (2.7, 3.3),
        (0.0, 3.3), (0.0, 2.7), (2.7, 2.7),
    ]
}

/// Case 5: 井-shape d=0.3 → 4-vertex inner hole (center room).
fn case_5_hole() -> Vec<(f64, f64)> {
    vec![(3.3, 3.3), (3.3, 6.7), (6.7, 6.7), (6.7, 3.3)]
}

/// Case 6: 井-shape d=0.8 → 28-vertex outer boundary.
fn case_6_outer() -> Vec<(f64, f64)> {
    vec![
        (2.2, 0.0), (3.8, 0.0), (3.8, 2.2), (6.2, 2.2), (6.2, 0.0),
        (7.8, 0.0), (7.8, 2.2), (10.0, 2.2), (10.0, 3.8), (7.8, 3.8),
        (7.8, 6.2), (10.0, 6.2), (10.0, 7.8), (7.8, 7.8), (7.8, 10.0),
        (6.2, 10.0), (6.2, 7.8), (3.8, 7.8), (3.8, 10.0), (2.2, 10.0),
        (2.2, 7.8), (0.0, 7.8), (0.0, 6.2), (2.2, 6.2), (2.2, 3.8),
        (0.0, 3.8), (0.0, 2.2), (2.2, 2.2),
    ]
}

/// Case 6: 井-shape d=0.8 → 4-vertex inner hole.
fn case_6_hole() -> Vec<(f64, f64)> {
    vec![(3.8, 3.8), (3.8, 6.2), (6.2, 6.2), (6.2, 3.8)]
}

/// Case 7: Y-fork d=0.5 → 9-vertex outline.
/// Center (5,5), arms at 120° intervals, arm length 5.
fn case_7_outer() -> Vec<(f64, f64)> {
    // Derived from: d=0.5, sin60=√3/2, cos60=1/2
    // Junction corners: (4.5, 5+√3/6), (5, 5−√3/3), (5.5, 5+√3/6)
    vec![
        (0.9199, 2.0670),
        (5.0, 4.4226),
        (9.0801, 2.0670),
        (9.5801, 2.9330),
        (5.5, 5.2887),
        (5.5, 10.0),
        (4.5, 10.0),
        (4.5, 5.2887),
        (0.4199, 2.9330),
    ]
}

/// Case 8: H-shape d=0.3 → 12-vertex outer boundary.
/// Two vertical columns (x=3, x=7) connected by horizontal bar at y=5.
fn case_8_outer() -> Vec<(f64, f64)> {
    vec![
        (2.7, 0.0), (3.3, 0.0), (3.3, 4.7), (6.7, 4.7), (6.7, 0.0),
        (7.3, 0.0), (7.3, 10.0), (6.7, 10.0), (6.7, 5.3), (3.3, 5.3),
        (3.3, 10.0), (2.7, 10.0),
    ]
}

/// Case 9: E-shape d=0.3 → 16-vertex outer boundary.
/// Horizontal spine y=3 with 3 vertical prongs at x=2, x=6, x=10.
fn case_9_outer() -> Vec<(f64, f64)> {
    vec![
        (0.0, 2.7), (12.0, 2.7), (12.0, 3.3), (10.3, 3.3),
        (10.3, 8.0), (9.7, 8.0), (9.7, 3.3), (6.3, 3.3),
        (6.3, 8.0), (5.7, 8.0), (5.7, 3.3), (2.3, 3.3),
        (2.3, 8.0), (1.7, 8.0), (1.7, 3.3), (0.0, 3.3),
    ]
}

/// Case 10: Angled cross 60° d=0.3 → 12-vertex outer boundary.
/// Horizontal line y=5 crossed by 60° diagonal through (5,5).
fn case_10_outer() -> Vec<(f64, f64)> {
    vec![
        (2.7598, 0.5199), (5.1732, 4.7), (10.0, 4.7), (10.0, 5.3),
        (5.5196, 5.3), (7.7598, 9.1801), (7.2402, 9.4801), (4.8268, 5.3),
        (0.0, 5.3), (0.0, 4.7), (4.4804, 4.7), (2.2402, 0.8199),
    ]
}

/// Case 11: Square room d=0.3 → outer boundary (CW, miter corners).
fn case_11_outer() -> Vec<(f64, f64)> {
    vec![(-0.3, -0.3), (10.3, -0.3), (10.3, 10.3), (-0.3, 10.3)]
}

/// Case 11: Square room d=0.3 → inner boundary.
fn case_11_inner() -> Vec<(f64, f64)> {
    vec![(0.3, 0.3), (9.7, 0.3), (9.7, 9.7), (0.3, 9.7)]
}

/// Case 12: Rectangle room (12×8) d=0.3 → outer boundary.
fn case_12_outer() -> Vec<(f64, f64)> {
    vec![(-0.3, -0.3), (12.3, -0.3), (12.3, 8.3), (-0.3, 8.3)]
}

/// Case 12: Rectangle room d=0.3 → inner boundary.
fn case_12_inner() -> Vec<(f64, f64)> {
    vec![(0.3, 0.3), (11.7, 0.3), (11.7, 7.7), (0.3, 7.7)]
}

/// Case 13: L-room d=0.3 → outer boundary.
fn case_13_outer() -> Vec<(f64, f64)> {
    vec![
        (-0.3, -0.3), (5.3, -0.3), (5.3, 3.3),
        (3.3, 3.3), (3.3, 5.3), (-0.3, 5.3),
    ]
}

/// Case 13: L-room d=0.3 → inner boundary.
fn case_13_inner() -> Vec<(f64, f64)> {
    vec![
        (0.3, 0.3), (4.7, 0.3), (4.7, 2.7),
        (2.7, 2.7), (2.7, 4.7), (0.3, 4.7),
    ]
}

/// Case 14: Room + corridor d=0.3 → 8-vertex outer.
fn case_14_outer() -> Vec<(f64, f64)> {
    vec![
        (4.7, -5.0), (5.3, -5.0), (5.3, -0.3), (10.3, -0.3),
        (10.3, 10.3), (-0.3, 10.3), (-0.3, -0.3), (4.7, -0.3),
    ]
}

/// Case 14: Room + corridor d=0.3 → 4-vertex inner room.
fn case_14_inner() -> Vec<(f64, f64)> {
    vec![(0.3, 0.3), (0.3, 9.7), (9.7, 9.7), (9.7, 0.3)]
}

/// Case 15: Room + partition d=0.3 → 4-vertex outer.
fn case_15_outer() -> Vec<(f64, f64)> {
    vec![(-0.3, -0.3), (10.3, -0.3), (10.3, 10.3), (-0.3, 10.3)]
}

/// Case 15: Room + partition d=0.3 → bottom room.
fn case_15_inner_bottom() -> Vec<(f64, f64)> {
    vec![(0.3, 0.3), (0.3, 4.7), (9.7, 4.7), (9.7, 0.3)]
}

/// Case 15: Room + partition d=0.3 → top room.
fn case_15_inner_top() -> Vec<(f64, f64)> {
    vec![(0.3, 5.3), (0.3, 9.7), (9.7, 9.7), (9.7, 5.3)]
}

/// Case 17: Room + diagonal wall d=0.3 → 12-vertex outer.
/// Diagonal from (-5,0) to (15,10), slope 1/2, direction (2,1)/sqrt(5).
/// At vertical wall x=X, diagonal offset y = 2.5 + X/2 ± d*sqrt(5)/2.
fn case_17_outer() -> Vec<(f64, f64)> {
    let dn = 0.3 / 5.0_f64.sqrt();
    let ds = 0.3 * 5.0_f64.sqrt() / 2.0; // d*sqrt(5)/2, y-shift at vertical wall
    vec![
        (-0.3, -0.3), (10.3, -0.3),
        (10.3, 2.5 + 10.3 / 2.0 - ds),  // right wall meets diagonal right side
        (15.0 + dn, 10.0 - 2.0 * dn),   // dead end cap right-bottom
        (15.0 - dn, 10.0 + 2.0 * dn),   // dead end cap right-top
        (10.3, 2.5 + 10.3 / 2.0 + ds),  // right wall meets diagonal left side
        (10.3, 10.3), (-0.3, 10.3),
        (-0.3, 2.5 - 0.3 / 2.0 + ds),   // left wall meets diagonal left side
        (-5.0 - dn, 2.0 * dn),           // dead end cap left-top
        (-5.0 + dn, -2.0 * dn),          // dead end cap left-bottom
        (-0.3, 2.5 - 0.3 / 2.0 - ds),   // left wall meets diagonal right side
    ]
}

/// Case 17: Room + diagonal wall d=0.3 → bottom-right inner room.
fn case_17_inner_bottom() -> Vec<(f64, f64)> {
    let ds = 0.3 * 5.0_f64.sqrt() / 2.0;
    vec![
        (0.3, 0.3),
        (0.3, 2.5 + 0.3 / 2.0 - ds),
        (9.7, 2.5 + 9.7 / 2.0 - ds),
        (9.7, 0.3),
    ]
}

/// Case 17: Room + diagonal wall d=0.3 → top-left inner room.
fn case_17_inner_top() -> Vec<(f64, f64)> {
    let ds = 0.3 * 5.0_f64.sqrt() / 2.0;
    vec![
        (0.3, 2.5 + 0.3 / 2.0 + ds),
        (0.3, 9.7),
        (9.7, 9.7),
        (9.7, 2.5 + 9.7 / 2.0 + ds),
    ]
}

/// Case 16: Room + penetrating wall d=0.3 → 12-vertex outer.
fn case_16_outer() -> Vec<(f64, f64)> {
    vec![
        (-0.3, -0.3), (10.3, -0.3), (10.3, 4.7), (13.0, 4.7),
        (13.0, 5.3), (10.3, 5.3), (10.3, 10.3), (-0.3, 10.3),
        (-0.3, 5.3), (-3.0, 5.3), (-3.0, 4.7), (-0.3, 4.7),
    ]
}

/// Case 16: Room + penetrating wall d=0.3 → bottom inner room.
fn case_16_inner_bottom() -> Vec<(f64, f64)> {
    vec![(0.3, 0.3), (0.3, 4.7), (9.7, 4.7), (9.7, 0.3)]
}

/// Case 16: Room + penetrating wall d=0.3 → top inner room.
fn case_16_inner_top() -> Vec<(f64, f64)> {
    vec![(0.3, 5.3), (0.3, 9.7), (9.7, 9.7), (9.7, 5.3)]
}

// ── Drawing helper ──────────────────────────────────────────────────

fn draw_ground_truth(
    storage: &MeshStorage,
    centerline: &[(f64, f64)],
    outer: &[(f64, f64)],
    holes: &[&[(f64, f64)]],
    closed: bool,
    bx: f64,
    by: f64,
) {
    // Centerline in gray.
    let cl: Vec<Point3> = centerline
        .iter()
        .map(|&(x, y)| Point3::new(x + bx, y + by, 0.0))
        .collect();
    if let Ok(s) = StrokeStyle::new(0.05) {
        register_stroke(storage, &cl, s, closed, GRAY);
    }

    // Expected outer boundary in green.
    let op: Vec<Point3> = outer
        .iter()
        .map(|&(x, y)| Point3::new(x + bx, y + by, 0.0))
        .collect();
    if let Ok(s) = StrokeStyle::new(0.08) {
        register_stroke(storage, &op, s, true, GREEN);
    }

    // Expected holes in blue.
    for hole in holes {
        let hp: Vec<Point3> = hole
            .iter()
            .map(|&(x, y)| Point3::new(x + bx, y + by, 0.0))
            .collect();
        if let Ok(s) = StrokeStyle::new(0.08) {
            register_stroke(storage, &hp, s, true, BLUE);
        }
    }
}

// ── Registration ────────────────────────────────────────────────────

/// Register `wall_offset` ground truth meshes.
#[allow(clippy::too_many_lines, clippy::type_complexity)]
pub fn register(storage: &MeshStorage) {
    let c5_hole = case_5_hole();
    let c6_hole = case_6_hole();
    let c11_inner = case_11_inner();
    let c12_inner = case_12_inner();
    let c13_inner = case_13_inner();
    let c14_inner = case_14_inner();
    let c15_bottom = case_15_inner_bottom();
    let c15_top = case_15_inner_top();
    let c16_bottom = case_16_inner_bottom();
    let c16_top = case_16_inner_top();
    let c17_bottom = case_17_inner_bottom();
    let c17_top = case_17_inner_top();

    // (centerline, outer, holes, closed, bx, by, lx, ly)
    #[allow(clippy::type_complexity)]
    let cases: Vec<(Vec<(f64, f64)>, Vec<(f64, f64)>, Vec<&[(f64, f64)]>, bool, f64, f64, f64, f64)> =
        vec![
            (single_line(), case_1_outer(), vec![], false, 0.0, 0.0, -1.5, 1.5),
            (l_shape(), case_2_outer(), vec![], false, 10.0, 0.0, 8.5, 6.0),
            (t_shape(), case_3_outer(), vec![], false, 22.0, 0.0, 20.5, 4.5),
            (cross_shape(), case_4_outer(), vec![], false, 36.0, 0.0, 34.5, 11.5),
            (double_cross(), case_5_outer(), vec![&c5_hole], false, 0.0, -16.0, -1.5, -4.5),
            (double_cross(), case_6_outer(), vec![&c6_hole], false, 16.0, -16.0, 14.5, -4.5),
            (y_fork(), case_7_outer(), vec![], false, 32.0, -16.0, 30.5, -4.5),
            (h_shape(), case_8_outer(), vec![], false, 0.0, -32.0, -1.5, -20.5),
            (e_shape(), case_9_outer(), vec![], false, 16.0, -32.0, 14.5, -22.5),
            (angled_cross(), case_10_outer(), vec![], false, 32.0, -32.0, 30.5, -20.5),
            (square_room(), case_11_outer(), vec![&c11_inner], true, 0.0, -48.0, -1.5, -36.5),
            (rectangle_room(), case_12_outer(), vec![&c12_inner], true, 16.0, -48.0, 14.5, -36.5),
            (l_room(), case_13_outer(), vec![&c13_inner], true, 32.0, -48.0, 30.5, -41.5),
            (room_with_corridor(), case_14_outer(), vec![&c14_inner], true, 0.0, -64.0, -1.5, -52.5),
            (room_with_partition(), case_15_outer(), vec![&c15_bottom, &c15_top], true, 16.0, -64.0, 14.5, -52.5),
            (room_with_penetrating_wall(), case_16_outer(), vec![&c16_bottom, &c16_top], true, 32.0, -64.0, 30.5, -52.5),
            (room_with_diagonal_wall(), case_17_outer(), vec![&c17_bottom, &c17_top], true, 48.0, -64.0, 46.5, -52.5),
        ];

    for (i, (cl, outer, holes, closed, bx, by, lx, ly)) in cases.iter().enumerate() {
        register_label(
            storage,
            *lx,
            *ly,
            &format!("{}", i + 1),
            LABEL_SIZE,
            LABEL_COLOR,
        );
        draw_ground_truth(storage, cl, outer, holes, *closed, *bx, *by);
    }
}
