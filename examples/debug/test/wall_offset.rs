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

/// Case 23: Open L-shape at 45° angle.
fn l_shape_45() -> Vec<(f64, f64)> {
    vec![(0.0, 0.0), (5.0, 0.0), (8.0, 3.0)]
}

/// Case 24: Open T-shape with 45° upward-right branch.
fn t_diagonal_branch() -> Vec<(f64, f64)> {
    vec![(0.0, 3.0), (10.0, 3.0), (5.0, 3.0), (7.0, 5.0)]
}

/// Case 25: Open Y-junction: L-shape + diagonal from corner.
fn y_mixed_junction() -> Vec<(f64, f64)> {
    vec![(0.0, 0.0), (5.0, 0.0), (5.0, 5.0), (5.0, 0.0), (8.0, 3.0)]
}

/// Case 26: Closed room with diagonal stub from corner (0,0).
fn room_with_corner_stub() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0), (-3.0, -3.0), (0.0, 0.0),
        (8.0, 0.0), (8.0, 8.0), (0.0, 8.0),
    ]
}

/// Case 27: Closed room with diagonal partition through corner (0,0).
fn room_with_corner_diagonal() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0), (10.0, 0.0), (10.0, 8.0), (8.0, 8.0),
        (11.0, 11.0), (-3.0, -3.0),
        (0.0, 0.0), (0.0, 8.0), (8.0, 8.0),
    ]
}

/// Case 28: Closed room + diagonal stub near corner, junction at (0, 0.5).
fn room_with_near_corner_stub() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.5), (-3.0, -2.5), (0.0, 0.5),
        (0.0, 8.0), (8.0, 8.0), (8.0, 0.0), (0.0, 0.0),
    ]
}

/// Case 29: Closed room + diagonal partition near corner, junctions at (0, 0.5) and (7.5, 8).
fn room_with_near_corner_diagonal() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0), (10.0, 0.0), (10.0, 8.0), (7.5, 8.0),
        (10.5, 11.0), (-3.0, -2.5),
        (0.0, 0.5), (0.0, 8.0), (7.5, 8.0),
        (0.0, 0.5),
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

/// Case 18: T-shape, arm shorter than d → concave notch.
/// d=1.0, arm length 0.5 (half of d). Cap at y=3.5 below spine wall y=4.
fn case_18_outer() -> Vec<(f64, f64)> {
    vec![
        (0.0, 2.0), (4.0, 2.0), (8.0, 2.0), (8.0, 4.0),
        (5.0, 4.0), (5.0, 3.5), (3.0, 3.5), (3.0, 4.0), (0.0, 4.0),
    ]
}

/// Case 19: T-shape, arm = 2d → clean T outline.
/// d=1.0, arm length 2 (tip extends d beyond spine wall).
fn case_19_outer() -> Vec<(f64, f64)> {
    vec![
        (0.0, 2.0), (4.0, 2.0), (8.0, 2.0), (8.0, 4.0),
        (5.0, 4.0), (5.0, 5.0), (3.0, 5.0), (3.0, 4.0), (0.0, 4.0),
    ]
}

/// Case 20: Cross, arm length = d — degenerates to square.
/// d=2.0, arm length 2. All junction corners reach arm tips.
fn case_20_outer() -> Vec<(f64, f64)> {
    vec![(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0)]
}

/// Case 21: L-shape, d > horizontal leg → miter extends past original.
/// d=2.5, horizontal leg = 2. Miter at (2,0) pushes to (4.5,-2.5).
fn case_21_outer() -> Vec<(f64, f64)> {
    vec![
        (0.0, -2.5), (4.5, -2.5), (4.5, 4.0),
        (-0.5, 4.0), (-0.5, 2.5), (0.0, 2.5),
    ]
}

/// Case 22: T-shape, arm = d → exact degeneration boundary.
/// d=1.0, arm length 1. Arm side edges degenerate, cap on spine wall.
fn case_22_outer() -> Vec<(f64, f64)> {
    vec![
        (0.0, 2.0), (4.0, 2.0), (8.0, 2.0), (8.0, 4.0),
        (5.0, 4.0), (3.0, 4.0), (0.0, 4.0),
    ]
}

/// Helper constants for diagonal offset at d=0.3.
/// `ds` = d/sqrt(2), `dm` = d*(sqrt(2)-1).
fn diag_offsets(d: f64) -> (f64, f64) {
    let ds = d / std::f64::consts::SQRT_2;
    let dm = d * (std::f64::consts::SQRT_2 - 1.0);
    (ds, dm)
}

/// Case 23: Open L at 45°, d=0.3 → 6-vertex outline.
fn case_23_outer() -> Vec<(f64, f64)> {
    let d = 0.3;
    let (ds, dm) = diag_offsets(d);
    vec![
        (0.0, -d), (5.0 + dm, -d), (8.0 + ds, 3.0 - ds),
        (8.0 - ds, 3.0 + ds), (5.0 - dm, d), (0.0, d),
    ]
}

/// Case 24: Open T with 45° upward branch, d=0.3 → 9-vertex outline.
fn case_24_outer() -> Vec<(f64, f64)> {
    let d = 0.3;
    let (ds, _dm) = diag_offsets(d);
    let s2 = std::f64::consts::SQRT_2;
    vec![
        (0.0, 3.0 - d), (5.0, 3.0 - d), (10.0, 3.0 - d),
        (10.0, 3.0 + d), (5.0 + d * (1.0 + s2), 3.0 + d),
        (7.0 + ds, 5.0 - ds), (7.0 - ds, 5.0 + ds),
        (5.0 + d * (1.0 - s2), 3.0 + d), (0.0, 3.0 + d),
    ]
}

/// Case 25: Open Y-junction (2 orthogonal + 1 diagonal), d=0.3 → 9 vertices.
fn case_25_outer() -> Vec<(f64, f64)> {
    let d = 0.3;
    let (ds, dm) = diag_offsets(d);
    let s2 = std::f64::consts::SQRT_2;
    vec![
        (0.0, -d), (5.0 + dm, -d),
        (8.0 + ds, 3.0 - ds), (8.0 - ds, 3.0 + ds),
        (5.0 + d, d * (1.0 + s2)),
        (5.0 + d, 5.0), (5.0 - d, 5.0),
        (5.0 - d, d), (0.0, d),
    ]
}

/// Case 26: Closed room + diagonal stub from corner, d=0.3 → 7 outer.
fn case_26_outer() -> Vec<(f64, f64)> {
    let d = 0.3;
    let (ds, dm) = diag_offsets(d);
    vec![
        (8.0 + d, -d), (8.0 + d, 8.0 + d), (-d, 8.0 + d),
        (-d, dm),
        (-3.0 - ds, -3.0 + ds), (-3.0 + ds, -3.0 - ds),
        (dm, -d),
    ]
}

/// Case 26: inner room boundary → 4 vertices.
fn case_26_inner() -> Vec<(f64, f64)> {
    let d = 0.3;
    vec![(d, d), (8.0 - d, d), (8.0 - d, 8.0 - d), (d, 8.0 - d)]
}

/// Case 27: Closed room + diagonal partition through corner, d=0.3 → 11 outer.
fn case_27_outer() -> Vec<(f64, f64)> {
    let d = 0.3;
    let (ds, dm) = diag_offsets(d);
    let s2 = std::f64::consts::SQRT_2;
    vec![
        (10.0 + d, -d), (10.0 + d, 8.0 + d),
        (8.0 + d * (1.0 + s2), 8.0 + d),
        (11.0 + ds, 11.0 - ds), (11.0 - ds, 11.0 + ds),
        (8.0 + d * (1.0 - s2), 8.0 + d),
        (-d, 8.0 + d), (-d, dm),
        (-3.0 - ds, -3.0 + ds), (-3.0 + ds, -3.0 - ds),
        (dm, -d),
    ]
}

/// Case 27: inner bottom-right triangle → 4 vertices.
fn case_27_inner_br() -> Vec<(f64, f64)> {
    let d = 0.3;
    let s2 = std::f64::consts::SQRT_2;
    let dm = d * (s2 - 1.0);
    vec![
        (d * (1.0 + s2), d),
        (10.0 - d, d),
        (10.0 - d, 8.0 - d),
        (8.0 + dm, 8.0 - d),
    ]
}

/// Case 27: inner top-left triangle → 3 vertices.
fn case_27_inner_tl() -> Vec<(f64, f64)> {
    let d = 0.3;
    let s2 = std::f64::consts::SQRT_2;
    vec![
        (d, d * (1.0 + s2)),
        (d, 8.0 - d),
        (8.0 - d * (1.0 + s2), 8.0 - d),
    ]
}

/// Case 28: Closed room + near-corner stub, d=0.3 → 8 outer vertices.
/// Junction at (0, 0.5) on left wall, diagonal stub to (-3, -2.5).
fn case_28_outer() -> Vec<(f64, f64)> {
    let d = 0.3;
    let (ds, dm) = diag_offsets(d);
    let s2 = std::f64::consts::SQRT_2;
    vec![
        (8.0 + d, -d), (8.0 + d, 8.0 + d), (-d, 8.0 + d),
        (-d, 0.5 + dm),
        (-3.0 - ds, -2.5 + ds), (-3.0 + ds, -2.5 - ds),
        (-d, 0.5 - d * (1.0 + s2)),
        (-d, -d),
    ]
}

/// Case 28: inner room boundary → 5 vertices (includes collinear junction corner).
fn case_28_inner() -> Vec<(f64, f64)> {
    let d = 0.3;
    vec![(d, d), (d, 0.5), (d, 8.0 - d), (8.0 - d, 8.0 - d), (8.0 - d, d)]
}

/// Case 29: Closed room + near-corner diagonal, d=0.3 → 12 outer vertices.
/// Junctions at (0, 0.5) and (7.5, 8), diagonal y=x+0.5.
fn case_29_outer() -> Vec<(f64, f64)> {
    let d = 0.3;
    let (ds, dm) = diag_offsets(d);
    let s2 = std::f64::consts::SQRT_2;
    vec![
        (10.0 + d, -d), (10.0 + d, 8.0 + d),
        (7.5 + d * (1.0 + s2), 8.0 + d),
        (10.5 + ds, 11.0 - ds), (10.5 - ds, 11.0 + ds),
        (7.5 + d * (1.0 - s2), 8.0 + d),
        (-d, 8.0 + d), (-d, 0.5 + dm),
        (-3.0 - ds, -2.5 + ds), (-3.0 + ds, -2.5 - ds),
        (-d, 0.5 - d * (1.0 + s2)),
        (-d, -d),
    ]
}

/// Case 29: inner bottom-right region → 5 vertices.
fn case_29_inner_br() -> Vec<(f64, f64)> {
    let d = 0.3;
    let dm = d * (std::f64::consts::SQRT_2 - 1.0);
    vec![
        (d, d),
        (d, 0.5 - dm),
        (7.5 + dm, 8.0 - d),
        (10.0 - d, 8.0 - d),
        (10.0 - d, d),
    ]
}

/// Case 29: inner top-left region → 3 vertices.
fn case_29_inner_tl() -> Vec<(f64, f64)> {
    let d = 0.3;
    let s2 = std::f64::consts::SQRT_2;
    vec![
        (d, 0.5 + d * (1.0 + s2)),
        (d, 8.0 - d),
        (7.5 - d * (1.0 + s2), 8.0 - d),
    ]
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
    let c26_inner = case_26_inner();
    let c27_br = case_27_inner_br();
    let c27_tl = case_27_inner_tl();
    let c28_inner = case_28_inner();
    let c29_br = case_29_inner_br();
    let c29_tl = case_29_inner_tl();

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
            (t_very_short_arm(), case_18_outer(), vec![], false, 0.0, -80.0, -1.5, -74.5),
            (t_arm_2d(), case_19_outer(), vec![], false, 16.0, -80.0, 14.5, -73.5),
            (cross_short(), case_20_outer(), vec![], false, 32.0, -80.0, 30.5, -72.5),
            (l_large_d(), case_21_outer(), vec![], false, 48.0, -80.0, 46.5, -74.5),
            (t_arm_eq_d(), case_22_outer(), vec![], false, 64.0, -80.0, 62.5, -74.5),
            (l_shape_45(), case_23_outer(), vec![], false, 0.0, -96.0, -1.5, -92.0),
            (t_diagonal_branch(), case_24_outer(), vec![], false, 16.0, -96.0, 14.5, -90.0),
            (y_mixed_junction(), case_25_outer(), vec![], false, 32.0, -96.0, 30.5, -90.0),
            (room_with_corner_stub(), case_26_outer(), vec![&c26_inner], true, 48.0, -96.0, 46.5, -85.0),
            (room_with_corner_diagonal(), case_27_outer(), vec![&c27_br, &c27_tl], true, 68.0, -96.0, 66.5, -82.0),
            (room_with_near_corner_stub(), case_28_outer(), vec![&c28_inner], true, 0.0, -112.0, -1.5, -100.5),
            (room_with_near_corner_diagonal(), case_29_outer(), vec![&c29_br, &c29_tl], true, 20.0, -112.0, 18.5, -98.0),
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
