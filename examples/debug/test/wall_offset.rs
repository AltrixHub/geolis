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

// ── Drawing helper ──────────────────────────────────────────────────

fn draw_ground_truth(
    storage: &MeshStorage,
    centerline: &[(f64, f64)],
    outer: &[(f64, f64)],
    holes: &[&[(f64, f64)]],
    bx: f64,
    by: f64,
) {
    // Centerline in gray.
    let cl: Vec<Point3> = centerline
        .iter()
        .map(|&(x, y)| Point3::new(x + bx, y + by, 0.0))
        .collect();
    if let Ok(s) = StrokeStyle::new(0.05) {
        register_stroke(storage, &cl, s, false, GRAY);
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

    // (centerline, outer, holes, bx, by, lx, ly)
    let cases: Vec<(Vec<(f64, f64)>, Vec<(f64, f64)>, Vec<&[(f64, f64)]>, f64, f64, f64, f64)> =
        vec![
            (single_line(), case_1_outer(), vec![], 0.0, 0.0, -1.5, 1.5),
            (l_shape(), case_2_outer(), vec![], 10.0, 0.0, 8.5, 6.0),
            (t_shape(), case_3_outer(), vec![], 22.0, 0.0, 20.5, 4.5),
            (cross_shape(), case_4_outer(), vec![], 36.0, 0.0, 34.5, 11.5),
            (double_cross(), case_5_outer(), vec![&c5_hole], 0.0, -16.0, -1.5, -4.5),
            (double_cross(), case_6_outer(), vec![&c6_hole], 16.0, -16.0, 14.5, -4.5),
            (y_fork(), case_7_outer(), vec![], 32.0, -16.0, 30.5, -4.5),
            (h_shape(), case_8_outer(), vec![], 0.0, -32.0, -1.5, -20.5),
            (e_shape(), case_9_outer(), vec![], 16.0, -32.0, 14.5, -22.5),
            (angled_cross(), case_10_outer(), vec![], 32.0, -32.0, 30.5, -20.5),
        ];

    for (i, (cl, outer, holes, bx, by, lx, ly)) in cases.iter().enumerate() {
        register_label(
            storage,
            *lx,
            *ly,
            &format!("{}", i + 1),
            LABEL_SIZE,
            LABEL_COLOR,
        );
        draw_ground_truth(storage, cl, outer, holes, *bx, *by);
    }
}
