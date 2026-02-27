//! Wall self-intersection debug visualization.
//!
//! Shows wall outlines that self-intersect due to sharp miter angles,
//! along with the centrelines and the result after self-intersection
//! resolution (distance-based miter vertex filtering).

use geolis::geometry::pline::{Pline, PlineVertex};
use geolis::math::distance_2d::point_to_segment_dist;
use geolis::math::Point3;
use geolis::operations::offset::WallOutline2D;
use geolis::tessellation::StrokeStyle;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_label, register_stroke};

const LABEL_SIZE: f64 = 1.0;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const GRAY: Color = Color::rgb(120, 120, 120);
const GREEN: Color = Color::rgb(100, 220, 100);
const RED: Color = Color::rgb(255, 80, 80);

// ── Self-intersection resolution (same logic as wall_group.rs) ──────

fn min_distance_to_centerlines(px: f64, py: f64, centerlines: &[Pline]) -> f64 {
    let mut min_d = f64::MAX;
    for cl in centerlines {
        let n = cl.vertices.len();
        let seg_count = if cl.closed { n } else { n.saturating_sub(1) };
        for i in 0..seg_count {
            let a = &cl.vertices[i];
            let b = &cl.vertices[(i + 1) % n];
            let d = point_to_segment_dist(px, py, a.x, a.y, b.x, b.y);
            if d < min_d {
                min_d = d;
            }
        }
    }
    min_d
}

fn find_clip_point(
    inside_x: f64, inside_y: f64,
    outside_x: f64, outside_y: f64,
    centerlines: &[Pline],
    threshold: f64,
) -> (f64, f64) {
    let mut lo = 0.0_f64;
    let mut hi = 1.0_f64;
    for _ in 0..50 {
        let mid = (lo + hi) * 0.5;
        let mx = inside_x + mid * (outside_x - inside_x);
        let my = inside_y + mid * (outside_y - inside_y);
        let d = min_distance_to_centerlines(mx, my, centerlines);
        if d <= threshold {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let t = (lo + hi) * 0.5;
    (
        inside_x + t * (outside_x - inside_x),
        inside_y + t * (outside_y - inside_y),
    )
}

fn clip_miter_vertices(
    verts: &mut Vec<PlineVertex>,
    centerlines: &[Pline],
    max_allowed: f64,
) {
    let n = verts.len();
    if n < 3 {
        return;
    }
    let mut result: Vec<PlineVertex> = Vec::with_capacity(n);
    for i in 0..n {
        let d = min_distance_to_centerlines(verts[i].x, verts[i].y, centerlines);
        if d > max_allowed {
            let prev = if i > 0 { i - 1 } else { n - 1 };
            let next = (i + 1) % n;
            let (cx1, cy1) = find_clip_point(
                verts[prev].x, verts[prev].y,
                verts[i].x, verts[i].y,
                centerlines, max_allowed,
            );
            let (cx2, cy2) = find_clip_point(
                verts[next].x, verts[next].y,
                verts[i].x, verts[i].y,
                centerlines, max_allowed,
            );
            result.push(PlineVertex::line(cx1, cy1));
            result.push(PlineVertex::line(cx2, cy2));
        } else {
            result.push(PlineVertex::line(verts[i].x, verts[i].y));
        }
    }
    *verts = result;
}

fn segment_segment_intersection(
    ax: f64, ay: f64, bx: f64, by: f64,
    cx: f64, cy: f64, dx: f64, dy: f64,
) -> Option<(f64, f64)> {
    let d1x = bx - ax;
    let d1y = by - ay;
    let d2x = dx - cx;
    let d2y = dy - cy;
    let cross = d1x * d2y - d1y * d2x;
    if cross.abs() < 1e-12 {
        return None;
    }
    let d3x = cx - ax;
    let d3y = cy - ay;
    let t = (d3x * d2y - d3y * d2x) / cross;
    let u = (d3x * d1y - d3y * d1x) / cross;
    let eps = 1e-9;
    if t > eps && t < 1.0 - eps && u > eps && u < 1.0 - eps {
        Some((t, u))
    } else {
        None
    }
}

fn find_self_intersection(vertices: &[PlineVertex]) -> Option<(usize, usize, f64, f64)> {
    let n = vertices.len();
    for i in 0..n {
        let i_next = (i + 1) % n;
        let (ax, ay) = (vertices[i].x, vertices[i].y);
        let (bx, by) = (vertices[i_next].x, vertices[i_next].y);
        for j in (i + 2)..n {
            let j_next = (j + 1) % n;
            if j_next == i {
                continue;
            }
            let (cx, cy) = (vertices[j].x, vertices[j].y);
            let (dx, dy) = (vertices[j_next].x, vertices[j_next].y);
            if let Some((t, _)) = segment_segment_intersection(ax, ay, bx, by, cx, cy, dx, dy) {
                return Some((i, j, ax + t * (bx - ax), ay + t * (by - ay)));
            }
        }
    }
    None
}

fn polygon_area_pv(vertices: &[PlineVertex]) -> f64 {
    let n = vertices.len();
    if n < 3 {
        return 0.0;
    }
    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        area += vertices[i].x * vertices[j].y;
        area -= vertices[j].x * vertices[i].y;
    }
    area / 2.0
}

/// Resolve to a single outer-contour polygon per boundary (ground truth).
/// Miter clip → split at self-intersection → keep largest area loop.
fn resolve_outer_contour(
    boundaries: &mut Vec<Pline>,
    centerlines: &[Pline],
    half_thickness: f64,
) {
    let max_allowed = half_thickness * 4.0;
    let mut result: Vec<Pline> = Vec::new();

    for boundary in boundaries.drain(..) {
        if !boundary.closed || boundary.vertices.len() < 4 {
            result.push(boundary);
            continue;
        }

        let mut verts = boundary.vertices.clone();
        clip_miter_vertices(&mut verts, centerlines, max_allowed);
        if verts.len() < 3 {
            continue;
        }

        for _iter in 0..100 {
            if verts.len() < 4 {
                break;
            }
            match find_self_intersection(&verts) {
                Some((edge_i, edge_j, ix, iy)) => {
                    let n = verts.len();
                    let mut loop_a = vec![PlineVertex::line(ix, iy)];
                    for k in (edge_i + 1)..=edge_j {
                        loop_a.push(PlineVertex::line(verts[k].x, verts[k].y));
                    }
                    let mut loop_b = vec![PlineVertex::line(ix, iy)];
                    for k in (edge_j + 1)..n {
                        loop_b.push(PlineVertex::line(verts[k].x, verts[k].y));
                    }
                    for k in 0..=edge_i {
                        loop_b.push(PlineVertex::line(verts[k].x, verts[k].y));
                    }
                    let area_a = polygon_area_pv(&loop_a).abs();
                    let area_b = polygon_area_pv(&loop_b).abs();
                    verts = if area_a >= area_b { loop_a } else { loop_b };
                }
                None => break,
            }
        }

        if verts.len() >= 3 {
            result.push(Pline { vertices: verts, closed: true });
        }
    }

    *boundaries = result;
}

/// Resolve keeping both loops (algorithm output for wall_group.rs).
/// Miter clip → split at self-intersection → keep both loops.
fn resolve_keep_both(
    boundaries: &mut Vec<Pline>,
    centerlines: &[Pline],
    half_thickness: f64,
) {
    let max_allowed = half_thickness * 4.0;
    let mut result: Vec<Pline> = Vec::new();

    for boundary in boundaries.drain(..) {
        if !boundary.closed || boundary.vertices.len() < 4 {
            result.push(boundary);
            continue;
        }

        let mut verts = boundary.vertices.clone();
        clip_miter_vertices(&mut verts, centerlines, max_allowed);
        if verts.len() < 3 {
            continue;
        }

        let mut pending = vec![verts];

        for _iter in 0..100 {
            let mut next_pending: Vec<Vec<PlineVertex>> = Vec::new();
            let mut made_progress = false;

            for verts in pending.drain(..) {
                if verts.len() < 4 {
                    if verts.len() >= 3 {
                        result.push(Pline { vertices: verts, closed: true });
                    }
                    continue;
                }

                match find_self_intersection(&verts) {
                    Some((edge_i, edge_j, ix, iy)) => {
                        let n = verts.len();
                        let mut loop_a = vec![PlineVertex::line(ix, iy)];
                        for k in (edge_i + 1)..=edge_j {
                            loop_a.push(PlineVertex::line(verts[k].x, verts[k].y));
                        }
                        let mut loop_b = vec![PlineVertex::line(ix, iy)];
                        for k in (edge_j + 1)..n {
                            loop_b.push(PlineVertex::line(verts[k].x, verts[k].y));
                        }
                        for k in 0..=edge_i {
                            loop_b.push(PlineVertex::line(verts[k].x, verts[k].y));
                        }

                        if loop_a.len() >= 3 {
                            next_pending.push(loop_a);
                        }
                        if loop_b.len() >= 3 {
                            next_pending.push(loop_b);
                        }
                        made_progress = true;
                    }
                    None => {
                        result.push(Pline { vertices: verts, closed: true });
                    }
                }
            }

            pending = next_pending;
            if !made_progress || pending.is_empty() {
                break;
            }
        }

        for verts in pending {
            if verts.len() >= 3 {
                result.push(Pline { vertices: verts, closed: true });
            }
        }
    }

    *boundaries = result;
}

// ── Test case definitions ───────────────────────────────────────────

/// 3-point V: last point's wall barely overlaps the 1st segment.
/// Centerline: (0,0)→(3,5)→(0.2,1) — 3rd point is ~0.15 from segment 1.
fn v_slight_overlap() -> Vec<(f64, f64)> {
    vec![(0.0, 0.0), (3.0, 5.0), (0.2, 1.0)]
}

/// 3-point V: last point's wall clearly overlaps the 1st segment.
/// Centerline: (0,0)→(3,5)→(-0.1,1) — 3rd point crosses segment 1.
fn v_clear_overlap() -> Vec<(f64, f64)> {
    vec![(0.0, 0.0), (3.0, 5.0), (-0.1, 1.0)]
}

/// 3-point V: same shape but thicker wall.
fn v_thick_overlap() -> Vec<(f64, f64)> {
    vec![(0.0, 0.0), (3.0, 5.0), (0.2, 1.0)]
}

/// Near-reversal (3 points, almost 180° turn).
fn near_reversal() -> Vec<(f64, f64)> {
    vec![(0.0, 0.0), (5.0, 0.0), (0.0, 0.5)]
}

/// 4-point zigzag: 4th point's wall slightly crosses 1st segment's wall.
fn zigzag_4_overlap() -> Vec<(f64, f64)> {
    vec![
        (0.0, 0.0),
        (3.0, 4.0),
        (0.0, 3.0),
        (0.1, 0.5),
    ]
}

// ── Drawing helper ──────────────────────────────────────────────────

fn draw_case(
    storage: &MeshStorage,
    pts: &[(f64, f64)],
    half_w: f64,
    bx: f64,
    by: f64,
    label: &str,
) {
    let thin = StrokeStyle::new(0.03).unwrap_or_else(|_| unreachable!());
    let medium = StrokeStyle::new(0.06).unwrap_or_else(|_| unreachable!());

    let pline = Pline {
        vertices: pts.iter().map(|&(x, y)| PlineVertex::line(x, y)).collect(),
        closed: false,
    };
    let centerlines = vec![pline.clone()];

    let raw_boundaries = WallOutline2D::new(vec![pline], half_w)
        .execute()
        .unwrap_or_default();

    // LEFT: RED = raw WallOutline2D output (BEFORE any processing)
    // Shows self-intersecting boundaries — the problem we're fixing.

    // Centreline (gray, thin)
    let center: Vec<Point3> = pts
        .iter()
        .map(|&(x, y)| Point3::new(x + bx, y + by, 0.0))
        .collect();
    register_stroke(storage, &center, thin, false, GRAY);

    for ol in &raw_boundaries {
        let p: Vec<Point3> = ol
            .vertices
            .iter()
            .map(|v| Point3::new(v.x + bx, v.y + by, 0.0))
            .collect();
        register_stroke(storage, &p, medium, ol.closed, RED);
    }

    // RIGHT: GREEN = algorithm output (keep both loops)
    let mut resolved = raw_boundaries.clone();
    resolve_keep_both(&mut resolved, &centerlines, half_w);

    let rx = bx + 15.0;
    let center2: Vec<Point3> = pts
        .iter()
        .map(|&(x, y)| Point3::new(x + rx, y + by, 0.0))
        .collect();
    register_stroke(storage, &center2, thin, false, GRAY);

    for ol in &resolved {
        let p: Vec<Point3> = ol
            .vertices
            .iter()
            .map(|v| Point3::new(v.x + rx, v.y + by, 0.0))
            .collect();
        register_stroke(storage, &p, medium, ol.closed, GREEN);
    }

    // Labels
    register_label(storage, bx - 1.0, by - 2.5, label, LABEL_SIZE, LABEL_COLOR);
    register_label(storage, rx - 1.0, by - 2.5, label, LABEL_SIZE, LABEL_COLOR);
}

// ── Registration ────────────────────────────────────────────────────

/// Register `wall_self_intersect` pattern meshes.
///
/// Left column: RED = raw WallOutline2D output (BEFORE, may self-intersect)
/// Right column: GREEN = resolved output (AFTER, miter clip + split + keep both)
pub fn register(storage: &MeshStorage) {
    // Labels for columns
    register_label(storage, 0.0, 5.0, "0", LABEL_SIZE * 1.5, RED);      // BEFORE
    register_label(storage, 15.0, 5.0, "0", LABEL_SIZE * 1.5, GREEN);   // AFTER

    let cases: Vec<(&str, Vec<(f64, f64)>, f64)> = vec![
        ("1", v_slight_overlap(), 0.15),       // 3pt V, slight overlap, thin
        ("2", v_clear_overlap(), 0.15),        // 3pt V, clear overlap, thin
        ("3", v_thick_overlap(), 0.3),         // 3pt V, thicker wall
        ("4", near_reversal(), 0.3),           // 3pt near-reversal
        ("5", zigzag_4_overlap(), 0.15),       // 4pt, last overlaps 1st seg
    ];

    for (i, (label, pts, hw)) in cases.iter().enumerate() {
        let by = -(i as f64) * 12.0;
        draw_case(storage, pts, *hw, 0.0, by, label);
    }
}
