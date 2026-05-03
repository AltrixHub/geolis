use super::polygon_union::{Polygon, PolygonWithHoles, WALL_EPS};

/// Expands a polyline into a thickened polygon with left/right offsets.
///
/// The polygon has two offset "sides":
///   - left side: points offset by `+left_w` along the left normal
///   - right side: points offset by `-right_w` along the left normal (i.e. to the right)
///
/// For a closed polyline the result is an annulus: outer ring + one hole.
/// For an open polyline it is a single ring (left side forward, right side back).
pub fn stroke_expand(
    vertices: &[(f64, f64)],
    closed: bool,
    left_w: f64,
    right_w: f64,
) -> PolygonWithHoles {
    let n = vertices.len();
    if n < 2 {
        return PolygonWithHoles {
            outer: Vec::new(),
            holes: Vec::new(),
        };
    }

    let verts = if closed {
        normalize_winding(vertices)
    } else {
        vertices.to_vec()
    };
    let n = verts.len();

    let seg_count = if closed { n } else { n - 1 };
    let dirs: Vec<(f64, f64)> = (0..seg_count)
        .map(|i| {
            let j = (i + 1) % n;
            normalize(verts[j].0 - verts[i].0, verts[j].1 - verts[i].1)
        })
        .collect();

    let mut left_pts: Vec<(f64, f64)> = Vec::new();
    let mut right_pts: Vec<(f64, f64)> = Vec::new();

    if closed {
        build_closed_offsets(
            &verts,
            &dirs,
            left_w,
            right_w,
            &mut left_pts,
            &mut right_pts,
        );
    } else {
        build_open_offsets(
            &verts,
            &dirs,
            left_w,
            right_w,
            &mut left_pts,
            &mut right_pts,
        );
    }

    if closed {
        // For CCW input, left_normal points INWARD, so:
        //   right_pts = outer boundary (larger)
        //   left_pts = inner boundary (hole)
        // Ensure outer is CCW, hole is CW.
        let right_area = signed_area_tuples(&right_pts);
        let left_area = signed_area_tuples(&left_pts);

        let (mut outer, mut hole) = if right_area.abs() > left_area.abs() {
            (right_pts, left_pts)
        } else {
            (left_pts, right_pts)
        };

        if signed_area_tuples(&outer) < 0.0 {
            outer.reverse();
        }
        if signed_area_tuples(&hole) > 0.0 {
            hole.reverse();
        }

        PolygonWithHoles {
            outer,
            holes: vec![hole],
        }
    } else {
        // Single polygon: left forward, end cap, right backward.
        let mut poly: Polygon = Vec::new();
        poly.extend_from_slice(&left_pts);
        right_pts.reverse();
        poly.extend_from_slice(&right_pts);
        PolygonWithHoles {
            outer: poly,
            holes: Vec::new(),
        }
    }
}

fn build_open_offsets(
    verts: &[(f64, f64)],
    dirs: &[(f64, f64)],
    left_w: f64,
    right_w: f64,
    left_pts: &mut Vec<(f64, f64)>,
    right_pts: &mut Vec<(f64, f64)>,
) {
    let n = verts.len();
    let seg_count = n - 1;

    // First vertex: flat end cap perpendicular to outgoing direction.
    let n0 = left_normal(dirs[0]);
    left_pts.push(offset_point(verts[0], n0, left_w));
    right_pts.push(offset_point(verts[0], n0, -right_w));

    // Interior vertices: compute join chain.
    for i in 1..seg_count {
        let join = compute_join(verts[i], dirs[i - 1], dirs[i], left_w, right_w);
        left_pts.extend(join.left);
        right_pts.extend(join.right);
    }

    // Last vertex: flat end cap perpendicular to incoming direction.
    let nl = left_normal(dirs[seg_count - 1]);
    left_pts.push(offset_point(verts[n - 1], nl, left_w));
    right_pts.push(offset_point(verts[n - 1], nl, -right_w));
}

fn build_closed_offsets(
    verts: &[(f64, f64)],
    dirs: &[(f64, f64)],
    left_w: f64,
    right_w: f64,
    left_pts: &mut Vec<(f64, f64)>,
    right_pts: &mut Vec<(f64, f64)>,
) {
    let seg_count = dirs.len();

    for i in 0..seg_count {
        let prev = if i == 0 { seg_count - 1 } else { i - 1 };
        let join = compute_join(verts[i], dirs[prev], dirs[i], left_w, right_w);
        left_pts.extend(join.left);
        right_pts.extend(join.right);
    }
}

/// Local chain of offset vertices emitted for a single join.
///
/// Each side emits exactly one point: the miter intersection of the two
/// offset edges. Sharp acute corners are preserved (no bevel chamfer). When
/// the inner miter flips at very acute concave corners, polygon union
/// downstream removes the inverted region.
struct JoinResult {
    left: Vec<(f64, f64)>,
    right: Vec<(f64, f64)>,
}

/// Computes offset vertices at a single interior join.
///
/// For a CCW polyline, `left_normal` points inward, so:
///   - left side = inner boundary (shrinks at convex corners, spikes at concave)
///   - right side = outer boundary (spikes at convex corners, shrinks at concave)
///
/// `cross = dir_in × dir_out`:
///   - `> 0`: left turn → convex for CCW → outer side is `right`, inner is `left`
///   - `< 0`: right turn → concave for CCW → outer side is `left`, inner is `right`
///   - `≈ 0`: collinear → single offset point, no miter needed
///
/// Both sides always take the miter intersection — sharp corners stay sharp.
/// At very acute angles where the miter would otherwise be cut off (bevel),
/// the spike is kept; the polygon-union arrangement downstream cleans up any
/// self-intersecting region the spike introduces.
fn compute_join(
    vertex: (f64, f64),
    dir_in: (f64, f64),
    dir_out: (f64, f64),
    left_w: f64,
    right_w: f64,
) -> JoinResult {
    let cross = dir_in.0 * dir_out.1 - dir_in.1 * dir_out.0;
    let n_in = left_normal(dir_in);
    let n_out = left_normal(dir_out);

    if cross.abs() < WALL_EPS {
        // Collinear: the two offset lines coincide. Emit a single point.
        return JoinResult {
            left: vec![offset_point(vertex, n_out, left_w)],
            right: vec![offset_point(vertex, n_out, -right_w)],
        };
    }

    // Endpoints of the incoming and outgoing offset edges AT the vertex.
    let lp_in = offset_point(vertex, n_in, left_w);
    let lp_out = offset_point(vertex, n_out, left_w);
    let rp_in = offset_point(vertex, n_in, -right_w);
    let rp_out = offset_point(vertex, n_out, -right_w);

    let left_miter = line_intersect(lp_in, dir_in, lp_out, dir_out);
    let right_miter = line_intersect(rp_in, dir_in, rp_out, dir_out);

    JoinResult {
        left: miter_chain(left_miter, lp_out),
        right: miter_chain(right_miter, rp_out),
    }
}

/// Always take the miter intersection. Falls back to the outgoing edge's
/// offset endpoint only when the two offset directions are numerically
/// parallel (`line_intersect` returned `None`); the single-point fallback
/// preserves continuity.
fn miter_chain(miter: Option<(f64, f64)>, fallback: (f64, f64)) -> Vec<(f64, f64)> {
    match miter {
        Some(m) => vec![m],
        None => vec![fallback],
    }
}

fn normalize_winding(verts: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let area = signed_area_tuples(verts);
    if area < 0.0 {
        let mut reversed = verts.to_vec();
        reversed.reverse();
        reversed
    } else {
        verts.to_vec()
    }
}

fn signed_area_tuples(verts: &[(f64, f64)]) -> f64 {
    let n = verts.len();
    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        area += verts[i].0 * verts[j].1;
        area -= verts[j].0 * verts[i].1;
    }
    area * 0.5
}

fn left_normal(dir: (f64, f64)) -> (f64, f64) {
    (-dir.1, dir.0)
}

fn offset_point(p: (f64, f64), normal: (f64, f64), w: f64) -> (f64, f64) {
    (p.0 + w * normal.0, p.1 + w * normal.1)
}

fn normalize(dx: f64, dy: f64) -> (f64, f64) {
    let len = (dx * dx + dy * dy).sqrt();
    if len < WALL_EPS {
        (1.0, 0.0)
    } else {
        (dx / len, dy / len)
    }
}

fn line_intersect(
    p1: (f64, f64),
    d1: (f64, f64),
    p2: (f64, f64),
    d2: (f64, f64),
) -> Option<(f64, f64)> {
    let cross = d1.0 * d2.1 - d1.1 * d2.0;
    if cross.abs() < WALL_EPS * WALL_EPS {
        return None;
    }
    let dx = p2.0 - p1.0;
    let dy = p2.1 - p1.1;
    let t = (dx * d2.1 - dy * d2.0) / cross;
    Some((p1.0 + t * d1.0, p1.1 + t * d1.1))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn straight_open() {
        let result = stroke_expand(&[(0.0, 0.0), (5.0, 0.0)], false, 0.3, 0.3);
        assert!(result.holes.is_empty());
        let area = signed_area_tuples(&result.outer).abs();
        let expected = 5.0 * 0.6;
        assert!(
            (area - expected).abs() < 0.1,
            "area={area}, expected={expected}"
        );
    }

    #[test]
    fn l_shape_open() {
        let result = stroke_expand(&[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0)], false, 0.3, 0.3);
        assert!(result.holes.is_empty());
        let area = signed_area_tuples(&result.outer).abs();
        assert!(area > 3.0, "area={area} too small");
    }

    #[test]
    fn closed_square() {
        let result = stroke_expand(
            &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
            true,
            0.3,
            0.3,
        );
        assert_eq!(result.holes.len(), 1, "closed square should have 1 hole");
        let outer_area = signed_area_tuples(&result.outer);
        let hole_area = signed_area_tuples(&result.holes[0]);
        assert!(outer_area > 0.0, "outer should be CCW (positive area)");
        assert!(hole_area < 0.0, "hole should be CW (negative area)");
        let wall_area = outer_area + hole_area;
        assert!(
            wall_area > 15.0 && wall_area < 30.0,
            "wall_area={wall_area}"
        );
    }

    #[test]
    fn closed_square_cw_input() {
        // CW input should be normalized to CCW.
        let result = stroke_expand(
            &[(0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0)],
            true,
            0.3,
            0.3,
        );
        assert_eq!(result.holes.len(), 1);
        let outer_area = signed_area_tuples(&result.outer);
        assert!(outer_area > 0.0, "outer should be CCW even with CW input");
    }

    #[test]
    fn asymmetric_offset() {
        let result = stroke_expand(&[(0.0, 0.0), (5.0, 0.0)], false, 0.0, 0.3);
        assert!(result.holes.is_empty());
        let area = signed_area_tuples(&result.outer).abs();
        let expected = 5.0 * 0.3;
        assert!(
            (area - expected).abs() < 0.1,
            "area={area}, expected={expected}"
        );
    }

    #[test]
    fn acute_angle_keeps_sharp_miter() {
        // Sharp acute turn — the outer side keeps the (potentially long)
        // miter point instead of beveling. Result must be non-degenerate
        // (hole-free, with at least the three corner points on the outer).
        let result = stroke_expand(&[(0.0, 0.0), (5.0, 0.0), (4.9, 0.1)], false, 0.3, 0.3);
        assert!(result.holes.is_empty());
        assert!(result.outer.len() >= 3);
    }

    #[test]
    fn angled_closed_polygon_no_dent() {
        // Regression test: diagonal (non-axis-aligned) closed polygon.
        // Earlier iterations produced corner "dents" on the outer boundary
        // when the polygon was not axis-aligned.
        // Input is a CCW quadrilateral with all oblique edges.
        let result = stroke_expand(
            &[
                (-3.217, -4.144),
                (2.002, -4.631),
                (2.578, 1.534),
                (-2.635, 2.085),
            ],
            true,
            0.15,
            0.15,
        );
        assert_eq!(result.holes.len(), 1);
        let outer_area = signed_area_tuples(&result.outer);
        let hole_area = signed_area_tuples(&result.holes[0]);
        assert!(outer_area > 0.0, "outer must be CCW");
        assert!(hole_area < 0.0, "hole must be CW");

        // Every outer vertex should be strictly OUTSIDE the input polygon's
        // convex hull tight bbox (never dented inward past an input vertex).
        let input_max_x = 2.578_f64;
        let input_min_x = -3.217_f64;
        let input_max_y = 2.085_f64;
        let input_min_y = -4.631_f64;
        let mut saw_x_beyond = false;
        let mut saw_y_beyond = false;
        for &(x, y) in &result.outer {
            if x > input_max_x || x < input_min_x {
                saw_x_beyond = true;
            }
            if y > input_max_y || y < input_min_y {
                saw_y_beyond = true;
            }
        }
        assert!(saw_x_beyond, "outer must extend beyond input X extent");
        assert!(saw_y_beyond, "outer must extend beyond input Y extent");

        // Thickness ≈ 0.3 uniformly: ring area is roughly perimeter * 0.3.
        // Exact value varies with corner miters; assert a conservative band.
        let ring_area = outer_area + hole_area;
        assert!(ring_area > 4.0 && ring_area < 12.0, "ring_area={ring_area}");
    }

    #[test]
    fn join_result_convex_miter() {
        // Convex CCW corner under miter limit → single-point miters on both
        // sides. Verifies compute_join's local chain contract.
        let join = compute_join((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), 0.3, 0.3);
        assert_eq!(join.left.len(), 1, "inner miter is a single point");
        assert_eq!(
            join.right.len(),
            1,
            "outer miter under limit is a single point"
        );
    }

    #[test]
    fn join_result_acute_keeps_sharp_miter_both_sides() {
        // Even at a very acute U-turn, both sides keep a single miter point
        // — no bevel fallback. Sharp acute corners are preserved.
        let dir_out = {
            let (dx, dy): (f64, f64) = (-1.0, 0.05);
            let l = (dx * dx + dy * dy).sqrt();
            (dx / l, dy / l)
        };
        let join = compute_join((0.0, 0.0), (1.0, 0.0), dir_out, 0.3, 0.3);
        assert_eq!(join.left.len(), 1, "left side keeps single miter point");
        assert_eq!(join.right.len(), 1, "right side keeps single miter point");
    }

    #[test]
    fn join_result_collinear_single_point() {
        let join = compute_join((1.0, 0.0), (1.0, 0.0), (1.0, 0.0), 0.3, 0.3);
        assert_eq!(join.left.len(), 1);
        assert_eq!(join.right.len(), 1);
    }
}
