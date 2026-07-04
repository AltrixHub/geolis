//! Self-intersection detection primitives for closed plines.
//!
//! Provides [`segment_segment_intersection_2d`] and
//! [`find_self_intersection`], used by the `WallOutline2D` test oracle
//! (P3.1 S2) and by per-fixture self-intersection assertions.
//!
//! **Scope: closed plines only.** Open plines are out of scope because
//! `WallOutline2D` always produces closed boundaries.

use super::Pline;

/// Tolerance for cross-product parallelism / degeneracy checks.
pub(crate) const CROSS_EPS: f64 = 1e-12;

/// Tolerance for segment parameter interior bounds.
pub(crate) const PARAM_EPS: f64 = 1e-9;

// ---------------------------------------------------------------------------
// Segment-segment intersection (2D)
// ---------------------------------------------------------------------------

/// Segment-segment intersection in 2D. Returns parameters `(t, u)` where
/// `t` is the position along the first segment AB and `u` along the
/// second segment CD, both strictly interior:
/// `PARAM_EPS < t, u < 1 - PARAM_EPS`.
///
/// Returns `None` for:
/// - Parallel or collinear segments (cross product magnitude below `CROSS_EPS`).
/// - Segments that touch only at an endpoint (parameter outside the
///   interior range on either segment).
/// - Tangent grazes that fall within either ε.
///
/// **This intentionally detects only transverse crossings.** Collinear
/// overlap and endpoint-touching are not flagged because they do not
/// trigger `spade::cdt` panics — those are the failure mode this module
/// is designed to prevent.
#[allow(
    clippy::many_single_char_names,
    reason = "a/b/c/d are the standard segment-endpoint convention in computational geometry"
)]
pub(crate) fn segment_segment_intersection_2d(
    a: (f64, f64),
    b: (f64, f64),
    c: (f64, f64),
    d: (f64, f64),
) -> Option<(f64, f64)> {
    let d1x = b.0 - a.0;
    let d1y = b.1 - a.1;
    let d2x = d.0 - c.0;
    let d2y = d.1 - c.1;

    let cross = d1x * d2y - d1y * d2x;
    if cross.abs() < CROSS_EPS {
        return None; // parallel or degenerate
    }

    let d3x = c.0 - a.0;
    let d3y = c.1 - a.1;

    let t = (d3x * d2y - d3y * d2x) / cross;
    let u = (d3x * d1y - d3y * d1x) / cross;

    if t > PARAM_EPS && t < 1.0 - PARAM_EPS && u > PARAM_EPS && u < 1.0 - PARAM_EPS {
        Some((t, u))
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Pline self-intersection scan
// ---------------------------------------------------------------------------

/// Find the first pair of non-adjacent edges in a **closed** `pline.vertices`
/// that cross at strictly-interior parameters.
///
/// Edges are indexed `0..n` where edge `k` connects `vertices[k]` to
/// `vertices[(k + 1) % n]` (the wrap-around closes the pline).
/// "Non-adjacent" means `i != j`, `(i + 1) % n != j`, and `(j + 1) % n != i`.
///
/// Returns `Some((i, j, x, y))` with `i < j`. The crossing point `(x, y)`
/// is the 2D location at parameter `t` on edge `i`.
/// Returns `None` when no such transverse crossing exists.
///
/// **Definition of "simple" used here.** This function detects only
/// transverse crossings (see [`segment_segment_intersection_2d`] for the
/// precise condition). Collinear-overlapping edges, endpoint-touching
/// edges, and tangent grazes are NOT flagged because they do not trigger
/// `spade::cdt` panics.
///
/// **Precondition:** `pline.closed == true`. Behavior on `closed == false`
/// is `None` (open input is out of scope; see module-level doc).
///
/// **Bulge handling:** This function reads only `vertex.x, vertex.y` and
/// treats every edge as a straight segment. Bulged (arc) edges are
/// approximated by their endpoint-to-endpoint chord here. For wall
/// outline boundaries returned by `polygon_union::union_all_with_holes`
/// — the only caller in this crate — bulges are always 0, so the
/// approximation is exact for that pipeline.
#[allow(
    clippy::many_single_char_names,
    reason = "a/b/c/d (segment endpoints), i/j (edge indices), n (vertex count), x/y (intersection coords) are domain-standard names in 2D segment intersection geometry"
)]
pub(crate) fn find_self_intersection(pline: &Pline) -> Option<(usize, usize, f64, f64)> {
    if !pline.closed {
        return None;
    }
    let n = pline.vertices.len();
    if n < 4 {
        return None; // < 4 vertices cannot have non-adjacent edge pairs in a closed pline
    }

    for i in 0..n {
        let a = (pline.vertices[i].x, pline.vertices[i].y);
        let b = (pline.vertices[(i + 1) % n].x, pline.vertices[(i + 1) % n].y);
        for j in (i + 2)..n {
            // Skip the wrap-around adjacency: when i = 0, j = n-1 the edges share v_0.
            if i == 0 && j == n - 1 {
                continue;
            }
            let c = (pline.vertices[j].x, pline.vertices[j].y);
            let d = (pline.vertices[(j + 1) % n].x, pline.vertices[(j + 1) % n].y);
            if let Some((t, _u)) = segment_segment_intersection_2d(a, b, c, d) {
                let x = a.0 + t * (b.0 - a.0);
                let y = a.1 + t * (b.1 - a.1);
                return Some((i, j, x, y));
            }
        }
    }
    None
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
#[allow(clippy::expect_used, clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::pline::PlineVertex;

    fn closed_pline_xy(pts: &[(f64, f64)]) -> Pline {
        Pline {
            vertices: pts.iter().map(|&(x, y)| PlineVertex::line(x, y)).collect(),
            closed: true,
        }
    }

    // --- segment_segment_intersection_2d ---

    #[test]
    fn segment_segment_intersection_crossing_returns_params() {
        // X-shape: AB = (0,0)→(2,2), CD = (0,2)→(2,0). Cross at (1,1) at t=u=0.5.
        let result =
            segment_segment_intersection_2d((0.0, 0.0), (2.0, 2.0), (0.0, 2.0), (2.0, 0.0));
        let (t, u) = result.expect("should detect transverse crossing");
        assert!((t - 0.5).abs() < 1e-9, "t expected 0.5, got {t}");
        assert!((u - 0.5).abs() < 1e-9, "u expected 0.5, got {u}");
    }

    #[test]
    fn segment_segment_intersection_parallel_returns_none() {
        // Two horizontal segments at y=0 and y=1.
        assert!(
            segment_segment_intersection_2d((0.0, 0.0), (2.0, 0.0), (0.0, 1.0), (2.0, 1.0))
                .is_none()
        );
    }

    #[test]
    fn segment_segment_intersection_endpoint_touch_returns_none() {
        // T-junction: AB = (0,0)→(2,0), CD = (2,0)→(2,2). Touch at AB endpoint
        // (t = 1.0 exactly), which falls outside (PARAM_EPS, 1 - PARAM_EPS).
        assert!(
            segment_segment_intersection_2d((0.0, 0.0), (2.0, 0.0), (2.0, 0.0), (2.0, 2.0))
                .is_none()
        );
    }

    // --- find_self_intersection ---

    #[test]
    fn find_self_intersection_simple_quad_returns_none() {
        // Unit square — simple closed pline, no self-intersection.
        let p = closed_pline_xy(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]);
        assert!(find_self_intersection(&p).is_none());
    }

    #[test]
    // Intersection coordinates read clearest in i/j/x/y notation.
    #[allow(clippy::many_single_char_names)]
    fn find_self_intersection_figure_eight_returns_crossing() {
        // Figure-8: vertices 0..3 = (0,0), (2,2), (0,2), (2,0).
        // Edges:
        //   0: (0,0)→(2,2)
        //   1: (2,2)→(0,2)
        //   2: (0,2)→(2,0)  ← crosses edge 0 at (1,1)
        //   3: (2,0)→(0,0)
        // Non-adjacent pair (0, 2) crosses at (1,1).
        let p = closed_pline_xy(&[(0.0, 0.0), (2.0, 2.0), (0.0, 2.0), (2.0, 0.0)]);
        let (i, j, x, y) = find_self_intersection(&p).expect("figure-8 should self-intersect");
        assert_eq!(i, 0);
        assert_eq!(j, 2);
        assert!((x - 1.0).abs() < 1e-9, "x expected 1.0, got {x}");
        assert!((y - 1.0).abs() < 1e-9, "y expected 1.0, got {y}");
    }

    #[test]
    fn find_self_intersection_open_pline_returns_none() {
        // Defensive: open pline is out of scope; function returns None.
        let p = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(2.0, 2.0),
                PlineVertex::line(0.0, 2.0),
                PlineVertex::line(2.0, 0.0),
            ],
            closed: false,
        };
        assert!(find_self_intersection(&p).is_none());
    }
}
