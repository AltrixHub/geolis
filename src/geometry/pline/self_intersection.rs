//! Self-intersection detection and resolution for closed plines.
//!
//! Used by [`WallOutline2D::execute`](crate::operations::offset::WallOutline2D::execute)
//! to guarantee its output boundaries are tessellation-safe (no
//! transversely-crossing non-adjacent edges, which would otherwise cause
//! `spade::cdt` to panic).
//!
//! **Scope: closed plines only.** Open plines are out of scope because
//! `WallOutline2D` always produces closed boundaries.

use super::{Pline, PlineVertex};
use crate::error::{OperationError, Result};

/// Tolerance for cross-product parallelism / degeneracy checks.
pub(crate) const CROSS_EPS: f64 = 1e-12;

/// Tolerance for segment parameter interior bounds.
pub(crate) const PARAM_EPS: f64 = 1e-9;

/// Work-queue safety bound for `split_at_self_intersections`.
pub(crate) const MAX_SPLIT_ITERATIONS: usize = 100;

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
// Recursive split at self-intersections
// ---------------------------------------------------------------------------

/// Recursively split a **closed** `pline` at every transverse self-intersection
/// until each output pline is simple. Returns one or more simple closed plines.
///
/// **Returning `Result` is intentional**, not stylistic: the caller
/// (`WallOutline2D::execute`) cannot honor its "simplicity guaranteed"
/// contract if the resolver returned a partial result containing
/// unresolved Plines. Bailing must be a hard failure that propagates up.
///
/// **Precondition:** `pline.closed == true`. The function asserts this
/// (`debug_assert!`) and returns `Ok(vec![pline])` unchanged in release
/// builds when called on an open pline (defensive no-op).
///
/// **Split contract** (formal). Given the first crossing at edges `i` and
/// `j` (`i < j`) with intersection point `P = (x, y)`, and the input
/// `pline.vertices = [v_0, ..., v_{n-1}]`:
///
/// - Loop A (length `j - i + 1`, closed): `[v_{i+1}, ..., v_j, P]`
/// - Loop B (length `n - j + i + 1`, closed):
///   `[v_{j+1}, ..., v_{n-1}, v_0, ..., v_i, P]`
///
/// (Sum of loop sizes is `n + 2`: every original vertex appears in exactly
/// one loop, P appears once per loop.) Both loops share P; they touch
/// geometrically at P but otherwise have disjoint interiors.
///
/// **Winding is NOT preserved.** Worked example: figure-8
/// `[(0,0), (2,2), (0,2), (2,0)]` with crossing at (1,1) splits into Loop A
/// `[(2,2), (0,2), (1,1)]` (signed area +1, CCW) and Loop B
/// `[(2,0), (0,0), (1,1)]` (signed area −1, CW). The original is a
/// degenerate figure-8 with signed area 0; both signs appear in the
/// children regardless of any "parent winding" notion. Callers needing
/// signed orientation must re-derive via shoelace area on the output.
///
/// **Termination guard.** Bails after `MAX_SPLIT_ITERATIONS = 100`
/// work-queue iterations and returns
/// `Err(OperationError::Failed("split_at_self_intersections: ..."))`.
/// The unresolved work queue is dropped — never returned, never silently
/// merged into output — so the simplicity guarantee is never violated.
///
/// **Output filter.** Loops with `< 3` vertices after split are dropped
/// (degenerate slivers). If, after termination, no simple loops were
/// produced (e.g. every split produced only degenerate children), returns
/// `Ok(vec![])`; the caller (`WallOutline2D::execute`) handles the empty
/// case via its existing `if outlines.is_empty()` check.
///
/// **Bulge handling.** Output `PlineVertex`s are constructed via
/// `PlineVertex::line(x, y)` (no bulge). Input vertices are read for x/y
/// only; any bulge on the input is dropped. For wall-outline boundaries
/// from `polygon_union::union_all_with_holes` (the only consumer in this
/// crate) bulges are always 0 on input, so this is a no-op for that path.
#[allow(
    dead_code,
    reason = "consumed by WallOutline2D::execute in plan-13k T3"
)]
#[allow(
    clippy::many_single_char_names,
    reason = "i/j (edge indices), x/y (intersection coords), n (vertex count), k (loop index), p (intersection vertex) match the formal split contract in this fn's doc-comment"
)]
pub(crate) fn split_at_self_intersections(pline: Pline) -> Result<Vec<Pline>> {
    debug_assert!(
        pline.closed,
        "split_at_self_intersections: open plines are out of scope"
    );
    if !pline.closed {
        return Ok(vec![pline]);
    }

    let mut work: Vec<Pline> = vec![pline];
    let mut output: Vec<Pline> = Vec::new();
    let mut iter: usize = 0;

    while let Some(current) = work.pop() {
        if iter >= MAX_SPLIT_ITERATIONS {
            return Err(OperationError::Failed(format!(
                "split_at_self_intersections: bailed at MAX_SPLIT_ITERATIONS={MAX_SPLIT_ITERATIONS}; \
                 input may be pathologically self-intersecting"
            ))
            .into());
        }
        iter += 1;

        match find_self_intersection(&current) {
            None => {
                if current.vertices.len() >= 3 {
                    output.push(current);
                }
            }
            Some((i, j, x, y)) => {
                let n = current.vertices.len();
                let p = PlineVertex::line(x, y);

                // Loop A: [v_{i+1}, ..., v_j, P]  (length j - i + 1)
                let mut loop_a: Vec<PlineVertex> = Vec::with_capacity(j - i + 1);
                for k in (i + 1)..=j {
                    loop_a.push(current.vertices[k]);
                }
                loop_a.push(p);

                // Loop B: [v_{j+1}, ..., v_{n-1}, v_0, ..., v_i, P]
                //         (length n - j + i + 1)
                let mut loop_b: Vec<PlineVertex> = Vec::with_capacity(n - j + i + 1);
                for k in (j + 1)..n {
                    loop_b.push(current.vertices[k]);
                }
                for k in 0..=i {
                    loop_b.push(current.vertices[k]);
                }
                loop_b.push(p);

                // Both loops are pushed back into the work queue; only loops with
                // < 3 vertices are dropped (degenerate slivers).
                if loop_a.len() >= 3 {
                    work.push(Pline {
                        vertices: loop_a,
                        closed: true,
                    });
                }
                if loop_b.len() >= 3 {
                    work.push(Pline {
                        vertices: loop_b,
                        closed: true,
                    });
                }
            }
        }
    }

    Ok(output)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
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

    // --- split_at_self_intersections ---

    /// Helper: assert all output plines are simple, closed, and have ≥3 vertices.
    fn assert_all_simple_and_well_formed(output: &[Pline]) {
        for (idx, p) in output.iter().enumerate() {
            assert!(
                p.closed,
                "output[{idx}] should be closed; got open pline with {} vertices",
                p.vertices.len()
            );
            assert!(
                p.vertices.len() >= 3,
                "output[{idx}] should have >=3 vertices; got {}",
                p.vertices.len()
            );
            assert!(
                find_self_intersection(p).is_none(),
                "output[{idx}] should be simple; found a self-intersection"
            );
        }
    }

    #[test]
    fn split_simple_closed_pline_returns_input_unchanged() {
        // Unit square. Already simple → output is one pline of length 4.
        let p = closed_pline_xy(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]);
        let original_xy: Vec<(f64, f64)> = p.vertices.iter().map(|v| (v.x, v.y)).collect();
        let out = split_at_self_intersections(p).expect("should not bail");
        assert_eq!(out.len(), 1, "simple pline should not be split");
        assert!(out[0].closed);
        let result_xy: Vec<(f64, f64)> = out[0].vertices.iter().map(|v| (v.x, v.y)).collect();
        assert_eq!(result_xy, original_xy, "vertex order should be preserved");
    }

    #[test]
    fn split_figure_eight_yields_two_simple_loops_touching_at_crossing() {
        // Figure-8 with crossing at (1, 1).
        let p = closed_pline_xy(&[(0.0, 0.0), (2.0, 2.0), (0.0, 2.0), (2.0, 0.0)]);
        let out = split_at_self_intersections(p).expect("should not bail");

        assert_eq!(out.len(), 2, "figure-8 should split into 2 loops");
        assert_all_simple_and_well_formed(&out);

        // Both output plines must contain the crossing point (1, 1).
        let count_with_crossing = out
            .iter()
            .filter(|p| {
                p.vertices
                    .iter()
                    .any(|v| (v.x - 1.0).abs() < 1e-9 && (v.y - 1.0).abs() < 1e-9)
            })
            .count();
        assert_eq!(
            count_with_crossing, 2,
            "both loops should share the crossing vertex (1, 1)"
        );
    }

    #[test]
    fn split_double_zigzag_yields_only_simple_outputs() {
        // 6-vertex closed pline with two non-overlapping crossings.
        // Layout: alternating-up zigzag designed so edges 0–2 cross AND
        // edges 1–4 cross (when read with closed wrap).
        //
        //   v0 = (0, 0)
        //   v1 = (4, 4)
        //   v2 = (1, 4)   ← edge 1 (v1→v2) is the top horizontal-ish
        //   v3 = (4, 0)   ← edge 2 (v2→v3) crosses edge 0 (v0→v1) at (2, 2)
        //   v4 = (3, 4)
        //   v5 = (0, 2)   ← edge 5 (v5→v0) and edge 4 (v4→v5) close the loop
        let p = closed_pline_xy(&[
            (0.0, 0.0),
            (4.0, 4.0),
            (1.0, 4.0),
            (4.0, 0.0),
            (3.0, 4.0),
            (0.0, 2.0),
        ]);
        let out = split_at_self_intersections(p).expect("should not bail");

        assert!(
            !out.is_empty(),
            "double-crossing pline should yield at least one simple loop"
        );
        assert!(
            out.len() <= 4,
            "double-crossing pline should yield at most 4 loops; got {}",
            out.len()
        );
        assert_all_simple_and_well_formed(&out);
    }
}
