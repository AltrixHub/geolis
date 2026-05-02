//! Tessellation-safety normalization for `WallOutline2D` boundary sets.
//!
//! `polygon_union::union_all_with_holes` is a "candidate boundary producer" —
//! it strokes input polylines into thick bands, unions the bands, carves
//! holes, and traces the resulting graph into closed polygons. It does NOT
//! guarantee its output is safe for direct ingestion by spade-CDT-based
//! tessellation. Specifically, the output set may contain:
//!
//! - Near-duplicate consecutive vertices (within `WALL_EPS`).
//! - Zero-length / sub-tolerance edges produced by trace artifacts.
//! - Collinear vertex chains where `polygon_union` inserted split points.
//! - Open traces (start ≠ end) that are forced into a `Pline { closed: true }`
//!   by the caller, producing an implicit closing edge that may transversely
//!   cross other edges in the set.
//! - **Inter-boundary transverse crossings** where a constraint edge in one
//!   pline crosses a constraint edge in another (outer vs hole, two outers,
//!   or nested layers). `TessellateFace` feeds outer + hole loops into the
//!   SAME CDT, so this is a CDT-rejected configuration even when each
//!   boundary is individually simple.
//!
//! This module enforces the final CDT-safe contract via
//! [`make_tessellation_safe`], which:
//!
//! 1. Per boundary: collapse near-duplicate vertices, drop zero-length
//!    edges (T7).
//! 2. Per boundary: simplify collinear chains (T7).
//! 3. Per boundary: validate / repair closure (T7).
//! 4. **Across the entire boundary set**: detect transverse crossings
//!    (intra- and inter-boundary) and split at every crossing point (T8).
//! 5. Optional CDT dry-run: insert all vertices + constraint edges into a
//!    fresh `spade::ConstrainedDelaunayTriangulation` and `Err` if any
//!    insertion is rejected (T9).
//!
//! The intra-boundary primitives [`find_self_intersection`](super::self_intersection::find_self_intersection)
//! and [`split_at_self_intersections`](super::self_intersection::split_at_self_intersections)
//! remain useful within `make_tessellation_safe` as the same-boundary
//! case of step 4. They are NOT removed; they are reused.

use super::self_intersection::{segment_segment_intersection_2d, split_at_self_intersections};
use super::{Pline, PlineVertex};
use crate::error::{OperationError, Result};
use crate::operations::offset::wall_outline::polygon_union::WALL_EPS;

/// Maximum work-loop iterations before bailing. Mirrors
/// `MAX_SPLIT_ITERATIONS` in `self_intersection` but at the boundary-set
/// scope. The set-scope bound is 10x the per-boundary bound: a real wall
/// network can produce many crossings legitimately, so we want headroom
/// while still defending against pathological inputs.
pub(crate) const MAX_SET_SPLIT_ITERATIONS: usize = 1000;

/// Final-contract enforcement pass for `WallOutline2D::execute` output.
///
/// Transforms a set of polygon-union candidate boundaries into a set that
/// is safe for direct ingestion by `spade::cdt`-based tessellation. See
/// the module-level doc-comment for the full sequence of normalization
/// steps.
///
/// **Returns** `Ok(Vec<Pline>)` where every boundary is closed and the
/// boundary set as a whole has zero transverse constraint-edge crossings.
/// **Returns `Err`** when:
///
/// - Boundary cleanup fails (e.g. degenerate trace that cannot be repaired).
/// - The cross-boundary split work-loop bails at its safety bound.
/// - The CDT dry-run rejects any insertion (defense-in-depth: catches
///   spade-reject cases not covered by transverse-crossing detection).
///
/// **Empty input** returns `Ok(vec![])` unchanged. Empty output (every
/// boundary collapsed to a degenerate sliver) also returns `Ok(vec![])` —
/// the caller (`WallOutline2D::execute`) handles empty via its existing
/// `if outlines.is_empty()` check.
#[allow(
    dead_code,
    reason = "Wired into WallOutline2D::execute in plan-13k T10."
)]
pub(crate) fn make_tessellation_safe(boundaries: Vec<Pline>) -> Result<Vec<Pline>> {
    // T7 — per-boundary cleanup. Each phase mutates the boundary in place
    // (or replaces it). Boundaries that collapse to < 3 vertices are
    // dropped silently — `WallOutline2D::execute`'s existing
    // `if outlines.is_empty()` guard handles the all-dropped case as Err.
    let mut out: Vec<Pline> = Vec::with_capacity(boundaries.len());
    for mut b in boundaries {
        collapse_near_duplicate_vertices(&mut b, WALL_EPS);
        simplify_collinear_chain(&mut b, WALL_EPS);
        if !validate_and_close(&mut b, WALL_EPS) {
            // Degenerate trace that closure repair could not fix — drop.
            continue;
        }
        if b.vertices.len() >= 3 {
            out.push(b);
        }
    }

    // T8 — cross-boundary transverse crossing split. Iterate until the
    // entire set has no remaining transverse crossings (intra- or inter-
    // boundary). Bails after MAX_SET_SPLIT_ITERATIONS to defend against
    // pathological inputs.
    let mut iter = 0usize;
    loop {
        if iter >= MAX_SET_SPLIT_ITERATIONS {
            return Err(OperationError::Failed(format!(
                "tessellation_safety: cross-boundary split bailed at \
                 MAX_SET_SPLIT_ITERATIONS={MAX_SET_SPLIT_ITERATIONS}; \
                 input may be pathologically self-intersecting"
            ))
            .into());
        }
        iter += 1;
        let Some(crossing) = find_first_set_crossing(&out) else {
            break;
        };
        out = split_at_crossing(out, &crossing)?;
    }

    // T9 — CDT dry-run verification (defense-in-depth). Catches
    // spade-reject configurations that our transverse-crossing detection
    // missed (e.g. floating-point noise putting a near-crossing just
    // outside the PARAM_EPS interior bound).
    verify_cdt_safe(&out)?;

    Ok(out)
}

// ---------------------------------------------------------------------------
// T7: per-boundary cleanup
// ---------------------------------------------------------------------------

/// Collapse runs of consecutive vertices that lie within `eps` of each other
/// in 2D, including the wrap-around for closed plines. Drops zero-length
/// edges as a side effect (a zero-length edge is exactly a duplicate vertex
/// pair).
///
/// **In-place.** May reduce `pline.vertices.len()` to below 3 — caller must
/// re-check after calling.
pub(crate) fn collapse_near_duplicate_vertices(pline: &mut Pline, eps: f64) {
    let eps_sq = eps * eps;
    let close = |a: &PlineVertex, b: &PlineVertex| {
        let dx = a.x - b.x;
        let dy = a.y - b.y;
        dx * dx + dy * dy < eps_sq
    };

    if pline.vertices.is_empty() {
        return;
    }

    // Sequential collapse: keep first, then keep each next vertex only if
    // it is not within eps of the previous kept vertex.
    let mut kept: Vec<PlineVertex> = Vec::with_capacity(pline.vertices.len());
    kept.push(pline.vertices[0]);
    for v in pline.vertices.iter().skip(1) {
        // SAFETY: kept is non-empty (we just pushed vertices[0]).
        let Some(prev) = kept.last() else { continue };
        if !close(prev, v) {
            kept.push(*v);
        }
    }

    // Wrap-around dedup for closed plines: drop the trailing vertex if it
    // collapses to the first.
    if pline.closed && kept.len() >= 2 {
        let (Some(last), Some(first)) = (kept.last(), kept.first()) else {
            pline.vertices = kept;
            return;
        };
        if close(last, first) {
            kept.pop();
        }
    }

    pline.vertices = kept;
}

/// Drop interior vertices that lie collinearly between their neighbours
/// (cross product magnitude below `eps`). For closed plines, the wrap-around
/// neighbours are checked too.
///
/// **In-place.** May reduce `pline.vertices.len()` to below 3 — caller must
/// re-check after calling. Iterates until stable to handle cascading
/// collinear chains (`A → B → C → D` where all four are collinear collapses
/// to `A → D`).
#[allow(
    clippy::many_single_char_names,
    reason = "i/k/n are standard polygon-loop indexing in computational geometry"
)]
pub(crate) fn simplify_collinear_chain(pline: &mut Pline, eps: f64) {
    if pline.vertices.len() < 3 {
        return;
    }
    loop {
        let n = pline.vertices.len();
        if n < 3 {
            return;
        }
        let mut next: Vec<PlineVertex> = Vec::with_capacity(n);
        let mut removed = false;
        for i in 0..n {
            // For closed plines, every vertex has wrap-around neighbours.
            // For open plines, the first and last vertices are kept
            // unconditionally (they are endpoints, not interior to a chain).
            if !pline.closed && (i == 0 || i == n - 1) {
                next.push(pline.vertices[i]);
                continue;
            }
            let prev = pline.vertices[(i + n - 1) % n];
            let curr = pline.vertices[i];
            let succ = pline.vertices[(i + 1) % n];
            let cross =
                (curr.x - prev.x) * (succ.y - curr.y) - (curr.y - prev.y) * (succ.x - curr.x);
            if cross.abs() > eps {
                next.push(curr);
            } else {
                removed = true;
            }
        }
        if !removed || next.len() == pline.vertices.len() {
            pline.vertices = next;
            return;
        }
        pline.vertices = next;
    }
}

/// For a closed pline, validate that the trace is consistent (start and
/// end are within `eps`). `polygon_union::trace` does NOT check closure
/// — it accepts any `len >= 3` traced sequence — so the caller may have
/// forced `closed: true` on an open trace, producing an implicit closing
/// edge from `vertices.last()` back to `vertices.first()` that may not
/// reflect the geometry.
///
/// Behavior:
/// - For open plines (`!closed`): no-op, returns `true`.
/// - For closed plines with `len < 3`: returns `false` (degenerate).
/// - For closed plines where `first ≈ last` within `eps`: pop the
///   trailing duplicate and return `true`. (The implicit wrap edge is
///   then identical to the explicit last edge — consistent.)
/// - For closed plines where `first ≢ last`: returns `true`. The implicit
///   wrap edge is part of the boundary by design (caller marked it
///   `closed: true`); we do NOT insert a second copy of the start
///   point as that would create a zero-length tail edge.
///
/// Returns `false` only for the degenerate-length case; in all other
/// cases the caller can assume the pline's vertex sequence is a valid
/// closed loop.
pub(crate) fn validate_and_close(pline: &mut Pline, eps: f64) -> bool {
    if !pline.closed {
        return true;
    }
    if pline.vertices.len() < 3 {
        return false;
    }
    let eps_sq = eps * eps;
    let first = pline.vertices[0];
    let last = pline.vertices[pline.vertices.len() - 1];
    let dx = first.x - last.x;
    let dy = first.y - last.y;
    if dx * dx + dy * dy < eps_sq {
        pline.vertices.pop();
    }
    pline.vertices.len() >= 3
}

// ---------------------------------------------------------------------------
// T8: cross-boundary transverse crossing split
// ---------------------------------------------------------------------------

/// A transverse crossing between two edges in a boundary set, identified
/// by the boundary index, the edge index within that boundary, and the
/// 2D crossing point.
///
/// Invariant: `boundary_a <= boundary_b`. When `boundary_a == boundary_b`
/// (intra-boundary same-loop crossing), `edge_a < edge_b`.
#[derive(Debug, Clone, Copy, PartialEq)]
struct SetCrossing {
    boundary_a: usize,
    edge_a: usize,
    boundary_b: usize,
    edge_b: usize,
    x: f64,
    y: f64,
}

/// Find the first transverse crossing between any pair of edges in the
/// boundary set. Returns `None` when the set is CDT-safe with respect
/// to transverse crossings.
///
/// Scans:
/// - All non-adjacent edge pairs WITHIN each boundary (intra-boundary;
///   reuses the same conditions as
///   [`super::self_intersection::find_self_intersection`]).
/// - All edge pairs ACROSS distinct boundaries (inter-boundary).
///
/// Determinism: outer loop over `boundary_a` ascending; inner loop over
/// `boundary_b >= boundary_a` ascending; then edge indices ascending.
/// The first crossing found in this scan order is returned.
#[allow(
    clippy::many_single_char_names,
    reason = "i/j/n/k are domain-standard polygon-loop indexing"
)]
#[allow(
    clippy::needless_range_loop,
    reason = "ba/bb are paired indices into a common slice; constructing a SetCrossing requires both indices, which iterator+enumerate would obscure"
)]
fn find_first_set_crossing(boundaries: &[Pline]) -> Option<SetCrossing> {
    for ba in 0..boundaries.len() {
        let pa = &boundaries[ba];
        if !pa.closed {
            continue;
        }
        let na = pa.vertices.len();
        if na < 3 {
            continue;
        }

        for bb in ba..boundaries.len() {
            let pb = &boundaries[bb];
            if !pb.closed {
                continue;
            }
            let nb = pb.vertices.len();
            if nb < 3 {
                continue;
            }

            if ba == bb {
                // Intra-boundary scan: non-adjacent pairs in the same
                // closed loop (mirrors find_self_intersection).
                if na < 4 {
                    continue;
                }
                for i in 0..na {
                    let a = (pa.vertices[i].x, pa.vertices[i].y);
                    let b = (pa.vertices[(i + 1) % na].x, pa.vertices[(i + 1) % na].y);
                    for j in (i + 2)..na {
                        if i == 0 && j == na - 1 {
                            continue; // wrap-around adjacency
                        }
                        let c = (pa.vertices[j].x, pa.vertices[j].y);
                        let d = (pa.vertices[(j + 1) % na].x, pa.vertices[(j + 1) % na].y);
                        if let Some((t, _u)) = segment_segment_intersection_2d(a, b, c, d) {
                            let x = a.0 + t * (b.0 - a.0);
                            let y = a.1 + t * (b.1 - a.1);
                            return Some(SetCrossing {
                                boundary_a: ba,
                                edge_a: i,
                                boundary_b: bb,
                                edge_b: j,
                                x,
                                y,
                            });
                        }
                    }
                }
            } else {
                // Inter-boundary scan: every edge of pa vs every edge of pb.
                // No adjacency to skip — distinct boundaries share no
                // vertices by index.
                for i in 0..na {
                    let a = (pa.vertices[i].x, pa.vertices[i].y);
                    let b = (pa.vertices[(i + 1) % na].x, pa.vertices[(i + 1) % na].y);
                    for j in 0..nb {
                        let c = (pb.vertices[j].x, pb.vertices[j].y);
                        let d = (pb.vertices[(j + 1) % nb].x, pb.vertices[(j + 1) % nb].y);
                        if let Some((t, _u)) = segment_segment_intersection_2d(a, b, c, d) {
                            let x = a.0 + t * (b.0 - a.0);
                            let y = a.1 + t * (b.1 - a.1);
                            return Some(SetCrossing {
                                boundary_a: ba,
                                edge_a: i,
                                boundary_b: bb,
                                edge_b: j,
                                x,
                                y,
                            });
                        }
                    }
                }
            }
        }
    }
    None
}

/// Split the affected boundaries at the crossing point.
///
/// - **Same-boundary case** (`boundary_a == boundary_b`): replace that
///   single boundary with the two simple loops produced by
///   [`super::self_intersection::split_at_self_intersections`] on it.
///   That function already handles the formal split contract correctly
///   for the same-boundary case.
/// - **Different-boundary case** (`boundary_a != boundary_b`): split each
///   boundary by INSERTING the crossing vertex P at the appropriate
///   position. This converts the transverse crossing into a T-junction
///   at P, which `spade::cdt` handles non-fatally.
///
/// Returns the rebuilt boundary set. The crossing is consumed from the
/// set in the sense that the geometric configuration that produced it
/// is replaced by either a split or a T-junction; the next call to
/// [`find_first_set_crossing`] will see fewer transverse crossings (or
/// none, if this was the only one).
fn split_at_crossing(boundaries: Vec<Pline>, crossing: &SetCrossing) -> Result<Vec<Pline>> {
    if crossing.boundary_a == crossing.boundary_b {
        // Same-boundary: delegate to the per-boundary splitter, which
        // produces two simple loops (or more, if the loop had multiple
        // self-intersections — though we resolve the first crossing
        // each work-loop iteration, so at most two).
        let mut out = Vec::with_capacity(boundaries.len() + 1);
        for (idx, b) in boundaries.into_iter().enumerate() {
            if idx == crossing.boundary_a {
                let split = split_at_self_intersections(b)?;
                out.extend(split.into_iter().filter(|p| p.vertices.len() >= 3));
            } else {
                out.push(b);
            }
        }
        Ok(out)
    } else {
        // Different boundaries: insert the crossing vertex into both.
        // Preserve the order of all other boundaries.
        let p = PlineVertex::line(crossing.x, crossing.y);
        let mut out = Vec::with_capacity(boundaries.len());
        for (idx, mut b) in boundaries.into_iter().enumerate() {
            if idx == crossing.boundary_a {
                insert_vertex_after_edge(&mut b, crossing.edge_a, p);
            } else if idx == crossing.boundary_b {
                insert_vertex_after_edge(&mut b, crossing.edge_b, p);
            }
            out.push(b);
        }
        Ok(out)
    }
}

/// Insert vertex `p` into `pline.vertices` immediately after the start
/// of edge `k`. Edge `k` connects `vertices[k]` to `vertices[(k+1) % n]`,
/// so the inserted vertex lands at index `k + 1` (i.e. between the
/// edge's two endpoints).
///
/// For closed plines (which is always the case in `make_tessellation_safe`),
/// this works for every edge index `0..n` because `Vec::insert` accepts
/// indices up to `len`.
fn insert_vertex_after_edge(pline: &mut Pline, edge: usize, p: PlineVertex) {
    let insert_at = edge + 1;
    debug_assert!(
        insert_at <= pline.vertices.len(),
        "insert_vertex_after_edge: edge index out of range"
    );
    pline.vertices.insert(insert_at, p);
}

// ---------------------------------------------------------------------------
// T9: CDT dry-run verification
// ---------------------------------------------------------------------------

/// Defense-in-depth: build a fresh
/// [`spade::ConstrainedDelaunayTriangulation`], insert every vertex of
/// every boundary, then attempt to add each boundary edge as a constraint
/// via `try_add_constraint`. Return `Err` on the first vertex insertion
/// or constraint-edge rejection.
///
/// This catches spade-reject cases that our transverse-crossing detection
/// missed (e.g. due to floating-point noise putting a near-crossing just
/// outside the `PARAM_EPS` interior bound). Cheap relative to a full
/// tessellation; expensive enough that we run it only once at the end of
/// [`make_tessellation_safe`] after all cleanup + splitting.
///
/// Uses spade's [`try_add_constraint`](spade::ConstrainedDelaunayTriangulation::try_add_constraint)
/// rather than the panicking `add_constraint`. On a constraint conflict
/// `try_add_constraint` returns an empty `Vec<FixedDirectedEdgeHandle>`
/// without modifying the triangulation, which we surface as
/// `OperationError::Failed`.
fn verify_cdt_safe(boundaries: &[Pline]) -> Result<()> {
    use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};

    let mut cdt: ConstrainedDelaunayTriangulation<Point2<f64>> =
        ConstrainedDelaunayTriangulation::new();

    for (bi, b) in boundaries.iter().enumerate() {
        let n = b.vertices.len();
        if n < 3 {
            continue;
        }

        // Insert every vertex first. Spade's insert tolerates duplicates
        // (returns the existing handle for a coincident point), so this
        // is safe even when separate boundaries share a T-junction
        // vertex inserted by T8.
        let mut handles = Vec::with_capacity(n);
        for (vi, v) in b.vertices.iter().enumerate() {
            match cdt.insert(Point2::new(v.x, v.y)) {
                Ok(h) => handles.push(h),
                Err(e) => {
                    return Err(OperationError::Failed(format!(
                        "tessellation_safety: CDT dry-run rejected vertex \
                         insert (boundary {bi}, vertex {vi}): {e:?}"
                    ))
                    .into());
                }
            }
        }

        // Add each boundary edge as a constraint. Edge k connects
        // vertices[k] to vertices[(k+1) % n].
        for k in 0..n {
            let from = handles[k];
            let to = handles[(k + 1) % n];
            if from == to {
                // Coincident vertices map to the same handle — skip
                // (zero-length edge, not a constraint).
                continue;
            }
            let added = cdt.try_add_constraint(from, to);
            if added.is_empty() {
                return Err(OperationError::Failed(format!(
                    "tessellation_safety: CDT dry-run rejected constraint edge \
                     (boundary {bi}, edge {k}): would intersect an existing \
                     constraint edge"
                ))
                .into());
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::pline::PlineVertex;

    fn closed_unit_square() -> Pline {
        Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(1.0, 0.0),
                PlineVertex::line(1.0, 1.0),
                PlineVertex::line(0.0, 1.0),
            ],
            closed: true,
        }
    }

    #[test]
    fn make_tessellation_safe_empty_input_returns_empty_output() {
        let result = make_tessellation_safe(vec![]).expect("empty input should not error");
        assert!(result.is_empty());
    }

    #[test]
    fn make_tessellation_safe_passthrough_simple_input() {
        // Single simple boundary passes through unchanged after T7 cleanup
        // (no near-dups, no collinear chains, no closure issue).
        let input = vec![closed_unit_square()];
        let expected_len = input[0].vertices.len();
        let result = make_tessellation_safe(input).expect("simple input should not error");
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vertices.len(), expected_len);
    }

    // --- T7: collapse_near_duplicate_vertices ---

    #[test]
    fn collapse_drops_consecutive_duplicates() {
        // A unit square with a near-duplicate of v0 inserted right after v0
        // (within WALL_EPS).
        let mut p = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(1e-9, 1e-9), // < WALL_EPS = 1e-6 from v0
                PlineVertex::line(1.0, 0.0),
                PlineVertex::line(1.0, 1.0),
                PlineVertex::line(0.0, 1.0),
            ],
            closed: true,
        };
        collapse_near_duplicate_vertices(&mut p, WALL_EPS);
        assert_eq!(
            p.vertices.len(),
            4,
            "near-dup of v0 should be dropped, leaving {{v0, v2, v3, v4}}"
        );
    }

    #[test]
    fn collapse_drops_wraparound_duplicate_for_closed_pline() {
        // A unit square whose last vertex is a near-dup of the first.
        // Closed-pline wrap-around dedup must drop the trailing dup.
        let mut p = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(1.0, 0.0),
                PlineVertex::line(1.0, 1.0),
                PlineVertex::line(0.0, 1.0),
                PlineVertex::line(1e-9, 1e-9), // wraps back to v0
            ],
            closed: true,
        };
        collapse_near_duplicate_vertices(&mut p, WALL_EPS);
        assert_eq!(
            p.vertices.len(),
            4,
            "trailing wrap-around dup should be dropped"
        );
    }

    // --- T7: simplify_collinear_chain ---

    #[test]
    fn simplify_drops_three_in_a_row_collinear_midpoint() {
        // Closed pline with a midpoint inserted on the bottom edge.
        let mut p = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(0.5, 0.0), // collinear midpoint on bottom edge
                PlineVertex::line(1.0, 0.0),
                PlineVertex::line(1.0, 1.0),
                PlineVertex::line(0.0, 1.0),
            ],
            closed: true,
        };
        simplify_collinear_chain(&mut p, WALL_EPS);
        assert_eq!(
            p.vertices.len(),
            4,
            "collinear midpoint on bottom edge should be dropped"
        );
    }

    #[test]
    fn simplify_keeps_real_corners_of_unit_square() {
        let mut p = closed_unit_square();
        simplify_collinear_chain(&mut p, WALL_EPS);
        assert_eq!(p.vertices.len(), 4, "all four corners should be kept");
    }

    // --- T7: validate_and_close ---

    #[test]
    fn validate_pops_trailing_duplicate_of_first_for_closed_pline() {
        let mut p = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(1.0, 0.0),
                PlineVertex::line(1.0, 1.0),
                PlineVertex::line(0.0, 1.0),
                PlineVertex::line(1e-9, 1e-9), // near-dup of first
            ],
            closed: true,
        };
        let ok = validate_and_close(&mut p, WALL_EPS);
        assert!(ok);
        assert_eq!(
            p.vertices.len(),
            4,
            "trailing near-duplicate of first should be popped"
        );
    }

    #[test]
    fn validate_returns_false_for_degenerate_closed_pline() {
        let mut p = Pline {
            vertices: vec![PlineVertex::line(0.0, 0.0), PlineVertex::line(1.0, 0.0)],
            closed: true,
        };
        let ok = validate_and_close(&mut p, WALL_EPS);
        assert!(
            !ok,
            "closed pline with < 3 vertices is degenerate; validate should return false"
        );
    }

    // --- T7 integration ---

    #[test]
    fn make_tessellation_safe_drops_degenerate_boundary_silently() {
        // 2-vertex closed Pline is degenerate — should be filtered out, not
        // cause an Err.
        let degenerate = Pline {
            vertices: vec![PlineVertex::line(0.0, 0.0), PlineVertex::line(1.0, 0.0)],
            closed: true,
        };
        let valid = closed_unit_square();
        let result =
            make_tessellation_safe(vec![degenerate, valid.clone()]).expect("should not error");
        assert_eq!(result.len(), 1, "only the valid boundary should remain");
        assert_eq!(result[0].vertices.len(), valid.vertices.len());
    }

    #[test]
    fn make_tessellation_safe_runs_full_t7_cleanup_pipeline() {
        // Single boundary that exercises all three cleanup phases:
        // - one near-duplicate vertex (collapse drops it)
        // - one collinear midpoint (simplify drops it)
        // - one trailing near-dup-of-first (validate pops it)
        let messy = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(1e-9, 1e-9), // near-dup of v0 → collapse
                PlineVertex::line(0.5, 0.0),   // collinear midpoint on bottom → simplify
                PlineVertex::line(1.0, 0.0),
                PlineVertex::line(1.0, 1.0),
                PlineVertex::line(0.0, 1.0),
                PlineVertex::line(1e-9, 1e-9), // wrap-back-to-v0 → collapse OR validate
            ],
            closed: true,
        };
        let result = make_tessellation_safe(vec![messy]).expect("should not error");
        assert_eq!(result.len(), 1);
        assert_eq!(
            result[0].vertices.len(),
            4,
            "all three cleanup phases should reduce to 4 unit-square corners"
        );
    }

    // --- T8: find_first_set_crossing ---

    fn closed_pline_xy(pts: &[(f64, f64)]) -> Pline {
        Pline {
            vertices: pts.iter().map(|&(x, y)| PlineVertex::line(x, y)).collect(),
            closed: true,
        }
    }

    /// Run `find_first_set_crossing` until it returns `None`, counting
    /// the iterations needed. Used by T8 integration tests to verify the
    /// output of `make_tessellation_safe` is fully crossing-free.
    fn count_remaining_crossings(boundaries: &[Pline]) -> usize {
        let mut count = 0;
        let mut bs: Vec<Pline> = boundaries.to_vec();
        while let Some(c) = find_first_set_crossing(&bs) {
            count += 1;
            // Sanity guard so a buggy splitter cannot loop forever in a
            // test helper.
            if count > 10_000 {
                panic!("count_remaining_crossings: runaway");
            }
            bs = split_at_crossing(bs, &c).expect("split should not error in helper");
        }
        count
    }

    #[test]
    fn find_first_set_crossing_simple_set_returns_none() {
        // Two disjoint unit squares, well separated.
        let s0 = closed_pline_xy(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]);
        let s1 = closed_pline_xy(&[(5.0, 5.0), (6.0, 5.0), (6.0, 6.0), (5.0, 6.0)]);
        assert!(find_first_set_crossing(&[s0, s1]).is_none());
    }

    #[test]
    fn find_first_set_crossing_intra_boundary_figure_eight() {
        // Single closed figure-8 pline: edges 0 and 2 cross at (1,1).
        let p = closed_pline_xy(&[(0.0, 0.0), (2.0, 2.0), (0.0, 2.0), (2.0, 0.0)]);
        let c = find_first_set_crossing(&[p]).expect("figure-8 should self-cross");
        assert_eq!(c.boundary_a, 0);
        assert_eq!(c.boundary_b, 0);
        assert_eq!(c.edge_a, 0);
        assert_eq!(c.edge_b, 2);
        assert!((c.x - 1.0).abs() < 1e-9, "x expected 1.0, got {}", c.x);
        assert!((c.y - 1.0).abs() < 1e-9, "y expected 1.0, got {}", c.y);
    }

    #[test]
    fn find_first_set_crossing_inter_boundary_two_overlapping_squares() {
        // Two unit squares: square 0 at (0,0)-(2,2), square 1 at (1,1)-(3,3).
        // Edge 1 of square 0 (right side, (2,0)→(2,2)) crosses edge 0 of
        // square 1 (bottom, (1,1)→(3,1)) at (2,1).
        // Edge 2 of square 0 (top, (2,2)→(0,2)) crosses edge 3 of square 1
        // (left, (1,3)→(1,1)) at (1,2).
        // Scan order: ba=0, bb=1, then i ascending, then j ascending.
        // For i=1 (the right side of square 0), the first j hit is j=0
        // (bottom of square 1) → crossing at (2, 1).
        let s0 = closed_pline_xy(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]);
        let s1 = closed_pline_xy(&[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)]);
        let c = find_first_set_crossing(&[s0, s1]).expect("squares should cross");
        assert_eq!(c.boundary_a, 0);
        assert_eq!(c.boundary_b, 1);
        // Earlier i (right side, i=1) is hit before later i (top, i=2).
        assert_eq!(c.edge_a, 1);
        assert_eq!(c.edge_b, 0);
        assert!((c.x - 2.0).abs() < 1e-9, "x expected 2.0, got {}", c.x);
        assert!((c.y - 1.0).abs() < 1e-9, "y expected 1.0, got {}", c.y);
    }

    #[test]
    fn make_tessellation_safe_resolves_intra_figure_eight() {
        // Figure-8 alone — the per-boundary case.
        let p = closed_pline_xy(&[(0.0, 0.0), (2.0, 2.0), (0.0, 2.0), (2.0, 0.0)]);
        // Sanity: input has at least one crossing.
        assert!(
            find_first_set_crossing(&[p.clone()]).is_some(),
            "figure-8 input should have a detectable crossing"
        );

        let result = make_tessellation_safe(vec![p]).expect("figure-8 should resolve");
        assert!(
            result.len() > 1,
            "figure-8 should split into more than one boundary; got {}",
            result.len()
        );
        for (idx, b) in result.iter().enumerate() {
            assert!(b.closed, "output[{idx}] should be closed");
        }
        assert_eq!(
            count_remaining_crossings(&result),
            0,
            "output should have zero remaining crossings"
        );
    }

    #[test]
    fn make_tessellation_safe_resolves_inter_overlapping_squares() {
        // Two overlapping unit squares (same as inter-boundary test above).
        let s0 = closed_pline_xy(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]);
        let s1 = closed_pline_xy(&[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)]);
        // Sanity: input has at least one crossing.
        assert!(
            find_first_set_crossing(&[s0.clone(), s1.clone()]).is_some(),
            "overlapping squares should have a detectable crossing"
        );

        let result =
            make_tessellation_safe(vec![s0, s1]).expect("overlapping squares should resolve");
        assert!(!result.is_empty());
        assert_eq!(
            count_remaining_crossings(&result),
            0,
            "output should have zero remaining crossings"
        );
    }

    // --- T9: verify_cdt_safe ---

    #[test]
    fn verify_cdt_safe_passes_simple_unit_square() {
        // A single simple unit square — must pass the CDT dry-run.
        let s = closed_unit_square();
        verify_cdt_safe(&[s]).expect("simple unit square should pass CDT dry-run");
    }

    #[test]
    fn verify_cdt_safe_rejects_unresolved_crossing() {
        // Manually constructed figure-8 — has a transverse self-crossing.
        // Bypass T7+T8 by calling verify_cdt_safe directly to confirm the
        // CDT dry-run rejects the unresolved crossing.
        let figure_eight = closed_pline_xy(&[(0.0, 0.0), (2.0, 2.0), (0.0, 2.0), (2.0, 0.0)]);
        let err = verify_cdt_safe(&[figure_eight])
            .expect_err("figure-8 with unresolved crossing should be rejected");
        let msg = format!("{err:?}");
        assert!(
            msg.contains("CDT dry-run rejected"),
            "error message should mention CDT dry-run rejection; got {msg}"
        );
    }
}
