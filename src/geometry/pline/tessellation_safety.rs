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

use super::{Pline, PlineVertex};
use crate::error::Result;
use crate::operations::offset::wall_outline::polygon_union::WALL_EPS;

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
#[allow(
    clippy::unnecessary_wraps,
    reason = "T7 cleanup never errors today (only filters degenerate boundaries); T8-T9 add real Err paths (work-loop bail, CDT dry-run rejection)."
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

    // T8 (cross-boundary crossing split) and T9 (CDT dry-run) land in
    // subsequent commits.

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
}
