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

use super::Pline;
use crate::error::Result;

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
    reason = "Wired into WallOutline2D::execute in plan-13k T10; T6 lands the scaffold first to keep commits small."
)]
#[allow(
    clippy::unnecessary_wraps,
    reason = "T6 scaffold returns Ok unconditionally; T7-T9 add real Err paths (cleanup failure, work-loop bail, CDT dry-run rejection)."
)]
pub(crate) fn make_tessellation_safe(boundaries: Vec<Pline>) -> Result<Vec<Pline>> {
    // T6 is scaffolding only. Each subsequent task fills in one phase:
    //   T7: per-boundary cleanup
    //   T8: cross-boundary crossing split
    //   T9: CDT dry-run verification
    //   T10: integration into WallOutline2D::execute
    //
    // Until T7 lands, this function is a no-op pass-through: every input
    // boundary appears in the output unchanged. This is safe because the
    // existing `split_at_self_intersections` integration in
    // `WallOutline2D::execute` (T3) is still in place and will be
    // replaced by `make_tessellation_safe` only at T10. Between T6 and
    // T10, this module compiles but is unused (#[allow(dead_code)]).
    Ok(boundaries)
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
        // T6 scaffold: single simple boundary passes through unchanged.
        // This assertion is intentionally weak and will be tightened as
        // T7-T9 add cleanup steps.
        let input = vec![closed_unit_square()];
        let result = make_tessellation_safe(input.clone()).expect("simple input should not error");
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vertices.len(), input[0].vertices.len());
    }
}
