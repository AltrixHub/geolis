//! 2D polygon boolean operations on [`PolygonWithHoles`] inputs.
//!
//! All operations share a single planar-arrangement engine (segment
//! split → vertex snap → bilateral half-edge classification →
//! polar-angle face walk → containment-matrix face assembly) and
//! differ only in the **fill oracle** that decides which result-region
//! a probe point lies in.
//!
//! # Operations
//!
//! | Function | Result region |
//! |---|---|
//! | [`union_all_with_holes`] | `⋃ inputs[i]` (OR-of-PWH-filled) |
//! | [`subtract_all_with_holes`] | `base ∩ (¬⋃ subtracts)` |
//!
//! Both return typed face topology (zero or more
//! [`PolygonWithHoles`]) where every output is guaranteed:
//! - CCW outer with `signed_area > 0`, CW holes with
//!   `signed_area < 0`.
//! - Every hole is fully contained in its outer.
//! - Sibling holes do not overlap.
//! - Boundaries are simple (no self-intersection).
//! - Outputs are CDT-safe (verified in every build; a violation from a
//!   near-degenerate input is returned as [`crate::error::OperationError`],
//!   never a panic).
//! - Determinism: outputs are topologically identical regardless of
//!   input order.
//!
//! # Adding a new operation
//!
//! 1. Implement [`engine::FillOracle`] for the new operation's fill
//!    rule.
//! 2. Add a public entry point that builds the segment-input list and
//!    calls [`engine::run_arrangement`] with the oracle.
//! 3. Add fixtures exercising the new fill rule.

mod engine;
mod intersect;
mod subtract;
mod types;
mod union;

pub use intersect::intersect_all_with_holes;
pub use subtract::subtract_all_with_holes;
pub use types::{
    point_in_polygon_class, signed_area, PointClass, Polygon, PolygonWithHoles, UnionResult,
    WALL_EPS, WALL_EPS_SQ,
};
pub use union::union_all_with_holes;

/// Crate-internal traced union: per-edge [`engine::SegmentSite`] source
/// tracking for callers that derive stable per-segment provenance
/// (`wall_outline::WallOutline2D::execute_faces_with_provenance`).
pub(crate) use engine::{RingRef, SegmentSite, TracedFace};
pub(crate) use union::union_all_with_holes_traced;

/// Crate-internal re-export of the engine's segment-segment intersection
/// primitive. Used by `wall_outline::try_from_parts` ring-validation
/// helpers to share a single `WALL_EPS`-tolerant implementation with
/// the arrangement engine.
pub(crate) use engine::seg_seg_intersect;

/// Crate-internal re-export of the planar arrangement engine.
/// Used by `boolean::merge::merge_component` to compute the union of
/// coplanar face fragments via the same vetted segment-split /
/// vertex-snap / half-edge-classification / face-walk pipeline that
/// drives the 2D polygon booleans, instead of the ad-hoc edge-cancel
/// loop the merge step used previously.
pub(crate) use engine::{run_arrangement, UnionOracle};
