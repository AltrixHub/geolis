//! Shared types and geometric primitives for the 2D boolean engine.
//!
//! All consumers of `boolean_2d` (union, subtract, future ops) share the
//! same [`PolygonWithHoles`] input/output shape, the same global epsilon
//! [`WALL_EPS`], and the same low-level point-in-polygon classifier.

use crate::math::distance_2d::point_to_segment_dist;

/// Single epsilon for all geometric decisions in the 2D boolean pipeline.
///
/// Inherited verbatim from the original `wall_outline::polygon_union`
/// engine so that all existing wall-outline regression fixtures remain
/// bit-identical after the move into `boolean_2d`.
pub const WALL_EPS: f64 = 1e-6;

/// `WALL_EPS` squared, used wherever a squared distance is compared
/// against a squared tolerance (avoids unnecessary `sqrt` calls in hot
/// loops).
pub const WALL_EPS_SQ: f64 = WALL_EPS * WALL_EPS;

/// A simple closed 2D polygon described by its ordered vertex list
/// (no explicit closing duplicate — the last vertex is implicitly
/// connected back to the first).
pub type Polygon = Vec<(f64, f64)>;

/// A planar face described by an outer boundary and zero or more holes.
///
/// # Winding contract
/// - `outer` is CCW with `signed_area(outer) > 0`.
/// - Each `holes[i]` is CW with `signed_area(holes[i]) < 0`.
/// - Every hole is fully contained in `outer`.
/// - Sibling holes are non-overlapping.
///
/// The contract is enforced by the engine's `assemble_faces` step for
/// every output of [`crate::operations::boolean_2d::union_all_with_holes`]
/// and [`crate::operations::boolean_2d::subtract_all_with_holes`].
/// External callers constructing a `PolygonWithHoles` directly are
/// responsible for upholding the invariants (use the higher-level
/// validator types — e.g. `WallFootprint2D::try_from_parts` — when
/// crossing crate boundaries with untrusted input).
#[derive(Clone, Debug, PartialEq)]
pub struct PolygonWithHoles {
    pub outer: Polygon,
    pub holes: Vec<Polygon>,
}

impl PolygonWithHoles {
    pub fn into_parts(self) -> (Polygon, Vec<Polygon>) {
        (self.outer, self.holes)
    }
}

/// Result of a boolean union: typed face topology.
pub struct UnionResult {
    pub faces: Vec<PolygonWithHoles>,
}

/// Three-valued classification of a point relative to a single ring.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointClass {
    Inside,
    Outside,
    Boundary,
}

/// Robust point-in-polygon classifier with a `WALL_EPS` boundary band.
///
/// Returns [`PointClass::Boundary`] when `p` lies within `WALL_EPS` of
/// any polygon edge. Otherwise performs a winding-number test and
/// reports [`PointClass::Inside`] (winding != 0) or
/// [`PointClass::Outside`].
#[must_use]
pub fn point_in_polygon_class(p: (f64, f64), poly: &Polygon) -> PointClass {
    let n = poly.len();
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        let dist = point_to_segment_dist(p.0, p.1, a.0, a.1, b.0, b.1);
        if dist < WALL_EPS {
            return PointClass::Boundary;
        }
    }
    let mut winding = 0i32;
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        if a.1 <= p.1 {
            if b.1 > p.1 && cross_2d(a, b, p) > 0.0 {
                winding += 1;
            }
        } else if b.1 <= p.1 && cross_2d(a, b, p) < 0.0 {
            winding -= 1;
        }
    }
    if winding != 0 {
        PointClass::Inside
    } else {
        PointClass::Outside
    }
}

fn cross_2d(a: (f64, f64), b: (f64, f64), p: (f64, f64)) -> f64 {
    (b.0 - a.0) * (p.1 - a.1) - (b.1 - a.1) * (p.0 - a.0)
}

/// Shoelace signed area. CCW > 0, CW < 0.
#[must_use]
pub fn signed_area(poly: &Polygon) -> f64 {
    let n = poly.len();
    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        area += poly[i].0 * poly[j].1;
        area -= poly[j].0 * poly[i].1;
    }
    area * 0.5
}
