//! 2D polygon boolean-intersection via the shared arrangement engine.
//!
//! [`intersect_all_with_holes`] computes `base ∩ (⋃others)` — the parts
//! of `base` that are covered by at least one of `others[i]`. Output is
//! typed face topology (zero, one, or many
//! [`PolygonWithHoles`]). The engine is identical to union / subtract;
//! only the fill oracle differs (see [`super::engine::IntersectOracle`]).

use crate::error::Result;

use super::engine::{IntersectOracle, run_arrangement};
use super::types::PolygonWithHoles;

/// Intersect a base region with the union of `others`. Returns the
/// covered parts of `base` as typed face topology.
///
/// Semantics: `result = base ∩ (⋃others)`.
///
/// Special cases:
/// - `others.is_empty()` returns an empty `Vec` (intersection with
///   nothing is empty).
/// - `others` fully covering `base` returns `vec![base]`.
/// - `others` disjoint from `base` returns an empty `Vec`.
///
/// # Errors
///
/// Propagates [`crate::error::OperationError::Failed`] from the
/// arrangement engine on the same degenerate-input cases as
/// [`super::union_all_with_holes`].
pub fn intersect_all_with_holes(
    base: PolygonWithHoles,
    others: &[PolygonWithHoles],
) -> Result<Vec<PolygonWithHoles>> {
    if others.is_empty() {
        return Ok(Vec::new());
    }
    let mut segment_inputs: Vec<PolygonWithHoles> = Vec::with_capacity(1 + others.len());
    segment_inputs.push(base.clone());
    segment_inputs.extend(others.iter().cloned());

    let oracle = IntersectOracle {
        base: &base,
        others,
    };
    run_arrangement(&segment_inputs, &oracle)
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::super::types::{Polygon, signed_area};
    use super::*;

    fn rect(x: f64, y: f64, w: f64, h: f64) -> Polygon {
        vec![(x, y), (x + w, y), (x + w, y + h), (x, y + h)]
    }

    fn pwh(outer: Polygon) -> PolygonWithHoles {
        PolygonWithHoles {
            outer,
            holes: Vec::new(),
        }
    }

    fn total_area(faces: &[PolygonWithHoles]) -> f64 {
        faces
            .iter()
            .map(|f| {
                signed_area(&f.outer).abs()
                    - f.holes.iter().map(|h| signed_area(h).abs()).sum::<f64>()
            })
            .sum()
    }

    #[test]
    fn overlapping_rects_intersect_to_the_shared_area() {
        let base = pwh(rect(0.0, 0.0, 10.0, 10.0));
        let other = pwh(rect(5.0, 5.0, 10.0, 10.0));
        let result = intersect_all_with_holes(base, &[other]).expect("intersect");
        assert_eq!(result.len(), 1);
        assert!((total_area(&result) - 25.0).abs() < 1e-9);
    }

    #[test]
    fn disjoint_rects_intersect_to_nothing() {
        let base = pwh(rect(0.0, 0.0, 4.0, 4.0));
        let other = pwh(rect(10.0, 10.0, 4.0, 4.0));
        let result = intersect_all_with_holes(base, &[other]).expect("intersect");
        assert!(result.is_empty());
    }

    #[test]
    fn empty_others_yield_empty_intersection() {
        let base = pwh(rect(0.0, 0.0, 4.0, 4.0));
        let result = intersect_all_with_holes(base, &[]).expect("intersect");
        assert!(result.is_empty());
    }

    #[test]
    fn full_cover_returns_base_area() {
        let base = pwh(rect(2.0, 2.0, 4.0, 4.0));
        let other = pwh(rect(0.0, 0.0, 10.0, 10.0));
        let result = intersect_all_with_holes(base.clone(), &[other]).expect("intersect");
        assert_eq!(result.len(), 1);
        assert!((total_area(&result) - 16.0).abs() < 1e-9);
    }

    #[test]
    fn union_of_others_drives_coverage() {
        // Two others jointly cover the base's left and right halves.
        let base = pwh(rect(0.0, 0.0, 10.0, 4.0));
        let left = pwh(rect(-1.0, -1.0, 6.0, 6.0));
        let right = pwh(rect(5.0, -1.0, 6.0, 6.0));
        let result = intersect_all_with_holes(base, &[left, right]).expect("intersect");
        assert!((total_area(&result) - 40.0).abs() < 1e-9, "full base covered");
    }
}
