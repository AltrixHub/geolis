//! 2D curve×curve intersection.
//!
//! Seeds come from bounding-box subdivision ([`super::bbox`]); each seed feeds a
//! 2×2 Newton iteration on `F(t_a, t_b) = a(t_a) - b(t_b)` whose Jacobian is
//! `[a'(t_a), -b'(t_b)]`. Results are clamped to the parameter domains and
//! deduplicated in parameter space.

use nalgebra::{Matrix2, Vector2};

use crate::error::Result;
use crate::geometry::nurbs::NurbsCurve2D;
use crate::math::Point2;

use super::bbox::seed_curve_curve;
use super::types::{CurveCurveIntersection2D, IntersectionOptions};

/// Computes the transversal intersection points of two 2D NURBS curves.
///
/// Tangential contacts (e.g. a line touching a circle) are reported as a single
/// point; the parameter-space dedup prevents duplicate or near-duplicate
/// outputs. Disjoint curves yield an empty vector.
///
/// # Errors
///
/// Returns an error only if a sub-curve construction during seeding or a curve
/// evaluation fails (e.g. a vanishing rational denominator).
pub fn intersect_curves_2d(
    a: &NurbsCurve2D,
    b: &NurbsCurve2D,
    options: &IntersectionOptions,
) -> Result<Vec<CurveCurveIntersection2D>> {
    let (a_min, a_max) = a.parameter_domain();
    let (b_min, b_max) = b.parameter_domain();

    // Leaf extent: a small fraction of the combined hull diagonal keeps seeds
    // tight without over-subdividing. The pad is a tiny absolute epsilon so
    // tangential hulls survive without inflating boxes into spurious overlaps.
    let leaf = seed_leaf_extent(a, b);
    let pad = 1e-7;
    let seeds = seed_curve_curve(a, b, leaf, pad, 40)?;

    let mut results: Vec<CurveCurveIntersection2D> = Vec::new();
    for (ba, bb) in seeds {
        let mut ta = 0.5 * (ba.t0 + ba.t1);
        let mut tb = 0.5 * (bb.t0 + bb.t1);
        let mut converged = false;
        for _ in 0..options.max_iterations {
            let da = a.derivatives(ta, 1)?;
            let db = b.derivatives(tb, 1)?;
            let f = Vector2::new(da[0].x - db[0].x, da[0].y - db[0].y);
            if f.norm() < options.tolerance {
                converged = true;
                break;
            }
            // J = [a'(ta), -b'(tb)] (columns).
            let j = Matrix2::new(da[1].x, -db[1].x, da[1].y, -db[1].y);
            let det = j.determinant();
            if det.abs() < 1e-14 {
                // Singular Jacobian: parallel tangents — abandon this seed.
                break;
            }
            let Some(jinv) = j.try_inverse() else {
                break;
            };
            let delta = jinv * f;
            ta = (ta - delta.x).clamp(a_min, a_max);
            tb = (tb - delta.y).clamp(b_min, b_max);
            if delta.norm() < options.tolerance {
                // Re-check residual after the clamped step.
                let pa = a.point_at(ta)?;
                let pb = b.point_at(tb)?;
                converged = (pa - pb).norm() < options.tolerance.max(1e-7);
                break;
            }
        }
        if !converged {
            continue;
        }
        let pa = a.point_at(ta)?;
        let pb = b.point_at(tb)?;
        if (pa - pb).norm() > options.tolerance.max(1e-7) {
            continue;
        }
        let point = Point2::from((pa.coords + pb.coords) * 0.5);
        push_unique(
            &mut results,
            CurveCurveIntersection2D {
                point,
                t_a: ta,
                t_b: tb,
            },
        );
    }
    Ok(results)
}

/// A leaf-extent heuristic: 1/16 of the larger control-hull diagonal. Coarse
/// enough to keep the candidate count bounded (Newton has a wide basin for the
/// transversal cases) yet fine enough to separate distinct crossings.
fn seed_leaf_extent(a: &NurbsCurve2D, b: &NurbsCurve2D) -> f64 {
    let (a_lo, a_hi) = a.bounding_box();
    let (b_lo, b_hi) = b.bounding_box();
    let da = (a_hi - a_lo).norm();
    let db = (b_hi - b_lo).norm();
    (da.max(db) / 16.0).max(1e-6)
}

/// Inserts a result unless an existing one coincides.
///
/// Two dedup scales are used: a tight parameter-space test for exact duplicate
/// convergence, and a looser geometric test that collapses the symmetric pair of
/// near-solutions a tangency produces (where Newton is ill-conditioned and lands
/// just off the true tangent point on either side) into a single reported hit.
fn push_unique(results: &mut Vec<CurveCurveIntersection2D>, candidate: CurveCurveIntersection2D) {
    const PARAM_DEDUP: f64 = 1e-6;
    const GEOM_DEDUP: f64 = 1e-4;
    for r in results.iter() {
        if (r.t_a - candidate.t_a).abs() < PARAM_DEDUP
            && (r.t_b - candidate.t_b).abs() < PARAM_DEDUP
        {
            return;
        }
        if (r.point - candidate.point).norm() < GEOM_DEDUP {
            return;
        }
    }
    results.push(candidate);
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::nurbs::{KnotVector, NurbsCurve, NurbsCurve3D};
    use crate::math::{Point2, Point3, Vector3};

    fn line_2d(p0: Point2, p1: Point2) -> NurbsCurve2D {
        NurbsCurve2D::from_unweighted(
            vec![p0, p1],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
        )
        .unwrap()
    }

    /// A full 2D rational circle, built from the exact 3D nine-point pattern's
    /// XY data (z is identically 0 for an XY-plane circle).
    fn circle_2d(center: Point2, radius: f64) -> NurbsCurve2D {
        let c3 = NurbsCurve3D::circle(
            Point3::new(center.x, center.y, 0.0),
            radius,
            Vector3::z(),
            Vector3::x(),
        )
        .unwrap();
        let pts: Vec<Point2> = c3
            .control_points()
            .iter()
            .map(|p| Point2::new(p.x, p.y))
            .collect();
        NurbsCurve::<2>::new(pts, c3.weights().to_vec(), c3.knots().clone(), c3.degree()).unwrap()
    }

    #[test]
    fn two_lines_cross_at_known_point() {
        let a = line_2d(Point2::new(0.0, 0.0), Point2::new(2.0, 2.0));
        let b = line_2d(Point2::new(0.0, 2.0), Point2::new(2.0, 0.0));
        let hits = intersect_curves_2d(&a, &b, &IntersectionOptions::default()).unwrap();
        assert_eq!(hits.len(), 1, "expected one crossing");
        assert!((hits[0].point - Point2::new(1.0, 1.0)).norm() < 1e-9);
        assert!((hits[0].t_a - 0.5).abs() < 1e-9);
        assert!((hits[0].t_b - 0.5).abs() < 1e-9);
    }

    #[test]
    fn line_through_circle_hits_two_points() {
        // Unit circle at origin, horizontal line y = 0 through it: hits (-1,0),(1,0).
        let circle = circle_2d(Point2::new(0.0, 0.0), 1.0);
        let line = line_2d(Point2::new(-2.0, 0.0), Point2::new(2.0, 0.0));
        let hits = intersect_curves_2d(&circle, &line, &IntersectionOptions::default()).unwrap();
        assert_eq!(hits.len(), 2, "expected two crossings, got {}", hits.len());
        for h in &hits {
            assert!(
                (h.point.coords.norm() - 1.0).abs() < 1e-8,
                "off circle: {:?}",
                h.point
            );
            assert!(h.point.y.abs() < 1e-8, "off line: {:?}", h.point);
        }
        // One near (1,0), one near (-1,0).
        assert!(hits
            .iter()
            .any(|h| (h.point - Point2::new(1.0, 0.0)).norm() < 1e-7));
        assert!(hits
            .iter()
            .any(|h| (h.point - Point2::new(-1.0, 0.0)).norm() < 1e-7));
    }

    #[test]
    fn disjoint_curves_yield_empty() {
        let a = line_2d(Point2::new(0.0, 0.0), Point2::new(1.0, 0.0));
        let b = line_2d(Point2::new(0.0, 5.0), Point2::new(1.0, 5.0));
        let hits = intersect_curves_2d(&a, &b, &IntersectionOptions::default()).unwrap();
        assert!(hits.is_empty());
    }

    #[test]
    fn separated_circle_and_line_yield_empty() {
        let circle = circle_2d(Point2::new(0.0, 0.0), 1.0);
        let line = line_2d(Point2::new(-2.0, 3.0), Point2::new(2.0, 3.0));
        let hits = intersect_curves_2d(&circle, &line, &IntersectionOptions::default()).unwrap();
        assert!(hits.is_empty());
    }

    #[test]
    fn tangent_line_touches_circle_once() {
        // Line y = 1 is tangent to the unit circle at (0, 1).
        //
        // Documented behavior: a tangency is an ill-conditioned root of
        // F = a - b (the curve tangents are parallel, so the Newton Jacobian is
        // near-singular). The solver still terminates and reports a *single*
        // hit (the symmetric near-solutions on either side of the true tangent
        // collapse via the geometric dedup), close to but not bit-exact at the
        // tangent point. Transversal crossings remain exact; only tangencies
        // carry this reduced precision.
        let circle = circle_2d(Point2::new(0.0, 0.0), 1.0);
        let line = line_2d(Point2::new(-2.0, 1.0), Point2::new(2.0, 1.0));
        let hits = intersect_curves_2d(&circle, &line, &IntersectionOptions::default()).unwrap();
        assert_eq!(
            hits.len(),
            1,
            "tangency must report exactly one point, got {}",
            hits.len()
        );
        assert!(
            (hits[0].point - Point2::new(0.0, 1.0)).norm() < 1e-3,
            "tangent point {:?} not near (0, 1)",
            hits[0].point
        );
    }

    #[test]
    fn full_circle_parameter_wraps_are_deduped() {
        // A line through the circle center crosses twice; the seam at t=0/1 must
        // not produce a spurious third hit.
        let circle = circle_2d(Point2::new(0.0, 0.0), 1.0);
        let line = line_2d(Point2::new(-2.0, 0.0), Point2::new(2.0, 0.0));
        let opts = IntersectionOptions {
            max_iterations: 80,
            ..IntersectionOptions::default()
        };
        let hits = intersect_curves_2d(&circle, &line, &opts).unwrap();
        assert_eq!(hits.len(), 2);
    }
}
