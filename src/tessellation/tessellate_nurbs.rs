//! Adaptive tessellation of NURBS curves and surfaces.
//!
//! Curves are sampled by recursive chord-deviation bisection, seeded at
//! distinct interior knots so degree-reducing discontinuities are honored.
//! Surfaces are meshed onto a uniformly refined grid whose division count is
//! doubled per direction until adjacent-sample normal deviation falls under
//! tolerance, also seeded with interior knot lines.

use crate::error::{Result, TessellationError};
use crate::geometry::nurbs::NurbsCurve3D;
use crate::math::Point3;

/// Options for adaptive NURBS curve tessellation.
#[derive(Debug, Clone, Copy)]
pub struct CurveTessellationOptions {
    /// Maximum chord deviation from the true curve, in model units.
    pub chord_tolerance: f64,
    /// Maximum recursion depth per seed interval.
    pub max_depth: usize,
}

impl Default for CurveTessellationOptions {
    fn default() -> Self {
        Self {
            chord_tolerance: 1e-3,
            max_depth: 16,
        }
    }
}

/// Adaptively samples a NURBS curve into a polyline.
///
/// The parameter domain is first split at every distinct interior knot, then
/// each sub-interval is recursively bisected while the curve's midpoint
/// deviates from the chord midpoint by more than `chord_tolerance` (capped at
/// `max_depth`). Endpoints are always emitted exactly once; the returned points
/// are ordered and contain no consecutive duplicates.
///
/// # Errors
///
/// Returns [`TessellationError::InvalidParameters`] if `chord_tolerance` is not
/// strictly positive, or propagates any evaluation error from the curve.
pub fn tessellate_nurbs_curve(
    curve: &NurbsCurve3D,
    options: &CurveTessellationOptions,
) -> Result<Vec<Point3>> {
    if !(options.chord_tolerance > 0.0) {
        return Err(TessellationError::InvalidParameters(
            "chord_tolerance must be strictly positive".to_owned(),
        )
        .into());
    }

    let (t_min, t_max) = curve.parameter_domain();
    let seeds = seed_parameters(curve, t_min, t_max);

    // Start with the first seed point, then append each interval's interior +
    // end points so endpoints are shared exactly once.
    let mut points = Vec::new();
    points.push(curve.point_at(seeds[0])?);
    for window in seeds.windows(2) {
        let (a, b) = (window[0], window[1]);
        let pa = *points.last().expect("points seeded with first sample");
        let pb = curve.point_at(b)?;
        bisect_curve(curve, a, b, pa, pb, options, 0, &mut points)?;
        points.push(pb);
    }

    Ok(points)
}

/// Distinct parameter seeds: domain endpoints plus every interior knot that
/// lies strictly inside the domain. Returned sorted and de-duplicated, always
/// containing at least the two endpoints.
fn seed_parameters(curve: &NurbsCurve3D, t_min: f64, t_max: f64) -> Vec<f64> {
    let mut seeds = vec![t_min, t_max];
    let span = (t_max - t_min).abs();
    let eps = span * 1e-9;
    for &k in curve.knots().as_slice() {
        if k > t_min + eps && k < t_max - eps {
            seeds.push(k);
        }
    }
    seeds.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    seeds.dedup_by(|a, b| (*a - *b).abs() <= eps);
    seeds
}

/// Recursively appends interior sample points for the open interval `(a, b)`.
///
/// `pa` / `pb` are the already-evaluated endpoints. The midpoint is appended
/// only when the curve deviates from the chord beyond tolerance and depth
/// remains; the endpoint `pb` is appended by the caller.
fn bisect_curve(
    curve: &NurbsCurve3D,
    a: f64,
    b: f64,
    pa: Point3,
    pb: Point3,
    options: &CurveTessellationOptions,
    depth: usize,
    out: &mut Vec<Point3>,
) -> Result<()> {
    let mid_t = 0.5 * (a + b);
    let pm = curve.point_at(mid_t)?;

    let chord_mid = Point3::from(0.5 * (pa.coords + pb.coords));
    let deviation = (pm - chord_mid).norm();

    if depth >= options.max_depth || deviation <= options.chord_tolerance {
        return Ok(());
    }

    bisect_curve(curve, a, mid_t, pa, pm, options, depth + 1, out)?;
    out.push(pm);
    bisect_curve(curve, mid_t, b, pm, pb, options, depth + 1, out)?;
    Ok(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod curve_tests {
    use super::*;
    use crate::geometry::nurbs::{KnotVector, NurbsCurve3D};
    use crate::math::Vector3;

    /// Straight degree-1 curve from (0,0,0) to (5,0,0).
    fn line_curve() -> NurbsCurve3D {
        NurbsCurve3D::from_unweighted(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(5.0, 0.0, 0.0)],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
        )
        .unwrap()
    }

    fn unit_circle() -> NurbsCurve3D {
        NurbsCurve3D::circle(
            Point3::new(0.0, 0.0, 0.0),
            1.0,
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 0.0),
        )
        .unwrap()
    }

    #[test]
    fn straight_line_tessellates_to_two_endpoints() {
        let curve = line_curve();
        let pts = tessellate_nurbs_curve(&curve, &CurveTessellationOptions::default()).unwrap();
        assert_eq!(pts.len(), 2, "degree-1 line must not oversample");
        assert!((pts[0] - Point3::new(0.0, 0.0, 0.0)).norm() < 1e-12);
        assert!((pts[1] - Point3::new(5.0, 0.0, 0.0)).norm() < 1e-12);
    }

    #[test]
    fn rational_circle_within_chord_tolerance() {
        let curve = unit_circle();
        let options = CurveTessellationOptions {
            chord_tolerance: 1e-3,
            max_depth: 16,
        };
        let pts = tessellate_nurbs_curve(&curve, &options).unwrap();

        // Every emitted point lies on the unit circle.
        for p in &pts {
            let r = (p.x * p.x + p.y * p.y).sqrt();
            assert!((r - 1.0).abs() < 1e-9, "point off circle: r = {r}");
            assert!(p.z.abs() < 1e-9, "point off plane: z = {}", p.z);
        }

        // Adjacent-chord deviation: sample the true curve midpoint of each
        // emitted segment and measure its distance to the chord midpoint.
        for w in pts.windows(2) {
            let chord_mid = Point3::from(0.5 * (w[0].coords + w[1].coords));
            // Invert each endpoint's angle to find the bracketing parameters,
            // then evaluate the true midpoint by angle.
            let a0 = w[0].y.atan2(w[0].x);
            let a1 = w[1].y.atan2(w[1].x);
            let mut da = a1 - a0;
            // Normalize to the short arc.
            while da > std::f64::consts::PI {
                da -= std::f64::consts::TAU;
            }
            while da < -std::f64::consts::PI {
                da += std::f64::consts::TAU;
            }
            let mid_angle = a0 + 0.5 * da;
            let true_mid = Point3::new(mid_angle.cos(), mid_angle.sin(), 0.0);
            let deviation = (true_mid - chord_mid).norm();
            assert!(
                deviation < options.chord_tolerance + 1e-9,
                "segment deviation {deviation} exceeds tolerance"
            );
        }

        assert!(pts.len() < 200, "unreasonable point count: {}", pts.len());
        assert!(pts.len() > 4, "circle should be refined: {}", pts.len());
    }

    #[test]
    fn tightening_tolerance_increases_point_count() {
        let curve = unit_circle();
        let coarse = tessellate_nurbs_curve(
            &curve,
            &CurveTessellationOptions {
                chord_tolerance: 1e-2,
                max_depth: 20,
            },
        )
        .unwrap();
        let fine = tessellate_nurbs_curve(
            &curve,
            &CurveTessellationOptions {
                chord_tolerance: 1e-3,
                max_depth: 20,
            },
        )
        .unwrap();
        assert!(
            fine.len() > coarse.len(),
            "tighter tolerance must add points: {} vs {}",
            fine.len(),
            coarse.len()
        );
    }
}
