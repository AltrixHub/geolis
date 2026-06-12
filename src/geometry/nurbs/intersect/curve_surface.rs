//! 3D curve×surface intersection.
//!
//! Seeds come from bounding-box subdivision ([`super::bbox::seed_curve_surface`]);
//! each seed feeds a 3×3 Newton iteration on
//! `F(t, u, v) = C(t) - S(u, v) = 0` whose Jacobian columns are
//! `[C'(t), -S_u, -S_v]`. Parameters are clamped to their domains and results are
//! deduplicated in parameter space.

use nalgebra::{Matrix3, Vector3 as NaVector3};

use crate::error::Result;
use crate::geometry::nurbs::{NurbsCurve3D, NurbsSurface};
use crate::math::Point3;

use super::bbox::seed_curve_surface;
use super::types::{CurveSurfaceIntersection, IntersectionOptions};

/// Computes the transversal intersection points of a 3D NURBS curve and a NURBS
/// surface.
///
/// Each result carries the 3D point and both the curve parameter `t` and the
/// surface parameters `(u, v)`. A curve lying tangent to or grazing the surface
/// is an unsupported degeneracy (P4 scope): the solver still terminates, but the
/// reported parameters may carry reduced precision. Disjoint inputs yield an
/// empty vector.
///
/// # Errors
///
/// Returns an error only if a sub-box construction during seeding or a
/// curve/surface evaluation fails (e.g. a vanishing rational denominator).
pub fn intersect_curve_surface(
    curve: &NurbsCurve3D,
    surface: &NurbsSurface,
    options: &IntersectionOptions,
) -> Result<Vec<CurveSurfaceIntersection>> {
    let (t_min, t_max) = curve.parameter_domain();
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_domain();

    let leaf = seed_leaf_extent(curve, surface);
    let pad = 1e-7;
    let seeds = seed_curve_surface(curve, surface, leaf, pad, 40)?;

    let mut results: Vec<CurveSurfaceIntersection> = Vec::new();
    for (cb, sb) in seeds {
        let mut t = 0.5 * (cb.t0 + cb.t1);
        let mut u = 0.5 * (sb.u0 + sb.u1);
        let mut v = 0.5 * (sb.v0 + sb.v1);
        let mut converged = false;
        for _ in 0..options.max_iterations {
            let cd = curve.derivatives(t, 1)?;
            let skl = surface.derivatives(u, v, 1)?;
            let c = cd[0];
            let s = skl[0][0];
            let f = NaVector3::new(c.x - s.x, c.y - s.y, c.z - s.z);
            if f.norm() < options.tolerance {
                converged = true;
                break;
            }
            // J = [C'(t), -S_u, -S_v] (columns).
            let ct = cd[1];
            let su = skl[1][0];
            let sv = skl[0][1];
            let j = Matrix3::new(
                ct.x, -su.x, -sv.x, //
                ct.y, -su.y, -sv.y, //
                ct.z, -su.z, -sv.z,
            );
            if j.determinant().abs() < 1e-14 {
                // Singular Jacobian: curve direction lies in the tangent plane
                // (grazing/tangential) — abandon this seed.
                break;
            }
            let Some(jinv) = j.try_inverse() else {
                break;
            };
            let delta = jinv * f;
            t = (t - delta.x).clamp(t_min, t_max);
            u = (u - delta.y).clamp(u_min, u_max);
            v = (v - delta.z).clamp(v_min, v_max);
            if delta.norm() < options.tolerance {
                let pc = curve.point_at(t)?;
                let ps = surface.point_at(u, v)?;
                converged = (pc - ps).norm() < options.tolerance.max(1e-7);
                break;
            }
        }
        if !converged {
            continue;
        }
        let pc = curve.point_at(t)?;
        let ps = surface.point_at(u, v)?;
        if (pc - ps).norm() > options.tolerance.max(1e-7) {
            continue;
        }
        let point = Point3::from((pc.coords + ps.coords) * 0.5);
        push_unique(&mut results, CurveSurfaceIntersection { point, t, u, v });
    }
    Ok(results)
}

/// A leaf-extent heuristic: 1/16 of the larger control-hull diagonal, mirroring
/// the curve×curve seeder. Coarse enough to keep the candidate count bounded yet
/// fine enough to separate distinct crossings.
fn seed_leaf_extent(c: &NurbsCurve3D, s: &NurbsSurface) -> f64 {
    let (c_lo, c_hi) = c.bounding_box();
    let (s_lo, s_hi) = s.bounding_box();
    let dc = (c_hi - c_lo).norm();
    let ds = (s_hi - s_lo).norm();
    (dc.max(ds) / 16.0).max(1e-6)
}

/// Inserts a result unless an existing one coincides in parameter space (tight)
/// or 3D position (loose, to collapse the symmetric near-solutions a grazing
/// contact produces).
fn push_unique(results: &mut Vec<CurveSurfaceIntersection>, candidate: CurveSurfaceIntersection) {
    const PARAM_DEDUP: f64 = 1e-6;
    const GEOM_DEDUP: f64 = 1e-4;
    for r in results.iter() {
        if (r.t - candidate.t).abs() < PARAM_DEDUP
            && (r.u - candidate.u).abs() < PARAM_DEDUP
            && (r.v - candidate.v).abs() < PARAM_DEDUP
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
    use crate::geometry::nurbs::{KnotVector, NurbsCurve3D, NurbsSurface};
    use crate::math::Point3;

    /// 2x2 bilinear patch spanning [0,2]x[0,2] in the z=0 plane.
    fn bilinear_patch() -> NurbsSurface {
        NurbsSurface::from_unweighted(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 2.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
                Point3::new(2.0, 2.0, 0.0),
            ],
            2,
            2,
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
            1,
        )
        .unwrap()
    }

    /// Quarter-cylinder shell: rational quadratic quarter circle in u (XY plane,
    /// radius 1) extruded along +Z by 2 in v.
    fn quarter_cylinder_patch() -> NurbsSurface {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        NurbsSurface::new(
            vec![
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 2.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(1.0, 1.0, 2.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 2.0),
            ],
            vec![1.0, 1.0, w, w, 1.0, 1.0],
            3,
            2,
            KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]).unwrap(),
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            2,
            1,
        )
        .unwrap()
    }

    fn line(p0: Point3, p1: Point3) -> NurbsCurve3D {
        NurbsCurve3D::from_unweighted(
            vec![p0, p1],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
        )
        .unwrap()
    }

    #[test]
    fn line_through_bilinear_patch_hits_known_point() {
        // Vertical line through the patch center (1, 1) crosses z=0 at (1, 1, 0).
        let s = bilinear_patch();
        let c = line(Point3::new(1.0, 1.0, -1.0), Point3::new(1.0, 1.0, 1.0));
        let hits = intersect_curve_surface(&c, &s, &IntersectionOptions::default()).unwrap();
        assert_eq!(hits.len(), 1, "expected one crossing, got {}", hits.len());
        assert!(
            (hits[0].point - Point3::new(1.0, 1.0, 0.0)).norm() < 1e-9,
            "hit {:?} not at (1,1,0)",
            hits[0].point
        );
        // Surface params: (1,1) is the patch center -> (u, v) = (0.5, 0.5).
        assert!((hits[0].u - 0.5).abs() < 1e-9);
        assert!((hits[0].v - 0.5).abs() < 1e-9);
        // Curve param: midpoint of z in [-1, 1] -> t = 0.5.
        assert!((hits[0].t - 0.5).abs() < 1e-9);
    }

    #[test]
    fn line_off_center_through_bilinear_patch() {
        let s = bilinear_patch();
        let c = line(Point3::new(0.5, 1.5, -3.0), Point3::new(0.5, 1.5, 2.0));
        let hits = intersect_curve_surface(&c, &s, &IntersectionOptions::default()).unwrap();
        assert_eq!(hits.len(), 1);
        assert!((hits[0].point - Point3::new(0.5, 1.5, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn line_through_quarter_cylinder_satisfies_cylinder_equation() {
        // A horizontal line at z = 1 (in the v-domain) pointed radially inward.
        // It must hit the cylinder x^2 + y^2 = 1.
        let s = quarter_cylinder_patch();
        // Line from outside, crossing the shell. Direction toward origin.
        let c = line(Point3::new(2.0, 0.6, 1.0), Point3::new(-2.0, 0.6, 1.0));
        let hits = intersect_curve_surface(&c, &s, &IntersectionOptions::default()).unwrap();
        assert!(!hits.is_empty(), "expected at least one hit");
        for h in &hits {
            let radial = (h.point.x * h.point.x + h.point.y * h.point.y).sqrt();
            assert!((radial - 1.0).abs() < 1e-8, "off cylinder: {:?}", h.point);
            assert!(
                (h.point.z - 1.0).abs() < 1e-8,
                "off z=1 plane: {:?}",
                h.point
            );
        }
    }

    #[test]
    fn parallel_line_missing_patch_yields_empty() {
        // A line in the z = 1 plane, parallel to the z=0 bilinear patch and
        // offset above it: no crossing.
        let s = bilinear_patch();
        let c = line(Point3::new(-1.0, 1.0, 1.0), Point3::new(3.0, 1.0, 1.0));
        let hits = intersect_curve_surface(&c, &s, &IntersectionOptions::default()).unwrap();
        assert!(hits.is_empty(), "expected no hits, got {:?}", hits);
    }

    #[test]
    fn line_outside_patch_extent_yields_empty() {
        // Vertical line crossing z=0 but at (5, 5) — far outside the [0,2]^2 patch.
        let s = bilinear_patch();
        let c = line(Point3::new(5.0, 5.0, -1.0), Point3::new(5.0, 5.0, 1.0));
        let hits = intersect_curve_surface(&c, &s, &IntersectionOptions::default()).unwrap();
        assert!(hits.is_empty(), "expected no hits, got {:?}", hits);
    }
}
