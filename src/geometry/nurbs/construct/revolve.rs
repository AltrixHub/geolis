use std::f64::consts::{FRAC_PI_2, TAU};

use crate::error::{GeometryError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};

use crate::geometry::nurbs::{KnotVector, NurbsCurve3D, NurbsSurface};

impl NurbsSurface {
    /// Builds a surface of revolution by sweeping a profile curve around an
    /// axis (The NURBS Book, A8.1).
    ///
    /// The v direction carries the exact rational-arc structure (same segment
    /// count, weights and knot pattern as [`NurbsCurve3D::arc`]): each profile
    /// control point is revolved into a circular arc around the axis through
    /// `axis_origin` along `axis_dir`, with weight products
    /// `w_profile * w_arc`. The u direction inherits the profile's degree and
    /// knot vector. `angle` is in `(0, 2*pi]`.
    ///
    /// # Errors
    ///
    /// Returns an error if `axis_dir` is zero-length, `angle` is not in
    /// `(0, 2*pi]`, the profile is empty, or the surface fails construction.
    // A8.1: segment-count and index conversions are exact small-integer casts;
    // single-char bindings (narcs, dtheta, w1, ...) follow The NURBS Book.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss,
        clippy::many_single_char_names
    )]
    pub fn revolve(
        profile: &NurbsCurve3D,
        axis_origin: Point3,
        axis_dir: Vector3,
        angle: f64,
    ) -> Result<Self> {
        let axis_len = axis_dir.norm();
        if axis_len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        let axis = axis_dir / axis_len;

        if !(TOLERANCE..=TAU + TOLERANCE).contains(&angle) {
            return Err(GeometryError::Degenerate(format!(
                "revolve angle {angle} must be in (0, 2*pi]"
            ))
            .into());
        }

        if profile.control_points().is_empty() {
            return Err(GeometryError::Degenerate("profile has no control points".into()).into());
        }

        // Rational-arc structure in v (matches A7.1 / NurbsCurve3D::arc).
        let narcs = ((angle / FRAC_PI_2).ceil() as usize).clamp(1, 4);
        let dtheta = angle / narcs as f64;
        let w1 = (dtheta / 2.0).cos();
        let nv = 2 * narcs + 1;

        // Per-segment v-weight multipliers: 1 at the endpoints, w1 at the
        // interior (odd) control points of each quadratic segment.
        let mut arc_weights = vec![1.0; nv];
        for k in 0..narcs {
            arc_weights[2 * k + 1] = w1;
        }

        let nu = profile.control_points().len();

        let mut control_points = Vec::with_capacity(nu * nv);
        let mut weights = Vec::with_capacity(nu * nv);

        for (point, &wp) in profile.control_points().iter().zip(profile.weights()) {
            // Decompose the profile point into its axial component and the
            // radial vector from the axis.
            let rel = *point - axis_origin;
            let axial = axis.dot(&rel);
            let foot = axis_origin + axis * axial;
            let radial = *point - foot;
            let r = radial.norm();

            if r < TOLERANCE {
                // Point lies on the axis: it is fixed under revolution. Every
                // v control point coincides with it; weights still follow the
                // arc pattern so the rational structure stays consistent.
                for &aw in &arc_weights {
                    control_points.push(*point);
                    weights.push(wp * aw);
                }
                continue;
            }

            // In-plane orthonormal frame for this point's circle.
            let x_axis = radial / r;
            let y_axis = axis.cross(&x_axis);

            let circle_point =
                |theta: f64| -> Point3 { foot + (x_axis * theta.cos() + y_axis * theta.sin()) * r };
            let circle_tangent =
                |theta: f64| -> Vector3 { x_axis * (-theta.sin()) + y_axis * theta.cos() };

            // Endpoint control points P0 (theta=0) and the per-segment mid
            // control points from tangent intersection (A7.1).
            let mut p0 = circle_point(0.0);
            let mut t0 = circle_tangent(0.0);
            control_points.push(p0);
            weights.push(wp);
            for s in 1..=narcs {
                let theta = dtheta * s as f64;
                let p2 = circle_point(theta);
                let t2 = circle_tangent(theta);
                let p1 = intersect_circle_tangents(&p0, &t0, &p2, &t2, &foot, &x_axis, &y_axis)?;
                control_points.push(p1);
                weights.push(wp * w1);
                control_points.push(p2);
                weights.push(wp);
                p0 = p2;
                t0 = t2;
            }
        }

        let knots_v = arc_knot_vector(narcs)?;

        NurbsSurface::new(
            control_points,
            weights,
            nu,
            nv,
            profile.knots().clone(),
            knots_v,
            profile.degree(),
            2,
        )
    }
}

/// Knot vector for the rational-arc v direction (A7.1 segment pattern).
fn arc_knot_vector(narcs: usize) -> Result<KnotVector> {
    let mut knots = vec![0.0; 3];
    match narcs {
        2 => knots.extend_from_slice(&[0.5, 0.5]),
        3 => knots.extend_from_slice(&[1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0]),
        4 => knots.extend_from_slice(&[0.25, 0.25, 0.5, 0.5, 0.75, 0.75]),
        _ => {}
    }
    knots.extend_from_slice(&[1.0; 3]);
    KnotVector::new(knots)
}

/// Intersects the two circle tangent lines for one quadratic arc segment,
/// projected into the plane spanned by `(x_axis, y_axis)` through `foot`.
#[allow(clippy::many_single_char_names)]
fn intersect_circle_tangents(
    p0: &Point3,
    t0: &Vector3,
    p2: &Point3,
    t2: &Vector3,
    foot: &Point3,
    x_axis: &Vector3,
    y_axis: &Vector3,
) -> Result<Point3> {
    let to_2d = |p: &Point3| -> (f64, f64) {
        let d = p - foot;
        (d.dot(x_axis), d.dot(y_axis))
    };
    let dir_2d = |v: &Vector3| -> (f64, f64) { (v.dot(x_axis), v.dot(y_axis)) };

    let (ax, ay) = to_2d(p0);
    let (bx, by) = to_2d(p2);
    let (ux, uy) = dir_2d(t0);
    let (vx, vy) = dir_2d(t2);

    // dtheta <= pi/2 guarantees the segment tangents are never parallel; the
    // determinant guard is defensive.
    let det = ux * vy - uy * vx;
    if det.abs() < TOLERANCE {
        return Err(GeometryError::Degenerate("parallel arc tangents".into()).into());
    }
    let s = ((bx - ax) * vy - (by - ay) * vx) / det;
    let ix = ax + s * ux;
    let iy = ay + s * uy;
    Ok(foot + *x_axis * ix + *y_axis * iy)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    /// Distance from `p` to the line through `origin` along unit `dir`.
    fn distance_to_axis(p: &Point3, origin: &Point3, dir: &Vector3) -> f64 {
        let rel = p - origin;
        let axial = dir.dot(&rel);
        (rel - dir * axial).norm()
    }

    /// Axial coordinate of `p` along the axis.
    fn axial_coord(p: &Point3, origin: &Point3, dir: &Vector3) -> f64 {
        dir.dot(&(p - origin))
    }

    fn vertical_line(x: f64) -> NurbsCurve3D {
        // Degree-1 line parallel to the Z axis at radius x, from z=0 to z=2.
        NurbsCurve3D::from_unweighted(
            vec![Point3::new(x, 0.0, 0.0), Point3::new(x, 0.0, 2.0)],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
        )
        .unwrap()
    }

    #[test]
    fn full_revolve_of_parallel_line_is_cylinder() {
        let r = 1.5;
        let profile = vertical_line(r);
        let s = NurbsSurface::revolve(&profile, Point3::origin(), Vector3::z(), TAU).unwrap();
        let ((u0, u1), (v0, v1)) = s.parameter_domain();
        for i in 0..=8 {
            let u = u0 + (u1 - u0) * f64::from(i) / 8.0;
            for j in 0..=16 {
                let v = v0 + (v1 - v0) * f64::from(j) / 16.0;
                let p = s.point_at(u, v).unwrap();
                let d = distance_to_axis(&p, &Point3::origin(), &Vector3::z());
                assert!((d - r).abs() < 1e-12, "radius off at ({u},{v}): {d}");
            }
        }
    }

    #[test]
    fn full_revolve_of_tilted_line_is_cone() {
        // Profile from (1, 0, 0) to (2, 0, 2): radius grows linearly with z.
        let profile = NurbsCurve3D::from_unweighted(
            vec![Point3::new(1.0, 0.0, 0.0), Point3::new(2.0, 0.0, 2.0)],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
        )
        .unwrap();
        let s = NurbsSurface::revolve(&profile, Point3::origin(), Vector3::z(), TAU).unwrap();
        let ((u0, u1), (v0, v1)) = s.parameter_domain();
        for i in 0..=8 {
            let u = u0 + (u1 - u0) * f64::from(i) / 8.0;
            for j in 0..=16 {
                let v = v0 + (v1 - v0) * f64::from(j) / 16.0;
                let p = s.point_at(u, v).unwrap();
                let h = axial_coord(&p, &Point3::origin(), &Vector3::z());
                let d = distance_to_axis(&p, &Point3::origin(), &Vector3::z());
                // radius = 1 + h/2 (linear interpolation along the profile).
                let expected = 1.0 + h / 2.0;
                assert!((d - expected).abs() < 1e-12, "cone radius off at ({u},{v})");
            }
        }
    }

    #[test]
    fn quarter_revolve_boundary_isocurves_at_correct_angles() {
        let r = 1.0;
        let profile = vertical_line(r);
        let s = NurbsSurface::revolve(&profile, Point3::origin(), Vector3::z(), FRAC_PI_2).unwrap();
        let ((_, _), (v0, v1)) = s.parameter_domain();
        // At v_min the swept point sits at angle 0 (along +X); at v_max at
        // angle pi/2 (along +Y).
        let start = s.point_at(0.0, v0).unwrap();
        let end = s.point_at(0.0, v1).unwrap();
        assert!((start - Point3::new(r, 0.0, 0.0)).norm() < 1e-12);
        assert!((end - Point3::new(0.0, r, 0.0)).norm() < 1e-12);
    }

    #[test]
    fn full_revolve_closes() {
        let profile = vertical_line(1.0);
        let s = NurbsSurface::revolve(&profile, Point3::origin(), Vector3::z(), TAU).unwrap();
        let ((u0, u1), (v0, v1)) = s.parameter_domain();
        for i in 0..=8 {
            let u = u0 + (u1 - u0) * f64::from(i) / 8.0;
            let a = s.point_at(u, v0).unwrap();
            let b = s.point_at(u, v1).unwrap();
            assert!((a - b).norm() < 1e-12, "not closed at u={u}");
        }
    }

    #[test]
    fn rejects_zero_axis() {
        let profile = vertical_line(1.0);
        assert!(NurbsSurface::revolve(&profile, Point3::origin(), Vector3::zeros(), PI).is_err());
    }

    #[test]
    fn rejects_zero_angle() {
        let profile = vertical_line(1.0);
        assert!(NurbsSurface::revolve(&profile, Point3::origin(), Vector3::z(), 0.0).is_err());
    }
}
