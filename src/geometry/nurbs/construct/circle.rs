use std::f64::consts::{FRAC_PI_2, TAU};

use crate::error::{GeometryError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};

use crate::geometry::nurbs::{KnotVector, NurbsCurve3D};

/// Builds an exact rational quadratic NURBS arc (The NURBS Book, A7.1).
///
/// The arc lies in the plane spanned by `x_axis` and `y_axis` through
/// `center`, sweeping from `start_angle` to `end_angle` (radians, measured
/// from `x_axis` toward `y_axis`). The sweep must be in `(0, 2*pi]`.
///
/// # Errors
///
/// Returns an error if the radius is non-positive, the axes are degenerate
/// or not perpendicular, or the sweep is not in `(0, 2*pi]`.
pub fn nurbs_arc(
    center: Point3,
    x_axis: Vector3,
    y_axis: Vector3,
    radius: f64,
    start_angle: f64,
    end_angle: f64,
) -> Result<NurbsCurve3D> {
    if radius < TOLERANCE {
        return Err(GeometryError::Degenerate("arc radius must be positive".into()).into());
    }
    let x_len = x_axis.norm();
    let y_len = y_axis.norm();
    if x_len < TOLERANCE || y_len < TOLERANCE {
        return Err(GeometryError::ZeroVector.into());
    }
    let x_axis = x_axis / x_len;
    let y_axis = y_axis / y_len;
    if x_axis.dot(&y_axis).abs() > TOLERANCE {
        return Err(GeometryError::Degenerate("arc axes must be perpendicular".into()).into());
    }
    let theta = end_angle - start_angle;
    if theta < TOLERANCE || theta > TAU + TOLERANCE {
        return Err(
            GeometryError::Degenerate(format!("arc sweep {theta} must be in (0, 2*pi]")).into(),
        );
    }

    // Number of quadratic segments: one per quarter turn.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let narcs = ((theta / FRAC_PI_2).ceil() as usize).clamp(1, 4);
    let dtheta = theta / narcs as f64;
    let w1 = (dtheta / 2.0).cos();

    let point_at_angle =
        |angle: f64| -> Point3 { center + (x_axis * angle.cos() + y_axis * angle.sin()) * radius };
    let tangent_at_angle =
        |angle: f64| -> Vector3 { x_axis * (-angle.sin()) + y_axis * angle.cos() };

    let n = 2 * narcs;
    let mut points = vec![Point3::origin(); n + 1];
    let mut weights = vec![1.0; n + 1];

    let mut angle = start_angle;
    points[0] = point_at_angle(angle);
    let mut p0 = points[0];
    let mut t0 = tangent_at_angle(angle);
    for i in 1..=narcs {
        angle = start_angle + dtheta * i as f64;
        let p2 = point_at_angle(angle);
        let t2 = tangent_at_angle(angle);
        // Mid control point: intersection of the two tangent lines.
        let p1 = intersect_tangents(&p0, &t0, &p2, &t2, &x_axis, &y_axis, &center)?;
        points[2 * i - 1] = p1;
        weights[2 * i - 1] = w1;
        points[2 * i] = p2;
        p0 = p2;
        t0 = t2;
    }

    // Knot vector by segment count (A7.1)
    let mut knots = vec![0.0; 3];
    match narcs {
        1 => {}
        2 => knots.extend_from_slice(&[0.5, 0.5]),
        3 => knots.extend_from_slice(&[1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0]),
        _ => knots.extend_from_slice(&[0.25, 0.25, 0.5, 0.5, 0.75, 0.75]),
    }
    knots.extend_from_slice(&[1.0; 3]);

    NurbsCurve3D::new(points, weights, KnotVector::new(knots)?, 2)
}

/// Builds an exact full NURBS circle.
///
/// # Errors
///
/// Same validation as [`nurbs_arc`].
pub fn nurbs_circle(
    center: Point3,
    x_axis: Vector3,
    y_axis: Vector3,
    radius: f64,
) -> Result<NurbsCurve3D> {
    nurbs_arc(center, x_axis, y_axis, radius, 0.0, TAU)
}

/// Intersects two coplanar tangent lines in the plane spanned by
/// `(x_axis, y_axis)` through `center`.
#[allow(clippy::too_many_arguments)]
fn intersect_tangents(
    p0: &Point3,
    t0: &Vector3,
    p2: &Point3,
    t2: &Vector3,
    x_axis: &Vector3,
    y_axis: &Vector3,
    center: &Point3,
) -> Result<Point3> {
    // Project to 2D plane coordinates
    let to_2d = |p: &Point3| -> (f64, f64) {
        let d = p - center;
        (d.dot(x_axis), d.dot(y_axis))
    };
    let dir_2d = |v: &Vector3| -> (f64, f64) { (v.dot(x_axis), v.dot(y_axis)) };

    let (ax, ay) = to_2d(p0);
    let (bx, by) = to_2d(p2);
    let (ux, uy) = dir_2d(t0);
    let (vx, vy) = dir_2d(t2);

    let det = ux * vy - uy * vx;
    if det.abs() < TOLERANCE {
        return Err(GeometryError::Degenerate("parallel arc tangents".into()).into());
    }
    let s = ((bx - ax) * vy - (by - ay) * vx) / det;
    let ix = ax + s * ux;
    let iy = ay + s * uy;
    Ok(center + *x_axis * ix + *y_axis * iy)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::TOLERANCE;
    use std::f64::consts::{FRAC_PI_2, PI, TAU};

    #[test]
    fn full_circle_is_exact() {
        let c = nurbs_circle(Point3::origin(), Vector3::x(), Vector3::y(), 2.0).unwrap();
        let (t0, t1) = c.parameter_domain();
        for i in 0..=100 {
            let t = t0 + (t1 - t0) * f64::from(i) / 100.0;
            let p = c.point_at(t).unwrap();
            assert!((p.coords.norm() - 2.0).abs() < 1e-12, "radius off at t={t}");
            assert!(p.z.abs() < 1e-12, "out of plane at t={t}");
        }
        assert!(c.is_endpoint_closed());
    }

    #[test]
    fn quarter_arc_endpoints() {
        let c = nurbs_arc(
            Point3::origin(),
            Vector3::x(),
            Vector3::y(),
            1.0,
            0.0,
            FRAC_PI_2,
        )
        .unwrap();
        let (t0, t1) = c.parameter_domain();
        let start = c.point_at(t0).unwrap();
        let end = c.point_at(t1).unwrap();
        assert!((start - Point3::new(1.0, 0.0, 0.0)).norm() < TOLERANCE);
        assert!((end - Point3::new(0.0, 1.0, 0.0)).norm() < TOLERANCE);
    }

    #[test]
    fn three_quarter_arc_stays_on_circle() {
        let c = nurbs_arc(
            Point3::new(1.0, 2.0, 3.0),
            Vector3::x(),
            Vector3::z(),
            1.5,
            0.3,
            0.3 + 1.5 * PI,
        )
        .unwrap();
        let (t0, t1) = c.parameter_domain();
        for i in 0..=60 {
            let t = t0 + (t1 - t0) * f64::from(i) / 60.0;
            let p = c.point_at(t).unwrap();
            let radial = p - Point3::new(1.0, 2.0, 3.0);
            assert!((radial.norm() - 1.5).abs() < 1e-12, "radius off at t={t}");
            assert!(radial.y.abs() < 1e-12, "out of XZ plane at t={t}");
        }
    }

    #[test]
    fn rejects_zero_sweep() {
        let r = nurbs_arc(Point3::origin(), Vector3::x(), Vector3::y(), 1.0, 1.0, 1.0);
        assert!(r.is_err());
    }

    #[test]
    fn rejects_sweep_beyond_full_circle() {
        let r = nurbs_arc(
            Point3::origin(),
            Vector3::x(),
            Vector3::y(),
            1.0,
            0.0,
            TAU + 0.1,
        );
        assert!(r.is_err());
    }
}
