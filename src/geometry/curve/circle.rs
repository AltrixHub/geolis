use crate::error::{GeometryError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};

use super::{Curve, CurveDomain};

/// A full circle in 3D space.
///
/// Defined by a center, radius, normal axis, and a reference direction
/// for the zero-angle. The parametric domain is `[0, 2*pi)` and the
/// curve is always closed.
///
/// `P(t) = center + radius * cos(t) * ref_dir + radius * sin(t) * binormal`
/// where `binormal = normal x ref_dir`.
#[derive(Debug, Clone)]
pub struct Circle {
    center: Point3,
    radius: f64,
    normal: Vector3,
    ref_dir: Vector3,
}

impl Circle {
    /// Creates a new circle.
    ///
    /// # Arguments
    ///
    /// * `center` - Center of the circle
    /// * `radius` - Radius (must be positive)
    /// * `normal` - Normal vector defining the circle plane
    /// * `ref_dir` - Reference direction for angle = 0 (must be perpendicular to normal)
    ///
    /// # Errors
    ///
    /// Returns an error if the radius is non-positive, the normal is zero-length,
    /// or the reference direction is not perpendicular to the normal.
    pub fn new(
        center: Point3,
        radius: f64,
        normal: Vector3,
        ref_dir: Vector3,
    ) -> Result<Self> {
        if radius < TOLERANCE {
            return Err(
                GeometryError::Degenerate("circle radius must be positive".into()).into(),
            );
        }

        let normal_len = normal.norm();
        if normal_len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        let normal = normal / normal_len;

        let ref_len = ref_dir.norm();
        if ref_len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        let ref_dir = ref_dir / ref_len;

        if normal.dot(&ref_dir).abs() > TOLERANCE {
            return Err(GeometryError::Degenerate(
                "reference direction must be perpendicular to normal".into(),
            )
            .into());
        }

        Ok(Self {
            center,
            radius,
            normal,
            ref_dir,
        })
    }

    /// Returns the center of the circle.
    #[must_use]
    pub fn center(&self) -> &Point3 {
        &self.center
    }

    /// Returns the radius of the circle.
    #[must_use]
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Returns the normal vector of the circle plane.
    #[must_use]
    pub fn normal(&self) -> &Vector3 {
        &self.normal
    }

    /// Returns the reference direction (t=0 direction).
    #[must_use]
    pub fn ref_dir(&self) -> &Vector3 {
        &self.ref_dir
    }

    /// Computes the binormal direction (`normal x ref_dir`).
    fn binormal(&self) -> Vector3 {
        self.normal.cross(&self.ref_dir)
    }
}

impl Curve for Circle {
    fn evaluate(&self, t: f64) -> Result<Point3> {
        let binormal = self.binormal();
        let x = self.radius * t.cos();
        let y = self.radius * t.sin();
        Ok(self.center + self.ref_dir * x + binormal * y)
    }

    fn tangent(&self, t: f64) -> Result<Vector3> {
        let binormal = self.binormal();
        let dx = -self.radius * t.sin();
        let dy = self.radius * t.cos();
        let tangent = self.ref_dir * dx + binormal * dy;
        let len = tangent.norm();
        if len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        Ok(tangent / len)
    }

    fn domain(&self) -> CurveDomain {
        CurveDomain::new(0.0, std::f64::consts::TAU)
    }

    fn is_closed(&self) -> bool {
        true
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, TAU};

    fn xy_circle(radius: f64) -> Circle {
        Circle::new(
            Point3::origin(),
            radius,
            Vector3::z(),
            Vector3::x(),
        )
        .unwrap()
    }

    #[test]
    fn evaluate_at_zero() {
        let c = xy_circle(2.0);
        let p = c.evaluate(0.0).unwrap();
        assert!((p - Point3::new(2.0, 0.0, 0.0)).norm() < TOLERANCE);
    }

    #[test]
    fn evaluate_at_pi_over_2() {
        let c = xy_circle(3.0);
        let p = c.evaluate(FRAC_PI_2).unwrap();
        assert!((p - Point3::new(0.0, 3.0, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn tangent_at_zero() {
        let c = xy_circle(1.0);
        let t = c.tangent(0.0).unwrap();
        // At t=0, tangent should be +Y direction
        assert!((t - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn is_always_closed() {
        let c = xy_circle(1.0);
        assert!(c.is_closed());
    }

    #[test]
    fn domain_is_full_circle() {
        let c = xy_circle(1.0);
        let d = c.domain();
        assert!((d.t_min).abs() < TOLERANCE);
        assert!((d.t_max - TAU).abs() < TOLERANCE);
    }

    #[test]
    fn offset_center() {
        let c = Circle::new(
            Point3::new(1.0, 2.0, 3.0),
            1.0,
            Vector3::z(),
            Vector3::x(),
        )
        .unwrap();
        let p = c.evaluate(0.0).unwrap();
        assert!((p - Point3::new(2.0, 2.0, 3.0)).norm() < TOLERANCE);
    }

    #[test]
    fn invalid_radius() {
        let r = Circle::new(Point3::origin(), 0.0, Vector3::z(), Vector3::x());
        assert!(r.is_err());
    }

    #[test]
    fn non_perpendicular_ref_dir() {
        let r = Circle::new(
            Point3::origin(),
            1.0,
            Vector3::z(),
            Vector3::new(1.0, 0.0, 1.0),
        );
        assert!(r.is_err());
    }
}
