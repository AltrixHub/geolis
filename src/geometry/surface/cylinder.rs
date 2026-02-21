use crate::error::{GeometryError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};

use super::{Surface, SurfaceDomain};

/// A cylindrical surface in 3D space.
///
/// Defined by a center point on the axis, radius, axis direction, and
/// a reference direction for u=0.
///
/// `P(u, v) = center + radius * cos(u) * ref_dir + radius * sin(u) * binormal + v * axis`
/// where `binormal = axis x ref_dir`.
///
/// The outward normal is `cos(u) * ref_dir + sin(u) * binormal`.
#[derive(Debug, Clone)]
pub struct Cylinder {
    center: Point3,
    radius: f64,
    axis: Vector3,
    ref_dir: Vector3,
}

impl Cylinder {
    /// Creates a new cylinder.
    ///
    /// # Arguments
    ///
    /// * `center` - A point on the cylinder axis
    /// * `radius` - Radius (must be positive)
    /// * `axis` - Axis direction (will be normalized)
    /// * `ref_dir` - Reference direction for u=0 (must be perpendicular to axis)
    ///
    /// # Errors
    ///
    /// Returns an error if the radius is non-positive, axis is zero-length,
    /// or the reference direction is not perpendicular to the axis.
    pub fn new(
        center: Point3,
        radius: f64,
        axis: Vector3,
        ref_dir: Vector3,
    ) -> Result<Self> {
        if radius < TOLERANCE {
            return Err(
                GeometryError::Degenerate("cylinder radius must be positive".into()).into(),
            );
        }

        let axis_len = axis.norm();
        if axis_len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        let axis = axis / axis_len;

        let ref_len = ref_dir.norm();
        if ref_len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        let ref_dir = ref_dir / ref_len;

        if axis.dot(&ref_dir).abs() > TOLERANCE {
            return Err(GeometryError::Degenerate(
                "reference direction must be perpendicular to axis".into(),
            )
            .into());
        }

        Ok(Self {
            center,
            radius,
            axis,
            ref_dir,
        })
    }

    /// Returns the center point on the axis.
    #[must_use]
    pub fn center(&self) -> &Point3 {
        &self.center
    }

    /// Returns the radius.
    #[must_use]
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Returns the axis direction (unit vector).
    #[must_use]
    pub fn axis(&self) -> &Vector3 {
        &self.axis
    }

    /// Returns the reference direction (u=0).
    #[must_use]
    pub fn ref_dir(&self) -> &Vector3 {
        &self.ref_dir
    }

    /// Computes the binormal direction (`axis x ref_dir`).
    fn binormal(&self) -> Vector3 {
        self.axis.cross(&self.ref_dir)
    }

    /// Computes the (u, v) parameters for a given 3D point on the surface.
    ///
    /// - `u` = angle in `(-pi, pi]` (atan2-based)
    /// - `v` = signed distance along the axis from the center
    #[must_use]
    pub fn inverse(&self, point: &Point3) -> (f64, f64) {
        let dp = point - self.center;
        let v = dp.dot(&self.axis);
        let binormal = self.binormal();
        let u = dp.dot(&binormal).atan2(dp.dot(&self.ref_dir));
        (u, v)
    }
}

impl Surface for Cylinder {
    fn evaluate(&self, u: f64, v: f64) -> Result<Point3> {
        let binormal = self.binormal();
        let x = self.radius * u.cos();
        let y = self.radius * u.sin();
        Ok(self.center + self.ref_dir * x + binormal * y + self.axis * v)
    }

    fn normal(&self, u: f64, _v: f64) -> Result<Vector3> {
        let binormal = self.binormal();
        let n = self.ref_dir * u.cos() + binormal * u.sin();
        let len = n.norm();
        if len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        Ok(n / len)
    }

    fn domain(&self) -> SurfaceDomain {
        SurfaceDomain::new(0.0, std::f64::consts::TAU, f64::NEG_INFINITY, f64::INFINITY)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, TAU};

    fn z_cylinder(radius: f64) -> Cylinder {
        Cylinder::new(
            Point3::origin(),
            radius,
            Vector3::z(),
            Vector3::x(),
        )
        .unwrap()
    }

    #[test]
    fn evaluate_at_origin() {
        let c = z_cylinder(2.0);
        let p = c.evaluate(0.0, 0.0).unwrap();
        assert!((p - Point3::new(2.0, 0.0, 0.0)).norm() < TOLERANCE);
    }

    #[test]
    fn evaluate_at_pi_over_2() {
        let c = z_cylinder(2.0);
        let p = c.evaluate(FRAC_PI_2, 0.0).unwrap();
        assert!((p - Point3::new(0.0, 2.0, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn evaluate_with_height() {
        let c = z_cylinder(1.0);
        let p = c.evaluate(0.0, 5.0).unwrap();
        assert!((p - Point3::new(1.0, 0.0, 5.0)).norm() < TOLERANCE);
    }

    #[test]
    fn normal_outward_at_zero() {
        let c = z_cylinder(1.0);
        let n = c.normal(0.0, 0.0).unwrap();
        assert!((n - Vector3::x()).norm() < TOLERANCE);
    }

    #[test]
    fn normal_outward_at_pi_over_2() {
        let c = z_cylinder(1.0);
        let n = c.normal(FRAC_PI_2, 0.0).unwrap();
        assert!((n - Vector3::y()).norm() < 1e-9);
    }

    #[test]
    fn domain_u_is_full_circle() {
        let c = z_cylinder(1.0);
        let d = c.domain();
        assert!((d.u_min).abs() < TOLERANCE);
        assert!((d.u_max - TAU).abs() < TOLERANCE);
        assert!(d.v_min.is_infinite());
        assert!(d.v_max.is_infinite());
    }

    #[test]
    fn invalid_radius() {
        let r = Cylinder::new(Point3::origin(), 0.0, Vector3::z(), Vector3::x());
        assert!(r.is_err());
    }

    #[test]
    fn inverse_roundtrip() {
        let c = z_cylinder(2.0);
        for &(u, v) in &[(0.0, 0.0), (FRAC_PI_2, 3.0), (1.0, -2.5), (TAU * 0.75, 1.0)] {
            let p = c.evaluate(u, v).unwrap();
            let (u2, v2) = c.inverse(&p);
            let p2 = c.evaluate(u2, v2).unwrap();
            assert!((p - p2).norm() < 1e-9, "roundtrip failed for u={u}, v={v}");
        }
    }
}
