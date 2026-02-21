use crate::error::{GeometryError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};

use super::{Surface, SurfaceDomain};

/// A spherical surface in 3D space.
///
/// Defined by a center, radius, axis (north pole direction), and a
/// reference direction for the equator at u=0.
///
/// `P(u, v) = center + r * cos(v) * cos(u) * ref_dir + r * cos(v) * sin(u) * binormal + r * sin(v) * axis`
/// where `binormal = axis x ref_dir`.
///
/// Parameters: `u` = longitude `[0, 2*pi)`, `v` = latitude `[-pi/2, pi/2]`.
/// The outward normal is `(P - center) / radius`.
#[derive(Debug, Clone)]
pub struct Sphere {
    center: Point3,
    radius: f64,
    axis: Vector3,
    ref_dir: Vector3,
}

impl Sphere {
    /// Creates a new sphere.
    ///
    /// # Arguments
    ///
    /// * `center` - Center of the sphere
    /// * `radius` - Radius (must be positive)
    /// * `axis` - North pole direction (will be normalized)
    /// * `ref_dir` - Equatorial reference direction for u=0 (must be perpendicular to axis)
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
                GeometryError::Degenerate("sphere radius must be positive".into()).into(),
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

    /// Returns the center of the sphere.
    #[must_use]
    pub fn center(&self) -> &Point3 {
        &self.center
    }

    /// Returns the radius.
    #[must_use]
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Returns the axis direction (north pole, unit vector).
    #[must_use]
    pub fn axis(&self) -> &Vector3 {
        &self.axis
    }

    /// Returns the reference direction (u=0 on equator).
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
    /// - `u` = longitude in `(-pi, pi]`
    /// - `v` = latitude in `[-pi/2, pi/2]`
    #[must_use]
    pub fn inverse(&self, point: &Point3) -> (f64, f64) {
        let dp = (point - self.center) / self.radius;
        let binormal = self.binormal();
        let v = dp.dot(&self.axis).clamp(-1.0, 1.0).asin();
        let u = dp.dot(&binormal).atan2(dp.dot(&self.ref_dir));
        (u, v)
    }
}

impl Surface for Sphere {
    fn evaluate(&self, u: f64, v: f64) -> Result<Point3> {
        let binormal = self.binormal();
        let cv = v.cos();
        let sv = v.sin();
        let cu = u.cos();
        let su = u.sin();
        Ok(self.center
            + self.ref_dir * (self.radius * cv * cu)
            + binormal * (self.radius * cv * su)
            + self.axis * (self.radius * sv))
    }

    fn normal(&self, u: f64, v: f64) -> Result<Vector3> {
        let binormal = self.binormal();
        let cv = v.cos();
        let sv = v.sin();
        let cu = u.cos();
        let su = u.sin();
        let n = self.ref_dir * (cv * cu) + binormal * (cv * su) + self.axis * sv;
        let len = n.norm();
        if len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        Ok(n / len)
    }

    fn domain(&self) -> SurfaceDomain {
        SurfaceDomain::new(
            0.0,
            std::f64::consts::TAU,
            -std::f64::consts::FRAC_PI_2,
            std::f64::consts::FRAC_PI_2,
        )
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, TAU};

    fn unit_sphere() -> Sphere {
        Sphere::new(
            Point3::origin(),
            1.0,
            Vector3::z(),
            Vector3::x(),
        )
        .unwrap()
    }

    #[test]
    fn evaluate_equator_at_zero() {
        let s = unit_sphere();
        let p = s.evaluate(0.0, 0.0).unwrap();
        assert!((p - Point3::new(1.0, 0.0, 0.0)).norm() < TOLERANCE);
    }

    #[test]
    fn evaluate_equator_at_pi_over_2() {
        let s = unit_sphere();
        let p = s.evaluate(FRAC_PI_2, 0.0).unwrap();
        assert!((p - Point3::new(0.0, 1.0, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn evaluate_north_pole() {
        let s = unit_sphere();
        let p = s.evaluate(0.0, FRAC_PI_2).unwrap();
        assert!((p - Point3::new(0.0, 0.0, 1.0)).norm() < 1e-9);
    }

    #[test]
    fn evaluate_south_pole() {
        let s = unit_sphere();
        let p = s.evaluate(0.0, -FRAC_PI_2).unwrap();
        assert!((p - Point3::new(0.0, 0.0, -1.0)).norm() < 1e-9);
    }

    #[test]
    fn normal_outward_equator() {
        let s = unit_sphere();
        let n = s.normal(0.0, 0.0).unwrap();
        assert!((n - Vector3::x()).norm() < TOLERANCE);
    }

    #[test]
    fn normal_at_north_pole() {
        let s = unit_sphere();
        let n = s.normal(0.0, FRAC_PI_2).unwrap();
        assert!((n - Vector3::z()).norm() < 1e-9);
    }

    #[test]
    fn normal_at_south_pole() {
        let s = unit_sphere();
        let n = s.normal(0.0, -FRAC_PI_2).unwrap();
        assert!((n - Vector3::new(0.0, 0.0, -1.0)).norm() < 1e-9);
    }

    #[test]
    fn domain_ranges() {
        let s = unit_sphere();
        let d = s.domain();
        assert!((d.u_min).abs() < TOLERANCE);
        assert!((d.u_max - TAU).abs() < TOLERANCE);
        assert!((d.v_min + FRAC_PI_2).abs() < TOLERANCE);
        assert!((d.v_max - FRAC_PI_2).abs() < TOLERANCE);
    }

    #[test]
    fn offset_center() {
        let s = Sphere::new(
            Point3::new(1.0, 2.0, 3.0),
            2.0,
            Vector3::z(),
            Vector3::x(),
        )
        .unwrap();
        let p = s.evaluate(0.0, 0.0).unwrap();
        assert!((p - Point3::new(3.0, 2.0, 3.0)).norm() < TOLERANCE);
    }

    #[test]
    fn invalid_radius() {
        let r = Sphere::new(Point3::origin(), 0.0, Vector3::z(), Vector3::x());
        assert!(r.is_err());
    }

    #[test]
    fn inverse_roundtrip() {
        let s = unit_sphere();
        for &(u, v) in &[
            (0.0, 0.0),
            (FRAC_PI_2, 0.0),
            (1.0, 0.5),
            (0.0, FRAC_PI_2 * 0.9),
            (0.0, -FRAC_PI_2 * 0.9),
        ] {
            let p = s.evaluate(u, v).unwrap();
            let (u2, v2) = s.inverse(&p);
            let p2 = s.evaluate(u2, v2).unwrap();
            assert!((p - p2).norm() < 1e-9, "roundtrip failed for u={u}, v={v}");
        }
    }
}
