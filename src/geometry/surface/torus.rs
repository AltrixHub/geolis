use crate::error::{GeometryError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};

use super::{Surface, SurfaceDomain};

/// A toroidal surface in 3D space.
///
/// Defined by a center, major radius (center to tube center), minor radius
/// (tube radius), axis of symmetry, and a reference direction for u=0.
///
/// `P(u, v) = center + (R + r*cos(v)) * (cos(u)*ref_dir + sin(u)*binormal) + r*sin(v)*axis`
/// where `binormal = axis x ref_dir`.
///
/// Parameters: `u, v` in `[0, 2*pi)`.
#[derive(Debug, Clone)]
pub struct Torus {
    center: Point3,
    major_radius: f64,
    minor_radius: f64,
    axis: Vector3,
    ref_dir: Vector3,
}

impl Torus {
    /// Creates a new torus.
    ///
    /// # Arguments
    ///
    /// * `center` - Center of the torus
    /// * `major_radius` - Distance from center to tube center (must be positive)
    /// * `minor_radius` - Tube radius (must be positive, must be less than major radius)
    /// * `axis` - Symmetry axis direction (will be normalized)
    /// * `ref_dir` - Reference direction for u=0 (must be perpendicular to axis)
    ///
    /// # Errors
    ///
    /// Returns an error if either radius is non-positive, minor >= major,
    /// axis is zero-length, or the reference direction is not perpendicular to the axis.
    pub fn new(
        center: Point3,
        major_radius: f64,
        minor_radius: f64,
        axis: Vector3,
        ref_dir: Vector3,
    ) -> Result<Self> {
        if major_radius < TOLERANCE {
            return Err(
                GeometryError::Degenerate("torus major radius must be positive".into()).into(),
            );
        }
        if minor_radius < TOLERANCE {
            return Err(
                GeometryError::Degenerate("torus minor radius must be positive".into()).into(),
            );
        }
        if minor_radius >= major_radius {
            return Err(GeometryError::Degenerate(
                "torus minor radius must be less than major radius".into(),
            )
            .into());
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
            major_radius,
            minor_radius,
            axis,
            ref_dir,
        })
    }

    /// Returns the center of the torus.
    #[must_use]
    pub fn center(&self) -> &Point3 {
        &self.center
    }

    /// Returns the major radius (center to tube center).
    #[must_use]
    pub fn major_radius(&self) -> f64 {
        self.major_radius
    }

    /// Returns the minor radius (tube radius).
    #[must_use]
    pub fn minor_radius(&self) -> f64 {
        self.minor_radius
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
    /// - `u` = angle around the major axis in `(-pi, pi]`
    /// - `v` = angle around the tube cross-section in `(-pi, pi]`
    #[must_use]
    pub fn inverse(&self, point: &Point3) -> (f64, f64) {
        let dp = point - self.center;
        let binormal = self.binormal();
        // u: angle of the tube center around the major axis
        let u = dp.dot(&binormal).atan2(dp.dot(&self.ref_dir));
        // Compute tube center position
        let radial = self.ref_dir * u.cos() + binormal * u.sin();
        let tube_center = self.center + radial * self.major_radius;
        let to_tube = point - tube_center;
        // v: angle around the tube cross-section
        let v = to_tube.dot(&self.axis).atan2(to_tube.dot(&radial));
        (u, v)
    }
}

impl Surface for Torus {
    fn evaluate(&self, u: f64, v: f64) -> Result<Point3> {
        let binormal = self.binormal();
        let cu = u.cos();
        let su = u.sin();
        let cv = v.cos();
        let sv = v.sin();
        let radial = self.ref_dir * cu + binormal * su;
        let r = self.major_radius + self.minor_radius * cv;
        Ok(self.center + radial * r + self.axis * (self.minor_radius * sv))
    }

    fn normal(&self, u: f64, v: f64) -> Result<Vector3> {
        let binormal = self.binormal();
        let cu = u.cos();
        let su = u.sin();
        let cv = v.cos();
        let sv = v.sin();
        let radial = self.ref_dir * cu + binormal * su;
        let n = radial * cv + self.axis * sv;
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
            0.0,
            std::f64::consts::TAU,
        )
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, TAU};

    fn xy_torus() -> Torus {
        Torus::new(
            Point3::origin(),
            3.0,
            1.0,
            Vector3::z(),
            Vector3::x(),
        )
        .unwrap()
    }

    #[test]
    fn evaluate_outer_equator() {
        let t = xy_torus();
        // u=0, v=0: outer point at (R+r, 0, 0) = (4, 0, 0)
        let p = t.evaluate(0.0, 0.0).unwrap();
        assert!((p - Point3::new(4.0, 0.0, 0.0)).norm() < TOLERANCE);
    }

    #[test]
    fn evaluate_inner_equator() {
        let t = xy_torus();
        // u=0, v=pi: inner point at (R-r, 0, 0) = (2, 0, 0)
        let p = t.evaluate(0.0, std::f64::consts::PI).unwrap();
        assert!((p - Point3::new(2.0, 0.0, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn evaluate_top() {
        let t = xy_torus();
        // u=0, v=pi/2: top point at (R, 0, r) = (3, 0, 1)
        let p = t.evaluate(0.0, FRAC_PI_2).unwrap();
        assert!((p - Point3::new(3.0, 0.0, 1.0)).norm() < 1e-9);
    }

    #[test]
    fn evaluate_rotated_u() {
        let t = xy_torus();
        // u=pi/2, v=0: outer point at (0, R+r, 0) = (0, 4, 0)
        let p = t.evaluate(FRAC_PI_2, 0.0).unwrap();
        assert!((p - Point3::new(0.0, 4.0, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn normal_outward_at_outer() {
        let t = xy_torus();
        let n = t.normal(0.0, 0.0).unwrap();
        // At (4,0,0), outward normal should be +X
        assert!((n - Vector3::x()).norm() < TOLERANCE);
    }

    #[test]
    fn normal_at_top() {
        let t = xy_torus();
        let n = t.normal(0.0, FRAC_PI_2).unwrap();
        // At (3,0,1), normal should be +Z
        assert!((n - Vector3::z()).norm() < 1e-9);
    }

    #[test]
    fn domain_full_torus() {
        let t = xy_torus();
        let d = t.domain();
        assert!((d.u_min).abs() < TOLERANCE);
        assert!((d.u_max - TAU).abs() < TOLERANCE);
        assert!((d.v_min).abs() < TOLERANCE);
        assert!((d.v_max - TAU).abs() < TOLERANCE);
    }

    #[test]
    fn invalid_major_radius() {
        let r = Torus::new(Point3::origin(), 0.0, 1.0, Vector3::z(), Vector3::x());
        assert!(r.is_err());
    }

    #[test]
    fn invalid_minor_radius() {
        let r = Torus::new(Point3::origin(), 3.0, 0.0, Vector3::z(), Vector3::x());
        assert!(r.is_err());
    }

    #[test]
    fn minor_exceeds_major() {
        let r = Torus::new(Point3::origin(), 1.0, 2.0, Vector3::z(), Vector3::x());
        assert!(r.is_err());
    }

    #[test]
    fn inverse_roundtrip() {
        let t = xy_torus();
        for &(u, v) in &[
            (0.0, 0.0),
            (FRAC_PI_2, 0.0),
            (1.0, 0.5),
            (0.0, FRAC_PI_2),
            (TAU * 0.75, TAU * 0.25),
        ] {
            let p = t.evaluate(u, v).unwrap();
            let (u2, v2) = t.inverse(&p);
            let p2 = t.evaluate(u2, v2).unwrap();
            assert!((p - p2).norm() < 1e-9, "roundtrip failed for u={u}, v={v}");
        }
    }
}
