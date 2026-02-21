use crate::error::{GeometryError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};

use super::{Curve, CurveDomain};

/// An ellipse (or elliptical arc) in 3D space.
///
/// Defined by a center, semi-major and semi-minor axes, a normal,
/// a major axis direction, and angular range.
///
/// `P(t) = center + a * cos(t) * major_dir + b * sin(t) * minor_dir`
/// where `minor_dir = normal x major_dir`.
#[derive(Debug, Clone)]
pub struct Ellipse {
    center: Point3,
    semi_major: f64,
    semi_minor: f64,
    normal: Vector3,
    major_dir: Vector3,
    start_angle: f64,
    end_angle: f64,
}

impl Ellipse {
    /// Creates a new ellipse.
    ///
    /// # Arguments
    ///
    /// * `center` - Center of the ellipse
    /// * `semi_major` - Semi-major axis length (must be positive)
    /// * `semi_minor` - Semi-minor axis length (must be positive)
    /// * `normal` - Normal vector defining the ellipse plane
    /// * `major_dir` - Major axis direction (must be perpendicular to normal)
    /// * `start_angle` - Start angle in radians
    /// * `end_angle` - End angle in radians
    ///
    /// # Errors
    ///
    /// Returns an error if either axis length is non-positive, the normal is
    /// zero-length, or the major direction is not perpendicular to the normal.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        center: Point3,
        semi_major: f64,
        semi_minor: f64,
        normal: Vector3,
        major_dir: Vector3,
        start_angle: f64,
        end_angle: f64,
    ) -> Result<Self> {
        if semi_major < TOLERANCE {
            return Err(
                GeometryError::Degenerate("semi-major axis must be positive".into()).into(),
            );
        }
        if semi_minor < TOLERANCE {
            return Err(
                GeometryError::Degenerate("semi-minor axis must be positive".into()).into(),
            );
        }

        let normal_len = normal.norm();
        if normal_len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        let normal = normal / normal_len;

        let major_len = major_dir.norm();
        if major_len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        let major_dir = major_dir / major_len;

        if normal.dot(&major_dir).abs() > TOLERANCE {
            return Err(GeometryError::Degenerate(
                "major direction must be perpendicular to normal".into(),
            )
            .into());
        }

        Ok(Self {
            center,
            semi_major,
            semi_minor,
            normal,
            major_dir,
            start_angle,
            end_angle,
        })
    }

    /// Returns the center of the ellipse.
    #[must_use]
    pub fn center(&self) -> &Point3 {
        &self.center
    }

    /// Returns the semi-major axis length.
    #[must_use]
    pub fn semi_major(&self) -> f64 {
        self.semi_major
    }

    /// Returns the semi-minor axis length.
    #[must_use]
    pub fn semi_minor(&self) -> f64 {
        self.semi_minor
    }

    /// Returns the normal vector of the ellipse plane.
    #[must_use]
    pub fn normal(&self) -> &Vector3 {
        &self.normal
    }

    /// Returns the major axis direction.
    #[must_use]
    pub fn major_dir(&self) -> &Vector3 {
        &self.major_dir
    }

    /// Returns the start angle.
    #[must_use]
    pub fn start_angle(&self) -> f64 {
        self.start_angle
    }

    /// Returns the end angle.
    #[must_use]
    pub fn end_angle(&self) -> f64 {
        self.end_angle
    }

    /// Computes the minor axis direction (`normal x major_dir`).
    fn minor_dir(&self) -> Vector3 {
        self.normal.cross(&self.major_dir)
    }
}

impl Curve for Ellipse {
    fn evaluate(&self, t: f64) -> Result<Point3> {
        let minor = self.minor_dir();
        let x = self.semi_major * t.cos();
        let y = self.semi_minor * t.sin();
        Ok(self.center + self.major_dir * x + minor * y)
    }

    fn tangent(&self, t: f64) -> Result<Vector3> {
        let minor = self.minor_dir();
        let dx = -self.semi_major * t.sin();
        let dy = self.semi_minor * t.cos();
        let tangent = self.major_dir * dx + minor * dy;
        let len = tangent.norm();
        if len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        Ok(tangent / len)
    }

    fn domain(&self) -> CurveDomain {
        CurveDomain::new(self.start_angle, self.end_angle)
    }

    fn is_closed(&self) -> bool {
        (self.end_angle - self.start_angle - std::f64::consts::TAU).abs() < TOLERANCE
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, TAU};

    fn xy_ellipse(a: f64, b: f64) -> Ellipse {
        Ellipse::new(
            Point3::origin(),
            a,
            b,
            Vector3::z(),
            Vector3::x(),
            0.0,
            TAU,
        )
        .unwrap()
    }

    fn xy_ellipse_arc(a: f64, b: f64, start: f64, end: f64) -> Ellipse {
        Ellipse::new(
            Point3::origin(),
            a,
            b,
            Vector3::z(),
            Vector3::x(),
            start,
            end,
        )
        .unwrap()
    }

    #[test]
    fn evaluate_at_zero() {
        let e = xy_ellipse(3.0, 2.0);
        let p = e.evaluate(0.0).unwrap();
        assert!((p - Point3::new(3.0, 0.0, 0.0)).norm() < TOLERANCE);
    }

    #[test]
    fn evaluate_at_pi_over_2() {
        let e = xy_ellipse(3.0, 2.0);
        let p = e.evaluate(FRAC_PI_2).unwrap();
        assert!((p - Point3::new(0.0, 2.0, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn tangent_at_zero() {
        let e = xy_ellipse(3.0, 2.0);
        let t = e.tangent(0.0).unwrap();
        // At t=0: dx = 0, dy = b => tangent is +Y
        assert!((t - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn full_ellipse_is_closed() {
        let e = xy_ellipse(3.0, 2.0);
        assert!(e.is_closed());
    }

    #[test]
    fn partial_ellipse_is_not_closed() {
        let e = xy_ellipse_arc(3.0, 2.0, 0.0, std::f64::consts::PI);
        assert!(!e.is_closed());
    }

    #[test]
    fn domain_matches_angles() {
        let e = xy_ellipse_arc(3.0, 2.0, 0.5, 2.0);
        let d = e.domain();
        assert!((d.t_min - 0.5).abs() < TOLERANCE);
        assert!((d.t_max - 2.0).abs() < TOLERANCE);
    }

    #[test]
    fn circle_is_special_case() {
        // When semi_major == semi_minor, ellipse degenerates to circle
        let e = xy_ellipse(2.0, 2.0);
        let p = e.evaluate(FRAC_PI_2).unwrap();
        assert!((p - Point3::new(0.0, 2.0, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn invalid_semi_major() {
        let r = Ellipse::new(
            Point3::origin(), 0.0, 1.0, Vector3::z(), Vector3::x(), 0.0, TAU,
        );
        assert!(r.is_err());
    }

    #[test]
    fn invalid_semi_minor() {
        let r = Ellipse::new(
            Point3::origin(), 1.0, 0.0, Vector3::z(), Vector3::x(), 0.0, TAU,
        );
        assert!(r.is_err());
    }
}
