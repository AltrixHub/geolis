use crate::error::{GeometryError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};

use super::{Curve, CurveDomain};

/// A circular arc in 3D space.
///
/// Defined by a center, radius, normal axis, and a reference direction
/// for the zero-angle. The parametric form sweeps from `start_angle`
/// to `end_angle` (in radians) around the normal axis.
#[derive(Debug, Clone)]
pub struct Arc {
    center: Point3,
    radius: f64,
    normal: Vector3,
    ref_dir: Vector3,
    start_angle: f64,
    end_angle: f64,
}

impl Arc {
    /// Creates a new arc.
    ///
    /// # Arguments
    ///
    /// * `center` - Center of the arc circle
    /// * `radius` - Radius (must be positive)
    /// * `normal` - Normal vector defining the arc plane
    /// * `ref_dir` - Reference direction for angle = 0 (must be perpendicular to normal)
    /// * `start_angle` - Start angle in radians
    /// * `end_angle` - End angle in radians
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
        start_angle: f64,
        end_angle: f64,
    ) -> Result<Self> {
        if radius < TOLERANCE {
            return Err(GeometryError::Degenerate("arc radius must be positive".into()).into());
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
            start_angle,
            end_angle,
        })
    }

    /// Returns the center of the arc.
    #[must_use]
    pub fn center(&self) -> &Point3 {
        &self.center
    }

    /// Returns the radius of the arc.
    #[must_use]
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Returns the normal vector of the arc plane.
    #[must_use]
    pub fn normal(&self) -> &Vector3 {
        &self.normal
    }

    /// Computes the second axis direction (perpendicular to both normal and `ref_dir`).
    fn binormal(&self) -> Vector3 {
        self.normal.cross(&self.ref_dir)
    }
}

impl Curve for Arc {
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
        CurveDomain::new(self.start_angle, self.end_angle)
    }

    fn is_closed(&self) -> bool {
        (self.end_angle - self.start_angle - std::f64::consts::TAU).abs() < TOLERANCE
    }
}
