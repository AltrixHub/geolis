use crate::error::Result;
use crate::math::{Point3, Vector3};

use super::{Curve, CurveDomain};

/// An infinite line defined by an origin point and a direction vector.
///
/// The parametric form is: `P(t) = origin + t * direction`.
#[derive(Debug, Clone)]
pub struct Line {
    origin: Point3,
    direction: Vector3,
}

impl Line {
    /// Creates a new line from an origin and direction.
    ///
    /// # Errors
    ///
    /// Returns an error if the direction vector is zero-length.
    pub fn new(origin: Point3, direction: Vector3) -> Result<Self> {
        let len = direction.norm();
        if len < crate::math::TOLERANCE {
            return Err(crate::error::GeometryError::ZeroVector.into());
        }
        Ok(Self {
            origin,
            direction: direction / len,
        })
    }

    /// Returns the origin point of the line.
    #[must_use]
    pub fn origin(&self) -> &Point3 {
        &self.origin
    }

    /// Returns the unit direction vector of the line.
    #[must_use]
    pub fn direction(&self) -> &Vector3 {
        &self.direction
    }
}

impl Curve for Line {
    fn evaluate(&self, t: f64) -> Result<Point3> {
        Ok(self.origin + self.direction * t)
    }

    fn tangent(&self, _t: f64) -> Result<Vector3> {
        Ok(self.direction)
    }

    fn domain(&self) -> CurveDomain {
        CurveDomain::new(f64::NEG_INFINITY, f64::INFINITY)
    }

    fn is_closed(&self) -> bool {
        false
    }
}
