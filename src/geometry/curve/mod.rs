mod arc;
mod line;

pub use arc::Arc;
pub use line::Line;

use crate::error::Result;
use crate::math::{Point3, Vector3};

/// Parameter domain for a curve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveDomain {
    /// Start of the parameter range.
    pub t_min: f64,
    /// End of the parameter range.
    pub t_max: f64,
}

impl CurveDomain {
    /// Creates a new curve domain.
    #[must_use]
    pub fn new(t_min: f64, t_max: f64) -> Self {
        Self { t_min, t_max }
    }
}

/// Trait for parametric curves in 3D space.
pub trait Curve {
    /// Evaluates the curve at parameter `t`, returning the 3D point.
    ///
    /// # Errors
    ///
    /// Returns an error if the parameter is out of range or evaluation fails.
    fn evaluate(&self, t: f64) -> Result<Point3>;

    /// Computes the tangent vector at parameter `t`.
    ///
    /// # Errors
    ///
    /// Returns an error if the parameter is out of range or the tangent is degenerate.
    fn tangent(&self, t: f64) -> Result<Vector3>;

    /// Returns the parameter domain of the curve.
    fn domain(&self) -> CurveDomain;

    /// Returns whether the curve is closed.
    fn is_closed(&self) -> bool;
}
