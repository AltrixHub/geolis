mod cone;
mod cylinder;
mod plane;
mod sphere;
mod torus;

pub use cone::Cone;
pub use cylinder::Cylinder;
pub use plane::Plane;
pub use sphere::Sphere;
pub use torus::Torus;

use crate::error::Result;
use crate::math::{Point3, Vector3};

/// Parameter domain for a surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceDomain {
    /// Start of the U parameter range.
    pub u_min: f64,
    /// End of the U parameter range.
    pub u_max: f64,
    /// Start of the V parameter range.
    pub v_min: f64,
    /// End of the V parameter range.
    pub v_max: f64,
}

impl SurfaceDomain {
    /// Creates a new surface domain.
    #[must_use]
    pub fn new(u_min: f64, u_max: f64, v_min: f64, v_max: f64) -> Self {
        Self {
            u_min,
            u_max,
            v_min,
            v_max,
        }
    }
}

/// Trait for parametric surfaces in 3D space.
pub trait Surface {
    /// Evaluates the surface at parameters `(u, v)`, returning the 3D point.
    ///
    /// # Errors
    ///
    /// Returns an error if the parameters are out of range or evaluation fails.
    fn evaluate(&self, u: f64, v: f64) -> Result<Point3>;

    /// Computes the surface normal at parameters `(u, v)`.
    ///
    /// # Errors
    ///
    /// Returns an error if the parameters are out of range or the normal is degenerate.
    fn normal(&self, u: f64, v: f64) -> Result<Vector3>;

    /// Returns the parameter domain of the surface.
    fn domain(&self) -> SurfaceDomain;
}
