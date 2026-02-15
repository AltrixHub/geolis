use crate::error::{GeometryError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};

use super::{Surface, SurfaceDomain};

/// An infinite plane in 3D space.
///
/// Defined by an origin point, and two orthogonal direction vectors
/// (`u_dir`, `v_dir`). The normal is `u_dir × v_dir`.
///
/// Parametric form: `P(u, v) = origin + u * u_dir + v * v_dir`.
#[derive(Debug, Clone)]
pub struct Plane {
    origin: Point3,
    u_dir: Vector3,
    v_dir: Vector3,
    normal: Vector3,
}

impl Plane {
    /// Creates a new plane from an origin and two direction vectors.
    ///
    /// # Errors
    ///
    /// Returns an error if the direction vectors are zero-length
    /// or parallel (degenerate plane).
    pub fn new(origin: Point3, u_dir: Vector3, v_dir: Vector3) -> Result<Self> {
        let u_len = u_dir.norm();
        if u_len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        let v_len = v_dir.norm();
        if v_len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }

        let u_dir = u_dir / u_len;
        let v_dir = v_dir / v_len;

        let normal = u_dir.cross(&v_dir);
        let normal_len = normal.norm();
        if normal_len < TOLERANCE {
            return Err(
                GeometryError::Degenerate("plane directions are parallel".into()).into(),
            );
        }
        let normal = normal / normal_len;

        Ok(Self {
            origin,
            u_dir,
            v_dir,
            normal,
        })
    }

    /// Creates a plane from an origin and a normal vector.
    ///
    /// The U and V directions are computed automatically.
    ///
    /// # Errors
    ///
    /// Returns an error if the normal vector is zero-length.
    pub fn from_normal(origin: Point3, normal: Vector3) -> Result<Self> {
        let len = normal.norm();
        if len < TOLERANCE {
            return Err(GeometryError::ZeroVector.into());
        }
        let normal = normal / len;

        // Choose a reference vector not parallel to the normal
        let reference = if normal.x.abs() < 0.9 {
            Vector3::new(1.0, 0.0, 0.0)
        } else {
            Vector3::new(0.0, 1.0, 0.0)
        };

        let u_dir = normal.cross(&reference).normalize();
        let v_dir = normal.cross(&u_dir);

        Ok(Self {
            origin,
            u_dir,
            v_dir,
            normal,
        })
    }

    /// Returns the origin point of the plane.
    #[must_use]
    pub fn origin(&self) -> &Point3 {
        &self.origin
    }

    /// Returns the U direction vector.
    #[must_use]
    pub fn u_dir(&self) -> &Vector3 {
        &self.u_dir
    }

    /// Returns the V direction vector.
    #[must_use]
    pub fn v_dir(&self) -> &Vector3 {
        &self.v_dir
    }

    /// Returns the normal vector of the plane.
    #[must_use]
    pub fn plane_normal(&self) -> &Vector3 {
        &self.normal
    }
}

impl Surface for Plane {
    fn evaluate(&self, u: f64, v: f64) -> Result<Point3> {
        Ok(self.origin + self.u_dir * u + self.v_dir * v)
    }

    fn normal(&self, _u: f64, _v: f64) -> Result<Vector3> {
        Ok(self.normal)
    }

    fn domain(&self) -> SurfaceDomain {
        SurfaceDomain::new(f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY)
    }
}
