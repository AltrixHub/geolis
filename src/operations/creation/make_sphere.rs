use crate::error::{OperationError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};
use crate::operations::shaping::Revolve;
use crate::topology::{SolidId, TopologyStore};

use super::{MakeFace, MakeWire};

/// Creates a sphere solid from center and radius.
///
/// Internally constructs a semicircular triangle profile (south pole, equator,
/// north pole) and revolves it 360 degrees around the vertical axis.
///
/// Note: The resulting surface uses Cone approximation (one cone for each half),
/// not a true Sphere surface. This is because Revolve generates surfaces based
/// on the profile edge geometry.
pub struct MakeSphere {
    center: Point3,
    radius: f64,
}

impl MakeSphere {
    /// Creates a new `MakeSphere` operation.
    #[must_use]
    pub fn new(center: Point3, radius: f64) -> Self {
        Self { center, radius }
    }

    /// Executes the operation, creating the sphere in the topology store.
    ///
    /// Builds a triangular profile (south pole -> equator -> north pole) and
    /// revolves it 360 degrees around the Z axis passing through the center.
    ///
    /// The result has 2 Cone faces (upper and lower hemisphere). The edge
    /// connecting the poles (both on-axis) is degenerate and skipped.
    ///
    /// # Errors
    ///
    /// Returns an error if the radius is near zero.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        if self.radius < TOLERANCE {
            return Err(
                OperationError::InvalidInput("sphere radius must be positive".into()).into(),
            );
        }

        let r = self.radius;
        let c = self.center;

        // Triangle profile in the XZ half-plane (y = center.y):
        // South pole (on axis) -> Equator (off axis) -> North pole (on axis)
        let south = Point3::new(c.x, c.y, c.z - r);
        let equator = Point3::new(c.x + r, c.y, c.z);
        let north = Point3::new(c.x, c.y, c.z + r);

        let profile = vec![south, equator, north];

        let wire = MakeWire::new(profile, true).execute(store)?;
        let face = MakeFace::new(wire, vec![]).execute(store)?;
        Revolve::new(face, self.center, Vector3::z()).execute(store)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::operations::query::{BoundingBox, IsValid};
    use crate::tessellation::{TessellateSolid, TessellationParams};

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    #[test]
    fn sphere_has_2_faces() {
        let mut store = TopologyStore::new();
        let solid = MakeSphere::new(p(0.0, 0.0, 0.0), 3.0)
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        // Triangle profile: 3 edges
        // south->equator: on-axis to off-axis -> Cone
        // equator->north: off-axis to on-axis -> Cone
        // north->south: both on-axis -> degenerate (skipped)
        assert_eq!(shell.faces.len(), 2);
        assert!(shell.is_closed);
    }

    #[test]
    fn sphere_bounding_box() {
        let mut store = TopologyStore::new();
        let solid = MakeSphere::new(p(0.0, 0.0, 0.0), 3.0)
            .execute(&mut store)
            .unwrap();

        let aabb = BoundingBox::new(solid).execute(&store).unwrap();
        assert!((aabb.min.z - (-3.0)).abs() < 1e-6);
        assert!((aabb.max.z - 3.0).abs() < 1e-6);
    }

    #[test]
    fn sphere_is_valid() {
        let mut store = TopologyStore::new();
        let solid = MakeSphere::new(p(0.0, 0.0, 0.0), 3.0)
            .execute(&mut store)
            .unwrap();

        assert!(IsValid::new(solid).execute(&store));
    }

    #[test]
    fn sphere_tessellates() {
        let mut store = TopologyStore::new();
        let solid = MakeSphere::new(p(0.0, 0.0, 0.0), 3.0)
            .execute(&mut store)
            .unwrap();

        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty());
        assert_eq!(mesh.vertices.len(), mesh.normals.len());
    }

    #[test]
    fn zero_radius_fails() {
        let mut store = TopologyStore::new();
        let result = MakeSphere::new(p(0.0, 0.0, 0.0), 0.0).execute(&mut store);
        assert!(result.is_err());
    }
}
