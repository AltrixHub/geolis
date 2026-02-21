use crate::error::{OperationError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};
use crate::operations::shaping::Extrude;
use crate::topology::{SolidId, TopologyStore};

use super::{MakeFace, MakeWire};

/// Creates a box solid from two corner points.
pub struct MakeBox {
    min_corner: Point3,
    max_corner: Point3,
}

impl MakeBox {
    /// Creates a new `MakeBox` operation.
    #[must_use]
    pub fn new(min_corner: Point3, max_corner: Point3) -> Self {
        Self {
            min_corner,
            max_corner,
        }
    }

    /// Executes the operation, creating the box in the topology store.
    ///
    /// Creates a bottom face from the four XY corners, then extrudes
    /// upward by the Z extent.
    ///
    /// # Errors
    ///
    /// Returns an error if the box is degenerate (any dimension is near zero).
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        let dx = self.max_corner.x - self.min_corner.x;
        let dy = self.max_corner.y - self.min_corner.y;
        let dz = self.max_corner.z - self.min_corner.z;

        if dx.abs() < TOLERANCE || dy.abs() < TOLERANCE || dz.abs() < TOLERANCE {
            return Err(OperationError::InvalidInput(
                "box dimensions must be non-zero".into(),
            )
            .into());
        }

        let (x0, y0, z0) = (self.min_corner.x, self.min_corner.y, self.min_corner.z);
        let (x1, y1) = (self.max_corner.x, self.max_corner.y);

        // Bottom face (CCW when viewed from +Z for positive dz)
        let bottom_pts = if dz > 0.0 {
            vec![
                Point3::new(x0, y0, z0),
                Point3::new(x1, y0, z0),
                Point3::new(x1, y1, z0),
                Point3::new(x0, y1, z0),
            ]
        } else {
            // Reverse winding for negative dz so extrude works correctly
            vec![
                Point3::new(x0, y0, z0),
                Point3::new(x0, y1, z0),
                Point3::new(x1, y1, z0),
                Point3::new(x1, y0, z0),
            ]
        };

        let wire = MakeWire::new(bottom_pts, true).execute(store)?;
        let face = MakeFace::new(wire, vec![]).execute(store)?;
        Extrude::new(face, Vector3::new(0.0, 0.0, dz)).execute(store)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::operations::query::{BoundingBox, IsValid};
    use crate::topology::TopologyStore;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    #[test]
    fn unit_box_has_6_faces() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(1.0, 1.0, 1.0))
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        assert_eq!(shell.faces.len(), 6);
        assert!(shell.is_closed);
    }

    #[test]
    fn box_aabb_matches_corners() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(1.0, 2.0, 3.0), p(4.0, 6.0, 8.0))
            .execute(&mut store)
            .unwrap();

        let aabb = BoundingBox::new(solid).execute(&store).unwrap();
        assert!((aabb.min.x - 1.0).abs() < 1e-10);
        assert!((aabb.min.y - 2.0).abs() < 1e-10);
        assert!((aabb.min.z - 3.0).abs() < 1e-10);
        assert!((aabb.max.x - 4.0).abs() < 1e-10);
        assert!((aabb.max.y - 6.0).abs() < 1e-10);
        assert!((aabb.max.z - 8.0).abs() < 1e-10);
    }

    #[test]
    fn box_is_valid() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(2.0, 3.0, 4.0))
            .execute(&mut store)
            .unwrap();

        assert!(IsValid::new(solid).execute(&store));
    }

    #[test]
    fn degenerate_box_fails() {
        let mut store = TopologyStore::new();
        let result = MakeBox::new(p(0.0, 0.0, 0.0), p(1.0, 1.0, 0.0)).execute(&mut store);
        assert!(result.is_err());
    }
}
