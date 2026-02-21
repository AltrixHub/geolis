use crate::error::{OperationError, Result};
use crate::math::TOLERANCE;
use crate::operations::shaping::Extrude;
use crate::topology::{FaceId, FaceSurface, SolidId, TopologyStore};

/// Thickens a face into a solid by offsetting in the normal direction.
///
/// Creates a solid by extruding the face along its normal by the given
/// thickness. Positive thickness extrudes in the normal direction;
/// negative thickness extrudes opposite to the normal.
pub struct ThickenFace {
    face: FaceId,
    thickness: f64,
}

impl ThickenFace {
    /// Creates a new `ThickenFace` operation.
    #[must_use]
    pub fn new(face: FaceId, thickness: f64) -> Self {
        Self { face, thickness }
    }

    /// Executes the thickening, creating the result solid in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if the thickness is near zero or the face cannot be extruded.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        if self.thickness.abs() < TOLERANCE {
            return Err(
                OperationError::InvalidInput("thickness must be non-zero".into()).into(),
            );
        }

        let face = store.face(self.face)?;
        let normal = match &face.surface {
            FaceSurface::Plane(plane) => {
                let n = *plane.plane_normal();
                if face.same_sense { n } else { -n }
            }
        };

        let direction = normal * self.thickness;
        Extrude::new(self.face, direction).execute(store)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use crate::operations::creation::{MakeFace, MakeWire};
    use crate::operations::query::{BoundingBox, IsValid};
    use crate::topology::TopologyStore;

    fn p(x: f64, y: f64) -> Point3 {
        Point3::new(x, y, 0.0)
    }

    #[test]
    fn thicken_square_face() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![p(0.0, 0.0), p(4.0, 0.0), p(4.0, 3.0), p(0.0, 3.0)],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(wire, vec![]).execute(&mut store).unwrap();

        let solid = ThickenFace::new(face, 2.0)
            .execute(&mut store)
            .unwrap();

        let aabb = BoundingBox::new(solid).execute(&store).unwrap();
        assert!((aabb.max.z - 2.0).abs() < 1e-10 || (aabb.min.z + 2.0).abs() < 1e-10);
        assert!(IsValid::new(solid).execute(&store));
    }

    #[test]
    fn thicken_negative_direction() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![p(0.0, 0.0), p(2.0, 0.0), p(2.0, 2.0), p(0.0, 2.0)],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(wire, vec![]).execute(&mut store).unwrap();

        let solid = ThickenFace::new(face, -3.0)
            .execute(&mut store)
            .unwrap();

        let aabb = BoundingBox::new(solid).execute(&store).unwrap();
        // Should extend in negative Z direction
        assert!(aabb.min.z < -1.0);
        assert!(IsValid::new(solid).execute(&store));
    }

    #[test]
    fn zero_thickness_returns_error() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0)],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(wire, vec![]).execute(&mut store).unwrap();

        let result = ThickenFace::new(face, 0.0).execute(&mut store);
        assert!(result.is_err());
    }
}
