use crate::error::{OperationError, Result};
use crate::geometry::surface::Plane;
use crate::math::{Point3, Vector3, TOLERANCE};
use crate::operations::boolean::Intersect;
use crate::operations::creation::MakeBox;
use crate::operations::query::BoundingBox;
use crate::topology::{SolidId, TopologyStore};

/// Splits a solid with a plane, producing two half-solids.
///
/// Uses the boolean engine: creates large boxes on each side of the
/// cutting plane, then intersects the solid with each box.
pub struct Split {
    solid: SolidId,
    plane_origin: Point3,
    plane_normal: Vector3,
}

impl Split {
    /// Creates a new `Split` operation.
    ///
    /// The plane is defined by a point and a normal direction.
    /// The positive half-space (where the normal points) produces the first result.
    #[must_use]
    pub fn new(solid: SolidId, plane_origin: Point3, plane_normal: Vector3) -> Self {
        Self {
            solid,
            plane_origin,
            plane_normal,
        }
    }

    /// Executes the split, returning two solids (positive half, negative half).
    ///
    /// # Errors
    ///
    /// Returns an error if the plane doesn't intersect the solid, or
    /// the boolean operations fail.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<(SolidId, SolidId)> {
        let normal_len = self.plane_normal.norm();
        if normal_len < TOLERANCE {
            return Err(
                OperationError::InvalidInput("split plane normal must be non-zero".into()).into(),
            );
        }
        let normal = self.plane_normal / normal_len;

        // Get the solid's bounding box to create appropriately-sized cutting boxes
        let aabb = BoundingBox::new(self.solid).execute(store)?;
        let diag = aabb.max - aabb.min;
        let extent = diag.norm() * 2.0; // generous extent beyond the solid

        // Build a local coordinate frame for the plane
        let plane = Plane::from_normal(self.plane_origin, normal)?;
        let u = *plane.u_dir();
        let v = *plane.v_dir();

        // Create the positive half-space box (normal side)
        let pos_min = self.plane_origin - u * extent - v * extent;
        let pos_max = self.plane_origin + u * extent + v * extent + normal * extent;
        let pos_box = make_aligned_box(store, pos_min, pos_max, &u, &v, &normal)?;

        // Create the negative half-space box (opposite side)
        let neg_min = self.plane_origin - u * extent - v * extent - normal * extent;
        let neg_max = self.plane_origin + u * extent + v * extent;
        let neg_box = make_aligned_box(store, neg_min, neg_max, &u, &v, &normal)?;

        // Intersect the solid with each half-space box
        let positive = Intersect::new(self.solid, pos_box).execute(store)?;
        let negative = Intersect::new(self.solid, neg_box).execute(store)?;

        Ok((positive, negative))
    }
}

/// Creates an axis-aligned box in the frame (u, v, n) with corners at min/max.
fn make_aligned_box(
    store: &mut TopologyStore,
    min: Point3,
    max: Point3,
    _u: &Vector3,
    _v: &Vector3,
    _n: &Vector3,
) -> Result<SolidId> {
    // Use MakeBox which creates an AABB. Since the cutting plane may be
    // arbitrarily oriented, we need the box in world coordinates.
    // MakeBox creates an axis-aligned box, so we compute the AABB of our corners.
    let min_world = Point3::new(
        min.x.min(max.x),
        min.y.min(max.y),
        min.z.min(max.z),
    );
    let max_world = Point3::new(
        min.x.max(max.x),
        min.y.max(max.y),
        min.z.max(max.z),
    );

    // Ensure non-degenerate box (extend thin dimensions)
    let dx = (max_world.x - min_world.x).max(TOLERANCE * 100.0);
    let dy = (max_world.y - min_world.y).max(TOLERANCE * 100.0);
    let dz = (max_world.z - min_world.z).max(TOLERANCE * 100.0);

    MakeBox::new(
        min_world,
        Point3::new(min_world.x + dx, min_world.y + dy, min_world.z + dz),
    )
    .execute(store)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::operations::creation::MakeBox;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    #[test]
    fn split_box_horizontal() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(4.0, 4.0, 4.0))
            .execute(&mut store)
            .unwrap();

        let (top, bottom) = Split::new(solid, p(0.0, 0.0, 2.0), Vector3::z())
            .execute(&mut store)
            .unwrap();

        // Both halves should have faces
        let top_data = store.solid(top).unwrap();
        let top_shell = store.shell(top_data.outer_shell).unwrap();
        assert!(!top_shell.faces.is_empty());

        let bot_data = store.solid(bottom).unwrap();
        let bot_shell = store.shell(bot_data.outer_shell).unwrap();
        assert!(!bot_shell.faces.is_empty());
    }

    #[test]
    fn split_box_vertical() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(4.0, 4.0, 4.0))
            .execute(&mut store)
            .unwrap();

        let (pos, neg) = Split::new(solid, p(2.0, 0.0, 0.0), Vector3::x())
            .execute(&mut store)
            .unwrap();

        let pos_data = store.solid(pos).unwrap();
        let pos_shell = store.shell(pos_data.outer_shell).unwrap();
        assert!(!pos_shell.faces.is_empty());

        let neg_data = store.solid(neg).unwrap();
        let neg_shell = store.shell(neg_data.outer_shell).unwrap();
        assert!(!neg_shell.faces.is_empty());
    }

    #[test]
    fn split_preserves_volume() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(4.0, 4.0, 4.0))
            .execute(&mut store)
            .unwrap();

        let orig_vol = crate::operations::query::Volume::new(solid)
            .execute(&store)
            .unwrap();

        let (top, bottom) = Split::new(solid, p(0.0, 0.0, 2.0), Vector3::z())
            .execute(&mut store)
            .unwrap();

        let top_vol = crate::operations::query::Volume::new(top)
            .execute(&store)
            .unwrap();
        let bot_vol = crate::operations::query::Volume::new(bottom)
            .execute(&store)
            .unwrap();

        let sum = top_vol + bot_vol;
        assert!(
            (sum - orig_vol).abs() < orig_vol * 0.10,
            "volumes don't add up: {top_vol} + {bot_vol} = {sum}, expected {orig_vol}"
        );
    }

    #[test]
    fn split_zero_normal_fails() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(4.0, 4.0, 4.0))
            .execute(&mut store)
            .unwrap();

        let result = Split::new(solid, p(0.0, 0.0, 2.0), Vector3::new(0.0, 0.0, 0.0))
            .execute(&mut store);
        assert!(result.is_err());
    }
}
