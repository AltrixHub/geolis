use crate::error::Result;
use crate::geometry::surface::Surface;
use crate::math::Point3;
use crate::topology::{FaceId, FaceSurface, TopologyStore};

/// Evaluates a point on a surface at given parameters.
pub struct PointOnSurface {
    face: FaceId,
    u: f64,
    v: f64,
}

impl PointOnSurface {
    /// Creates a new `PointOnSurface` query.
    #[must_use]
    pub fn new(face: FaceId, u: f64, v: f64) -> Self {
        Self { face, u, v }
    }

    /// Executes the query, returning the 3D point.
    ///
    /// # Errors
    ///
    /// Returns an error if the face is not found or evaluation fails.
    pub fn execute(&self, store: &TopologyStore) -> Result<Point3> {
        let face = store.face(self.face)?;
        match &face.surface {
            FaceSurface::Plane(plane) => plane.evaluate(self.u, self.v),
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use crate::operations::creation::{MakeFace, MakeWire};
    use crate::topology::TopologyStore;

    #[test]
    fn point_on_xy_plane() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(4.0, 0.0, 0.0),
                Point3::new(4.0, 4.0, 0.0),
                Point3::new(0.0, 4.0, 0.0),
            ],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(wire, vec![]).execute(&mut store).unwrap();

        let pt = PointOnSurface::new(face, 1.0, 2.0)
            .execute(&store)
            .unwrap();
        // The plane is constructed from the wire centroid; exact coords depend
        // on the plane's u_dir/v_dir, but z should always be 0
        assert!(pt.z.abs() < 1e-10);
    }
}
