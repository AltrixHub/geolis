use crate::error::Result;
use crate::math::Point3;
use crate::topology::{SolidId, TopologyStore};

/// An axis-aligned bounding box.
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    /// Minimum corner of the bounding box.
    pub min: Point3,
    /// Maximum corner of the bounding box.
    pub max: Point3,
}

/// Computes the axis-aligned bounding box of a solid.
pub struct BoundingBox {
    solid: SolidId,
}

impl BoundingBox {
    /// Creates a new `BoundingBox` query.
    #[must_use]
    pub fn new(solid: SolidId) -> Self {
        Self { solid }
    }

    /// Executes the query, returning the AABB.
    ///
    /// Iterates over all vertices in the solid to compute min/max coordinates.
    ///
    /// # Errors
    ///
    /// Returns an error if the solid or any referenced entity is not found.
    pub fn execute(&self, store: &TopologyStore) -> Result<Aabb> {
        let solid = store.solid(self.solid)?;
        let outer_shell_id = solid.outer_shell;
        let inner_shell_ids = solid.inner_shells.clone();

        let mut min = Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut max = Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);

        let mut process_shell =
            |shell_id: crate::topology::ShellId| -> Result<()> {
                let shell = store.shell(shell_id)?;
                for &face_id in &shell.faces {
                    let face = store.face(face_id)?;
                    let wire_ids: Vec<_> = std::iter::once(face.outer_wire)
                        .chain(face.inner_wires.iter().copied())
                        .collect();

                    for wire_id in wire_ids {
                        let wire = store.wire(wire_id)?;
                        for oe in &wire.edges {
                            let edge = store.edge(oe.edge)?;
                            for &vid in &[edge.start, edge.end] {
                                let pt = store.vertex(vid)?.point;
                                min.x = min.x.min(pt.x);
                                min.y = min.y.min(pt.y);
                                min.z = min.z.min(pt.z);
                                max.x = max.x.max(pt.x);
                                max.y = max.y.max(pt.y);
                                max.z = max.z.max(pt.z);
                            }
                        }
                    }
                }
                Ok(())
            };

        process_shell(outer_shell_id)?;
        for &shell_id in &inner_shell_ids {
            process_shell(shell_id)?;
        }

        Ok(Aabb { min, max })
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Vector3;
    use crate::operations::creation::{MakeFace, MakeWire};
    use crate::operations::shaping::Extrude;
    use crate::topology::TopologyStore;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    #[test]
    fn unit_cube_aabb() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(0.0, 1.0, 0.0)],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(wire, vec![]).execute(&mut store).unwrap();
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 1.0))
            .execute(&mut store)
            .unwrap();

        let aabb = BoundingBox::new(solid).execute(&store).unwrap();
        assert!((aabb.min.x).abs() < 1e-10);
        assert!((aabb.min.y).abs() < 1e-10);
        assert!((aabb.min.z).abs() < 1e-10);
        assert!((aabb.max.x - 1.0).abs() < 1e-10);
        assert!((aabb.max.y - 1.0).abs() < 1e-10);
        assert!((aabb.max.z - 1.0).abs() < 1e-10);
    }

    #[test]
    fn offset_box_aabb() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![p(2.0, 3.0, 0.0), p(5.0, 3.0, 0.0), p(5.0, 7.0, 0.0), p(2.0, 7.0, 0.0)],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(wire, vec![]).execute(&mut store).unwrap();
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 4.0))
            .execute(&mut store)
            .unwrap();

        let aabb = BoundingBox::new(solid).execute(&store).unwrap();
        assert!((aabb.min.x - 2.0).abs() < 1e-10);
        assert!((aabb.min.y - 3.0).abs() < 1e-10);
        assert!((aabb.min.z).abs() < 1e-10);
        assert!((aabb.max.x - 5.0).abs() < 1e-10);
        assert!((aabb.max.y - 7.0).abs() < 1e-10);
        assert!((aabb.max.z - 4.0).abs() < 1e-10);
    }
}
