use crate::error::Result;
use crate::math::{Matrix4, Vector3};
use crate::topology::{SolidId, TopologyStore};

use super::GeneralTransform;

/// Translates a solid by a displacement vector.
pub struct Translate {
    solid: SolidId,
    displacement: Vector3,
}

impl Translate {
    /// Creates a new `Translate` operation.
    #[must_use]
    pub fn new(solid: SolidId, displacement: Vector3) -> Self {
        Self {
            solid,
            displacement,
        }
    }

    /// Executes the translation, modifying the solid in-place.
    ///
    /// # Errors
    ///
    /// Returns an error if the underlying transform fails.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<()> {
        let mut matrix = Matrix4::identity();
        matrix[(0, 3)] = self.displacement.x;
        matrix[(1, 3)] = self.displacement.y;
        matrix[(2, 3)] = self.displacement.z;
        GeneralTransform::new(self.solid, matrix).execute(store)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use crate::operations::creation::{MakeFace, MakeWire};
    use crate::operations::shaping::Extrude;
    use crate::topology::TopologyStore;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    #[test]
    fn translate_shifts_vertices() {
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

        Translate::new(solid, Vector3::new(10.0, 20.0, 30.0))
            .execute(&mut store)
            .unwrap();

        // All vertices should be in [10, 11] x [20, 21] x [30, 31]
        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        for &fid in &shell.faces {
            let face = store.face(fid).unwrap();
            let wire = store.wire(face.outer_wire).unwrap();
            for oe in &wire.edges {
                let edge = store.edge(oe.edge).unwrap();
                let pt = store.vertex(edge.start).unwrap().point;
                assert!(pt.x >= 10.0 - 1e-10 && pt.x <= 11.0 + 1e-10);
                assert!(pt.y >= 20.0 - 1e-10 && pt.y <= 21.0 + 1e-10);
                assert!(pt.z >= 30.0 - 1e-10 && pt.z <= 31.0 + 1e-10);
            }
        }
    }
}
