use crate::error::{OperationError, Result};
use crate::math::{Matrix4, Point3, TOLERANCE};
use crate::topology::{SolidId, TopologyStore};

use super::GeneralTransform;

/// Scales a solid uniformly from a center point.
pub struct Scale {
    solid: SolidId,
    center: Point3,
    factor: f64,
}

impl Scale {
    /// Creates a new `Scale` operation.
    #[must_use]
    pub fn new(solid: SolidId, center: Point3, factor: f64) -> Self {
        Self {
            solid,
            center,
            factor,
        }
    }

    /// Executes the scaling, modifying the solid in-place.
    ///
    /// # Errors
    ///
    /// Returns an error if the scale factor is near zero.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<()> {
        if self.factor.abs() < TOLERANCE {
            return Err(
                OperationError::InvalidInput("scale factor must be non-zero".into()).into(),
            );
        }

        // Translate to origin, scale, translate back
        let t_neg = Matrix4::new_translation(&(-self.center.coords));
        let s = Matrix4::new_scaling(self.factor);
        let t_pos = Matrix4::new_translation(&self.center.coords);
        let matrix = t_pos * s * t_neg;

        GeneralTransform::new(self.solid, matrix).execute(store)
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
    fn scale_from_origin_doubles() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![p(1.0, 1.0, 0.0), p(2.0, 1.0, 0.0), p(2.0, 2.0, 0.0), p(1.0, 2.0, 0.0)],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(wire, vec![]).execute(&mut store).unwrap();
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 1.0))
            .execute(&mut store)
            .unwrap();

        Scale::new(solid, p(0.0, 0.0, 0.0), 2.0)
            .execute(&mut store)
            .unwrap();

        // Vertices should now be in [2, 4] x [2, 4] x [0, 2]
        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        for &fid in &shell.faces {
            let face = store.face(fid).unwrap();
            let wire = store.wire(face.outer_wire).unwrap();
            for oe in &wire.edges {
                let edge = store.edge(oe.edge).unwrap();
                let pt = store.vertex(edge.start).unwrap().point;
                assert!(pt.x >= 2.0 - 1e-6 && pt.x <= 4.0 + 1e-6);
                assert!(pt.y >= 2.0 - 1e-6 && pt.y <= 4.0 + 1e-6);
                assert!(pt.z >= 0.0 - 1e-6 && pt.z <= 2.0 + 1e-6);
            }
        }
    }

    #[test]
    fn zero_factor_returns_error() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0)],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(wire, vec![]).execute(&mut store).unwrap();
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 1.0))
            .execute(&mut store)
            .unwrap();

        let result = Scale::new(solid, p(0.0, 0.0, 0.0), 0.0).execute(&mut store);
        assert!(result.is_err());
    }
}
