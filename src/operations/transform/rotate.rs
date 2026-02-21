use crate::error::{OperationError, Result};
use crate::math::{Matrix4, Point3, Vector3, TOLERANCE};
use crate::topology::{SolidId, TopologyStore};

use super::GeneralTransform;

/// Rotates a solid around an axis.
pub struct Rotate {
    solid: SolidId,
    axis_origin: Point3,
    axis_direction: Vector3,
    angle: f64,
}

impl Rotate {
    /// Creates a new `Rotate` operation.
    ///
    /// * `angle` - Rotation angle in radians.
    #[must_use]
    pub fn new(solid: SolidId, axis_origin: Point3, axis_direction: Vector3, angle: f64) -> Self {
        Self {
            solid,
            axis_origin,
            axis_direction,
            angle,
        }
    }

    /// Executes the rotation, modifying the solid in-place.
    ///
    /// # Errors
    ///
    /// Returns an error if the axis direction is zero-length.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<()> {
        let len = self.axis_direction.norm();
        if len < TOLERANCE {
            return Err(
                OperationError::InvalidInput("rotation axis must be non-zero".into()).into(),
            );
        }
        let axis = self.axis_direction / len;

        // Build rotation matrix: Translate to origin, rotate, translate back
        let t_neg = Matrix4::new_translation(&(-self.axis_origin.coords));
        let rot = rotation_matrix(&axis, self.angle);
        let t_pos = Matrix4::new_translation(&self.axis_origin.coords);
        let matrix = t_pos * rot * t_neg;

        GeneralTransform::new(self.solid, matrix).execute(store)
    }
}

/// Builds a 4x4 rotation matrix around a unit axis by an angle (Rodrigues).
#[allow(clippy::many_single_char_names)]
fn rotation_matrix(axis: &Vector3, angle: f64) -> Matrix4 {
    let c = angle.cos();
    let s = angle.sin();
    let t = 1.0 - c;
    let (x, y, z) = (axis.x, axis.y, axis.z);

    #[allow(clippy::suspicious_operation_groupings)]
    Matrix4::new(
        t * x * x + c,     t * x * y - s * z, t * x * z + s * y, 0.0,
        t * x * y + s * z, t * y * y + c,     t * y * z - s * x, 0.0,
        t * x * z - s * y, t * y * z + s * x, t * z * z + c,     0.0,
        0.0,               0.0,               0.0,               1.0,
    )
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use std::f64::consts::FRAC_PI_2;

    use super::*;
    use crate::operations::creation::{MakeFace, MakeWire};
    use crate::operations::shaping::Extrude;
    use crate::topology::TopologyStore;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    #[test]
    fn rotate_90_around_z() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![p(1.0, 0.0, 0.0), p(2.0, 0.0, 0.0), p(2.0, 1.0, 0.0), p(1.0, 1.0, 0.0)],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(wire, vec![]).execute(&mut store).unwrap();
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 1.0))
            .execute(&mut store)
            .unwrap();

        // Rotate 90 degrees around Z axis at origin
        Rotate::new(
            solid,
            p(0.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            FRAC_PI_2,
        )
        .execute(&mut store)
        .unwrap();

        // Point (1, 0, 0) -> (0, 1, 0), (2, 0, 0) -> (0, 2, 0), etc.
        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        for &fid in &shell.faces {
            let face = store.face(fid).unwrap();
            let wire = store.wire(face.outer_wire).unwrap();
            for oe in &wire.edges {
                let edge = store.edge(oe.edge).unwrap();
                let pt = store.vertex(edge.start).unwrap().point;
                // After 90° Z rotation: x ∈ [-1, 0], y ∈ [0, 2], z ∈ [0, 1]
                assert!(
                    pt.x >= -1.0 - 1e-6 && pt.x <= 0.0 + 1e-6,
                    "x={} out of range",
                    pt.x
                );
                assert!(
                    pt.y >= 0.0 - 1e-6 && pt.y <= 2.0 + 1e-6,
                    "y={} out of range",
                    pt.y
                );
            }
        }
    }

    #[test]
    fn zero_axis_returns_error() {
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

        let result = Rotate::new(
            solid,
            p(0.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 0.0),
            1.0,
        )
        .execute(&mut store);
        assert!(result.is_err());
    }
}
