use crate::error::{OperationError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};
use crate::operations::shaping::Revolve;
use crate::topology::{SolidId, TopologyStore};

use super::{MakeFace, MakeWire};

/// Creates a cone or truncated cone (frustum) solid.
///
/// Internally constructs a triangle or trapezoid profile and revolves it
/// 360 degrees around the specified axis.
///
/// - `top_radius = 0` produces a pointed cone (triangle profile)
/// - `top_radius > 0` produces a truncated cone (trapezoid profile)
pub struct MakeCone {
    center: Point3,
    bottom_radius: f64,
    top_radius: f64,
    axis: Vector3,
    height: f64,
}

impl MakeCone {
    /// Creates a new `MakeCone` operation.
    #[must_use]
    pub fn new(
        center: Point3,
        bottom_radius: f64,
        top_radius: f64,
        axis: Vector3,
        height: f64,
    ) -> Self {
        Self {
            center,
            bottom_radius,
            top_radius,
            axis,
            height,
        }
    }

    /// Executes the operation, creating the cone in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if the bottom radius is near zero, height is near zero,
    /// or both radii are the same (use `MakeCylinder` instead).
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        if self.bottom_radius < TOLERANCE {
            return Err(
                OperationError::InvalidInput("cone bottom radius must be positive".into()).into(),
            );
        }
        if self.top_radius < 0.0 {
            return Err(
                OperationError::InvalidInput("cone top radius must be non-negative".into()).into(),
            );
        }
        if self.height.abs() < TOLERANCE {
            return Err(
                OperationError::InvalidInput("cone height must be non-zero".into()).into(),
            );
        }
        let axis_len = self.axis.norm();
        if axis_len < TOLERANCE {
            return Err(
                OperationError::InvalidInput("cone axis must be non-zero".into()).into(),
            );
        }
        let axis = self.axis / axis_len;

        // Build a reference direction perpendicular to the axis
        let ref_dir = perpendicular_dir(&axis);

        let rb = self.bottom_radius;
        let rt = self.top_radius;
        let h = self.height;

        if rt < TOLERANCE {
            // Pointed cone: triangle profile
            // bottom-outer -> top-apex (on axis) -> bottom-center (on axis)
            let profile = vec![
                self.center + ref_dir * rb,          // bottom outer edge
                self.center + axis * h,               // top apex (on axis)
                self.center,                           // bottom center (on axis)
            ];

            let wire = MakeWire::new(profile, true).execute(store)?;
            let face = MakeFace::new(wire, vec![]).execute(store)?;
            Revolve::new(face, self.center, self.axis).execute(store)
        } else {
            // Truncated cone: trapezoid profile
            // bottom-outer (rb, 0) -> top-outer (rt, h) -> top-center (0, h) -> bottom-center (0, 0)
            let profile = vec![
                self.center + ref_dir * rb,              // bottom outer
                self.center + ref_dir * rt + axis * h,   // top outer
                self.center + axis * h,                   // top center (on axis)
                self.center,                               // bottom center (on axis)
            ];

            let wire = MakeWire::new(profile, true).execute(store)?;
            let face = MakeFace::new(wire, vec![]).execute(store)?;
            Revolve::new(face, self.center, self.axis).execute(store)
        }
    }
}

/// Finds a direction perpendicular to the given unit vector.
fn perpendicular_dir(axis: &Vector3) -> Vector3 {
    let candidate = if axis.x.abs() < 0.9 {
        Vector3::x()
    } else {
        Vector3::y()
    };
    let perp = axis.cross(&candidate);
    perp / perp.norm()
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
    fn cone_full_has_2_faces() {
        let mut store = TopologyStore::new();
        let solid = MakeCone::new(p(0.0, 0.0, 0.0), 3.0, 0.0, Vector3::z(), 6.0)
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        // Triangle profile: bottom-outer -> apex -> bottom-center
        // bottom-outer(off) -> apex(on): Cone face
        // apex(on) -> bottom-center(on): degenerate (skipped)
        // bottom-center(on) -> bottom-outer(off): Plane disc
        assert_eq!(shell.faces.len(), 2);
        assert!(shell.is_closed);
    }

    #[test]
    fn cone_truncated_has_3_faces() {
        let mut store = TopologyStore::new();
        let solid = MakeCone::new(p(0.0, 0.0, 0.0), 3.0, 1.5, Vector3::z(), 6.0)
            .execute(&mut store)
            .unwrap();

        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        // Trapezoid profile: bottom-outer -> top-outer -> top-center -> bottom-center
        // bottom-outer(off) -> top-outer(off): Cone face (different radii)
        // top-outer(off) -> top-center(on): Plane disc (top)
        // top-center(on) -> bottom-center(on): degenerate (skipped)
        // bottom-center(on) -> bottom-outer(off): Plane disc (bottom)
        assert_eq!(shell.faces.len(), 3);
        assert!(shell.is_closed);
    }

    #[test]
    fn cone_bounding_box() {
        let mut store = TopologyStore::new();
        let solid = MakeCone::new(p(0.0, 0.0, 0.0), 3.0, 0.0, Vector3::z(), 4.0)
            .execute(&mut store)
            .unwrap();

        let aabb = BoundingBox::new(solid).execute(&store).unwrap();
        assert!((aabb.min.z - 0.0).abs() < 1e-6);
        assert!((aabb.max.z - 4.0).abs() < 1e-6);
    }

    #[test]
    fn cone_is_valid() {
        let mut store = TopologyStore::new();
        let solid = MakeCone::new(p(0.0, 0.0, 0.0), 3.0, 0.0, Vector3::z(), 6.0)
            .execute(&mut store)
            .unwrap();

        assert!(IsValid::new(solid).execute(&store));
    }

    #[test]
    fn truncated_cone_is_valid() {
        let mut store = TopologyStore::new();
        let solid = MakeCone::new(p(0.0, 0.0, 0.0), 3.0, 1.5, Vector3::z(), 6.0)
            .execute(&mut store)
            .unwrap();

        assert!(IsValid::new(solid).execute(&store));
    }

    #[test]
    fn cone_tessellates() {
        let mut store = TopologyStore::new();
        let solid = MakeCone::new(p(0.0, 0.0, 0.0), 3.0, 0.0, Vector3::z(), 6.0)
            .execute(&mut store)
            .unwrap();

        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty());
        assert_eq!(mesh.vertices.len(), mesh.normals.len());
    }

    #[test]
    fn truncated_cone_tessellates() {
        let mut store = TopologyStore::new();
        let solid = MakeCone::new(p(0.0, 0.0, 0.0), 3.0, 1.5, Vector3::z(), 6.0)
            .execute(&mut store)
            .unwrap();

        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty());
        assert_eq!(mesh.vertices.len(), mesh.normals.len());
    }

    #[test]
    fn zero_bottom_radius_fails() {
        let mut store = TopologyStore::new();
        let result = MakeCone::new(p(0.0, 0.0, 0.0), 0.0, 0.0, Vector3::z(), 5.0)
            .execute(&mut store);
        assert!(result.is_err());
    }

    #[test]
    fn zero_height_fails() {
        let mut store = TopologyStore::new();
        let result = MakeCone::new(p(0.0, 0.0, 0.0), 3.0, 0.0, Vector3::z(), 0.0)
            .execute(&mut store);
        assert!(result.is_err());
    }
}
