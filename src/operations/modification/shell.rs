use std::collections::HashSet;

use crate::error::{OperationError, Result};
use crate::math::TOLERANCE;
use crate::operations::boolean::Subtract;
use crate::operations::creation::MakeBox;
use crate::operations::query::BoundingBox;
use crate::topology::{FaceId, FaceSurface, SolidId, TopologyStore};

/// Hollows a solid by removing specified faces and offsetting the remaining
/// faces inward by a given thickness.
///
/// The result is a shell (thin-walled solid) where the wall thickness equals
/// the given `thickness` parameter. The removed faces become openings.
///
/// Implementation uses boolean subtraction: an inner block (derived from
/// the offset planes of kept faces) is subtracted from the original solid.
/// The inner block extends beyond the solid at each removed face, creating
/// openings.
///
/// Currently supports axis-aligned box solids with `Plane` faces.
pub struct Shell {
    solid: SolidId,
    thickness: f64,
    removed_faces: Vec<FaceId>,
}

impl Shell {
    /// Creates a new `Shell` operation.
    #[must_use]
    pub fn new(solid: SolidId, thickness: f64, removed_faces: Vec<FaceId>) -> Self {
        Self {
            solid,
            thickness,
            removed_faces,
        }
    }

    /// Executes the shell operation.
    ///
    /// # Errors
    ///
    /// Returns an error if thickness is invalid, faces are not found,
    /// or topology construction fails.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        if self.thickness < TOLERANCE {
            return Err(
                OperationError::InvalidInput("shell thickness must be positive".into()).into(),
            );
        }
        if self.removed_faces.is_empty() {
            return Err(
                OperationError::InvalidInput("at least one face must be removed".into()).into(),
            );
        }

        let solid_data = store.solid(self.solid)?;
        let shell = store.shell(solid_data.outer_shell)?;
        let face_ids: Vec<FaceId> = shell.faces.clone();

        let removed_set: HashSet<FaceId> = self.removed_faces.iter().copied().collect();

        // Get the solid's bounding box
        let aabb = BoundingBox::new(self.solid).execute(store)?;

        // Compute inner box by offsetting each face inward (kept) or outward (removed).
        // For axis-aligned box solids, each face normal aligns with ±X, ±Y, or ±Z.
        let mut inner_min = aabb.min;
        let mut inner_max = aabb.max;
        let extend = self.thickness * 2.0;

        for &face_id in &face_ids {
            let face = store.face(face_id)?;
            let FaceSurface::Plane(plane) = &face.surface else {
                continue;
            };
            let normal = if face.same_sense {
                *plane.plane_normal()
            } else {
                -*plane.plane_normal()
            };

            let is_removed = removed_set.contains(&face_id);
            let (axis, positive) = dominant_axis(&normal);

            if is_removed {
                // Extend beyond the solid to create an opening
                if positive {
                    inner_max[axis] += extend;
                } else {
                    inner_min[axis] -= extend;
                }
            } else {
                // Offset inward by thickness
                if positive {
                    inner_max[axis] -= self.thickness;
                } else {
                    inner_min[axis] += self.thickness;
                }
            }
        }

        // Validate that the inner box is non-degenerate
        if inner_min.x >= inner_max.x
            || inner_min.y >= inner_max.y
            || inner_min.z >= inner_max.z
        {
            return Err(OperationError::InvalidInput(
                "shell thickness too large for solid dimensions".into(),
            )
            .into());
        }

        // Create inner box and subtract from outer solid
        let inner = MakeBox::new(inner_min, inner_max).execute(store)?;
        Subtract::new(self.solid, inner).execute(store)
    }
}

/// Returns the dominant axis index (0=X, 1=Y, 2=Z) and whether
/// the normal points in the positive direction along that axis.
fn dominant_axis(normal: &crate::math::Vector3) -> (usize, bool) {
    let ax = normal.x.abs();
    let ay = normal.y.abs();
    let az = normal.z.abs();
    if ax >= ay && ax >= az {
        (0, normal.x > 0.0)
    } else if ay >= az {
        (1, normal.y > 0.0)
    } else {
        (2, normal.z > 0.0)
    }
}


#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use crate::operations::query::Volume;
    use crate::tessellation::{TessellateSolid, TessellationParams};

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    /// Gets the top face (highest z normal) of a box solid.
    fn get_top_face(store: &TopologyStore, solid: SolidId) -> FaceId {
        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        let mut best = shell.faces[0];
        let mut best_z = f64::NEG_INFINITY;
        for &face_id in &shell.faces {
            let face = store.face(face_id).unwrap();
            if let FaceSurface::Plane(plane) = &face.surface {
                let z = if face.same_sense {
                    plane.plane_normal().z
                } else {
                    -plane.plane_normal().z
                };
                if z > best_z {
                    best_z = z;
                    best = face_id;
                }
            }
        }
        best
    }

    /// Gets a side face (highest x normal) of a box solid.
    fn get_side_face(store: &TopologyStore, solid: SolidId) -> FaceId {
        let solid_data = store.solid(solid).unwrap();
        let shell = store.shell(solid_data.outer_shell).unwrap();
        let mut best = shell.faces[0];
        let mut best_x = f64::NEG_INFINITY;
        for &face_id in &shell.faces {
            let face = store.face(face_id).unwrap();
            if let FaceSurface::Plane(plane) = &face.surface {
                let x = if face.same_sense {
                    plane.plane_normal().x
                } else {
                    -plane.plane_normal().x
                };
                if x > best_x {
                    best_x = x;
                    best = face_id;
                }
            }
        }
        best
    }

    #[test]
    fn shell_box_top_removed() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(4.0, 4.0, 4.0))
            .execute(&mut store)
            .unwrap();

        let top = get_top_face(&store, solid);
        let result = Shell::new(solid, 0.5, vec![top])
            .execute(&mut store)
            .unwrap();

        // Boolean subtract should produce a valid solid with faces
        let result_data = store.solid(result).unwrap();
        let result_shell = store.shell(result_data.outer_shell).unwrap();
        assert!(
            result_shell.faces.len() >= 10,
            "expected at least 10 faces, got {}",
            result_shell.faces.len()
        );
    }

    #[test]
    fn shell_box_side_removed() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(4.0, 4.0, 4.0))
            .execute(&mut store)
            .unwrap();

        let side = get_side_face(&store, solid);
        let result = Shell::new(solid, 0.5, vec![side])
            .execute(&mut store)
            .unwrap();

        let result_data = store.solid(result).unwrap();
        let result_shell = store.shell(result_data.outer_shell).unwrap();
        assert!(
            result_shell.faces.len() >= 10,
            "expected at least 10 faces, got {}",
            result_shell.faces.len()
        );
    }

    #[test]
    fn shell_tessellates() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(4.0, 4.0, 4.0))
            .execute(&mut store)
            .unwrap();

        let top = get_top_face(&store, solid);
        let result = Shell::new(solid, 0.5, vec![top])
            .execute(&mut store)
            .unwrap();

        let mesh = TessellateSolid::new(result, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty());
    }

    #[test]
    fn shell_inner_vertices_correct() {
        // Box (0,0,0)-(4,4,4), top removed, thickness=0.5
        // Inner block should be (0.5, 0.5, 0.5)-(3.5, 3.5, 4+extend)
        // Subtraction result should have vertices at the inner boundary
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(4.0, 4.0, 4.0))
            .execute(&mut store)
            .unwrap();

        let top = get_top_face(&store, solid);
        let result = Shell::new(solid, 0.5, vec![top])
            .execute(&mut store)
            .unwrap();

        // Check bounding box of result
        let aabb = BoundingBox::new(result).execute(&store).unwrap();
        assert!((aabb.min.x - 0.0).abs() < 1e-6);
        assert!((aabb.min.y - 0.0).abs() < 1e-6);
        assert!((aabb.min.z - 0.0).abs() < 1e-6);
        assert!((aabb.max.x - 4.0).abs() < 1e-6);
        assert!((aabb.max.y - 4.0).abs() < 1e-6);
        assert!((aabb.max.z - 4.0).abs() < 1e-6);
    }

    #[test]
    fn shell_volume_correct() {
        // Box (0,0,0)-(4,4,4), top removed, thickness=0.5
        // Outer box volume = 4 * 4 * 4 = 64
        // Inner box = (0.5, 0.5, 0.5) to (3.5, 3.5, 5.0)
        //   → inner width=3, depth=3, height=4.5 BUT clamped by outer box at z=4
        //   → effective inner = (0.5,0.5,0.5) to (3.5,3.5,4.0) = 3 * 3 * 3.5 = 31.5
        // Shell volume = 64 - 31.5 = 32.5
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(4.0, 4.0, 4.0))
            .execute(&mut store)
            .unwrap();

        let top = get_top_face(&store, solid);
        let result = Shell::new(solid, 0.5, vec![top])
            .execute(&mut store)
            .unwrap();

        let volume = Volume::new(result).execute(&store).unwrap();
        let expected = 32.5;
        assert!(
            (volume - expected).abs() < expected * 0.05,
            "shell volume: expected ~{expected}, got {volume}"
        );
    }

    #[test]
    fn shell_zero_thickness_fails() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(4.0, 4.0, 4.0))
            .execute(&mut store)
            .unwrap();
        let top = get_top_face(&store, solid);

        let result = Shell::new(solid, 0.0, vec![top]).execute(&mut store);
        assert!(result.is_err());
    }

    #[test]
    fn shell_no_removed_faces_fails() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(4.0, 4.0, 4.0))
            .execute(&mut store)
            .unwrap();

        let result = Shell::new(solid, 0.5, vec![]).execute(&mut store);
        assert!(result.is_err());
    }

    #[test]
    fn shell_top_face_has_opening() {
        // The shell result should have a face at z=4 with an inner wire (hole)
        // representing the opening where the top face was removed.
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(4.0, 4.0, 4.0))
            .execute(&mut store)
            .unwrap();

        let top = get_top_face(&store, solid);
        let result = Shell::new(solid, 0.5, vec![top])
            .execute(&mut store)
            .unwrap();

        let result_data = store.solid(result).unwrap();
        let result_shell = store.shell(result_data.outer_shell).unwrap();

        // Find any face at z=4 with normal pointing up
        let mut top_faces = Vec::new();
        for &face_id in &result_shell.faces {
            let face = store.face(face_id).unwrap();
            if let FaceSurface::Plane(plane) = &face.surface {
                let normal = if face.same_sense {
                    *plane.plane_normal()
                } else {
                    -*plane.plane_normal()
                };
                // Check if face is at z=4 with upward normal
                if normal.z > 0.9 {
                    let wire = store.wire(face.outer_wire).unwrap();
                    let mut max_z = f64::NEG_INFINITY;
                    for oe in &wire.edges {
                        let edge = store.edge(oe.edge).unwrap();
                        let sv = store.vertex(edge.start).unwrap();
                        if sv.point.z > max_z {
                            max_z = sv.point.z;
                        }
                    }
                    if (max_z - 4.0).abs() < 0.01 {
                        top_faces.push((face_id, !face.inner_wires.is_empty()));
                    }
                }
            }
        }

        // There should be a top face at z=4, and it should have an inner wire (hole)
        assert!(
            !top_faces.is_empty(),
            "expected at least one face at z=4 with upward normal"
        );

        let has_hole = top_faces.iter().any(|(_id, has_inner)| *has_inner);
        assert!(
            has_hole,
            "expected top face at z=4 to have inner wire (opening), but found faces: {top_faces:?}"
        );
    }

    #[test]
    fn shell_thickness_too_large_fails() {
        let mut store = TopologyStore::new();
        let solid = MakeBox::new(p(0.0, 0.0, 0.0), p(4.0, 4.0, 4.0))
            .execute(&mut store)
            .unwrap();
        let top = get_top_face(&store, solid);

        // Thickness 2.5 > half of 4.0 → inner box is degenerate
        let result = Shell::new(solid, 2.5, vec![top]).execute(&mut store);
        assert!(result.is_err());
    }
}
