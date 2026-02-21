use crate::error::{OperationError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};
use crate::operations::creation::{MakeFace, MakeSolid, MakeWire};
use crate::topology::{
    FaceId, FaceSurface, ShellData, SolidId, TopologyStore,
};

/// Hollows a solid by removing specified faces and offsetting the remaining
/// faces inward by a given thickness.
///
/// The result is a shell (thin-walled solid) with:
/// - The original outer faces (unchanged)
/// - Offset inner faces (moved inward by thickness)
/// - Side faces connecting the boundaries of removed faces
///
/// Currently supports Plane faces only.
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

        let solid = store.solid(self.solid)?;
        let shell = store.shell(solid.outer_shell)?;
        let face_ids: Vec<FaceId> = shell.faces.clone();

        // Separate faces into kept and removed
        let removed_set: std::collections::HashSet<FaceId> =
            self.removed_faces.iter().copied().collect();

        let mut all_result_faces: Vec<FaceId> = Vec::new();

        // For each kept face, create the outer face (copy) and inner face (offset inward)
        for &face_id in &face_ids {
            if removed_set.contains(&face_id) {
                continue;
            }

            let face = store.face(face_id)?;
            match &face.surface {
                FaceSurface::Plane(plane) => {
                    let plane = plane.clone();
                    let same_sense = face.same_sense;
                    let outer_wire_id = face.outer_wire;

                    // Get outer boundary points
                    let outer_points = collect_wire_points(store, outer_wire_id)?;

                    // Create outer face (original)
                    let outer_wire = MakeWire::new(outer_points.clone(), true).execute(store)?;
                    let outer_face = MakeFace::new(outer_wire, vec![]).execute(store)?;
                    all_result_faces.push(outer_face);

                    // Create inner face (offset inward)
                    let normal = if same_sense {
                        *plane.plane_normal()
                    } else {
                        -*plane.plane_normal()
                    };
                    // Inward = opposite of outward normal
                    let offset_dir = -normal * self.thickness;
                    let inner_points: Vec<Point3> =
                        outer_points.iter().map(|p| p + offset_dir).collect();

                    // Inner face has reversed winding (points inward)
                    let inner_reversed: Vec<Point3> = inner_points.iter().copied().rev().collect();
                    let inner_wire = MakeWire::new(inner_reversed, true).execute(store)?;
                    let inner_face = MakeFace::new(inner_wire, vec![]).execute(store)?;
                    all_result_faces.push(inner_face);
                }
                _ => {
                    // For non-plane faces, just copy the outer face for now
                    let outer_wire_id = face.outer_wire;
                    let outer_points = collect_wire_points(store, outer_wire_id)?;
                    let outer_wire = MakeWire::new(outer_points, true).execute(store)?;
                    let outer_face = MakeFace::new(outer_wire, vec![]).execute(store)?;
                    all_result_faces.push(outer_face);
                }
            }
        }

        // For each removed face, create side faces connecting outer and inner edges
        for &face_id in &self.removed_faces {
            let face = store.face(face_id)?;
            let outer_wire_id = face.outer_wire;
            let outer_points = collect_wire_points(store, outer_wire_id)?;

            let same_sense = face.same_sense;
            let normal = match &face.surface {
                FaceSurface::Plane(plane) => {
                    if same_sense {
                        *plane.plane_normal()
                    } else {
                        -*plane.plane_normal()
                    }
                }
                _ => Vector3::z(), // fallback
            };

            let offset_dir = -normal * self.thickness;
            let inner_points: Vec<Point3> =
                outer_points.iter().map(|p| p + offset_dir).collect();

            // Create side faces (quads connecting outer[i]→outer[i+1]→inner[i+1]→inner[i])
            let n = outer_points.len();
            for i in 0..n {
                let j = (i + 1) % n;
                let quad = vec![
                    outer_points[i],
                    outer_points[j],
                    inner_points[j],
                    inner_points[i],
                ];
                let wire = MakeWire::new(quad, true).execute(store)?;
                let side_face = MakeFace::new(wire, vec![]).execute(store)?;
                all_result_faces.push(side_face);
            }
        }

        let shell_id = store.add_shell(ShellData {
            faces: all_result_faces,
            is_closed: true,
        });
        MakeSolid::new(shell_id, vec![]).execute(store)
    }
}

/// Collects vertex positions from a wire in traversal order.
fn collect_wire_points(store: &TopologyStore, wire_id: crate::topology::WireId) -> Result<Vec<Point3>> {
    let edges = store.wire(wire_id)?.edges.clone();
    let mut points = Vec::with_capacity(edges.len());
    for oe in &edges {
        let edge = store.edge(oe.edge)?;
        let vertex_id = if oe.forward { edge.start } else { edge.end };
        points.push(store.vertex(vertex_id)?.point);
    }
    Ok(points)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::operations::creation::MakeBox;
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

        let result_data = store.solid(result).unwrap();
        let result_shell = store.shell(result_data.outer_shell).unwrap();
        // 5 kept faces × 2 (outer + inner) + 4 side faces = 14
        assert_eq!(result_shell.faces.len(), 14);
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
        assert_eq!(result_shell.faces.len(), 14);
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
}
