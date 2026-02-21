use crate::error::{OperationError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};
use crate::operations::creation::{MakeFace, MakeWire};
use crate::topology::{FaceId, FaceSurface, TopologyStore};

/// Offsets a face along its normal direction by a given distance.
///
/// Creates a new face that is a copy of the original, displaced by `distance`
/// along the face's outward normal direction.
///
/// - **Plane**: translates origin and boundary by `normal * distance`
/// - **Cylinder**: changes radius by `distance` (positive = outward)
/// - **Sphere**: changes radius by `distance`
/// - **Cone**: adjusts apex position along axis
///
/// Currently fully supports Plane faces. Curved surface support is
/// limited to translating boundary vertices.
pub struct FaceOffset {
    face: FaceId,
    distance: f64,
}

impl FaceOffset {
    /// Creates a new `FaceOffset` operation.
    #[must_use]
    pub fn new(face: FaceId, distance: f64) -> Self {
        Self { face, distance }
    }

    /// Executes the offset, returning a new face ID.
    ///
    /// # Errors
    ///
    /// Returns an error if the face is not found or the offset would
    /// create a degenerate face (e.g., negative radius).
    pub fn execute(&self, store: &mut TopologyStore) -> Result<FaceId> {
        let face = store.face(self.face)?;
        let same_sense = face.same_sense;
        let outer_wire_id = face.outer_wire;

        match &face.surface {
            FaceSurface::Plane(plane) => {
                let normal = if same_sense {
                    *plane.plane_normal()
                } else {
                    -*plane.plane_normal()
                };
                let offset = normal * self.distance;

                let outer_points = collect_wire_points(store, outer_wire_id)?;
                let offset_points: Vec<Point3> =
                    outer_points.iter().map(|p| p + offset).collect();

                let wire = MakeWire::new(offset_points, true).execute(store)?;
                MakeFace::new(wire, vec![]).execute(store)
            }
            FaceSurface::Cylinder(cyl) => {
                let new_radius = cyl.radius() + self.distance;
                if new_radius < TOLERANCE {
                    return Err(OperationError::InvalidInput(
                        "cylinder offset would produce zero or negative radius".into(),
                    )
                    .into());
                }

                // Offset boundary vertices radially
                let outer_points = collect_wire_points(store, outer_wire_id)?;
                let offset_points: Vec<Point3> = outer_points
                    .iter()
                    .map(|p| offset_radially(p, cyl.center(), cyl.axis(), self.distance))
                    .collect();

                let wire = MakeWire::new(offset_points, true).execute(store)?;
                MakeFace::new(wire, vec![]).execute(store)
            }
            FaceSurface::Sphere(sph) => {
                let new_radius = sph.radius() + self.distance;
                if new_radius < TOLERANCE {
                    return Err(OperationError::InvalidInput(
                        "sphere offset would produce zero or negative radius".into(),
                    )
                    .into());
                }

                // Offset boundary vertices radially from center
                let outer_points = collect_wire_points(store, outer_wire_id)?;
                let offset_points: Vec<Point3> = outer_points
                    .iter()
                    .map(|p| {
                        let dp = p - sph.center();
                        let len = dp.norm();
                        if len < TOLERANCE {
                            *p
                        } else {
                            *sph.center() + dp * (new_radius / len)
                        }
                    })
                    .collect();

                let wire = MakeWire::new(offset_points, true).execute(store)?;
                MakeFace::new(wire, vec![]).execute(store)
            }
            FaceSurface::Cone(cone) => {
                // For cone, offset boundary vertices radially
                let outer_points = collect_wire_points(store, outer_wire_id)?;
                let offset_points: Vec<Point3> = outer_points
                    .iter()
                    .map(|p| offset_radially(p, cone.apex(), cone.axis(), self.distance))
                    .collect();

                let wire = MakeWire::new(offset_points, true).execute(store)?;
                MakeFace::new(wire, vec![]).execute(store)
            }
            FaceSurface::Torus(_) => {
                // For torus, use simple normal offset on boundary vertices
                let normal = face.surface.normal_at_centroid(store, outer_wire_id);
                let offset = normal * self.distance;

                let outer_points = collect_wire_points(store, outer_wire_id)?;
                let offset_points: Vec<Point3> =
                    outer_points.iter().map(|p| p + offset).collect();

                let wire = MakeWire::new(offset_points, true).execute(store)?;
                MakeFace::new(wire, vec![]).execute(store)
            }
        }
    }
}

/// Offsets a point radially from an axis.
fn offset_radially(point: &Point3, axis_point: &Point3, axis: &Vector3, distance: f64) -> Point3 {
    let dp = point - axis_point;
    let axis_proj = dp.dot(axis);
    let foot = axis_point + axis * axis_proj;
    let radial = point - foot;
    let radial_len = radial.norm();

    if radial_len < TOLERANCE {
        // Point is on axis; offset along ref direction (can't determine radial)
        *point
    } else {
        let radial_dir = radial / radial_len;
        *point + radial_dir * distance
    }
}

/// Collects vertex positions from a wire in traversal order.
fn collect_wire_points(
    store: &TopologyStore,
    wire_id: crate::topology::WireId,
) -> Result<Vec<Point3>> {
    let edges = store.wire(wire_id)?.edges.clone();
    let mut points = Vec::with_capacity(edges.len());
    for oe in &edges {
        let edge = store.edge(oe.edge)?;
        let vertex_id = if oe.forward { edge.start } else { edge.end };
        points.push(store.vertex(vertex_id)?.point);
    }
    Ok(points)
}

/// Extension trait for `FaceSurface` to compute an approximate normal at centroid.
trait FaceSurfaceExt {
    fn normal_at_centroid(&self, store: &TopologyStore, wire_id: crate::topology::WireId)
        -> Vector3;
}

impl FaceSurfaceExt for FaceSurface {
    fn normal_at_centroid(
        &self,
        store: &TopologyStore,
        wire_id: crate::topology::WireId,
    ) -> Vector3 {
        use crate::geometry::surface::Surface;
        let points = collect_wire_points(store, wire_id).unwrap_or_default();
        if points.is_empty() {
            return Vector3::z();
        }
        let centroid = points.iter().fold(Point3::origin(), |acc, p| {
            Point3::new(acc.x + p.x, acc.y + p.y, acc.z + p.z)
        });
        #[allow(clippy::cast_precision_loss)]
        let n = points.len() as f64;
        let centroid = Point3::new(centroid.x / n, centroid.y / n, centroid.z / n);

        match self {
            FaceSurface::Torus(t) => {
                let (u, v) = t.inverse(&centroid);
                t.normal(u, v).unwrap_or(Vector3::z())
            }
            _ => Vector3::z(),
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::operations::creation::{MakeFace, MakeWire};

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    fn make_xy_face(store: &mut TopologyStore) -> FaceId {
        let wire = MakeWire::new(
            vec![p(0.0, 0.0, 0.0), p(4.0, 0.0, 0.0), p(4.0, 4.0, 0.0), p(0.0, 4.0, 0.0)],
            true,
        )
        .execute(store)
        .unwrap();
        MakeFace::new(wire, vec![]).execute(store).unwrap()
    }

    #[test]
    fn plane_offset_positive() {
        let mut store = TopologyStore::new();
        let face = make_xy_face(&mut store);

        let new_face = FaceOffset::new(face, 1.0).execute(&mut store).unwrap();

        let new_face_data = store.face(new_face).unwrap();
        let wire = store.wire(new_face_data.outer_wire).unwrap();
        // All vertices should be at z = 1.0
        for oe in &wire.edges {
            let edge = store.edge(oe.edge).unwrap();
            let start = store.vertex(edge.start).unwrap();
            assert!(
                (start.point.z - 1.0).abs() < 1e-10,
                "expected z=1.0, got z={}",
                start.point.z
            );
        }
    }

    #[test]
    fn plane_offset_negative() {
        let mut store = TopologyStore::new();
        let face = make_xy_face(&mut store);

        let new_face = FaceOffset::new(face, -2.0).execute(&mut store).unwrap();

        let new_face_data = store.face(new_face).unwrap();
        let wire = store.wire(new_face_data.outer_wire).unwrap();
        for oe in &wire.edges {
            let edge = store.edge(oe.edge).unwrap();
            let start = store.vertex(edge.start).unwrap();
            assert!(
                (start.point.z - (-2.0)).abs() < 1e-10,
                "expected z=-2.0, got z={}",
                start.point.z
            );
        }
    }
}
