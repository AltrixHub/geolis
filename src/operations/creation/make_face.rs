use crate::error::{OperationError, Result, TopologyError};
use crate::geometry::surface::Plane;
use crate::math::{Point3, Vector3, TOLERANCE};
use crate::topology::{FaceData, FaceId, FaceSurface, TopologyStore, WireId};

/// Creates a face from a wire boundary and a surface.
pub struct MakeFace {
    outer_wire: WireId,
    inner_wires: Vec<WireId>,
}

impl MakeFace {
    /// Creates a new `MakeFace` operation.
    #[must_use]
    pub fn new(outer_wire: WireId, inner_wires: Vec<WireId>) -> Self {
        Self {
            outer_wire,
            inner_wires,
        }
    }

    /// Executes the operation, creating the face in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if any wire is not closed, points are not coplanar,
    /// or all points are collinear.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<FaceId> {
        // Validate outer wire is closed
        validate_wire_closed(store, self.outer_wire)?;

        // Validate inner wires are closed
        for &wire_id in &self.inner_wires {
            validate_wire_closed(store, wire_id)?;
        }

        // Collect points from outer wire
        let outer_points = collect_wire_points(store, self.outer_wire)?;

        // Compute plane from outer wire points
        let plane = compute_plane_from_points(&outer_points)?;

        // Validate outer wire coplanarity
        validate_coplanarity(&plane, &outer_points)?;

        // Validate inner wire coplanarity
        for &wire_id in &self.inner_wires {
            let inner_points = collect_wire_points(store, wire_id)?;
            validate_coplanarity(&plane, &inner_points)?;
        }

        // Create face (same_sense = true because we construct the plane from the wire)
        let face_id = store.add_face(FaceData {
            surface: FaceSurface::Plane(plane),
            outer_wire: self.outer_wire,
            inner_wires: self.inner_wires.clone(),
            same_sense: true,
        });

        Ok(face_id)
    }
}

/// Validates that a wire exists and is closed.
fn validate_wire_closed(store: &TopologyStore, wire_id: WireId) -> Result<()> {
    let wire = store.wire(wire_id)?;
    if !wire.is_closed {
        return Err(TopologyError::WireNotClosed.into());
    }
    Ok(())
}

/// Collects vertex positions from a wire in traversal order.
fn collect_wire_points(store: &TopologyStore, wire_id: WireId) -> Result<Vec<Point3>> {
    let edges = store.wire(wire_id)?.edges.clone();
    let mut points = Vec::with_capacity(edges.len());

    for oe in &edges {
        let edge = store.edge(oe.edge)?;
        let vertex_id = if oe.forward { edge.start } else { edge.end };
        let vertex = store.vertex(vertex_id)?;
        points.push(vertex.point);
    }

    Ok(points)
}

/// Computes a plane from a set of points using Newell's method.
fn compute_plane_from_points(points: &[Point3]) -> Result<Plane> {
    let n = points.len();
    if n < 3 {
        return Err(OperationError::Failed(
            "at least 3 points are required to define a plane".into(),
        )
        .into());
    }

    // Newell's method for normal estimation
    let mut normal = Vector3::new(0.0, 0.0, 0.0);
    for i in 0..n {
        let curr = &points[i];
        let next = &points[(i + 1) % n];
        normal.x += (curr.y - next.y) * (curr.z + next.z);
        normal.y += (curr.z - next.z) * (curr.x + next.x);
        normal.z += (curr.x - next.x) * (curr.y + next.y);
    }

    let len = normal.norm();
    if len < TOLERANCE {
        return Err(
            OperationError::Failed("all points are collinear, cannot define a plane".into())
                .into(),
        );
    }

    // Centroid as origin
    #[allow(clippy::cast_precision_loss)]
    let inv_n = 1.0 / n as f64;
    let centroid = Point3::new(
        points.iter().map(|p| p.x).sum::<f64>() * inv_n,
        points.iter().map(|p| p.y).sum::<f64>() * inv_n,
        points.iter().map(|p| p.z).sum::<f64>() * inv_n,
    );

    Plane::from_normal(centroid, normal)
}

/// Validates that all points lie within `TOLERANCE` of the plane.
fn validate_coplanarity(plane: &Plane, points: &[Point3]) -> Result<()> {
    let origin = plane.origin();
    let normal = plane.plane_normal();
    for (i, p) in points.iter().enumerate() {
        let dist = (p - origin).dot(normal).abs();
        if dist > TOLERANCE {
            return Err(OperationError::InvalidInput(format!(
                "point {i} is not coplanar (distance = {dist})"
            ))
            .into());
        }
    }
    Ok(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::operations::creation::MakeWire;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    fn make_closed_wire(store: &mut TopologyStore, points: Vec<Point3>) -> crate::topology::WireId {
        MakeWire::new(points, true).execute(store).unwrap()
    }

    #[test]
    fn square_xy_plane_normal_is_z() {
        let mut store = TopologyStore::new();
        let pts = vec![p(0.0, 0.0, 0.0), p(4.0, 0.0, 0.0), p(4.0, 4.0, 0.0), p(0.0, 4.0, 0.0)];
        let wire = make_closed_wire(&mut store, pts);
        let face_id = MakeFace::new(wire, vec![]).execute(&mut store).unwrap();
        let face = store.face(face_id).unwrap();
        let FaceSurface::Plane(plane) = &face.surface;
        let n = plane.plane_normal();
        assert!(n.z.abs() > 0.99);
    }

    #[test]
    fn triangle_creates_face() {
        let mut store = TopologyStore::new();
        let pts = vec![p(0.0, 0.0, 0.0), p(3.0, 0.0, 0.0), p(1.5, 2.0, 0.0)];
        let wire = make_closed_wire(&mut store, pts);
        let face_id = MakeFace::new(wire, vec![]).execute(&mut store).unwrap();
        let face = store.face(face_id).unwrap();
        assert!(face.same_sense);
        assert!(face.inner_wires.is_empty());
    }

    #[test]
    fn face_with_inner_wire() {
        let mut store = TopologyStore::new();
        let outer = vec![
            p(0.0, 0.0, 0.0), p(10.0, 0.0, 0.0),
            p(10.0, 10.0, 0.0), p(0.0, 10.0, 0.0),
        ];
        let inner = vec![
            p(2.0, 2.0, 0.0), p(8.0, 2.0, 0.0),
            p(8.0, 8.0, 0.0), p(2.0, 8.0, 0.0),
        ];
        let outer_wire = make_closed_wire(&mut store, outer);
        let inner_wire = make_closed_wire(&mut store, inner);
        let face_id = MakeFace::new(outer_wire, vec![inner_wire])
            .execute(&mut store)
            .unwrap();
        let face = store.face(face_id).unwrap();
        assert_eq!(face.inner_wires.len(), 1);
    }

    #[test]
    fn open_wire_fails() {
        let mut store = TopologyStore::new();
        let pts = vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0)];
        let wire = MakeWire::new(pts, false).execute(&mut store).unwrap();
        let result = MakeFace::new(wire, vec![]).execute(&mut store);
        assert!(result.is_err());
    }

    #[test]
    fn collinear_points_fail() {
        let mut store = TopologyStore::new();
        let pts = vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(2.0, 0.0, 0.0)];
        let wire = make_closed_wire(&mut store, pts);
        let result = MakeFace::new(wire, vec![]).execute(&mut store);
        assert!(result.is_err());
    }

    #[test]
    fn non_coplanar_points_fail() {
        let mut store = TopologyStore::new();
        let pts = vec![
            p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0),
            p(1.0, 1.0, 0.0), p(0.0, 1.0, 1.0),
        ];
        let wire = make_closed_wire(&mut store, pts);
        let result = MakeFace::new(wire, vec![]).execute(&mut store);
        assert!(result.is_err());
    }
}
