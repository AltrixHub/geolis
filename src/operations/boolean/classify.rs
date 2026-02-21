use crate::error::Result;
use crate::math::intersect_3d::{line_plane_intersect, LinePlaneRelation};
use crate::math::polygon_3d::point_in_polygon_3d;
use crate::math::{Point3, Vector3, TOLERANCE};
use crate::topology::{FaceSurface, SolidId, TopologyStore};

/// Classification of a point relative to a solid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointClassification {
    Inside,
    Outside,
    OnBoundary,
}

/// Classifies a point as inside, outside, or on the boundary of a solid.
///
/// Uses ray casting: shoots a ray from the point and counts face crossings.
/// Odd crossings = inside, even = outside. If the ray is degenerate
/// (hits an edge/vertex), retries with alternative directions.
///
/// # Errors
///
/// Returns an error if the solid or its topology cannot be read.
pub fn classify_point_in_solid(
    point: &Point3,
    solid_id: SolidId,
    store: &TopologyStore,
) -> Result<PointClassification> {
    let solid = store.solid(solid_id)?;
    let shell = store.shell(solid.outer_shell)?;

    // Collect face polygons and planes
    let faces = collect_face_data(store, &shell.faces)?;

    // Try up to 3 ray directions
    let directions = [
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    ];

    for dir in &directions {
        if let RayCastResult::Clear(classification) = ray_cast_classify(point, dir, &faces) {
            return Ok(classification);
        }
    }

    // All directions degenerate — very unlikely, treat as outside
    Ok(PointClassification::Outside)
}

struct FaceInfo {
    polygon: Vec<Point3>,
    plane: crate::geometry::surface::Plane,
    #[allow(dead_code)]
    outward_normal: Vector3,
}

fn collect_face_data(
    store: &TopologyStore,
    face_ids: &[crate::topology::FaceId],
) -> Result<Vec<FaceInfo>> {
    let mut faces = Vec::with_capacity(face_ids.len());
    for &face_id in face_ids {
        let face = store.face(face_id)?;
        let FaceSurface::Plane(ref plane) = face.surface else {
            todo!("Point classification for non-planar faces")
        };

        let wire = store.wire(face.outer_wire)?;
        let mut polygon = Vec::with_capacity(wire.edges.len());
        for oe in &wire.edges {
            let edge = store.edge(oe.edge)?;
            let vid = if oe.forward { edge.start } else { edge.end };
            let vertex = store.vertex(vid)?;
            polygon.push(vertex.point);
        }

        let outward_normal = if face.same_sense {
            *plane.plane_normal()
        } else {
            -plane.plane_normal()
        };

        faces.push(FaceInfo {
            polygon,
            plane: plane.clone(),
            outward_normal,
        });
    }
    Ok(faces)
}

enum RayCastResult {
    Clear(PointClassification),
    Degenerate,
}

fn ray_cast_classify(point: &Point3, dir: &Vector3, faces: &[FaceInfo]) -> RayCastResult {
    let mut crossings = 0u32;
    let boundary_tol = TOLERANCE * 10.0;

    for face in faces {
        match line_plane_intersect(point, dir, &face.plane) {
            LinePlaneRelation::Point { point: hit, t } => {
                // Only count forward intersections (t > 0)
                if t < boundary_tol {
                    if t > -boundary_tol {
                        // Point is very close to a face plane — check if on boundary
                        if point_in_polygon_3d(point, &face.polygon, &face.plane) {
                            return RayCastResult::Clear(PointClassification::OnBoundary);
                        }
                    }
                    continue;
                }

                // Check if hit point is inside the face polygon
                if !point_in_polygon_3d(&hit, &face.polygon, &face.plane) {
                    continue;
                }

                // Check if hit point is near a polygon edge (degenerate)
                if is_near_polygon_edge(&hit, &face.polygon, &face.plane) {
                    return RayCastResult::Degenerate;
                }

                crossings += 1;
            }
            LinePlaneRelation::OnPlane => {
                // Ray lies in the face plane — degenerate
                return RayCastResult::Degenerate;
            }
            LinePlaneRelation::Parallel => {
                // No intersection
            }
        }
    }

    if crossings % 2 == 1 {
        RayCastResult::Clear(PointClassification::Inside)
    } else {
        RayCastResult::Clear(PointClassification::Outside)
    }
}

/// Check if a point is near any edge of the polygon (within tolerance).
fn is_near_polygon_edge(point: &Point3, polygon: &[Point3], _plane: &crate::geometry::surface::Plane) -> bool {
    let n = polygon.len();
    let edge_tol = TOLERANCE * 100.0;
    for i in 0..n {
        let a = &polygon[i];
        let b = &polygon[(i + 1) % n];
        let ab = b - a;
        let ab_len_sq = ab.dot(&ab);
        if ab_len_sq < TOLERANCE * TOLERANCE {
            continue;
        }
        let ap = point - a;
        let t = ap.dot(&ab) / ab_len_sq;
        if t < -edge_tol || t > 1.0 + edge_tol {
            continue;
        }
        let closest = a + ab * t.clamp(0.0, 1.0);
        let dist = (point - closest).norm();
        if dist < edge_tol {
            return true;
        }
    }
    false
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::operations::creation::{MakeFace, MakeWire};
    use crate::operations::shaping::Extrude;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    fn make_box(store: &mut TopologyStore) -> SolidId {
        let pts = vec![
            p(0.0, 0.0, 0.0),
            p(2.0, 0.0, 0.0),
            p(2.0, 2.0, 0.0),
            p(0.0, 2.0, 0.0),
        ];
        let wire = MakeWire::new(pts, true).execute(store).unwrap();
        let face = MakeFace::new(wire, vec![]).execute(store).unwrap();
        Extrude::new(face, Vector3::new(0.0, 0.0, 2.0))
            .execute(store)
            .unwrap()
    }

    #[test]
    fn center_is_inside() {
        let mut store = TopologyStore::new();
        let solid = make_box(&mut store);
        let result = classify_point_in_solid(&p(1.0, 1.0, 1.0), solid, &store).unwrap();
        assert_eq!(result, PointClassification::Inside);
    }

    #[test]
    fn far_point_is_outside() {
        let mut store = TopologyStore::new();
        let solid = make_box(&mut store);
        let result = classify_point_in_solid(&p(10.0, 10.0, 10.0), solid, &store).unwrap();
        assert_eq!(result, PointClassification::Outside);
    }

    #[test]
    fn point_near_face_is_outside() {
        let mut store = TopologyStore::new();
        let solid = make_box(&mut store);
        let result = classify_point_in_solid(&p(1.0, 1.0, -1.0), solid, &store).unwrap();
        assert_eq!(result, PointClassification::Outside);
    }

    #[test]
    fn point_on_face_is_boundary() {
        let mut store = TopologyStore::new();
        let solid = make_box(&mut store);
        let result = classify_point_in_solid(&p(1.0, 1.0, 0.0), solid, &store).unwrap();
        assert_eq!(result, PointClassification::OnBoundary);
    }

    #[test]
    fn point_just_inside() {
        let mut store = TopologyStore::new();
        let solid = make_box(&mut store);
        let result =
            classify_point_in_solid(&p(0.001, 0.001, 0.001), solid, &store).unwrap();
        assert_eq!(result, PointClassification::Inside);
    }
}
