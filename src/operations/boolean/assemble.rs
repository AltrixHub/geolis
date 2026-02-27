use std::collections::HashMap;

use crate::error::{OperationError, Result};
use crate::geometry::curve::Line;
use crate::math::{Point3, TOLERANCE};
use crate::operations::creation::{MakeFace, MakeSolid};
use crate::topology::{
    EdgeCurve, EdgeData, FaceId, OrientedEdge, ShellData, SolidId, TopologyStore, VertexData,
    VertexId, WireData,
};

use super::select::KeepDecision;
use super::split::FaceFragment;

/// Assembles a `BRep` solid from a set of face fragments.
///
/// Creates new topology (vertices, edges, wires, faces, shell, solid) from
/// the polygon boundaries of the fragments. Uses spatial hashing to merge
/// coincident vertices.
pub fn assemble_result(
    store: &mut TopologyStore,
    fragments: &[(FaceFragment, KeepDecision)],
) -> Result<SolidId> {
    let kept: Vec<&FaceFragment> = fragments
        .iter()
        .filter(|(_, decision)| *decision != KeepDecision::Discard)
        .map(|(frag, _)| frag)
        .collect();

    let decisions: Vec<KeepDecision> = fragments
        .iter()
        .filter(|(_, decision)| *decision != KeepDecision::Discard)
        .map(|(_, decision)| *decision)
        .collect();

    if kept.is_empty() {
        return Err(OperationError::Failed("boolean operation produced no faces".into()).into());
    }

    // Build vertex merger
    let mut merger = VertexMerger::new(TOLERANCE * 1000.0);

    // Create faces from each fragment
    let mut all_faces: Vec<FaceId> = Vec::with_capacity(kept.len());

    for (frag, &decision) in kept.iter().zip(decisions.iter()) {
        let boundary = if decision == KeepDecision::KeepFlipped {
            // Reverse winding order to flip the normal
            frag.boundary.iter().rev().copied().collect::<Vec<_>>()
        } else {
            frag.boundary.clone()
        };

        let inner_boundaries: Vec<Vec<Point3>> = if decision == KeepDecision::KeepFlipped {
            frag.inner_boundaries
                .iter()
                .map(|ib| ib.iter().rev().copied().collect())
                .collect()
        } else {
            frag.inner_boundaries.clone()
        };

        let face_id =
            create_face_from_polygon(store, &boundary, &inner_boundaries, &mut merger)?;
        all_faces.push(face_id);
    }

    // Create shell
    let shell_id = store.add_shell(ShellData {
        faces: all_faces,
        is_closed: true,
    });

    // Create solid
    MakeSolid::new(shell_id, vec![]).execute(store)
}

/// Creates a face from a polygon boundary with optional inner boundaries, reusing merged vertices.
fn create_face_from_polygon(
    store: &mut TopologyStore,
    boundary: &[Point3],
    inner_boundaries: &[Vec<Point3>],
    merger: &mut VertexMerger,
) -> Result<FaceId> {
    let n = boundary.len();
    if n < 3 {
        return Err(
            OperationError::Failed("face fragment has fewer than 3 vertices".into()).into(),
        );
    }

    // Get or create vertices
    let vertex_ids: Vec<VertexId> = boundary
        .iter()
        .map(|p| merger.get_or_create(store, p))
        .collect();

    // Create edges and oriented edges for the outer wire
    let mut oriented_edges = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let start = vertex_ids[i];
        let end = vertex_ids[j];
        let start_pt = boundary[i];
        let end_pt = boundary[j];

        let edge_id = create_line_edge(store, start, end, start_pt, end_pt)?;
        oriented_edges.push(OrientedEdge::new(edge_id, true));
    }

    // Create closed outer wire
    let wire_id = store.add_wire(WireData {
        edges: oriented_edges,
        is_closed: true,
    });

    // Create inner wires
    let mut inner_wire_ids = Vec::with_capacity(inner_boundaries.len());
    for inner in inner_boundaries {
        let m = inner.len();
        if m < 3 {
            continue;
        }

        let inner_vids: Vec<VertexId> = inner
            .iter()
            .map(|p| merger.get_or_create(store, p))
            .collect();

        let mut inner_edges = Vec::with_capacity(m);
        for i in 0..m {
            let j = (i + 1) % m;
            let edge_id =
                create_line_edge(store, inner_vids[i], inner_vids[j], inner[i], inner[j])?;
            inner_edges.push(OrientedEdge::new(edge_id, true));
        }

        let inner_wire_id = store.add_wire(WireData {
            edges: inner_edges,
            is_closed: true,
        });
        inner_wire_ids.push(inner_wire_id);
    }

    // Create face via MakeFace
    MakeFace::new(wire_id, inner_wire_ids).execute(store)
}

/// Creates a line edge between two vertices.
fn create_line_edge(
    store: &mut TopologyStore,
    start: VertexId,
    end: VertexId,
    start_point: Point3,
    end_point: Point3,
) -> Result<crate::topology::EdgeId> {
    let direction = end_point - start_point;
    let t_end = direction.norm();
    let line = Line::new(start_point, direction)?;
    Ok(store.add_edge(EdgeData {
        start,
        end,
        curve: EdgeCurve::Line(line),
        t_start: 0.0,
        t_end,
    }))
}

/// Spatial hash-based vertex merger.
///
/// Groups points by grid cell and merges vertices that are within `tolerance`
/// of each other.
struct VertexMerger {
    cell_size: f64,
    map: HashMap<(i64, i64, i64), Vec<(VertexId, Point3)>>,
}

impl VertexMerger {
    fn new(cell_size: f64) -> Self {
        Self {
            cell_size,
            map: HashMap::new(),
        }
    }

    #[allow(clippy::cast_possible_truncation)]
    fn cell_key(&self, p: &Point3) -> (i64, i64, i64) {
        let inv = 1.0 / self.cell_size;
        (
            (p.x * inv).floor() as i64,
            (p.y * inv).floor() as i64,
            (p.z * inv).floor() as i64,
        )
    }

    fn get_or_create(&mut self, store: &mut TopologyStore, point: &Point3) -> VertexId {
        let key = self.cell_key(point);

        // Search in neighboring cells (3x3x3) for a match
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let neighbor = (key.0 + dx, key.1 + dy, key.2 + dz);
                    if let Some(entries) = self.map.get(&neighbor) {
                        for &(vid, ref existing) in entries {
                            if (point - existing).norm() < self.cell_size {
                                return vid;
                            }
                        }
                    }
                }
            }
        }

        // No match found — create new vertex
        let vid = store.add_vertex(VertexData::new(*point));
        self.map.entry(key).or_default().push((vid, *point));
        vid
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::surface::Plane;
    use crate::math::Vector3;
    use crate::operations::boolean::split::SolidSource;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    fn make_fragment(boundary: Vec<Point3>, flip: bool) -> (FaceFragment, KeepDecision) {
        let plane = Plane::from_normal(p(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 1.0)).unwrap();
        let frag = FaceFragment {
            boundary,
            inner_boundaries: vec![],
            plane,
            same_sense: true,
            source_face: crate::topology::FaceId::default(),
            source: SolidSource::A,
        };
        let decision = if flip {
            KeepDecision::KeepFlipped
        } else {
            KeepDecision::Keep
        };
        (frag, decision)
    }

    #[test]
    fn assemble_single_box_fragments() {
        let mut store = TopologyStore::new();

        // Build 6 faces of a unit cube as fragments
        let fragments = vec![
            // Bottom (z=0)
            make_fragment(
                vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(0.0, 1.0, 0.0)],
                false,
            ),
            // Top (z=1)
            make_fragment(
                vec![p(0.0, 0.0, 1.0), p(1.0, 0.0, 1.0), p(1.0, 1.0, 1.0), p(0.0, 1.0, 1.0)],
                false,
            ),
            // Front (y=0)
            make_fragment(
                vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(1.0, 0.0, 1.0), p(0.0, 0.0, 1.0)],
                false,
            ),
            // Back (y=1)
            make_fragment(
                vec![p(0.0, 1.0, 0.0), p(1.0, 1.0, 0.0), p(1.0, 1.0, 1.0), p(0.0, 1.0, 1.0)],
                false,
            ),
            // Left (x=0)
            make_fragment(
                vec![p(0.0, 0.0, 0.0), p(0.0, 1.0, 0.0), p(0.0, 1.0, 1.0), p(0.0, 0.0, 1.0)],
                false,
            ),
            // Right (x=1)
            make_fragment(
                vec![p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(1.0, 1.0, 1.0), p(1.0, 0.0, 1.0)],
                false,
            ),
        ];

        let solid_id = assemble_result(&mut store, &fragments).unwrap();
        let solid = store.solid(solid_id).unwrap();
        let shell = store.shell(solid.outer_shell).unwrap();
        assert_eq!(shell.faces.len(), 6);
        assert!(shell.is_closed);
    }

    #[test]
    fn assemble_discards_discarded_fragments() {
        let mut store = TopologyStore::new();

        let fragments = vec![
            make_fragment(
                vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(0.0, 1.0, 0.0)],
                false,
            ),
            (
                FaceFragment {
                    boundary: vec![
                        p(2.0, 2.0, 0.0),
                        p(3.0, 2.0, 0.0),
                        p(3.0, 3.0, 0.0),
                        p(2.0, 3.0, 0.0),
                    ],
                    inner_boundaries: vec![],
                    plane: Plane::from_normal(p(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 1.0))
                        .unwrap(),
                    same_sense: true,
                    source_face: crate::topology::FaceId::default(),
                    source: SolidSource::B,
                },
                KeepDecision::Discard,
            ),
        ];

        let solid_id = assemble_result(&mut store, &fragments).unwrap();
        let solid = store.solid(solid_id).unwrap();
        let shell = store.shell(solid.outer_shell).unwrap();
        // Only the non-discarded fragment should be present
        assert_eq!(shell.faces.len(), 1);
    }
}
