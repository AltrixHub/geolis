use std::collections::{HashMap, HashSet};

use crate::error::{OperationError, Result};
use crate::geometry::curve::Line;
use crate::math::{Point3, TOLERANCE};
use crate::operations::creation::{MakeFace, MakeSolid};
use crate::topology::{
    EdgeCurve, EdgeData, FaceId, OrientedEdge, ShellData, SolidId, TopologyStore, VertexData,
    VertexId, WireData,
};

use super::select::KeepDecision;
use super::split::{newell_normal_3d, FaceFragment};

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

        let fragment_faces =
            create_face_from_polygon(store, &boundary, &inner_boundaries, &mut merger)?;
        all_faces.extend(fragment_faces);
    }

    if all_faces.is_empty() {
        return Err(OperationError::Failed(
            "boolean operation produced no non-degenerate faces".into(),
        )
        .into());
    }

    // Create shell
    let shell_id = store.add_shell(ShellData {
        faces: all_faces,
        is_closed: true,
    });

    // Create solid
    MakeSolid::new(shell_id, vec![]).execute(store)
}

/// Creates one or more faces from a polygon boundary with optional
/// inner boundaries.
///
/// Returns multiple faces (instead of one) when the input boundary is
/// self-intersecting in the face plane: a single boolean fragment can
/// come out with a long boundary segment that overlaps several smaller
/// collinear segments (the multi-door-at-z=0 case), and packing it
/// into a single `MakeFace` call panics `spade::cdt` downstream. To
/// keep those fragments rendering, run the input through the shared
/// 2D planar arrangement engine — for simple inputs it returns the
/// same polygon unchanged; for self-intersecting inputs it decomposes
/// the input into the simple sub-polygons that fill the same area.
fn create_face_from_polygon(
    store: &mut TopologyStore,
    boundary: &[Point3],
    inner_boundaries: &[Vec<Point3>],
    merger: &mut VertexMerger,
) -> Result<Vec<FaceId>> {
    let n = boundary.len();
    if n < 3 {
        return Ok(Vec::new());
    }

    // Reject 3D-collinear fragments before they reach MakeFace, which would
    // otherwise fail with "all points are collinear, cannot define a plane".
    // The split stage filters its own slivers, but vertex merging here can
    // collapse a previously-valid fragment into a degenerate one.
    if newell_normal_3d(boundary).norm() <= TOLERANCE {
        return Ok(Vec::new());
    }

    // Get or create vertices
    let vertex_ids: Vec<VertexId> = boundary
        .iter()
        .map(|p| merger.get_or_create(store, p))
        .collect();

    // Post-snap validation. `VertexMerger` snaps points within its cell
    // tolerance into the same `VertexId`, so a fragment whose edges sit
    // below that tolerance (sub-snap slivers emitted by upstream boolean
    // ops) can collapse into a degenerate boundary. Drop consecutive
    // duplicates so `create_line_edge` never sees a zero-length
    // direction, then re-validate Newell normal / unique vertex count.
    //
    // Note: we do NOT additionally reject fragments whose inner wire
    // shares a vertex with the outer wire (the door-touches-floor
    // pattern). The merge stage emits such "inner touches outer"
    // representations that look invalid in isolation but are how the
    // current `merge_component` encodes opening cuts that connect to
    // the wall boundary. Rejecting them strips legitimate faces and
    // leaves spurious mesh holes. Downstream `BRepSolidToMesh` wraps
    // tessellation in `catch_unwind` so a residual spade panic on
    // those still-invalid PSLG inputs only loses its single face
    // instead of taking the whole mesh build down. Replacing the
    // touching-inner pattern with proper notched outer / split faces
    // is tracked as follow-up work on `merge_component` (planar
    // arrangement refactor).
    let raw_effective: Vec<Point3> = vertex_ids
        .iter()
        .map(|&vid| store.vertex(vid).map(|v| v.point).map_err(Into::into))
        .collect::<Result<Vec<_>>>()?;
    let (_deduped_vertex_ids, effective_boundary) =
        dedupe_consecutive_pairs(&vertex_ids, &raw_effective);
    if effective_boundary.len() < 3 || newell_normal_3d(&effective_boundary).norm() <= TOLERANCE {
        return Ok(Vec::new());
    }

    // Inner wires — apply the same pre/post-snap dedupe.
    let mut snapped_inners: Vec<Vec<Point3>> = Vec::with_capacity(inner_boundaries.len());
    for inner in inner_boundaries {
        if inner.len() < 3 || newell_normal_3d(inner).norm() <= TOLERANCE {
            continue;
        }
        let inner_vids_raw: Vec<VertexId> = inner
            .iter()
            .map(|p| merger.get_or_create(store, p))
            .collect();
        let inner_effective_raw: Vec<Point3> = inner_vids_raw
            .iter()
            .map(|&vid| store.vertex(vid).map(|v| v.point).map_err(Into::into))
            .collect::<Result<Vec<_>>>()?;
        let (_inner_vids, inner_effective) =
            dedupe_consecutive_pairs(&inner_vids_raw, &inner_effective_raw);
        if inner_effective.len() < 3 || newell_normal_3d(&inner_effective).norm() <= TOLERANCE {
            continue;
        }
        snapped_inners.push(inner_effective);
    }

    // Fast path: if every wire is already simple in the face plane,
    // build the face directly. This is the common case and avoids the
    // planar arrangement engine's tendency to insert extra vertices
    // at T-junctions even on simple inputs (which would show up as
    // spurious horizontal seams across the wall mesh).
    //
    // Slow path: if any wire is self-intersecting (proper crossing or
    // collinear overlap of non-adjacent edges), launder the input
    // through the 2D arrangement engine so the result is a set of
    // simple sub-polygons spade can tessellate.
    let normal = newell_normal_3d(&effective_boundary);
    let needs_arrangement =
        !loops_are_simple_on_plane(&normal, &effective_boundary, &snapped_inners);
    if needs_arrangement {
        return build_faces_via_planar_arrangement(
            store,
            merger,
            &effective_boundary,
            &snapped_inners,
        );
    }

    // Direct path — build outer wire + inner wires + face from the
    // snapped (effective) polygon as-is.
    let Some(outer_wire) = build_wire_via_merger(store, merger, &effective_boundary)? else {
        return Ok(Vec::new());
    };
    let mut inner_wires = Vec::with_capacity(snapped_inners.len());
    for inner in &snapped_inners {
        if let Some(w) = build_wire_via_merger(store, merger, inner)? {
            inner_wires.push(w);
        }
    }
    let face_id = MakeFace::new(outer_wire, inner_wires).execute(store)?;
    Ok(vec![face_id])
}

/// Returns `true` when every supplied loop (outer + each inner) is
/// simple in the face plane: no two non-adjacent edges cross at
/// strictly-interior parameters, and no parallel pair overlaps on an
/// interior segment. End-to-end touches are ignored — they're how
/// adjacent boolean fragments meet at shared boundaries, not a defect.
fn loops_are_simple_on_plane(
    normal: &crate::math::Vector3,
    outer: &[Point3],
    inners: &[Vec<Point3>],
) -> bool {
    use crate::math::Vector3;

    let norm = normal.norm();
    if norm <= TOLERANCE {
        return false;
    }
    let n_unit = normal / norm;
    let seed = if n_unit.x.abs() < 0.9 {
        Vector3::new(1.0, 0.0, 0.0)
    } else {
        Vector3::new(0.0, 1.0, 0.0)
    };
    let u_axis = n_unit.cross(&seed).normalize();
    let v_axis = n_unit.cross(&u_axis);
    let origin = outer
        .first()
        .copied()
        .unwrap_or_else(|| Point3::new(0.0, 0.0, 0.0));
    let project = |p: &Point3| -> (f64, f64) {
        let d = *p - origin;
        (d.dot(&u_axis), d.dot(&v_axis))
    };

    let outer_uv: Vec<(f64, f64)> = outer.iter().map(project).collect();
    if !loop_uv_is_simple(&outer_uv) {
        return false;
    }
    for inner in inners {
        let inner_uv: Vec<(f64, f64)> = inner.iter().map(project).collect();
        if !loop_uv_is_simple(&inner_uv) {
            return false;
        }
    }
    true
}

/// Detects proper-crossings and collinear-overlap between non-adjacent
/// edges of a single closed loop projected to UV.
// Segment-intersection math reads clearest in textbook a/b/c/d notation.
#[allow(clippy::many_single_char_names)]
fn loop_uv_is_simple(uv: &[(f64, f64)]) -> bool {
    const CROSS_EPS: f64 = 1e-12;
    const PARAM_EPS: f64 = 1e-9;

    let n = uv.len();
    if n < 4 {
        return true;
    }
    for i in 0..n {
        let a = uv[i];
        let b = uv[(i + 1) % n];
        for j in (i + 2)..n {
            // Skip the wrap-around adjacency (edge n-1 meets edge 0).
            if i == 0 && j == n - 1 {
                continue;
            }
            let c = uv[j];
            let d = uv[(j + 1) % n];
            let d1x = b.0 - a.0;
            let d1y = b.1 - a.1;
            let d2x = d.0 - c.0;
            let d2y = d.1 - c.1;
            let cross = d1x * d2y - d1y * d2x;
            if cross.abs() >= CROSS_EPS {
                // Non-parallel — classic proper crossing test.
                let d3x = c.0 - a.0;
                let d3y = c.1 - a.1;
                let t = (d3x * d2y - d3y * d2x) / cross;
                let u_p = (d3x * d1y - d3y * d1x) / cross;
                if t > PARAM_EPS && t < 1.0 - PARAM_EPS && u_p > PARAM_EPS && u_p < 1.0 - PARAM_EPS
                {
                    return false;
                }
            } else {
                // Parallel — collinear iff c lies on the a→b line.
                let acx = c.0 - a.0;
                let acy = c.1 - a.1;
                let cross_abc = d1x * acy - d1y * acx;
                if cross_abc.abs() >= CROSS_EPS {
                    continue;
                }
                let len2 = d1x * d1x + d1y * d1y;
                if len2 < CROSS_EPS {
                    continue;
                }
                let tc = (acx * d1x + acy * d1y) / len2;
                let adx = d.0 - a.0;
                let ady = d.1 - a.1;
                let td = (adx * d1x + ady * d1y) / len2;
                let (lo, hi) = if tc <= td { (tc, td) } else { (td, tc) };
                if hi > PARAM_EPS && lo < 1.0 - PARAM_EPS {
                    return false;
                }
            }
        }
    }
    true
}

/// Materialises planar arrangement output back into `BRep` faces. Used
/// by `create_face_from_polygon` to launder boolean fragments through
/// the 2D engine so self-intersecting inputs become a set of simple
/// sub-faces.
fn build_faces_via_planar_arrangement(
    store: &mut TopologyStore,
    merger: &mut VertexMerger,
    outer_3d: &[Point3],
    inner_3ds: &[Vec<Point3>],
) -> Result<Vec<FaceId>> {
    use crate::math::Vector3;
    use crate::operations::boolean_2d::{
        run_arrangement, Polygon as B2dPolygon, PolygonWithHoles as B2dPwh, UnionOracle,
    };

    // Plane frame: use the outer Newell normal as the face plane.
    let normal_raw = newell_normal_3d(outer_3d);
    let norm = normal_raw.norm();
    if norm <= TOLERANCE {
        return Ok(Vec::new());
    }
    let n_unit = normal_raw / norm;
    let seed = if n_unit.x.abs() < 0.9 {
        Vector3::new(1.0, 0.0, 0.0)
    } else {
        Vector3::new(0.0, 1.0, 0.0)
    };
    let u_axis = n_unit.cross(&seed).normalize();
    let v_axis = n_unit.cross(&u_axis);
    let origin = outer_3d[0];

    let project = |p: &Point3| -> (f64, f64) {
        let d = *p - origin;
        (d.dot(&u_axis), d.dot(&v_axis))
    };
    let unproject = |uv: (f64, f64)| -> Point3 { origin + u_axis * uv.0 + v_axis * uv.1 };

    let outer_uv: B2dPolygon = outer_3d.iter().map(project).collect();
    let holes_uv: Vec<B2dPolygon> = inner_3ds
        .iter()
        .map(|inner| inner.iter().map(project).collect())
        .collect();
    let pwh = B2dPwh {
        outer: outer_uv,
        holes: holes_uv,
    };
    let inputs = [pwh];
    let oracle = UnionOracle { inputs: &inputs };
    let Ok(arranged) = run_arrangement(&inputs, &oracle) else {
        return Ok(Vec::new());
    };

    let mut face_ids = Vec::with_capacity(arranged.len());
    for pwh in arranged {
        let outer_pts: Vec<Point3> = pwh.outer.iter().copied().map(unproject).collect();
        if outer_pts.len() < 3 {
            continue;
        }
        let outer_wire = build_wire_via_merger(store, merger, &outer_pts)?;
        let Some(outer_wire) = outer_wire else {
            continue;
        };
        let mut inner_wires = Vec::with_capacity(pwh.holes.len());
        for hole in pwh.holes {
            let hole_pts: Vec<Point3> = hole.iter().copied().map(unproject).collect();
            if hole_pts.len() < 3 {
                continue;
            }
            if let Some(wire) = build_wire_via_merger(store, merger, &hole_pts)? {
                inner_wires.push(wire);
            }
        }
        let face_id = MakeFace::new(outer_wire, inner_wires).execute(store)?;
        face_ids.push(face_id);
    }
    Ok(face_ids)
}

/// Builds a closed wire from 3D points, snapping each point through
/// `merger` so neighbouring fragments share `VertexId`s, then dropping
/// consecutive-coincident pairs so `Line::new` never sees a zero-length
/// direction. Returns `None` if the wire collapses to fewer than 3
/// unique vertices.
fn build_wire_via_merger(
    store: &mut TopologyStore,
    merger: &mut VertexMerger,
    points: &[Point3],
) -> Result<Option<crate::topology::WireId>> {
    let vids_raw: Vec<VertexId> = points
        .iter()
        .map(|p| merger.get_or_create(store, p))
        .collect();
    let effective_raw: Vec<Point3> = vids_raw
        .iter()
        .map(|&vid| store.vertex(vid).map(|v| v.point).map_err(Into::into))
        .collect::<Result<Vec<_>>>()?;
    let (vids, effective) = dedupe_consecutive_pairs(&vids_raw, &effective_raw);
    let unique = vids.iter().copied().collect::<HashSet<_>>().len();
    if effective.len() < 3 || unique < 3 {
        return Ok(None);
    }
    let dn = effective.len();
    let mut oriented_edges = Vec::with_capacity(dn);
    for i in 0..dn {
        let j = (i + 1) % dn;
        let edge_id = create_line_edge(store, vids[i], vids[j], effective[i], effective[j])?;
        oriented_edges.push(OrientedEdge::new(edge_id, true));
    }
    let wire_id = store.add_wire(WireData {
        edges: oriented_edges,
        is_closed: true,
    });
    Ok(Some(wire_id))
}

/// Walks the parallel `(VertexId, Point3)` arrays produced by snapping a
/// polygon's corners through `VertexMerger` and drops any entry whose
/// point coincides with its predecessor (within `TOLERANCE`), including
/// the wrap-around pair `[n-1, 0]`. Required before `create_line_edge`
/// because `Line::new` rejects zero-length direction vectors.
fn dedupe_consecutive_pairs(
    vertex_ids: &[VertexId],
    points: &[Point3],
) -> (Vec<VertexId>, Vec<Point3>) {
    debug_assert_eq!(vertex_ids.len(), points.len());
    if points.is_empty() {
        return (Vec::new(), Vec::new());
    }
    let mut out_vids: Vec<VertexId> = Vec::with_capacity(points.len());
    let mut out_pts: Vec<Point3> = Vec::with_capacity(points.len());
    for (&vid, &pt) in vertex_ids.iter().zip(points.iter()) {
        if let Some(prev) = out_pts.last() {
            if (pt - *prev).norm() <= TOLERANCE {
                continue;
            }
        }
        out_vids.push(vid);
        out_pts.push(pt);
    }
    // Drop wrap-around duplicate (last == first). The loop guard
    // ensures `last()` always returns `Some`.
    while out_pts.len() >= 2 {
        let Some(&last) = out_pts.last() else { break };
        let first = out_pts[0];
        if (last - first).norm() <= TOLERANCE {
            out_pts.pop();
            out_vids.pop();
        } else {
            break;
        }
    }
    (out_vids, out_pts)
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
                vec![
                    p(0.0, 0.0, 0.0),
                    p(1.0, 0.0, 0.0),
                    p(1.0, 1.0, 0.0),
                    p(0.0, 1.0, 0.0),
                ],
                false,
            ),
            // Top (z=1)
            make_fragment(
                vec![
                    p(0.0, 0.0, 1.0),
                    p(1.0, 0.0, 1.0),
                    p(1.0, 1.0, 1.0),
                    p(0.0, 1.0, 1.0),
                ],
                false,
            ),
            // Front (y=0)
            make_fragment(
                vec![
                    p(0.0, 0.0, 0.0),
                    p(1.0, 0.0, 0.0),
                    p(1.0, 0.0, 1.0),
                    p(0.0, 0.0, 1.0),
                ],
                false,
            ),
            // Back (y=1)
            make_fragment(
                vec![
                    p(0.0, 1.0, 0.0),
                    p(1.0, 1.0, 0.0),
                    p(1.0, 1.0, 1.0),
                    p(0.0, 1.0, 1.0),
                ],
                false,
            ),
            // Left (x=0)
            make_fragment(
                vec![
                    p(0.0, 0.0, 0.0),
                    p(0.0, 1.0, 0.0),
                    p(0.0, 1.0, 1.0),
                    p(0.0, 0.0, 1.0),
                ],
                false,
            ),
            // Right (x=1)
            make_fragment(
                vec![
                    p(1.0, 0.0, 0.0),
                    p(1.0, 1.0, 0.0),
                    p(1.0, 1.0, 1.0),
                    p(1.0, 0.0, 1.0),
                ],
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
                vec![
                    p(0.0, 0.0, 0.0),
                    p(1.0, 0.0, 0.0),
                    p(1.0, 1.0, 0.0),
                    p(0.0, 1.0, 0.0),
                ],
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
