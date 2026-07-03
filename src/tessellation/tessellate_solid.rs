use crate::error::Result;
use crate::topology::{SolidId, TopologyStore};

use super::edge_samples::EdgeSampleCache;
use super::{TessellateFace, TessellationParams, TriangleMesh};

/// Tessellates all faces of a solid into a combined triangle mesh.
pub struct TessellateSolid {
    solid: SolidId,
    params: TessellationParams,
}

impl TessellateSolid {
    /// Creates a new `TessellateSolid` operation.
    #[must_use]
    pub fn new(solid: SolidId, params: TessellationParams) -> Self {
        Self { solid, params }
    }

    /// Executes the tessellation, returning a combined triangle mesh.
    ///
    /// # Errors
    ///
    /// Returns an error if the solid or any of its faces cannot be tessellated.
    pub fn execute(&self, store: &TopologyStore) -> Result<TriangleMesh> {
        let solid = store.solid(self.solid)?;
        let shell = store.shell(solid.outer_shell)?;

        // One edge-sample cache for the whole shell: faces sharing an edge
        // consume the identical boundary polyline (structural conformance).
        let mut cache = EdgeSampleCache::new(self.params);

        let mut combined = TriangleMesh::default();
        for &face_id in &shell.faces {
            let face_mesh =
                TessellateFace::new(face_id, self.params).execute_with_cache(store, &mut cache)?;
            combined.merge(&face_mesh);
        }

        Ok(combined)
    }
}

/// Squared distance from point `p` to segment `[a, b]`.
#[cfg(test)]
fn point_segment_dist_sq(
    p: crate::math::Point3,
    a: crate::math::Point3,
    b: crate::math::Point3,
) -> f64 {
    let ab = b - a;
    let len_sq = ab.norm_squared();
    if len_sq < 1e-30 {
        return (p - a).norm_squared();
    }
    let t = ((p - a).dot(&ab) / len_sq).clamp(0.0, 1.0);
    (p - (a + ab * t)).norm_squared()
}

/// A boundary edge of a face mesh: mesh vertex indices plus 3D endpoints.
#[cfg(test)]
#[derive(Clone, Copy)]
struct BoundarySeg {
    a_idx: u32,
    b_idx: u32,
    a: crate::math::Point3,
    b: crate::math::Point3,
}

/// Boundary edges of a single face mesh: undirected triangle edges used by
/// exactly one triangle.
#[cfg(test)]
fn face_boundary_edges(mesh: &TriangleMesh) -> Vec<BoundarySeg> {
    use std::collections::HashMap;
    let mut counts: HashMap<(u32, u32), usize> = HashMap::new();
    for tri in &mesh.indices {
        for &(a, b) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
            let key = if a < b { (a, b) } else { (b, a) };
            *counts.entry(key).or_insert(0) += 1;
        }
    }
    counts
        .iter()
        .filter(|(_, &c)| c == 1)
        .map(|(&(a, b), _)| BoundarySeg {
            a_idx: a,
            b_idx: b,
            a: mesh.vertices[a as usize],
            b: mesh.vertices[b as usize],
        })
        .collect()
}

/// Maximum 3D deviation between boundary polylines: for every boundary vertex
/// of every face, the distance to the nearest boundary segment on any OTHER
/// face — or on the SAME face when the segment is not incident to the vertex
/// (by mesh index). On a conforming closed solid this is ~0; where adjacent
/// faces tessellate a shared boundary curve with disagreeing chords it equals
/// the chord sagitta — the visible silhouette sliver.
///
/// Same-face non-incident segments are legitimate conforming targets because
/// of parametric seams: a face on a geometrically closed surface has two
/// boundary columns (`u0` / `u1`) that coincide in 3D — an internal seam, not
/// a crack. Only index-incident segments are excluded so a vertex never
/// trivially conforms to its own edges.
///
/// Used as the regression metric for the boundary-conforming tessellation,
/// here and in the NURBS boolean tests.
#[cfg(test)]
pub(crate) fn max_adjacent_boundary_deviation(store: &TopologyStore, solid: SolidId) -> f64 {
    #[allow(clippy::unwrap_used)]
    let shell = store
        .shell(store.solid(solid).unwrap().outer_shell)
        .unwrap();
    let per_face: Vec<Vec<BoundarySeg>> = shell
        .faces
        .iter()
        .map(|&f| {
            #[allow(clippy::unwrap_used)]
            let mesh = TessellateFace::new(f, TessellationParams::default())
                .execute(store)
                .unwrap();
            face_boundary_edges(&mesh)
        })
        .collect();

    let mut max_dev = 0.0_f64;
    for (i, edges_i) in per_face.iter().enumerate() {
        for seg in edges_i {
            for (v_idx, v) in [(seg.a_idx, seg.a), (seg.b_idx, seg.b)] {
                let mut best = f64::INFINITY;
                for (j, edges_j) in per_face.iter().enumerate() {
                    for other in edges_j {
                        // Same face: skip segments incident to the vertex.
                        if i == j && (other.a_idx == v_idx || other.b_idx == v_idx) {
                            continue;
                        }
                        best = best.min(point_segment_dist_sq(v, other.a, other.b));
                    }
                }
                if best.is_finite() {
                    max_dev = max_dev.max(best.sqrt());
                }
            }
        }
    }
    max_dev
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::operations::creation::MakeCurvedSlab;

    /// The plain curved slab's adjacent faces (curved top/bottom vs ruled side
    /// walls) now tessellate their shared boundary curves at identical
    /// parameters, so the silhouette slivers are gone: the max adjacent-boundary
    /// deviation drops from the chord sagitta (~2e-2) to floating-point noise.
    #[test]
    fn plain_slab_boundaries_conform() {
        let mut store = TopologyStore::new();
        let solid = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut store)
            .unwrap();
        let dev = max_adjacent_boundary_deviation(&store, solid);
        assert!(
            dev < 1e-6,
            "plain slab adjacent-boundary deviation {dev} exceeds 1e-6"
        );
    }

    /// F2 target: the tube's planar caps and NURBS side wall share true ring
    /// edges sampled once per edge, so the cap rim and the side boundary emit
    /// identical vertices.
    #[test]
    fn tube_boundaries_conform() {
        use crate::operations::creation::MakeNurbsTube;
        let mut store = TopologyStore::new();
        let solid = MakeNurbsTube::new(crate::math::Point3::new(0.0, 0.0, 0.0), 0.8, 3.0)
            .execute(&mut store)
            .unwrap();
        let dev = max_adjacent_boundary_deviation(&store, solid);
        assert!(dev < 1e-6, "tube cap/side deviation {dev} exceeds 1e-6");
    }

    /// F2 target: the revolved solid's disk caps reference the wall's true
    /// boundary circles instead of independent 48-gons.
    #[test]
    fn revolved_solid_boundaries_conform() {
        use crate::operations::creation::MakeRevolvedSolid;
        let mut store = TopologyStore::new();
        let solid = MakeRevolvedSolid::new(vec![(2.0, 0.0), (2.4, 1.2), (2.1, 2.4)])
            .execute(&mut store)
            .unwrap();
        let dev = max_adjacent_boundary_deviation(&store, solid);
        assert!(dev < 1e-6, "revolved cap/wall deviation {dev} exceeds 1e-6");
    }
}
