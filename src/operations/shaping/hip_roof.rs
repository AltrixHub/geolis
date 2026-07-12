//! Hip roof solid from a baseline polygon via the straight skeleton.
//!
//! The baseline is the eave line footprint drawn by the user. Every roof
//! face rises from the eave at the uniform slope implied by `rise` (the
//! height difference between the baseline and the ridge):
//! `slope = rise / max_inset(baseline skeleton)`, so `ridge_z =
//! baseline_z + rise` exactly. An optional eave overhang extends every
//! face outward past the baseline (mitred corners), lowering the eave
//! edge to `baseline_z - slope * overhang`.
//!
//! The solid is bounded by one planar sloped face per baseline edge (its
//! straight-skeleton cell plus the overhang band) and a flat bottom cap
//! at the eave polygon, so it is watertight without vertical sides.

use std::collections::HashMap;

use super::extrude::{create_closed_wire, create_line_edge};
use crate::error::{OperationError, Result};
use crate::math::straight_skeleton::{
    compute_straight_skeleton, ring_self_intersection, SkeletonCell,
};
use crate::math::{Point2, Point3, Vector2, TOLERANCE};
use crate::operations::creation::{MakeFace, MakeSolid};
use crate::topology::{
    EdgeId, OrientedEdge, ShellData, SolidId, TopologyStore, VertexData, VertexId,
};

/// Longest allowed mitre extension relative to the overhang: corners whose
/// mitre point would land farther than `MITER_CAP * overhang` from the
/// baseline vertex are rejected (interior angle sharper than ~6 degrees).
const MITER_CAP: f64 = 20.0;

/// Builds a hip roof solid from a baseline polygon.
///
/// The baseline is interpreted as a closed ring in the XY plane (z is
/// ignored); non-convex simple polygons are supported.
pub struct MakeHipRoof {
    baseline: Vec<Point3>,
    rise: f64,
    overhang: f64,
    baseline_z: f64,
}

impl MakeHipRoof {
    /// Creates a hip roof operation with no overhang, based at z = 0.
    #[must_use]
    pub fn new(baseline: Vec<Point3>, rise: f64) -> Self {
        Self {
            baseline,
            rise,
            overhang: 0.0,
            baseline_z: 0.0,
        }
    }

    /// Sets the eave overhang: every roof face is extended outward past
    /// the baseline by this distance (measured in plan).
    #[must_use]
    pub fn with_overhang(mut self, overhang: f64) -> Self {
        self.overhang = overhang;
        self
    }

    /// Sets the z at which the baseline (eave line) sits.
    #[must_use]
    pub fn with_baseline_z(mut self, z: f64) -> Self {
        self.baseline_z = z;
        self
    }

    /// Executes the operation, creating the roof solid in the store.
    ///
    /// # Errors
    ///
    /// Returns [`OperationError::InvalidInput`] when the baseline is not a
    /// simple polygon with at least 3 distinct vertices, `rise` is not
    /// strictly positive, `overhang` is negative, any parameter is
    /// non-finite, or the overhang is too large for the footprint (the
    /// offset eave polygon self-intersects or a mitre explodes). Returns
    /// [`OperationError::Failed`] on numerically degenerate input.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        let (skeleton, slope, eave) = self.prepare()?;
        let lift = |p: Point2, inset: f64| -> Point3 {
            Point3::new(p.x, p.y, self.baseline_z + slope * inset)
        };

        let mut mesh = SharedTopology::default();
        let mut faces = Vec::with_capacity(skeleton.cells.len() + 1);
        for cell in &skeleton.cells {
            let polygon = face_polygon(cell, &eave, self.overhang, &lift);
            faces.push(mesh.add_face(store, &polygon)?);
        }

        // Bottom cap at the eave plane, wound to face downward.
        let cap: Vec<Point3> = eave
            .iter()
            .rev()
            .map(|&p| lift(p, -self.overhang))
            .collect();
        faces.push(mesh.add_face(store, &cap)?);

        let shell = store.add_shell(ShellData {
            faces,
            is_closed: true,
        });
        MakeSolid::new(shell, vec![]).execute(store)
    }

    /// Returns the eave polygon of the roof: the baseline offset outward
    /// by the overhang (mitred corners), lifted to the eave height
    /// `baseline_z - slope * overhang`. This is the plan-view outline of
    /// the roof.
    ///
    /// # Errors
    ///
    /// Same conditions as [`MakeHipRoof::execute`].
    pub fn eave_ring(&self) -> Result<Vec<Point3>> {
        let (_, slope, eave) = self.prepare()?;
        let eave_z = self.baseline_z - slope * self.overhang;
        Ok(eave.iter().map(|p| Point3::new(p.x, p.y, eave_z)).collect())
    }

    /// Validates parameters and computes the skeleton, slope, and eave
    /// corners shared by [`MakeHipRoof::execute`] and
    /// [`MakeHipRoof::eave_ring`].
    fn prepare(
        &self,
    ) -> Result<(
        crate::math::straight_skeleton::StraightSkeleton,
        f64,
        Vec<Point2>,
    )> {
        if !self.rise.is_finite() || self.rise <= 0.0 {
            return Err(OperationError::InvalidInput(format!(
                "hip roof rise must be strictly positive, got {}",
                self.rise
            ))
            .into());
        }
        if !self.overhang.is_finite() || self.overhang < 0.0 {
            return Err(OperationError::InvalidInput(format!(
                "hip roof overhang must be non-negative, got {}",
                self.overhang
            ))
            .into());
        }
        if !self.baseline_z.is_finite() {
            return Err(
                OperationError::InvalidInput("hip roof baseline z must be finite".into()).into(),
            );
        }
        let skeleton = compute_straight_skeleton(&self.baseline)?;
        let ring: Vec<Point2> = skeleton
            .polygon
            .iter()
            .map(|p| Point2::new(p.x, p.y))
            .collect();
        if skeleton.max_inset < TOLERANCE {
            return Err(OperationError::Failed(
                "hip roof: degenerate skeleton (zero inset)".into(),
            )
            .into());
        }
        let slope = self.rise / skeleton.max_inset;
        let eave = self.eave_corners(&ring)?;
        Ok((skeleton, slope, eave))
    }

    /// Computes the mitred outward offset of the baseline ring by the
    /// overhang distance. With zero overhang this is the ring itself.
    fn eave_corners(&self, ring: &[Point2]) -> Result<Vec<Point2>> {
        if self.overhang == 0.0 {
            return Ok(ring.to_vec());
        }
        let n = ring.len();
        let mut corners = Vec::with_capacity(n);
        for i in 0..n {
            let prev = ring[(i + n - 1) % n];
            let curr = ring[i];
            let next = ring[(i + 1) % n];
            let dir_in = (curr - prev).normalize();
            let dir_out = (next - curr).normalize();
            let normal_in = Vector2::new(-dir_in.y, dir_in.x);
            let normal_out = Vector2::new(-dir_out.y, dir_out.x);
            let denom = 1.0 + normal_in.dot(&normal_out);
            if denom.abs() < 1e-9 {
                return Err(OperationError::InvalidInput(
                    "hip roof: overhang cannot be mitred at a spike corner".into(),
                )
                .into());
            }
            let velocity = (normal_in + normal_out) / denom;
            if velocity.norm() > MITER_CAP {
                return Err(OperationError::InvalidInput(
                    "hip roof: overhang too large for this footprint (sharp corner)".into(),
                )
                .into());
            }
            corners.push(curr - self.overhang * velocity);
        }
        for i in 0..n {
            let j = (i + 1) % n;
            // Each eave edge must keep the direction of its baseline edge:
            // a reversed or collapsed edge means opposing wavefronts have
            // crossed (e.g. the overhang bridges a cavity of the footprint).
            let edge_dir = (ring[j] - ring[i]).normalize();
            if (corners[j] - corners[i]).dot(&edge_dir) < TOLERANCE {
                return Err(OperationError::InvalidInput(
                    "hip roof: overhang too large for this footprint (eave edge collapses)".into(),
                )
                .into());
            }
        }
        if ring_self_intersection(&corners).is_some() {
            return Err(OperationError::InvalidInput(
                "hip roof: overhang too large for this footprint (eave self-intersects)".into(),
            )
            .into());
        }
        Ok(corners)
    }
}

/// Builds the 3D boundary polygon (CCW seen from above) of the roof face
/// belonging to one baseline edge: the mitred overhang band followed by
/// the lifted skeleton cell chain.
fn face_polygon(
    cell: &SkeletonCell,
    eave: &[Point2],
    overhang: f64,
    lift: &impl Fn(Point2, f64) -> Point3,
) -> Vec<Point3> {
    let edge = cell.edge_index;
    let mut polygon = Vec::with_capacity(cell.vertices.len() + 2);
    if overhang > 0.0 {
        let next = (edge + 1) % eave.len();
        polygon.push(lift(eave[edge], -overhang));
        polygon.push(lift(eave[next], -overhang));
        // Cell vertices are [start, end, chain from end back to start]:
        // with the eave band in front, traversal continues at the edge end
        // vertex and walks the chain back to the edge start vertex.
        for v in cell.vertices.iter().skip(1) {
            polygon.push(lift(Point2::new(v.position.x, v.position.y), v.inset));
        }
        polygon.push(lift(
            Point2::new(cell.vertices[0].position.x, cell.vertices[0].position.y),
            cell.vertices[0].inset,
        ));
    } else {
        for v in &cell.vertices {
            polygon.push(lift(Point2::new(v.position.x, v.position.y), v.inset));
        }
    }
    polygon
}

/// Vertex- and edge-deduplicating face builder: faces created through it
/// share vertices and edges, so the resulting shell has every edge used
/// exactly twice. Positions are keyed exactly (by bit pattern), which is
/// sound here because coincident face corners are copies of the same
/// skeleton node / eave corner values.
#[derive(Default)]
struct SharedTopology {
    vertices: HashMap<[u64; 3], VertexId>,
    edges: HashMap<(VertexId, VertexId), EdgeId>,
}

impl SharedTopology {
    fn vertex(&mut self, store: &mut TopologyStore, p: Point3) -> VertexId {
        let key = [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
        if let Some(&id) = self.vertices.get(&key) {
            return id;
        }
        let id = store.add_vertex(VertexData::new(p));
        self.vertices.insert(key, id);
        id
    }

    fn add_face(
        &mut self,
        store: &mut TopologyStore,
        polygon: &[Point3],
    ) -> Result<crate::topology::FaceId> {
        if polygon.len() < 3 {
            return Err(OperationError::Failed("hip roof: degenerate face polygon".into()).into());
        }
        let ids: Vec<VertexId> = polygon.iter().map(|&p| self.vertex(store, p)).collect();
        let mut oriented = Vec::with_capacity(ids.len());
        for i in 0..ids.len() {
            let j = (i + 1) % ids.len();
            let (a, b) = (ids[i], ids[j]);
            if a == b {
                return Err(OperationError::Failed(
                    "hip roof: zero-length edge in face polygon".into(),
                )
                .into());
            }
            let key = (a.min(b), a.max(b));
            let forward = a < b;
            let edge_id = if let Some(&id) = self.edges.get(&key) {
                id
            } else {
                let (start, end) = if forward { (i, j) } else { (j, i) };
                let id =
                    create_line_edge(store, ids[start], ids[end], polygon[start], polygon[end])?;
                self.edges.insert(key, id);
                id
            };
            oriented.push(OrientedEdge::new(edge_id, forward));
        }
        let wire = create_closed_wire(store, oriented);
        MakeFace::new(wire, vec![]).execute(store)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use std::collections::HashMap;

    use super::*;
    use crate::operations::query::Volume;
    use crate::tessellation::{TessellateSolid, TessellationParams};

    fn ring(pts: &[(f64, f64)]) -> Vec<Point3> {
        pts.iter().map(|&(x, y)| Point3::new(x, y, 0.0)).collect()
    }

    fn square() -> Vec<Point3> {
        ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)])
    }

    fn rect() -> Vec<Point3> {
        ring(&[(0.0, 0.0), (6.0, 0.0), (6.0, 4.0), (0.0, 4.0)])
    }

    fn l_shape() -> Vec<Point3> {
        ring(&[
            (0.0, 0.0),
            (6.0, 0.0),
            (6.0, 3.0),
            (3.0, 3.0),
            (3.0, 6.0),
            (0.0, 6.0),
        ])
    }

    fn edge_usage_is_two_everywhere(store: &TopologyStore, solid: SolidId) -> bool {
        let shell_id = store.solid(solid).unwrap().outer_shell;
        let shell = store.shell(shell_id).unwrap();
        let mut usage: HashMap<crate::topology::EdgeId, usize> = HashMap::new();
        for &face_id in &shell.faces {
            let face = store.face(face_id).unwrap();
            let wire = store.wire(face.outer_wire).unwrap();
            for oe in &wire.edges {
                *usage.entry(oe.edge).or_insert(0) += 1;
            }
        }
        usage.values().all(|&count| count == 2)
    }

    /// Quantizes tessellated vertex positions and asserts every undirected
    /// triangle edge is shared by exactly two triangles (no boundary).
    #[allow(clippy::cast_possible_truncation)]
    fn assert_position_weld_watertight(mesh: &crate::tessellation::TriangleMesh) {
        let quantize = |value: f64| (value * 1e6).round() as i64;
        let key = |p: &Point3| -> (i64, i64, i64) { (quantize(p.x), quantize(p.y), quantize(p.z)) };
        let mut ids: HashMap<(i64, i64, i64), usize> = HashMap::new();
        let mut vertex_ids = Vec::with_capacity(mesh.vertices.len());
        for p in &mesh.vertices {
            let next = ids.len();
            let id = *ids.entry(key(p)).or_insert(next);
            vertex_ids.push(id);
        }
        let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
        for tri in &mesh.indices {
            for k in 0..3 {
                let a = vertex_ids[tri[k] as usize];
                let b = vertex_ids[tri[(k + 1) % 3] as usize];
                if a == b {
                    continue;
                }
                *edge_count.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
        }
        let boundary = edge_count.values().filter(|&&count| count != 2).count();
        assert_eq!(
            boundary, 0,
            "mesh has {boundary} non-manifold/boundary edges"
        );
    }

    fn face_count(store: &TopologyStore, solid: SolidId) -> usize {
        let shell_id = store.solid(solid).unwrap().outer_shell;
        store.shell(shell_id).unwrap().faces.len()
    }

    fn mesh_z_range(mesh: &crate::tessellation::TriangleMesh) -> (f64, f64) {
        let mut min = f64::INFINITY;
        let mut max = f64::NEG_INFINITY;
        for p in &mesh.vertices {
            min = min.min(p.z);
            max = max.max(p.z);
        }
        (min, max)
    }

    #[test]
    fn square_pyramid() {
        let mut store = TopologyStore::new();
        let solid = MakeHipRoof::new(square(), 2.0).execute(&mut store).unwrap();
        // 4 sloped faces + bottom cap.
        assert_eq!(face_count(&store, solid), 5);
        assert!(edge_usage_is_two_everywhere(&store, solid));
        let volume = Volume::new(solid).execute(&store).unwrap();
        assert!((volume - 32.0 / 3.0).abs() < 1e-9, "volume {volume}");
        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert_eq!(mesh.vertices.len(), mesh.normals.len());
        let (min_z, max_z) = mesh_z_range(&mesh);
        assert!(min_z.abs() < 1e-9);
        assert!((max_z - 2.0).abs() < 1e-9, "apex at {max_z}");
        assert_position_weld_watertight(&mesh);
    }

    #[test]
    fn rectangle_ridge_volume() {
        let mut store = TopologyStore::new();
        let solid = MakeHipRoof::new(rect(), 1.0).execute(&mut store).unwrap();
        assert_eq!(face_count(&store, solid), 5);
        assert!(edge_usage_is_two_everywhere(&store, solid));
        // V = integral of (6 - 4z)(4 - 4z) dz over [0, 1] = 28/3.
        let volume = Volume::new(solid).execute(&store).unwrap();
        assert!((volume - 28.0 / 3.0).abs() < 1e-9, "volume {volume}");
        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();
        let (min_z, max_z) = mesh_z_range(&mesh);
        assert!(min_z.abs() < 1e-9);
        assert!((max_z - 1.0).abs() < 1e-9, "ridge at {max_z}");
    }

    #[test]
    fn overhang_lowers_eave_and_matches_rebased_equivalent() {
        // rect 6x4, rise 1 -> max_inset 2, slope 0.5. With overhang 0.5 the
        // eave polygon is the 7x5 rectangle at z = -0.25.
        let mut store = TopologyStore::new();
        let solid = MakeHipRoof::new(rect(), 1.0)
            .with_overhang(0.5)
            .execute(&mut store)
            .unwrap();
        assert!(edge_usage_is_two_everywhere(&store, solid));
        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();
        let (min_z, max_z) = mesh_z_range(&mesh);
        assert!((min_z + 0.25).abs() < 1e-9, "eave at {min_z}");
        assert!((max_z - 1.0).abs() < 1e-9, "ridge at {max_z}");
        assert_position_weld_watertight(&mesh);
        let volume = Volume::new(solid).execute(&store).unwrap();

        // The same solid expressed with the eave rectangle as baseline:
        // 7x5 has max_inset 2.5, rise 1.25 gives the same slope 0.5.
        let mut store2 = TopologyStore::new();
        let eave_rect = ring(&[(-0.5, -0.5), (6.5, -0.5), (6.5, 4.5), (-0.5, 4.5)]);
        let equivalent = MakeHipRoof::new(eave_rect, 1.25)
            .with_baseline_z(-0.25)
            .execute(&mut store2)
            .unwrap();
        let volume2 = Volume::new(equivalent).execute(&store2).unwrap();
        assert!(
            (volume - volume2).abs() < 1e-9,
            "overhang form {volume} vs rebased form {volume2}"
        );
    }

    #[test]
    fn l_shape_with_overhang_is_watertight() {
        let mut store = TopologyStore::new();
        let solid = MakeHipRoof::new(l_shape(), 1.5)
            .with_overhang(0.3)
            .execute(&mut store)
            .unwrap();
        // 6 sloped faces + bottom cap.
        assert_eq!(face_count(&store, solid), 7);
        assert!(edge_usage_is_two_everywhere(&store, solid));
        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert_position_weld_watertight(&mesh);
        // slope = 1.5 / 1.5 = 1.0; eave sits at -0.3, ridge at 1.5.
        let (min_z, max_z) = mesh_z_range(&mesh);
        assert!((min_z + 0.3).abs() < 1e-9, "eave at {min_z}");
        assert!((max_z - 1.5).abs() < 1e-9, "ridge at {max_z}");
        for p in &mesh.vertices {
            assert!(p.x.is_finite() && p.y.is_finite() && p.z.is_finite());
        }
    }

    #[test]
    fn baseline_z_lifts_the_solid() {
        let mut store = TopologyStore::new();
        let solid = MakeHipRoof::new(square(), 2.0)
            .with_overhang(0.5)
            .with_baseline_z(5.0)
            .execute(&mut store)
            .unwrap();
        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();
        // slope = 2 / 2 = 1: eave at 5 - 0.5, ridge at 5 + 2.
        let (min_z, max_z) = mesh_z_range(&mesh);
        assert!((min_z - 4.5).abs() < 1e-9, "eave at {min_z}");
        assert!((max_z - 7.0).abs() < 1e-9, "ridge at {max_z}");
    }

    #[test]
    fn eave_ring_is_the_mitred_offset_at_eave_height() {
        // square 4x4, rise 2 -> slope 1; overhang 0.5 -> eave z = -0.5.
        let op = MakeHipRoof::new(square(), 2.0).with_overhang(0.5);
        let eave = op.eave_ring().unwrap();
        assert_eq!(eave.len(), 4);
        for p in &eave {
            assert!((p.z + 0.5).abs() < 1e-9, "eave z {}", p.z);
            assert!((p.x + 0.5).abs() < 1e-9 || (p.x - 4.5).abs() < 1e-9);
            assert!((p.y + 0.5).abs() < 1e-9 || (p.y - 4.5).abs() < 1e-9);
        }
        // Zero overhang: the eave ring is the baseline at baseline_z.
        let flat = MakeHipRoof::new(square(), 2.0)
            .with_baseline_z(3.0)
            .eave_ring()
            .unwrap();
        assert_eq!(flat.len(), 4);
        for p in &flat {
            assert!((p.z - 3.0).abs() < 1e-9);
        }
    }

    #[test]
    fn rejects_invalid_parameters() {
        let mut store = TopologyStore::new();
        assert!(MakeHipRoof::new(square(), 0.0).execute(&mut store).is_err());
        assert!(MakeHipRoof::new(square(), -1.0)
            .execute(&mut store)
            .is_err());
        assert!(MakeHipRoof::new(square(), f64::NAN)
            .execute(&mut store)
            .is_err());
        assert!(MakeHipRoof::new(square(), 2.0)
            .with_overhang(-0.1)
            .execute(&mut store)
            .is_err());
        assert!(MakeHipRoof::new(square(), 2.0)
            .with_baseline_z(f64::INFINITY)
            .execute(&mut store)
            .is_err());
    }

    #[test]
    fn rejects_degenerate_baselines() {
        let mut store = TopologyStore::new();
        // Bowtie.
        assert!(
            MakeHipRoof::new(ring(&[(0.0, 0.0), (4.0, 4.0), (4.0, 0.0), (0.0, 4.0)]), 1.0)
                .execute(&mut store)
                .is_err()
        );
        // Too few vertices.
        assert!(MakeHipRoof::new(ring(&[(0.0, 0.0), (4.0, 0.0)]), 1.0)
            .execute(&mut store)
            .is_err());
    }

    #[test]
    fn rejects_overhang_bridging_a_cavity() {
        // U-shape with a 3-wide cavity: an overhang of 5 makes the two
        // cavity-wall wavefronts cross, reversing the cavity-bottom eave
        // edge.
        let mut store = TopologyStore::new();
        let u_shape = ring(&[
            (0.0, 0.0),
            (9.0, 0.0),
            (9.0, 6.0),
            (6.0, 6.0),
            (6.0, 2.0),
            (3.0, 2.0),
            (3.0, 6.0),
            (0.0, 6.0),
        ]);
        assert!(MakeHipRoof::new(u_shape, 1.5)
            .with_overhang(5.0)
            .execute(&mut store)
            .is_err());
    }

    #[test]
    fn large_overhang_on_l_shape_stays_consistent() {
        // The L-shape's mitred outward offset remains a simple ring even
        // for a large overhang; the roof must stay watertight.
        let mut store = TopologyStore::new();
        let solid = MakeHipRoof::new(l_shape(), 1.5)
            .with_overhang(5.0)
            .execute(&mut store)
            .unwrap();
        assert!(edge_usage_is_two_everywhere(&store, solid));
        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert_position_weld_watertight(&mesh);
    }
}
