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
//! By default the solid is a wedge, bounded by one planar sloped face per
//! baseline edge (its straight-skeleton cell plus the overhang band) and a
//! flat bottom cap at the eave polygon, so it is watertight without
//! vertical sides.
//!
//! # Thickness
//!
//! [`MakeHipRoof::with_thickness`] turns the wedge into a roof *shell*:
//! the sloped faces become the top sheet, a second copy of the same
//! piecewise-planar height field lowered by `thickness` becomes the
//! bottom sheet, and a vertical fascia band of height `thickness` closes
//! the two along the outer eave ring. The flat bottom cap only exists in
//! the wedge form (`thickness == 0`, the default).
//!
//! The offset is **vertical, not normal**: the bottom sheet is the top
//! sheet translated by `-thickness` in z, so both sheets share the exact
//! same plan-view decomposition into skeleton cells and overhang bands.
//! That makes the shell watertight by construction — every interior
//! ridge/hip edge is shared by the two adjacent cells on its own sheet,
//! and every eave-ring edge is shared with the fascia quad above or below
//! it — and it makes the enclosed volume exactly
//! `plan_area(eave_polygon) * thickness`.
//!
//! The consequence of the vertical (rather than normal) offset is that the
//! shell's thickness measured perpendicular to the roof plane is
//! `thickness * cos(pitch_angle)`, i.e. `thickness / sqrt(1 + slope^2)`,
//! not `thickness`. That is the standard simplified-BIM roof
//! representation: the parameter is the vertical build-up a section cut
//! shows, and it keeps the eave fascia a clean vertical band instead of
//! the mitred, self-intersection-prone geometry a true normal offset
//! would need at hips and valleys.
//!
//! # Slope clearance
//!
//! Inside the overhang band the sloped surface passes *below* the
//! baseline: `surface_z = baseline_z - slope * distance_outside_baseline`.
//! A caller that drew the baseline on a wall centreline therefore needs
//! the whole solid lifted so the surface clears the wall's outer face.
//! [`MakeHipRoof::with_slope_clearance`] expresses that as a plan
//! distance: the solid is raised by `slope * clearance`, so the surface at
//! `clearance` outside the baseline sits back at `baseline_z`.
//!
//! The lift belongs to this operation because `slope` is derived from the
//! skeleton's `max_inset` — it is not knowable from the constructor
//! arguments, so a caller could only pre-compute the lift by running the
//! skeleton itself.

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
    slope_clearance: f64,
    thickness: f64,
}

/// The skeleton-derived quantities every public entry point shares.
struct RoofPlan {
    skeleton: crate::math::straight_skeleton::StraightSkeleton,
    /// Rise per unit plan run — the slope of every sloped face.
    slope: f64,
    /// Baseline ring offset outward by the overhang (mitred corners).
    eave: Vec<Point2>,
    /// The baseline z after the slope-clearance lift.
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
            slope_clearance: 0.0,
            thickness: 0.0,
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

    /// Lifts the whole solid so the sloped surface at `clearance` (a plan
    /// distance, measured outward from the baseline) sits at the baseline
    /// z instead of `slope * clearance` below it.
    ///
    /// The lift is `slope * clearance`: a pure z translation, so it changes
    /// neither the plan-view outline nor the volume nor any surface area.
    /// A clearance of `0` reproduces the unlifted solid exactly. See the
    /// module docs for why the lift belongs to this operation.
    #[must_use]
    pub fn with_slope_clearance(mut self, clearance: f64) -> Self {
        self.slope_clearance = clearance;
        self
    }

    /// Builds the roof as a shell of this vertical thickness instead of a
    /// solid wedge: the sloped faces become the top sheet, the same height
    /// field lowered by `thickness` becomes the bottom sheet, and a
    /// vertical fascia band closes the two along the eave ring.
    ///
    /// The offset is vertical, so the plan-view outline is unchanged and
    /// the enclosed volume is exactly `plan_area(eave_polygon) *
    /// thickness`; the thickness measured perpendicular to the roof plane
    /// is `thickness / sqrt(1 + slope^2)`. See the module docs.
    ///
    /// A thickness of `0` (the default) reproduces the solid wedge —
    /// sloped faces plus a flat bottom cap — exactly. The whole shell
    /// participates in [`MakeHipRoof::with_slope_clearance`]'s lift, and
    /// [`MakeHipRoof::slope`] / [`MakeHipRoof::eave_ring`] keep describing
    /// the top surface.
    #[must_use]
    pub fn with_thickness(mut self, thickness: f64) -> Self {
        self.thickness = thickness;
        self
    }

    /// Returns the roof pitch: the rise per unit plan run, i.e. the slope
    /// of every sloped face (`rise / max_inset` of the baseline skeleton).
    ///
    /// # Errors
    ///
    /// Same conditions as [`MakeHipRoof::execute`].
    pub fn slope(&self) -> Result<f64> {
        Ok(self.plan()?.slope)
    }

    /// Executes the operation, creating the roof solid in the store.
    ///
    /// # Errors
    ///
    /// Returns [`OperationError::InvalidInput`] when the baseline is not a
    /// simple polygon with at least 3 distinct vertices, `rise` is not
    /// strictly positive, `overhang`, the slope clearance or the thickness
    /// is negative, any parameter is non-finite, or the overhang is too
    /// large for the footprint (the offset eave polygon self-intersects or
    /// a mitre explodes). Returns [`OperationError::Failed`] on
    /// numerically degenerate input.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        let RoofPlan {
            skeleton,
            slope,
            eave,
            baseline_z,
        } = self.plan()?;
        let top =
            |p: Point2, inset: f64| -> Point3 { Point3::new(p.x, p.y, baseline_z + slope * inset) };
        // The bottom sheet subtracts from the *evaluated* top height rather
        // than re-deriving the height field from a lowered base: that keeps
        // the drop exactly `thickness` at every vertex in floating point,
        // which is what makes the enclosed volume exactly plan area times
        // thickness.
        let thickness = self.thickness;
        let bottom = |p: Point2, inset: f64| -> Point3 {
            let t = top(p, inset);
            Point3::new(t.x, t.y, t.z - thickness)
        };

        let mut mesh = SharedTopology::default();
        let mut faces = Vec::with_capacity(2 * skeleton.cells.len() + eave.len());
        // Top sheet: one sloped face per baseline edge.
        for cell in &skeleton.cells {
            let polygon = face_polygon(cell, &eave, self.overhang, &top);
            faces.push(mesh.add_face(store, &polygon)?);
        }

        if thickness > 0.0 {
            // Bottom sheet: the same height field lowered by `thickness`,
            // each face reversed so it faces downward.
            for cell in &skeleton.cells {
                let mut polygon = face_polygon(cell, &eave, self.overhang, &bottom);
                polygon.reverse();
                faces.push(mesh.add_face(store, &polygon)?);
            }
            // Eave fascia: one outward-facing vertical quad per eave-ring
            // edge, closing the two sheets. The eave ring is CCW in plan,
            // so `start_top -> start_bottom -> end_bottom -> end_top`
            // winds outward.
            for i in 0..eave.len() {
                let j = (i + 1) % eave.len();
                let quad = [
                    top(eave[i], -self.overhang),
                    bottom(eave[i], -self.overhang),
                    bottom(eave[j], -self.overhang),
                    top(eave[j], -self.overhang),
                ];
                faces.push(mesh.add_face(store, &quad)?);
            }
        } else {
            // Wedge form: flat bottom cap at the eave plane, wound to face
            // downward.
            let cap: Vec<Point3> = eave.iter().rev().map(|&p| top(p, -self.overhang)).collect();
            faces.push(mesh.add_face(store, &cap)?);
        }

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
        let plan = self.plan()?;
        let eave_z = plan.baseline_z - plan.slope * self.overhang;
        Ok(plan
            .eave
            .iter()
            .map(|p| Point3::new(p.x, p.y, eave_z))
            .collect())
    }

    /// Validates parameters and computes the skeleton, slope, lifted
    /// baseline z, and eave corners shared by every public entry point.
    fn plan(&self) -> Result<RoofPlan> {
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
        if !self.slope_clearance.is_finite() || self.slope_clearance < 0.0 {
            return Err(OperationError::InvalidInput(format!(
                "hip roof slope clearance must be non-negative, got {}",
                self.slope_clearance
            ))
            .into());
        }
        if !self.thickness.is_finite() || self.thickness < 0.0 {
            return Err(OperationError::InvalidInput(format!(
                "hip roof thickness must be non-negative, got {}",
                self.thickness
            ))
            .into());
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
        Ok(RoofPlan {
            skeleton,
            slope,
            eave,
            baseline_z: self.baseline_z + slope * self.slope_clearance,
        })
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
/// the skeleton cell chain, with `height` mapping each plan point and its
/// inset onto the sheet being built (top or thickness-lowered bottom).
fn face_polygon(
    cell: &SkeletonCell,
    eave: &[Point2],
    overhang: f64,
    height: &impl Fn(Point2, f64) -> Point3,
) -> Vec<Point3> {
    let edge = cell.edge_index;
    let mut polygon = Vec::with_capacity(cell.vertices.len() + 2);
    if overhang > 0.0 {
        let next = (edge + 1) % eave.len();
        polygon.push(height(eave[edge], -overhang));
        polygon.push(height(eave[next], -overhang));
        // Cell vertices are [start, end, chain from end back to start]:
        // with the eave band in front, traversal continues at the edge end
        // vertex and walks the chain back to the edge start vertex.
        for v in cell.vertices.iter().skip(1) {
            polygon.push(height(Point2::new(v.position.x, v.position.y), v.inset));
        }
        polygon.push(height(
            Point2::new(cell.vertices[0].position.x, cell.vertices[0].position.y),
            cell.vertices[0].inset,
        ));
    } else {
        for v in &cell.vertices {
            polygon.push(height(Point2::new(v.position.x, v.position.y), v.inset));
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
    use crate::error::GeolisError;
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

    /// The published pitch is `rise / max_inset` of the baseline skeleton
    /// — the same slope every sloped face rises at. A 4x4 square insets to
    /// 2, a 6x4 rectangle also to 2.
    #[test]
    fn slope_is_rise_over_max_inset() {
        let slope = MakeHipRoof::new(square(), 2.0).slope().unwrap();
        assert!((slope - 1.0).abs() < 1e-12, "square slope {slope}");
        let slope = MakeHipRoof::new(rect(), 1.0).slope().unwrap();
        assert!((slope - 0.5).abs() < 1e-12, "rect slope {slope}");
        // The overhang and the baseline z are pure post-skeleton offsets —
        // neither changes the pitch.
        let slope = MakeHipRoof::new(square(), 2.0)
            .with_overhang(0.7)
            .with_baseline_z(9.0)
            .slope()
            .unwrap();
        assert!((slope - 1.0).abs() < 1e-12, "offsets changed slope {slope}");
    }

    #[test]
    fn slope_rejects_the_same_input_execute_rejects() {
        assert!(MakeHipRoof::new(square(), 0.0).slope().is_err());
        assert!(MakeHipRoof::new(ring(&[(0.0, 0.0), (4.0, 0.0)]), 1.0)
            .slope()
            .is_err());
    }

    /// A slope clearance lifts the whole solid by `slope * clearance` — a
    /// pure z translation of both the eave and the ridge.
    #[test]
    fn slope_clearance_lifts_the_whole_solid() {
        // square 4x4, rise 2 -> slope 1. Clearance 0.5 -> lift 0.5.
        let mut store = TopologyStore::new();
        let solid = MakeHipRoof::new(square(), 2.0)
            .with_slope_clearance(0.5)
            .execute(&mut store)
            .unwrap();
        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();
        let (min_z, max_z) = mesh_z_range(&mesh);
        assert!((min_z - 0.5).abs() < 1e-9, "eave at {min_z}");
        assert!((max_z - 2.5).abs() < 1e-9, "ridge at {max_z}");
        assert_position_weld_watertight(&mesh);
    }

    /// The guarantee the parameter exists for, pinned exactly: with the
    /// overhang equal to the clearance, the eave ring lands back on the
    /// baseline z — i.e. the sloped surface at `clearance` outside the
    /// baseline is level with the baseline, which is what lets a wall top
    /// at the baseline z stay under the roof.
    #[test]
    fn surface_at_the_clearance_offset_returns_to_the_baseline_z() {
        let clearance = 0.09;
        let eave = MakeHipRoof::new(square(), 1.5)
            .with_overhang(clearance)
            .with_baseline_z(7.8)
            .with_slope_clearance(clearance)
            .eave_ring()
            .unwrap();
        for p in &eave {
            assert!(
                (p.z - 7.8).abs() < 1e-12,
                "eave must return to the baseline z, got {}",
                p.z
            );
        }

        // Control: without the lift the same surface dips below it.
        let flat = MakeHipRoof::new(square(), 1.5)
            .with_overhang(clearance)
            .with_baseline_z(7.8)
            .eave_ring()
            .unwrap();
        assert!(flat[0].z < 7.8 - 1e-6, "control eave at {}", flat[0].z);
    }

    /// A zero clearance must reproduce the unlifted solid bit for bit, so
    /// the parameter is a pure opt-in.
    #[test]
    fn zero_slope_clearance_is_the_unlifted_solid() {
        let mut store = TopologyStore::new();
        let plain = MakeHipRoof::new(l_shape(), 1.5)
            .with_overhang(0.3)
            .execute(&mut store)
            .unwrap();
        let plain_volume = Volume::new(plain).execute(&store).unwrap();
        let plain_mesh = TessellateSolid::new(plain, TessellationParams::default())
            .execute(&store)
            .unwrap();

        let mut store2 = TopologyStore::new();
        let zeroed = MakeHipRoof::new(l_shape(), 1.5)
            .with_overhang(0.3)
            .with_slope_clearance(0.0)
            .execute(&mut store2)
            .unwrap();
        let zeroed_volume = Volume::new(zeroed).execute(&store2).unwrap();
        let zeroed_mesh = TessellateSolid::new(zeroed, TessellationParams::default())
            .execute(&store2)
            .unwrap();

        assert!((plain_volume - zeroed_volume).abs() < 1e-12);
        assert_eq!(mesh_z_range(&plain_mesh), mesh_z_range(&zeroed_mesh));
    }

    /// The lift is a z translation, so it must leave the volume and the
    /// sloped surface area untouched.
    #[test]
    fn slope_clearance_preserves_volume_and_sloped_area() {
        let sloped_area = |clearance: f64| -> (f64, f64) {
            let mut store = TopologyStore::new();
            let solid = MakeHipRoof::new(l_shape(), 1.5)
                .with_overhang(0.3)
                .with_slope_clearance(clearance)
                .execute(&mut store)
                .unwrap();
            (
                Volume::new(solid).execute(&store).unwrap(),
                crate::operations::query::Area::new(solid)
                    .with_facing(crate::math::Vector3::z())
                    .execute(&store)
                    .unwrap(),
            )
        };
        let (v0, a0) = sloped_area(0.0);
        let (v1, a1) = sloped_area(0.4);
        assert!((v0 - v1).abs() < 1e-9, "volume {v0} vs {v1}");
        assert!((a0 - a1).abs() < 1e-9, "sloped area {a0} vs {a1}");
    }

    /// The up-facing area of a hip roof is exactly its roofing surface —
    /// the sloped faces, with the downward bottom cap excluded. A 4x4
    /// square at rise 2 (slope 1, i.e. 45 degrees) roofs 16 m2 in plan, so
    /// the sloped surface is `16 * sqrt(2)`.
    #[test]
    fn up_facing_area_is_the_roofing_surface() {
        let mut store = TopologyStore::new();
        let solid = MakeHipRoof::new(square(), 2.0).execute(&mut store).unwrap();
        let sloped = crate::operations::query::Area::new(solid)
            .with_facing(crate::math::Vector3::z())
            .execute(&store)
            .unwrap();
        let expected = 16.0 * std::f64::consts::SQRT_2;
        assert!(
            (sloped - expected).abs() < 1e-9,
            "sloped area {sloped} vs {expected}"
        );
        // The bottom cap is the 16 m2 plan square, so the total is the sum.
        let total = crate::operations::query::Area::new(solid)
            .execute(&store)
            .unwrap();
        assert!(
            (total - (expected + 16.0)).abs() < 1e-9,
            "total area {total}"
        );
    }

    /// A zero thickness must reproduce the solid wedge exactly, so the
    /// parameter is a pure opt-in.
    #[test]
    fn thickness_zero_matches_the_wedge() {
        let mut store = TopologyStore::new();
        let solid = MakeHipRoof::new(square(), 2.0)
            .with_thickness(0.0)
            .execute(&mut store)
            .unwrap();
        // The `square_pyramid` wedge: 4 sloped faces + bottom cap.
        assert_eq!(face_count(&store, solid), 5);
        assert!(edge_usage_is_two_everywhere(&store, solid));
        let volume = Volume::new(solid).execute(&store).unwrap();
        assert!((volume - 32.0 / 3.0).abs() < 1e-9, "volume {volume}");

        // And identical to the wedge built without the parameter at all.
        let mut store2 = TopologyStore::new();
        let plain = MakeHipRoof::new(square(), 2.0)
            .execute(&mut store2)
            .unwrap();
        assert_eq!(face_count(&store2, plain), face_count(&store, solid));
        let plain_volume = Volume::new(plain).execute(&store2).unwrap();
        assert!((volume - plain_volume).abs() < 1e-12);
        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();
        let plain_mesh = TessellateSolid::new(plain, TessellationParams::default())
            .execute(&store2)
            .unwrap();
        assert_eq!(mesh_z_range(&mesh), mesh_z_range(&plain_mesh));
    }

    /// The invariant the vertical offset buys: both sheets share the plan
    /// decomposition, so the shell encloses exactly `plan area * thickness`
    /// regardless of pitch or footprint convexity.
    #[test]
    fn shell_volume_is_plan_area_times_thickness() {
        let thickness = 0.12;

        // square 4x4, rise 2 -> slope 1; overhang 0.5 makes the eave
        // polygon the 5x5 square.
        let mut store = TopologyStore::new();
        let solid = MakeHipRoof::new(square(), 2.0)
            .with_overhang(0.5)
            .with_thickness(thickness)
            .execute(&mut store)
            .unwrap();
        let volume = Volume::new(solid).execute(&store).unwrap();
        let expected = 5.0 * 5.0 * thickness;
        assert!(
            (volume - expected).abs() < 1e-9,
            "square shell {volume} vs {expected}"
        );

        // Non-convex control: the L-shape's plan area is 6*6 - 3*3 = 27,
        // and with no overhang the eave polygon is the baseline itself.
        let mut store2 = TopologyStore::new();
        let l_solid = MakeHipRoof::new(l_shape(), 1.5)
            .with_thickness(thickness)
            .execute(&mut store2)
            .unwrap();
        let l_volume = Volume::new(l_solid).execute(&store2).unwrap();
        let l_expected = 27.0 * thickness;
        assert!(
            (l_volume - l_expected).abs() < 1e-9,
            "L shell {l_volume} vs {l_expected}"
        );
    }

    /// Top sheet + bottom sheet + eave fascia must close the solid: every
    /// topological edge used exactly twice, and the tessellated mesh free
    /// of boundary edges after a position weld.
    #[test]
    fn shell_is_watertight() {
        let mut store = TopologyStore::new();
        let solid = MakeHipRoof::new(l_shape(), 1.5)
            .with_overhang(0.3)
            .with_thickness(0.25)
            .execute(&mut store)
            .unwrap();
        // 6 top faces + 6 bottom faces + 6 fascia quads, and no bottom cap.
        assert_eq!(face_count(&store, solid), 18);
        assert!(edge_usage_is_two_everywhere(&store, solid));
        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert_position_weld_watertight(&mesh);
        // slope = 1.5 / 1.5 = 1: the eave is at -0.3, so the underside of
        // the fascia sits at -0.55 while the ridge stays at 1.5.
        let (min_z, max_z) = mesh_z_range(&mesh);
        assert!((min_z + 0.55).abs() < 1e-9, "underside at {min_z}");
        assert!((max_z - 1.5).abs() < 1e-9, "ridge at {max_z}");
    }

    #[test]
    fn negative_or_nan_thickness_is_a_typed_error() {
        let mut store = TopologyStore::new();
        for bad in [-0.1, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err = MakeHipRoof::new(square(), 2.0)
                .with_thickness(bad)
                .execute(&mut store)
                .unwrap_err();
            assert!(
                matches!(err, GeolisError::Operation(OperationError::InvalidInput(_))),
                "thickness {bad} gave {err:?}"
            );
        }
        // The check lives in the shared plan, so the pitch query and the
        // eave ring reject the same input.
        assert!(MakeHipRoof::new(square(), 2.0)
            .with_thickness(-1.0)
            .slope()
            .is_err());
        assert!(MakeHipRoof::new(square(), 2.0)
            .with_thickness(f64::NAN)
            .eave_ring()
            .is_err());
    }

    /// The slope-clearance lift is a pure z translation of the *whole*
    /// shell — top sheet, bottom sheet and fascia move together.
    #[test]
    fn clearance_lifts_the_whole_shell() {
        // square 4x4, rise 2 -> slope 1. Clearance 0.5 lifts the baseline
        // to z = 0.5, so the overhang of 0.25 puts the eave at 0.25 and the
        // fascia's underside at 0.15; the ridge sits at 0.5 + 2 = 2.5.
        let mut store = TopologyStore::new();
        let solid = MakeHipRoof::new(square(), 2.0)
            .with_overhang(0.25)
            .with_slope_clearance(0.5)
            .with_thickness(0.1)
            .execute(&mut store)
            .unwrap();
        assert!(edge_usage_is_two_everywhere(&store, solid));
        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert_position_weld_watertight(&mesh);
        let (min_z, max_z) = mesh_z_range(&mesh);
        assert!((min_z - 0.15).abs() < 1e-9, "underside at {min_z}");
        assert!((max_z - 2.5).abs() < 1e-9, "ridge at {max_z}");
        // A translation cannot change the enclosed volume: still the
        // 4.5x4.5 eave polygon times the thickness.
        let volume = Volume::new(solid).execute(&store).unwrap();
        let expected = 4.5 * 4.5 * 0.1;
        assert!(
            (volume - expected).abs() < 1e-9,
            "volume {volume} vs {expected}"
        );
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
        assert!(MakeHipRoof::new(square(), 2.0)
            .with_slope_clearance(-0.1)
            .execute(&mut store)
            .is_err());
        assert!(MakeHipRoof::new(square(), 2.0)
            .with_slope_clearance(f64::NAN)
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
