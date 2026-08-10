//! Vertical-sided prism whose caps are ruled (non-planar) surfaces.
//!
//! [`MakeRuledLoft`] is the sibling of [`MakeLoft`](super::MakeLoft) for
//! sloped bands: a ramp slab, or the raking rail of a stair that follows
//! a curved path. Both rings share one plan projection and the sides are
//! vertical; only the caps tilt, and along a CURVED band that tilt is a
//! warped surface, not a plane — three ring vertices fix a plane and the
//! fourth already leaves it. `MakeLoft` cannot build those bands because
//! it caps each ring with a single [`MakeFace`], which validates
//! coplanarity and rejects them outright.
//!
//! The cap here is therefore not one face but its triangulation: the
//! shared plan projection is ear-clipped ONCE
//! ([`triangulate_polygon_xy`]) and the resulting index triples are
//! lifted twice — each corner to its own bottom z and its own top z.
//! Every triangle is planar by construction, so the caps carry no
//! planarity requirement at all while the solid stays a plain `BRep` of
//! planar faces. The two caps share a triangulation, which is what keeps
//! a constant-thickness band's volume exactly `thickness × plan area`.
//!
//! Topology matches `MakeLoft`'s: `2n` side triangles (each vertical
//! quad split along a diagonal), `n − 2` triangles per cap, every edge
//! referenced by exactly two faces. Interior cap diagonals are created
//! once and shared by the two triangles that name them, and the ring
//! edges are shared with the side walls.
//!
//! # Faces carry the surface they are a piece OF
//!
//! Almost every face this op creates is a fragment of a larger surface
//! the caller actually designed, and the split is the BUILDER's choice,
//! not the design's:
//!
//! | Face | Is a piece of | Who chose the split |
//! |---|---|---|
//! | A cap triangle | that whole cap, one ruled surface | [`triangulate_polygon_xy`]'s ear order |
//! | A side triangle | the vertical quad over ring segment `i` | this op, splitting the quad on `b_i → t_j` |
//! | Side quad `i` itself | the source curve the caller's [`Self::with_segment_tags`] names | the caller's chord sampling |
//!
//! An edge whose existence and position an ALGORITHM picked can never be
//! a feature edge of the model, so every such face is bound as
//! [`FaceRole::Tagged`]`("{source}:p{ordinal}")` — one shared source
//! string per surface, one dense ordinal per piece to keep the names
//! unique. A consumer that strips the trailing `:p{ordinal}` recovers
//! the surface, which is how a wireframe extractor drops the interior
//! seams and keeps only the boundaries between genuinely different
//! surfaces.
//!
//! Naming is opt-in through [`Self::with_op_id`], exactly as
//! [`MakeSegmentedPrism`](crate::operations::creation::MakeSegmentedPrism)
//! does it: geolis never invents identity.

use std::collections::HashMap;

use crate::error::{OperationError, Result};
use crate::math::polygon_2d::{signed_area_2d, triangulate_polygon_xy};
use crate::math::{Point3, TOLERANCE};
use crate::operations::boolean_2d::WALL_EPS;
use crate::operations::creation::bind_created_face;
use crate::operations::creation::{MakeFace, MakeSolid};
use crate::topology::{
    EdgeId, FaceId, FaceRole, OpId, OrientedEdge, SegmentTag, ShellData, SolidId, TopologyStore,
    VertexData, VertexId,
};

use super::extrude::{create_closed_wire, create_line_edge, create_loop_edges};

/// Source string the BOTTOM cap's triangles share.
///
/// Namespaced under `loft:` so it cannot alias a caller segment tag —
/// see [`MakeRuledLoft::with_segment_tags`].
const BOTTOM_CAP_SOURCE: &str = "loft:cap:bottom";
/// Source string the TOP cap's triangles share.
const TOP_CAP_SOURCE: &str = "loft:cap:top";

/// The source string an UNTAGGED side quad's two triangles share: the
/// quad is still one planar surface the op split on a diagonal, even
/// when the caller cannot say which curve the quad came from.
fn untagged_side_source(segment: usize) -> String {
    format!("loft:side:{segment}")
}

/// Lofts two index-matched rings that share a plan projection into a
/// faceted solid with vertical sides and ruled (possibly non-planar)
/// caps.
///
/// See the module docs for why the caps are triangulated rather than
/// faced.
pub struct MakeRuledLoft {
    bottom: Vec<Point3>,
    top: Vec<Point3>,
    op_id: Option<OpId>,
    segment_tags: Option<Vec<SegmentTag>>,
}

impl MakeRuledLoft {
    /// Creates a ruled loft between `bottom` and `top`. The rings
    /// correspond by index; every invariant is checked in
    /// [`Self::execute`].
    #[must_use]
    pub fn new(bottom: Vec<Point3>, top: Vec<Point3>) -> Self {
        Self {
            bottom,
            top,
            op_id: None,
            segment_tags: None,
        }
    }

    /// Binds persistent [`FaceName::Created`](crate::topology::FaceName)
    /// names under `op`. Without it the solid carries no names at all —
    /// geolis never invents identity.
    #[must_use]
    pub fn with_op_id(mut self, op: OpId) -> Self {
        self.op_id = Some(op);
        self
    }

    /// Names the SOURCE CURVE behind each ring segment: `tags[i]` covers
    /// the segment from ring vertex `i` to vertex `i + 1`, so the list is
    /// exactly as long as either ring.
    ///
    /// Consecutive segments that chord ONE source curve must repeat one
    /// tag — that repetition is the whole point. The two triangles of
    /// every quad those segments raise then share a source string, and a
    /// wireframe extractor drops both the quad diagonals and the chord
    /// seams BETWEEN the quads, leaving only the boundary where the tag
    /// actually changes. Without tags each quad is its own surface and
    /// every chord seam survives as an edge.
    ///
    /// Tags are namespaced away from this op's intrinsic sources by the
    /// `loft:` prefix it reserves ([`BOTTOM_CAP_SOURCE`],
    /// [`TOP_CAP_SOURCE`], [`untagged_side_source`]); a caller tag that
    /// collides with one of those strings would merge unrelated faces.
    ///
    /// The length is validated in [`Self::execute`].
    #[must_use]
    pub fn with_segment_tags(mut self, tags: Vec<SegmentTag>) -> Self {
        self.segment_tags = Some(tags);
        self
    }

    /// Executes the loft, creating a closed solid in the topology store.
    ///
    /// # Errors
    ///
    /// Returns [`OperationError::InvalidInput`] when the rings have
    /// mismatched counts or fewer than 3 vertices, carry a non-finite
    /// coordinate, disagree in plan at any index (the sides are vertical
    /// by construction), lack strictly positive vertical separation at
    /// any index, repeat a vertex consecutively, or enclose no area in
    /// plan. Propagates [`OperationError::Failed`] from
    /// [`triangulate_polygon_xy`] when the plan projection is not a
    /// simple polygon, and any error raised while building the faces.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        let (bottom, top, reversed) = self.validated_rings()?;
        let n = bottom.len();
        let segment_tags = self.validated_segment_tags(n, reversed)?;

        // Triangulated before any store mutation, so a ring the ear
        // clipper rejects leaves no partial topology behind. One
        // triangulation serves both caps: the rings share a plan
        // projection, so the triples index both.
        let triangles = triangulate_polygon_xy(&bottom)?;

        let bottom_verts: Vec<VertexId> = bottom
            .iter()
            .map(|p| store.add_vertex(VertexData::new(*p)))
            .collect();
        let top_verts: Vec<VertexId> = top
            .iter()
            .map(|p| store.add_vertex(VertexData::new(*p)))
            .collect();

        let bottom_ring = create_loop_edges(store, &bottom_verts, &bottom)?;
        let top_ring = create_loop_edges(store, &top_verts, &top)?;
        // Vertical edges i → i and the diagonal splitting side quad i.
        let mut vert_edges = Vec::with_capacity(n);
        let mut diag_edges = Vec::with_capacity(n);
        for i in 0..n {
            let j = (i + 1) % n;
            vert_edges.push(create_line_edge(
                store,
                bottom_verts[i],
                top_verts[i],
                bottom[i],
                top[i],
            )?);
            diag_edges.push(create_line_edge(
                store,
                bottom_verts[i],
                top_verts[j],
                bottom[i],
                top[j],
            )?);
        }

        let mut faces = Vec::with_capacity(4 * n - 4);
        let bottom_cap = build_cap(
            store,
            &triangles,
            &bottom_verts,
            &bottom,
            &bottom_ring,
            Facing::Down,
        )?;
        let top_cap = build_cap(store, &triangles, &top_verts, &top, &top_ring, Facing::Up)?;
        let cap_faces = bottom_cap.len();
        faces.extend(bottom_cap);
        faces.extend(top_cap);

        // Side triangles: quad (b_i, b_j, t_j, t_i) split along b_i → t_j,
        // outward because the rings wind counter-clockwise in plan.
        for i in 0..n {
            let j = (i + 1) % n;
            let lower = create_closed_wire(
                store,
                vec![
                    OrientedEdge::new(bottom_ring[i], true),
                    OrientedEdge::new(vert_edges[j], true),
                    OrientedEdge::new(diag_edges[i], false),
                ],
            );
            faces.push(MakeFace::new(lower, vec![]).execute(store)?);

            let upper = create_closed_wire(
                store,
                vec![
                    OrientedEdge::new(diag_edges[i], true),
                    OrientedEdge::new(top_ring[i], false),
                    OrientedEdge::new(vert_edges[i], false),
                ],
            );
            faces.push(MakeFace::new(upper, vec![]).execute(store)?);
        }

        self.bind_names(store, &faces, cap_faces, segment_tags.as_deref());

        let shell = store.add_shell(ShellData {
            faces,
            is_closed: true,
        });
        MakeSolid::new(shell, vec![]).execute(store)
    }

    /// Binds one [`FaceRole::Tagged`] name per face, grouping the pieces
    /// of each surface under a shared source string (see the module
    /// docs). `faces` is `[bottom cap.., top cap.., side pairs..]` with
    /// `cap_faces` triangles in each cap; side quad `i` owns
    /// `faces[2 * cap_faces + 2 * i]` and the face after it.
    fn bind_names(
        &self,
        store: &mut TopologyStore,
        faces: &[FaceId],
        cap_faces: usize,
        segment_tags: Option<&[SegmentTag]>,
    ) {
        let Some(op) = &self.op_id else {
            return;
        };
        // One dense ordinal per source, so the names stay unique while
        // the shared prefix keeps saying which surface each piece is of.
        let mut pieces: HashMap<String, usize> = HashMap::new();
        let mut bind = |store: &mut TopologyStore, face: FaceId, source: &str| {
            let ordinal = pieces.entry(source.to_string()).or_default();
            let tag = SegmentTag::new(format!("{source}:p{ordinal}"));
            *ordinal += 1;
            bind_created_face(store, face, op, FaceRole::Tagged(tag));
        };

        let sides = &faces[2 * cap_faces..];
        for &face in &faces[..cap_faces] {
            bind(store, face, BOTTOM_CAP_SOURCE);
        }
        for &face in &faces[cap_faces..2 * cap_faces] {
            bind(store, face, TOP_CAP_SOURCE);
        }
        for (segment, pair) in sides.chunks(2).enumerate() {
            let source = segment_tags.map_or_else(
                || untagged_side_source(segment),
                |tags| tags[segment].as_str().to_string(),
            );
            for &face in pair {
                bind(store, face, &source);
            }
        }
    }

    /// The caller's segment tags, re-indexed onto the ring
    /// [`Self::validated_rings`] actually built.
    ///
    /// A clockwise ring is reversed there to put the side normals
    /// outward, which renumbers its segments: reversed segment `k` runs
    /// between the same two vertices as original segment `n − 2 − k`
    /// (mod `n`), so the tags must follow or every side would name the
    /// wrong curve.
    ///
    /// # Errors
    ///
    /// The tag list length does not match the ring's.
    fn validated_segment_tags(&self, n: usize, reversed: bool) -> Result<Option<Vec<SegmentTag>>> {
        let Some(tags) = &self.segment_tags else {
            return Ok(None);
        };
        if tags.len() != n {
            return Err(OperationError::InvalidInput(format!(
                "ruled loft has {n} ring segments but {} segment tags — one tag \
                 names one segment",
                tags.len()
            ))
            .into());
        }
        if !reversed {
            return Ok(Some(tags.clone()));
        }
        Ok(Some(
            (0..n).map(|k| tags[(n + n - 2 - k) % n].clone()).collect(),
        ))
    }

    /// Validates every invariant and returns the two rings normalized to
    /// a counter-clockwise plan winding, index correspondence intact,
    /// plus whether normalizing REVERSED them (which renumbers the ring
    /// segments — see [`Self::validated_segment_tags`]).
    fn validated_rings(&self) -> Result<(Vec<Point3>, Vec<Point3>, bool)> {
        let n = self.bottom.len();
        if n < 3 || self.top.len() != n {
            return Err(OperationError::InvalidInput(format!(
                "ruled loft rings need matching vertex counts >= 3, got {n} and {}",
                self.top.len()
            ))
            .into());
        }

        for (i, (b, t)) in self.bottom.iter().zip(&self.top).enumerate() {
            if ![b.x, b.y, b.z, t.x, t.y, t.z].iter().all(|v| v.is_finite()) {
                return Err(OperationError::InvalidInput(format!(
                    "ruled loft vertex {i} has a non-finite coordinate"
                ))
                .into());
            }
            let plan_drift = (t.x - b.x).abs().max((t.y - b.y).abs());
            if plan_drift > WALL_EPS {
                return Err(OperationError::InvalidInput(format!(
                    "ruled loft vertex {i} moves in plan by {plan_drift}: top ({}, {}) \
                     against bottom ({}, {}) — the sides are vertical by construction",
                    t.x, t.y, b.x, b.y
                ))
                .into());
            }
            if t.z - b.z <= TOLERANCE {
                return Err(OperationError::InvalidInput(format!(
                    "ruled loft vertex {i} has no vertical separation: top z {} must \
                     clear bottom z {}",
                    t.z, b.z
                ))
                .into());
            }
            let next = self.bottom[(i + 1) % n];
            if (next.x - b.x).abs().max((next.y - b.y).abs()) <= TOLERANCE {
                return Err(OperationError::InvalidInput(format!(
                    "ruled loft ring repeats vertex {i} at ({}, {}) — a zero-length \
                     ring segment has no side wall",
                    b.x, b.y
                ))
                .into());
            }
        }

        let area = signed_area_2d(&self.bottom);
        if area.abs() <= TOLERANCE {
            return Err(OperationError::InvalidInput(
                "ruled loft ring encloses no area in plan".into(),
            )
            .into());
        }
        // Counter-clockwise in plan puts the side normals outward and the
        // cap triangles' own winding in agreement with the ring's.
        if area > 0.0 {
            Ok((self.bottom.clone(), self.top.clone(), false))
        } else {
            Ok((
                self.bottom.iter().rev().copied().collect(),
                self.top.iter().rev().copied().collect(),
                true,
            ))
        }
    }
}

/// Which way a cap's faces look. The rings wind counter-clockwise in
/// plan, so the top cap keeps the triangulation's winding and the bottom
/// cap reverses it.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Facing {
    Up,
    Down,
}

/// Builds one warped cap as `n − 2` planar triangles over the ring's own
/// vertices, reusing `ring` for the boundary edges (shared with the side
/// walls) and creating each interior diagonal once — the two triangles
/// meeting at a diagonal name it with the same index pair, so they share
/// the edge and the shell stays watertight.
fn build_cap(
    store: &mut TopologyStore,
    triangles: &[[usize; 3]],
    verts: &[VertexId],
    points: &[Point3],
    ring: &[EdgeId],
    facing: Facing,
) -> Result<Vec<FaceId>> {
    let n = verts.len();
    let mut edges: HashMap<(usize, usize), EdgeId> = ring
        .iter()
        .enumerate()
        .map(|(i, &edge)| ((i, (i + 1) % n), edge))
        .collect();

    let mut faces = Vec::with_capacity(triangles.len());
    for tri in triangles {
        let corners = match facing {
            Facing::Up => *tri,
            Facing::Down => [tri[0], tri[2], tri[1]],
        };
        let mut wire_edges = Vec::with_capacity(3);
        for k in 0..3 {
            wire_edges.push(shared_edge(
                store,
                &mut edges,
                verts,
                points,
                corners[k],
                corners[(k + 1) % 3],
            )?);
        }
        let wire = create_closed_wire(store, wire_edges);
        faces.push(MakeFace::new(wire, vec![]).execute(store)?);
    }
    Ok(faces)
}

/// Returns the edge between ring vertices `a` and `b`, oriented from `a`
/// to `b`, creating it on first use.
fn shared_edge(
    store: &mut TopologyStore,
    edges: &mut HashMap<(usize, usize), EdgeId>,
    verts: &[VertexId],
    points: &[Point3],
    a: usize,
    b: usize,
) -> Result<OrientedEdge> {
    if let Some(&edge) = edges.get(&(a, b)) {
        return Ok(OrientedEdge::new(edge, true));
    }
    if let Some(&edge) = edges.get(&(b, a)) {
        return Ok(OrientedEdge::new(edge, false));
    }
    let edge = create_line_edge(store, verts[a], verts[b], points[a], points[b])?;
    edges.insert((a, b), edge);
    Ok(OrientedEdge::new(edge, true))
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::operations::query::Volume;
    use crate::operations::shaping::MakeLoft;
    use crate::tessellation::{TessellateSolid, TessellationParams};
    use crate::topology::{FaceName, FaceSurface};
    use std::f64::consts::FRAC_PI_2;

    /// Chords per quarter turn of the curved test band.
    const CHORDS: usize = 16;
    /// Constant vertical thickness of the test bands.
    const THICKNESS: f64 = 0.3;

    fn lift(ring: &[Point3], dz: f64) -> Vec<Point3> {
        ring.iter()
            .map(|p| Point3::new(p.x, p.y, p.z + dz))
            .collect()
    }

    fn shell_faces(store: &TopologyStore, solid: SolidId) -> Vec<FaceId> {
        let shell = store.solid(solid).unwrap().outer_shell;
        store.shell(shell).unwrap().faces.clone()
    }

    /// Every edge of the solid must be referenced by exactly two face
    /// wires — the shell is watertight only if the cap triangulation
    /// shares its interior diagonals and its ring edges with the walls.
    fn assert_every_edge_used_twice(store: &TopologyStore, solid: SolidId) {
        let mut counts: HashMap<EdgeId, usize> = HashMap::new();
        for face in shell_faces(store, solid) {
            let face_data = store.face(face).unwrap();
            assert!(face_data.inner_wires.is_empty());
            for oe in &store.wire(face_data.outer_wire).unwrap().edges {
                *counts.entry(oe.edge).or_insert(0) += 1;
            }
        }
        for (&edge, &count) in &counts {
            assert_eq!(count, 2, "edge {edge:?} referenced by {count} faces");
        }
    }

    /// Rectangular band sloping linearly in x: both caps stay planar, so
    /// `MakeLoft` accepts the same input and the two ops must agree.
    fn sloped_rectangle() -> Vec<Point3> {
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(4.0, 0.0, 2.0),
            Point3::new(4.0, 2.0, 2.0),
            Point3::new(0.0, 2.0, 0.0),
        ]
    }

    /// Annular sector band (quarter turn, `CHORDS` chords per arc) whose
    /// z rises linearly with the sweep angle: the cap is a faceted
    /// helicoid, non-planar at every interior vertex.
    fn helical_band() -> Vec<Point3> {
        let rise = 3.0;
        let point = |radius: f64, i: usize| {
            #[allow(clippy::cast_precision_loss)]
            let frac = i as f64 / CHORDS as f64;
            let angle = frac * FRAC_PI_2;
            Point3::new(radius * angle.cos(), radius * angle.sin(), frac * rise)
        };
        let mut ring: Vec<Point3> = (0..=CHORDS).map(|i| point(3.0, i)).collect();
        ring.extend((0..=CHORDS).rev().map(|i| point(2.0, i)));
        ring
    }

    // ── Planar-cap equivalence ─────────────────────────────────

    /// With planar (merely sloped) caps the ruled loft is the ordinary
    /// loft: `MakeLoft` accepts the input, and both solids enclose the
    /// same volume — plan area times the vertical separation.
    #[test]
    fn planar_sloped_caps_agree_with_the_planar_loft() {
        const SEPARATION: f64 = 1.0;
        let bottom = sloped_rectangle();
        let top = lift(&bottom, SEPARATION);

        let mut ruled_store = TopologyStore::new();
        let ruled = MakeRuledLoft::new(bottom.clone(), top.clone())
            .execute(&mut ruled_store)
            .unwrap();
        let ruled_volume = Volume::new(ruled).execute(&ruled_store).unwrap();

        let mut loft_store = TopologyStore::new();
        let lofted = MakeLoft::new(bottom.clone(), top)
            .execute(&mut loft_store)
            .unwrap();
        let loft_volume = Volume::new(lofted).execute(&loft_store).unwrap();

        let expected = signed_area_2d(&bottom).abs() * SEPARATION;
        assert!(
            (ruled_volume - expected).abs() < 1e-9,
            "ruled volume {ruled_volume} != {expected}"
        );
        assert!(
            (ruled_volume - loft_volume).abs() < 1e-9,
            "ruled volume {ruled_volume} != loft volume {loft_volume}"
        );

        // Same topology class: planar faces, every edge shared by two.
        let faces = shell_faces(&ruled_store, ruled);
        assert_eq!(faces.len(), 4 * bottom.len() - 4, "2 caps + 2n sides");
        for face in faces {
            assert!(matches!(
                ruled_store.face(face).unwrap().surface,
                FaceSurface::Plane(_)
            ));
        }
        assert_every_edge_used_twice(&ruled_store, ruled);
    }

    // ── Warped caps ────────────────────────────────────────────

    /// The curved band's caps are genuinely non-planar: `MakeLoft`
    /// rejects them, which is the whole reason this op exists.
    #[test]
    fn the_planar_loft_rejects_the_warped_band() {
        let bottom = helical_band();
        let top = lift(&bottom, THICKNESS);
        let mut store = TopologyStore::new();
        assert!(MakeLoft::new(bottom, top).execute(&mut store).is_err());
    }

    #[test]
    fn warped_band_builds_a_watertight_solid() {
        let bottom = helical_band();
        let top = lift(&bottom, THICKNESS);
        let n = bottom.len();
        assert_eq!(n, 2 * (CHORDS + 1), "outer arc + inner arc");

        let mut store = TopologyStore::new();
        let solid = MakeRuledLoft::new(bottom.clone(), top.clone())
            .execute(&mut store)
            .unwrap();

        let faces = shell_faces(&store, solid);
        assert_eq!(
            faces.len(),
            4 * n - 4,
            "2n side triangles + (n - 2) per cap"
        );
        assert_every_edge_used_twice(&store, solid);

        // Constant thickness over the ear-clipped plan: the two caps
        // share a triangulation, so the volume is exact.
        let volume = Volume::new(solid).execute(&store).unwrap();
        let expected = signed_area_2d(&bottom).abs() * THICKNESS;
        assert!(
            (volume - expected).abs() < 1e-9,
            "volume {volume} != {expected}"
        );
    }

    /// Every cap vertex keeps the z it was given — the caps are not
    /// flattened onto a fitted plane anywhere.
    #[test]
    fn warped_band_cap_vertices_keep_their_input_z() {
        let bottom = helical_band();
        let top = lift(&bottom, THICKNESS);

        let mut store = TopologyStore::new();
        let solid = MakeRuledLoft::new(bottom.clone(), top.clone())
            .execute(&mut store)
            .unwrap();

        // Collect the solid's distinct vertex positions.
        let mut positions: Vec<Point3> = Vec::new();
        for face in shell_faces(&store, solid) {
            let wire = store.face(face).unwrap().outer_wire;
            for oe in &store.wire(wire).unwrap().edges {
                let edge = store.edge(oe.edge).unwrap();
                for vertex in [edge.start, edge.end] {
                    positions.push(store.vertex(vertex).unwrap().point);
                }
            }
        }
        for expected in bottom.iter().chain(&top) {
            assert!(
                positions.iter().any(|p| (p - expected).norm() < 1e-12),
                "cap vertex {expected:?} is missing from the solid"
            );
        }

        // Spot-check the interior samples that make the cap non-planar:
        // three of them already fix a plane the rest leave.
        for i in [1, CHORDS / 2, CHORDS - 1] {
            let sample = bottom[i];
            let z = positions
                .iter()
                .find(|p| (p.x - sample.x).abs() < 1e-12 && (p.y - sample.y).abs() < 1e-12)
                .map(|p| p.z)
                .unwrap();
            assert!((z - sample.z).abs() < 1e-12, "vertex {i} moved to z {z}");
        }
    }

    /// Interns a mesh vertex by quantized position, welding the separate
    /// per-face vertices that sit at the same point.
    #[allow(clippy::cast_possible_truncation)]
    fn weld(ids: &mut HashMap<(i64, i64, i64), u32>, p: &Point3) -> u32 {
        const Q: f64 = 1e6;
        let key = (
            (p.x * Q).round() as i64,
            (p.y * Q).round() as i64,
            (p.z * Q).round() as i64,
        );
        let next = ids.len() as u32;
        *ids.entry(key).or_insert(next)
    }

    /// The band must also mesh watertight: after welding mesh vertices by
    /// position, every undirected triangle edge is used exactly twice.
    #[test]
    fn warped_band_tessellates_watertight() {
        let bottom = helical_band();
        let top = lift(&bottom, THICKNESS);
        let mut store = TopologyStore::new();
        let solid = MakeRuledLoft::new(bottom, top).execute(&mut store).unwrap();

        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty());

        let mut ids: HashMap<(i64, i64, i64), u32> = HashMap::new();
        let mut counts: HashMap<(u32, u32), usize> = HashMap::new();
        for tri in &mesh.indices {
            let a = weld(&mut ids, &mesh.vertices[tri[0] as usize]);
            let b = weld(&mut ids, &mesh.vertices[tri[1] as usize]);
            let c = weld(&mut ids, &mesh.vertices[tri[2] as usize]);
            for (x, y) in [(a, b), (b, c), (c, a)] {
                *counts
                    .entry(if x < y { (x, y) } else { (y, x) })
                    .or_insert(0) += 1;
            }
        }
        let boundary = counts.values().filter(|&&c| c != 2).count();
        assert_eq!(boundary, 0, "found {boundary} boundary edges in the mesh");
    }

    // ── Orientation ────────────────────────────────────────────

    /// The traversal points of a face's outer wire, in order.
    fn wire_points(store: &TopologyStore, face: FaceId) -> Vec<Point3> {
        let wire = store.face(face).unwrap().outer_wire;
        store
            .wire(wire)
            .unwrap()
            .edges
            .iter()
            .map(|oe| {
                let edge = store.edge(oe.edge).unwrap();
                let vertex = if oe.forward { edge.start } else { edge.end };
                store.vertex(vertex).unwrap().point
            })
            .collect()
    }

    /// Bottom cap down, top cap up, sides outward — pinned two ways:
    /// each cap face's plane normal by its z sign, and the whole shell by
    /// the divergence theorem, whose signed volume comes out POSITIVE
    /// only when every face winds outward. (The `Volume` query takes an
    /// absolute value, so it cannot see a flipped face.)
    #[test]
    fn caps_face_out_and_every_face_winds_outward() {
        for bottom in [sloped_rectangle(), helical_band()] {
            let n = bottom.len();
            let top = lift(&bottom, THICKNESS);
            let mut store = TopologyStore::new();
            let solid = MakeRuledLoft::new(bottom.clone(), top)
                .execute(&mut store)
                .unwrap();
            let faces = shell_faces(&store, solid);

            // Caps come first: bottom triangles, then top triangles.
            for (i, &face) in faces[..2 * (n - 2)].iter().enumerate() {
                let FaceSurface::Plane(plane) = &store.face(face).unwrap().surface else {
                    panic!("cap faces are planar");
                };
                let normal_z = plane.plane_normal().z;
                if i < n - 2 {
                    assert!(normal_z < 0.0, "bottom cap face {i} looks up ({normal_z})");
                } else {
                    assert!(normal_z > 0.0, "top cap face {i} looks down ({normal_z})");
                }
            }

            let mut signed_volume = 0.0;
            for face in faces {
                let pts = wire_points(&store, face);
                assert_eq!(pts.len(), 3, "every face of a ruled loft is a triangle");
                signed_volume += pts[0].coords.dot(&pts[1].coords.cross(&pts[2].coords));
            }
            signed_volume /= 6.0;
            let expected = signed_area_2d(&bottom).abs() * THICKNESS;
            assert!(
                (signed_volume - expected).abs() < 1e-9,
                "outward-wound volume {signed_volume} != {expected}"
            );
        }
    }

    // ── Rejections ─────────────────────────────────────────────

    fn is_invalid_input(result: &Result<SolidId>) -> bool {
        matches!(
            result,
            Err(crate::error::GeolisError::Operation(
                OperationError::InvalidInput(_)
            ))
        )
    }

    #[test]
    fn mismatched_vertex_counts_are_rejected() {
        let bottom = sloped_rectangle();
        let mut top = lift(&bottom, 1.0);
        top.pop();
        let mut store = TopologyStore::new();
        assert!(is_invalid_input(
            &MakeRuledLoft::new(bottom, top).execute(&mut store)
        ));
    }

    #[test]
    fn rings_shorter_than_three_vertices_are_rejected() {
        let bottom = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
        let top = lift(&bottom, 1.0);
        let mut store = TopologyStore::new();
        assert!(is_invalid_input(
            &MakeRuledLoft::new(bottom, top).execute(&mut store)
        ));
    }

    #[test]
    fn a_vertex_without_vertical_separation_is_rejected() {
        let bottom = sloped_rectangle();
        let mut top = lift(&bottom, 1.0);
        top[2].z = bottom[2].z;
        let mut store = TopologyStore::new();
        assert!(is_invalid_input(
            &MakeRuledLoft::new(bottom.clone(), top).execute(&mut store)
        ));

        // A top BELOW the bottom is the same rejection.
        let mut inverted = lift(&bottom, 1.0);
        inverted[1].z = bottom[1].z - 0.5;
        let mut store = TopologyStore::new();
        assert!(is_invalid_input(
            &MakeRuledLoft::new(bottom, inverted).execute(&mut store)
        ));
    }

    #[test]
    fn a_vertex_that_moves_in_plan_is_rejected() {
        let bottom = sloped_rectangle();
        let mut top = lift(&bottom, 1.0);
        top[1].x += 0.25;
        let mut store = TopologyStore::new();
        assert!(is_invalid_input(
            &MakeRuledLoft::new(bottom, top).execute(&mut store)
        ));
    }

    #[test]
    fn a_non_finite_coordinate_is_rejected() {
        let bottom = sloped_rectangle();
        let mut top = lift(&bottom, 1.0);
        top[0].x = f64::NAN;
        let mut store = TopologyStore::new();
        assert!(is_invalid_input(
            &MakeRuledLoft::new(bottom, top).execute(&mut store)
        ));
    }

    // ── Degenerate rings ───────────────────────────────────────

    /// `MakeLoft` fails on a repeated consecutive vertex (its ring edge
    /// has zero length); so does this op, with a named error.
    #[test]
    fn a_repeated_consecutive_vertex_is_rejected() {
        let mut bottom = sloped_rectangle();
        bottom.insert(1, bottom[0]);
        let top = lift(&bottom, 1.0);

        let mut store = TopologyStore::new();
        assert!(is_invalid_input(
            &MakeRuledLoft::new(bottom.clone(), top.clone()).execute(&mut store)
        ));

        let mut store = TopologyStore::new();
        assert!(MakeLoft::new(bottom, top).execute(&mut store).is_err());
    }

    #[test]
    fn a_ring_enclosing_no_plan_area_is_rejected() {
        let bottom = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.5),
            Point3::new(2.0, 0.0, 1.0),
        ];
        let top = lift(&bottom, 1.0);
        let mut store = TopologyStore::new();
        assert!(is_invalid_input(
            &MakeRuledLoft::new(bottom, top).execute(&mut store)
        ));
    }

    // ── Face naming: every face says which surface it is a piece of ──

    /// The source string a face's `Tagged` name carries, with the
    /// trailing `:p{ordinal}` piece suffix stripped.
    fn face_source(store: &TopologyStore, face: FaceId) -> String {
        match store.names().name_of_face(face) {
            Some(FaceName::Created {
                role: FaceRole::Tagged(tag),
                ..
            }) => {
                let tag = tag.as_str();
                let Some(at) = tag.rfind(":p") else {
                    panic!("`{tag}` carries no piece ordinal")
                };
                assert!(
                    tag[at + 2..].chars().all(|c| c.is_ascii_digit()) && at + 2 < tag.len(),
                    "`{tag}` must end in a numeric piece ordinal",
                );
                tag[..at].to_string()
            }
            other => panic!("face {face:?} carries no tagged name: {other:?}"),
        }
    }

    /// Without an op id the solid stays anonymous — geolis never invents
    /// identity, so an unnamed caller keeps the pre-naming behaviour.
    #[test]
    fn no_op_id_binds_no_names() {
        let bottom = helical_band();
        let top = lift(&bottom, THICKNESS);
        let mut store = TopologyStore::new();
        let solid = MakeRuledLoft::new(bottom, top).execute(&mut store).unwrap();
        for face in shell_faces(&store, solid) {
            assert!(
                store.names().name_of_face(face).is_none(),
                "face {face:?} was named without an op id",
            );
        }
    }

    /// Every face is a piece of some surface, and the pieces of one
    /// surface share a source: the two caps take one source each, and an
    /// untagged side quad's two triangles take one between them. Every
    /// bound name is still unique.
    #[test]
    fn every_face_names_the_surface_it_is_a_piece_of() {
        let bottom = helical_band();
        let n = bottom.len();
        let top = lift(&bottom, THICKNESS);
        let mut store = TopologyStore::new();
        let solid = MakeRuledLoft::new(bottom, top)
            .with_op_id(OpId::new("ramp"))
            .execute(&mut store)
            .unwrap();

        let faces = shell_faces(&store, solid);
        let mut per_source: HashMap<String, usize> = HashMap::new();
        let mut names: Vec<String> = Vec::new();
        for &face in &faces {
            *per_source.entry(face_source(&store, face)).or_default() += 1;
            let Some(FaceName::Created {
                role: FaceRole::Tagged(tag),
                ..
            }) = store.names().name_of_face(face)
            else {
                unreachable!("face_source already asserted the shape")
            };
            names.push(tag.as_str().to_string());
        }
        names.sort_unstable();
        let unique = names.len();
        names.dedup();
        assert_eq!(names.len(), unique, "piece ordinals must keep names unique");

        // Each cap is `n - 2` triangles of one surface.
        assert_eq!(per_source.get(BOTTOM_CAP_SOURCE), Some(&(n - 2)));
        assert_eq!(per_source.get(TOP_CAP_SOURCE), Some(&(n - 2)));
        // Each of the `n` side quads is two triangles of one surface.
        for segment in 0..n {
            assert_eq!(
                per_source.get(&untagged_side_source(segment)),
                Some(&2),
                "side quad {segment} must own exactly its two triangles",
            );
        }
        assert_eq!(per_source.len(), n + 2, "no other surface exists");
    }

    /// Consecutive ring segments that chord ONE source curve repeat one
    /// tag, and every quad they raise then lands under that one source —
    /// which is what lets a wireframe drop the chord seams between them.
    #[test]
    fn segments_sharing_a_tag_land_under_one_source() {
        let bottom = helical_band();
        let n = bottom.len();
        let top = lift(&bottom, THICKNESS);
        // The band is [outer arc, end cap, inner arc, start cap]: the two
        // arcs are one curve each, the two caps one segment each.
        let tags: Vec<SegmentTag> = (0..n)
            .map(|i| {
                SegmentTag::new(if i < CHORDS {
                    "outer-arc"
                } else if i == CHORDS {
                    "end-cap"
                } else if i < 2 * CHORDS + 1 {
                    "inner-arc"
                } else {
                    "start-cap"
                })
            })
            .collect();

        let mut store = TopologyStore::new();
        let solid = MakeRuledLoft::new(bottom, top)
            .with_op_id(OpId::new("ramp"))
            .with_segment_tags(tags)
            .execute(&mut store)
            .unwrap();

        let mut per_source: HashMap<String, usize> = HashMap::new();
        for face in shell_faces(&store, solid) {
            *per_source.entry(face_source(&store, face)).or_default() += 1;
        }
        // Two triangles per quad, and each arc raises CHORDS quads.
        assert_eq!(per_source.get("outer-arc"), Some(&(2 * CHORDS)));
        assert_eq!(per_source.get("inner-arc"), Some(&(2 * CHORDS)));
        assert_eq!(per_source.get("end-cap"), Some(&2));
        assert_eq!(per_source.get("start-cap"), Some(&2));
        assert_eq!(
            per_source.len(),
            6,
            "four side sources plus the two caps: {per_source:?}",
        );
    }

    /// Normalizing a CLOCKWISE ring reverses it, which renumbers the
    /// segments. The tags must follow, so a face still names the curve
    /// its own geometry came from — never the one that happens to sit at
    /// its new index.
    #[test]
    fn a_clockwise_ring_keeps_every_tag_on_its_own_segment() {
        // A trapezoid, so every side has a distinct plan direction.
        let ccw = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(4.0, 0.0, 1.0),
            Point3::new(3.0, 2.0, 1.0),
            Point3::new(1.0, 2.0, 0.0),
        ];
        let tags: Vec<SegmentTag> = ["south", "east", "north", "west"]
            .iter()
            .map(|t| SegmentTag::new(*t))
            .collect();

        // The same four sides, wound the other way: segment (v_i, v_i+1)
        // of the CCW ring is segment (n-2-k) of the reversed one.
        let n = ccw.len();
        let cw: Vec<Point3> = ccw.iter().rev().copied().collect();
        let cw_tags: Vec<SegmentTag> = (0..n).map(|k| tags[(n + n - 2 - k) % n].clone()).collect();

        // Each run maps "plan midpoint of a side quad" → its source.
        let sides_by_midpoint = |ring: Vec<Point3>, tags: Vec<SegmentTag>| {
            let top = lift(&ring, THICKNESS);
            let mut store = TopologyStore::new();
            let solid = MakeRuledLoft::new(ring, top)
                .with_op_id(OpId::new("trapezoid"))
                .with_segment_tags(tags)
                .execute(&mut store)
                .unwrap();
            let mut out: Vec<(String, String)> = Vec::new();
            for face in shell_faces(&store, solid) {
                let source = face_source(&store, face);
                if source.starts_with("loft:cap:") {
                    continue;
                }
                // A side triangle spans exactly two plan positions.
                let mut plan: Vec<(i64, i64)> = Vec::new();
                let wire = store.face(face).unwrap().outer_wire;
                for oe in &store.wire(wire).unwrap().edges {
                    let edge = store.edge(oe.edge).unwrap();
                    for v in [edge.start, edge.end] {
                        let p = store.vertex(v).unwrap().point;
                        #[allow(clippy::cast_possible_truncation)]
                        let key = ((p.x * 1e6) as i64, (p.y * 1e6) as i64);
                        if !plan.contains(&key) {
                            plan.push(key);
                        }
                    }
                }
                plan.sort_unstable();
                assert_eq!(plan.len(), 2, "a side triangle spans one plan segment");
                out.push((source, format!("{plan:?}")));
            }
            out.sort();
            out
        };

        assert_eq!(
            sides_by_midpoint(ccw, tags),
            sides_by_midpoint(cw, cw_tags),
            "reversing the ring must carry each tag to the same geometry",
        );
    }

    /// One tag names one segment: a list of any other length is refused
    /// rather than silently shifting every later side onto the wrong
    /// curve.
    #[test]
    fn a_tag_list_that_does_not_match_the_ring_is_refused() {
        let bottom = sloped_rectangle();
        let top = lift(&bottom, THICKNESS);
        let mut store = TopologyStore::new();
        assert!(is_invalid_input(
            &MakeRuledLoft::new(bottom, top)
                .with_op_id(OpId::new("short"))
                .with_segment_tags(vec![SegmentTag::new("only-one")])
                .execute(&mut store)
        ));
    }
}
