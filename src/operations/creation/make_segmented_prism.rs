//! Segmented prism — one side face per profile segment (F5 Phase A).
//!
//! [`MakeSegmentedPrism`] extrudes a closed chain of line / arc segments
//! along a direction, producing **one side face per segment** (arcs become
//! exact cylindrical patches). Adjacent side faces share a vertical kink
//! edge at each segment joint, and every side face shares its bottom / top
//! boundary edge with the corresponding planar cap (F2 shared-edge
//! topology + per-face pcurves), so the solid tessellates watertight by
//! construction: every boundary edge is sampled once by the
//! [`EdgeSampleCache`](crate::tessellation) and consumed identically by
//! both incident faces.
//!
//! Naming (op id present): side face `k` binds
//! [`FaceRole::Tagged`]`(tag_k)` when the caller supplies segment tags
//! (junction-stable identity — geolis never invents it, the `OpId`
//! precedent), or positional [`FaceRole::Side`]`(k)` otherwise; the caps
//! bind `CapStart` / `CapEnd`.
//!
//! **No edge names are bound in Phase A**: the prism has N bottom and N top
//! boundary edges (one per segment), so the whole-ring
//! `EdgeRole::RingStart` / `RingEnd` names would be ambiguous. Phase B
//! decides the per-segment edge-naming scheme together with the multi-face
//! band names.

use crate::error::{OperationError, Result};
use crate::geometry::nurbs::{KnotVector, NurbsCurve2D, NurbsCurve3D, NurbsSurface};
use crate::math::{Point2, Point3, Vector3, TOLERANCE};
use crate::topology::{
    EdgeCurve, EdgeData, EdgeId, FaceId, FacePcurve, FaceRole, FaceSurface, OpId, OrientedEdge,
    SegmentTag, SolidId, TopologyStore, VertexData, VertexId, WireData,
};

use super::make_nurbs_solid::{bind_created_face, finish_solid, iso_pcurve_u};
use super::{MakeFace, MakeNurbsFace};

/// Interior samples per segment when computing the profile's signed area
/// (winding / orientation test only — never used for boundary geometry).
const WINDING_SAMPLES: usize = 16;

/// One segment of a closed planar profile chain.
#[derive(Debug, Clone)]
pub enum ProfileSegment {
    /// A straight segment from `start` to `end`.
    Line {
        /// Segment start point.
        start: Point3,
        /// Segment end point.
        end: Point3,
    },
    /// An exact rational-quadratic arc in center / normal / angles form,
    /// matching [`NurbsCurve3D::arc`] (angles in radians measured from
    /// `ref_dir` toward `normal × ref_dir`; the sweep must be positive).
    Arc {
        /// Arc center.
        center: Point3,
        /// Arc radius.
        radius: f64,
        /// Plane normal of the arc.
        normal: Vector3,
        /// Zero-angle direction (perpendicular to `normal`).
        ref_dir: Vector3,
        /// Start angle in radians.
        start_angle: f64,
        /// End angle in radians (must exceed `start_angle`).
        end_angle: f64,
    },
}

impl ProfileSegment {
    /// Builds the exact 3D NURBS curve of this segment (degree-1 line or
    /// rational-quadratic arc).
    fn curve(&self) -> Result<NurbsCurve3D> {
        match self {
            Self::Line { start, end } => NurbsCurve3D::polyline(&[*start, *end]),
            Self::Arc {
                center,
                radius,
                normal,
                ref_dir,
                start_angle,
                end_angle,
            } => NurbsCurve3D::arc(
                *center,
                *radius,
                *normal,
                *ref_dir,
                *start_angle,
                *end_angle,
            ),
        }
    }
}

/// A prism over a closed segment chain with one side face per segment.
///
/// See the module docs for the shared-edge topology and naming contract.
pub struct MakeSegmentedPrism {
    profile: Vec<ProfileSegment>,
    direction: Vector3,
    op_id: Option<OpId>,
    segment_tags: Option<Vec<SegmentTag>>,
}

impl MakeSegmentedPrism {
    /// Creates a segmented prism extruding `profile` along `direction`.
    #[must_use]
    pub fn new(profile: Vec<ProfileSegment>, direction: Vector3) -> Self {
        Self {
            profile,
            direction,
            op_id: None,
            segment_tags: None,
        }
    }

    /// Registers persistent names for the prism's faces under the
    /// caller-supplied operation identity.
    #[must_use]
    pub fn with_op_id(mut self, op: OpId) -> Self {
        self.op_id = Some(op);
        self
    }

    /// Supplies one junction-stable tag per profile segment (parallel to the
    /// profile). Side face `k` then binds [`FaceRole::Tagged`]`(tags[k])`
    /// instead of the positional [`FaceRole::Side`]`(k)`.
    #[must_use]
    pub fn with_segment_tags(mut self, tags: Vec<SegmentTag>) -> Self {
        self.segment_tags = Some(tags);
        self
    }

    /// Builds the prism solid in the store.
    ///
    /// # Errors
    ///
    /// Returns a typed error if `direction` is zero-length, the profile has
    /// fewer than 3 segments, the segment chain is not endpoint-closed
    /// (consecutive endpoints must coincide within [`TOLERANCE`]), the tag
    /// count does not match the segment count, an untagged profile exceeds
    /// the positional `Side(u8)` range, the profile encloses no area
    /// transverse to `direction`, or any curve / surface / face construction
    /// fails.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        let curves = self.validate()?;
        let n = curves.len();

        // Extruded-surface normal (∂u × ∂v = tangent × direction) points
        // outward exactly when the profile winds counter-clockwise about the
        // extrusion direction; the winding is a global property, so one flip
        // decision covers every side face (non-convex profiles included).
        // Computed before any store mutation so a degenerate profile fails
        // without leaving partial topology behind.
        let ccw = profile_signed_area(&curves, &self.direction)? > 0.0;

        // One extruded surface per segment; exact for rational profiles
        // (v = 0 row is the segment curve, v = 1 its translate).
        let mut surfaces = Vec::with_capacity(n);
        for curve in &curves {
            surfaces.push(NurbsSurface::extrude(curve, self.direction)?);
        }

        let edges = build_shared_edges(store, &curves, &surfaces, self.direction)?;

        // Side faces: one per segment, each a full-domain extruded surface
        // whose 4-edge wire references only shared edges, with pcurves for
        // all four boundary edges.
        let mut side_faces = Vec::with_capacity(n);
        for (k, surface) in surfaces.iter().enumerate() {
            side_faces.push(build_side_face(store, surface, &edges, k, ccw)?);
        }

        // Caps: one planar face per end, each with a single wire around ALL
        // per-segment boundary edges (the same edges the side faces use).
        let bottom_cap = cap_over_edges(store, &edges.bottom, -self.direction)?;
        let top_cap = cap_over_edges(store, &edges.top, self.direction)?;

        if let Some(op) = &self.op_id {
            for (k, &face) in side_faces.iter().enumerate() {
                let role = match &self.segment_tags {
                    Some(tags) => FaceRole::Tagged(tags[k].clone()),
                    // Range validated up front in `validate`.
                    #[allow(clippy::cast_possible_truncation)]
                    None => FaceRole::Side(k as u8),
                };
                bind_created_face(store, face, op, role);
            }
            bind_created_face(store, bottom_cap, op, FaceRole::CapStart);
            bind_created_face(store, top_cap, op, FaceRole::CapEnd);
        }

        let mut faces = side_faces;
        faces.push(bottom_cap);
        faces.push(top_cap);
        Ok(finish_solid(store, faces))
    }

    /// Validates the inputs and builds the per-segment curves.
    fn validate(&self) -> Result<Vec<NurbsCurve3D>> {
        if self.direction.norm() <= TOLERANCE {
            return Err(OperationError::InvalidInput(
                "segmented prism direction must be non-zero".into(),
            )
            .into());
        }
        let n = self.profile.len();
        if n < 3 {
            return Err(OperationError::InvalidInput(format!(
                "segmented prism profile needs at least 3 segments, got {n}"
            ))
            .into());
        }
        if let Some(tags) = &self.segment_tags {
            if tags.len() != n {
                return Err(OperationError::InvalidInput(format!(
                    "segment tag count {} does not match profile segment count {n}",
                    tags.len()
                ))
                .into());
            }
        } else if self.op_id.is_some() && n > usize::from(u8::MAX) + 1 {
            return Err(OperationError::InvalidInput(format!(
                "untagged segmented prism supports at most 256 named segments \
                 (positional Side(u8)); got {n} — supply segment tags"
            ))
            .into());
        }

        let curves: Vec<NurbsCurve3D> = self
            .profile
            .iter()
            .map(ProfileSegment::curve)
            .collect::<Result<_>>()?;

        // The chain must be endpoint-closed: each segment's end coincides
        // with the next segment's start (cyclically) within TOLERANCE.
        for k in 0..n {
            let next = (k + 1) % n;
            let (_, t1) = curves[k].parameter_domain();
            let end = curves[k].point_at(t1)?;
            let (t0, _) = curves[next].parameter_domain();
            let start = curves[next].point_at(t0)?;
            if (end - start).norm() > TOLERANCE {
                return Err(OperationError::InvalidInput(format!(
                    "profile segment chain is not closed between segment {k} and {next}"
                ))
                .into());
            }
        }
        Ok(curves)
    }
}

/// The prism's shared boundary edges, indexed per segment: `bottom[k]` /
/// `top[k]` run along segment `k` (shared between side face `k` and the
/// corresponding cap), `kink[k]` rises at joint `k` (shared between side
/// faces `k − 1` and `k`).
struct SharedEdges {
    bottom: Vec<EdgeId>,
    top: Vec<EdgeId>,
    kink: Vec<EdgeId>,
}

/// Builds the shared joint vertices and the `3n` shared edges of the prism.
fn build_shared_edges(
    store: &mut TopologyStore,
    curves: &[NurbsCurve3D],
    surfaces: &[NurbsSurface],
    direction: Vector3,
) -> Result<SharedEdges> {
    let n = curves.len();

    // Shared joint vertices: bottom / top ring vertices at each segment
    // start point.
    let mut bottom_joints: Vec<VertexId> = Vec::with_capacity(n);
    let mut top_joints: Vec<VertexId> = Vec::with_capacity(n);
    let mut joint_points: Vec<Point3> = Vec::with_capacity(n);
    for curve in curves {
        let (t0, _) = curve.parameter_domain();
        let point = curve.point_at(t0)?;
        bottom_joints.push(store.add_vertex(VertexData::new(point)));
        top_joints.push(store.add_vertex(VertexData::new(point + direction)));
        joint_points.push(point);
    }

    // Per-segment bottom / top boundary edges (shared with the caps) and one
    // vertical kink edge per joint (shared by the two adjacent side faces).
    let mut edges = SharedEdges {
        bottom: Vec::with_capacity(n),
        top: Vec::with_capacity(n),
        kink: Vec::with_capacity(n),
    };
    for k in 0..n {
        let next = (k + 1) % n;
        let (t0, t1) = curves[k].parameter_domain();
        edges.bottom.push(store.add_edge(EdgeData {
            start: bottom_joints[k],
            end: bottom_joints[next],
            curve: EdgeCurve::Nurbs(curves[k].clone()),
            t_start: t0,
            t_end: t1,
        }));
        // The surface's exact v = 1 isocurve is the segment curve translated
        // by `direction`, with identical knots.
        let (_, (_, v1)) = surfaces[k].parameter_domain();
        let top_curve = surfaces[k].isocurve_v(v1)?;
        let (tt0, tt1) = top_curve.parameter_domain();
        edges.top.push(store.add_edge(EdgeData {
            start: top_joints[k],
            end: top_joints[next],
            curve: EdgeCurve::Nurbs(top_curve),
            t_start: tt0,
            t_end: tt1,
        }));
        let kink_curve = NurbsCurve3D::polyline(&[joint_points[k], joint_points[k] + direction])?;
        let (kt0, kt1) = kink_curve.parameter_domain();
        edges.kink.push(store.add_edge(EdgeData {
            start: bottom_joints[k],
            end: top_joints[k],
            curve: EdgeCurve::Nurbs(kink_curve),
            t_start: kt0,
            t_end: kt1,
        }));
    }
    Ok(edges)
}

/// Builds side face `k`: a full-domain extruded surface whose 4-edge wire
/// references only shared edges, with pcurves for all four boundary edges.
/// Flips the face sense for clockwise profiles so the normal points outward.
fn build_side_face(
    store: &mut TopologyStore,
    surface: &NurbsSurface,
    edges: &SharedEdges,
    k: usize,
    ccw: bool,
) -> Result<FaceId> {
    let next = (k + 1) % edges.bottom.len();
    let ((u0, u1), (v0, v1)) = surface.parameter_domain();
    let wire = store.add_wire(WireData {
        edges: vec![
            OrientedEdge::new(edges.bottom[k], true),
            OrientedEdge::new(edges.kink[next], true),
            OrientedEdge::new(edges.top[k], false),
            OrientedEdge::new(edges.kink[k], false),
        ],
        is_closed: true,
    });
    let pcurves = vec![
        FacePcurve {
            edge: edges.bottom[k],
            curve: iso_pcurve_u(u0, u1, v0)?,
        },
        FacePcurve {
            edge: edges.top[k],
            curve: iso_pcurve_u(u0, u1, v1)?,
        },
        FacePcurve {
            edge: edges.kink[k],
            curve: iso_pcurve_v_unit(u0, v0, v1)?,
        },
        FacePcurve {
            edge: edges.kink[next],
            curve: iso_pcurve_v_unit(u1, v0, v1)?,
        },
    ];
    let face = MakeNurbsFace::new(surface.clone())
        .with_boundary(wire, pcurves)
        .execute(store)?;
    if !ccw {
        store.face_mut(face)?.same_sense = false;
    }
    Ok(face)
}

/// Builds a planar cap over the given shared boundary edges (one wire around
/// all segments), oriented so its stored normal points along `outward`.
fn cap_over_edges(store: &mut TopologyStore, edges: &[EdgeId], outward: Vector3) -> Result<FaceId> {
    let wire = store.add_wire(WireData {
        edges: edges.iter().map(|&e| OrientedEdge::new(e, true)).collect(),
        is_closed: true,
    });
    cap_from_wire(store, wire, outward)
}

/// Degree-1 pcurve for a vertical kink edge: maps the edge's `[0, 1]`
/// parameter onto the `v` axis at a fixed `u` (`t → (u, v0 + t·(v1 − v0))`).
/// Exact same-parameter by construction: the extruded surface is degree-1
/// linear in `v` and the kink edge is the chord-parameterized segment
/// `p + t·direction`.
fn iso_pcurve_v_unit(u: f64, v0: f64, v1: f64) -> Result<NurbsCurve2D> {
    NurbsCurve2D::from_unweighted(
        vec![Point2::new(u, v0), Point2::new(u, v1)],
        KnotVector::new(vec![0.0, 0.0, 1.0, 1.0])?,
        1,
    )
}

/// Signed area of the closed profile projected onto the plane transverse to
/// `direction` (positive = counter-clockwise about `direction`).
///
/// # Errors
///
/// Returns a typed error when the projected area is below [`TOLERANCE`]
/// (profile plane contains the extrusion direction — orientation would be
/// ambiguous and the prism degenerate).
fn profile_signed_area(curves: &[NurbsCurve3D], direction: &Vector3) -> Result<f64> {
    let axis = direction / direction.norm();
    let (t0, _) = curves[0].parameter_domain();
    let origin = curves[0].point_at(t0)?;

    let mut poly: Vec<Vector3> = Vec::with_capacity(curves.len() * WINDING_SAMPLES);
    for curve in curves {
        let (t0, t1) = curve.parameter_domain();
        for i in 0..WINDING_SAMPLES {
            #[allow(clippy::cast_precision_loss)]
            let frac = i as f64 / WINDING_SAMPLES as f64;
            poly.push(curve.point_at(t0 + (t1 - t0) * frac)? - origin);
        }
    }
    let mut area2 = 0.0;
    for i in 0..poly.len() {
        let a = &poly[i];
        let b = &poly[(i + 1) % poly.len()];
        area2 += a.cross(b).dot(&axis);
    }
    let area = 0.5 * area2;
    if area.abs() <= TOLERANCE {
        return Err(OperationError::InvalidInput(
            "segmented prism profile encloses no area transverse to the extrusion direction".into(),
        )
        .into());
    }
    Ok(area)
}

/// Planar cap face on a pre-built shared-edge wire, oriented so its stored
/// normal points along `outward`.
fn cap_from_wire(
    store: &mut TopologyStore,
    wire: crate::topology::WireId,
    outward: Vector3,
) -> Result<FaceId> {
    let face = MakeFace::new(wire, vec![]).execute(store)?;
    let flip = match &store.face(face)?.surface {
        FaceSurface::Plane(plane) => plane.plane_normal().dot(&outward) < 0.0,
        _ => false,
    };
    if flip {
        store.face_mut(face)?.same_sense = false;
    }
    Ok(face)
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use crate::tessellation::{TessellateSolid, TessellationParams};
    use crate::topology::{FaceName, ShellData};
    use std::collections::HashMap;
    use std::f64::consts::{FRAC_PI_2, PI};

    fn p(x: f64, y: f64) -> Point3 {
        Point3::new(x, y, 0.0)
    }

    /// L-shaped 5-segment closed profile (CCW): 4 lines + 1 concave arc
    /// rounding the L's inner corner. Non-convex, so the caps exercise the
    /// CDT interior classification too.
    fn l_profile() -> Vec<ProfileSegment> {
        vec![
            ProfileSegment::Line {
                start: p(0.0, 0.0),
                end: p(4.0, 0.0),
            },
            ProfileSegment::Line {
                start: p(4.0, 0.0),
                end: p(4.0, 2.0),
            },
            ProfileSegment::Line {
                start: p(4.0, 2.0),
                end: p(3.0, 2.0),
            },
            // Concave arc (3,2) → (2,3) about (3,3), radius 1: with the −Z
            // normal the sweep π/2 → π traverses toward the L's interior.
            ProfileSegment::Arc {
                center: Point3::new(3.0, 3.0, 0.0),
                radius: 1.0,
                normal: -Vector3::z(),
                ref_dir: Vector3::x(),
                start_angle: FRAC_PI_2,
                end_angle: PI,
            },
            ProfileSegment::Line {
                start: p(2.0, 3.0),
                end: p(0.0, 0.0),
            },
        ]
    }

    const HEIGHT: f64 = 2.5;

    fn direction() -> Vector3 {
        Vector3::new(0.0, 0.0, HEIGHT)
    }

    fn tags() -> Vec<SegmentTag> {
        ["south", "east", "notch", "fillet", "hypotenuse"]
            .iter()
            .map(|t| SegmentTag::new(*t))
            .collect()
    }

    fn build_l_prism(store: &mut TopologyStore) -> SolidId {
        MakeSegmentedPrism::new(l_profile(), direction())
            .execute(store)
            .unwrap()
    }

    fn shell_of(store: &TopologyStore, solid: SolidId) -> &ShellData {
        store
            .shell(store.solid(solid).unwrap().outer_shell)
            .unwrap()
    }

    #[test]
    fn builds_one_side_face_per_segment_plus_two_caps() {
        let mut store = TopologyStore::new();
        let solid = build_l_prism(&mut store);
        let shell = shell_of(&store, solid);
        assert_eq!(shell.faces.len(), 7, "5 side faces + 2 caps");
        for &face in &shell.faces[..5] {
            assert!(
                matches!(store.face(face).unwrap().surface, FaceSurface::Nurbs(_)),
                "side faces are extruded NURBS surfaces"
            );
        }
        for &face in &shell.faces[5..] {
            assert!(
                matches!(store.face(face).unwrap().surface, FaceSurface::Plane(_)),
                "caps are planar"
            );
        }
    }

    /// Every edge of the solid is referenced by exactly two faces' wires:
    /// bottom/top segment edges by one side face + one cap, kink edges by
    /// the two adjacent side faces. 5 bottom + 5 top + 5 kink = 15 edges.
    #[test]
    fn every_shared_edge_appears_in_exactly_two_face_wires() {
        let mut store = TopologyStore::new();
        let solid = build_l_prism(&mut store);
        let shell = shell_of(&store, solid);

        let mut counts: HashMap<EdgeId, usize> = HashMap::new();
        for &face in &shell.faces {
            let face_data = store.face(face).unwrap();
            let wire = store.wire(face_data.outer_wire).unwrap();
            for oe in &wire.edges {
                *counts.entry(oe.edge).or_insert(0) += 1;
            }
            assert!(face_data.inner_wires.is_empty());
        }
        assert_eq!(counts.len(), 15, "5 bottom + 5 top + 5 kink edges");
        for (&edge, &count) in &counts {
            assert_eq!(count, 2, "edge {edge:?} referenced by {count} faces");
        }
    }

    /// Each side face records a pcurve for all 4 boundary edges, and every
    /// pcurve satisfies the same-parameter convention:
    /// `surface(pcurve(t)) == edge_curve(t)`.
    #[test]
    fn side_faces_carry_same_parameter_pcurves_for_all_boundary_edges() {
        let mut store = TopologyStore::new();
        let solid = build_l_prism(&mut store);
        let shell = shell_of(&store, solid);

        for &face in &shell.faces[..5] {
            let face_data = store.face(face).unwrap();
            let FaceSurface::Nurbs(surface) = &face_data.surface else {
                panic!("side face must be NURBS");
            };
            let wire = store.wire(face_data.outer_wire).unwrap();
            assert_eq!(wire.edges.len(), 4);
            for oe in &wire.edges {
                let pcurve = face_data
                    .pcurve_for(oe.edge)
                    .expect("every boundary edge has a pcurve");
                let edge = store.edge(oe.edge).unwrap();
                let EdgeCurve::Nurbs(edge_curve) = &edge.curve else {
                    panic!("boundary edges are NURBS curves");
                };
                for i in 0..=8 {
                    let t = edge.t_start + (edge.t_end - edge.t_start) * f64::from(i) / 8.0;
                    let uv = pcurve.point_at(t).unwrap();
                    let on_surface = surface.point_at(uv.x, uv.y).unwrap();
                    let on_edge = edge_curve.point_at(t).unwrap();
                    assert!(
                        (on_surface - on_edge).norm() < 1e-9,
                        "pcurve breaks same-parameter convention at t={t}"
                    );
                }
            }
        }
    }

    /// Position-weld watertightness (F2 pattern): after deduplicating mesh
    /// vertices by quantized position, every undirected triangle edge is
    /// used exactly twice — no boundary edges anywhere.
    #[test]
    fn tessellation_is_watertight() {
        #[allow(clippy::cast_possible_truncation)]
        fn canon_id(canon: &mut HashMap<(i64, i64, i64), u32>, p: &Point3) -> u32 {
            const Q: f64 = 1e6;
            let k = (
                (p.x * Q).round() as i64,
                (p.y * Q).round() as i64,
                (p.z * Q).round() as i64,
            );
            let next = canon.len() as u32;
            *canon.entry(k).or_insert(next)
        }

        let mut store = TopologyStore::new();
        let solid = build_l_prism(&mut store);
        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty());

        let mut canon: HashMap<(i64, i64, i64), u32> = HashMap::new();
        let mut counts: HashMap<(u32, u32), usize> = HashMap::new();
        for tri in &mesh.indices {
            let a = canon_id(&mut canon, &mesh.vertices[tri[0] as usize]);
            let b = canon_id(&mut canon, &mesh.vertices[tri[1] as usize]);
            let c = canon_id(&mut canon, &mesh.vertices[tri[2] as usize]);
            for &(x, y) in &[(a, b), (b, c), (c, a)] {
                let key = if x < y { (x, y) } else { (y, x) };
                *counts.entry(key).or_insert(0) += 1;
            }
        }
        let boundary = counts.values().filter(|&&c| c != 2).count();
        assert_eq!(
            boundary, 0,
            "segmented prism must position-weld watertight (found {boundary} boundary edges)"
        );
    }

    /// The arc segment's side face is an exact cylindrical patch: every
    /// surface sample sits at the arc radius from the vertical axis through
    /// the arc center.
    #[test]
    fn arc_side_face_is_exact_cylindrical_patch() {
        let mut store = TopologyStore::new();
        let solid = build_l_prism(&mut store);
        let shell = shell_of(&store, solid);

        // Side faces are stored in segment order; the arc is segment 3.
        let FaceSurface::Nurbs(surface) = &store.face(shell.faces[3]).unwrap().surface else {
            panic!("arc side face must be NURBS");
        };
        let ((u0, u1), (v0, v1)) = surface.parameter_domain();
        for i in 0..=16 {
            for j in 0..=4 {
                let u = u0 + (u1 - u0) * f64::from(i) / 16.0;
                let v = v0 + (v1 - v0) * f64::from(j) / 4.0;
                let sample = surface.point_at(u, v).unwrap();
                let radial = ((sample.x - 3.0).powi(2) + (sample.y - 3.0).powi(2)).sqrt();
                assert!(
                    (radial - 1.0).abs() < 1e-9,
                    "arc face sample off the cylinder: distance {radial} vs radius 1"
                );
            }
        }
    }

    /// F4-style rebuild stability: tagged names resolve across two builds
    /// into fresh stores and address identical geometry; without an op id
    /// nothing is registered.
    #[test]
    fn tagged_names_are_rebuild_stable() {
        let build = || {
            let mut store = TopologyStore::new();
            let solid = MakeSegmentedPrism::new(l_profile(), direction())
                .with_op_id(OpId::new("wall1"))
                .with_segment_tags(tags())
                .execute(&mut store)
                .unwrap();
            (store, solid)
        };
        let (store_a, _) = build();
        let (store_b, _) = build();

        for tag in tags() {
            let name = FaceName::Created {
                op: OpId::new("wall1"),
                role: FaceRole::Tagged(tag.clone()),
            };
            let fa = store_a.names().face(&name).expect("resolves in build A");
            let fb = store_b.names().face(&name).expect("resolves in build B");
            let sample = |store: &TopologyStore, f| match &store.face(f).unwrap().surface {
                FaceSurface::Nurbs(s) => s.point_at(0.3, 0.6).unwrap(),
                other => panic!("unexpected surface {other:?}"),
            };
            let pa = sample(&store_a, fa);
            let pb = sample(&store_b, fb);
            assert!(
                (pa - pb).norm() < 1e-12,
                "rebuilt face for tag {tag} moved: {pa:?} vs {pb:?}"
            );
        }
        for role in [FaceRole::CapStart, FaceRole::CapEnd] {
            let name = FaceName::Created {
                op: OpId::new("wall1"),
                role: role.clone(),
            };
            assert!(store_a.names().face(&name).is_some(), "{role:?} named");
            assert!(store_b.names().face(&name).is_some(), "{role:?} named");
        }

        // Without an op id nothing is registered.
        let mut store = TopologyStore::new();
        MakeSegmentedPrism::new(l_profile(), direction())
            .execute(&mut store)
            .unwrap();
        let name = FaceName::Created {
            op: OpId::new("wall1"),
            role: FaceRole::Tagged(SegmentTag::new("south")),
        };
        assert!(store.names().face(&name).is_none());
    }

    /// Junction stability (the D2 rationale): rebuilding with a REDUCED
    /// profile (the notch line + fillet arc replaced by one straight
    /// segment) keeps every surviving tag resolving; only the dropped
    /// segment's tag disappears. Positional Side(k) would have shifted.
    #[test]
    fn reduced_profile_keeps_surviving_tags_resolving() {
        let reduced_profile = vec![
            ProfileSegment::Line {
                start: p(0.0, 0.0),
                end: p(4.0, 0.0),
            },
            ProfileSegment::Line {
                start: p(4.0, 0.0),
                end: p(4.0, 2.0),
            },
            ProfileSegment::Line {
                start: p(4.0, 2.0),
                end: p(2.0, 3.0),
            },
            ProfileSegment::Line {
                start: p(2.0, 3.0),
                end: p(0.0, 0.0),
            },
        ];
        let reduced_tags: Vec<SegmentTag> = ["south", "east", "notch", "hypotenuse"]
            .iter()
            .map(|t| SegmentTag::new(*t))
            .collect();

        let mut store = TopologyStore::new();
        MakeSegmentedPrism::new(reduced_profile, direction())
            .with_op_id(OpId::new("wall1"))
            .with_segment_tags(reduced_tags)
            .execute(&mut store)
            .unwrap();

        for surviving in ["south", "east", "notch", "hypotenuse"] {
            let name = FaceName::Created {
                op: OpId::new("wall1"),
                role: FaceRole::Tagged(SegmentTag::new(surviving)),
            };
            assert!(
                store.names().face(&name).is_some(),
                "surviving tag {surviving} must still resolve"
            );
        }
        let dropped = FaceName::Created {
            op: OpId::new("wall1"),
            role: FaceRole::Tagged(SegmentTag::new("fillet")),
        };
        assert!(
            store.names().face(&dropped).is_none(),
            "dropped segment's tag must not resolve"
        );
    }

    /// Untagged fallback: with an op id but no tags, side faces bind
    /// positional Side(k) in segment order.
    #[test]
    fn untagged_sides_bind_positional_roles() {
        let mut store = TopologyStore::new();
        MakeSegmentedPrism::new(l_profile(), direction())
            .with_op_id(OpId::new("wall1"))
            .execute(&mut store)
            .unwrap();
        for k in 0..5 {
            let name = FaceName::Created {
                op: OpId::new("wall1"),
                role: FaceRole::Side(k),
            };
            assert!(store.names().face(&name).is_some(), "Side({k}) named");
        }
    }

    #[test]
    fn tag_count_mismatch_is_a_typed_error() {
        let mut store = TopologyStore::new();
        let result = MakeSegmentedPrism::new(l_profile(), direction())
            .with_op_id(OpId::new("wall1"))
            .with_segment_tags(vec![SegmentTag::new("only-one")])
            .execute(&mut store);
        assert!(
            matches!(
                result,
                Err(crate::error::GeolisError::Operation(
                    OperationError::InvalidInput(_)
                ))
            ),
            "tag/profile length mismatch must be a typed InvalidInput error"
        );
    }

    #[test]
    fn rejects_degenerate_inputs() {
        // Open chain: last segment does not return to the first start.
        let mut open = l_profile();
        open[4] = ProfileSegment::Line {
            start: p(2.0, 3.0),
            end: p(0.5, 0.0),
        };
        let mut store = TopologyStore::new();
        assert!(MakeSegmentedPrism::new(open, direction())
            .execute(&mut store)
            .is_err());

        // Fewer than 3 segments.
        let two = vec![
            ProfileSegment::Line {
                start: p(0.0, 0.0),
                end: p(1.0, 0.0),
            },
            ProfileSegment::Line {
                start: p(1.0, 0.0),
                end: p(0.0, 0.0),
            },
        ];
        let mut store = TopologyStore::new();
        assert!(MakeSegmentedPrism::new(two, direction())
            .execute(&mut store)
            .is_err());

        // Zero-length extrusion direction.
        let mut store = TopologyStore::new();
        assert!(
            MakeSegmentedPrism::new(l_profile(), Vector3::new(0.0, 0.0, 0.0))
                .execute(&mut store)
                .is_err()
        );

        // Profile plane contains the extrusion direction (no transverse area).
        let mut store = TopologyStore::new();
        assert!(
            MakeSegmentedPrism::new(l_profile(), Vector3::new(1.0, 0.0, 0.0))
                .execute(&mut store)
                .is_err()
        );
    }
}
