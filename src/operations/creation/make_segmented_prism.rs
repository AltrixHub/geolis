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
//! **Holes** ([`MakeSegmentedPrism::with_holes`]): each hole ring is a
//! closed segment chain of its own and receives the exact same F2
//! shared-edge treatment as the outer ring — per-segment NURBS side faces,
//! kink edges at every joint, and per-segment top / bottom ring edges
//! shared with the caps. The caps become planar faces whose
//! `inner_wires` hold one wire per hole (annulus caps), watertight by
//! construction because they reference the same `EdgeId`s as the hole
//! side faces. Ring orientation is decided per ring by the same
//! signed-area logic: an outer ring's faces point away from the enclosed
//! area, a hole ring's faces point INTO its enclosed area (the courtyard
//! side of a room wall), whichever way the caller winds the chains. Cap
//! inner wires are normalized to wind opposite the cap's outer wire (the
//! planar hole-loop convention the cap-notch rebuild classifies by).
//!
//! Naming (op id present): side face `k` binds
//! [`FaceRole::Tagged`]`(tag_k)` when the caller supplies segment tags
//! (junction-stable identity — geolis never invents it, the `OpId`
//! precedent), or positional [`FaceRole::Side`]`(k)` otherwise; the caps
//! bind `CapStart` / `CapEnd`. Hole ring faces tag via
//! [`MakeSegmentedPrism::with_hole_tags`]; the positional fallback
//! continues the `Side(k)` numbering across rings in ring order (outer
//! ring segments first, then each hole ring in the order supplied), so the
//! `FaceName` grammar is unchanged and untagged naming stays
//! deterministic.
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
    holes: Vec<Vec<ProfileSegment>>,
    direction: Vector3,
    op_id: Option<OpId>,
    segment_tags: Option<Vec<SegmentTag>>,
    hole_tags: Option<Vec<Vec<SegmentTag>>>,
}

impl MakeSegmentedPrism {
    /// Creates a segmented prism extruding `profile` along `direction`.
    #[must_use]
    pub fn new(profile: Vec<ProfileSegment>, direction: Vector3) -> Self {
        Self {
            profile,
            holes: Vec::new(),
            direction,
            op_id: None,
            segment_tags: None,
            hole_tags: None,
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

    /// Supplies hole rings: each is a closed segment chain enclosing an
    /// area inside the outer profile (validated by the same rules — at
    /// least 3 segments, endpoint-closed). Every hole ring builds its own
    /// per-segment side faces (oriented into the hole's interior), and the
    /// caps carry one inner wire per hole. Winding is free: orientation is
    /// decided per ring by its signed area.
    #[must_use]
    pub fn with_holes(mut self, holes: Vec<Vec<ProfileSegment>>) -> Self {
        self.holes = holes;
        self
    }

    /// Supplies one tag ring per hole (parallel to
    /// [`Self::with_holes`], each ring parallel to its hole's segments).
    /// Hole side faces then bind [`FaceRole::Tagged`] instead of the
    /// positional cross-ring [`FaceRole::Side`] numbering.
    #[must_use]
    pub fn with_hole_tags(mut self, tags: Vec<Vec<SegmentTag>>) -> Self {
        self.hole_tags = Some(tags);
        self
    }

    /// Builds the prism solid in the store.
    ///
    /// # Errors
    ///
    /// Returns a typed error if `direction` is zero-length, any ring (the
    /// profile or a hole) has fewer than 3 segments or is not
    /// endpoint-closed (consecutive endpoints must coincide within
    /// [`TOLERANCE`]), a tag ring does not match its segment ring, an
    /// untagged prism exceeds the positional `Side(u8)` range (counted
    /// across all rings), a ring encloses no area transverse to
    /// `direction`, or any curve / surface / face construction fails.
    pub fn execute(&self, store: &mut TopologyStore) -> Result<SolidId> {
        let rings = self.validate()?;

        // Extruded-surface normal (∂u × ∂v = tangent × direction) points
        // away from the ring's enclosed area exactly when the ring winds
        // counter-clockwise about the extrusion direction. The outer ring's
        // faces must point away from its enclosed area (outward from the
        // solid); a hole ring's faces must point INTO its enclosed area
        // (away from the material surrounding the hole). The winding is a
        // global per-ring property, so one flip decision covers every side
        // face of a ring (non-convex rings included). Computed before any
        // store mutation so a degenerate ring fails without leaving partial
        // topology behind.
        let mut areas = Vec::with_capacity(rings.len());
        for (ri, curves) in rings.iter().enumerate() {
            areas.push(ring_signed_area(curves, &self.direction, ri)?);
        }

        let mut side_faces: Vec<FaceId> = Vec::new();
        let mut bottom_rings: Vec<Vec<EdgeId>> = Vec::with_capacity(rings.len());
        let mut top_rings: Vec<Vec<EdgeId>> = Vec::with_capacity(rings.len());
        for (ri, curves) in rings.iter().enumerate() {
            let is_hole = ri > 0;
            let flip = (areas[ri] > 0.0) == is_hole;

            // One extruded surface per segment; exact for rational profiles
            // (v = 0 row is the segment curve, v = 1 its translate).
            let mut surfaces = Vec::with_capacity(curves.len());
            for curve in curves {
                surfaces.push(NurbsSurface::extrude(curve, self.direction)?);
            }

            let edges = build_shared_edges(store, curves, &surfaces, self.direction)?;

            // Side faces: one per segment, each a full-domain extruded
            // surface whose 4-edge wire references only shared edges, with
            // pcurves for all four boundary edges.
            for (k, surface) in surfaces.iter().enumerate() {
                side_faces.push(build_side_face(store, surface, &edges, k, flip)?);
            }
            bottom_rings.push(edges.bottom);
            top_rings.push(edges.top);
        }

        // Caps: one planar face per end whose outer wire runs around ALL of
        // the outer ring's per-segment boundary edges and whose inner wires
        // hold one hole ring each (the same edges the side faces use). A
        // hole wire is reversed when its winding matches the outer ring's,
        // so cap inner wires always wind opposite the cap's outer wire.
        let outer_ccw = areas[0] > 0.0;
        let reverse_holes: Vec<bool> = areas[1..]
            .iter()
            .map(|&area| (area > 0.0) == outer_ccw)
            .collect();
        let bottom_cap = cap_over_rings(store, &bottom_rings, &reverse_holes, -self.direction)?;
        let top_cap = cap_over_rings(store, &top_rings, &reverse_holes, self.direction)?;

        self.bind_names(store, &rings, &side_faces, bottom_cap, top_cap);

        let mut faces = side_faces;
        faces.push(bottom_cap);
        faces.push(top_cap);
        Ok(finish_solid(store, faces))
    }

    /// Binds the persistent face names when an op id is present: per-ring
    /// tags where supplied, positional `Side(k)` continuing across rings
    /// otherwise, and `CapStart` / `CapEnd` for the caps.
    fn bind_names(
        &self,
        store: &mut TopologyStore,
        rings: &[Vec<NurbsCurve3D>],
        side_faces: &[FaceId],
        bottom_cap: FaceId,
        top_cap: FaceId,
    ) {
        let Some(op) = &self.op_id else {
            return;
        };
        let mut global = 0usize;
        for (ri, curves) in rings.iter().enumerate() {
            let ring_tags: Option<&Vec<SegmentTag>> = if ri == 0 {
                self.segment_tags.as_ref()
            } else {
                self.hole_tags.as_ref().map(|rings| &rings[ri - 1])
            };
            for k in 0..curves.len() {
                let role = match ring_tags {
                    Some(tags) => FaceRole::Tagged(tags[k].clone()),
                    // Range validated up front in `validate`.
                    #[allow(clippy::cast_possible_truncation)]
                    None => FaceRole::Side(global as u8),
                };
                bind_created_face(store, side_faces[global], op, role);
                global += 1;
            }
        }
        bind_created_face(store, bottom_cap, op, FaceRole::CapStart);
        bind_created_face(store, top_cap, op, FaceRole::CapEnd);
    }

    /// Validates the inputs and builds the per-segment curves of every
    /// ring (`rings[0]` is the outer profile, `rings[1..]` the holes).
    fn validate(&self) -> Result<Vec<Vec<NurbsCurve3D>>> {
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
        for (hi, hole) in self.holes.iter().enumerate() {
            let hn = hole.len();
            if hn < 3 {
                return Err(OperationError::InvalidInput(format!(
                    "segmented prism hole {hi} needs at least 3 segments, got {hn}"
                ))
                .into());
            }
        }
        if let Some(tags) = &self.segment_tags {
            if tags.len() != n {
                return Err(OperationError::InvalidInput(format!(
                    "segment tag count {} does not match profile segment count {n}",
                    tags.len()
                ))
                .into());
            }
        }
        if let Some(hole_tags) = &self.hole_tags {
            if hole_tags.len() != self.holes.len() {
                return Err(OperationError::InvalidInput(format!(
                    "hole tag ring count {} does not match hole count {}",
                    hole_tags.len(),
                    self.holes.len()
                ))
                .into());
            }
            for (hi, (tags, hole)) in hole_tags.iter().zip(&self.holes).enumerate() {
                if tags.len() != hole.len() {
                    return Err(OperationError::InvalidInput(format!(
                        "hole {hi} tag count {} does not match its segment count {}",
                        tags.len(),
                        hole.len()
                    ))
                    .into());
                }
            }
        }
        // Positional Side(u8) numbering continues across rings, so ANY
        // untagged ring requires the cross-ring total to stay in range.
        let total = n + self.holes.iter().map(Vec::len).sum::<usize>();
        let any_untagged =
            self.segment_tags.is_none() || (!self.holes.is_empty() && self.hole_tags.is_none());
        if self.op_id.is_some() && any_untagged && total > usize::from(u8::MAX) + 1 {
            return Err(OperationError::InvalidInput(format!(
                "untagged segmented prism supports at most 256 named segments \
                 (positional Side(u8)); got {total} across all rings — supply \
                 segment tags"
            ))
            .into());
        }

        let mut rings = Vec::with_capacity(1 + self.holes.len());
        for (ri, ring) in std::iter::once(&self.profile)
            .chain(self.holes.iter())
            .enumerate()
        {
            let curves: Vec<NurbsCurve3D> = ring
                .iter()
                .map(ProfileSegment::curve)
                .collect::<Result<_>>()?;

            // The chain must be endpoint-closed: each segment's end
            // coincides with the next segment's start (cyclically) within
            // TOLERANCE.
            let rn = curves.len();
            for k in 0..rn {
                let next = (k + 1) % rn;
                let (_, t1) = curves[k].parameter_domain();
                let end = curves[k].point_at(t1)?;
                let (t0, _) = curves[next].parameter_domain();
                let start = curves[next].point_at(t0)?;
                if (end - start).norm() > TOLERANCE {
                    return Err(OperationError::InvalidInput(format!(
                        "{} segment chain is not closed between segment {k} and {next}",
                        ring_label(ri)
                    ))
                    .into());
                }
            }
            rings.push(curves);
        }
        Ok(rings)
    }
}

/// Human-readable ring identifier for error messages (`ri == 0` is the
/// outer profile; holes are numbered from 0).
fn ring_label(ri: usize) -> String {
    if ri == 0 {
        "profile".to_owned()
    } else {
        format!("hole {}", ri - 1)
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
/// `flip` reverses the face sense when the ring winds against its required
/// orientation (outer rings point away from the enclosed area, hole rings
/// into it).
fn build_side_face(
    store: &mut TopologyStore,
    surface: &NurbsSurface,
    edges: &SharedEdges,
    k: usize,
    flip: bool,
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
    if flip {
        store.face_mut(face)?.same_sense = false;
    }
    Ok(face)
}

/// Builds a planar cap over the given per-ring shared boundary edges
/// (`rings[0]` is the outer ring, `rings[1..]` the holes), oriented so its
/// stored normal points along `outward`. Hole wires are reversed when
/// `reverse_holes` says their winding matches the outer ring's, so every
/// cap inner wire winds opposite the cap's outer wire.
fn cap_over_rings(
    store: &mut TopologyStore,
    rings: &[Vec<EdgeId>],
    reverse_holes: &[bool],
    outward: Vector3,
) -> Result<FaceId> {
    let ring_wire = |store: &mut TopologyStore, edges: &[EdgeId], reverse: bool| {
        let oriented: Vec<OrientedEdge> = if reverse {
            edges
                .iter()
                .rev()
                .map(|&e| OrientedEdge::new(e, false))
                .collect()
        } else {
            edges.iter().map(|&e| OrientedEdge::new(e, true)).collect()
        };
        store.add_wire(WireData {
            edges: oriented,
            is_closed: true,
        })
    };
    let outer = ring_wire(store, &rings[0], false);
    let inners: Vec<crate::topology::WireId> = rings[1..]
        .iter()
        .zip(reverse_holes)
        .map(|(edges, &reverse)| ring_wire(store, edges, reverse))
        .collect();
    cap_from_wires(store, outer, inners, outward)
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

/// Signed area of a closed ring projected onto the plane transverse to
/// `direction` (positive = counter-clockwise about `direction`). `ri` is
/// the ring index for error labelling (0 = outer profile).
///
/// # Errors
///
/// Returns a typed error when the projected area is below [`TOLERANCE`]
/// (ring plane contains the extrusion direction — orientation would be
/// ambiguous and the prism degenerate).
fn ring_signed_area(curves: &[NurbsCurve3D], direction: &Vector3, ri: usize) -> Result<f64> {
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
        return Err(OperationError::InvalidInput(format!(
            "segmented prism {} encloses no area transverse to the extrusion direction",
            ring_label(ri)
        ))
        .into());
    }
    Ok(area)
}

/// Planar cap face on pre-built shared-edge wires, oriented so its stored
/// normal points along `outward`.
fn cap_from_wires(
    store: &mut TopologyStore,
    outer: crate::topology::WireId,
    inners: Vec<crate::topology::WireId>,
    outward: Vector3,
) -> Result<FaceId> {
    let face = MakeFace::new(outer, inners).execute(store)?;
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

    /// Asserts each given side face records a pcurve for all 4 boundary
    /// edges, and every pcurve satisfies the same-parameter convention:
    /// `surface(pcurve(t)) == edge_curve(t)`.
    fn assert_same_parameter_pcurves(store: &TopologyStore, faces: &[FaceId]) {
        for &face in faces {
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

    /// Each side face records a pcurve for all 4 boundary edges, and every
    /// pcurve satisfies the same-parameter convention.
    #[test]
    fn side_faces_carry_same_parameter_pcurves_for_all_boundary_edges() {
        let mut store = TopologyStore::new();
        let solid = build_l_prism(&mut store);
        let shell = shell_of(&store, solid);
        assert_same_parameter_pcurves(&store, &shell.faces[..5]);
    }

    /// Position-weld watertightness helper (F2 pattern): after
    /// deduplicating mesh vertices by quantized position, every undirected
    /// triangle edge must be used exactly twice — no boundary edges
    /// anywhere.
    fn assert_position_weld_watertight(store: &TopologyStore, solid: SolidId) {
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

        let mesh = TessellateSolid::new(solid, TessellationParams::default())
            .execute(store)
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

    #[test]
    fn tessellation_is_watertight() {
        let mut store = TopologyStore::new();
        let solid = build_l_prism(&mut store);
        assert_position_weld_watertight(&store, solid);
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

    // ── Hole (annulus) support ──────────────────────────────────────────

    /// Axis-aligned rectangle chain. `ccw = true` winds counter-clockwise
    /// about `+Z`, `false` clockwise.
    fn rect_ring(x0: f64, y0: f64, x1: f64, y1: f64, ccw: bool) -> Vec<ProfileSegment> {
        let corners = if ccw {
            [p(x0, y0), p(x1, y0), p(x1, y1), p(x0, y1)]
        } else {
            [p(x0, y0), p(x0, y1), p(x1, y1), p(x1, y0)]
        };
        (0..4)
            .map(|i| ProfileSegment::Line {
                start: corners[i],
                end: corners[(i + 1) % 4],
            })
            .collect()
    }

    /// A closed square room: outer 6 × 6 footprint with a 5 × 5 courtyard
    /// hole (wall thickness 0.5). `hole_ccw` picks the hole winding — the
    /// builder auto-orients either way.
    fn room_prism(hole_ccw: bool) -> MakeSegmentedPrism {
        MakeSegmentedPrism::new(rect_ring(0.0, 0.0, 6.0, 6.0, true), direction())
            .with_holes(vec![rect_ring(0.5, 0.5, 5.5, 5.5, hole_ccw)])
    }

    fn outer_room_tags() -> Vec<SegmentTag> {
        ["outer-s", "outer-e", "outer-n", "outer-w"]
            .iter()
            .map(|t| SegmentTag::new(*t))
            .collect()
    }

    fn hole_room_tags() -> Vec<Vec<SegmentTag>> {
        vec![["inner-w", "inner-n", "inner-e", "inner-s"]
            .iter()
            .map(|t| SegmentTag::new(*t))
            .collect()]
    }

    #[test]
    fn annulus_builds_hole_faces_and_cap_inner_wires() {
        let mut store = TopologyStore::new();
        let solid = room_prism(false).execute(&mut store).unwrap();
        let shell = shell_of(&store, solid);
        assert_eq!(shell.faces.len(), 10, "4 outer + 4 hole sides + 2 caps");
        for &face in &shell.faces[..8] {
            assert!(
                matches!(store.face(face).unwrap().surface, FaceSurface::Nurbs(_)),
                "all side faces (both rings) are extruded NURBS surfaces"
            );
        }
        for &face in &shell.faces[8..] {
            let face_data = store.face(face).unwrap();
            assert!(
                matches!(face_data.surface, FaceSurface::Plane(_)),
                "caps are planar"
            );
            assert_eq!(
                face_data.inner_wires.len(),
                1,
                "each cap carries one inner wire per hole"
            );
            let inner = store.wire(face_data.inner_wires[0]).unwrap();
            assert_eq!(inner.edges.len(), 4, "hole ring has 4 boundary edges");
        }
    }

    /// Every edge of the annulus solid is referenced by exactly two faces'
    /// wires (inner cap wires included): 2 rings × (4 bottom + 4 top +
    /// 4 kink) = 24 edges.
    #[test]
    fn annulus_shared_edges_appear_in_exactly_two_face_wires() {
        let mut store = TopologyStore::new();
        let solid = room_prism(false).execute(&mut store).unwrap();
        let shell = shell_of(&store, solid);

        let mut counts: HashMap<EdgeId, usize> = HashMap::new();
        for &face in &shell.faces {
            let face_data = store.face(face).unwrap();
            let wires =
                std::iter::once(face_data.outer_wire).chain(face_data.inner_wires.iter().copied());
            for wire in wires {
                for oe in &store.wire(wire).unwrap().edges {
                    *counts.entry(oe.edge).or_insert(0) += 1;
                }
            }
        }
        assert_eq!(counts.len(), 24, "2 rings x (4 bottom + 4 top + 4 kink)");
        for (&edge, &count) in &counts {
            assert_eq!(count, 2, "edge {edge:?} referenced by {count} faces");
        }
    }

    #[test]
    fn annulus_tessellation_is_watertight_for_both_hole_windings() {
        for hole_ccw in [false, true] {
            let mut store = TopologyStore::new();
            let solid = room_prism(hole_ccw).execute(&mut store).unwrap();
            assert_position_weld_watertight(&store, solid);
        }
    }

    #[test]
    fn annulus_side_faces_carry_same_parameter_pcurves() {
        let mut store = TopologyStore::new();
        let solid = room_prism(false).execute(&mut store).unwrap();
        let shell = shell_of(&store, solid);
        assert_same_parameter_pcurves(&store, &shell.faces[..8]);
    }

    /// Hole side faces point INTO the courtyard (away from the wall
    /// material), whichever way the caller winds the hole chain.
    #[test]
    fn annulus_hole_faces_point_into_courtyard() {
        for hole_ccw in [false, true] {
            let mut store = TopologyStore::new();
            let solid = room_prism(hole_ccw).execute(&mut store).unwrap();
            let shell = shell_of(&store, solid);
            let courtyard = Point3::new(3.0, 3.0, HEIGHT / 2.0);

            for &face in &shell.faces[4..8] {
                let face_data = store.face(face).unwrap();
                let FaceSurface::Nurbs(surface) = &face_data.surface else {
                    panic!("hole side face must be NURBS");
                };
                let ((u0, u1), (v0, v1)) = surface.parameter_domain();
                let (um, vm) = (0.5 * (u0 + u1), 0.5 * (v0 + v1));
                let eps = 1e-6;
                let p0 = surface.point_at(um, vm).unwrap();
                let du = surface.point_at(um + eps, vm).unwrap() - p0;
                let dv = surface.point_at(um, vm + eps).unwrap() - p0;
                let mut normal = du.cross(&dv);
                if !face_data.same_sense {
                    normal = -normal;
                }
                assert!(
                    normal.dot(&(courtyard - p0)) > 0.0,
                    "hole face normal must point into the courtyard \
                     (hole_ccw = {hole_ccw}, sample {p0:?})"
                );
            }
        }
    }

    /// Cap inner wires are normalized to wind opposite the cap's outer
    /// wire, whichever way the caller wound the hole chain.
    #[test]
    fn annulus_cap_inner_wires_wind_opposite_the_outer_wire() {
        /// Shoelace area over a wire's oriented traversal points in the
        /// XY plane (the caps are horizontal here).
        fn wire_area_xy(store: &TopologyStore, wire: crate::topology::WireId) -> f64 {
            let wire = store.wire(wire).unwrap();
            let mut poly: Vec<Point3> = Vec::new();
            for oe in &wire.edges {
                let edge = store.edge(oe.edge).unwrap();
                let v = if oe.forward { edge.start } else { edge.end };
                poly.push(store.vertex(v).unwrap().point);
            }
            let mut area2 = 0.0;
            for i in 0..poly.len() {
                let a = &poly[i];
                let b = &poly[(i + 1) % poly.len()];
                area2 += a.x * b.y - b.x * a.y;
            }
            0.5 * area2
        }

        for hole_ccw in [false, true] {
            let mut store = TopologyStore::new();
            let solid = room_prism(hole_ccw).execute(&mut store).unwrap();
            let shell = shell_of(&store, solid);
            for &cap in &shell.faces[8..] {
                let face_data = store.face(cap).unwrap();
                let outer = wire_area_xy(&store, face_data.outer_wire);
                let inner = wire_area_xy(&store, face_data.inner_wires[0]);
                assert!(
                    outer * inner < 0.0,
                    "cap inner wire must wind opposite the outer wire \
                     (hole_ccw = {hole_ccw}, outer {outer}, inner {inner})"
                );
            }
        }
    }

    /// Hole tags bind hole ring faces, resolve across rebuilds into
    /// identical geometry, and the outer ring's tags stay untouched.
    #[test]
    fn annulus_hole_tags_are_rebuild_stable() {
        let build = || {
            let mut store = TopologyStore::new();
            room_prism(false)
                .with_op_id(OpId::new("room1"))
                .with_segment_tags(outer_room_tags())
                .with_hole_tags(hole_room_tags())
                .execute(&mut store)
                .unwrap();
            store
        };
        let store_a = build();
        let store_b = build();

        let all_tags: Vec<SegmentTag> = outer_room_tags()
            .into_iter()
            .chain(hole_room_tags().remove(0))
            .collect();
        for tag in all_tags {
            let name = FaceName::Created {
                op: OpId::new("room1"),
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
    }

    /// Positional fallback: `Side(k)` numbering continues across rings in
    /// ring order — outer segments 0..3, hole segments 4..7.
    #[test]
    fn annulus_untagged_continues_side_numbering_across_rings() {
        let mut store = TopologyStore::new();
        let solid = room_prism(false)
            .with_op_id(OpId::new("room1"))
            .execute(&mut store)
            .unwrap();
        let shell = shell_of(&store, solid);
        for k in 0..8u8 {
            let name = FaceName::Created {
                op: OpId::new("room1"),
                role: FaceRole::Side(k),
            };
            let face = store
                .names()
                .face(&name)
                .unwrap_or_else(|| panic!("Side({k}) must resolve"));
            assert_eq!(
                face,
                shell.faces[usize::from(k)],
                "Side({k}) must bind ring-order face {k}"
            );
        }
    }

    #[test]
    fn annulus_rejects_invalid_hole_inputs() {
        // A hole with fewer than 3 segments.
        let two_segment_hole = vec![
            ProfileSegment::Line {
                start: p(1.0, 1.0),
                end: p(2.0, 1.0),
            },
            ProfileSegment::Line {
                start: p(2.0, 1.0),
                end: p(1.0, 1.0),
            },
        ];
        let mut store = TopologyStore::new();
        assert!(
            MakeSegmentedPrism::new(rect_ring(0.0, 0.0, 6.0, 6.0, true), direction())
                .with_holes(vec![two_segment_hole])
                .execute(&mut store)
                .is_err()
        );

        // An open hole chain.
        let mut open_hole = rect_ring(0.5, 0.5, 5.5, 5.5, false);
        open_hole[3] = ProfileSegment::Line {
            start: p(5.5, 0.5),
            end: p(1.0, 0.5),
        };
        let mut store = TopologyStore::new();
        assert!(
            MakeSegmentedPrism::new(rect_ring(0.0, 0.0, 6.0, 6.0, true), direction())
                .with_holes(vec![open_hole])
                .execute(&mut store)
                .is_err()
        );

        // Hole tag ring count mismatch (1 hole, 2 tag rings).
        let mut store = TopologyStore::new();
        assert!(matches!(
            room_prism(false)
                .with_op_id(OpId::new("room1"))
                .with_segment_tags(outer_room_tags())
                .with_hole_tags(vec![
                    hole_room_tags().remove(0),
                    vec![SegmentTag::new("extra")]
                ])
                .execute(&mut store),
            Err(crate::error::GeolisError::Operation(
                OperationError::InvalidInput(_)
            ))
        ));

        // Hole tag length mismatch inside a ring.
        let mut store = TopologyStore::new();
        assert!(matches!(
            room_prism(false)
                .with_op_id(OpId::new("room1"))
                .with_segment_tags(outer_room_tags())
                .with_hole_tags(vec![vec![SegmentTag::new("only-one")]])
                .execute(&mut store),
            Err(crate::error::GeolisError::Operation(
                OperationError::InvalidInput(_)
            ))
        ));
    }
}
