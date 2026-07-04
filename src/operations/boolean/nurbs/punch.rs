//! Punching SSI loops as trim holes on target NURBS faces.
//!
//! Each cut loop's `uv_a` trace (in the target face's parameter space) becomes a
//! clockwise trim hole on that face, and its 3D `points` become a closed inner
//! wire so the 3D topology stays consistent with the UV trim.
//!
//! ## Documented approximations (v1)
//!
//! - The hole pcurve is a **degree-1 polyline** through the SSI trace points
//!   (not a fitted rational arc); the trimmed CDT tessellator samples loops into
//!   polylines anyway, so this is lossless for tessellation. Chained rings
//!   ([`punch_chain`]) also build their 3D ring edges as degree-1 polylines
//!   through the same samples, with knot-compatible target/tool pcurves —
//!   the edge sample cache then serves ONE canonical 3D polyline per rim
//!   edge to every face that references it (F6 R3 rim weld).
//! - Loops on a geometrically closed tool arrive genuinely `closed` from the
//!   SSI marcher (periodic-domain wrapping) with exact seam samples at every
//!   crossing, shared by the punch (`uv_a`) and band (`uv_b`) rings. The ring's
//!   wrap-around segment therefore closes a genuinely small sub-step arc. If a
//!   seam sample did not converge in the marcher, the wrap segment degrades to
//!   a single straight chord across the sub-step gap — bounded and documented,
//!   never fabricated geometry.

use crate::error::{OperationError, Result};
use crate::geometry::nurbs::{KnotVector, NurbsCurve2D, NurbsCurve3D, NurbsSurface};
use crate::math::{Point2, Point3, TOLERANCE};
use crate::topology::{
    EdgeCurve, EdgeData, FaceId, FacePcurve, FaceSurface, FaceTrim, OrientedEdge, TopologyStore,
    TrimLoop, VertexData, WireData, WireId,
};

use super::loops::CutLoop;
use super::stitch::CutChain;

/// The punched hole ring of a chained (multi-face-tool) loop: the closed ring
/// wire plus its per-segment edges and junction vertices, all in chain order.
/// `edges[i]` runs chain segment `i`; `junctions[i]` is the shared vertex at
/// segment `i`'s START (the kink crossing between segments `i - 1` and `i`).
#[derive(Debug, Clone)]
pub(crate) struct ChainRing {
    /// The closed hole ring wire (one edge per chain segment).
    pub wire: WireId,
    /// Per-segment ring edges, in chain order.
    pub edges: Vec<crate::topology::EdgeId>,
    /// Per-junction shared vertices, in chain order.
    pub junctions: Vec<crate::topology::VertexId>,
    /// Tool-UV pcurve of `edges[i]` (in segment `i`'s tool face parameter
    /// space), knot-compatible with the ring edge — the band fragments
    /// register these so edge-driven tessellation pins their rim vertices
    /// to the same canonical edge samples as the punched target face.
    pub tool_pcurves: Vec<NurbsCurve2D>,
}

/// Punches a chained loop onto its target face: the concatenated target-UV
/// traces become one CW trim hole, and the 3D trace becomes a closed inner
/// wire of one degree-1 polyline edge PER chain segment (through the exact
/// SSI samples) with shared junction vertices at the kink crossings — so the
/// per-tool-face band fragments can each reference exactly their segment's
/// ring edge (F2 shared-edge topology).
///
/// Each ring edge carries knot-compatible pcurves on BOTH sides: the
/// target-UV pcurve is registered on the punched face here, the tool-UV
/// pcurve is returned in the [`ChainRing`] for the band fragments. Both
/// faces therefore take the edge-driven tessellation path and pin their rim
/// vertices to the ONE canonical per-edge 3D polyline (exact cross-face
/// coincidence — the F6 R3 rim weld).
///
/// # Errors
///
/// Returns an error if the target face is not a NURBS face, the concatenated
/// trace degenerates, or curve / wire construction fails.
pub(crate) fn punch_chain(store: &mut TopologyStore, chain: &CutChain) -> Result<ChainRing> {
    let target_face = chain.single_target_face().ok_or_else(|| {
        OperationError::Failed(
            "punch_chain requires a chained loop on a single target face \
             (target-crossing chains go through the face splitter)"
                .into(),
        )
    })?;
    let surface = nurbs_surface_of(store, target_face)?;

    // 1. CW hole pcurve from the concatenated target-UV traces (junction
    //    samples are welded, so consecutive duplicates collapse in dedup).
    let mut uv: Vec<Point2> = Vec::new();
    for seg in &chain.segments {
        uv.extend_from_slice(&seg.branch.uv_a);
    }
    let hole_loop = ssi_trim_loop(&uv, true)?;

    // 2. Shared junction vertices + one degree-1 trace edge per segment with
    //    lockstep target/tool pcurves (the shared trace-edge construction).
    let n = chain.segments.len();
    let mut junctions = Vec::with_capacity(n);
    for seg in &chain.segments {
        let start = *seg
            .branch
            .points
            .first()
            .ok_or_else(|| OperationError::Failed("empty chain segment trace".into()))?;
        junctions.push(store.add_vertex(VertexData::new(start)));
    }
    let mut edges = Vec::with_capacity(n);
    let mut target_pcurves = Vec::with_capacity(n);
    let mut tool_pcurves = Vec::with_capacity(n);
    for (i, seg) in chain.segments.iter().enumerate() {
        let (edge, pcurve, tool_pcurve) = super::split::trace_edge(
            store,
            &seg.branch.points,
            &seg.branch.uv_a,
            &seg.branch.uv_b,
            junctions[i],
            junctions[(i + 1) % n],
        )?;
        edges.push(edge);
        target_pcurves.push(pcurve);
        tool_pcurves.push(tool_pcurve);
    }
    let wire = store.add_wire(WireData {
        edges: edges.iter().map(|&e| OrientedEdge::new(e, true)).collect(),
        is_closed: true,
    });

    // 3. Attach to the face: full-domain outer trim if absent, push the hole,
    //    append the inner wire, and register the target-UV ring pcurves so
    //    the hole tessellates edge-driven.
    let face = store.face_mut(target_face)?;
    let trim = face
        .trim
        .get_or_insert_with(|| FaceTrim::new(full_domain_outer_loop(&surface), Vec::new()));
    trim.holes.push(hole_loop);
    face.inner_wires.push(wire);
    for (&edge, curve) in edges.iter().zip(target_pcurves) {
        face.pcurves.push(FacePcurve { edge, curve });
    }

    Ok(ChainRing {
        wire,
        edges,
        junctions,
        tool_pcurves,
    })
}

/// Punches a single cut loop onto its target face: adds a CW trim hole and a 3D
/// inner wire. Returns the [`WireId`] of the hole ring wire it created so callers
/// can share the same wire with the band (hole-wall) face.
///
/// # Errors
///
/// Returns an error if the target face is not a NURBS face, the loop degenerates
/// to fewer than 3 distinct UV points, or curve / wire construction fails.
pub(crate) fn punch_loop(store: &mut TopologyStore, cut: &CutLoop) -> Result<WireId> {
    let surface = nurbs_surface_of(store, cut.target_face)?;

    // 1. Build the CW hole pcurve from the target-UV trace.
    let hole_loop = hole_loop_from_trace(&cut.branch.uv_a)?;

    // 2. Build the 3D inner wire from the SSI 3D points.
    let inner_wire = build_ring_wire(store, &cut.branch.points)?;

    // 3. Attach to the face: create a full-domain outer trim if absent, push the
    //    hole, and append the inner wire.
    let face = store.face_mut(cut.target_face)?;
    let trim = face
        .trim
        .get_or_insert_with(|| FaceTrim::new(full_domain_outer_loop(&surface), Vec::new()));
    trim.holes.push(hole_loop);
    face.inner_wires.push(inner_wire);
    Ok(inner_wire)
}

/// Returns a clone of the NURBS surface backing `face`, or an error if `face` is
/// not a NURBS face (planar target splitting is out of scope).
fn nurbs_surface_of(store: &TopologyStore, face: FaceId) -> Result<NurbsSurface> {
    match &store.face(face)?.surface {
        FaceSurface::Nurbs(s) => Ok(s.clone()),
        _ => Err(OperationError::Failed(
            "through-cut subtract can only punch NURBS target faces".into(),
        )
        .into()),
    }
}

/// Builds the full-domain outer trim loop (CCW rectangle over the surface's
/// parameter domain) for a NURBS surface.
///
/// Shared with the trimming pipeline: when a target face has no prior trim, the
/// hole must sit inside an explicit outer boundary, which for an untrimmed face
/// is the whole parameter rectangle.
pub(crate) fn full_domain_outer_loop(surface: &NurbsSurface) -> TrimLoop {
    let ((u0, u1), (v0, v1)) = surface.parameter_domain();
    let corners = [(u0, v0), (u1, v0), (u1, v1), (u0, v1)];
    let mut curves = Vec::with_capacity(4);
    for i in 0..4 {
        let a = corners[i];
        let b = corners[(i + 1) % 4];
        curves.push(uv_segment(Point2::new(a.0, a.1), Point2::new(b.0, b.1)));
    }
    TrimLoop::new(curves)
}

/// A degree-1 two-point UV line segment.
fn uv_segment(a: Point2, b: Point2) -> NurbsCurve2D {
    // Construction is infallible for two distinct points with a valid clamped
    // knot vector; the domain-loop builder only ever passes distinct corners.
    NurbsCurve2D::from_unweighted(
        vec![a, b],
        KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap_or_else(|_| unreachable!()),
        1,
    )
    .unwrap_or_else(|_| unreachable!())
}

/// Converts a target-UV SSI trace into a clockwise (hole) `TrimLoop` of degree-1
/// segments. The marcher emits genuinely closed traces, so the wrap-around
/// segment closes a small sub-step arc (it only degrades to a straight chord if
/// a marcher seam sample did not converge).
fn hole_loop_from_trace(uv_a: &[Point2]) -> Result<TrimLoop> {
    ssi_trim_loop(uv_a, true)
}

/// Converts an SSI UV trace into a closed `TrimLoop` of degree-1 segments with
/// the requested winding (`clockwise` for a subtract hole, counter-clockwise for
/// an intersect keep-inside outer boundary). Shared by the subtract and intersect
/// paths so the ring geometry is identical in both modes.
pub(crate) fn ssi_trim_loop(uv: &[Point2], clockwise: bool) -> Result<TrimLoop> {
    let mut pts = dedup_uv(uv);
    if pts.len() < 3 {
        return Err(OperationError::Failed(
            "SSI loop degenerated to fewer than 3 distinct UV points".into(),
        )
        .into());
    }
    // Reverse only when the current winding disagrees with the requested one.
    let is_clockwise = signed_area(&pts) < 0.0;
    if is_clockwise != clockwise {
        pts.reverse();
    }

    let n = pts.len();
    let mut curves = Vec::with_capacity(n);
    for i in 0..n {
        let a = pts[i];
        let b = pts[(i + 1) % n];
        // The wrap-around segment (i = n-1) closes the ring; with the marcher's
        // exact seam samples it is a short sub-step arc, not a seam chord.
        curves.push(uv_segment(a, b));
    }
    Ok(TrimLoop::new(curves))
}

/// Removes consecutive near-duplicate UV points and a coincident wrap point.
fn dedup_uv(pts: &[Point2]) -> Vec<Point2> {
    let mut out: Vec<Point2> = Vec::with_capacity(pts.len());
    for &p in pts {
        if out.last().is_none_or(|q| (p - q).norm() > 1e-9) {
            out.push(p);
        }
    }
    while out.len() >= 2 && (out[0] - out[out.len() - 1]).norm() < 1e-9 {
        out.pop();
    }
    out
}

/// Shoelace signed area of a UV polygon. Positive = counter-clockwise.
fn signed_area(pts: &[Point2]) -> f64 {
    let n = pts.len();
    let mut a2 = 0.0;
    for i in 0..n {
        let p = pts[i];
        let q = pts[(i + 1) % n];
        a2 += p.x * q.y - q.x * p.y;
    }
    0.5 * a2
}

/// Builds a closed 3D ring wire from the SSI 3D trace via a single degree-3
/// interpolated NURBS edge (closing the seam gap). Shared by the subtract (hole
/// inner wire) and intersect (kept-disc outer wire) paths.
pub(crate) fn build_ring_wire(store: &mut TopologyStore, points: &[Point3]) -> Result<WireId> {
    let mut pts = dedup3d(points);
    if pts.len() < 4 {
        return Err(OperationError::Failed(
            "inner hole wire needs at least 4 distinct 3D points".into(),
        )
        .into());
    }
    // Close the ring explicitly so the interpolated curve returns to its start.
    if (pts[0] - pts[pts.len() - 1]).norm() > TOLERANCE {
        pts.push(pts[0]);
    }
    let degree = 3.min(pts.len() - 1);
    let (curve, _) = NurbsCurve3D::interpolate(&pts, degree)?;

    let start = *pts.first().unwrap_or_else(|| unreachable!());
    let end = *pts.last().unwrap_or_else(|| unreachable!());
    let v_start = store.add_vertex(VertexData::new(start));
    // Start and end coincide (closed ring); reuse the start vertex.
    let v_end = if (end - start).norm() < TOLERANCE {
        v_start
    } else {
        store.add_vertex(VertexData::new(end))
    };
    let (t0, t1) = curve.parameter_domain();
    let edge = store.add_edge(EdgeData {
        start: v_start,
        end: v_end,
        curve: EdgeCurve::Nurbs(curve),
        t_start: t0,
        t_end: t1,
    });
    Ok(store.add_wire(WireData {
        edges: vec![OrientedEdge::new(edge, true)],
        is_closed: true,
    }))
}

/// Removes consecutive near-duplicate 3D points.
fn dedup3d(pts: &[Point3]) -> Vec<Point3> {
    let mut out: Vec<Point3> = Vec::with_capacity(pts.len());
    for &p in pts {
        if out.last().is_none_or(|q| (p - q).norm() > TOLERANCE) {
            out.push(p);
        }
    }
    out
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::cast_precision_loss)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use crate::operations::boolean::nurbs::loops::{
        collect_nurbs_faces, extract_cut_loops, ToolFaceCut,
    };
    use crate::operations::creation::{MakeCurvedSlab, MakeNurbsTube};
    use crate::tessellation::{TessellateFace, TessellationParams};
    use crate::topology::SolidId;

    fn solid_faces(store: &TopologyStore, solid: SolidId) -> Vec<FaceId> {
        let shell = store
            .shell(store.solid(solid).unwrap().outer_shell)
            .unwrap();
        shell.faces.clone()
    }

    /// Punches the slab×tube front-face loop and returns (store, that face id).
    fn punch_front(radius: f64) -> (TopologyStore, FaceId, Vec<Point2>) {
        let mut store = TopologyStore::new();
        let slab = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut store)
            .unwrap();
        let tube = MakeNurbsTube::new(Point3::new(3.0, 3.0, -1.5), radius, 5.0)
            .execute(&mut store)
            .unwrap();
        let target = collect_nurbs_faces(&store, &solid_faces(&store, slab));
        let tool = collect_nurbs_faces(&store, &solid_faces(&store, tube));
        let cuts = extract_cut_loops(&target, &tool).unwrap();
        // Punch both loops; return the higher-z (front) face for inspection.
        let mut front_face = None;
        let mut front_trace = Vec::new();
        for loop_pair in &cuts {
            let ToolFaceCut::Through { loops, .. } = loop_pair else {
                panic!("expected a through cut");
            };
            for cut in loops {
                punch_loop(&mut store, cut).unwrap();
            }
            // Identify the front (top) loop: higher mean 3D z.
            let lo = &loops[0];
            let hi = &loops[1];
            let mean_z = |c: &CutLoop| {
                c.branch.points.iter().map(|p| p.z).sum::<f64>() / c.branch.points.len() as f64
            };
            let (front, trace) = if mean_z(hi) >= mean_z(lo) {
                (hi.target_face, hi.branch.uv_a.clone())
            } else {
                (lo.target_face, lo.branch.uv_a.clone())
            };
            front_face = Some(front);
            front_trace = trace;
        }
        (store, front_face.unwrap(), front_trace)
    }

    #[test]
    fn punched_face_has_one_hole_and_inner_wire() {
        let (store, face_id, _) = punch_front(0.7);
        let face = store.face(face_id).unwrap();
        let trim = face.trim.as_ref().expect("trim added");
        assert_eq!(trim.holes.len(), 1, "one hole on the front face");
        assert_eq!(face.inner_wires.len(), 1, "one inner wire");
    }

    #[test]
    fn hole_trim_winds_clockwise() {
        let (store, face_id, _) = punch_front(0.7);
        let face = store.face(face_id).unwrap();
        let hole = &face.trim.as_ref().unwrap().holes[0];
        // Sample the hole pcurves into a polygon and check the signed area is
        // negative (clockwise).
        let mut poly = Vec::new();
        for c in &hole.curves {
            let (t0, _t1) = c.parameter_domain();
            poly.push(c.point_at(t0).unwrap());
        }
        assert!(signed_area(&poly) < 0.0, "hole must wind clockwise");
    }

    #[test]
    fn punched_tessellation_excludes_hole_centroid() {
        let (store, face_id, trace) = punch_front(0.7);
        let mesh = TessellateFace::new(face_id, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty());

        // Hole centroid in UV.
        let cu = trace.iter().map(|p| p.x).sum::<f64>() / trace.len() as f64;
        let cv = trace.iter().map(|p| p.y).sum::<f64>() / trace.len() as f64;
        // Approximate hole radius in UV (max trace deviation).
        let hole_r_uv = trace
            .iter()
            .map(|p| ((p.x - cu).powi(2) + (p.y - cv).powi(2)).sqrt())
            .fold(0.0_f64, f64::max);

        let surface = nurbs_surface_of(&store, face_id).unwrap();
        let hole_center_3d = surface.point_at(cu, cv).unwrap();
        // No triangle centroid lies near the hole center in 3D (well inside the
        // hole). Use 40% of the hole's 3D radius as the exclusion band.
        let hole_r_3d = {
            let edge = surface.point_at((cu + hole_r_uv).min(1.0), cv).unwrap();
            (edge - hole_center_3d).norm()
        };
        for tri in &mesh.indices {
            let pa = mesh.vertices[tri[0] as usize];
            let pb = mesh.vertices[tri[1] as usize];
            let pc = mesh.vertices[tri[2] as usize];
            let cen = Point3::new(
                (pa.x + pb.x + pc.x) / 3.0,
                (pa.y + pb.y + pc.y) / 3.0,
                (pa.z + pb.z + pc.z) / 3.0,
            );
            let dist = (cen - hole_center_3d).norm();
            assert!(
                dist > hole_r_3d * 0.4,
                "triangle centroid {dist} too close to hole center (r3d={hole_r_3d})"
            );
        }
    }

    #[test]
    fn punching_twice_accumulates_two_holes() {
        // Build a slab and punch two separate tube loops onto the SAME face.
        let mut store = TopologyStore::new();
        let slab = MakeCurvedSlab::new(8.0, 0.0, 1.5, 1.0)
            .execute(&mut store)
            .unwrap();
        let slab_faces = solid_faces(&store, slab);
        let target = collect_nurbs_faces(&store, &slab_faces);

        // Two tubes at different XY positions, both through the slab.
        let punch_one = |store: &mut TopologyStore, cx: f64, cy: f64| -> FaceId {
            let tube = MakeNurbsTube::new(Point3::new(cx, cy, -1.5), 0.6, 5.0)
                .execute(store)
                .unwrap();
            let tool = collect_nurbs_faces(store, &solid_faces(store, tube));
            let cuts = extract_cut_loops(&target, &tool).unwrap();
            // Front loop = higher mean z.
            let ToolFaceCut::Through { loops, .. } = &cuts[0] else {
                panic!("expected a through cut");
            };
            let mean_z = |c: &CutLoop| {
                c.branch.points.iter().map(|p| p.z).sum::<f64>() / c.branch.points.len() as f64
            };
            let front = if mean_z(&loops[1]) >= mean_z(&loops[0]) {
                loops[1].clone()
            } else {
                loops[0].clone()
            };
            let face = front.target_face;
            punch_loop(store, &front).unwrap();
            face
        };

        let f1 = punch_one(&mut store, 2.5, 2.5);
        let f2 = punch_one(&mut store, 5.5, 5.5);
        assert_eq!(f1, f2, "both punches land on the same front face");
        let face = store.face(f1).unwrap();
        assert_eq!(
            face.trim.as_ref().unwrap().holes.len(),
            2,
            "two accumulated holes"
        );
        assert_eq!(face.inner_wires.len(), 2, "two inner wires");
    }
}
