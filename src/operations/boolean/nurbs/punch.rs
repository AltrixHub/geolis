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
//!   polylines anyway, so this is lossless for tessellation.
//! - The SSI trace omits the short arc at the tool's parametric seam (the
//!   marcher terminates there). The loop is closed with a single straight
//!   segment across that gap. The gap is one marching step wide, so the induced
//!   error is sub-step and bounded.

use crate::error::{OperationError, Result};
use crate::geometry::nurbs::{KnotVector, NurbsCurve2D, NurbsCurve3D, NurbsSurface};
use crate::math::{Point2, Point3, TOLERANCE};
use crate::topology::{
    EdgeCurve, EdgeData, FaceId, FaceSurface, FaceTrim, OrientedEdge, TopologyStore, TrimLoop,
    VertexData, WireData, WireId,
};

use super::loops::CutLoop;

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
    let inner_wire = build_inner_wire(store, &cut.branch.points)?;

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
/// segments, closing the seam gap with a straight segment.
fn hole_loop_from_trace(uv_a: &[Point2]) -> Result<TrimLoop> {
    let mut pts = dedup_uv(uv_a);
    if pts.len() < 3 {
        return Err(OperationError::Failed(
            "punched loop degenerated to fewer than 3 distinct UV points".into(),
        )
        .into());
    }
    // Hole loops must wind clockwise (negative signed area).
    if signed_area(&pts) > 0.0 {
        pts.reverse();
    }

    let n = pts.len();
    let mut curves = Vec::with_capacity(n);
    for i in 0..n {
        let a = pts[i];
        let b = pts[(i + 1) % n];
        // The wrap-around segment (i = n-1) closes the seam gap.
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

/// Builds a closed 3D inner wire from the SSI 3D trace via a single degree-3
/// interpolated NURBS edge (closing the seam gap).
fn build_inner_wire(store: &mut TopologyStore, points: &[Point3]) -> Result<WireId> {
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
    use crate::operations::boolean::nurbs::loops::{collect_nurbs_faces, extract_cut_loops};
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
            for cut in &loop_pair.loops {
                punch_loop(&mut store, cut).unwrap();
            }
            // Identify the front (top) loop: higher mean 3D z.
            let lo = &loop_pair.loops[0];
            let hi = &loop_pair.loops[1];
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
        let r = trace
            .iter()
            .map(|p| ((p.x - cu).powi(2) + (p.y - cv).powi(2)).sqrt())
            .fold(0.0_f64, f64::max);

        let surface = nurbs_surface_of(&store, face_id).unwrap();
        let hole_center_3d = surface.point_at(cu, cv).unwrap();
        // No triangle centroid lies near the hole center in 3D (well inside the
        // hole). Use 40% of the hole's 3D radius as the exclusion band.
        let hole_r_3d = {
            let edge = surface.point_at((cu + r).min(1.0), cv).unwrap();
            (edge - hole_center_3d).norm()
        };
        for t in &mesh.indices {
            let a = mesh.vertices[t[0] as usize];
            let b = mesh.vertices[t[1] as usize];
            let c = mesh.vertices[t[2] as usize];
            let cen = Point3::new(
                (a.x + b.x + c.x) / 3.0,
                (a.y + b.y + c.y) / 3.0,
                (a.z + b.z + c.z) / 3.0,
            );
            let d = (cen - hole_center_3d).norm();
            assert!(
                d > hole_r_3d * 0.4,
                "triangle centroid {d} too close to hole center (r3d={hole_r_3d})"
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
            let pair = &cuts[0];
            let mean_z = |c: &CutLoop| {
                c.branch.points.iter().map(|p| p.z).sum::<f64>() / c.branch.points.len() as f64
            };
            let front = if mean_z(&pair.loops[1]) >= mean_z(&pair.loops[0]) {
                pair.loops[1].clone()
            } else {
                pair.loops[0].clone()
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
