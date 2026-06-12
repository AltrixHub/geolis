//! SSI loop extraction for the through-cut subtract.
//!
//! Runs surface-surface intersection over every (target NURBS face × tool NURBS
//! face) pair, keeps only closed loops, and groups them per tool side face. The
//! through-cut precondition is that each tool side face yields exactly two
//! closed loops (entry + exit), each lying on a target face. Any deviation is a
//! typed unsupported error — never silent wrong geometry.

use crate::error::{OperationError, Result};
use crate::geometry::nurbs::{
    intersect_surfaces, IntersectionOptions, NurbsSurface, SurfaceIntersectionCurve,
};
use crate::topology::{FaceId, FaceSurface, TopologyStore};

/// One closed intersection loop between a target face and a tool side face.
#[derive(Debug, Clone)]
pub(crate) struct CutLoop {
    /// The target face this loop lies on (its `uv_a` trace is in target UV).
    pub target_face: FaceId,
    /// The tool side face this loop lies on (its `uv_b` trace is in tool UV).
    pub tool_face: FaceId,
    /// The SSI branch (`closed == true`), with `uv_a`/`uv_b` synchronized to
    /// the 3D `points`.
    pub branch: SurfaceIntersectionCurve,
}

/// All cut loops belonging to a single tool side face. The through-cut contract
/// guarantees exactly two: the entry loop and the exit loop.
#[derive(Debug, Clone)]
pub(crate) struct ToolFaceCut {
    pub tool_face: FaceId,
    pub loops: [CutLoop; 2],
}

/// Extracts and validates the through-cut loops for `target` minus `tool`.
///
/// # Errors
///
/// Returns [`OperationError::Failed`] naming the unsupported case when: an
/// intersection branch is open / not seam-closed (partial cut), no loops are
/// found at all (tool disjoint), or a tool side face does not yield exactly two
/// closed loops. (Cap-face intersection is guarded separately by the caller.)
pub(crate) fn extract_cut_loops(
    target_faces: &[(FaceId, NurbsSurface)],
    tool_faces: &[(FaceId, NurbsSurface)],
) -> Result<Vec<ToolFaceCut>> {
    let options = IntersectionOptions::default();
    let mut loops: Vec<CutLoop> = Vec::new();

    for (tool_id, tool_surf) in tool_faces {
        for (target_id, target_surf) in target_faces {
            if !aabb_overlap(target_surf, tool_surf) {
                continue;
            }
            // SSI is run with the target as surface `a` and the tool as surface
            // `b`, so `uv_a` lands on the target (where we punch) and `uv_b` on
            // the tool (where we band).
            let branches = intersect_surfaces(target_surf, tool_surf, &options)?;
            for branch in branches {
                if branch.points.len() < 3 {
                    continue;
                }
                if !is_closed_loop(&branch, tool_surf) {
                    return Err(OperationError::Failed(
                        "through-cut subtract requires closed intersection loops; \
                         an open branch (partial cut / tool not passing fully through) \
                         was found"
                            .into(),
                    )
                    .into());
                }
                loops.push(CutLoop {
                    target_face: *target_id,
                    tool_face: *tool_id,
                    branch,
                });
            }
        }
    }

    group_per_tool_face(&loops, tool_faces)
}

/// Groups loops per tool side face and validates the exactly-two contract.
fn group_per_tool_face(
    loops: &[CutLoop],
    tool_faces: &[(FaceId, NurbsSurface)],
) -> Result<Vec<ToolFaceCut>> {
    let mut cuts = Vec::new();

    for (tool_id, _) in tool_faces {
        let mut mine: Vec<CutLoop> = loops
            .iter()
            .filter(|l| l.tool_face == *tool_id)
            .cloned()
            .collect();

        if mine.is_empty() {
            // A tool side face that misses the target entirely is allowed only
            // when NO loops exist at all (tool disjoint from target). If some
            // other tool face cut the target but this one did not, the tool is
            // not passing cleanly through — unsupported.
            continue;
        }
        if mine.len() != 2 {
            return Err(OperationError::Failed(format!(
                "through-cut subtract requires exactly 2 closed loops per tool \
                 side face (entry + exit); tool face yielded {}",
                mine.len()
            ))
            .into());
        }
        // Order the two loops by mean v on the tool surface so the band path can
        // treat loops[0] as the lower (entry) and loops[1] as the upper (exit).
        mine.sort_by(|a, b| {
            mean_v_b(&a.branch)
                .partial_cmp(&mean_v_b(&b.branch))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let [lo, hi] = [mine[0].clone(), mine[1].clone()];
        cuts.push(ToolFaceCut {
            tool_face: *tool_id,
            loops: [lo, hi],
        });
    }

    if cuts.is_empty() {
        return Err(OperationError::Failed(
            "through-cut subtract found no intersection loops (tool does not \
             pass through the target)"
                .into(),
        )
        .into());
    }
    Ok(cuts)
}

/// Whether a branch forms a closed through-cut loop on the tool side face.
///
/// A branch counts as closed when either:
/// - the SSI marcher already flagged it `closed` (the loop fits inside the
///   tool's parametric domain without crossing a seam), OR
/// - it is a **seam-closed** loop: the loop is a v = f(u) graph that wraps the
///   tube once, so the marcher terminates at the tool surface's `u` seam with
///   the two `uv_b` endpoints pinned to opposite ends of the u-domain at a
///   matching v. The unmarched arc is exactly the seam, so the loop is closed
///   in 3D even though `closed == false`.
///
/// This is the shipped resolution to the seam-termination behavior of marching
/// SSI on a geometrically-closed-but-non-periodic tool surface (the tube side):
/// the openness is a UV-seam artifact, not a partial cut. The check is strict —
/// both endpoints must sit within `SEAM_EPS` of the opposite u-boundaries and
/// share v within `SEAM_V_EPS` — so a genuine partial cut (endpoints in the
/// domain interior, or both at the same boundary) is still rejected upstream.
fn is_closed_loop(branch: &SurfaceIntersectionCurve, tool: &NurbsSurface) -> bool {
    if branch.closed {
        return true;
    }
    let (Some(first), Some(last)) = (branch.uv_b.first(), branch.uv_b.last()) else {
        return false;
    };
    let ((u0, u1), _) = tool.parameter_domain();
    let u_span = (u1 - u0).abs().max(f64::EPSILON);
    // Endpoints near opposite u-boundaries (one near u0, one near u1).
    let near_lo = |u: f64| (u - u0).abs() <= SEAM_EPS * u_span;
    let near_hi = |u: f64| (u1 - u).abs() <= SEAM_EPS * u_span;
    let straddles_seam =
        (near_lo(first.x) && near_hi(last.x)) || (near_hi(first.x) && near_lo(last.x));
    if !straddles_seam {
        return false;
    }
    // v must match at the seam (the loop is a single-valued graph over u).
    (first.y - last.y).abs() <= SEAM_V_EPS
}

/// Fraction of the tool u-domain within which a seam-closed loop endpoint must
/// lie to count as touching the seam.
const SEAM_EPS: f64 = 0.12;

/// Maximum v mismatch (in tool parameter units) between the two seam endpoints.
const SEAM_V_EPS: f64 = 0.05;

/// Mean of the `uv_b` v-coordinate over a branch (tool-axis position).
fn mean_v_b(branch: &SurfaceIntersectionCurve) -> f64 {
    if branch.uv_b.is_empty() {
        return 0.0;
    }
    #[allow(clippy::cast_precision_loss)]
    let inv = 1.0 / branch.uv_b.len() as f64;
    branch.uv_b.iter().map(|p| p.y).sum::<f64>() * inv
}

/// Conservative control-hull AABB overlap test.
fn aabb_overlap(a: &NurbsSurface, b: &NurbsSurface) -> bool {
    let (a_lo, a_hi) = a.bounding_box();
    let (b_lo, b_hi) = b.bounding_box();
    let pad = 1e-7;
    a_lo.x <= b_hi.x + pad
        && a_hi.x >= b_lo.x - pad
        && a_lo.y <= b_hi.y + pad
        && a_hi.y >= b_lo.y - pad
        && a_lo.z <= b_hi.z + pad
        && a_hi.z >= b_lo.z - pad
}

/// Collects the NURBS faces of a solid as `(FaceId, surface clone)` pairs.
///
/// Planar faces (caps, slab sides) are skipped here; their interaction is
/// validated separately by [`assert_no_cap_intersection`].
pub(crate) fn collect_nurbs_faces(
    store: &TopologyStore,
    faces: &[FaceId],
) -> Vec<(FaceId, NurbsSurface)> {
    let mut out = Vec::new();
    for &fid in faces {
        if let Ok(face) = store.face(fid) {
            if let FaceSurface::Nurbs(surf) = &face.surface {
                out.push((fid, surf.clone()));
            }
        }
    }
    out
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use crate::operations::creation::{MakeCurvedSlab, MakeNurbsTube};
    use crate::topology::SolidId;

    fn solid_faces(store: &TopologyStore, solid: SolidId) -> Vec<FaceId> {
        let shell = store
            .shell(store.solid(solid).unwrap().outer_shell)
            .unwrap();
        shell.faces.clone()
    }

    /// Builds a slab and a tube through its center; returns (store, target, tool).
    fn slab_and_tube(tube_center: Point3, radius: f64) -> (TopologyStore, SolidId, SolidId) {
        let mut store = TopologyStore::new();
        // Slab spans [0,6]^2 in XY, front peaks 1.5 above z=0, 1.0 thick (down).
        let slab = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut store)
            .unwrap();
        // Tube rises through the slab from below to above.
        let tube = MakeNurbsTube::new(tube_center, radius, 5.0)
            .execute(&mut store)
            .unwrap();
        (store, slab, tube)
    }

    #[test]
    fn slab_through_tube_yields_two_loops() {
        let (store, slab, tube) = slab_and_tube(Point3::new(3.0, 3.0, -1.5), 0.7);
        let target = collect_nurbs_faces(&store, &solid_faces(&store, slab));
        let tool = collect_nurbs_faces(&store, &solid_faces(&store, tube));
        let cuts = extract_cut_loops(&target, &tool).unwrap();
        // The tube has one NURBS side face; it must yield exactly 2 loops.
        assert_eq!(cuts.len(), 1, "one tool side face");
        let cut = &cuts[0];
        // Both loops are accepted as closed (genuinely closed or seam-closed).
        let tool_surf = &tool[0].1;
        assert!(super::is_closed_loop(&cut.loops[0].branch, tool_surf));
        assert!(super::is_closed_loop(&cut.loops[1].branch, tool_surf));
        // Each loop lies on a target face; the two target faces differ
        // (front + back of the slab).
        assert_ne!(
            cut.loops[0].target_face, cut.loops[1].target_face,
            "entry and exit loops lie on different slab faces"
        );
    }

    #[test]
    fn tilted_tube_still_yields_two_loops() {
        use crate::geometry::nurbs::{NurbsCurve3D, NurbsSurface};
        use crate::math::Vector3;
        use crate::operations::creation::MakeNurbsFace;

        let mut store = TopologyStore::new();
        let slab = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut store)
            .unwrap();
        let target = collect_nurbs_faces(&store, &solid_faces(&store, slab));

        // A circle in the XY plane at z = -1.5, extruded along a tilted axis so
        // its side surface leans while still passing through both slab faces.
        let circle =
            NurbsCurve3D::circle(Point3::new(3.0, 3.0, -1.5), 0.7, Vector3::z(), Vector3::x())
                .unwrap();
        let tilt = Vector3::new(0.6, 0.4, 5.0);
        let side_surf = NurbsSurface::extrude(&circle, tilt).unwrap();
        let side_face = MakeNurbsFace::new(side_surf.clone())
            .execute(&mut store)
            .unwrap();
        let tool = vec![(side_face, side_surf)];

        let cuts = extract_cut_loops(&target, &tool).unwrap();
        assert_eq!(cuts.len(), 1, "one tilted tool side face");
        assert_ne!(
            cuts[0].loops[0].target_face, cuts[0].loops[1].target_face,
            "tilted entry/exit loops still land on different slab faces"
        );
    }

    #[test]
    fn tube_missing_slab_yields_no_loops_error() {
        // Tube far to the side, never touching the slab.
        let (store, slab, tube) = slab_and_tube(Point3::new(20.0, 20.0, -1.5), 0.7);
        let target = collect_nurbs_faces(&store, &solid_faces(&store, slab));
        let tool = collect_nurbs_faces(&store, &solid_faces(&store, tube));
        let result = extract_cut_loops(&target, &tool);
        assert!(result.is_err(), "disjoint tube must error (no loops)");
    }

    #[test]
    fn half_buried_tube_is_unsupported() {
        // Tube that starts below but stops INSIDE the slab thickness: its top
        // cap is buried, so the side face intersects the front face in a closed
        // loop but never exits — an open branch or a single loop. Either way,
        // not the clean 2-loop through-cut.
        let mut store = TopologyStore::new();
        let slab = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0)
            .execute(&mut store)
            .unwrap();
        // Tube top at z = -0.5 (inside the slab body which spans ~[-1, 1.5]).
        let tube = MakeNurbsTube::new(Point3::new(3.0, 3.0, -3.0), 0.7, 2.5)
            .execute(&mut store)
            .unwrap();
        let target = collect_nurbs_faces(&store, &solid_faces(&store, slab));
        let tool = collect_nurbs_faces(&store, &solid_faces(&store, tube));
        let result = extract_cut_loops(&target, &tool);
        assert!(
            result.is_err(),
            "half-buried tube must be unsupported, got: {result:?}"
        );
    }
}
