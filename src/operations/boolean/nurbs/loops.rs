//! SSI loop extraction for the through-cut subtract.
//!
//! Runs surface-surface intersection over every (target NURBS face × tool NURBS
//! face) pair and classifies the cut loops. Closed branches (periodic tool
//! faces) are grouped per tool side face: exactly two closed loops (entry +
//! exit) form a through cut, a single loop a pocket. Open branches are
//! acceptable only when BOTH endpoints land on a tool kink edge; the
//! [`super::stitch`] module chains them across adjacent tool side faces into
//! closed loops, which are then grouped per chained-loop face set with the
//! same 2-loop (through) / 1-loop (pocket) classification. Any other
//! deviation is a typed unsupported error — never silent wrong geometry.
//!
//! Loops on geometrically closed surfaces arrive genuinely `closed` from the
//! SSI marcher (periodic-domain wrapping), with traces that reach both
//! parametric boundaries exactly at every seam crossing — no reclassification
//! or gap filling happens here. A loop whose trace crosses the TARGET face's
//! own seam cannot be punched as a simple polygon in the unrolled UV rectangle
//! and is rejected with a typed error (general boolean face splitting will
//! lift this).

use crate::error::{OperationError, Result};
use crate::geometry::nurbs::{
    intersect_surfaces, IntersectionOptions, NurbsSurface, SurfaceIntersectionCurve,
};
use crate::topology::{FaceId, FaceSurface, TopologyStore};

use super::stitch::{self, CutChain};

/// One intersection loop (or chained-loop segment) between a target face and
/// a tool side face.
#[derive(Debug, Clone)]
pub(crate) struct CutLoop {
    /// The target face this loop lies on (its `uv_a` trace is in target UV).
    pub target_face: FaceId,
    /// The tool side face this loop lies on (its `uv_b` trace is in tool UV).
    pub tool_face: FaceId,
    /// The SSI branch (`closed == true` for a whole periodic-face loop; an
    /// open kink-to-kink segment of a chained loop otherwise), with
    /// `uv_a`/`uv_b` synchronized to the 3D `points`.
    pub branch: SurfaceIntersectionCurve,
}

/// The cut class of one tool side face (periodic single-face loops) or of one
/// chained loop group (multi-face tools whose open branches were stitched
/// across kink edges).
#[derive(Debug, Clone)]
pub(crate) enum ToolFaceCut {
    /// Entry + exit loops, sorted by mean tool-`v` — the tool passes fully
    /// through the target.
    Through {
        tool_face: FaceId,
        loops: [CutLoop; 2],
    },
    /// A single entry loop — the tool enters the target and ends inside it
    /// (pocket / blind cut).
    Pocket { tool_face: FaceId, entry: CutLoop },
    /// Entry + exit chained loops crossing the same set of tool side faces,
    /// sorted by mean tool-`v` — a multi-face (box-like) tool passing fully
    /// through the target.
    MultiFaceThrough { chains: [CutChain; 2] },
    /// A single chained entry loop — a multi-face tool ending inside the
    /// target.
    MultiFacePocket { entry: CutChain },
}

/// The typed error for an open branch that cannot participate in kink-edge
/// chaining (message kept verbatim from the pre-chaining guard).
pub(crate) fn open_branch_error() -> crate::error::GeolisError {
    OperationError::Failed(
        "through-cut subtract requires closed intersection loops; \
         an open branch (partial cut / tool not passing fully through) \
         was found"
            .into(),
    )
    .into()
}

/// Extracts and validates the through-cut loops for `target` minus `tool`.
///
/// # Errors
///
/// Returns [`OperationError::Failed`] naming the unsupported case when: an
/// intersection branch is open without both endpoints on tool kink edges
/// (partial cut), a loop crosses the target face's parametric seam, no loops
/// are found at all (tool disjoint), a tool side face does not yield exactly
/// two closed loops, or open branches cannot be chained into closed loops.
/// (Cap-face intersection is guarded separately by the caller.)
pub(crate) fn extract_cut_loops(
    target_faces: &[(FaceId, NurbsSurface)],
    tool_faces: &[(FaceId, NurbsSurface)],
) -> Result<Vec<ToolFaceCut>> {
    let options = IntersectionOptions::default();
    let mut loops: Vec<CutLoop> = Vec::new();
    let mut open_segments: Vec<CutLoop> = Vec::new();

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
                // An open branch is acceptable only when BOTH endpoints sit
                // on a tool kink edge or a target face boundary (chained
                // across adjacent faces below); every other open branch keeps
                // the pre-chaining typed error.
                if !branch.closed
                    && !stitch::open_branch_admissible(&branch, target_surf, tool_surf)
                {
                    return Err(open_branch_error());
                }
                if crosses_target_seam(&branch, target_surf) {
                    return Err(OperationError::Failed(
                        "through-cut loop crosses the target face's parametric seam \
                         (unsupported until general boolean face splitting)"
                            .into(),
                    )
                    .into());
                }
                let cut = CutLoop {
                    target_face: *target_id,
                    tool_face: *tool_id,
                    branch,
                };
                if cut.branch.closed {
                    loops.push(cut);
                } else {
                    open_segments.push(cut);
                }
            }
        }
    }

    let mut cuts = group_per_tool_face(&loops, tool_faces)?;
    let chains = stitch::chain_open_segments(&open_segments, target_faces, tool_faces)?;
    cuts.extend(group_chains(chains, tool_faces)?);

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

/// Groups loops per tool side face and classifies each face's cut: two loops
/// form a through cut, a single loop a pocket. Three or more are unsupported.
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

        match mine.len() {
            // A tool side face that misses the target entirely is allowed only
            // when NO loops exist at all (tool disjoint from target, caught by
            // the caller). If some other tool face cut the target but this one
            // did not, the tool is not passing cleanly through — unsupported.
            0 => {}
            1 => {
                cuts.push(ToolFaceCut::Pocket {
                    tool_face: *tool_id,
                    entry: mine[0].clone(),
                });
            }
            2 => {
                // Order the two loops by mean v on the tool surface so the band
                // path can treat loops[0] as the lower (entry) and loops[1] as
                // the upper (exit).
                mine.sort_by(|a, b| {
                    mean_v_b(&a.branch)
                        .partial_cmp(&mean_v_b(&b.branch))
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
                let [lo, hi] = [mine[0].clone(), mine[1].clone()];
                cuts.push(ToolFaceCut::Through {
                    tool_face: *tool_id,
                    loops: [lo, hi],
                });
            }
            n => {
                return Err(OperationError::Failed(format!(
                    "NURBS subtract supports 1 (pocket) or 2 (through) closed \
                     loops per tool side face; tool face yielded {n}"
                ))
                .into());
            }
        }
    }
    Ok(cuts)
}

/// Groups chained loops by the SET of tool side faces they cross and
/// classifies each group: two chains form a through cut (sorted by mean
/// tool-`v`, the same convention as the single-face path), a single chain a
/// pocket. Groups are ordered deterministically by their lowest tool-face
/// index.
fn group_chains(
    chains: Vec<CutChain>,
    tool_faces: &[(FaceId, NurbsSurface)],
) -> Result<Vec<ToolFaceCut>> {
    use std::collections::BTreeMap;

    let index_of = |face: FaceId| -> usize {
        tool_faces
            .iter()
            .position(|(id, _)| *id == face)
            .unwrap_or(usize::MAX)
    };

    let mut groups: BTreeMap<Vec<usize>, Vec<CutChain>> = BTreeMap::new();
    for chain in chains {
        let mut key: Vec<usize> = chain
            .segments
            .iter()
            .map(|s| index_of(s.tool_face))
            .collect();
        key.sort_unstable();
        // A chain may cross one tool face in several segments (split by
        // target-face boundaries); the group identity is the face SET.
        key.dedup();
        groups.entry(key).or_default().push(chain);
    }

    let mut cuts = Vec::new();
    for (_, mut group) in groups {
        match group.len() {
            1 => {
                let entry = group.pop().unwrap_or_else(|| unreachable!());
                cuts.push(ToolFaceCut::MultiFacePocket { entry });
            }
            2 => {
                group.sort_by(|a, b| {
                    a.mean_v()
                        .partial_cmp(&b.mean_v())
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
                let hi = group.pop().unwrap_or_else(|| unreachable!());
                let lo = group.pop().unwrap_or_else(|| unreachable!());
                cuts.push(ToolFaceCut::MultiFaceThrough { chains: [lo, hi] });
            }
            n => {
                return Err(OperationError::Failed(format!(
                    "NURBS subtract supports 1 (pocket) or 2 (through) chained \
                     loops per tool face group; got {n}"
                ))
                .into());
            }
        }
    }
    Ok(cuts)
}

/// Whether a closed loop's target-UV trace crosses the target face's own
/// parametric seam: a geometrically closed target direction whose consecutive
/// trace samples jump more than half the period.
///
/// Such a loop is genuinely closed in 3D, but its punch polygon would have to
/// cross the edge of the target's unrolled UV rectangle — not a simple polygon.
/// The caller rejects it with a typed error until general boolean face
/// splitting can split the target at the seam.
fn crosses_target_seam(branch: &SurfaceIntersectionCurve, target: &NurbsSurface) -> bool {
    let ((u0, u1), (v0, v1)) = target.parameter_domain();
    let u_closed = target.is_closed_in_u();
    let v_closed = target.is_closed_in_v();
    if !u_closed && !v_closed {
        return false;
    }
    branch.uv_a.windows(2).any(|w| {
        (u_closed && (w[1].x - w[0].x).abs() > 0.5 * (u1 - u0))
            || (v_closed && (w[1].y - w[0].y).abs() > 0.5 * (v1 - v0))
    })
}

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
#[allow(clippy::unwrap_used, clippy::expect_used)]
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
        let ToolFaceCut::Through {
            loops: cut_loops, ..
        } = &cuts[0]
        else {
            panic!("expected a through cut");
        };
        // Both loops arrive genuinely closed from the marcher and carry exact
        // seam samples, so each `uv_b` trace spans the full tool u domain (it
        // reaches both the u0 and u1 boundaries across the seam).
        let tool_surf = &tool[0].1;
        let ((u0, u1), _) = tool_surf.parameter_domain();
        let u_span = u1 - u0;
        for loop_ in cut_loops {
            let umin = loop_
                .branch
                .uv_b
                .iter()
                .map(|p| p.x)
                .fold(f64::INFINITY, f64::min);
            let umax = loop_
                .branch
                .uv_b
                .iter()
                .map(|p| p.x)
                .fold(f64::NEG_INFINITY, f64::max);
            assert!(
                (umin - u0).abs() <= 1e-6 * u_span && (u1 - umax).abs() <= 1e-6 * u_span,
                "closed trace must reach both u boundaries: [{umin}, {umax}] \
                 vs domain [{u0}, {u1}]"
            );
        }
        // Each loop lies on a target face; the two target faces differ
        // (front + back of the slab).
        assert_ne!(
            cut_loops[0].target_face, cut_loops[1].target_face,
            "entry and exit loops lie on different slab faces"
        );
    }

    #[test]
    fn loop_crossing_target_seam_is_rejected() {
        use crate::geometry::nurbs::NurbsCurve3D;
        use crate::math::Vector3;
        use crate::operations::creation::{MakeNurbsPrism, MakeRevolvedSolid};

        let mut store = TopologyStore::new();
        // Plain cylindrical revolved wall (closed in u, seam azimuth at +X).
        let vase = MakeRevolvedSolid::new(vec![(2.0, 0.0), (2.0, 3.0)])
            .execute(&mut store)
            .unwrap();
        // Tube along +X: its entry hole straddles the revolved wall's
        // parametric seam. Punching that hole would need a trim polygon that
        // crosses the unrolled UV rectangle's edge — unsupported until general
        // boolean face splitting.
        let circle =
            NurbsCurve3D::circle(Point3::new(-4.0, 0.0, 1.5), 0.4, Vector3::x(), Vector3::y())
                .unwrap();
        let tube = MakeNurbsPrism::new(circle, Vector3::new(8.0, 0.0, 0.0))
            .execute(&mut store)
            .unwrap();
        let target = collect_nurbs_faces(&store, &solid_faces(&store, vase));
        let tool = collect_nurbs_faces(&store, &solid_faces(&store, tube));
        let result = extract_cut_loops(&target, &tool);
        let err = result.expect_err("seam-straddling hole must be a typed error");
        assert!(
            err.to_string().contains("parametric seam"),
            "error must name the target-seam limitation, got: {err}"
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
        let ToolFaceCut::Through { loops, .. } = &cuts[0] else {
            panic!("expected a through cut");
        };
        assert_ne!(
            loops[0].target_face, loops[1].target_face,
            "tilted entry/exit loops still land on different slab faces"
        );
    }

    /// A genuine 4-side-face box cutter through a segmented-prism wall: every
    /// (wall face × box face) SSI branch is OPEN (it ends on the box's kink
    /// edges), and the stitcher chains them into exactly two closed loops —
    /// the entry and exit windows — classified as one multi-face through cut.
    #[test]
    fn box_cutter_open_branches_chain_into_through_loops() {
        use crate::math::Vector3;
        use crate::operations::creation::{MakeSegmentedPrism, ProfileSegment};

        let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
        let line = |a: Point3, b: Point3| ProfileSegment::Line { start: a, end: b };

        let mut store = TopologyStore::new();
        // Wall: 6 x 0.4 footprint extruded 3 up.
        let wall = MakeSegmentedPrism::new(
            vec![
                line(p(0.0, 0.0, 0.0), p(6.0, 0.0, 0.0)),
                line(p(6.0, 0.0, 0.0), p(6.0, 0.4, 0.0)),
                line(p(6.0, 0.4, 0.0), p(0.0, 0.4, 0.0)),
                line(p(0.0, 0.4, 0.0), p(0.0, 0.0, 0.0)),
            ],
            Vector3::new(0.0, 0.0, 3.0),
        )
        .execute(&mut store)
        .unwrap();
        // Box cutter: window rectangle in the XZ plane at y = -1, extruded
        // horizontally through the wall.
        let cutter = MakeSegmentedPrism::new(
            vec![
                line(p(2.0, -1.0, 1.0), p(3.5, -1.0, 1.0)),
                line(p(3.5, -1.0, 1.0), p(3.5, -1.0, 2.0)),
                line(p(3.5, -1.0, 2.0), p(2.0, -1.0, 2.0)),
                line(p(2.0, -1.0, 2.0), p(2.0, -1.0, 1.0)),
            ],
            Vector3::new(0.0, 2.4, 0.0),
        )
        .execute(&mut store)
        .unwrap();

        let target = collect_nurbs_faces(&store, &solid_faces(&store, wall));
        let tool = collect_nurbs_faces(&store, &solid_faces(&store, cutter));
        let cuts = extract_cut_loops(&target, &tool).unwrap();
        assert_eq!(cuts.len(), 1, "one multi-face through cut");
        let ToolFaceCut::MultiFaceThrough { chains } = &cuts[0] else {
            panic!("expected a multi-face through cut, got {:?}", cuts[0]);
        };

        for chain in chains {
            assert_eq!(chain.segments.len(), 4, "one segment per box side face");
            // All four box side faces are distinct.
            for i in 0..4 {
                for j in (i + 1)..4 {
                    assert_ne!(chain.segments[i].tool_face, chain.segments[j].tool_face);
                }
            }
            // Junctions are welded EXACTLY: each segment's last sample is the
            // next segment's first sample (same 3D point, same target UV).
            for i in 0..4 {
                let a = &chain.segments[i];
                let b = &chain.segments[(i + 1) % 4];
                let pa = *a.branch.points.last().unwrap();
                let pb = *b.branch.points.first().unwrap();
                assert!(
                    (pa - pb).norm() == 0.0,
                    "junction {i} not welded: {pa:?} vs {pb:?}"
                );
                let ua = *a.branch.uv_a.last().unwrap();
                let ub = *b.branch.uv_a.first().unwrap();
                assert!((ua - ub).norm() == 0.0, "junction {i} target UV not welded");
            }
            // Every segment lies on the chain's single target face.
            assert!(
                chain.single_target_face().is_some(),
                "mid-wall box chain stays on one target face"
            );
        }

        // Entry (lower mean tool v) and exit land on the two opposite wall
        // faces, ordered by the same mean-v convention as the tube path.
        assert_ne!(
            chains[0].single_target_face(),
            chains[1].single_target_face()
        );
        assert!(chains[0].mean_v() < chains[1].mean_v());
    }

    /// F5 Phase C: a box window straddling a TARGET kink edge (two collinear
    /// wall segments sharing a vertical joint). The entry-side SSI branches
    /// end on the target faces' shared boundary and are chained ACROSS the
    /// target faces into one closed loop with exactly welded junctions.
    #[test]
    fn window_across_target_kink_chains_across_target_faces() {
        use crate::math::Vector3;
        use crate::operations::creation::{MakeSegmentedPrism, ProfileSegment};

        let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
        let line = |a: Point3, b: Point3| ProfileSegment::Line { start: a, end: b };

        let mut store = TopologyStore::new();
        // Wall: 6 x 0.4 footprint, height 3; the OUTER side is segmented into
        // two collinear pieces joined at x = 3 (a vertical target kink edge).
        let wall = MakeSegmentedPrism::new(
            vec![
                line(p(0.0, 0.0, 0.0), p(3.0, 0.0, 0.0)), // outer-a
                line(p(3.0, 0.0, 0.0), p(6.0, 0.0, 0.0)), // outer-b
                line(p(6.0, 0.0, 0.0), p(6.0, 0.4, 0.0)), // end-east
                line(p(6.0, 0.4, 0.0), p(0.0, 0.4, 0.0)), // inner
                line(p(0.0, 0.4, 0.0), p(0.0, 0.0, 0.0)), // end-west
            ],
            Vector3::new(0.0, 0.0, 3.0),
        )
        .execute(&mut store)
        .unwrap();
        // Box window x in [2, 3.5], z in [1, 2]: straddles the x = 3 joint.
        let cutter = MakeSegmentedPrism::new(
            vec![
                line(p(2.0, -1.0, 1.0), p(3.5, -1.0, 1.0)),
                line(p(3.5, -1.0, 1.0), p(3.5, -1.0, 2.0)),
                line(p(3.5, -1.0, 2.0), p(2.0, -1.0, 2.0)),
                line(p(2.0, -1.0, 2.0), p(2.0, -1.0, 1.0)),
            ],
            Vector3::new(0.0, 2.4, 0.0),
        )
        .execute(&mut store)
        .unwrap();

        let target = collect_nurbs_faces(&store, &solid_faces(&store, wall));
        let tool = collect_nurbs_faces(&store, &solid_faces(&store, cutter));
        let cuts = extract_cut_loops(&target, &tool).unwrap();
        assert_eq!(cuts.len(), 1, "one multi-face through cut");
        let ToolFaceCut::MultiFaceThrough { chains } = &cuts[0] else {
            panic!("expected a multi-face through cut, got {:?}", cuts[0]);
        };

        // Entry (outer side, lower mean tool v): crosses BOTH outer faces —
        // 6 segments (sill and head are split at the target kink).
        let entry = &chains[0];
        let exit = &chains[1];
        assert!(entry.mean_v() < exit.mean_v());
        assert!(
            entry.crosses_target_faces(),
            "entry chain crosses the two outer wall faces"
        );
        assert_eq!(entry.segments.len(), 6, "4 tool faces + 2 kink crossings");
        let entry_targets: std::collections::HashSet<FaceId> =
            entry.segments.iter().map(|s| s.target_face).collect();
        assert_eq!(entry_targets.len(), 2, "entry spans exactly 2 target faces");

        // Exit (inner side): a single unsegmented face — a plain 4-segment
        // chained loop as in Phase B.
        assert!(exit.single_target_face().is_some());
        assert_eq!(exit.segments.len(), 4);

        // Every junction (both kinds) is welded EXACTLY: same 3D point on
        // both sides, and target-boundary junctions pin the target UV on its
        // domain bound.
        for chain in [entry, exit] {
            let n = chain.segments.len();
            for i in 0..n {
                let a = &chain.segments[i];
                let b = &chain.segments[(i + 1) % n];
                let pa = *a.branch.points.last().unwrap();
                let pb = *b.branch.points.first().unwrap();
                assert!(
                    (pa - pb).norm() == 0.0,
                    "junction {i} not welded: {pa:?} vs {pb:?}"
                );
                if a.target_face != b.target_face {
                    // Target-boundary junction: both target UVs pinned
                    // exactly on their domain bounds.
                    for (seg, uv) in [
                        (a, *a.branch.uv_a.last().unwrap()),
                        (b, *b.branch.uv_a.first().unwrap()),
                    ] {
                        let surf = &target
                            .iter()
                            .find(|(id, _)| *id == seg.target_face)
                            .unwrap()
                            .1;
                        let ((u0, u1), _) = surf.parameter_domain();
                        assert!(
                            (uv.x - u0).abs() == 0.0 || (u1 - uv.x).abs() == 0.0,
                            "target UV not pinned on its u bound: {uv:?}"
                        );
                    }
                }
            }
        }
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
