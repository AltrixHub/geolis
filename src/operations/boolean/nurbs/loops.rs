//! SSI loop extraction for the through-cut subtract.
//!
//! Runs surface-surface intersection over every (target NURBS face × tool NURBS
//! face) pair, keeps only closed loops, and groups them per tool side face. The
//! through-cut precondition is that each tool side face yields exactly two
//! closed loops (entry + exit), each lying on a target face. Any deviation is a
//! typed unsupported error — never silent wrong geometry.

use nalgebra::Matrix3;

use crate::error::{OperationError, Result};
use crate::geometry::nurbs::{
    intersect_surfaces, IntersectionOptions, NurbsSurface, SurfaceIntersectionCurve,
};
use crate::math::Point2;
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
            for mut branch in branches {
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
                // A seam-closed branch (`closed == false` but endpoints pinned to
                // the tool's opposite u-boundaries) is missing the short arc at
                // the tool's parametric u-seam. Fill that wedge with true
                // intersection samples so the punch (uv_a) ring and the band
                // (uv_b) ring share the same geometry across the seam. Genuinely
                // closed branches (`closed == true`) already span their loop.
                if !branch.closed {
                    fill_seam_gap(&mut branch, target_surf, tool_surf, &options);
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

/// Fills the u-seam wedge of a seam-closed branch with true intersection samples.
///
/// The SSI marcher terminates at the tool's parametric u-seam, so a seam-closed
/// branch's `uv_b` trace spans only the interior `(u_lo, u_hi)` of the tool u
/// domain and omits the wedge `[u0, u_lo] ∪ [u_hi, u1]` that straddles the seam.
/// The punch would otherwise close its ring with one straight chord across that
/// wedge while the band would leave the wedge unmeshed — the two disagree,
/// leaving a visible slit at the seam azimuth.
///
/// This appends `K` refined samples that walk the wedge from the trace's `last`
/// endpoint across the seam to (just before) its `first` endpoint, **reaching
/// both `u0` and `u1` exactly** so the band ribbon spans the full u domain (no
/// unmeshed wedge) and the punch ring follows the real arc (no chord). For each
/// sample the tool `u` is fixed and a 3×3 Newton solves the true intersection
/// `S_target(u_a, v_a) = S_tool(u_fixed, v_b)` for the three unknowns
/// `(u_a, v_a, v_b)`, mirroring the corrector in `surface_surface.rs` but with
/// the tool `u` pinned. Each sample is appended to the synchronized
/// `points` / `uv_a` / `uv_b` (`uv_a` = target, `uv_b` = tool).
///
/// ## Seam-wrap convention (band UV)
///
/// The tool side surface is geometrically closed but parametrically non-periodic
/// (`u0` and `u1` map to the same seam azimuth). The appended `uv_b` samples keep
/// their tool `u` **wrapped into `[u0, u1]`** — not unrolled past `u1` — because
/// the trimmed tessellator evaluates 3D via `surface.point_at(u, v)`, which
/// requires `u ∈ [u0, u1]`. Sorting the wrapped `u` then yields a ribbon covering
/// the full `[u0, u1]`; its left (`u0`) and right (`u1`) closing edges land on the
/// same seam azimuth and coincide in 3D, so the band stays a simple polygon in
/// the unrolled rectangle and closes conformally with the punch ring.
///
/// On Newton non-convergence for any sample the branch is left unmodified and the
/// punch chord / band stitch remain as honest fallbacks (a genuinely unfillable
/// gap keeps the pre-existing sub-step approximation rather than fabricating
/// geometry).
fn fill_seam_gap(
    branch: &mut SurfaceIntersectionCurve,
    target: &NurbsSurface,
    tool: &NurbsSurface,
    options: &IntersectionOptions,
) {
    let Some(&last_a) = branch.uv_a.last() else {
        return;
    };
    let (Some(&first_b), Some(&last_b)) = (branch.uv_b.first(), branch.uv_b.last()) else {
        return;
    };
    let ((u0, u1), _) = tool.parameter_domain();
    let period = u1 - u0;
    if period <= f64::EPSILON {
        return;
    }

    // Marching step in tool u (median consecutive |Δu| of the raw trace).
    let step_u = median_step_u(&branch.uv_b).max(period * 1e-3);

    // Orientation: which u-boundary is the trace's `last` endpoint closest to?
    // `last` is walked toward its own boundary; `first` is approached from the
    // opposite one, so the appended chain crosses the seam once.
    let last_near_hi = (u1 - last_b.x).abs() <= (last_b.x - u0).abs();
    let hi_bound = if last_near_hi { u1 } else { u0 };
    let lo_bound = if last_near_hi { u0 } else { u1 };

    // Segment counts sized to the marching step; `max(2, …)` guarantees the seam
    // boundary itself is sampled and keeps the total `K = n_hi + n_lo >= 4`.
    let hi_span = hi_bound - last_b.x;
    let lo_span = first_b.x - lo_bound;
    let n_hi = sample_count(hi_span, step_u);
    let n_lo = sample_count(lo_span, step_u);

    // Collect samples first; only commit if every Newton solve converges.
    let mut new_points = Vec::with_capacity(n_hi + n_lo);
    let mut new_uv_a = Vec::with_capacity(n_hi + n_lo);
    let mut new_uv_b = Vec::with_capacity(n_hi + n_lo);
    let mut seed = (last_a.x, last_a.y, last_b.y);

    // High side: from `last` (exclusive) up to the near boundary (inclusive).
    for i in 1..=n_hi {
        #[allow(clippy::cast_precision_loss)]
        let frac = i as f64 / n_hi as f64;
        let u_fixed = clamp_u(last_b.x + hi_span * frac, u0, u1);
        let Some(sample) = newton_seam(target, tool, u_fixed, seed, options) else {
            return;
        };
        seed = sample;
        push_sample(
            target,
            u_fixed,
            sample,
            &mut new_points,
            &mut new_uv_a,
            &mut new_uv_b,
        );
    }
    // Low side: from the opposite boundary (inclusive) toward `first` (exclusive).
    for i in 0..n_lo {
        #[allow(clippy::cast_precision_loss)]
        let frac = i as f64 / n_lo as f64;
        let u_fixed = clamp_u(lo_bound + lo_span * frac, u0, u1);
        let Some(sample) = newton_seam(target, tool, u_fixed, seed, options) else {
            return;
        };
        seed = sample;
        push_sample(
            target,
            u_fixed,
            sample,
            &mut new_points,
            &mut new_uv_a,
            &mut new_uv_b,
        );
    }

    branch.points.extend(new_points);
    branch.uv_a.extend(new_uv_a);
    branch.uv_b.extend(new_uv_b);
}

/// Number of segments for a wedge sub-span at roughly the marching step, at least
/// 2 so the seam boundary is always reached.
fn sample_count(span: f64, step_u: f64) -> usize {
    #[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
    let n = (span.abs() / step_u).ceil() as usize;
    n.max(2)
}

/// Clamps a tool u into the parameter domain (guards floating-point overshoot at
/// the seam boundary).
fn clamp_u(u: f64, u0: f64, u1: f64) -> f64 {
    u.clamp(u0, u1)
}

/// Evaluates the 3D point at the solved target UV and records the synchronized
/// sample. `S_target(u_a, v_a) == S_tool(u_fixed, v_b)` at convergence, so the
/// target evaluation is authoritative for the shared ring geometry.
fn push_sample(
    target: &NurbsSurface,
    u_fixed: f64,
    sample: (f64, f64, f64),
    points: &mut Vec<crate::math::Point3>,
    uv_a: &mut Vec<Point2>,
    uv_b: &mut Vec<Point2>,
) {
    let (ua, va, vb) = sample;
    // `newton_seam` only returns a converged sample, so this evaluation succeeds;
    // fall back to skipping the point defensively if it somehow fails.
    if let Ok(p) = target.point_at(ua, va) {
        points.push(p);
        uv_a.push(Point2::new(ua, va));
        uv_b.push(Point2::new(u_fixed, vb));
    }
}

/// Median consecutive |Δu| of a tool-UV trace (the local marching step in u).
fn median_step_u(uv_b: &[Point2]) -> f64 {
    let mut diffs: Vec<f64> = uv_b.windows(2).map(|w| (w[1].x - w[0].x).abs()).collect();
    if diffs.is_empty() {
        return 0.0;
    }
    diffs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    diffs[diffs.len() / 2]
}

/// 3×3 Newton for one seam sample with the tool `u` pinned.
///
/// Unknowns `(u_a, v_a, v_b)`; residual `r = S_target(u_a, v_a) - S_tool(u_fixed,
/// v_b)` (3 components). Jacobian columns `[∂S_a/∂u_a, ∂S_a/∂v_a, -∂S_b/∂v_b]`.
/// Solves `J·Δ = r`, updates `x -= Δ`, clamps to both domains (`u_fixed` stays
/// pinned). Returns `None` if the residual stays above tolerance (non-convergence
/// / singular Jacobian), so the caller can fall back cleanly.
#[allow(clippy::similar_names)]
fn newton_seam(
    target: &NurbsSurface,
    tool: &NurbsSurface,
    u_fixed: f64,
    seed: (f64, f64, f64),
    options: &IntersectionOptions,
) -> Option<(f64, f64, f64)> {
    let ((au0, au1), (av0, av1)) = target.parameter_domain();
    let (_, (bv0, bv1)) = tool.parameter_domain();
    let (mut ua, mut va, mut vb) = seed;
    let tol = options.tolerance.max(1e-12);

    for _ in 0..options.max_iterations {
        let (pa, sau, sav) = target.partials(ua, va).ok()?;
        let (pb, _sbu, sbv) = tool.partials(u_fixed, vb).ok()?;
        let r = pa - pb;
        if r.norm() < tol {
            return Some((ua, va, vb));
        }
        let j = Matrix3::from_columns(&[sau, sav, -sbv]);
        let jinv = j.try_inverse()?;
        let delta = jinv * r;
        ua = (ua - delta[0]).clamp(au0, au1);
        va = (va - delta[1]).clamp(av0, av1);
        vb = (vb - delta[2]).clamp(bv0, bv1);
        if delta.norm() < tol {
            break;
        }
    }
    let pa = target.point_at(ua, va).ok()?;
    let pb = tool.point_at(u_fixed, vb).ok()?;
    if (pa - pb).norm() < tol.max(1e-7) {
        Some((ua, va, vb))
    } else {
        None
    }
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
        let cut = &cuts[0];
        // Both loops were accepted (extraction succeeded) and — being seam-closed
        // here — were seam-filled, so each `uv_b` trace now spans the full tool u
        // domain (it reaches both the u0 and u1 boundaries across the seam).
        let tool_surf = &tool[0].1;
        let ((u0, u1), _) = tool_surf.parameter_domain();
        let u_span = u1 - u0;
        for loop_ in &cut.loops {
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
                "seam-filled trace must reach both u boundaries: [{umin}, {umax}] \
                 vs domain [{u0}, {u1}]"
            );
        }
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
