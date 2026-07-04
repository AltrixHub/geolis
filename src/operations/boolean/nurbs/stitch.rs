//! Open-branch stitcher: chains open SSI branches across adjacent tool side
//! faces into closed cut loops (F5 Phase B).
//!
//! A multi-face (box-like) tool has planar side faces meeting at kink edges.
//! Each side face's intersection with a target face is an OPEN SSI branch
//! whose endpoints land exactly on the tool's kink edges (the marcher pins
//! open-boundary endpoints — see the SSI boundary-pinned termination). This
//! module accepts such branches, chains them end-to-end across adjacent tool
//! faces into closed [`CutChain`]s, and welds the junction samples so the
//! punched hole ring, the per-face band fragments, and the kink-crossing
//! edges all share exactly the same junction geometry.
//!
//! Acceptance is structural, not tolerance-driven: an open branch is
//! acceptable iff BOTH endpoints sit on an open `u` boundary of the tool face
//! (within the marcher's own [`BOUNDARY_EPS`]) while staying strictly inside
//! the target face's domain, and every endpoint must find exactly one partner
//! endpoint from an ADJACENT (different) tool face at the same 3D point. Any
//! violation keeps the pre-chaining typed errors verbatim.
//!
//! F3b (target-side face splitting, Phase C) reuses this machinery with the
//! roles of target and tool boundaries swapped.

use crate::error::Result;
use crate::geometry::nurbs::{NurbsSurface, SurfaceIntersectionCurve, BOUNDARY_EPS};
use crate::math::Point3;
use crate::operations::boolean::nurbs::loops::open_branch_error;
use crate::topology::FaceId;

use super::loops::CutLoop;

/// 3D junction coincidence bound: the SSI marcher's point-acceptance bound
/// (`IntersectionOptions.tolerance.max(1e-7)` in the corrector's final
/// check). Two branch endpoints converged onto the same tool-kink crossing
/// agree within this — no new tolerance is introduced, and Newton typically
/// converges orders of magnitude tighter.
pub(crate) const JUNCTION_TOLERANCE: f64 = 1e-7;

/// A closed cut loop chained from open SSI branches across adjacent tool
/// side faces. `segments[i]` runs on one tool face; segment `i`'s last sample
/// coincides exactly (welded) with segment `i + 1`'s first sample (cyclic).
#[derive(Debug, Clone)]
pub(crate) struct CutChain {
    /// The single target face every segment of this loop lies on.
    pub target_face: FaceId,
    /// The chained per-tool-face segments, in cyclic order.
    pub segments: Vec<CutLoop>,
}

impl CutChain {
    /// Mean tool-`v` over all segment samples (the loop's position along the
    /// tool axis) — the same ordering key the single-face path uses.
    pub(crate) fn mean_v(&self) -> f64 {
        let mut sum = 0.0;
        let mut count = 0usize;
        for seg in &self.segments {
            for p in &seg.branch.uv_b {
                sum += p.y;
                count += 1;
            }
        }
        if count == 0 {
            0.0
        } else {
            #[allow(clippy::cast_precision_loss)]
            let inv = 1.0 / count as f64;
            sum * inv
        }
    }
}

/// Whether an open SSI branch qualifies for kink-edge chaining: both
/// endpoints terminate on an open `u` boundary of the TOOL face (a kink-edge
/// candidate), away from the tool's `v` boundaries (cap rings) and strictly
/// inside the TARGET face's domain (endpoints on target boundaries are the
/// F3b face-splitting case, still unsupported).
pub(crate) fn open_branch_on_tool_kinks(
    branch: &SurfaceIntersectionCurve,
    target: &NurbsSurface,
    tool: &NurbsSurface,
) -> bool {
    let ((tu0, tu1), (tv0, tv1)) = tool.parameter_domain();
    let ((au0, au1), (av0, av1)) = target.parameter_domain();
    if tool.is_closed_in_u() {
        // A periodic tool direction has no kink boundary to chain across.
        return false;
    }
    let ends = [
        (branch.uv_b.first(), branch.uv_a.first()),
        (branch.uv_b.last(), branch.uv_a.last()),
    ];
    for (uv_b, uv_a) in ends {
        let (Some(b), Some(a)) = (uv_b, uv_a) else {
            return false;
        };
        let on_u_bound = (b.x - tu0).abs() < BOUNDARY_EPS || (tu1 - b.x).abs() < BOUNDARY_EPS;
        if !on_u_bound {
            return false;
        }
        // A corner endpoint (also on a tool cap ring) is degenerate.
        if !tool.is_closed_in_v()
            && ((b.y - tv0).abs() < BOUNDARY_EPS || (tv1 - b.y).abs() < BOUNDARY_EPS)
        {
            return false;
        }
        // Endpoints on the target's own open boundaries are F3b territory.
        if !target.is_closed_in_u()
            && ((a.x - au0).abs() < BOUNDARY_EPS || (au1 - a.x).abs() < BOUNDARY_EPS)
        {
            return false;
        }
        if !target.is_closed_in_v()
            && ((a.y - av0).abs() < BOUNDARY_EPS || (av1 - a.y).abs() < BOUNDARY_EPS)
        {
            return false;
        }
    }
    true
}

/// Chains accepted open segments into closed [`CutChain`]s.
///
/// Deterministic: chains start at the earliest unused segment in the input
/// order (the SSI extraction iterates tool faces then target faces, so the
/// input order is stable) and extend from that segment's natural direction.
///
/// # Errors
///
/// Returns the verbatim open-branch typed error when an endpoint finds no
/// partner (partial cut) or the partner lies on the same tool face; a typed
/// error when a junction is ambiguous, a chain mixes target faces, or a chain
/// crosses one tool face twice.
pub(crate) fn chain_open_segments(
    segments: Vec<CutLoop>,
    tool_faces: &[(FaceId, NurbsSurface)],
) -> Result<Vec<CutChain>> {
    use crate::error::OperationError;

    let n = segments.len();
    let mut used = vec![false; n];
    let mut chains = Vec::new();

    for start in 0..n {
        if used[start] {
            continue;
        }
        used[start] = true;
        let mut chain: Vec<CutLoop> = vec![segments[start].clone()];
        let chain_start = first_point(&chain[0]);

        loop {
            let tail = chain.last().unwrap_or_else(|| unreachable!());
            let current_end = last_point(tail);

            // Closure: back at the chain's first sample.
            if chain.len() >= 2 && (current_end - chain_start).norm() <= JUNCTION_TOLERANCE {
                break;
            }

            // Find the unique unused partner endpoint at the current end.
            let mut candidates: Vec<(usize, bool)> = Vec::new();
            for (idx, seg) in segments.iter().enumerate() {
                if used[idx] {
                    continue;
                }
                if (first_point(seg) - current_end).norm() <= JUNCTION_TOLERANCE {
                    candidates.push((idx, false));
                }
                if (last_point(seg) - current_end).norm() <= JUNCTION_TOLERANCE {
                    candidates.push((idx, true));
                }
            }
            match candidates.as_slice() {
                [] => return Err(open_branch_error()),
                [(idx, flip)] => {
                    let mut next = segments[*idx].clone();
                    if *flip {
                        reverse_segment(&mut next);
                    }
                    if next.tool_face == tail.tool_face {
                        // A kink junction joins two ADJACENT tool faces.
                        return Err(open_branch_error());
                    }
                    if next.target_face != tail.target_face {
                        return Err(OperationError::Failed(
                            "chained cut loop crosses target faces (unsupported \
                             until general boolean face splitting)"
                                .into(),
                        )
                        .into());
                    }
                    used[*idx] = true;
                    chain.push(next);
                }
                _ => {
                    return Err(OperationError::Failed(
                        "ambiguous kink-edge junction: multiple open branch \
                         endpoints coincide at one tool kink crossing"
                            .into(),
                    )
                    .into());
                }
            }
        }

        // Each tool face contributes at most one segment per chained loop.
        for i in 0..chain.len() {
            for j in (i + 1)..chain.len() {
                if chain[i].tool_face == chain[j].tool_face {
                    return Err(OperationError::Failed(
                        "chained cut loop crosses one tool side face twice \
                         (unsupported)"
                            .into(),
                    )
                    .into());
                }
            }
        }

        let mut cut_chain = CutChain {
            target_face: chain[0].target_face,
            segments: chain,
        };
        weld_chain(&mut cut_chain, tool_faces)?;
        chains.push(cut_chain);
    }
    Ok(chains)
}

/// Welds every chain junction so adjacent segments share EXACT geometry: the
/// outgoing segment's last sample and the incoming segment's first sample get
/// the same 3D point, the same target UV, and tool UVs pinned exactly on
/// their respective `u` bounds at a common `v`.
///
/// The common 3D point is evaluated on the outgoing face at its pinned bound,
/// which lies on the shared kink edge; for extruded tools adjacent side faces
/// parameterize the kink identically in `v`, which is verified (within
/// [`JUNCTION_TOLERANCE`]) rather than assumed.
fn weld_chain(chain: &mut CutChain, tool_faces: &[(FaceId, NurbsSurface)]) -> Result<()> {
    use crate::error::OperationError;
    use crate::math::Point2;

    let surface_of = |face: FaceId| -> Option<&NurbsSurface> {
        tool_faces
            .iter()
            .find(|(id, _)| *id == face)
            .map(|(_, s)| s)
    };

    let n = chain.segments.len();
    for j in 0..n {
        let prev = (j + n - 1) % n;
        let (prev_face, prev_uv) = {
            let seg = &chain.segments[prev];
            (
                seg.tool_face,
                *seg.branch.uv_b.last().unwrap_or_else(|| unreachable!()),
            )
        };
        let (next_face, next_uv) = {
            let seg = &chain.segments[j];
            (
                seg.tool_face,
                *seg.branch.uv_b.first().unwrap_or_else(|| unreachable!()),
            )
        };
        let (Some(prev_surf), Some(next_surf)) = (surface_of(prev_face), surface_of(next_face))
        else {
            return Err(OperationError::Failed(
                "chained cut loop references an unknown tool face".into(),
            )
            .into());
        };

        let bound_prev = nearest_u_bound(prev_surf, prev_uv.x);
        let bound_next = nearest_u_bound(next_surf, next_uv.x);
        let v_weld = prev_uv.y;
        let point = prev_surf.point_at(bound_prev, v_weld)?;
        let next_point = next_surf.point_at(bound_next, v_weld)?;
        if (next_point - point).norm() > JUNCTION_TOLERANCE {
            return Err(OperationError::Failed(
                "adjacent tool side faces do not share the kink-edge \
                 parameterization at a chain junction"
                    .into(),
            )
            .into());
        }

        let target_uv = {
            let seg = &chain.segments[prev];
            *seg.branch.uv_a.last().unwrap_or_else(|| unreachable!())
        };

        {
            let seg = &mut chain.segments[prev];
            let last = seg.branch.points.len() - 1;
            seg.branch.points[last] = point;
            seg.branch.uv_b[last] = Point2::new(bound_prev, v_weld);
        }
        {
            let seg = &mut chain.segments[j];
            seg.branch.points[0] = point;
            seg.branch.uv_b[0] = Point2::new(bound_next, v_weld);
            seg.branch.uv_a[0] = target_uv;
        }
    }
    Ok(())
}

/// The nearest `u` domain bound to a value already known to sit within
/// [`BOUNDARY_EPS`] of one of them.
fn nearest_u_bound(surface: &NurbsSurface, u: f64) -> f64 {
    let ((u0, u1), _) = surface.parameter_domain();
    if (u - u0).abs() <= (u1 - u).abs() {
        u0
    } else {
        u1
    }
}

fn first_point(seg: &CutLoop) -> Point3 {
    *seg.branch.points.first().unwrap_or_else(|| unreachable!())
}

fn last_point(seg: &CutLoop) -> Point3 {
    *seg.branch.points.last().unwrap_or_else(|| unreachable!())
}

/// Reverses a segment's traversal direction in place (all three synchronized
/// traces).
fn reverse_segment(seg: &mut CutLoop) {
    seg.branch.points.reverse();
    seg.branch.uv_a.reverse();
    seg.branch.uv_b.reverse();
}
