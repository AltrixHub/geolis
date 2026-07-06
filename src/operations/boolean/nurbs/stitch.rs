//! Open-branch stitcher: chains open SSI branches across adjacent tool side
//! faces AND adjacent target faces into closed cut loops (F5 Phases B + C).
//!
//! A multi-face (box-like) tool has planar side faces meeting at kink edges;
//! each side face's intersection with a target face is an OPEN SSI branch
//! whose endpoints land exactly on the tool's kink edges. Symmetrically, a
//! cut that crosses a TARGET kink edge (a window sliding across a
//! segmented-prism joint) produces open branches whose endpoints land exactly
//! on the target face's own domain boundary (the marcher pins open-boundary
//! endpoints — see the SSI boundary-pinned termination). This module accepts
//! both endpoint kinds, chains the branches end-to-end into closed
//! [`CutChain`]s, and welds every junction so the punched hole rings, the
//! per-face band fragments, and the crossing edges all share exactly the same
//! junction geometry.
//!
//! Acceptance is structural, not tolerance-driven: an open-branch endpoint is
//! acceptable iff it sits on an open `u` boundary of the tool face (a tool
//! kink) while staying strictly inside the target face's domain, OR on an
//! open boundary of the target face (a target kink / ring boundary) while
//! staying strictly inside the tool face's domain — always within the
//! marcher's own [`BOUNDARY_EPS`], never a new tolerance. A corner endpoint
//! (pinned on both — a tool kink coinciding with a target boundary, e.g. a
//! cutter sill flush with the target's cap) is admissible as a
//! target-boundary terminal (F6 R3): the chain continues through the kink
//! when a partner branch exists and terminates on the target boundary
//! otherwise. Every endpoint must find exactly one partner endpoint at the
//! same 3D point; a junction must change the tool face (kink crossing, same
//! target) or the target face (target boundary crossing, same tool) —
//! changing both at once is ambiguous. Any violation keeps the pre-chaining
//! typed errors verbatim.

use crate::error::Result;
use crate::geometry::nurbs::{NurbsSurface, SurfaceIntersectionCurve, BOUNDARY_EPS};
use crate::math::Point3;
use crate::operations::boolean::nurbs::loops::open_branch_error;
use crate::topology::FaceId;

use super::loops::CutLoop;

/// 3D junction coincidence bound: the SSI marcher's point-acceptance bound
/// (`IntersectionOptions.tolerance.max(1e-7)` in the corrector's final
/// check). Two branch endpoints converged onto the same kink crossing agree
/// within this — no new tolerance is introduced, and Newton typically
/// converges orders of magnitude tighter.
pub(crate) const JUNCTION_TOLERANCE: f64 = 1e-7;

/// A cut loop chained from open SSI branches. `segments[i]` runs on one
/// (tool face × target face) pair; segment `i`'s last sample coincides
/// exactly (welded) with segment `i + 1`'s first sample (cyclic when
/// `closed`).
///
/// An OPEN chain (`closed == false`) is a cap-touching cut: both terminal
/// endpoints are pinned exactly on an open boundary of their target face —
/// the ring edge shared with a planar cap. The cut's circuit is closed by
/// the cap-plane closure edges the band assembly creates between the entry
/// and exit chains' matching terminals (F6 R2).
#[derive(Debug, Clone)]
pub(crate) struct CutChain {
    /// The chained per-face-pair segments, in cyclic (closed) or
    /// terminal-to-terminal (open) order.
    pub segments: Vec<CutLoop>,
    /// Whether the chain closes back on itself.
    pub closed: bool,
}

impl CutChain {
    /// The single target face every segment lies on, or `None` when the
    /// chain crosses target faces (the F3b splitting case).
    pub(crate) fn single_target_face(&self) -> Option<FaceId> {
        let first = self.segments.first()?.target_face;
        self.segments
            .iter()
            .all(|s| s.target_face == first)
            .then_some(first)
    }

    /// Whether the chain crosses more than one target face.
    pub(crate) fn crosses_target_faces(&self) -> bool {
        self.single_target_face().is_none() && !self.segments.is_empty()
    }

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

/// How one endpoint of an open SSI branch terminates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum EndpointKind {
    /// Pinned on an open `u` boundary of the TOOL face, strictly inside the
    /// target face — a tool kink-edge crossing (Phase B).
    ToolKink,
    /// Pinned on an open boundary of the TARGET face, strictly inside the
    /// tool face — a target boundary crossing (Phase C, face splitting).
    TargetBoundary,
}

/// Classifies one open-branch endpoint, or `None` when it is inadmissible
/// (interior endpoint = partial cut, or a degenerate corner pinned on both
/// the tool and the target boundary).
fn classify_endpoint(
    uv_b: crate::math::Point2,
    uv_a: crate::math::Point2,
    target: &NurbsSurface,
    tool: &NurbsSurface,
) -> Option<EndpointKind> {
    let ((tu0, tu1), (tv0, tv1)) = tool.parameter_domain();
    let ((au0, au1), (av0, av1)) = target.parameter_domain();

    let on_tool_u = !tool.is_closed_in_u()
        && ((uv_b.x - tu0).abs() < BOUNDARY_EPS || (tu1 - uv_b.x).abs() < BOUNDARY_EPS);
    let on_tool_v = !tool.is_closed_in_v()
        && ((uv_b.y - tv0).abs() < BOUNDARY_EPS || (tv1 - uv_b.y).abs() < BOUNDARY_EPS);
    let on_target_u = !target.is_closed_in_u()
        && ((uv_a.x - au0).abs() < BOUNDARY_EPS || (au1 - uv_a.x).abs() < BOUNDARY_EPS);
    let on_target_v = !target.is_closed_in_v()
        && ((uv_a.y - av0).abs() < BOUNDARY_EPS || (av1 - uv_a.y).abs() < BOUNDARY_EPS);

    let on_target = on_target_u || on_target_v;
    // A tool `v` boundary is a cap ring, never a kink; an endpoint there is
    // degenerate unless it is a clean target-boundary crossing.
    let on_tool_kink = on_tool_u && !on_tool_v;

    match (on_tool_kink, on_target) {
        (true, false) => Some(EndpointKind::ToolKink),
        (false, true) if !on_tool_u && !on_tool_v => Some(EndpointKind::TargetBoundary),
        // A corner endpoint — tool kink coinciding with a target boundary
        // (e.g. a cutter sill flush with the target's bottom cap, F6 R3).
        // Admissible as a target-boundary terminal: chain extension still
        // continues through the kink when a partner branch exists, and
        // terminates on the target boundary when the flush face's
        // tangential branches were dropped by the loop extraction.
        (true, true) => Some(EndpointKind::TargetBoundary),
        _ => None,
    }
}

/// Whether an open SSI branch qualifies for chaining: BOTH endpoints must be
/// admissible ([`classify_endpoint`]) — each either on a tool kink edge or on
/// a target face boundary.
pub(crate) fn open_branch_admissible(
    branch: &SurfaceIntersectionCurve,
    target: &NurbsSurface,
    tool: &NurbsSurface,
) -> bool {
    let ends = [
        (branch.uv_b.first(), branch.uv_a.first()),
        (branch.uv_b.last(), branch.uv_a.last()),
    ];
    for (uv_b, uv_a) in ends {
        let (Some(&b), Some(&a)) = (uv_b, uv_a) else {
            return false;
        };
        if classify_endpoint(b, a, target, tool).is_none() {
            return false;
        }
    }
    true
}

/// Chains accepted open segments into [`CutChain`]s.
///
/// Deterministic: chains start at the earliest unused segment in the input
/// order (the SSI extraction iterates tool faces then target faces, so the
/// input order is stable) and extend from that segment's natural direction.
/// A chain that returns to its first sample is CLOSED. A chain whose forward
/// extension exhausts on a target-boundary endpoint is extended backward
/// from its head; when both terminals rest on target-boundary crossings the
/// chain is accepted OPEN (a cap-touching cut), oriented so its
/// lexicographically smaller terminal comes first, and its terminals are
/// pinned exactly onto their target boundary bounds.
///
/// # Errors
///
/// Returns the verbatim open-branch typed error when an exhausted endpoint
/// is not a clean target-boundary crossing (partial cut) or a junction
/// partner changes neither the tool nor the target face; a typed error when
/// a junction is ambiguous (multiple partners, or both faces change at
/// once) or a chain crosses one tool face in two separate runs.
pub(crate) fn chain_open_segments(
    segments: &[CutLoop],
    target_faces: &[(FaceId, NurbsSurface)],
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

        let closed = extend_chain_forward(&mut chain, segments, &mut used)?;
        if !closed {
            extend_chain_backward(&mut chain, segments, &mut used)?;
            // Both exhausted terminals must be clean target-boundary
            // crossings (the ring shared with a planar cap); anything else
            // keeps the verbatim partial-cut error.
            for (uv_b, uv_a, target, tool) in [
                terminal(&chain, false, target_faces, tool_faces)?,
                terminal(&chain, true, target_faces, tool_faces)?,
            ] {
                if classify_endpoint(uv_b, uv_a, &target, &tool)
                    != Some(EndpointKind::TargetBoundary)
                {
                    return Err(open_branch_error());
                }
            }
            normalize_open_chain(&mut chain);
        }

        // Each tool face contributes at most ONE contiguous run of segments
        // per chained loop (cyclic for closed chains): a loop entering the
        // same tool face in two separate runs is out of scope.
        if run_count_exceeds_one(&chain, closed, |s| s.tool_face) {
            return Err(OperationError::Failed(
                "chained cut loop crosses one tool side face twice \
                 (unsupported)"
                    .into(),
            )
            .into());
        }

        let mut cut_chain = CutChain {
            segments: chain,
            closed,
        };
        weld_chain(&mut cut_chain, target_faces, tool_faces)?;
        if !closed {
            pin_open_terminals(&mut cut_chain, target_faces)?;
        }
        chains.push(cut_chain);
    }
    Ok(chains)
}

/// Extends a chain forward from its tail until it closes (`Ok(true)`) or no
/// partner endpoint exists (`Ok(false)`).
fn extend_chain_forward(
    chain: &mut Vec<CutLoop>,
    segments: &[CutLoop],
    used: &mut [bool],
) -> Result<bool> {
    let chain_start = first_point(&chain[0]);
    loop {
        let tail = chain.last().unwrap_or_else(|| unreachable!());
        let current_end = last_point(tail);

        // Closure: back at the chain's first sample.
        if chain.len() >= 2 && (current_end - chain_start).norm() <= JUNCTION_TOLERANCE {
            return Ok(true);
        }

        match find_partner(segments, used, current_end)? {
            None => return Ok(false),
            Some((idx, flip)) => {
                let mut next = segments[idx].clone();
                if flip {
                    reverse_segment(&mut next);
                }
                check_junction(tail, &next)?;
                used[idx] = true;
                chain.push(next);
            }
        }
    }
}

/// Extends a chain backward from its head until no partner endpoint exists.
fn extend_chain_backward(
    chain: &mut Vec<CutLoop>,
    segments: &[CutLoop],
    used: &mut [bool],
) -> Result<()> {
    loop {
        let head = chain.first().unwrap_or_else(|| unreachable!());
        let current_start = first_point(head);
        match find_partner(segments, used, current_start)? {
            None => return Ok(()),
            Some((idx, flip)) => {
                let mut prev = segments[idx].clone();
                // A backward partner matched via its FIRST point must be
                // flipped so its last sample meets the chain head.
                if !flip {
                    reverse_segment(&mut prev);
                }
                check_junction(&prev, head)?;
                used[idx] = true;
                chain.insert(0, prev);
            }
        }
    }
}

/// Finds the unique unused partner endpoint at `point`. `Ok(None)` when no
/// endpoint coincides; the `bool` is `true` when the partner matched via its
/// LAST sample (it must be reversed for forward extension).
///
/// # Errors
///
/// A typed error when multiple endpoints coincide (ambiguous junction).
fn find_partner(
    segments: &[CutLoop],
    used: &[bool],
    point: Point3,
) -> Result<Option<(usize, bool)>> {
    use crate::error::OperationError;

    let mut candidates: Vec<(usize, bool)> = Vec::new();
    for (idx, seg) in segments.iter().enumerate() {
        if used[idx] {
            continue;
        }
        if (first_point(seg) - point).norm() <= JUNCTION_TOLERANCE {
            candidates.push((idx, false));
        }
        if (last_point(seg) - point).norm() <= JUNCTION_TOLERANCE {
            candidates.push((idx, true));
        }
    }
    match candidates.as_slice() {
        [] => Ok(None),
        [unique] => Ok(Some(*unique)),
        _ => Err(OperationError::Failed(
            "ambiguous kink-edge junction: multiple open branch \
             endpoints coincide at one tool kink crossing"
                .into(),
        )
        .into()),
    }
}

/// Validates a chain junction: it must cross a tool kink (adjacent tool
/// faces) or a target boundary (adjacent target faces) — never neither,
/// never both.
fn check_junction(from: &CutLoop, to: &CutLoop) -> Result<()> {
    use crate::error::OperationError;

    let tool_changes = to.tool_face != from.tool_face;
    let target_changes = to.target_face != from.target_face;
    match (tool_changes, target_changes) {
        (false, false) => Err(open_branch_error()),
        (true, true) => Err(OperationError::Failed(
            "ambiguous chain junction: tool face and target \
             face change at the same crossing point"
                .into(),
        )
        .into()),
        _ => Ok(()),
    }
}

/// The (tool UV, target UV, target surface, tool surface) of an open
/// chain's head (`tail == false`) or tail terminal.
fn terminal(
    chain: &[CutLoop],
    tail: bool,
    target_faces: &[(FaceId, NurbsSurface)],
    tool_faces: &[(FaceId, NurbsSurface)],
) -> Result<(
    crate::math::Point2,
    crate::math::Point2,
    NurbsSurface,
    NurbsSurface,
)> {
    use crate::error::OperationError;

    let seg = if tail {
        chain.last().unwrap_or_else(|| unreachable!())
    } else {
        &chain[0]
    };
    let (uv_b, uv_a) = if tail {
        (
            *seg.branch.uv_b.last().unwrap_or_else(|| unreachable!()),
            *seg.branch.uv_a.last().unwrap_or_else(|| unreachable!()),
        )
    } else {
        (seg.branch.uv_b[0], seg.branch.uv_a[0])
    };
    let target = target_faces
        .iter()
        .find(|(id, _)| *id == seg.target_face)
        .map(|(_, s)| s.clone())
        .ok_or_else(|| {
            OperationError::Failed("chained cut loop references an unknown target face".into())
        })?;
    let tool = tool_faces
        .iter()
        .find(|(id, _)| *id == seg.tool_face)
        .map(|(_, s)| s.clone())
        .ok_or_else(|| {
            OperationError::Failed("chained cut loop references an unknown tool face".into())
        })?;
    Ok((uv_b, uv_a, target, tool))
}

/// Orients an open chain deterministically: the lexicographically smaller
/// 3D terminal comes first, so the entry and exit chains of one
/// cap-touching cut traverse their shared tool faces in the same direction
/// and their terminals pair index-to-index.
///
/// The comparison treats coordinates within [`JUNCTION_TOLERANCE`] as equal
/// before moving to the next coordinate: the two terminals of one chain
/// nominally share coordinates that differ only by evaluation noise (a jamb
/// chain's terminals share x and y exactly up to the marcher's residual),
/// and a raw float comparison would flip the orientation on that noise.
fn normalize_open_chain(chain: &mut [CutLoop]) {
    let head = first_point(&chain[0]);
    let tail = last_point(chain.last().unwrap_or_else(|| unreachable!()));
    if approx_lex_less(tail, head) {
        chain.reverse();
        for seg in chain.iter_mut() {
            reverse_segment(seg);
        }
    }
}

/// Noise-robust lexicographic 3D comparison: coordinates within
/// [`JUNCTION_TOLERANCE`] compare equal, then the next coordinate decides.
fn approx_lex_less(a: Point3, b: Point3) -> bool {
    for (x, y) in [(a.x, b.x), (a.y, b.y), (a.z, b.z)] {
        if (x - y).abs() > JUNCTION_TOLERANCE {
            return x < y;
        }
    }
    false
}

/// Pins an open chain's two terminals exactly onto their target boundary
/// bounds: the target UV is pinned per [`pin_to_boundary`] and the 3D point
/// re-evaluated on the target surface at the pinned UV — exactly on the
/// boundary ring curve shared with the planar cap.
fn pin_open_terminals(chain: &mut CutChain, target_faces: &[(FaceId, NurbsSurface)]) -> Result<()> {
    use crate::error::OperationError;

    let surface_of = |face: FaceId| -> Result<NurbsSurface> {
        target_faces
            .iter()
            .find(|(id, _)| *id == face)
            .map(|(_, s)| s.clone())
            .ok_or_else(|| {
                OperationError::Failed("chained cut loop references an unknown target face".into())
                    .into()
            })
    };

    {
        let seg = &mut chain.segments[0];
        let surf = surface_of(seg.target_face)?;
        let pin = pin_to_boundary(&surf, seg.branch.uv_a[0]);
        seg.branch.points[0] = surf.point_at(pin.x, pin.y)?;
        seg.branch.uv_a[0] = pin;
    }
    {
        let seg = chain.segments.last_mut().unwrap_or_else(|| unreachable!());
        let surf = surface_of(seg.target_face)?;
        let last = seg.branch.points.len() - 1;
        let pin = pin_to_boundary(&surf, seg.branch.uv_a[last]);
        seg.branch.points[last] = surf.point_at(pin.x, pin.y)?;
        seg.branch.uv_a[last] = pin;
    }
    Ok(())
}

/// Whether some key (per-segment face) appears in more than one contiguous
/// run of the chain (cyclically contiguous for closed chains).
fn run_count_exceeds_one(
    chain: &[CutLoop],
    closed: bool,
    key: impl Fn(&CutLoop) -> FaceId,
) -> bool {
    use std::collections::HashMap;
    let n = chain.len();
    let mut runs: HashMap<FaceId, usize> = HashMap::new();
    for i in 0..n {
        let k = key(&chain[i]);
        let starts_run = if i == 0 {
            !closed || n == 1 || k != key(&chain[n - 1])
        } else {
            k != key(&chain[i - 1])
        };
        if starts_run {
            *runs.entry(k).or_insert(0) += 1;
        }
    }
    // A closed single-face chain (n >= 2, all same face) has one cyclic run.
    if closed && runs.is_empty() {
        return false;
    }
    runs.values().any(|&c| c > 1)
}

/// Welds every chain junction so adjacent segments share EXACT geometry.
///
/// Tool-kink junction (tool face changes, target face stays): the outgoing
/// and incoming tool UVs are pinned exactly on their `u` bounds at a common
/// `v`; the common 3D point is evaluated on the outgoing tool face at its
/// pinned bound (on the shared kink edge). Adjacent extruded tool faces
/// parameterize the kink identically in `v`, which is verified (within
/// [`JUNCTION_TOLERANCE`]) rather than assumed.
///
/// Target-boundary junction (target face changes, tool face stays): the
/// mirror image — target UVs are pinned exactly on their target-domain
/// bounds at the outgoing segment's free coordinate; the common 3D point is
/// evaluated on the outgoing target face at its pinned bound (on the shared
/// target kink / boundary edge), and the incoming target face's
/// parameterization of that boundary is verified the same way. The tool UV
/// is carried over unchanged (same tool face).
fn weld_chain(
    chain: &mut CutChain,
    target_faces: &[(FaceId, NurbsSurface)],
    tool_faces: &[(FaceId, NurbsSurface)],
) -> Result<()> {
    use crate::error::OperationError;

    let surface_in = |faces: &[(FaceId, NurbsSurface)], face: FaceId| -> Option<NurbsSurface> {
        faces
            .iter()
            .find(|(id, _)| *id == face)
            .map(|(_, s)| s.clone())
    };

    let n = chain.segments.len();
    // Open chains have no wrap-around junction: their terminals are pinned
    // on their target boundary bounds instead ([`pin_open_terminals`]).
    let first_junction = usize::from(!chain.closed);
    for j in first_junction..n {
        let prev = (j + n - 1) % n;
        let (prev_tool, prev_target, prev_uv_b, prev_uv_a) = {
            let seg = &chain.segments[prev];
            (
                seg.tool_face,
                seg.target_face,
                *seg.branch.uv_b.last().unwrap_or_else(|| unreachable!()),
                *seg.branch.uv_a.last().unwrap_or_else(|| unreachable!()),
            )
        };
        let (next_tool, next_target) = {
            let seg = &chain.segments[j];
            (seg.tool_face, seg.target_face)
        };

        if prev_tool == next_tool {
            // Target-boundary junction (same tool face on both sides).
            let (Some(prev_surf), Some(next_surf)) = (
                surface_in(target_faces, prev_target),
                surface_in(target_faces, next_target),
            ) else {
                return Err(OperationError::Failed(
                    "chained cut loop references an unknown target face".into(),
                )
                .into());
            };
            weld_target_junction(chain, prev, j, &prev_surf, &next_surf, prev_uv_b)?;
        } else {
            // Tool-kink junction (same target face on both sides).
            let (Some(prev_surf), Some(next_surf)) = (
                surface_in(tool_faces, prev_tool),
                surface_in(tool_faces, next_tool),
            ) else {
                return Err(OperationError::Failed(
                    "chained cut loop references an unknown tool face".into(),
                )
                .into());
            };
            weld_tool_junction(chain, prev, j, &prev_surf, &next_surf, prev_uv_a)?;
        }
    }
    Ok(())
}

/// Welds a tool-kink junction: tool UVs pinned exactly on their `u` bounds
/// at a common `v`; the shared 3D point evaluated on the outgoing tool face.
fn weld_tool_junction(
    chain: &mut CutChain,
    prev: usize,
    next: usize,
    prev_surf: &NurbsSurface,
    next_surf: &NurbsSurface,
    prev_uv_a: crate::math::Point2,
) -> Result<()> {
    use crate::error::OperationError;
    use crate::math::Point2;

    let prev_uv_b = *chain.segments[prev]
        .branch
        .uv_b
        .last()
        .unwrap_or_else(|| unreachable!());
    let next_uv_b = *chain.segments[next]
        .branch
        .uv_b
        .first()
        .unwrap_or_else(|| unreachable!());
    let bound_prev = nearest_u_bound(prev_surf, prev_uv_b.x);
    let bound_next = nearest_u_bound(next_surf, next_uv_b.x);
    let v_weld = prev_uv_b.y;
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
    {
        let seg = &mut chain.segments[prev];
        let last = seg.branch.points.len() - 1;
        seg.branch.points[last] = point;
        seg.branch.uv_b[last] = Point2::new(bound_prev, v_weld);
    }
    {
        let seg = &mut chain.segments[next];
        seg.branch.points[0] = point;
        seg.branch.uv_b[0] = Point2::new(bound_next, v_weld);
        seg.branch.uv_a[0] = prev_uv_a;
    }
    Ok(())
}

/// Welds a target-boundary junction: target UVs pinned exactly on their
/// domain bounds; the shared 3D point evaluated on the outgoing target face
/// and verified against the incoming face's parameterization.
fn weld_target_junction(
    chain: &mut CutChain,
    prev: usize,
    next: usize,
    prev_surf: &NurbsSurface,
    next_surf: &NurbsSurface,
    prev_uv_b: crate::math::Point2,
) -> Result<()> {
    use crate::error::OperationError;

    let prev_uv_a = *chain.segments[prev]
        .branch
        .uv_a
        .last()
        .unwrap_or_else(|| unreachable!());
    let next_uv_a = *chain.segments[next]
        .branch
        .uv_a
        .first()
        .unwrap_or_else(|| unreachable!());
    let prev_pin = pin_to_boundary(prev_surf, prev_uv_a);
    let next_pin = pin_to_boundary(next_surf, next_uv_a);
    let point = prev_surf.point_at(prev_pin.x, prev_pin.y)?;
    let next_point = next_surf.point_at(next_pin.x, next_pin.y)?;
    if (next_point - point).norm() > JUNCTION_TOLERANCE {
        return Err(OperationError::Failed(
            "adjacent target faces do not share the boundary-edge \
             parameterization at a chain junction"
                .into(),
        )
        .into());
    }
    {
        let seg = &mut chain.segments[prev];
        let last = seg.branch.points.len() - 1;
        seg.branch.points[last] = point;
        seg.branch.uv_a[last] = prev_pin;
    }
    {
        let seg = &mut chain.segments[next];
        seg.branch.points[0] = point;
        seg.branch.uv_a[0] = next_pin;
        seg.branch.uv_b[0] = prev_uv_b;
    }
    Ok(())
}

/// Pins a UV point already known to sit within [`BOUNDARY_EPS`] of an open
/// domain boundary exactly onto that boundary (the nearer coordinate wins
/// when both are near a bound — the caller rejects genuine corners earlier).
fn pin_to_boundary(surface: &NurbsSurface, uv: crate::math::Point2) -> crate::math::Point2 {
    use crate::math::Point2;
    let ((u0, u1), (v0, v1)) = surface.parameter_domain();
    let du = (uv.x - u0).abs().min((u1 - uv.x).abs());
    let dv = (uv.y - v0).abs().min((v1 - uv.y).abs());
    let u_open = !surface.is_closed_in_u();
    let v_open = !surface.is_closed_in_v();
    if u_open && (!v_open || du <= dv) {
        let u = if (uv.x - u0).abs() <= (u1 - uv.x).abs() {
            u0
        } else {
            u1
        };
        Point2::new(u, uv.y)
    } else {
        let v = if (uv.y - v0).abs() <= (v1 - uv.y).abs() {
            v0
        } else {
            v1
        };
        Point2::new(uv.x, v)
    }
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
