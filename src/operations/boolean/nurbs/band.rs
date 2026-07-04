//! Tool band (hole-wall) face construction for the through-cut subtract.
//!
//! ## Shipped path: generic stitched-CDT
//!
//! The plan offered a fallback ladder for the band's UV topology (generic
//! stitched-CDT first, dedicated quad-strip only if that proved infeasible).
//! The **generic stitched-CDT path shipped**: the band between the entry and
//! exit loops on the tool side surface is represented as an ordinary trimmed
//! NURBS face whose `FaceTrim.outer` is a single simple polygon in the tool's
//! unrolled UV rectangle, and the existing constrained-Delaunay trimmed
//! tessellator (P5) meshes it directly. No dedicated band tessellation was
//! needed.
//!
//! ### Why it works
//!
//! Each loop's `uv_b` trace is a single-valued graph `v = f(u)` over the tool's
//! `u` domain. The SSI marcher wraps the tool's periodic `u` direction, so the
//! loop arrives genuinely `closed` with its tool `u` kept wrapped into
//! `[u0, u1]` and with **exact seam samples** (the crossing point emitted at
//! both `u0` and `u1`), so each trace spans the **full** `u` domain at a
//! roughly constant `v`. The entry loop sits at a lower mean `v` than the exit
//! loop (the loops are pre-sorted by mean `v` in [`super::loops`]). Stitching
//!
//! ```text
//!   entry trace  (u increasing)
//!   -> exit trace (u decreasing)
//!   -> close
//! ```
//!
//! yields a ribbon polygon that is simple (non-self-intersecting) in the
//! unrolled rectangle, so the generic trimmed CDT meshes it without a seam cut.
//! Because the traces reach `u0` and `u1` exactly, the ribbon's left (`u0`)
//! and right (`u1`) closing edges land on the same seam azimuth and coincide in
//! 3D, covering the seam wedge. The two rings (entry/exit) share the exact seam
//! samples with the punched target faces, so the hole rim tessellates
//! conformally with no slit at the seam. If a marcher seam sample did not
//! converge, the ribbon degrades to the marched span (a sub-step gap at the
//! seam) — the honest fallback.
//!
//! ### Orientation
//!
//! Subtract pushes the band normals INTO the hole, so the band face is built
//! with `same_sense = false`.

use std::collections::HashMap;

use crate::error::{OperationError, Result};
use crate::geometry::nurbs::{KnotVector, NurbsCurve2D, NurbsCurve3D, NurbsSurface};
use crate::math::Point2;
use crate::topology::{
    EdgeCurve, EdgeData, EdgeId, FaceData, FaceId, FaceSurface, FaceTrim, OrientedEdge,
    TopologyStore, TrimLoop, VertexId, WireData, WireId,
};

use super::loops::CutLoop;
use super::punch::ChainRing;
use super::stitch::CutChain;

/// The two hole-ring wires shared with the punched target faces for one tool
/// side face: the entry ring (lower mean v) and the exit ring (upper mean v).
///
/// These are the exact [`WireId`]s returned by [`super::punch::punch_loop`] for
/// the same tool face's two loops, so the band face shares its boundary edges
/// with the punched target faces (correct `BRep` adjacency).
#[derive(Debug, Clone, Copy)]
pub(crate) struct BandRingWires {
    /// Entry ring wire (matches `cut.loops[0]`, the lower-v loop).
    pub entry: WireId,
    /// Exit ring wire (matches `cut.loops[1]`, the upper-v loop).
    pub exit: WireId,
}

/// Builds the band (hole-wall) face for one tool side face from its two cut
/// loops, and returns the new face's id.
///
/// `rings` carries the two hole-ring wires already created by the punch step for
/// this tool face's loops; the band face reuses them as its boundary so it
/// shares edges/wires with the punched target faces instead of fabricating a new
/// full-surface boundary.
///
/// # Errors
///
/// Returns an error if the tool face is not a NURBS face or the stitched band
/// polygon degenerates (fewer than 3 distinct UV points).
pub(crate) fn build_band_face(
    store: &mut TopologyStore,
    tool_face: FaceId,
    loops: &[CutLoop; 2],
    rings: BandRingWires,
) -> Result<FaceId> {
    // Subtract: band normals point INTO the hole (`same_sense = false`).
    build_band_face_oriented(store, tool_face, loops, rings, false)
}

/// Builds the band (hole-wall / plug-wall) face for one tool side face with an
/// explicit `same_sense`.
///
/// The band region (the tool side surface strip between the entry and exit
/// loops) is identical for the subtract and intersect through-cuts; only the
/// normal orientation differs. Subtract points the band normals into the hole
/// (`same_sense = false`); intersect points them outward from the kept plug
/// (`same_sense = true`).
///
/// # Errors
///
/// Returns an error if the tool face is not a NURBS face or the stitched band
/// polygon degenerates (fewer than 3 distinct UV points).
pub(crate) fn build_band_face_oriented(
    store: &mut TopologyStore,
    tool_face: FaceId,
    loops: &[CutLoop; 2],
    rings: BandRingWires,
    same_sense: bool,
) -> Result<FaceId> {
    let surface = match &store.face(tool_face)?.surface {
        FaceSurface::Nurbs(s) => s.clone(),
        _ => {
            return Err(OperationError::Failed(
                "through-cut band requires a NURBS tool side face".into(),
            )
            .into())
        }
    };

    let entry = clamp_trace(&loops[0].branch.uv_b, &surface);
    let exit = clamp_trace(&loops[1].branch.uv_b, &surface);

    let outer = stitch_band_loop(&entry, &exit)?;
    let trim = FaceTrim::new(outer, Vec::new());

    // The band's real boundary is the two SSI rings (entry + exit). Share the
    // exact wires the punch step attached to the target faces so the band face
    // has correct BRep adjacency (no fabricated full-surface seam wire).
    Ok(store.add_face(FaceData {
        surface: FaceSurface::Nurbs(surface),
        outer_wire: rings.entry,
        inner_wires: vec![rings.exit],
        same_sense,
        trim: Some(trim),
        pcurves: Vec::new(),
    }))
}

/// Builds the pocket band face: the tool side strip from the entry loop down
/// to the buried ring (`v = v_boundary` across the full `u` domain).
///
/// `entry_ring` is the punched entry hole's wire (shared with the entry
/// face); `buried_ring` is the tool's own shared ring wire at the buried end
/// (shared with the pocket floor). `buried_uv` carries the ring's UV samples
/// at the cache-identical chord parameters (see
/// [`super::pocket::buried_ring_uv`]).
///
/// # Errors
///
/// Returns an error if the tool face is not NURBS or the ribbon degenerates.
pub(crate) fn build_pocket_band_face(
    store: &mut TopologyStore,
    tool_face: FaceId,
    entry_loop: &CutLoop,
    buried_uv: &[Point2],
    entry_ring: WireId,
    buried_ring: WireId,
) -> Result<FaceId> {
    let surface = match &store.face(tool_face)?.surface {
        FaceSurface::Nurbs(s) => s.clone(),
        _ => {
            return Err(OperationError::Failed(
                "pocket band requires a NURBS tool side face".into(),
            )
            .into())
        }
    };

    let entry = clamp_trace(&entry_loop.branch.uv_b, &surface);
    let outer = stitch_band_loop(&entry, buried_uv)?;
    let trim = FaceTrim::new(outer, Vec::new());

    // Pocket band normals point INTO the cavity (`same_sense = false`), like
    // the through-cut hole wall.
    Ok(store.add_face(FaceData {
        surface: FaceSurface::Nurbs(surface),
        outer_wire: entry_ring,
        inner_wires: vec![buried_ring],
        same_sense: false,
        trim: Some(trim),
        pcurves: Vec::new(),
    }))
}

/// One band fragment of a multi-face through cut: the hole-wall face built on
/// one tool side face, trimmed to that face's chain segments.
#[derive(Debug, Clone, Copy)]
pub(crate) struct BandFragment {
    /// The tool side face this fragment lies on (its name seeds the band
    /// fragment's persistent name).
    pub tool_face: FaceId,
    /// The new band fragment face.
    pub face: FaceId,
}

/// One contiguous cyclic run of chain segments on a single tool side face.
#[derive(Debug, Clone)]
struct ToolRun {
    face: FaceId,
    /// Chain-order segment indices (cyclically contiguous).
    indices: Vec<usize>,
}

impl ToolRun {
    fn first(&self) -> usize {
        self.indices[0]
    }
    fn last(&self) -> usize {
        *self.indices.last().unwrap_or_else(|| unreachable!())
    }
}

/// Groups a chain's segments into cyclic contiguous runs per tool face
/// (the stitcher guarantees at most one run per face).
fn tool_runs(chain: &CutChain) -> Result<Vec<ToolRun>> {
    let n = chain.segments.len();
    let Some(first_change) =
        (0..n).find(|&i| chain.segments[i].tool_face != chain.segments[(i + n - 1) % n].tool_face)
    else {
        return Err(OperationError::Failed(
            "multi-face band requires a chain crossing tool side faces".into(),
        )
        .into());
    };
    let mut runs: Vec<ToolRun> = Vec::new();
    let mut i = 0usize;
    while i < n {
        let start = (first_change + i) % n;
        let face = chain.segments[start].tool_face;
        let mut len = 1usize;
        while len < n && chain.segments[(start + len) % n].tool_face == face {
            len += 1;
        }
        runs.push(ToolRun {
            face,
            indices: (0..len).map(|k| (start + k) % n).collect(),
        });
        i += len;
    }
    Ok(runs)
}

/// Groups an OPEN chain's segments into linear contiguous runs per tool
/// face, from the chain head (a single-face open chain yields one run).
fn open_tool_runs(chain: &CutChain) -> Vec<ToolRun> {
    let n = chain.segments.len();
    let mut runs: Vec<ToolRun> = Vec::new();
    let mut start = 0usize;
    while start < n {
        let face = chain.segments[start].tool_face;
        let mut len = 1usize;
        while start + len < n && chain.segments[start + len].tool_face == face {
            len += 1;
        }
        runs.push(ToolRun {
            face,
            indices: (start..start + len).collect(),
        });
        start += len;
    }
    runs
}

/// Builds the band fragments of a multi-face through cut: ONE face per tool
/// side face crossed by the chained loops, each trimmed to that face's chain
/// segments (possibly several, when the loop crosses target face boundaries
/// on that tool face), sharing its entry / exit ring edges with the punched
/// or split target faces and NEW kink-crossing edges with the adjacent
/// fragments (F2 shared-edge topology).
///
/// The kink-crossing edges are straight segments between the entry and exit
/// junction points — exact for extruded tools, whose kink edges are straight
/// lines; they map to the fragments' `u = const` trim closings, so adjacent
/// fragments emit identical rim vertices.
///
/// # Errors
///
/// Returns an error when the chains disagree on the crossed tool faces or
/// kink junctions, a tool face is not NURBS, or a fragment polygon
/// degenerates.
pub(crate) fn build_band_fragments(
    store: &mut TopologyStore,
    entry: &CutChain,
    exit: &CutChain,
    entry_ring: &ChainRing,
    exit_ring: &ChainRing,
) -> Result<Vec<BandFragment>> {
    let entry_runs = tool_runs(entry)?;
    let exit_runs = tool_runs(exit)?;
    let n_runs = entry_runs.len();
    if n_runs < 3 || exit_runs.len() != n_runs {
        return Err(OperationError::Failed(
            "multi-face through cut requires entry and exit chains crossing \
             the same three or more tool side faces"
                .into(),
        )
        .into());
    }

    // Run index per face within the ENTRY chain: normalizes junction pairs.
    let entry_index: HashMap<FaceId, usize> = entry_runs
        .iter()
        .enumerate()
        .map(|(i, r)| (r.face, i))
        .collect();
    let exit_run_of = |face: FaceId| -> Result<&ToolRun> {
        exit_runs.iter().find(|r| r.face == face).ok_or_else(|| {
            OperationError::Failed(
                "entry and exit chained loops cross different tool side \
                 faces"
                    .into(),
            )
            .into()
        })
    };
    let kink_edges = build_run_kink_edges(
        store,
        entry,
        exit,
        &entry_runs,
        &exit_runs,
        &entry_index,
        entry_ring,
        exit_ring,
    )?;

    // One fragment per entry run.
    let mut fragments = Vec::with_capacity(n_runs);
    for (i, run) in entry_runs.iter().enumerate() {
        let tool_face = run.face;
        let surface = match &store.face(tool_face)?.surface {
            FaceSurface::Nurbs(s) => s.clone(),
            _ => {
                return Err(OperationError::Failed(
                    "through-cut band requires a NURBS tool side face".into(),
                )
                .into())
            }
        };
        let exit_run = exit_run_of(tool_face)?;

        // Trim ribbon between the two concatenated tool-UV traces.
        let mut entry_uv: Vec<Point2> = Vec::new();
        for &k in &run.indices {
            entry_uv.extend_from_slice(&entry.segments[k].branch.uv_b);
        }
        let mut exit_uv: Vec<Point2> = Vec::new();
        for &k in &exit_run.indices {
            exit_uv.extend_from_slice(&exit.segments[k].branch.uv_b);
        }
        let entry_trace = clamp_trace(&entry_uv, &surface);
        let exit_trace = clamp_trace(&exit_uv, &surface);
        let outer = stitch_band_loop(&entry_trace, &exit_trace)?;
        let trim = FaceTrim::new(outer, Vec::new());

        let prev = &entry_runs[(i + n_runs - 1) % n_runs];
        let next = &entry_runs[(i + 1) % n_runs];
        let (kink_start, _) = kink_edges[&run_pair_key(&entry_index, prev, run)?];
        let (kink_end, end_exit_vertex) = kink_edges[&run_pair_key(&entry_index, run, next)?];

        let exit_edges =
            exit_edges_from(exit_run, exit_ring, exit.segments.len(), end_exit_vertex)?;

        let mut cycle: Vec<EdgeId> = run.indices.iter().map(|&k| entry_ring.edges[k]).collect();
        cycle.push(kink_end);
        cycle.extend(exit_edges);
        cycle.push(kink_start);
        let wire_edges = orient_cycle(store, &cycle)?;
        let wire = store.add_wire(WireData {
            edges: wire_edges,
            is_closed: true,
        });

        // Subtract: band normals point INTO the hole (`same_sense = false`).
        let face = store.add_face(FaceData {
            surface: FaceSurface::Nurbs(surface),
            outer_wire: wire,
            inner_wires: Vec::new(),
            same_sense: false,
            trim: Some(trim),
            pcurves: Vec::new(),
        });
        fragments.push(BandFragment { tool_face, face });
    }
    Ok(fragments)
}

/// Builds the band fragments of an OPEN (cap-touching) through cut: one
/// face per tool-face run, sharing the entry / exit trace edges with the
/// split target-face fragments and NEW kink-crossing edges with adjacent
/// band fragments — plus TWO cap-plane closure edges joining the entry and
/// exit chains' matching terminals. The closure edges are returned so the
/// cap-notch rebuild can share them (the notched cap gains exactly these
/// edges — watertight by construction).
///
/// The entry and exit chains are normalized to the same geometric direction
/// (their lexicographically smaller terminals come first), so their
/// tool-face runs align index-to-index and their terminals pair start-to-
/// start / end-to-end.
///
/// # Errors
///
/// Returns an error when the chains' tool-face runs disagree, a tool face
/// is not NURBS, or a fragment polygon degenerates.
pub(crate) fn build_open_band_fragments(
    store: &mut TopologyStore,
    entry: &CutChain,
    exit: &CutChain,
    entry_ring: &ChainRing,
    exit_ring: &ChainRing,
) -> Result<(Vec<BandFragment>, [EdgeId; 2])> {
    let entry_runs = open_tool_runs(entry);
    let exit_runs = open_tool_runs(exit);
    let n_runs = entry_runs.len();
    if n_runs == 0 || exit_runs.len() != n_runs {
        return Err(OperationError::Failed(
            "cap-touching through cut requires entry and exit chains \
             crossing the same tool side faces"
                .into(),
        )
        .into());
    }
    for (a, b) in entry_runs.iter().zip(&exit_runs) {
        if a.face != b.face {
            return Err(OperationError::Failed(
                "entry and exit chained loops cross different tool side \
                 faces"
                    .into(),
            )
            .into());
        }
    }

    // Terminal closure edges (start / end), lying in the cap planes.
    let start_closure = closure_edge(store, entry, exit, entry_ring, exit_ring, false)?;
    let end_closure = closure_edge(store, entry, exit, entry_ring, exit_ring, true)?;

    // One interior kink-crossing edge per run junction (i | i + 1).
    let mut kinks: Vec<EdgeId> = Vec::with_capacity(n_runs.saturating_sub(1));
    for i in 0..n_runs.saturating_sub(1) {
        let entry_idx = entry_runs[i + 1].first();
        let exit_idx = exit_runs[i + 1].first();
        kinks.push(straight_edge(
            store,
            entry_ring.junctions[entry_idx],
            exit_ring.junctions[exit_idx],
            entry.segments[entry_idx].branch.points[0],
            exit.segments[exit_idx].branch.points[0],
        )?);
    }

    let mut fragments = Vec::with_capacity(n_runs);
    for (i, run) in entry_runs.iter().enumerate() {
        let tool_face = run.face;
        let surface = match &store.face(tool_face)?.surface {
            FaceSurface::Nurbs(s) => s.clone(),
            _ => {
                return Err(OperationError::Failed(
                    "through-cut band requires a NURBS tool side face".into(),
                )
                .into())
            }
        };
        let exit_run = &exit_runs[i];

        // Trim ribbon between the two concatenated tool-UV traces.
        let mut entry_uv: Vec<Point2> = Vec::new();
        for &k in &run.indices {
            entry_uv.extend_from_slice(&entry.segments[k].branch.uv_b);
        }
        let mut exit_uv: Vec<Point2> = Vec::new();
        for &k in &exit_run.indices {
            exit_uv.extend_from_slice(&exit.segments[k].branch.uv_b);
        }
        let entry_trace = clamp_trace(&entry_uv, &surface);
        let exit_trace = clamp_trace(&exit_uv, &surface);
        let outer = stitch_band_loop(&entry_trace, &exit_trace)?;
        let trim = FaceTrim::new(outer, Vec::new());

        let left = if i == 0 { start_closure } else { kinks[i - 1] };
        let right = if i == n_runs - 1 {
            end_closure
        } else {
            kinks[i]
        };
        let mut cycle: Vec<EdgeId> = run.indices.iter().map(|&k| entry_ring.edges[k]).collect();
        cycle.push(right);
        cycle.extend(exit_run.indices.iter().rev().map(|&k| exit_ring.edges[k]));
        cycle.push(left);
        let wire_edges = orient_cycle(store, &cycle)?;
        let wire = store.add_wire(WireData {
            edges: wire_edges,
            is_closed: true,
        });

        // Subtract: band normals point INTO the doorway (`same_sense = false`).
        let face = store.add_face(FaceData {
            surface: FaceSurface::Nurbs(surface),
            outer_wire: wire,
            inner_wires: Vec::new(),
            same_sense: false,
            trim: Some(trim),
            pcurves: Vec::new(),
        });
        fragments.push(BandFragment { tool_face, face });
    }
    Ok((fragments, [start_closure, end_closure]))
}

/// Builds one terminal closure edge between the entry and exit chains'
/// matching terminals (`tail == false`: chain heads; `true`: chain tails).
/// The edge is a straight segment in the cap plane, exact for extruded
/// tools whose section at the cap is a straight line across the target
/// thickness.
fn closure_edge(
    store: &mut TopologyStore,
    entry: &CutChain,
    exit: &CutChain,
    entry_ring: &ChainRing,
    exit_ring: &ChainRing,
    tail: bool,
) -> Result<EdgeId> {
    let terminal = |chain: &CutChain, ring: &ChainRing| -> (VertexId, crate::math::Point3) {
        if tail {
            let seg = chain.segments.last().unwrap_or_else(|| unreachable!());
            let p = *seg.branch.points.last().unwrap_or_else(|| unreachable!());
            (*ring.junctions.last().unwrap_or_else(|| unreachable!()), p)
        } else {
            (ring.junctions[0], chain.segments[0].branch.points[0])
        }
    };
    let (entry_v, entry_p) = terminal(entry, entry_ring);
    let (exit_v, exit_p) = terminal(exit, exit_ring);
    straight_edge(store, entry_v, exit_v, entry_p, exit_p)
}

/// Builds the pocket band fragments of a multi-face blind cut: ONE face per
/// tool side face crossed by the entry chain, each running from its entry
/// chain segment down to that side face's shared buried ring edge (whose
/// other incident face is the buried cap — the pocket floor), with NEW
/// kink-crossing edges from the entry junctions to the tool's own buried ring
/// corner vertices.
///
/// # Errors
///
/// Returns an error when a tool face is not NURBS, adjacent buried ring
/// edges share no corner vertex, or a fragment polygon degenerates.
pub(crate) fn build_pocket_band_fragments(
    store: &mut TopologyStore,
    entry: &CutChain,
    entry_ring: &ChainRing,
    buried: &super::pocket::BuriedChainEnd,
) -> Result<Vec<BandFragment>> {
    let n = entry.segments.len();
    if n < 3 || buried.rings.len() != n {
        return Err(OperationError::Failed(
            "multi-face pocket cut requires an entry chain crossing three or \
             more tool side faces with one buried ring edge each"
                .into(),
        )
        .into());
    }

    // One kink-crossing edge per junction: from the entry junction vertex
    // down to the buried ring corner shared by the two adjacent ring edges.
    let mut kinks: Vec<EdgeId> = Vec::with_capacity(n);
    for i in 0..n {
        let prev_ring = buried.rings[(i + n - 1) % n].0;
        let next_ring = buried.rings[i].0;
        let corner = common_vertex(store, prev_ring, next_ring)?;
        let corner_point = store.vertex(corner)?.point;
        let entry_point = entry.segments[i].branch.points[0];
        kinks.push(straight_edge(
            store,
            entry_ring.junctions[i],
            corner,
            entry_point,
            corner_point,
        )?);
    }

    let mut fragments = Vec::with_capacity(n);
    for (i, seg) in entry.segments.iter().enumerate() {
        let tool_face = seg.tool_face;
        let surface = match &store.face(tool_face)?.surface {
            FaceSurface::Nurbs(s) => s.clone(),
            _ => {
                return Err(OperationError::Failed(
                    "pocket band requires a NURBS tool side face".into(),
                )
                .into())
            }
        };
        let (ring_edge, v_boundary) = buried.rings[i];
        let entry_trace = clamp_trace(&seg.branch.uv_b, &surface);
        let buried_uv = super::pocket::buried_edge_uv(store, ring_edge, v_boundary)?;
        let outer = stitch_band_loop(&entry_trace, &buried_uv)?;
        let trim = FaceTrim::new(outer, Vec::new());

        let wire_edges = orient_cycle(
            store,
            &[entry_ring.edges[i], kinks[(i + 1) % n], ring_edge, kinks[i]],
        )?;
        let wire = store.add_wire(WireData {
            edges: wire_edges,
            is_closed: true,
        });

        // Pocket band normals point INTO the cavity (`same_sense = false`).
        let face = store.add_face(FaceData {
            surface: FaceSurface::Nurbs(surface),
            outer_wire: wire,
            inner_wires: Vec::new(),
            same_sense: false,
            trim: Some(trim),
            pcurves: Vec::new(),
        });
        fragments.push(BandFragment { tool_face, face });
    }
    Ok(fragments)
}

/// The vertex shared by two edges (a tool ring corner).
fn common_vertex(store: &TopologyStore, a: EdgeId, b: EdgeId) -> Result<VertexId> {
    let ea = store.edge(a)?;
    let eb = store.edge(b)?;
    for va in [ea.start, ea.end] {
        if va == eb.start || va == eb.end {
            return Ok(va);
        }
    }
    Err(OperationError::Failed("adjacent buried ring edges share no corner vertex".into()).into())
}

/// The exit run's ring edges in cycle-adjacent order: after the end kink
/// edge the boundary walk stands on the exit junction vertex; the exit run
/// is traversed from whichever of its two ends that vertex is.
fn exit_edges_from(
    exit_run: &ToolRun,
    exit_ring: &ChainRing,
    n_exit: usize,
    end_exit_vertex: VertexId,
) -> Result<Vec<EdgeId>> {
    let exit_start_vertex = exit_ring.junctions[exit_run.first()];
    let exit_end_vertex = exit_ring.junctions[(exit_run.last() + 1) % n_exit];
    if end_exit_vertex == exit_start_vertex {
        Ok(exit_run
            .indices
            .iter()
            .map(|&k| exit_ring.edges[k])
            .collect())
    } else if end_exit_vertex == exit_end_vertex {
        Ok(exit_run
            .indices
            .iter()
            .rev()
            .map(|&k| exit_ring.edges[k])
            .collect())
    } else {
        Err(OperationError::Failed(
            "entry and exit chained loops disagree on tool kink crossings".into(),
        )
        .into())
    }
}

/// The unordered entry-run-index pair identifying the junction between two
/// adjacent tool-face runs.
fn run_pair_key(
    entry_index: &HashMap<FaceId, usize>,
    prev: &ToolRun,
    next: &ToolRun,
) -> Result<(usize, usize)> {
    let (Some(&a), Some(&b)) = (entry_index.get(&prev.face), entry_index.get(&next.face)) else {
        return Err(OperationError::Failed(
            "entry and exit chained loops cross different tool side faces".into(),
        )
        .into());
    };
    Ok((a.min(b), a.max(b)))
}

/// Builds one shared kink-crossing edge per junction face pair, recording
/// the exit junction vertex it lands on (for exit-run direction resolution).
#[allow(clippy::too_many_arguments)]
fn build_run_kink_edges(
    store: &mut TopologyStore,
    entry: &CutChain,
    exit: &CutChain,
    entry_runs: &[ToolRun],
    exit_runs: &[ToolRun],
    entry_index: &HashMap<FaceId, usize>,
    entry_ring: &ChainRing,
    exit_ring: &ChainRing,
) -> Result<HashMap<(usize, usize), (EdgeId, VertexId)>> {
    let n_runs = entry_runs.len();
    let mut kink_edges: HashMap<(usize, usize), (EdgeId, VertexId)> = HashMap::new();
    for i in 0..n_runs {
        let prev = &entry_runs[(i + n_runs - 1) % n_runs];
        let next = &entry_runs[i];
        let key = run_pair_key(entry_index, prev, next)?;

        // The matching exit junction: the exit run boundary between the same
        // two tool faces.
        let mut exit_vertex: Option<(usize, VertexId)> = None;
        for j in 0..n_runs {
            let e_prev = &exit_runs[(j + n_runs - 1) % n_runs];
            let e_next = &exit_runs[j];
            if run_pair_key(entry_index, e_prev, e_next)? == key {
                let idx = e_next.first();
                exit_vertex = Some((idx, exit_ring.junctions[idx]));
                break;
            }
        }
        let Some((exit_idx, exit_v)) = exit_vertex else {
            return Err(OperationError::Failed(
                "entry and exit chained loops disagree on tool kink crossings".into(),
            )
            .into());
        };

        let entry_idx = next.first();
        let entry_point = entry.segments[entry_idx].branch.points[0];
        let exit_point = exit.segments[exit_idx].branch.points[0];
        let edge = straight_edge(
            store,
            entry_ring.junctions[entry_idx],
            exit_v,
            entry_point,
            exit_point,
        )?;
        kink_edges.insert(key, (edge, exit_v));
    }
    Ok(kink_edges)
}

/// A straight (degree-1) 3D edge between two existing vertices.
fn straight_edge(
    store: &mut TopologyStore,
    start: VertexId,
    end: VertexId,
    start_point: crate::math::Point3,
    end_point: crate::math::Point3,
) -> Result<EdgeId> {
    let curve = NurbsCurve3D::polyline(&[start_point, end_point])?;
    let (t0, t1) = curve.parameter_domain();
    Ok(store.add_edge(EdgeData {
        start,
        end,
        curve: EdgeCurve::Nurbs(curve),
        t_start: t0,
        t_end: t1,
    }))
}

/// Orients a cyclic edge list into a closed wire: the first edge runs
/// forward, every following edge is flipped as needed to continue from the
/// previous end vertex, and the cycle must return to the first edge's start.
fn orient_cycle(store: &TopologyStore, edges: &[EdgeId]) -> Result<Vec<OrientedEdge>> {
    let first = store.edge(edges[0])?;
    let start_vertex = first.start;
    let mut current = first.end;
    let mut oriented = vec![OrientedEdge::new(edges[0], true)];
    for &e in &edges[1..] {
        let edge = store.edge(e)?;
        if edge.start == current {
            oriented.push(OrientedEdge::new(e, true));
            current = edge.end;
        } else if edge.end == current {
            oriented.push(OrientedEdge::new(e, false));
            current = edge.start;
        } else {
            return Err(OperationError::Failed(
                "band fragment boundary edges do not form a closed cycle".into(),
            )
            .into());
        }
    }
    if current != start_vertex {
        return Err(OperationError::Failed(
            "band fragment boundary edges do not close back to their start".into(),
        )
        .into());
    }
    Ok(oriented)
}

/// Clamps a UV trace into the surface's parameter domain (the SSI corrector may
/// land a hair outside on the seam side) and deduplicates.
fn clamp_trace(uv: &[Point2], surface: &NurbsSurface) -> Vec<Point2> {
    let ((u0, u1), (v0, v1)) = surface.parameter_domain();
    let mut out: Vec<Point2> = Vec::with_capacity(uv.len());
    for p in uv {
        let c = Point2::new(p.x.clamp(u0, u1), p.y.clamp(v0, v1));
        if out.last().is_none_or(|q| (c - q).norm() > 1e-9) {
            out.push(c);
        }
    }
    out
}

/// Stitches the entry and exit traces into a single simple band polygon (a
/// CCW outer trim loop of degree-1 segments).
///
/// The entry trace is walked in its natural (u-increasing) direction and the
/// exit trace reversed (u-decreasing), so the ribbon closes without crossing.
/// Winding is normalized to counter-clockwise (the trim outer convention).
fn stitch_band_loop(entry: &[Point2], exit: &[Point2]) -> Result<TrimLoop> {
    if entry.len() < 2 || exit.len() < 2 {
        return Err(
            OperationError::Failed("band trace degenerated to fewer than 2 points".into()).into(),
        );
    }

    // Order each trace by ascending u so the stitch direction is unambiguous.
    let mut e = entry.to_vec();
    let mut x = exit.to_vec();
    e.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal));
    x.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal));

    // Polygon: entry (u up) then exit (u down).
    let mut poly: Vec<Point2> = Vec::with_capacity(e.len() + x.len());
    poly.extend_from_slice(&e);
    poly.extend(x.iter().rev().copied());
    dedup_closed(&mut poly);

    if poly.len() < 3 {
        return Err(OperationError::Failed(
            "stitched band polygon degenerated to fewer than 3 points".into(),
        )
        .into());
    }

    // Normalize to CCW for the trim outer convention.
    if signed_area(&poly) < 0.0 {
        poly.reverse();
    }

    let n = poly.len();
    let mut curves = Vec::with_capacity(n);
    for i in 0..n {
        curves.push(uv_segment(poly[i], poly[(i + 1) % n]));
    }
    Ok(TrimLoop::new(curves))
}

/// A degree-1 two-point UV line segment.
fn uv_segment(a: Point2, b: Point2) -> NurbsCurve2D {
    NurbsCurve2D::from_unweighted(
        vec![a, b],
        KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap_or_else(|_| unreachable!()),
        1,
    )
    .unwrap_or_else(|_| unreachable!())
}

/// Removes consecutive near-duplicate points and a coincident wrap point.
fn dedup_closed(pts: &mut Vec<Point2>) {
    pts.dedup_by(|a, b| (*a - *b).norm() < 1e-9);
    while pts.len() >= 2 && (pts[0] - pts[pts.len() - 1]).norm() < 1e-9 {
        pts.pop();
    }
}

/// Shoelace signed area. Positive = counter-clockwise.
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

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::nurbs::InversionOptions;
    use crate::math::Point3;
    use crate::operations::boolean::nurbs::loops::{collect_nurbs_faces, extract_cut_loops};
    use crate::operations::boolean::nurbs::punch::punch_loop;
    use crate::operations::creation::{MakeCurvedSlab, MakeNurbsTube};
    use crate::tessellation::{TessellateFace, TessellationParams};
    use crate::topology::SolidId;
    use std::collections::HashMap;

    fn solid_faces(store: &TopologyStore, solid: SolidId) -> Vec<FaceId> {
        let shell = store
            .shell(store.solid(solid).unwrap().outer_shell)
            .unwrap();
        shell.faces.clone()
    }

    /// Builds the band face for a slab×tube and returns (store, band face,
    /// tool surface). The two hole loops are punched first (as the real pipeline
    /// does) so the band shares their ring wires.
    fn band_face(radius: f64) -> (TopologyStore, FaceId, NurbsSurface) {
        let (store, band, surf, _rings) = band_face_with_rings(radius);
        (store, band, surf)
    }

    /// Like [`band_face`] but also returns the entry/exit ring wires the band
    /// shares with the punched faces.
    fn band_face_with_rings(radius: f64) -> (TopologyStore, FaceId, NurbsSurface, BandRingWires) {
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
        let tool_surf = tool[0].1.clone();
        // Punch both loops (entry then exit) and share their ring wires.
        let crate::operations::boolean::nurbs::loops::ToolFaceCut::Through { tool_face, loops } =
            &cuts[0]
        else {
            panic!("expected a through cut");
        };
        let entry = punch_loop(&mut store, &loops[0]).unwrap();
        let exit = punch_loop(&mut store, &loops[1]).unwrap();
        let rings = BandRingWires { entry, exit };
        let band = build_band_face(&mut store, *tool_face, loops, rings).unwrap();
        (store, band, tool_surf, rings)
    }

    #[test]
    fn band_mesh_is_watertight() {
        let (store, band, _) = band_face(0.7);
        let mesh = TessellateFace::new(band, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty(), "band produced no triangles");

        // Edge-manifold check: every undirected edge used 1 or 2 times.
        let mut counts: HashMap<(u32, u32), usize> = HashMap::new();
        for tri in &mesh.indices {
            for &(a, b) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let key = if a < b { (a, b) } else { (b, a) };
                *counts.entry(key).or_insert(0) += 1;
            }
        }
        for (&(a, b), &c) in &counts {
            assert!(c == 1 || c == 2, "band edge ({a},{b}) used {c} times");
        }
    }

    #[test]
    fn band_vertices_lie_on_tool_surface() {
        let (store, band, tool_surf) = band_face(0.7);
        let mesh = TessellateFace::new(band, TessellationParams::default())
            .execute(&store)
            .unwrap();
        let opts = InversionOptions::default();
        for v in &mesh.vertices {
            let inv = tool_surf.closest_point(v, &opts).unwrap();
            assert!(
                inv.distance < 1e-6,
                "band vertex off tool surface: d = {}",
                inv.distance
            );
        }
    }

    #[test]
    fn band_z_extent_spans_slab_thickness() {
        let (store, band, _) = band_face(0.7);
        let mesh = TessellateFace::new(band, TessellationParams::default())
            .execute(&store)
            .unwrap();
        let zmin = mesh
            .vertices
            .iter()
            .map(|p| p.z)
            .fold(f64::INFINITY, f64::min);
        let zmax = mesh
            .vertices
            .iter()
            .map(|p| p.z)
            .fold(f64::NEG_INFINITY, f64::max);
        // The band runs from the back face (z ~ -1 .. 0.5) up to the front face
        // (z ~ 0 .. 1.5); its vertical extent should be a meaningful fraction of
        // the slab thickness (>= 0.5).
        assert!(
            zmax - zmin > 0.5,
            "band z-extent {} too small (zmin={zmin}, zmax={zmax})",
            zmax - zmin
        );
    }

    #[test]
    fn band_boundary_is_exactly_the_two_ring_wires() {
        let (store, band, _surf, rings) = band_face_with_rings(0.7);
        let face = store.face(band).unwrap();
        // The band's outer wire is the entry ring; its single inner wire is the
        // exit ring — the exact wires the punch step created.
        assert_eq!(face.outer_wire, rings.entry, "outer wire = entry ring");
        assert_eq!(face.inner_wires.len(), 1, "exactly one inner wire");
        assert_eq!(face.inner_wires[0], rings.exit, "inner wire = exit ring");
        // The two rings are distinct wires (entry != exit).
        assert_ne!(rings.entry, rings.exit, "entry and exit rings differ");
    }
}
