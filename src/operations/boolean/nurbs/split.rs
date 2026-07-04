//! F3b target face splitting (F5 Phase C).
//!
//! A chained cut loop that crosses target face boundaries (a window sliding
//! across a segmented-prism joint) or the target's own parametric seam (a
//! revolved wall cut across its seam) cannot be punched as an interior trim
//! hole: on each affected face the cut is an OPEN trace from boundary point
//! to boundary point. This module splits such faces along their traces:
//!
//! 1. The face's outer UV rectangle is subdivided by the traces (a UV trim
//!    boolean: outer loop ∪ trace → two loops per trace).
//! 2. Each region is classified against the tool: at a trace sample the
//!    region-side tangent direction is compared with the tool face's outward
//!    normal; the region lying inside the tool is removed material, the rest
//!    is kept. Classification is local and sign-based — no point-in-solid,
//!    no tolerances beyond the SSI marcher's own bounds.
//! 3. Every kept region becomes a fragment face sharing the trace edges
//!    (degree-1 polylines through the exact SSI samples, with matching
//!    pcurves on both the fragment and — via the shared [`EdgeId`]s — the
//!    band fragments), and sharing split boundary sub-edges with the
//!    neighboring fragments (F2 shared-edge topology).
//!
//! ## Name evolution
//!
//! - One kept fragment (a boundary notch — the common window case): the
//!   parent face's name TRANSFERS unchanged. A notched wall face is still
//!   the wall's face, exactly like a punched face keeps its name; a window
//!   sliding across a joint therefore never churns the wall face names.
//! - Two kept fragments (the cut severs the face): the parent name retires
//!   and the fragments bind [`FaceName::Split`] with [`SplitSide::Left`] /
//!   [`SplitSide::Right`] per the deterministic canonical-trace rule
//!   documented on [`SplitSide`].
//! - Three or more kept fragments: a typed unsupported error.
//!
//! Unnamed parents propagate unnamed fragments.

use std::collections::HashMap;

use crate::error::{OperationError, Result};
use crate::geometry::nurbs::{KnotVector, NurbsCurve2D, NurbsCurve3D, NurbsSurface};
use crate::math::{Point2, Point3, TOLERANCE};
use crate::topology::{
    EdgeCurve, EdgeData, EdgeId, FaceData, FaceId, FacePcurve, FaceSurface, FaceTrim, OpId,
    OrientedEdge, TopologyStore, TrimLoop, VertexData, VertexId, WireData,
};

use super::loops::CutLoop;
use super::stitch::CutChain;

/// Exactness bound for UV bookkeeping (side assignment, pinned boundary
/// coordinates). All compared values are constructed exactly (welded pins,
/// iso pcurves), so this only absorbs floating-point noise — it is NOT a
/// geometric tolerance.
const UV_EXACT: f64 = 1e-9;

/// The 3D topology of one chained cut loop: one degree-1 polyline edge per
/// chain segment through the exact SSI samples, shared junction vertices,
/// and a target-UV pcurve per edge (knot-compatible with the 3D edge, so
/// edge-driven tessellation reproduces exactly the SSI sample points).
#[derive(Debug, Clone)]
pub(crate) struct ChainTopology {
    /// One edge per chain segment, in chain order.
    pub edges: Vec<EdgeId>,
    /// `junctions[i]` is the shared vertex at segment `i`'s START.
    pub junctions: Vec<VertexId>,
    /// Target-UV pcurve of `edges[i]`, parameterized identically.
    pub pcurves: Vec<NurbsCurve2D>,
}

/// Builds the shared trace topology of a chained loop (see
/// [`ChainTopology`]).
///
/// # Errors
///
/// Returns an error when a segment degenerates to fewer than 2 distinct
/// samples or curve construction fails.
pub(crate) fn build_chain_topology(
    store: &mut TopologyStore,
    chain: &CutChain,
) -> Result<ChainTopology> {
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
    let mut pcurves = Vec::with_capacity(n);
    for (i, seg) in chain.segments.iter().enumerate() {
        let (edge, pcurve) = trace_edge(
            store,
            &seg.branch.points,
            &seg.branch.uv_a,
            junctions[i],
            junctions[(i + 1) % n],
        )?;
        edges.push(edge);
        pcurves.push(pcurve);
    }
    Ok(ChainTopology {
        edges,
        junctions,
        pcurves,
    })
}

/// Builds one trace edge: a degree-1 3D polyline through the SSI samples
/// plus its target-UV pcurve on the SAME knot vector, so the edge sample
/// cache (degree-1 breakpoint rule) reproduces exactly the SSI sample
/// points on both sides of the shared edge.
fn trace_edge(
    store: &mut TopologyStore,
    points: &[Point3],
    uv: &[Point2],
    start: VertexId,
    end: VertexId,
) -> Result<(EdgeId, NurbsCurve2D)> {
    // Deduplicate 3D and UV in lockstep so the pcurve stays synchronized.
    let mut pts: Vec<Point3> = Vec::with_capacity(points.len());
    let mut uvs: Vec<Point2> = Vec::with_capacity(uv.len());
    for (p, q) in points.iter().zip(uv) {
        if pts
            .last()
            .is_none_or(|last| (*p - *last).norm() > TOLERANCE)
        {
            pts.push(*p);
            uvs.push(*q);
        }
    }
    if pts.len() < 2 {
        return Err(OperationError::Failed(
            "chain segment degenerated to fewer than 2 distinct 3D points".into(),
        )
        .into());
    }
    let curve = NurbsCurve3D::polyline(&pts)?;
    let knots = KnotVector::new(curve.knots().as_slice().to_vec())?;
    let pcurve = NurbsCurve2D::from_unweighted(uvs, knots, 1)?;
    let (t0, t1) = curve.parameter_domain();
    let edge = store.add_edge(EdgeData {
        start,
        end,
        curve: EdgeCurve::Nurbs(curve),
        t_start: t0,
        t_end: t1,
    });
    Ok((edge, pcurve))
}

/// One contiguous run of chain segments on a single target face — an open
/// trace from face boundary to face boundary, ready for splitting.
#[derive(Debug, Clone)]
pub(crate) struct TraceRun {
    /// The target face the run lies on.
    pub target_face: FaceId,
    /// The run's segments in chain order (welded).
    pub segments: Vec<CutLoop>,
    /// One shared trace edge per segment (chain order).
    pub edges: Vec<EdgeId>,
    /// Target-UV pcurve per trace edge (chain order).
    pub pcurves: Vec<NurbsCurve2D>,
    /// Deterministic loop index of the owning chain (entry 0 / exit 1).
    pub loop_index: u32,
    /// The junction vertex at the run's start (a target boundary crossing).
    pub start_vertex: VertexId,
    /// The junction vertex at the run's end.
    pub end_vertex: VertexId,
}

impl TraceRun {
    /// The run's UV start point (pinned on the face boundary).
    fn uv_start(&self) -> Point2 {
        self.segments[0].branch.uv_a[0]
    }

    /// The run's UV end point (pinned on the face boundary).
    fn uv_end(&self) -> Point2 {
        *self
            .segments
            .last()
            .unwrap_or_else(|| unreachable!())
            .branch
            .uv_a
            .last()
            .unwrap_or_else(|| unreachable!())
    }

    /// The concatenated, deduplicated UV samples of the run in run order.
    fn uv_points(&self) -> Vec<Point2> {
        let mut out: Vec<Point2> = Vec::new();
        for seg in &self.segments {
            for &p in &seg.branch.uv_a {
                if out.last().is_none_or(|q| (p - *q).norm() > UV_EXACT) {
                    out.push(p);
                }
            }
        }
        out
    }
}

/// Extracts the per-target-face contiguous runs of a target-crossing chain.
///
/// The chain is rotated so it starts at a target-face change, then cut at
/// every change; each maximal cyclic run becomes one [`TraceRun`].
///
/// # Errors
///
/// Returns an error when the chain does not cross target faces (callers use
/// the interior-punch path for those).
pub(crate) fn trace_runs(
    chain: &CutChain,
    topo: &ChainTopology,
    loop_index: u32,
) -> Result<Vec<TraceRun>> {
    let n = chain.segments.len();
    let Some(first_change) = (0..n)
        .find(|&i| chain.segments[i].target_face != chain.segments[(i + n - 1) % n].target_face)
    else {
        return Err(OperationError::Failed(
            "trace_runs requires a chain crossing target face boundaries".into(),
        )
        .into());
    };

    let mut runs: Vec<TraceRun> = Vec::new();
    let mut i = 0usize;
    while i < n {
        let start = (first_change + i) % n;
        let face = chain.segments[start].target_face;
        let mut len = 1usize;
        while len < n && chain.segments[(start + len) % n].target_face == face {
            len += 1;
        }
        let indices: Vec<usize> = (0..len).map(|k| (start + k) % n).collect();
        runs.push(TraceRun {
            target_face: face,
            segments: indices.iter().map(|&k| chain.segments[k].clone()).collect(),
            edges: indices.iter().map(|&k| topo.edges[k]).collect(),
            pcurves: indices.iter().map(|&k| topo.pcurves[k].clone()).collect(),
            loop_index,
            start_vertex: topo.junctions[start],
            end_vertex: topo.junctions[(start + len) % n],
        });
        i += len;
    }
    Ok(runs)
}

/// One kept fragment of a split target face.
#[derive(Debug, Clone)]
pub(crate) struct Fragment {
    /// The new fragment face.
    pub face: FaceId,
    /// The fragment's outer UV polygon (for interior-hole assignment).
    pub polygon: Vec<Point2>,
    /// The first trace edge on the fragment's boundary (rim naming).
    pub first_trace_edge: EdgeId,
    /// The loop index of the run providing [`Self::first_trace_edge`].
    pub first_loop_index: u32,
}

/// Splits every affected target face along its trace runs and returns the
/// kept fragments per original face. Persistent names transfer (one kept
/// fragment) or split (two kept fragments) per the module-level rules.
///
/// `all_target_faces` is the full target face list, used to guard edge
/// splits: an edge cut by a run endpoint must only be referenced by faces
/// that are themselves being split.
///
/// # Errors
///
/// Returns typed errors for every unsupported configuration: missing
/// pcurves, non-rectangular trim boundaries, cuts landing on existing
/// vertices, ambiguous or tangential classifications, more than two kept
/// fragments, or a split edge shared with an unaffected face.
pub(crate) fn split_target_faces(
    store: &mut TopologyStore,
    affected: &[(FaceId, Vec<TraceRun>)],
    all_target_faces: &[FaceId],
    op_id: Option<&OpId>,
) -> Result<HashMap<FaceId, Vec<Fragment>>> {
    // Phase A: per-face planning (perimeter model, regions, classification,
    // boundary cut registration).
    let mut splitter = EdgeSplitter::default();
    let mut plans: Vec<FacePlan> = Vec::with_capacity(affected.len());
    for (face_id, runs) in affected {
        plans.push(plan_face(store, *face_id, runs, &mut splitter)?);
    }

    // Guard: every split edge is referenced only by affected faces.
    let affected_ids: Vec<FaceId> = affected.iter().map(|(f, _)| *f).collect();
    splitter.assert_only_affected_faces(store, all_target_faces, &affected_ids)?;

    // Phase B: materialize sub-edges.
    splitter.materialize(store)?;

    // Phase C: build fragment faces and evolve names.
    let mut out: HashMap<FaceId, Vec<Fragment>> = HashMap::new();
    for plan in plans {
        let fragments = build_fragments(store, &plan, &splitter)?;
        apply_names(store, plan.face, &fragments, &plan.runs, op_id)?;
        out.insert(plan.face, fragments);
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// Perimeter model
// ---------------------------------------------------------------------------

/// One outer-wire edge mapped onto a domain-rectangle side.
#[derive(Debug, Clone)]
struct SideSpan {
    edge: EdgeId,
    /// Side coordinate at the edge's `t_start` / `t_end`.
    c0: f64,
    c1: f64,
    /// Edge parameter range.
    t0: f64,
    t1: f64,
    /// The face's pcurve for this edge (edge-parameterized UV image).
    pcurve: NurbsCurve2D,
    /// Which UV coordinate is the side coordinate (0 = u, 1 = v).
    coord: usize,
}

impl SideSpan {
    fn c_min(&self) -> f64 {
        self.c0.min(self.c1)
    }
    fn c_max(&self) -> f64 {
        self.c0.max(self.c1)
    }
}

/// The UV domain rectangle of a face with its outer-wire edges assigned to
/// the four sides. Sides may be edge-less (a closed direction's seam), in
/// which case fragments leave a UV gap that the edge-driven tessellation
/// closes with its axis-aligned seam connectors.
#[derive(Debug)]
struct PerimeterModel {
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    /// Per side (0: v=v0, 1: u=u1, 2: v=v1, 3: u=u0), the edges on it.
    sides: [Vec<SideSpan>; 4],
}

impl PerimeterModel {
    /// Builds the model from the face's outer wire and pcurves.
    fn build(store: &TopologyStore, face: &FaceData, surface: &NurbsSurface) -> Result<Self> {
        let ((u0, u1), (v0, v1)) = surface.parameter_domain();
        let mut sides: [Vec<SideSpan>; 4] = Default::default();
        let wire = store.wire(face.outer_wire)?;
        for oe in &wire.edges {
            let Some(pcurve) = face.pcurve_for(oe.edge) else {
                return Err(OperationError::Failed(
                    "face splitting requires per-edge pcurves on the target \
                     face's outer wire"
                        .into(),
                )
                .into());
            };
            let edge = store.edge(oe.edge)?;
            let (t0, t1) = (edge.t_start, edge.t_end);
            let ps = pcurve.point_at(t0)?;
            let pe = pcurve.point_at(t1)?;
            let eps_u = UV_EXACT * (u1 - u0).abs().max(1.0);
            let eps_v = UV_EXACT * (v1 - v0).abs().max(1.0);
            let side = if (ps.y - v0).abs() < eps_v && (pe.y - v0).abs() < eps_v {
                0
            } else if (ps.x - u1).abs() < eps_u && (pe.x - u1).abs() < eps_u {
                1
            } else if (ps.y - v1).abs() < eps_v && (pe.y - v1).abs() < eps_v {
                2
            } else if (ps.x - u0).abs() < eps_u && (pe.x - u0).abs() < eps_u {
                3
            } else {
                return Err(OperationError::Failed(
                    "face splitting requires a rectangular outer boundary \
                     (an outer-wire edge does not lie on a UV domain side)"
                        .into(),
                )
                .into());
            };
            let (c0, c1, coord) = match side {
                0 | 2 => (ps.x, pe.x, 0),
                _ => (ps.y, pe.y, 1),
            };
            sides[side].push(SideSpan {
                edge: oe.edge,
                c0,
                c1,
                t0,
                t1,
                pcurve: pcurve.clone(),
                coord,
            });
        }
        for side in &mut sides {
            side.sort_by(|a, b| {
                a.c_min()
                    .partial_cmp(&b.c_min())
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        }
        Ok(Self {
            u0,
            u1,
            v0,
            v1,
            sides,
        })
    }

    /// The perimeter coordinate `s ∈ [0, 4)` of a UV point pinned on the
    /// rectangle boundary (side + normalized position, counter-clockwise).
    fn s_of(&self, uv: Point2) -> Result<f64> {
        let eps_u = UV_EXACT * (self.u1 - self.u0).abs().max(1.0);
        let eps_v = UV_EXACT * (self.v1 - self.v0).abs().max(1.0);
        let on_west = (uv.x - self.u0).abs() < eps_u;
        let on_east = (self.u1 - uv.x).abs() < eps_u;
        let on_south = (uv.y - self.v0).abs() < eps_v;
        let on_north = (self.v1 - uv.y).abs() < eps_v;
        let du = (self.u1 - self.u0).abs();
        let dv = (self.v1 - self.v0).abs();
        match (on_west, on_east, on_south, on_north) {
            (false, false, true, false) => Ok((uv.x - self.u0) / du),
            (false, true, false, false) => Ok(1.0 + (uv.y - self.v0) / dv),
            (false, false, false, true) => Ok(2.0 + (self.u1 - uv.x) / du),
            (true, false, false, false) => Ok(3.0 + (self.v1 - uv.y) / dv),
            _ => Err(OperationError::Failed(
                "trace endpoint lands on a UV domain corner (degenerate, \
                 unsupported)"
                    .into(),
            )
            .into()),
        }
    }

    /// The UV point of a perimeter coordinate.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    fn uv_of(&self, s: f64) -> Point2 {
        let s = s.rem_euclid(4.0);
        let side = (s.floor() as usize).min(3);
        let f = s - s.floor();
        match side {
            0 => Point2::new(self.u0 + f * (self.u1 - self.u0), self.v0),
            1 => Point2::new(self.u1, self.v0 + f * (self.v1 - self.v0)),
            2 => Point2::new(self.u1 - f * (self.u1 - self.u0), self.v1),
            _ => Point2::new(self.u0, self.v1 - f * (self.v1 - self.v0)),
        }
    }

    /// The side coordinate (`u` for horizontal sides, `v` for vertical) of a
    /// perimeter coordinate on side `side`.
    fn side_coord(&self, side: usize, s: f64) -> f64 {
        let uv = self.uv_of(s);
        match side {
            0 | 2 => uv.x,
            _ => uv.y,
        }
    }
}

// ---------------------------------------------------------------------------
// Edge splitting
// ---------------------------------------------------------------------------

/// One materialized sub-edge of a split boundary edge.
#[derive(Debug, Clone)]
struct SubEdge {
    id: EdgeId,
    /// Edge parameter range (ascending in the parent's parameterization).
    t0: f64,
    t1: f64,
    /// Side coordinates at `t0` / `t1`.
    c0: f64,
    c1: f64,
}

/// Registers boundary-edge cuts (run endpoints) and materializes shared
/// sub-edges once per parent edge.
#[derive(Debug, Default)]
struct EdgeSplitter {
    /// Parent edge → registered cuts `(t, side coordinate, shared vertex)`.
    cuts: HashMap<EdgeId, Vec<(f64, f64, VertexId)>>,
    /// Parent edge → materialized sub-edges, ascending in `t`.
    subs: HashMap<EdgeId, Vec<SubEdge>>,
}

impl EdgeSplitter {
    /// Registers a cut on `edge` at parameter `t` (side coordinate `c`),
    /// sharing `vertex`. Duplicate registrations (the same junction vertex
    /// reached from both adjacent faces) collapse.
    fn register(&mut self, edge: EdgeId, t: f64, c: f64, vertex: VertexId) {
        let cuts = self.cuts.entry(edge).or_default();
        if !cuts.iter().any(|&(_, _, v)| v == vertex) {
            cuts.push((t, c, vertex));
        }
    }

    /// Errors when a cut edge is referenced by a target face that is not
    /// itself being split (its wire would dangle on the removed edge).
    fn assert_only_affected_faces(
        &self,
        store: &TopologyStore,
        all_target_faces: &[FaceId],
        affected: &[FaceId],
    ) -> Result<()> {
        for &fid in all_target_faces {
            if affected.contains(&fid) {
                continue;
            }
            let face = store.face(fid)?;
            let mut wires = vec![face.outer_wire];
            wires.extend(face.inner_wires.iter().copied());
            for w in wires {
                for oe in &store.wire(w)?.edges {
                    if self.cuts.contains_key(&oe.edge) {
                        return Err(OperationError::Failed(
                            "cut trace splits a boundary edge shared with a \
                             face the cut does not cross (unsupported)"
                                .into(),
                        )
                        .into());
                    }
                }
            }
        }
        Ok(())
    }

    /// Materializes the sub-edges of every registered parent edge: the cuts
    /// are sorted along the parent parameterization and each interval
    /// becomes a sub-edge on the SAME curve (restricted parameter range),
    /// sharing the cut junction vertices and the parent's end vertices.
    fn materialize(&mut self, store: &mut TopologyStore) -> Result<()> {
        for (&edge, cuts) in &mut self.cuts {
            let parent = store.edge(edge)?.clone();
            let ascending = parent.t_end > parent.t_start;
            cuts.sort_by(|a, b| {
                let ord = a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal);
                if ascending {
                    ord
                } else {
                    ord.reverse()
                }
            });
            // Guard: cuts strictly inside the parameter range.
            let t_eps = 1e-12 * (parent.t_end - parent.t_start).abs().max(1.0);
            for &(t, _, _) in cuts.iter() {
                if (t - parent.t_start).abs() < t_eps || (t - parent.t_end).abs() < t_eps {
                    return Err(OperationError::Failed(
                        "cut trace lands exactly on a target boundary vertex \
                         (unsupported)"
                            .into(),
                    )
                    .into());
                }
            }

            // Side coordinates of the parent ends for the sub-edge records.
            let mut stations: Vec<(f64, Option<f64>, VertexId)> =
                vec![(parent.t_start, None, parent.start)];
            for &(t, c, v) in cuts.iter() {
                stations.push((t, Some(c), v));
            }
            stations.push((parent.t_end, None, parent.end));

            let mut subs = Vec::with_capacity(stations.len() - 1);
            for w in stations.windows(2) {
                let (t0, c0_opt, v0) = w[0];
                let (t1, c1_opt, v1) = w[1];
                let id = store.add_edge(EdgeData {
                    start: v0,
                    end: v1,
                    curve: parent.curve.clone(),
                    t_start: t0,
                    t_end: t1,
                });
                subs.push(SubEdge {
                    id,
                    t0,
                    t1,
                    // End coordinates are filled by the caller per face (the
                    // side coordinate is face-local); store cut coordinates
                    // when known and NaN placeholders for parent ends —
                    // resolved in `sub_edges_between` via the face's span.
                    c0: c0_opt.unwrap_or(f64::NAN),
                    c1: c1_opt.unwrap_or(f64::NAN),
                });
            }
            self.subs.insert(edge, subs);
        }
        Ok(())
    }

    /// The sub-edges of `span`'s parent edge covering the side-coordinate
    /// interval `[c_lo, c_hi]`, in ascending-coordinate order.
    fn sub_edges_between(&self, span: &SideSpan, c_lo: f64, c_hi: f64) -> Result<Vec<SubEdge>> {
        let subs = self.subs.get(&span.edge).ok_or_else(|| {
            OperationError::Failed("sub-edge lookup on an unsplit boundary edge".into())
        })?;
        // Resolve the parent-end placeholder coordinates from the span.
        let resolved: Vec<SubEdge> = subs
            .iter()
            .map(|s| {
                let mut r = s.clone();
                if r.c0.is_nan() {
                    r.c0 = if (r.t0 - span.t0).abs() <= (r.t0 - span.t1).abs() {
                        span.c0
                    } else {
                        span.c1
                    };
                }
                if r.c1.is_nan() {
                    r.c1 = if (r.t1 - span.t1).abs() <= (r.t1 - span.t0).abs() {
                        span.c1
                    } else {
                        span.c0
                    };
                }
                r
            })
            .collect();
        let eps = UV_EXACT * (c_hi - c_lo).abs().max(1.0);
        let mut picked: Vec<SubEdge> = resolved
            .into_iter()
            .filter(|s| {
                let lo = s.c0.min(s.c1);
                let hi = s.c0.max(s.c1);
                lo >= c_lo - eps && hi <= c_hi + eps
            })
            .collect();
        picked.sort_by(|a, b| {
            a.c0.min(a.c1)
                .partial_cmp(&b.c0.min(b.c1))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        // The picked sub-edges must exactly tile the interval.
        let total: f64 = picked.iter().map(|s| (s.c1 - s.c0).abs()).sum();
        if picked.is_empty() || (total - (c_hi - c_lo).abs()).abs() > 1e-6 {
            return Err(OperationError::Failed(
                "boundary sub-edges do not tile the requested interval \
                 (inconsistent cut bookkeeping)"
                    .into(),
            )
            .into());
        }
        Ok(picked)
    }

    fn is_split(&self, edge: EdgeId) -> bool {
        self.cuts.contains_key(&edge)
    }
}

// ---------------------------------------------------------------------------
// Region subdivision
// ---------------------------------------------------------------------------

/// One boundary piece of a region cycle.
#[derive(Debug, Clone)]
enum Piece {
    /// A counter-clockwise perimeter arc from `s0` to `s1` (`s0 < s1`, no
    /// wrap — cycles are pre-normalized).
    Boundary { s0: f64, s1: f64 },
    /// A trace run traversed forward (start → end) or backward.
    Trace { run: usize, forward: bool },
}

/// One face's split plan: the perimeter model, the region cycles that
/// survive classification, and the runs.
#[derive(Debug)]
struct FacePlan {
    face: FaceId,
    runs: Vec<TraceRun>,
    perimeter: PerimeterModel,
    /// Kept region cycles (counter-clockwise piece lists).
    kept: Vec<Vec<Piece>>,
}

/// Plans one face: subdivides its UV rectangle by the runs, classifies the
/// regions against the tool, and registers boundary-edge cuts for every run
/// endpoint that lands on a real boundary edge.
fn plan_face(
    store: &TopologyStore,
    face_id: FaceId,
    runs: &[TraceRun],
    splitter: &mut EdgeSplitter,
) -> Result<FacePlan> {
    let face = store.face(face_id)?;
    let FaceSurface::Nurbs(surface) = &face.surface else {
        return Err(
            OperationError::Failed("face splitting requires a NURBS target face".into()).into(),
        );
    };
    let surface = surface.clone();
    let perimeter = PerimeterModel::build(store, face, &surface)?;

    // Subdivide: one region initially, split once per run.
    let mut regions: Vec<Vec<Piece>> = vec![vec![Piece::Boundary { s0: 0.0, s1: 4.0 }]];
    for (run_idx, run) in runs.iter().enumerate() {
        let s_a = perimeter.s_of(run.uv_start())?;
        let s_b = perimeter.s_of(run.uv_end())?;
        let region_idx = regions
            .iter()
            .position(|cycle| contains_s(cycle, s_a) && contains_s(cycle, s_b))
            .ok_or_else(|| {
                OperationError::Failed(
                    "trace runs on one face are not nested consistently \
                     (crossing traces are unsupported)"
                        .into(),
                )
            })?;
        let cycle = regions.swap_remove(region_idx);
        let (with_fwd, with_rev) = split_cycle(&cycle, run_idx, s_a, s_b)?;
        regions.push(with_fwd);
        regions.push(with_rev);
    }

    // Classify each region: keep the regions outside the tool.
    let mut kept = Vec::new();
    for cycle in regions {
        if region_is_outside_tool(store, &surface, &cycle, runs)? {
            kept.push(cycle);
        }
    }
    if kept.is_empty() {
        return Err(OperationError::Failed(
            "face splitting classified every fragment as removed material \
             (inconsistent cut)"
                .into(),
        )
        .into());
    }

    // Register boundary cuts for run endpoints on real boundary edges.
    for run in runs {
        for (uv, vertex) in [
            (run.uv_start(), run.start_vertex),
            (run.uv_end(), run.end_vertex),
        ] {
            register_endpoint_cut(&perimeter, uv, vertex, splitter)?;
        }
    }

    Ok(FacePlan {
        face: face_id,
        runs: runs.to_vec(),
        perimeter,
        kept,
    })
}

/// Whether a region cycle's boundary pieces contain perimeter coordinate `s`
/// strictly inside one of them.
fn contains_s(cycle: &[Piece], s: f64) -> bool {
    cycle.iter().any(|p| match p {
        Piece::Boundary { s0, s1 } => s > *s0 + f64::EPSILON && s < *s1 - f64::EPSILON,
        Piece::Trace { .. } => false,
    })
}

/// Splits a region cycle along run `run_idx` (endpoints at perimeter
/// coordinates `s_a` → start, `s_b` → end). Returns the two child cycles;
/// both stay counter-clockwise.
fn split_cycle(
    cycle: &[Piece],
    run_idx: usize,
    s_a: f64,
    s_b: f64,
) -> Result<(Vec<Piece>, Vec<Piece>)> {
    // Explode the cycle into atomic items, splitting the boundary pieces
    // that contain s_a / s_b and remembering the insertion positions.
    #[derive(Debug, Clone)]
    enum Item {
        P(Piece),
        MarkA,
        MarkB,
    }
    let mut items: Vec<Item> = Vec::new();
    for piece in cycle {
        match piece {
            Piece::Boundary { s0, s1 } => {
                let a_in = s_a > *s0 + f64::EPSILON && s_a < *s1 - f64::EPSILON;
                let b_in = s_b > *s0 + f64::EPSILON && s_b < *s1 - f64::EPSILON;
                match (a_in, b_in) {
                    (false, false) => items.push(Item::P(piece.clone())),
                    (true, false) => {
                        items.push(Item::P(Piece::Boundary { s0: *s0, s1: s_a }));
                        items.push(Item::MarkA);
                        items.push(Item::P(Piece::Boundary { s0: s_a, s1: *s1 }));
                    }
                    (false, true) => {
                        items.push(Item::P(Piece::Boundary { s0: *s0, s1: s_b }));
                        items.push(Item::MarkB);
                        items.push(Item::P(Piece::Boundary { s0: s_b, s1: *s1 }));
                    }
                    (true, true) => {
                        let (first, second) = if s_a < s_b { (s_a, s_b) } else { (s_b, s_a) };
                        let mark_first = if s_a < s_b { Item::MarkA } else { Item::MarkB };
                        let mark_second = if s_a < s_b { Item::MarkB } else { Item::MarkA };
                        items.push(Item::P(Piece::Boundary { s0: *s0, s1: first }));
                        items.push(mark_first);
                        items.push(Item::P(Piece::Boundary {
                            s0: first,
                            s1: second,
                        }));
                        items.push(mark_second);
                        items.push(Item::P(Piece::Boundary {
                            s0: second,
                            s1: *s1,
                        }));
                    }
                }
            }
            Piece::Trace { .. } => items.push(Item::P(piece.clone())),
        }
    }

    let pos_a = items
        .iter()
        .position(|i| matches!(i, Item::MarkA))
        .ok_or_else(|| {
            OperationError::Failed("trace endpoint not on the region boundary".into())
        })?;
    // Walk from A: collect pieces until B → path A→B; the rest → path B→A.
    let mut path_ab: Vec<Piece> = Vec::new();
    let mut path_ba: Vec<Piece> = Vec::new();
    let mut in_ab = true;
    let len = items.len();
    for k in 1..len {
        match &items[(pos_a + k) % len] {
            Item::MarkB => in_ab = false,
            Item::MarkA => {}
            Item::P(p) => {
                if in_ab {
                    path_ab.push(p.clone());
                } else {
                    path_ba.push(p.clone());
                }
            }
        }
    }
    // Region closing A→B boundary path with the trace B→A (reversed), and
    // B→A boundary path with the trace A→B (forward).
    let mut region_rev = path_ab;
    region_rev.push(Piece::Trace {
        run: run_idx,
        forward: false,
    });
    let mut region_fwd = path_ba;
    region_fwd.push(Piece::Trace {
        run: run_idx,
        forward: true,
    });
    Ok((region_fwd, region_rev))
}

/// Classifies a region against the tool: `true` when the region lies
/// OUTSIDE the tool (kept material for subtract).
///
/// At each trace piece the region's interior lies on the LEFT of the piece
/// traversal (cycles are counter-clockwise). The left normal in UV, pushed
/// through the target surface partials, gives the 3D direction into the
/// region; the region is inside the tool iff that direction opposes the
/// tool face's outward normal. All trace pieces must agree.
fn region_is_outside_tool(
    store: &TopologyStore,
    surface: &NurbsSurface,
    cycle: &[Piece],
    runs: &[TraceRun],
) -> Result<bool> {
    let mut verdict: Option<bool> = None;
    for piece in cycle {
        let Piece::Trace { run, forward } = piece else {
            continue;
        };
        let run = &runs[*run];
        let uv = run.uv_points();
        if uv.len() < 2 {
            return Err(OperationError::Failed("degenerate trace run".into()).into());
        }
        let mid = uv.len() / 2;
        let lo = mid.saturating_sub(1);
        let hi = (mid + 1).min(uv.len() - 1);
        let mut tangent = uv[hi] - uv[lo];
        if !*forward {
            tangent = -tangent;
        }
        // Left normal in UV → 3D direction into the region.
        let left = Point2::new(-tangent.y, tangent.x);
        let sample_uv = uv[mid];
        let (_, su, sv) = surface.partials(sample_uv.x, sample_uv.y)?;
        let dir = su * left.x + sv * left.y;

        // Tool outward normal at the corresponding tool-UV sample.
        let seg = segment_at_sample(run, mid)?;
        let tool_face = store.face(seg.0)?;
        let FaceSurface::Nurbs(tool_surface) = &tool_face.surface else {
            return Err(OperationError::Failed(
                "face splitting requires NURBS tool side faces".into(),
            )
            .into());
        };
        let tool_uv = seg.1;
        let mut n_tool =
            crate::geometry::surface::Surface::normal(tool_surface, tool_uv.x, tool_uv.y)?;
        if !tool_face.same_sense {
            n_tool = -n_tool;
        }

        let dot = dir.dot(&n_tool);
        if dot.abs() <= f64::EPSILON * dir.norm().max(1.0) {
            return Err(OperationError::Failed(
                "tangential cut trace: cannot classify the fragment against \
                 the tool (grazing cuts are unsupported)"
                    .into(),
            )
            .into());
        }
        let outside = dot > 0.0;
        match verdict {
            None => verdict = Some(outside),
            Some(prev) if prev != outside => {
                return Err(OperationError::Failed(
                    "inconsistent fragment classification: trace pieces of one \
                     region disagree on the tool side"
                        .into(),
                )
                .into());
            }
            Some(_) => {}
        }
    }
    verdict.ok_or_else(|| {
        OperationError::Failed("region without any trace piece (splitting invariant broken)".into())
            .into()
    })
}

/// The (tool face, tool UV) of the run sample at concatenated index `idx`.
fn segment_at_sample(run: &TraceRun, idx: usize) -> Result<(FaceId, Point2)> {
    // Walk the segments with the same dedup rule as `uv_points`.
    let mut count = 0usize;
    let mut last: Option<Point2> = None;
    for seg in &run.segments {
        for (k, &p) in seg.branch.uv_a.iter().enumerate() {
            if last.is_none_or(|q| (p - q).norm() > UV_EXACT) {
                if count == idx {
                    return Ok((seg.tool_face, seg.branch.uv_b[k]));
                }
                count += 1;
                last = Some(p);
            }
        }
    }
    Err(OperationError::Failed("trace sample index out of range".into()).into())
}

/// Registers the boundary-edge cut of one run endpoint, when the endpoint
/// lands on a side that carries real edges (edge-less seam sides need no
/// cut).
fn register_endpoint_cut(
    perimeter: &PerimeterModel,
    uv: Point2,
    vertex: VertexId,
    splitter: &mut EdgeSplitter,
) -> Result<()> {
    let s = perimeter.s_of(uv)?;
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let side = (s.rem_euclid(4.0).floor() as usize).min(3);
    let c = perimeter.side_coord(side, s);
    let spans = &perimeter.sides[side];
    if spans.is_empty() {
        return Ok(()); // Edge-less side (parametric seam): nothing to split.
    }
    let eps = UV_EXACT
        * match side {
            0 | 2 => (perimeter.u1 - perimeter.u0).abs().max(1.0),
            _ => (perimeter.v1 - perimeter.v0).abs().max(1.0),
        };
    for span in spans {
        if c > span.c_min() + eps && c < span.c_max() - eps {
            let t = solve_edge_param(span, c)?;
            splitter.register(span.edge, t, c, vertex);
            return Ok(());
        }
        if (c - span.c_min()).abs() <= eps || (c - span.c_max()).abs() <= eps {
            return Err(OperationError::Failed(
                "cut trace lands exactly on a target boundary vertex \
                 (unsupported)"
                    .into(),
            )
            .into());
        }
    }
    Err(OperationError::Failed(
        "trace endpoint does not lie on any outer-wire edge of the target \
         face"
            .into(),
    )
    .into())
}

/// Solves the edge parameter whose side coordinate equals `c` by monotone
/// bisection on the span's ACTUAL pcurve coordinate (deterministic,
/// converges to f64 precision; no tolerance knob).
fn solve_edge_param(span: &SideSpan, c: f64) -> Result<f64> {
    let (mut t_lo, mut c_lo, mut t_hi, mut c_hi) = (span.t0, span.c0, span.t1, span.c1);
    if c_lo > c_hi {
        std::mem::swap(&mut t_lo, &mut t_hi);
        std::mem::swap(&mut c_lo, &mut c_hi);
    }
    if !(c_lo..=c_hi).contains(&c) {
        return Err(
            OperationError::Failed("cut coordinate outside the boundary edge span".into()).into(),
        );
    }
    let coord_at = |t: f64| -> Result<f64> {
        let p = span.pcurve.point_at(t)?;
        Ok(if span.coord == 0 { p.x } else { p.y })
    };
    let (mut lo, mut hi) = (t_lo, t_hi);
    for _ in 0..128 {
        let mid = 0.5 * (lo + hi);
        // Bit-exhausted interval: the midpoint no longer separates the ends.
        if !(lo.min(hi) < mid && mid < lo.max(hi)) {
            break;
        }
        if coord_at(mid)? < c {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

// ---------------------------------------------------------------------------
// Fragment construction
// ---------------------------------------------------------------------------

/// Builds the kept fragment faces of one planned split.
fn build_fragments(
    store: &mut TopologyStore,
    plan: &FacePlan,
    splitter: &EdgeSplitter,
) -> Result<Vec<Fragment>> {
    let parent = store.face(plan.face)?.clone();
    let FaceSurface::Nurbs(surface) = &parent.surface else {
        return Err(
            OperationError::Failed("face splitting requires a NURBS target face".into()).into(),
        );
    };
    let surface = surface.clone();

    let mut fragments = Vec::with_capacity(plan.kept.len());
    for cycle in &plan.kept {
        let polygon = cycle_polygon(&plan.perimeter, cycle, &plan.runs);
        let build = cycle_wire(&plan.perimeter, cycle, &plan.runs, splitter)?;
        let (first_trace_edge, first_loop_index) = build
            .first_trace
            .ok_or_else(|| OperationError::Failed("kept fragment without any trace edge".into()))?;

        let outer_wire = store.add_wire(WireData {
            edges: build.edges,
            is_closed: true,
        });
        let trim = FaceTrim::new(polygon_trim_loop(&polygon)?, Vec::new());
        let face = store.add_face(FaceData {
            surface: FaceSurface::Nurbs(surface.clone()),
            outer_wire,
            inner_wires: Vec::new(),
            same_sense: parent.same_sense,
            trim: Some(trim),
            pcurves: build.pcurves,
        });
        fragments.push(Fragment {
            face,
            polygon,
            first_trace_edge,
            first_loop_index,
        });
    }
    Ok(fragments)
}

/// The UV polygon of a region cycle (counter-clockwise, deduplicated).
fn cycle_polygon(perimeter: &PerimeterModel, cycle: &[Piece], runs: &[TraceRun]) -> Vec<Point2> {
    let mut poly: Vec<Point2> = Vec::new();
    let push = |p: Point2, poly: &mut Vec<Point2>| {
        if poly.last().is_none_or(|q| (p - *q).norm() > UV_EXACT) {
            poly.push(p);
        }
    };
    for piece in cycle {
        match piece {
            Piece::Boundary { s0, s1 } => {
                push(perimeter.uv_of(*s0), &mut poly);
                // Interior rectangle corners.
                let mut k = s0.floor() + 1.0;
                while k < *s1 - f64::EPSILON {
                    if k > *s0 + f64::EPSILON {
                        push(perimeter.uv_of(k), &mut poly);
                    }
                    k += 1.0;
                }
                push(perimeter.uv_of(*s1), &mut poly);
            }
            Piece::Trace { run, forward } => {
                let mut pts = runs[*run].uv_points();
                if !*forward {
                    pts.reverse();
                }
                for p in pts {
                    push(p, &mut poly);
                }
            }
        }
    }
    while poly.len() >= 2 && (poly[0] - poly[poly.len() - 1]).norm() < UV_EXACT {
        poly.pop();
    }
    poly
}

/// Converts a UV polygon into a degree-1 trim loop.
fn polygon_trim_loop(poly: &[Point2]) -> Result<TrimLoop> {
    if poly.len() < 3 {
        return Err(OperationError::Failed(
            "fragment polygon degenerated to fewer than 3 UV points".into(),
        )
        .into());
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

/// Adds a pcurve entry once per edge.
fn add_pcurve(
    pcurves: &mut Vec<FacePcurve>,
    seen: &mut Vec<EdgeId>,
    edge: EdgeId,
    curve: NurbsCurve2D,
) {
    if !seen.contains(&edge) {
        seen.push(edge);
        pcurves.push(FacePcurve { edge, curve });
    }
}

/// A fragment's assembled outer wire: oriented edges in cycle order, the
/// pcurve table, and the first trace edge encountered (rim naming).
struct WireBuild {
    edges: Vec<OrientedEdge>,
    pcurves: Vec<FacePcurve>,
    first_trace: Option<(EdgeId, u32)>,
}

/// Builds the fragment's outer wire (oriented edges in cycle order) plus its
/// pcurve table, and reports the first trace edge encountered.
fn cycle_wire(
    perimeter: &PerimeterModel,
    cycle: &[Piece],
    runs: &[TraceRun],
    splitter: &EdgeSplitter,
) -> Result<WireBuild> {
    let mut wire: Vec<OrientedEdge> = Vec::new();
    let mut pcurves: Vec<FacePcurve> = Vec::new();
    let mut seen: Vec<EdgeId> = Vec::new();
    let mut first_trace: Option<(EdgeId, u32)> = None;

    for piece in cycle {
        match piece {
            Piece::Trace { run, forward } => {
                let run = &runs[*run];
                if first_trace.is_none() {
                    first_trace = Some((run.edges[0], run.loop_index));
                }
                let indices: Vec<usize> = if *forward {
                    (0..run.edges.len()).collect()
                } else {
                    (0..run.edges.len()).rev().collect()
                };
                for i in indices {
                    wire.push(OrientedEdge::new(run.edges[i], *forward));
                    add_pcurve(
                        &mut pcurves,
                        &mut seen,
                        run.edges[i],
                        run.pcurves[i].clone(),
                    );
                }
            }
            Piece::Boundary { s0, s1 } => {
                append_boundary_edges(
                    perimeter,
                    *s0,
                    *s1,
                    splitter,
                    &mut wire,
                    &mut pcurves,
                    &mut seen,
                )?;
            }
        }
    }
    Ok(WireBuild {
        edges: wire,
        pcurves,
        first_trace,
    })
}

/// Appends the boundary edges (whole or split) covering perimeter arc
/// `[s0, s1]` in counter-clockwise traversal order.
fn append_boundary_edges(
    perimeter: &PerimeterModel,
    s0: f64,
    s1: f64,
    splitter: &EdgeSplitter,
    wire: &mut Vec<OrientedEdge>,
    pcurves: &mut Vec<FacePcurve>,
    seen: &mut Vec<EdgeId>,
) -> Result<()> {
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let first_side = (s0.floor() as usize).min(3);
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let last_side = ((s1 - 1e-12).floor() as usize).min(3);

    for side in first_side..=last_side {
        #[allow(clippy::cast_precision_loss)]
        let side_f = side as f64;
        let seg_s0 = s0.max(side_f);
        let seg_s1 = s1.min(side_f + 1.0);
        if seg_s1 - seg_s0 <= 1e-12 {
            continue;
        }
        let c_from = perimeter.side_coord(side, seg_s0);
        let c_to = perimeter.side_coord(side, seg_s1);
        let spans = &perimeter.sides[side];
        if spans.is_empty() {
            continue; // Edge-less seam side: UV gap, closed by seam connectors.
        }
        let (lo, hi) = (c_from.min(c_to), c_from.max(c_to));
        let ascending = c_to > c_from;

        // Collect (c_min, edge, natural ascending, pcurve) per span overlap.
        let mut collected: Vec<(f64, EdgeId, bool, NurbsCurve2D)> = Vec::new();
        let eps = UV_EXACT * (hi - lo).abs().max(1.0);
        for span in spans {
            let o_lo = span.c_min().max(lo);
            let o_hi = span.c_max().min(hi);
            if o_hi - o_lo <= eps {
                continue;
            }
            let natural_ascending = span.c1 > span.c0;
            let full = (o_lo - span.c_min()).abs() <= eps && (o_hi - span.c_max()).abs() <= eps;
            if full && !splitter.is_split(span.edge) {
                collected.push((o_lo, span.edge, natural_ascending, span.pcurve.clone()));
            } else {
                for sub in splitter.sub_edges_between(span, o_lo, o_hi)? {
                    let sub_ascending = sub.c1 > sub.c0;
                    collected.push((
                        sub.c0.min(sub.c1),
                        sub.id,
                        sub_ascending,
                        span.pcurve.clone(),
                    ));
                }
            }
        }
        collected.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        if !ascending {
            collected.reverse();
        }
        for (_, edge, natural_ascending, pcurve) in collected {
            wire.push(OrientedEdge::new(edge, natural_ascending == ascending));
            add_pcurve(pcurves, seen, edge, pcurve);
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Name evolution
// ---------------------------------------------------------------------------

/// Applies the split naming rules to the kept fragments of one face.
fn apply_names(
    store: &mut TopologyStore,
    parent: FaceId,
    fragments: &[Fragment],
    runs: &[TraceRun],
    op_id: Option<&OpId>,
) -> Result<()> {
    match fragments {
        [] => {
            Err(OperationError::Failed("face splitting produced no kept fragments".into()).into())
        }
        [single] => {
            // A boundary notch: the fragment IS the parent face (same
            // identity, new trim) — the name transfers unchanged.
            store.names_mut().transfer_face(parent, single.face);
            Ok(())
        }
        [a, b] => {
            let Some(op) = op_id else {
                return Ok(()); // Unnamed operation: fragments stay unnamed.
            };
            let (left, right) = order_fragments(a, b, runs)?;
            store.names_mut().split_face(parent, op, left, right);
            Ok(())
        }
        _ => Err(OperationError::Failed(
            "face splitting produced more than two kept fragments \
             (unsupported)"
                .into(),
        )
        .into()),
    }
}

/// Orders two fragments into (left, right) by the canonical-trace rule
/// documented on [`crate::topology::SplitSide`]: canonical trace = the run
/// whose canonically-oriented UV chord starts lexicographically first; a
/// fragment is Left when its polygon centroid lies on the positive
/// cross-product side of that chord.
fn order_fragments(a: &Fragment, b: &Fragment, runs: &[TraceRun]) -> Result<(FaceId, FaceId)> {
    let chord = runs
        .iter()
        .map(|run| canonical_chord(run.uv_start(), run.uv_end()))
        .min_by(|x, y| {
            let kx = (x.0.x, x.0.y, x.1.x, x.1.y);
            let ky = (y.0.x, y.0.y, y.1.x, y.1.y);
            kx.partial_cmp(&ky).unwrap_or(std::cmp::Ordering::Equal)
        })
        .ok_or_else(|| OperationError::Failed("split without trace runs".into()))?;

    let side = |frag: &Fragment| -> f64 {
        let center = centroid(&frag.polygon);
        let dir = chord.1 - chord.0;
        let rel = center - chord.0;
        dir.x * rel.y - dir.y * rel.x
    };
    let side_a = side(a);
    let side_b = side(b);
    if side_a > 0.0 && side_b < 0.0 {
        Ok((a.face, b.face))
    } else if side_a < 0.0 && side_b > 0.0 {
        Ok((b.face, a.face))
    } else {
        Err(OperationError::Failed(
            "split fragments do not lie on opposite sides of the canonical \
             trace (ambiguous SplitSide)"
                .into(),
        )
        .into())
    }
}

/// Orients a chord so `end - start` is lexicographically positive.
fn canonical_chord(a: Point2, b: Point2) -> (Point2, Point2) {
    let d = b - a;
    if d.x > 0.0 || (d.x == 0.0 && d.y > 0.0) {
        (a, b)
    } else {
        (b, a)
    }
}

/// Vertex-average centroid of a UV polygon (adequate for sidedness of
/// convex-ish fragment regions; documented with the [`SplitSide`] rule).
///
/// [`SplitSide`]: crate::topology::SplitSide
fn centroid(poly: &[Point2]) -> Point2 {
    if poly.is_empty() {
        return Point2::new(0.0, 0.0);
    }
    let mut sum = Point2::new(0.0, 0.0);
    for p in poly {
        sum = Point2::new(sum.x + p.x, sum.y + p.y);
    }
    #[allow(clippy::cast_precision_loss)]
    let inv = 1.0 / poly.len() as f64;
    Point2::new(sum.x * inv, sum.y * inv)
}

/// Ray-cast point-in-polygon test in UV (for interior-hole assignment).
pub(crate) fn polygon_contains(poly: &[Point2], point: Point2) -> bool {
    let count = poly.len();
    let mut inside = false;
    for i in 0..count {
        let from = poly[i];
        let to = poly[(i + 1) % count];
        if (from.y > point.y) != (to.y > point.y) {
            let crossing = from.x + (point.y - from.y) / (to.y - from.y) * (to.x - from.x);
            if crossing > point.x {
                inside = !inside;
            }
        }
    }
    inside
}
