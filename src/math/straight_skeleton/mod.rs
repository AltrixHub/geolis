//! Straight skeleton of simple polygons in the XY plane.
//!
//! The skeleton is computed with a shrinking-wavefront algorithm that
//! handles both edge-collapse events and reflex-vertex split events, so
//! non-convex footprints (L/T/U shapes and similar) are supported. The
//! result is returned as one *cell* per input edge: the region of the
//! polygon swept by that edge's wavefront, with every cell vertex carrying
//! its inset distance (the offset time at which the wavefront reached it).
//!
//! Cells are the natural building block for uniform-slope (hip) roofs:
//! lifting each cell vertex by `slope * inset` turns every cell into one
//! planar roof face, because the inset is the exact distance to the cell's
//! supporting edge line.
//!
//! Degenerate simultaneous collapses (parallel opposite edges producing
//! ridges, symmetric footprints collapsing onto a tree of segments) are
//! resolved by closing *flat* wavefront cycles: when a cycle's area reaches
//! zero at the current event time, the cycle has degenerated into a tree of
//! segments traversed twice (an Euler tour), and each segment becomes one
//! skeleton ridge bordered by the two faces that swept onto it.

use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;

use super::{Point2, Point3, Vector2, TOLERANCE};
use crate::error::{OperationError, Result};

/// Welding radius for skeleton nodes: events landing within this distance
/// (and matching inset) reuse the same node. Footprints are meter-scale, so
/// an absolute epsilon is appropriate.
const WELD_EPS: f64 = 1e-6;

/// Epsilon for event-time comparisons.
const TIME_EPS: f64 = 1e-9;

/// Epsilon guarding near-parallel denominators in bisector math.
const PARALLEL_EPS: f64 = 1e-9;

/// Containment slack when locating a point on a wavefront segment.
const FRONT_EPS: f64 = 1e-6;

/// A wavefront cycle whose swept area drops below this is considered
/// collapsed (flat) and is closed into ridge arcs.
const FLAT_AREA_EPS: f64 = 1e-9;

/// Consecutive input vertices closer than this are merged.
const DEDUP_EPS: f64 = 1e-7;

/// A vertex of a skeleton cell: its planar position (z = 0) and the inset
/// distance from the polygon boundary at which the wavefront created it.
#[derive(Clone, Copy, Debug)]
pub struct SkeletonVertex {
    /// Position in the XY plane (z is always 0).
    pub position: Point3,
    /// Offset distance from the polygon boundary (0 on the contour).
    pub inset: f64,
}

/// The skeleton cell of one polygon edge: the region swept by that edge's
/// wavefront, as a CCW simple polygon. The first two vertices are the edge
/// endpoints (inset 0) in edge order.
#[derive(Clone, Debug)]
pub struct SkeletonCell {
    /// Index of the polygon edge this cell belongs to. Edge `i` runs from
    /// `polygon[i]` to `polygon[(i + 1) % n]`.
    pub edge_index: usize,
    /// Cell boundary in CCW order, starting with the edge endpoints.
    pub vertices: Vec<SkeletonVertex>,
}

/// Straight skeleton of a simple polygon, decomposed into per-edge cells.
#[derive(Clone, Debug)]
pub struct StraightSkeleton {
    /// The normalized polygon the cells refer to: CCW, consecutive
    /// duplicates removed. Edge/cell indices are relative to this ring.
    pub polygon: Vec<Point3>,
    /// One cell per polygon edge, in edge order.
    pub cells: Vec<SkeletonCell>,
    /// Largest inset distance reached (the "ridge height" of the skeleton).
    pub max_inset: f64,
}

/// Computes the straight skeleton of a simple polygon in the XY plane.
///
/// The input is interpreted as a closed ring (do not repeat the first
/// point); z coordinates are ignored. Clockwise rings are reversed to CCW.
///
/// # Errors
///
/// Returns [`OperationError::InvalidInput`] when the ring has fewer than 3
/// distinct vertices, contains non-finite coordinates, is self-intersecting,
/// has (near-)zero area, or contains a zero-angle spike. Returns
/// [`OperationError::Failed`] if the wavefront propagation does not
/// converge (numerically degenerate input).
pub fn compute_straight_skeleton(points: &[Point3]) -> Result<StraightSkeleton> {
    let ring = normalize_polygon(points)?;
    let mut builder = Builder::new(&ring)?;
    builder.run()?;
    let cells = builder.extract_cells()?;
    let max_inset = builder
        .nodes
        .iter()
        .map(|n| n.inset)
        .fold(0.0_f64, f64::max);
    let polygon = ring
        .iter()
        .map(|p| Point3::new(p.x, p.y, 0.0))
        .collect::<Vec<_>>();
    Ok(StraightSkeleton {
        polygon,
        cells,
        max_inset,
    })
}

/// Validates and normalizes the input ring: dedup, finiteness, simplicity,
/// CCW orientation.
fn normalize_polygon(points: &[Point3]) -> Result<Vec<Point2>> {
    for p in points {
        if !p.x.is_finite() || !p.y.is_finite() {
            return Err(OperationError::InvalidInput(
                "polygon contains non-finite coordinates".into(),
            )
            .into());
        }
    }
    let mut ring: Vec<Point2> = Vec::with_capacity(points.len());
    for p in points {
        let q = Point2::new(p.x, p.y);
        if let Some(last) = ring.last() {
            if (q - last).norm() < DEDUP_EPS {
                continue;
            }
        }
        ring.push(q);
    }
    while ring.len() >= 2 {
        let first = ring[0];
        let last = ring[ring.len() - 1];
        if (first - last).norm() < DEDUP_EPS {
            ring.pop();
        } else {
            break;
        }
    }
    if ring.len() < 3 {
        return Err(OperationError::InvalidInput(format!(
            "polygon needs at least 3 distinct vertices, got {}",
            ring.len()
        ))
        .into());
    }
    let area = ring_area(&ring);
    if area.abs() < TOLERANCE {
        return Err(OperationError::InvalidInput("polygon has zero area".into()).into());
    }
    if area < 0.0 {
        ring.reverse();
    }
    if let Some((i, j)) = ring_self_intersection(&ring) {
        return Err(OperationError::InvalidInput(format!(
            "polygon is self-intersecting (edges {i} and {j} cross)"
        ))
        .into());
    }
    Ok(ring)
}

/// Signed area of a 2D ring (shoelace; positive for CCW).
pub(crate) fn ring_area(pts: &[Point2]) -> f64 {
    let n = pts.len();
    let mut sum = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        sum += pts[i].x * pts[j].y - pts[j].x * pts[i].y;
    }
    sum * 0.5
}

/// Returns the first pair of non-adjacent edges that cross transversely.
pub(crate) fn ring_self_intersection(pts: &[Point2]) -> Option<(usize, usize)> {
    let n = pts.len();
    for i in 0..n {
        let a0 = pts[i];
        let a1 = pts[(i + 1) % n];
        for j in (i + 2)..n {
            if i == 0 && j == n - 1 {
                continue;
            }
            let b0 = pts[j];
            let b1 = pts[(j + 1) % n];
            if segments_cross_interior(a0, a1, b0, b1) {
                return Some((i, j));
            }
        }
    }
    None
}

/// True when segments `(a0, a1)` and `(b0, b1)` cross strictly in their
/// interiors.
fn segments_cross_interior(a0: Point2, a1: Point2, b0: Point2, b1: Point2) -> bool {
    let da = a1 - a0;
    let db = b1 - b0;
    let denom = cross_2d(da, db);
    if denom.abs() < TOLERANCE {
        return false;
    }
    let d = b0 - a0;
    let t = cross_2d(d, db) / denom;
    let u = cross_2d(d, da) / denom;
    let eps = 1e-9;
    t > eps && t < 1.0 - eps && u > eps && u < 1.0 - eps
}

fn cross_2d(a: Vector2, b: Vector2) -> f64 {
    a.x * b.y - a.y * b.x
}

/// One original polygon edge with its unit direction and inward unit
/// normal. The wavefront line of edge `e` at time `t` is the supporting
/// line translated by `t * normal`.
#[derive(Clone, Copy)]
struct Edge {
    base: Point2,
    dir: Vector2,
    normal: Vector2,
}

impl Edge {
    /// Signed offset distance of `p` from this edge's supporting line
    /// (positive on the polygon interior side).
    fn offset_dist(&self, p: Point2) -> f64 {
        (p - self.base).dot(&self.normal)
    }
}

/// A wavefront vertex: the moving intersection point of the offset lines
/// of its two incident edges.
struct Vert {
    pos: Point2,
    time: f64,
    /// `pos(t) = pos + (t - time) * velocity`. `None` for vertices between
    /// antiparallel fronts (parallel opposite edges): those never move —
    /// their region closes at the same instant they are born.
    velocity: Option<Vector2>,
    left_edge: usize,
    right_edge: usize,
    reflex: bool,
    prev: usize,
    next: usize,
    alive: bool,
    origin_node: usize,
}

/// A skeleton graph node (an event point or a contour vertex).
struct Node {
    pos: Point2,
    inset: f64,
}

/// A skeleton arc between two nodes, bordering exactly two faces.
struct SkelArc {
    from: usize,
    to: usize,
    faces: [usize; 2],
}

#[derive(Clone, Copy)]
enum EventKind {
    Collapse { va: usize, vb: usize },
    Split { v: usize, opposite: usize },
}

struct Event {
    time: f64,
    rank: u8,
    order: (usize, usize),
    point: Point2,
    kind: EventKind,
}

impl PartialEq for Event {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Event {}

impl PartialOrd for Event {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Event {
    fn cmp(&self, other: &Self) -> Ordering {
        self.time
            .total_cmp(&other.time)
            .then(self.rank.cmp(&other.rank))
            .then(self.order.cmp(&other.order))
    }
}

struct Builder {
    edges: Vec<Edge>,
    verts: Vec<Vert>,
    nodes: Vec<Node>,
    arcs: Vec<SkelArc>,
    queue: BinaryHeap<Reverse<Event>>,
}

impl Builder {
    fn new(ring: &[Point2]) -> Result<Self> {
        let n = ring.len();
        let mut edges = Vec::with_capacity(n);
        for i in 0..n {
            let a = ring[i];
            let b = ring[(i + 1) % n];
            let d = b - a;
            let len = d.norm();
            if len < TOLERANCE {
                return Err(OperationError::InvalidInput("zero-length polygon edge".into()).into());
            }
            let dir = d / len;
            edges.push(Edge {
                base: a,
                dir,
                normal: Vector2::new(-dir.y, dir.x),
            });
        }
        let nodes = ring
            .iter()
            .map(|p| Node {
                pos: *p,
                inset: 0.0,
            })
            .collect::<Vec<_>>();
        let mut builder = Self {
            edges,
            verts: Vec::with_capacity(2 * n),
            nodes,
            arcs: Vec::new(),
            queue: BinaryHeap::new(),
        };
        for (i, p) in ring.iter().enumerate() {
            let left = (i + n - 1) % n;
            let right = i;
            let velocity = builder.bisector_velocity(left, right);
            if velocity.is_none() {
                return Err(OperationError::InvalidInput(
                    "polygon has a zero-angle spike vertex".into(),
                )
                .into());
            }
            let reflex = builder.is_reflex(left, right);
            builder.verts.push(Vert {
                pos: *p,
                time: 0.0,
                velocity,
                left_edge: left,
                right_edge: right,
                reflex,
                prev: (i + n - 1) % n,
                next: (i + 1) % n,
                alive: true,
                origin_node: i,
            });
        }
        for i in 0..n {
            builder.push_collapse_event(i, (i + 1) % n);
        }
        for i in 0..n {
            if builder.verts[i].reflex {
                builder.push_split_events(i);
            }
        }
        Ok(builder)
    }

    /// Vertex velocity for edge pair `(left, right)`: the unique `w` with
    /// `w . n_left = w . n_right = 1`, so the vertex stays on both offset
    /// lines. `None` when the normals are antiparallel.
    fn bisector_velocity(&self, left: usize, right: usize) -> Option<Vector2> {
        let nl = self.edges[left].normal;
        let nr = self.edges[right].normal;
        let denom = 1.0 + nl.dot(&nr);
        if denom.abs() < PARALLEL_EPS {
            return None;
        }
        Some((nl + nr) / denom)
    }

    fn is_reflex(&self, left: usize, right: usize) -> bool {
        cross_2d(self.edges[left].dir, self.edges[right].dir) < -PARALLEL_EPS
    }

    fn pos_at(&self, v: usize, t: f64) -> Point2 {
        let vert = &self.verts[v];
        match vert.velocity {
            Some(w) => vert.pos + (t - vert.time) * w,
            None => vert.pos,
        }
    }

    /// Finds or creates the skeleton node at `p` with inset `t`.
    fn weld(&mut self, p: Point2, t: f64) -> usize {
        for (i, node) in self.nodes.iter().enumerate() {
            if (node.pos - p).norm() < WELD_EPS && (node.inset - t).abs() < WELD_EPS {
                return i;
            }
        }
        self.nodes.push(Node { pos: p, inset: t });
        self.nodes.len() - 1
    }

    /// Emits the skeleton arc traced by vertex `v` dying at `node`.
    fn emit_vertex_arc(&mut self, v: usize, node: usize) {
        let (origin, faces) = {
            let vert = &self.verts[v];
            (vert.origin_node, [vert.left_edge, vert.right_edge])
        };
        if origin != node {
            self.arcs.push(SkelArc {
                from: origin,
                to: node,
                faces,
            });
        }
    }

    fn cycle_len(&self, start: usize) -> usize {
        let mut len = 1;
        let mut cur = self.verts[start].next;
        while cur != start && len <= self.verts.len() {
            len += 1;
            cur = self.verts[cur].next;
        }
        len
    }

    fn same_cycle(&self, a: usize, b: usize) -> bool {
        if a == b {
            return true;
        }
        let mut cur = self.verts[a].next;
        let mut steps = 0;
        while cur != a && steps <= self.verts.len() {
            if cur == b {
                return true;
            }
            cur = self.verts[cur].next;
            steps += 1;
        }
        false
    }

    /// Area of the wavefront cycle containing `start`, evaluated at time
    /// `t`. Zero means the cycle has fully collapsed.
    fn cycle_area_at(&self, start: usize, t: f64) -> f64 {
        let mut pts = Vec::new();
        let mut cur = start;
        loop {
            pts.push(self.pos_at(cur, t));
            cur = self.verts[cur].next;
            if cur == start || pts.len() > self.verts.len() {
                break;
            }
        }
        ring_area(&pts)
    }

    fn push_collapse_event(&mut self, va: usize, vb: usize) {
        let (pa, ta, wa) = {
            let v = &self.verts[va];
            let Some(w) = v.velocity else { return };
            (v.pos, v.time, w)
        };
        let (pb, tb, wb) = {
            let v = &self.verts[vb];
            let Some(w) = v.velocity else { return };
            (v.pos, v.time, w)
        };
        let denom = cross_2d(wa, wb);
        if denom.abs() < PARALLEL_EPS {
            return;
        }
        let delta = pb - pa;
        let param_a = cross_2d(delta, wb) / denom;
        let param_b = cross_2d(delta, wa) / denom;
        if param_a < -TIME_EPS || param_b < -TIME_EPS {
            return;
        }
        let point = pa + param_a * wa;
        if !point.x.is_finite() || !point.y.is_finite() {
            return;
        }
        let shared = self.verts[va].right_edge;
        let time = self.edges[shared].offset_dist(point);
        if !time.is_finite() || time < ta - TIME_EPS || time < tb - TIME_EPS {
            return;
        }
        self.queue.push(Reverse(Event {
            time,
            rank: 0,
            order: (va, vb),
            point,
            kind: EventKind::Collapse { va, vb },
        }));
    }

    /// Pushes split-event candidates for reflex vertex `v` against every
    /// non-incident original edge. Candidates are validated lazily at pop
    /// time against the then-current wavefront.
    fn push_split_events(&mut self, v: usize) {
        let (pos, time, w, left, right) = {
            let vert = &self.verts[v];
            let Some(w) = vert.velocity else { return };
            (vert.pos, vert.time, w, vert.left_edge, vert.right_edge)
        };
        for e in 0..self.edges.len() {
            if e == left || e == right {
                continue;
            }
            let edge = self.edges[e];
            let denom = 1.0 - w.dot(&edge.normal);
            if denom.abs() < PARALLEL_EPS {
                continue;
            }
            let d0 = edge.offset_dist(pos);
            let along = (d0 - time) / denom;
            if along <= TIME_EPS || !along.is_finite() {
                continue;
            }
            let hit_time = time + along;
            let hit_point = pos + along * w;
            if !hit_point.x.is_finite() || !hit_point.y.is_finite() || !hit_time.is_finite() {
                continue;
            }
            self.queue.push(Reverse(Event {
                time: hit_time,
                rank: 1,
                order: (v, e),
                point: hit_point,
                kind: EventKind::Split { v, opposite: e },
            }));
        }
    }

    /// Finds the alive wavefront segment derived from original edge `e`
    /// whose extent at time `t` contains `p`. Returns the segment's start
    /// vertex.
    fn find_front(&self, e: usize, p: Point2, t: f64) -> Option<usize> {
        let dir = self.edges[e].dir;
        for (i, vert) in self.verts.iter().enumerate() {
            if !vert.alive || vert.right_edge != e {
                continue;
            }
            let next_id = vert.next;
            if !self.verts[next_id].alive || self.verts[next_id].left_edge != e {
                continue;
            }
            let front_start = self.pos_at(i, t);
            let front_end = self.pos_at(next_id, t);
            if (p - front_start).dot(&dir) >= -FRONT_EPS && (p - front_end).dot(&dir) <= FRONT_EPS {
                return Some(i);
            }
        }
        None
    }

    /// Creates a wavefront vertex at node `node` (time `t`) between edges
    /// `(left, right)` and returns its index. Links are set by the caller.
    fn add_vertex(&mut self, node: usize, t: f64, left: usize, right: usize) -> usize {
        let velocity = self.bisector_velocity(left, right);
        let reflex = velocity.is_some() && self.is_reflex(left, right);
        self.verts.push(Vert {
            pos: self.nodes[node].pos,
            time: t,
            velocity,
            left_edge: left,
            right_edge: right,
            reflex,
            prev: usize::MAX,
            next: usize::MAX,
            alive: true,
            origin_node: node,
        });
        self.verts.len() - 1
    }

    /// Closes a fully collapsed (zero-area) wavefront cycle. The cycle has
    /// degenerated into a tree of segments traversed twice; every distinct
    /// segment becomes one ridge arc bordered by the two faces whose fronts
    /// swept onto it.
    fn close_flat_cycle(&mut self, start: usize, t: f64) -> Result<()> {
        let mut members = Vec::new();
        let mut cur = start;
        loop {
            members.push(cur);
            cur = self.verts[cur].next;
            if cur == start {
                break;
            }
            if members.len() > self.verts.len() {
                return Err(OperationError::Failed(
                    "straight skeleton: broken wavefront cycle".into(),
                )
                .into());
            }
        }
        let mut node_of = Vec::with_capacity(members.len());
        for &v in &members {
            let p = self.pos_at(v, t);
            let node = self.weld(p, t);
            node_of.push(node);
            self.emit_vertex_arc(v, node);
        }
        // Euler tour of the collapsed tree: pair up the two traversals of
        // each distinct segment into one ridge arc. A traversal may pass
        // over nodes welded mid-segment (collinear wavefront vertices), so
        // each half-arc is first subdivided at every cycle node lying on
        // its interior to make the two traversals structurally identical.
        let mut cycle_nodes = node_of.clone();
        cycle_nodes.sort_unstable();
        cycle_nodes.dedup();
        let mut half_arcs: Vec<(usize, usize, usize)> = Vec::new();
        for (k, &v) in members.iter().enumerate() {
            let a = node_of[k];
            let b = node_of[(k + 1) % members.len()];
            if a == b {
                continue;
            }
            let face = self.verts[v].right_edge;
            for (from, to) in self.subdivide_segment(a, b, &cycle_nodes) {
                half_arcs.push((from.min(to), from.max(to), face));
            }
        }
        half_arcs.sort_unstable();
        let mut idx = 0;
        while idx < half_arcs.len() {
            let (a, b, f1) = half_arcs[idx];
            if idx + 1 >= half_arcs.len() {
                return Err(OperationError::Failed(
                    "straight skeleton: unpaired ridge segment in flat cycle".into(),
                )
                .into());
            }
            let (a2, b2, f2) = half_arcs[idx + 1];
            if a2 != a || b2 != b {
                return Err(OperationError::Failed(
                    "straight skeleton: unpaired ridge segment in flat cycle".into(),
                )
                .into());
            }
            if idx + 2 < half_arcs.len() {
                let (a3, b3, _) = half_arcs[idx + 2];
                if a3 == a && b3 == b {
                    return Err(OperationError::Failed(
                        "straight skeleton: overtraversed ridge segment in flat cycle".into(),
                    )
                    .into());
                }
            }
            self.arcs.push(SkelArc {
                from: a,
                to: b,
                faces: [f1, f2],
            });
            idx += 2;
        }
        for &v in &members {
            self.verts[v].alive = false;
        }
        Ok(())
    }

    /// Splits the segment from node `a` to node `b` at every node in
    /// `candidates` lying strictly inside it, returning consecutive
    /// sub-segments in traversal order.
    fn subdivide_segment(&self, a: usize, b: usize, candidates: &[usize]) -> Vec<(usize, usize)> {
        let pa = self.nodes[a].pos;
        let pb = self.nodes[b].pos;
        let d = pb - pa;
        let len = d.norm();
        if len < WELD_EPS {
            return vec![(a, b)];
        }
        let dir = d / len;
        let mut interior: Vec<(f64, usize)> = Vec::new();
        for &c in candidates {
            if c == a || c == b {
                continue;
            }
            let pc = self.nodes[c].pos;
            let along = (pc - pa).dot(&dir);
            if along <= WELD_EPS || along >= len - WELD_EPS {
                continue;
            }
            let offside = (pc - pa - along * dir).norm();
            if offside < WELD_EPS {
                interior.push((along, c));
            }
        }
        if interior.is_empty() {
            return vec![(a, b)];
        }
        interior.sort_by(|x, y| x.0.total_cmp(&y.0));
        let mut result = Vec::with_capacity(interior.len() + 1);
        let mut prev = a;
        for (_, c) in interior {
            result.push((prev, c));
            prev = c;
        }
        result.push((prev, b));
        result
    }

    fn handle_collapse(&mut self, va: usize, vb: usize, p: Point2, t: f64) -> Result<()> {
        if self.cycle_len(va) == 3 || self.cycle_area_at(va, t).abs() < FLAT_AREA_EPS {
            return self.close_flat_cycle(va, t);
        }
        let node = self.weld(p, t);
        self.emit_vertex_arc(va, node);
        self.emit_vertex_arc(vb, node);
        let prev = self.verts[va].prev;
        let next = self.verts[vb].next;
        let left = self.verts[va].left_edge;
        let right = self.verts[vb].right_edge;
        self.verts[va].alive = false;
        self.verts[vb].alive = false;
        let nv = self.add_vertex(node, t, left, right);
        self.verts[nv].prev = prev;
        self.verts[nv].next = next;
        self.verts[prev].next = nv;
        self.verts[next].prev = nv;
        if self.cycle_area_at(nv, t).abs() < FLAT_AREA_EPS {
            return self.close_flat_cycle(nv, t);
        }
        self.push_collapse_event(prev, nv);
        self.push_collapse_event(nv, next);
        if self.verts[nv].reflex {
            self.push_split_events(nv);
        }
        Ok(())
    }

    fn handle_split(&mut self, v: usize, x: usize, p: Point2, t: f64) -> Result<()> {
        let node = self.weld(p, t);
        self.emit_vertex_arc(v, node);
        let front_end = self.verts[x].next;
        let opposite = self.verts[x].right_edge;
        let v_prev = self.verts[v].prev;
        let v_next = self.verts[v].next;
        let v_left = self.verts[v].left_edge;
        let v_right = self.verts[v].right_edge;
        self.verts[v].alive = false;
        // Cycle A: x -> vert_a -> v_next -> ...
        let vert_a = self.add_vertex(node, t, opposite, v_right);
        self.verts[vert_a].prev = x;
        self.verts[vert_a].next = v_next;
        self.verts[x].next = vert_a;
        self.verts[v_next].prev = vert_a;
        // Cycle B: v_prev -> vert_b -> front_end -> ...
        let vert_b = self.add_vertex(node, t, v_left, opposite);
        self.verts[vert_b].prev = v_prev;
        self.verts[vert_b].next = front_end;
        self.verts[v_prev].next = vert_b;
        self.verts[front_end].prev = vert_b;
        for nv in [vert_a, vert_b] {
            if !self.verts[nv].alive {
                continue;
            }
            if self.cycle_len(nv) == 2 || self.cycle_area_at(nv, t).abs() < FLAT_AREA_EPS {
                self.close_flat_cycle(nv, t)?;
                continue;
            }
            let prev = self.verts[nv].prev;
            let next = self.verts[nv].next;
            self.push_collapse_event(prev, nv);
            self.push_collapse_event(nv, next);
            if self.verts[nv].reflex {
                self.push_split_events(nv);
            }
        }
        Ok(())
    }

    fn run(&mut self) -> Result<()> {
        let n = self.edges.len();
        let cap = 64 * n * n + 256;
        let mut processed = 0;
        while let Some(Reverse(event)) = self.queue.pop() {
            processed += 1;
            if processed > cap {
                return Err(
                    OperationError::Failed("straight skeleton did not converge".into()).into(),
                );
            }
            match event.kind {
                EventKind::Collapse { va, vb } => {
                    if !self.verts[va].alive || !self.verts[vb].alive || self.verts[va].next != vb {
                        continue;
                    }
                    self.handle_collapse(va, vb, event.point, event.time)?;
                }
                EventKind::Split { v, opposite } => {
                    if !self.verts[v].alive {
                        continue;
                    }
                    let Some(x) = self.find_front(opposite, event.point, event.time) else {
                        continue;
                    };
                    if !self.same_cycle(v, x) {
                        continue;
                    }
                    self.handle_split(v, x, event.point, event.time)?;
                }
            }
        }
        if self.verts.iter().any(|v| v.alive) {
            return Err(OperationError::Failed(
                "straight skeleton: wavefront did not fully collapse".into(),
            )
            .into());
        }
        Ok(())
    }

    /// Traces the cell of every edge by chaining the arcs bordering it from
    /// the edge's end vertex back to its start vertex.
    fn extract_cells(&self) -> Result<Vec<SkeletonCell>> {
        let n = self.edges.len();
        let mut by_edge: Vec<Vec<usize>> = vec![Vec::new(); n];
        for (i, arc) in self.arcs.iter().enumerate() {
            by_edge[arc.faces[0]].push(i);
            if arc.faces[1] != arc.faces[0] {
                by_edge[arc.faces[1]].push(i);
            }
        }
        let mut cells = Vec::with_capacity(n);
        for (e, edge_arcs) in by_edge.iter().enumerate() {
            let goal = e;
            let start = (e + 1) % n;
            let mut chain = Vec::new();
            let mut used = vec![false; edge_arcs.len()];
            let mut current = start;
            let mut steps = 0;
            while current != goal {
                steps += 1;
                if steps > self.arcs.len() + 2 {
                    return Err(OperationError::Failed(format!(
                        "straight skeleton: cell trace for edge {e} did not close"
                    ))
                    .into());
                }
                let mut found = None;
                for (k, &ai) in edge_arcs.iter().enumerate() {
                    if used[k] {
                        continue;
                    }
                    let arc = &self.arcs[ai];
                    if arc.from == current {
                        found = Some((k, arc.to));
                        break;
                    }
                    if arc.to == current {
                        found = Some((k, arc.from));
                        break;
                    }
                }
                let Some((k, nxt)) = found else {
                    return Err(OperationError::Failed(format!(
                        "straight skeleton: cell trace for edge {e} is disconnected"
                    ))
                    .into());
                };
                used[k] = true;
                current = nxt;
                if current != goal {
                    chain.push(current);
                }
            }
            let mut vertices = Vec::with_capacity(chain.len() + 2);
            for node_id in [goal, start].into_iter().chain(chain) {
                let node = &self.nodes[node_id];
                vertices.push(SkeletonVertex {
                    position: Point3::new(node.pos.x, node.pos.y, 0.0),
                    inset: node.inset,
                });
            }
            cells.push(SkeletonCell {
                edge_index: e,
                vertices,
            });
        }
        Ok(cells)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn ring(pts: &[(f64, f64)]) -> Vec<Point3> {
        pts.iter().map(|&(x, y)| Point3::new(x, y, 0.0)).collect()
    }

    fn cell_area(cell: &SkeletonCell) -> f64 {
        let pts = cell
            .vertices
            .iter()
            .map(|v| Point2::new(v.position.x, v.position.y))
            .collect::<Vec<_>>();
        ring_area(&pts)
    }

    fn assert_common_invariants(skeleton: &StraightSkeleton) {
        let poly = skeleton
            .polygon
            .iter()
            .map(|p| Point2::new(p.x, p.y))
            .collect::<Vec<_>>();
        let polygon_area = ring_area(&poly);
        let mut cells_area = 0.0;
        for cell in &skeleton.cells {
            let area = cell_area(cell);
            assert!(
                area > -1e-9,
                "cell {} is not CCW (area {area})",
                cell.edge_index
            );
            cells_area += area;
            assert!(cell.vertices.len() >= 3);
            // First two vertices are the edge endpoints at inset 0.
            assert!(cell.vertices[0].inset.abs() < 1e-9);
            assert!(cell.vertices[1].inset.abs() < 1e-9);
            for v in &cell.vertices {
                assert!(v.inset >= -1e-9);
                assert!(v.inset <= skeleton.max_inset + 1e-9);
            }
        }
        assert!(
            (cells_area - polygon_area).abs() < 1e-6,
            "cells cover {cells_area}, polygon area {polygon_area}"
        );
    }

    fn find_cell(skeleton: &StraightSkeleton, edge: usize) -> &SkeletonCell {
        skeleton
            .cells
            .iter()
            .find(|c| c.edge_index == edge)
            .unwrap()
    }

    #[test]
    fn square_collapses_to_center_apex() {
        let skeleton =
            compute_straight_skeleton(&ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]))
                .unwrap();
        assert_eq!(skeleton.cells.len(), 4);
        assert!((skeleton.max_inset - 2.0).abs() < 1e-9);
        for cell in &skeleton.cells {
            assert_eq!(cell.vertices.len(), 3, "square cells are triangles");
            let apex = cell.vertices[2];
            assert!((apex.position.x - 2.0).abs() < 1e-9);
            assert!((apex.position.y - 2.0).abs() < 1e-9);
            assert!((apex.inset - 2.0).abs() < 1e-9);
        }
        assert_common_invariants(&skeleton);
    }

    #[test]
    fn rectangle_has_ridge() {
        let skeleton =
            compute_straight_skeleton(&ring(&[(0.0, 0.0), (6.0, 0.0), (6.0, 4.0), (0.0, 4.0)]))
                .unwrap();
        assert!((skeleton.max_inset - 2.0).abs() < 1e-9);
        // Long edges get trapezoids, short edges get triangles.
        let bottom = find_cell(&skeleton, 0);
        assert_eq!(bottom.vertices.len(), 4);
        let right = find_cell(&skeleton, 1);
        assert_eq!(right.vertices.len(), 3);
        // Ridge endpoints at (2,2) and (4,2), inset 2.
        let ridge_pts: Vec<_> = bottom
            .vertices
            .iter()
            .filter(|v| (v.inset - 2.0).abs() < 1e-9)
            .collect();
        assert_eq!(ridge_pts.len(), 2);
        let mut xs = [ridge_pts[0].position.x, ridge_pts[1].position.x];
        xs.sort_by(f64::total_cmp);
        assert!((xs[0] - 2.0).abs() < 1e-9);
        assert!((xs[1] - 4.0).abs() < 1e-9);
        assert_common_invariants(&skeleton);
    }

    #[test]
    fn triangle_apex_is_incenter() {
        // 3-4-5 right triangle: inradius r = (3 + 4 - 5) / 2 = 1,
        // incenter at (1, 1).
        let skeleton =
            compute_straight_skeleton(&ring(&[(0.0, 0.0), (4.0, 0.0), (0.0, 3.0)])).unwrap();
        assert_eq!(skeleton.cells.len(), 3);
        assert!((skeleton.max_inset - 1.0).abs() < 1e-9);
        for cell in &skeleton.cells {
            assert_eq!(cell.vertices.len(), 3);
            let apex = cell.vertices[2];
            assert!((apex.position.x - 1.0).abs() < 1e-9);
            assert!((apex.position.y - 1.0).abs() < 1e-9);
        }
        assert_common_invariants(&skeleton);
    }

    #[test]
    fn l_shape_split_event() {
        // 6x6 square with the top-right 3x3 notch removed; reflex at (3,3).
        let skeleton = compute_straight_skeleton(&ring(&[
            (0.0, 0.0),
            (6.0, 0.0),
            (6.0, 3.0),
            (3.0, 3.0),
            (3.0, 6.0),
            (0.0, 6.0),
        ]))
        .unwrap();
        assert_eq!(skeleton.cells.len(), 6);
        assert!((skeleton.max_inset - 1.5).abs() < 1e-9);
        assert_common_invariants(&skeleton);
        // The reflex vertex projects to the skeleton node (1.5, 1.5).
        let bottom = find_cell(&skeleton, 0);
        assert_eq!(bottom.vertices.len(), 4);
        assert!(bottom
            .vertices
            .iter()
            .any(|v| (v.position.x - 1.5).abs() < 1e-9
                && (v.position.y - 1.5).abs() < 1e-9
                && (v.inset - 1.5).abs() < 1e-9));
        // Cell of the notch bottom edge (6,3)->(3,3).
        let notch = find_cell(&skeleton, 2);
        assert_eq!(notch.vertices.len(), 4);
        assert!(notch
            .vertices
            .iter()
            .any(|v| (v.position.x - 4.5).abs() < 1e-9 && (v.position.y - 1.5).abs() < 1e-9));
    }

    #[test]
    fn t_shape_split_events() {
        // Horizontal bar 10x2 with a 4-wide stem up to y = 6.
        let skeleton = compute_straight_skeleton(&ring(&[
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 2.0),
            (7.0, 2.0),
            (7.0, 6.0),
            (3.0, 6.0),
            (3.0, 2.0),
            (0.0, 2.0),
        ]))
        .unwrap();
        assert_eq!(skeleton.cells.len(), 8);
        assert!((skeleton.max_inset - 2.0).abs() < 1e-9);
        assert_common_invariants(&skeleton);
        // Bottom cell walks the bar ridge and around the stem.
        let bottom = find_cell(&skeleton, 0);
        assert_eq!(bottom.vertices.len(), 7);
    }

    #[test]
    fn plus_shape_fully_symmetric() {
        // Arm width 2, everything collapses at inset 1 simultaneously.
        let skeleton = compute_straight_skeleton(&ring(&[
            (2.0, 0.0),
            (4.0, 0.0),
            (4.0, 2.0),
            (6.0, 2.0),
            (6.0, 4.0),
            (4.0, 4.0),
            (4.0, 6.0),
            (2.0, 6.0),
            (2.0, 4.0),
            (0.0, 4.0),
            (0.0, 2.0),
            (2.0, 2.0),
        ]))
        .unwrap();
        assert_eq!(skeleton.cells.len(), 12);
        assert!((skeleton.max_inset - 1.0).abs() < 1e-9);
        assert_common_invariants(&skeleton);
    }

    #[test]
    fn u_shape_two_reflex_vertices() {
        let skeleton = compute_straight_skeleton(&ring(&[
            (0.0, 0.0),
            (9.0, 0.0),
            (9.0, 6.0),
            (6.0, 6.0),
            (6.0, 2.0),
            (3.0, 2.0),
            (3.0, 6.0),
            (0.0, 6.0),
        ]))
        .unwrap();
        assert_eq!(skeleton.cells.len(), 8);
        assert!((skeleton.max_inset - 1.5).abs() < 1e-9);
        assert_common_invariants(&skeleton);
    }

    #[test]
    fn irregular_convex_pentagon() {
        let skeleton = compute_straight_skeleton(&ring(&[
            (0.0, 0.0),
            (4.0, -1.0),
            (7.0, 2.0),
            (3.5, 5.0),
            (-1.0, 3.0),
        ]))
        .unwrap();
        assert_eq!(skeleton.cells.len(), 5);
        assert!(skeleton.max_inset > 1.0);
        assert_common_invariants(&skeleton);
    }

    #[test]
    fn collinear_mid_edge_vertex() {
        let skeleton = compute_straight_skeleton(&ring(&[
            (0.0, 0.0),
            (3.0, 0.0),
            (6.0, 0.0),
            (6.0, 4.0),
            (0.0, 4.0),
        ]))
        .unwrap();
        assert_eq!(skeleton.cells.len(), 5);
        assert!((skeleton.max_inset - 2.0).abs() < 1e-9);
        assert_common_invariants(&skeleton);
    }

    #[test]
    fn clockwise_input_is_reversed() {
        let skeleton =
            compute_straight_skeleton(&ring(&[(0.0, 4.0), (4.0, 4.0), (4.0, 0.0), (0.0, 0.0)]))
                .unwrap();
        assert_eq!(skeleton.cells.len(), 4);
        assert!((skeleton.max_inset - 2.0).abs() < 1e-9);
        assert_common_invariants(&skeleton);
    }

    #[test]
    fn duplicate_and_closing_points_are_deduped() {
        let skeleton = compute_straight_skeleton(&ring(&[
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (0.0, 0.0),
        ]))
        .unwrap();
        assert_eq!(skeleton.cells.len(), 4);
        assert_common_invariants(&skeleton);
    }

    #[test]
    fn rejects_too_few_vertices() {
        assert!(compute_straight_skeleton(&ring(&[(0.0, 0.0), (1.0, 0.0)])).is_err());
        assert!(compute_straight_skeleton(&[]).is_err());
    }

    #[test]
    fn rejects_non_finite() {
        assert!(compute_straight_skeleton(&[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(f64::NAN, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ])
        .is_err());
    }

    #[test]
    fn rejects_bowtie() {
        assert!(compute_straight_skeleton(&ring(&[
            (0.0, 0.0),
            (4.0, 4.0),
            (4.0, 0.0),
            (0.0, 4.0),
        ]))
        .is_err());
    }

    #[test]
    fn rejects_zero_area() {
        assert!(compute_straight_skeleton(&ring(&[(0.0, 0.0), (2.0, 0.0), (4.0, 0.0),])).is_err());
    }
}
