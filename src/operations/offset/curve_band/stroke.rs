use super::polygon_union::{Polygon, WALL_EPS};
use super::provenance::{CapEnd, OffsetSide};
use crate::operations::offset::pline_offset::RawOffset;

/// Structural origin of one stroke-polygon edge, expressed in the
/// **caller's** vertex frame: `seg` indexes the segments of the
/// `vertices` slice exactly as passed to [`stroke_expand_labeled`]
/// (segment `k` connects `vertices[k]` to `vertices[(k + 1) % n]`),
/// and `side` is relative to the caller's traversal direction.
///
/// A join contributes a single shared point — and hence no edge of its
/// own — for as long as it miters. Past the miter limit it bevels (see
/// [`compute_join`]): it emits BOTH offset endpoints and the chamfer
/// edge between them, which is a `Join` origin. So a stroke edge that
/// bounds material is a `Side` offset, a flat end `Cap`, or a `Join`
/// chamfer. The fourth variant, `Slit`, is the internal cut a closed
/// band is opened along; it never bounds material (see
/// [`stroke_expand_labeled`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrokeOrigin {
    Side {
        seg: usize,
        side: OffsetSide,
    },
    Cap {
        end: CapEnd,
    },
    /// The chamfer that truncates a **bevelled** join, on the given side.
    ///
    /// `seg` names the segment LEAVING the join — the join sits at the
    /// START vertex of segment `seg` — which keys it uniquely for both
    /// open and closed inputs. (Every join has an outgoing segment: an
    /// open input's two free ends are caps, not joins.)
    ///
    /// A chamfer does NOT lie at the distance of the width named by
    /// `side` — it cuts across from segment `seg - 1`'s offset line to
    /// segment `seg`'s, and only its two endpoints are at the width.
    /// That is what separates it from `Side`.
    Join {
        seg: usize,
        side: OffsetSide,
    },
    /// The zero-width slit that opens a closed band's ring at vertex 0
    /// so both offset rings live in ONE polygon (see
    /// [`assemble_closed_keyhole`]). The stroke polygon traverses it
    /// once in each direction, so material lies on BOTH of its sides and
    /// the arrangement always drops it: a `Slit` edge can never reach an
    /// output boundary, and therefore never carries provenance.
    Slit,
}

/// One side's offset chain, grouped by the input vertex that produced
/// each group: `groups[v]` holds the points vertex `v` contributed to
/// this side — ONE for a miter join or a flat cap endpoint, TWO for a
/// bevelled join (the incoming offset endpoint, then the outgoing one).
///
/// Grouping is what lets the assemblers label edges by walking the real
/// per-join point counts instead of deriving them from the segment
/// count: a bevel adds a point, therefore an edge, therefore a label.
type SideGroups = Vec<Vec<(f64, f64)>>;

/// Expands a polyline into a thickened polygon with left/right offsets,
/// reporting per polygon edge the [`StrokeOrigin`] it came from:
/// `labels[e]` describes the edge `poly[e] → poly[(e + 1) % n]`.
///
/// The polygon has two offset "sides", both expressed in the
/// **caller's** traversal frame — the SAME frame the labels use:
///   - left side: points offset by `+left_w` along the left normal of
///     the caller's traversal direction
///   - right side: points offset by `-right_w` along the left normal
///     (i.e. to the caller's right)
///
/// An edge labeled `Side { seg, side }` always lies at the distance of
/// the width parameter named by `side`, for every input winding. A
/// [`StrokeOrigin::Join`] chamfer does not — it truncates a corner
/// whose miter ran past the limit, so only its endpoints are at the
/// width.
///
/// # The result is a fill specification, not a simple polygon
///
/// Both the open and the closed assembly return ONE ring that describes
/// the band as *"the points the ring winds around a non-zero number of
/// times"* — the rule
/// [`point_in_polygon_class`](super::polygon_union::point_in_polygon_class)
/// and hence the union oracle apply. The ring may self-intersect and may
/// touch itself; the arrangement stage flattens it into real topology
/// (dropping every edge with material on both sides). It is deliberately
/// NOT a `PolygonWithHoles` respecting the winding contract: a
/// self-crossing centerline has no such representation.
///
/// - open input: left side forward, end cap, right side back, start cap.
/// - closed input: left ring forward, slit, right ring backward, slit —
///   see [`assemble_closed_keyhole`] for why the band cannot be modelled
///   as an annulus.
pub fn stroke_expand_labeled(
    vertices: &[(f64, f64)],
    closed: bool,
    left_w: f64,
    right_w: f64,
) -> (Polygon, Vec<StrokeOrigin>) {
    let n = vertices.len();
    if n < 2 {
        return (Vec::new(), Vec::new());
    }

    let seg_count = if closed { n } else { n - 1 };
    let dirs: Vec<(f64, f64)> = (0..seg_count)
        .map(|i| {
            let j = (i + 1) % n;
            normalize(vertices[j].0 - vertices[i].0, vertices[j].1 - vertices[i].1)
        })
        .collect();

    if closed {
        let (left, right) = build_closed_offsets(vertices, &dirs, left_w, right_w);
        assemble_closed_keyhole(&left, &right, seg_count)
    } else {
        let (left, right) = build_open_offsets(vertices, &dirs, left_w, right_w);
        assemble_open_ring(&left, &right, seg_count)
    }
}

/// Flattens one side's per-vertex groups into a point run plus the
/// origin of every edge BETWEEN consecutive points of that run.
///
/// Two kinds of edge come out of the walk, and the walk — not the
/// segment count — is what decides how many of each there are:
///
/// - inside group `v`, an edge joins two points the same join emitted,
///   so it is that join's chamfer: [`StrokeOrigin::Join`] with
///   `seg = v` (the segment leaving the join);
/// - between group `v`'s last point and group `v + 1`'s first, the edge
///   runs along segment `v`'s offset line: [`StrokeOrigin::Side`].
///
/// `wrap` closes the run back onto its own first point (closed input);
/// an open run stops at the last group and the caller labels the caps.
///
/// The run stays a polyline either way — `wrap` repeats the first point
/// rather than closing the run — so `labels.len()` is always
/// `pts.len() - 1`, the run's own edge count. The caller supplies the
/// labels of the edges that JOIN the runs (the caps, or the slit).
fn flatten_side(
    groups: &SideGroups,
    side: OffsetSide,
    wrap: bool,
) -> (Vec<(f64, f64)>, Vec<StrokeOrigin>) {
    let total: usize = groups.iter().map(Vec::len).sum();
    let mut pts: Vec<(f64, f64)> = Vec::with_capacity(total + usize::from(wrap));
    let mut labels: Vec<StrokeOrigin> = Vec::with_capacity(total);

    for (v, group) in groups.iter().enumerate() {
        debug_assert!(!group.is_empty(), "every vertex emits at least one point");
        for (k, p) in group.iter().enumerate() {
            if k > 0 {
                labels.push(StrokeOrigin::Join { seg: v, side });
            }
            pts.push(*p);
        }
        if wrap || v + 1 < groups.len() {
            labels.push(StrokeOrigin::Side { seg: v, side });
        }
    }
    if wrap {
        pts.push(groups[0][0]);
    }

    // A run is a polyline either way — `wrap` repeats the first point
    // rather than closing the run — so it always has one more point than
    // it has edges.
    debug_assert_eq!(labels.len(), pts.len() - 1);
    (pts, labels)
}

/// Reverses a flattened run's traversal. Edge `i` of the reversed run is
/// edge `len - 1 - i` of the forward run walked the other way, so it
/// keeps the same origin — reversing the label list is exactly right.
fn reverse_run(
    pts: &[(f64, f64)],
    labels: &[StrokeOrigin],
) -> (Vec<(f64, f64)>, Vec<StrokeOrigin>) {
    (
        pts.iter().rev().copied().collect(),
        labels.iter().rev().copied().collect(),
    )
}

/// Assembles the closed-input band as a single **keyhole** ring: the
/// left offset ring forward, a zero-width slit across the band at
/// vertex 0, the right offset ring backward, and the slit again on the
/// way back.
///
/// # Why not an annulus
///
/// Modelling the band as `outer offset ring + other offset ring as a
/// HOLE` makes the filled set `inside(outer) ∧ ¬inside(hole)`, which is
/// only the band while the centerline ring is **simple**. Let a ring
/// cross itself and both offset rings self-intersect: their signed areas
/// partially cancel (they can be exactly zero), so the larger-area
/// outer/hole pick is arbitrary, and — worse — the "hole" subtracts real
/// band material wherever the two rings' non-zero-winding regions
/// overlap. Whole segments of the input silently lose their material.
///
/// The keyhole ring instead makes the filled set the winding difference
/// `w_left − w_right ≠ 0`, because the slit's two traversals cancel and
/// the two rings are wound in OPPOSITE directions. Every strand of band
/// covering a point contributes exactly `−1` to that difference — a
/// left-hand strand and a right-hand strand of the same crossing add
/// instead of cancelling — so the filled set is exactly the swept band,
/// for simple and self-crossing centerlines alike. It is the same rule
/// the open assembly already relies on, and it needs no winding
/// normalisation of the input: sides and widths stay in the caller's
/// frame throughout.
///
/// The slit is traversed once per direction, so material lies on both
/// its sides and the arrangement drops it; it is never part of an output
/// boundary. Its label is [`StrokeOrigin::Slit`], which carries no
/// provenance — a `Slit` reaching an output boundary is a topology bug
/// and is reported as one.
fn assemble_closed_keyhole(
    left: &SideGroups,
    right: &SideGroups,
    seg_count: usize,
) -> (Polygon, Vec<StrokeOrigin>) {
    // One group per join, and a closed input has a join at every vertex.
    // How many POINTS each group holds is the walk's business, not this
    // assertion's: a bevelled join holds two.
    debug_assert_eq!(left.len(), seg_count);
    debug_assert_eq!(right.len(), seg_count);

    // Left ring forward: group 0 … group n-1, back to group 0's first
    // point. Right ring backward over the same construction.
    let (left_pts, left_labels) = flatten_side(left, OffsetSide::Left, true);
    let (right_fwd_pts, right_fwd_labels) = flatten_side(right, OffsetSide::Right, true);
    let (right_pts, right_labels) = reverse_run(&right_fwd_pts, &right_fwd_labels);

    let mut poly: Polygon = left_pts;
    let mut labels: Vec<StrokeOrigin> = left_labels;
    // Slit L[0] → R[0], across the band at vertex 0.
    labels.push(StrokeOrigin::Slit);
    poly.extend_from_slice(&right_pts);
    labels.extend_from_slice(&right_labels);
    // Closing edge R[0] → L[0]: the slit, retraced.
    labels.push(StrokeOrigin::Slit);

    debug_assert_eq!(labels.len(), poly.len());
    (poly, labels)
}

/// Assembles the open-input single ring: left side forward, end cap,
/// right side backward, start cap — with matching per-edge labels.
fn assemble_open_ring(
    left: &SideGroups,
    right: &SideGroups,
    seg_count: usize,
) -> (Polygon, Vec<StrokeOrigin>) {
    // One group per vertex: a start cap point, one join per interior
    // vertex, an end cap point.
    debug_assert_eq!(left.len(), seg_count + 1);
    debug_assert_eq!(right.len(), seg_count + 1);

    let (left_pts, left_labels) = flatten_side(left, OffsetSide::Left, false);
    let (right_fwd_pts, right_fwd_labels) = flatten_side(right, OffsetSide::Right, false);
    let (right_pts, right_labels) = reverse_run(&right_fwd_pts, &right_fwd_labels);

    let mut poly: Polygon = left_pts;
    let mut labels: Vec<StrokeOrigin> = left_labels;
    labels.push(StrokeOrigin::Cap { end: CapEnd::End });
    poly.extend_from_slice(&right_pts);
    labels.extend_from_slice(&right_labels);
    labels.push(StrokeOrigin::Cap { end: CapEnd::Start });

    debug_assert_eq!(labels.len(), poly.len());
    (poly, labels)
}

fn build_open_offsets(
    verts: &[(f64, f64)],
    dirs: &[(f64, f64)],
    left_w: f64,
    right_w: f64,
) -> (SideGroups, SideGroups) {
    let n = verts.len();
    let seg_count = n - 1;
    let mut left: SideGroups = Vec::with_capacity(n);
    let mut right: SideGroups = Vec::with_capacity(n);

    // First vertex: flat end cap perpendicular to outgoing direction.
    let n0 = left_normal(dirs[0]);
    left.push(vec![offset_point(verts[0], n0, left_w)]);
    right.push(vec![offset_point(verts[0], n0, -right_w)]);

    // Interior vertices: compute join chain.
    for i in 1..seg_count {
        let join = compute_join(verts[i], dirs[i - 1], dirs[i], left_w, right_w);
        left.push(join.left);
        right.push(join.right);
    }

    // Last vertex: flat end cap perpendicular to incoming direction.
    let nl = left_normal(dirs[seg_count - 1]);
    left.push(vec![offset_point(verts[n - 1], nl, left_w)]);
    right.push(vec![offset_point(verts[n - 1], nl, -right_w)]);

    (left, right)
}

fn build_closed_offsets(
    verts: &[(f64, f64)],
    dirs: &[(f64, f64)],
    left_w: f64,
    right_w: f64,
) -> (SideGroups, SideGroups) {
    let seg_count = dirs.len();
    let mut left: SideGroups = Vec::with_capacity(seg_count);
    let mut right: SideGroups = Vec::with_capacity(seg_count);

    for i in 0..seg_count {
        let prev = if i == 0 { seg_count - 1 } else { i - 1 };
        let join = compute_join(verts[i], dirs[prev], dirs[i], left_w, right_w);
        left.push(join.left);
        right.push(join.right);
    }

    (left, right)
}

/// Local chain of offset vertices emitted for a single join.
///
/// A side emits ONE point while its miter stays within the limit — the
/// intersection of the two offset edges, shared by both. Past the limit
/// it emits TWO, the incoming and outgoing offset endpoints, and the
/// chamfer edge between them replaces the spike. When the inner miter
/// flips at very acute concave corners, polygon union downstream removes
/// the inverted region.
struct JoinResult {
    left: Vec<(f64, f64)>,
    right: Vec<(f64, f64)>,
}

/// Computes offset vertices at a single interior join.
///
/// For a CCW polyline, `left_normal` points inward, so:
///   - left side = inner boundary (shrinks at convex corners, reaches out at
///     concave — spiking only for as long as the miter holds)
///   - right side = outer boundary (reaches out at convex corners — spiking
///     only for as long as the miter holds — shrinks at concave)
///
/// `cross = dir_in × dir_out`:
///   - `> 0`: left turn → convex for CCW → outer side is `right`, inner is `left`
///   - `< 0`: right turn → concave for CCW → outer side is `left`, inner is `right`
///   - `≈ 0`: collinear → single offset point, no miter needed
///
/// Each side takes the miter intersection while it stays within
/// [`RawOffset::MITER_LIMIT`] times THAT side's own offset distance
/// (`left_w` for the left, `right_w` for the right), so an asymmetric
/// band can miter on one side and bevel on the other. Past the limit —
/// the outer side of a near-reversal drives the miter point arbitrarily
/// far from the vertex — the corner bevels instead of spiking.
///
/// The limit and its measurement are [`RawOffset`]'s, unchanged: a band
/// stroke and a pline offset of the same corner therefore agree on which
/// corners spike and which are cut.
fn compute_join(
    vertex: (f64, f64),
    dir_in: (f64, f64),
    dir_out: (f64, f64),
    left_w: f64,
    right_w: f64,
) -> JoinResult {
    let cross = dir_in.0 * dir_out.1 - dir_in.1 * dir_out.0;
    let n_in = left_normal(dir_in);
    let n_out = left_normal(dir_out);

    if cross.abs() < WALL_EPS {
        // Collinear: the two offset lines coincide. Emit a single point.
        return JoinResult {
            left: vec![offset_point(vertex, n_out, left_w)],
            right: vec![offset_point(vertex, n_out, -right_w)],
        };
    }

    // Endpoints of the incoming and outgoing offset edges AT the vertex.
    let lp_in = offset_point(vertex, n_in, left_w);
    let lp_out = offset_point(vertex, n_out, left_w);
    let rp_in = offset_point(vertex, n_in, -right_w);
    let rp_out = offset_point(vertex, n_out, -right_w);

    let left_miter = line_intersect(lp_in, dir_in, lp_out, dir_out);
    let right_miter = line_intersect(rp_in, dir_in, rp_out, dir_out);

    JoinResult {
        left: join_chain(left_miter, vertex, left_w, lp_in, lp_out),
        right: join_chain(right_miter, vertex, right_w, rp_in, rp_out),
    }
}

/// Resolves one side of a join to its point chain.
///
/// Takes the miter intersection while it lies within
/// `RawOffset::MITER_LIMIT * |w|` of the SOURCE vertex — the exact test
/// [`RawOffset`] applies to a pline corner, measured the same way, so
/// the two stages agree corner for corner. Past the limit the corner
/// bevels: both offset endpoints are emitted and the chamfer edge
/// between them replaces the spike, which is bounded by construction
/// (each endpoint is exactly `|w|` from the vertex).
///
/// Falls back to the outgoing edge's offset endpoint alone when the two
/// offset directions are numerically parallel (`line_intersect` returned
/// `None`); the single-point fallback preserves continuity.
fn join_chain(
    miter: Option<(f64, f64)>,
    vertex: (f64, f64),
    w: f64,
    p_in: (f64, f64),
    p_out: (f64, f64),
) -> Vec<(f64, f64)> {
    let Some(m) = miter else {
        return vec![p_out];
    };
    let dx = m.0 - vertex.0;
    let dy = m.1 - vertex.1;
    let miter_dist_sq = dx * dx + dy * dy;
    let limit = RawOffset::MITER_LIMIT * w.abs();
    if miter_dist_sq > limit * limit {
        vec![p_in, p_out]
    } else {
        vec![m]
    }
}

fn left_normal(dir: (f64, f64)) -> (f64, f64) {
    (-dir.1, dir.0)
}

fn offset_point(p: (f64, f64), normal: (f64, f64), w: f64) -> (f64, f64) {
    (p.0 + w * normal.0, p.1 + w * normal.1)
}

fn normalize(dx: f64, dy: f64) -> (f64, f64) {
    let len = (dx * dx + dy * dy).sqrt();
    if len < WALL_EPS {
        (1.0, 0.0)
    } else {
        (dx / len, dy / len)
    }
}

fn line_intersect(
    p1: (f64, f64),
    d1: (f64, f64),
    p2: (f64, f64),
    d2: (f64, f64),
) -> Option<(f64, f64)> {
    let cross = d1.0 * d2.1 - d1.1 * d2.0;
    if cross.abs() < WALL_EPS * WALL_EPS {
        return None;
    }
    let dx = p2.0 - p1.0;
    let dy = p2.1 - p1.1;
    let t = (dx * d2.1 - dy * d2.0) / cross;
    Some((p1.0 + t * d1.0, p1.1 + t * d1.1))
}

#[cfg(test)]
mod tests {
    use super::*;

    use super::super::polygon_union::{point_in_polygon_class, PointClass};
    use crate::geometry::pline::{Pline, PlineVertex};

    /// Geometry-only view for tests written against the pre-labels API.
    fn stroke_expand(vertices: &[(f64, f64)], closed: bool, left_w: f64, right_w: f64) -> Polygon {
        stroke_expand_labeled(vertices, closed, left_w, right_w).0
    }

    fn signed_area_tuples(verts: &[(f64, f64)]) -> f64 {
        let n = verts.len();
        let mut area = 0.0;
        for i in 0..n {
            let j = (i + 1) % n;
            area += verts[i].0 * verts[j].1;
            area -= verts[j].0 * verts[i].1;
        }
        area * 0.5
    }

    /// The stroke ring is a fill specification under the non-zero
    /// winding rule — the same rule the union oracle applies — so "is
    /// this point band material?" is exactly `Inside`.
    fn fills(ring: &Polygon, p: (f64, f64)) -> bool {
        point_in_polygon_class(p, ring) == PointClass::Inside
    }

    #[test]
    fn straight_open() {
        let result = stroke_expand(&[(0.0, 0.0), (5.0, 0.0)], false, 0.3, 0.3);
        let area = signed_area_tuples(&result).abs();
        let expected = 5.0 * 0.6;
        assert!(
            (area - expected).abs() < 0.1,
            "area={area}, expected={expected}"
        );
    }

    #[test]
    fn l_shape_open() {
        let result = stroke_expand(&[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0)], false, 0.3, 0.3);
        let area = signed_area_tuples(&result).abs();
        assert!(area > 3.0, "area={area} too small");
    }

    #[test]
    fn closed_square() {
        let square = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)];
        let result = stroke_expand(&square, true, 0.3, 0.3);
        // Keyhole ring: one edge per segment per side, plus the slit
        // traversed twice.
        assert_eq!(result.len(), 2 * square.len() + 2);
        // Band material on every wall, nothing in the room, nothing out.
        for p in [(5.0, 0.0), (10.0, 5.0), (5.0, 10.0), (0.0, 5.0)] {
            assert!(fills(&result, p), "wall centre {p:?} must be material");
        }
        assert!(!fills(&result, (5.0, 5.0)), "the room must stay empty");
        assert!(!fills(&result, (-1.0, 5.0)), "outside must stay empty");
        // Winding difference over the annulus is a full turn of area.
        let area = signed_area_tuples(&result).abs();
        assert!(area > 15.0 && area < 30.0, "band area={area}");
    }

    #[test]
    fn closed_square_cw_input() {
        // Winding of the input is irrelevant: the keyhole ring never
        // normalises it, and the band is the same set either way.
        let cw = [(0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0)];
        let result = stroke_expand(&cw, true, 0.3, 0.3);
        for p in [(5.0, 0.0), (10.0, 5.0), (5.0, 10.0), (0.0, 5.0)] {
            assert!(fills(&result, p), "wall centre {p:?} must be material");
        }
        assert!(!fills(&result, (5.0, 5.0)), "the room must stay empty");
        let ccw = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)];
        assert!(
            (signed_area_tuples(&result).abs()
                - signed_area_tuples(&stroke_expand(&ccw, true, 0.3, 0.3)).abs())
            .abs()
                < 1e-9,
            "CW and CCW inputs must sweep the same band"
        );
    }

    #[test]
    fn asymmetric_offset() {
        let result = stroke_expand(&[(0.0, 0.0), (5.0, 0.0)], false, 0.0, 0.3);
        let area = signed_area_tuples(&result).abs();
        let expected = 5.0 * 0.3;
        assert!(
            (area - expected).abs() < 0.1,
            "area={area}, expected={expected}"
        );
    }

    #[test]
    fn acute_angle_keeps_sharp_miter() {
        // Acute (45°) turn, still INSIDE the miter limit — its miter
        // reaches 0.784, under `MITER_LIMIT * 0.3 = 1.2` — so the corner
        // stays sharp. Result must be non-degenerate.
        let result = stroke_expand(&[(0.0, 0.0), (5.0, 0.0), (4.9, 0.1)], false, 0.3, 0.3);
        assert!(result.len() >= 3);
    }

    #[test]
    // The bbox extremes are clearest as input_{min,max}_{x,y}.
    #[allow(clippy::similar_names)]
    fn angled_closed_polygon_no_dent() {
        // Regression test: diagonal (non-axis-aligned) closed polygon.
        // Earlier iterations produced corner "dents" on the outer boundary
        // when the polygon was not axis-aligned.
        // Input is a CCW quadrilateral with all oblique edges.
        let quad = [
            (-3.217, -4.144),
            (2.002, -4.631),
            (2.578, 1.534),
            (-2.635, 2.085),
        ];
        let result = stroke_expand(&quad, true, 0.15, 0.15);

        // Every ring vertex except the inner offsets should reach beyond
        // the input extents (never dented inward past an input vertex).
        let input_max_x = 2.578_f64;
        let input_min_x = -3.217_f64;
        let input_max_y = 2.085_f64;
        let input_min_y = -4.631_f64;
        let mut saw_x_beyond = false;
        let mut saw_y_beyond = false;
        for &(x, y) in &result {
            if x > input_max_x || x < input_min_x {
                saw_x_beyond = true;
            }
            if y > input_max_y || y < input_min_y {
                saw_y_beyond = true;
            }
        }
        assert!(saw_x_beyond, "band must extend beyond input X extent");
        assert!(saw_y_beyond, "band must extend beyond input Y extent");

        // Every edge midpoint of the centerline carries material, and the
        // room interior does not.
        for k in 0..quad.len() {
            let a = quad[k];
            let b = quad[(k + 1) % quad.len()];
            let mid = ((a.0 + b.0) * 0.5, (a.1 + b.1) * 0.5);
            assert!(
                fills(&result, mid),
                "segment {k} midpoint carries no material"
            );
        }
        assert!(!fills(&result, (0.0, -1.0)), "the room must stay empty");

        // Thickness ≈ 0.3 uniformly: band area is roughly perimeter * 0.3.
        let ring_area = signed_area_tuples(&result).abs();
        assert!(ring_area > 4.0 && ring_area < 12.0, "ring_area={ring_area}");
    }

    /// The band of a SELF-CROSSING closed ring keeps material on every
    /// segment. Modelled as an annulus (outer offset ring minus the other
    /// as a hole) it does not: both offset rings self-intersect, their
    /// non-zero regions overlap, and the "hole" eats whole segments.
    #[test]
    fn closed_self_crossing_ring_fills_every_segment() {
        let ring = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (2.0, -2.0)];
        let result = stroke_expand(&ring, true, 0.09, 0.09);
        for k in 0..ring.len() {
            let a = ring[k];
            let b = ring[(k + 1) % ring.len()];
            let mid = ((a.0 + b.0) * 0.5, (a.1 + b.1) * 0.5);
            assert!(
                fills(&result, mid),
                "segment {k} midpoint {mid:?} carries no material"
            );
        }
        // The crossing itself (segment 2 over segment 0 at x ≈ 2.667) is
        // covered twice over — two strands ADD, they must not cancel.
        assert!(
            fills(&result, (2.6667, 0.0)),
            "the crossing must be material"
        );
    }

    /// A figure-8 ring's two lobes wind in OPPOSITE directions, so its
    /// signed area is exactly zero — the degenerate case that made the
    /// annulus model's larger-area outer/hole pick arbitrary.
    #[test]
    fn closed_figure_8_fills_every_segment() {
        let ring = [(0.0, 0.0), (2.0, 2.0), (0.0, 2.0), (2.0, 0.0)];
        assert!(
            signed_area_tuples(&ring).abs() < 1e-12,
            "fixture must be the zero-signed-area case"
        );
        let result = stroke_expand(&ring, true, 0.09, 0.09);
        for k in 0..ring.len() {
            let a = ring[k];
            let b = ring[(k + 1) % ring.len()];
            let mid = ((a.0 + b.0) * 0.5, (a.1 + b.1) * 0.5);
            assert!(
                fills(&result, mid),
                "segment {k} midpoint {mid:?} carries no material"
            );
        }
        assert!(fills(&result, (1.0, 1.0)), "the crossing must be material");
        // The lobes themselves stay empty.
        assert!(!fills(&result, (1.0, 1.6)), "upper lobe must be empty");
        assert!(!fills(&result, (1.0, 0.4)), "lower lobe must be empty");
    }

    /// The closed ring's slit is emitted twice, in opposite directions,
    /// over the same two points — the property that makes the
    /// arrangement drop it (material on both sides) and that keeps it
    /// out of every output boundary.
    #[test]
    fn closed_slit_is_traversed_once_per_direction() {
        let square = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)];
        let (ring, labels) = stroke_expand_labeled(&square, true, 0.3, 0.3);
        let slits: Vec<usize> = labels
            .iter()
            .enumerate()
            .filter(|(_, l)| matches!(l, StrokeOrigin::Slit))
            .map(|(e, _)| e)
            .collect();
        assert_eq!(slits.len(), 2, "exactly one slit, traversed twice");
        let edge = |e: usize| (ring[e], ring[(e + 1) % ring.len()]);
        let (a0, b0) = edge(slits[0]);
        let (a1, b1) = edge(slits[1]);
        assert_eq!((a0, b0), (b1, a1), "the two traversals must be opposite");
        // It spans the band at vertex 0, from the left miter to the right
        // one — for this CCW square, (0.3, 0.3) to (−0.3, −0.3).
        assert!(
            (a0.0 - 0.3).abs() < 1e-12 && (a0.1 - 0.3).abs() < 1e-12,
            "{a0:?}"
        );
        assert!(
            (b0.0 + 0.3).abs() < 1e-12 && (b0.1 + 0.3).abs() < 1e-12,
            "{b0:?}"
        );
    }

    #[test]
    fn join_result_convex_miter() {
        // Convex CCW corner under miter limit → single-point miters on both
        // sides. Verifies compute_join's local chain contract.
        let join = compute_join((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), 0.3, 0.3);
        assert_eq!(join.left.len(), 1, "inner miter is a single point");
        assert_eq!(
            join.right.len(),
            1,
            "outer miter under limit is a single point"
        );
    }

    #[test]
    fn join_result_collinear_single_point() {
        let join = compute_join((1.0, 0.0), (1.0, 0.0), (1.0, 0.0), 0.3, 0.3);
        assert_eq!(join.left.len(), 1);
        assert_eq!(join.right.len(), 1);
    }

    // ===== the miter limit =====

    fn dist_from(p: (f64, f64), origin: (f64, f64)) -> f64 {
        ((p.0 - origin.0).powi(2) + (p.1 - origin.1).powi(2)).sqrt()
    }

    fn assert_pt_near(actual: (f64, f64), expected: (f64, f64), what: &str) {
        assert!(
            (actual.0 - expected.0).abs() < 1e-12 && (actual.1 - expected.1).abs() < 1e-12,
            "{what}: got {actual:?}, expected {expected:?}"
        );
    }

    /// Unit direction from `a` to `b`.
    fn dir_between(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
        normalize(b.0 - a.0, b.1 - a.1)
    }

    /// A near-reversal (hairpin) join drove the miter point 12 m from the
    /// vertex of a 0.3 m band — the unbounded spike this limit exists to
    /// cut. Both sides now bevel to their two offset endpoints, which sit
    /// at exactly `w` from the vertex, so nothing the join emits can
    /// reach past `MITER_LIMIT * w`.
    #[test]
    fn join_near_reversal_bevels_to_the_two_offset_endpoints() {
        let vertex = (0.0, 0.0);
        let dir_in = (1.0, 0.0);
        let dir_out = normalize(-1.0, 0.05);
        let w = 0.3;

        let n_in = left_normal(dir_in);
        let n_out = left_normal(dir_out);
        let lp_in = offset_point(vertex, n_in, w);
        let lp_out = offset_point(vertex, n_out, w);
        let rp_in = offset_point(vertex, n_in, -w);
        let rp_out = offset_point(vertex, n_out, -w);

        // The fixture really is a runaway miter, not a merely sharp one.
        let Some(unbounded) = line_intersect(lp_in, dir_in, lp_out, dir_out) else {
            panic!("the hairpin's offset lines must still meet");
        };
        let spike = dist_from(unbounded, vertex);
        assert!(spike > 10.0, "fixture must run away; miter reached {spike}");

        let join = compute_join(vertex, dir_in, dir_out, w, w);

        assert_eq!(join.left.len(), 2, "left side must bevel");
        assert_pt_near(join.left[0], lp_in, "left bevel start");
        assert_pt_near(join.left[1], lp_out, "left bevel end");
        assert_eq!(join.right.len(), 2, "right side must bevel");
        assert_pt_near(join.right[0], rp_in, "right bevel start");
        assert_pt_near(join.right[1], rp_out, "right bevel end");

        let limit = RawOffset::MITER_LIMIT * w;
        for p in join.left.iter().chain(&join.right) {
            let d = dist_from(*p, vertex);
            assert!(
                d <= limit,
                "join emitted {p:?} at {d} from the vertex, past the {limit} limit"
            );
        }
    }

    /// Each side is measured against ITS OWN width, and the sides can
    /// therefore reach opposite verdicts at one corner.
    ///
    /// A line-line corner's miter reaches `w / sin(φ/2)` against a limit
    /// of `4 w`, so the ratio is `1 / (4 sin(φ/2))` — width-free, and the
    /// two sides normally agree. The exception is a side of width ZERO:
    /// its "offset" is the centerline itself, its corner IS the vertex,
    /// and it has no spike to cut. It must stay a single point while the
    /// wide side of the same hairpin bevels.
    #[test]
    fn join_zero_width_side_never_bevels() {
        let vertex = (0.0, 0.0);
        let dir_in = (1.0, 0.0);
        let dir_out = normalize(-1.0, 0.05);
        let join = compute_join(vertex, dir_in, dir_out, 0.0, 0.3);

        assert_eq!(join.left.len(), 1, "a zero-width side has no spike to cut");
        assert_pt_near(join.left[0], vertex, "zero-width miter");
        assert_eq!(join.right.len(), 2, "the 0.3-wide side still bevels");
    }

    /// Regression guard: a moderate corner is untouched. 60° between the
    /// legs puts the miter at `w / sin(30°) = 2 w`, half the limit, so
    /// both sides still emit the single shared point they always did.
    #[test]
    fn join_moderate_corner_still_miters() {
        let w = 0.3;
        let vertex = (0.0, 0.0);
        let dir_in = (1.0, 0.0);
        let dir_out = normalize(-0.5, 3.0_f64.sqrt() * 0.5);
        let join = compute_join(vertex, dir_in, dir_out, w, w);
        assert_eq!(join.left.len(), 1, "left side keeps one miter point");
        assert_eq!(join.right.len(), 1, "right side keeps one miter point");
        for p in join.left.iter().chain(&join.right) {
            let d = dist_from(*p, vertex);
            assert!(
                (d - 2.0 * w).abs() < 1e-9,
                "miter reached {d}, expected {}",
                2.0 * w
            );
        }
    }

    /// A 20° corner: past the `4 w` miter limit, and short of
    /// `RawOffset`'s own near-antiparallel flat-cap shortcut, so the
    /// LIMIT is what decides in both stages.
    const HAIRPIN: [(f64, f64); 3] = [
        (-3.0, 0.0),
        (0.0, 0.0),
        (-2.819_077_862_357_725, 1.026_060_429_977_006),
    ];
    /// The same shape opened out to 60°: miter `2 w`, well under the limit.
    const MODERATE: [(f64, f64); 3] = [(-3.0, 0.0), (0.0, 0.0), (-1.5, 2.598_076_211_353_316)];

    /// The limit and the way it is measured are [`RawOffset`]'s, so the
    /// band stroke and a pline offset of the SAME corner reach the same
    /// verdict — that is what "band strokes and pline offsets agree"
    /// means, and it is checked on both sides of the threshold.
    #[test]
    fn band_stroke_and_raw_offset_agree_on_bevelling() {
        let w = 0.3;
        for (verts, expect_bevel) in [(HAIRPIN, true), (MODERATE, false)] {
            let pline = Pline {
                vertices: verts
                    .iter()
                    .map(|&(x, y)| PlineVertex::line(x, y))
                    .collect(),
                closed: false,
            };
            // Both offset sides of the pline, which is what the band's
            // two sides correspond to.
            for signed in [w, -w] {
                let Ok(raw) = RawOffset::build(&pline, signed) else {
                    panic!("raw offset of {verts:?} at {signed} must build");
                };
                assert_eq!(
                    raw.segments[1].bevelled_entry, expect_bevel,
                    "RawOffset at distance {signed} disagrees for {verts:?}"
                );
            }

            let join = compute_join(
                verts[1],
                dir_between(verts[0], verts[1]),
                dir_between(verts[1], verts[2]),
                w,
                w,
            );
            let expected = if expect_bevel { 2 } else { 1 };
            assert_eq!(join.left.len(), expected, "band left side, {verts:?}");
            assert_eq!(join.right.len(), expected, "band right side, {verts:?}");
        }
    }

    /// A bevelled join adds a point, therefore an edge, therefore a
    /// label. The OPEN assembly builds its labels by walking the real
    /// per-join point counts, so the 1:1 edge↔label alignment survives —
    /// and the chamfer is a `Join`, keyed by the segment LEAVING the
    /// corner, never a `Side` (which would claim it lies at the width).
    #[test]
    fn open_assembly_labels_a_bevelled_join_without_losing_alignment() {
        let w = 0.3;
        let (ring, labels) = stroke_expand_labeled(&HAIRPIN, false, w, w);
        assert_eq!(labels.len(), ring.len(), "one label per ring edge");

        let l = |seg| StrokeOrigin::Side {
            seg,
            side: OffsetSide::Left,
        };
        let r = |seg| StrokeOrigin::Side {
            seg,
            side: OffsetSide::Right,
        };
        assert_eq!(
            labels,
            vec![
                l(0),
                StrokeOrigin::Join {
                    seg: 1,
                    side: OffsetSide::Left
                },
                l(1),
                StrokeOrigin::Cap { end: CapEnd::End },
                r(1),
                StrokeOrigin::Join {
                    seg: 1,
                    side: OffsetSide::Right
                },
                r(0),
                StrokeOrigin::Cap { end: CapEnd::Start },
            ]
        );
        // Every `Side` edge still lies on its own offset line — the
        // chamfers are exactly the edges that do not, which is why they
        // are labelled apart.
        assert_side_edges_on_offset_lines(&ring, &labels, &HAIRPIN, w, w);
        // The bevel's own points sit at exactly `w` from the corner —
        // the miter they replaced reached 1.73.
        for p in [ring[1], ring[2], ring[5], ring[6]] {
            let d = dist_from(p, HAIRPIN[1]);
            assert!(
                (d - w).abs() < 1e-12,
                "bevel point {p:?} is at {d}, not {w}"
            );
        }
    }

    /// The CLOSED keyhole assembly does the same walk, on both rings, and
    /// the right ring's labels stay aligned through its reversal.
    #[test]
    fn closed_assembly_labels_bevelled_joins_without_losing_alignment() {
        // A needle triangle: the two tip corners (3.4° between their
        // legs, miter 3.34) bevel; the blunt one (173°, miter 0.10) does
        // not.
        let needle = [(0.0, 0.0), (10.0, 0.0), (5.0, 0.3)];
        let w = 0.1;
        let (ring, labels) = stroke_expand_labeled(&needle, true, w, w);
        assert_eq!(labels.len(), ring.len(), "one label per ring edge");

        let jl = |seg| StrokeOrigin::Join {
            seg,
            side: OffsetSide::Left,
        };
        let jr = |seg| StrokeOrigin::Join {
            seg,
            side: OffsetSide::Right,
        };
        let l = |seg| StrokeOrigin::Side {
            seg,
            side: OffsetSide::Left,
        };
        let r = |seg| StrokeOrigin::Side {
            seg,
            side: OffsetSide::Right,
        };
        assert_eq!(
            labels,
            vec![
                jl(0),
                l(0),
                jl(1),
                l(1),
                l(2),
                StrokeOrigin::Slit,
                r(2),
                r(1),
                jr(1),
                r(0),
                jr(0),
                StrokeOrigin::Slit,
            ]
        );
        assert_side_edges_on_offset_lines(&ring, &labels, &needle, w, w);
        // Unbounded, the two tip miters reached 3.34 and the ring spanned
        // x ∈ [−3.24, 13.24]. Bounded, every ring point is a bevel
        // endpoint (exactly `w` from its vertex) or the blunt corner's
        // miter (0.1002 from its own), so the ring stays inside the
        // needle's own box grown by 0.11.
        let pad = 0.11;
        for p in &ring {
            assert!(
                p.0 >= -pad && p.0 <= 10.0 + pad && p.1 >= -pad && p.1 <= 0.3 + pad,
                "ring point {p:?} spiked outside the padded needle box"
            );
        }
    }

    // ===== stroke_expand_labeled — per-edge origin labels =====

    #[test]
    fn labels_open_l_shape_sequence() {
        let (ring, labels) =
            stroke_expand_labeled(&[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0)], false, 0.3, 0.3);
        assert_eq!(labels.len(), ring.len());
        let l = |seg| StrokeOrigin::Side {
            seg,
            side: OffsetSide::Left,
        };
        let r = |seg| StrokeOrigin::Side {
            seg,
            side: OffsetSide::Right,
        };
        assert_eq!(
            labels,
            vec![
                l(0),
                l(1),
                StrokeOrigin::Cap { end: CapEnd::End },
                r(1),
                r(0),
                StrokeOrigin::Cap { end: CapEnd::Start },
            ]
        );
    }

    /// The closed ring's label sequence: every segment's LEFT offset in
    /// order, the slit, every segment's RIGHT offset in reverse order,
    /// the slit again.
    #[test]
    fn labels_closed_square_sequence() {
        let square = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)];
        let (ring, labels) = stroke_expand_labeled(&square, true, 0.3, 0.3);
        assert_eq!(labels.len(), ring.len());
        let l = |seg| StrokeOrigin::Side {
            seg,
            side: OffsetSide::Left,
        };
        let r = |seg| StrokeOrigin::Side {
            seg,
            side: OffsetSide::Right,
        };
        assert_eq!(
            labels,
            vec![
                l(0),
                l(1),
                l(2),
                l(3),
                StrokeOrigin::Slit,
                r(3),
                r(2),
                r(1),
                r(0),
                StrokeOrigin::Slit,
            ]
        );
    }

    /// Asserts that ring edge `e` (labelled `Side { seg, side }`) lies on
    /// the offset supporting line of centerline segment `seg`.
    fn assert_side_edges_on_offset_lines(
        ring: &[(f64, f64)],
        labels: &[StrokeOrigin],
        centerline: &[(f64, f64)],
        left_w: f64,
        right_w: f64,
    ) {
        let vert_count = centerline.len();
        for (e, label) in labels.iter().enumerate() {
            let StrokeOrigin::Side { seg, side } = *label else {
                continue;
            };
            let seg_a = centerline[seg];
            let seg_b = centerline[(seg + 1) % vert_count];
            let dir = normalize(seg_b.0 - seg_a.0, seg_b.1 - seg_a.1);
            let nn = left_normal(dir);
            let width = match side {
                OffsetSide::Left => left_w,
                OffsetSide::Right => -right_w,
            };
            let base = (seg_a.0 + width * nn.0, seg_a.1 + width * nn.1);
            for p in [ring[e], ring[(e + 1) % ring.len()]] {
                let perp = (p.0 - base.0) * dir.1 - (p.1 - base.1) * dir.0;
                assert!(
                    perp.abs() < 1e-9,
                    "edge {e} labelled seg={seg} side={side:?} not on its \
                     offset line: perp={perp}"
                );
            }
        }
    }

    #[test]
    fn labels_closed_square_sides_are_on_their_offset_lines() {
        let square = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)];
        let (ring, labels) = stroke_expand_labeled(&square, true, 0.3, 0.3);
        assert_eq!(labels.len(), ring.len());
        // Every centerline segment appears exactly once per side; no caps.
        for side in [OffsetSide::Left, OffsetSide::Right] {
            let mut segs: Vec<usize> = labels
                .iter()
                .filter_map(|l| match l {
                    StrokeOrigin::Side { seg, side: s } if *s == side => Some(*seg),
                    StrokeOrigin::Cap { .. } => panic!("closed input must not emit caps"),
                    _ => None,
                })
                .collect();
            segs.sort_unstable();
            assert_eq!(segs, vec![0, 1, 2, 3], "{side:?}");
        }
        assert_side_edges_on_offset_lines(&ring, &labels, &square, 0.3, 0.3);
    }

    /// Sides are reported relative to the caller's traversal direction,
    /// so the same square traversed CW puts the LEFT side outward.
    #[test]
    fn labels_closed_square_cw_input_sides_follow_caller_direction() {
        let square_cw = [(0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0)];
        let (ring, labels) = stroke_expand_labeled(&square_cw, true, 0.3, 0.3);
        assert_side_edges_on_offset_lines(&ring, &labels, &square_cw, 0.3, 0.3);
        // CW ring: caller-left is the ring exterior.
        let leftmost = |side: OffsetSide| {
            labels
                .iter()
                .enumerate()
                .filter(|(_, l)| matches!(l, StrokeOrigin::Side { side: s, .. } if *s == side))
                .flat_map(|(e, _)| [ring[e].0, ring[(e + 1) % ring.len()].0])
                .fold(f64::INFINITY, f64::min)
        };
        assert!((leftmost(OffsetSide::Left) + 0.3).abs() < 1e-9);
        assert!((leftmost(OffsetSide::Right) - 0.3).abs() < 1e-9);
    }

    /// A CW closed ring with ASYMMETRIC widths: the widths follow the
    /// caller frame exactly like the labels do, so an edge labelled
    /// `Left` lies at `left_w` on the caller's left.
    #[test]
    fn labels_closed_square_cw_asymmetric_widths_follow_caller_sides() {
        let square_cw = [(0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0)];
        let (ring, labels) = stroke_expand_labeled(&square_cw, true, 0.4, 0.1);
        assert_side_edges_on_offset_lines(&ring, &labels, &square_cw, 0.4, 0.1);
        // CW ring: caller-left is the ring exterior, so the 0.4 band
        // grows outward — the band reaches x = −0.4.
        let min_x = ring.iter().map(|p| p.0).fold(f64::INFINITY, f64::min);
        assert!(
            (min_x + 0.4).abs() < 1e-9,
            "caller-left width 0.4 must expand the CW ring outward to \
             x = −0.4; got min x = {min_x}"
        );
    }

    #[test]
    fn labels_asymmetric_open_segment_on_offset_lines() {
        let line = [(0.0, 0.0), (5.0, 0.0)];
        let (ring, labels) = stroke_expand_labeled(&line, false, 0.0, 0.3);
        assert_eq!(labels.len(), ring.len());
        assert_side_edges_on_offset_lines(&ring, &labels, &line, 0.0, 0.3);
    }
}
