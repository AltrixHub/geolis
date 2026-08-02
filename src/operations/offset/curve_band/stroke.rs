use super::polygon_union::{Polygon, WALL_EPS};
use super::provenance::{CapEnd, OffsetSide};

/// Structural origin of one stroke-polygon edge, expressed in the
/// **caller's** vertex frame: `seg` indexes the segments of the
/// `vertices` slice exactly as passed to [`stroke_expand_labeled`]
/// (segment `k` connects `vertices[k]` to `vertices[(k + 1) % n]`),
/// and `side` is relative to the caller's traversal direction.
///
/// Joins are miters (each join contributes a single shared point, never
/// its own edge), so every stroke edge that bounds material is either a
/// `Side` offset or a flat end `Cap` — there is no join-arc origin. The
/// third variant, `Slit`, is the internal cut a closed band is opened
/// along; it never bounds material (see [`stroke_expand_labeled`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrokeOrigin {
    Side {
        seg: usize,
        side: OffsetSide,
    },
    Cap {
        end: CapEnd,
    },
    /// The zero-width slit that opens a closed band's ring at vertex 0
    /// so both offset rings live in ONE polygon (see
    /// [`assemble_closed_keyhole`]). The stroke polygon traverses it
    /// once in each direction, so material lies on BOTH of its sides and
    /// the arrangement always drops it: a `Slit` edge can never reach an
    /// output boundary, and therefore never carries provenance.
    Slit,
}

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
/// the width parameter named by `side`, for every input winding.
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

    let mut left_pts: Vec<(f64, f64)> = Vec::new();
    let mut right_pts: Vec<(f64, f64)> = Vec::new();

    if closed {
        build_closed_offsets(
            vertices,
            &dirs,
            left_w,
            right_w,
            &mut left_pts,
            &mut right_pts,
        );
        assemble_closed_keyhole(&left_pts, &right_pts, seg_count)
    } else {
        build_open_offsets(
            vertices,
            &dirs,
            left_w,
            right_w,
            &mut left_pts,
            &mut right_pts,
        );
        assemble_open_ring(&left_pts, right_pts, seg_count)
    }
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
    left_pts: &[(f64, f64)],
    right_pts: &[(f64, f64)],
    seg_count: usize,
) -> (Polygon, Vec<StrokeOrigin>) {
    // Each join contributes exactly one miter point per side, so both
    // rings have exactly seg_count points; ring edge k runs along the
    // offset line of segment k.
    debug_assert_eq!(left_pts.len(), seg_count);
    debug_assert_eq!(right_pts.len(), seg_count);

    let mut poly: Polygon = Vec::with_capacity(2 * seg_count + 2);
    let mut labels: Vec<StrokeOrigin> = Vec::with_capacity(2 * seg_count + 2);

    // Left ring forward: L[0] → … → L[n-1] → L[0]. Edge k lies on
    // segment k's left offset line.
    poly.extend_from_slice(left_pts);
    poly.push(left_pts[0]);
    for k in 0..seg_count {
        labels.push(StrokeOrigin::Side {
            seg: k,
            side: OffsetSide::Left,
        });
    }
    // Slit L[0] → R[0], across the band at vertex 0.
    labels.push(StrokeOrigin::Slit);

    // Right ring backward: R[0] → R[n-1] → … → R[1] → R[0]. The edge
    // leaving R[j] retraces segment `j - 1`'s right offset line.
    poly.push(right_pts[0]);
    for j in 1..seg_count {
        poly.push(right_pts[seg_count - j]);
    }
    poly.push(right_pts[0]);
    for j in 0..seg_count {
        labels.push(StrokeOrigin::Side {
            seg: seg_count - j - 1,
            side: OffsetSide::Right,
        });
    }
    // Closing edge R[0] → L[0]: the slit, retraced.
    labels.push(StrokeOrigin::Slit);

    debug_assert_eq!(poly.len(), 2 * seg_count + 2);
    debug_assert_eq!(labels.len(), poly.len());
    (poly, labels)
}

/// Assembles the open-input single ring: left side forward, end cap,
/// right side backward, start cap — with matching per-edge labels.
fn assemble_open_ring(
    left_pts: &[(f64, f64)],
    mut right_pts: Vec<(f64, f64)>,
    seg_count: usize,
) -> (Polygon, Vec<StrokeOrigin>) {
    // Each side has seg_count + 1 points (start cap point, one miter per
    // interior vertex, end point); side edge k lies on the offset line
    // of segment k.
    debug_assert_eq!(left_pts.len(), seg_count + 1);
    debug_assert_eq!(right_pts.len(), seg_count + 1);
    let mut poly: Polygon = Vec::new();
    poly.extend_from_slice(left_pts);
    right_pts.reverse();
    poly.extend_from_slice(&right_pts);

    let mut labels: Vec<StrokeOrigin> = Vec::with_capacity(poly.len());
    for k in 0..seg_count {
        labels.push(StrokeOrigin::Side {
            seg: k,
            side: OffsetSide::Left,
        });
    }
    labels.push(StrokeOrigin::Cap { end: CapEnd::End });
    for k in (0..seg_count).rev() {
        labels.push(StrokeOrigin::Side {
            seg: k,
            side: OffsetSide::Right,
        });
    }
    labels.push(StrokeOrigin::Cap { end: CapEnd::Start });
    debug_assert_eq!(labels.len(), poly.len());

    (poly, labels)
}

fn build_open_offsets(
    verts: &[(f64, f64)],
    dirs: &[(f64, f64)],
    left_w: f64,
    right_w: f64,
    left_pts: &mut Vec<(f64, f64)>,
    right_pts: &mut Vec<(f64, f64)>,
) {
    let n = verts.len();
    let seg_count = n - 1;

    // First vertex: flat end cap perpendicular to outgoing direction.
    let n0 = left_normal(dirs[0]);
    left_pts.push(offset_point(verts[0], n0, left_w));
    right_pts.push(offset_point(verts[0], n0, -right_w));

    // Interior vertices: compute join chain.
    for i in 1..seg_count {
        let join = compute_join(verts[i], dirs[i - 1], dirs[i], left_w, right_w);
        left_pts.extend(join.left);
        right_pts.extend(join.right);
    }

    // Last vertex: flat end cap perpendicular to incoming direction.
    let nl = left_normal(dirs[seg_count - 1]);
    left_pts.push(offset_point(verts[n - 1], nl, left_w));
    right_pts.push(offset_point(verts[n - 1], nl, -right_w));
}

fn build_closed_offsets(
    verts: &[(f64, f64)],
    dirs: &[(f64, f64)],
    left_w: f64,
    right_w: f64,
    left_pts: &mut Vec<(f64, f64)>,
    right_pts: &mut Vec<(f64, f64)>,
) {
    let seg_count = dirs.len();

    for i in 0..seg_count {
        let prev = if i == 0 { seg_count - 1 } else { i - 1 };
        let join = compute_join(verts[i], dirs[prev], dirs[i], left_w, right_w);
        left_pts.extend(join.left);
        right_pts.extend(join.right);
    }
}

/// Local chain of offset vertices emitted for a single join.
///
/// Each side emits exactly one point: the miter intersection of the two
/// offset edges. Sharp acute corners are preserved (no bevel chamfer). When
/// the inner miter flips at very acute concave corners, polygon union
/// downstream removes the inverted region.
struct JoinResult {
    left: Vec<(f64, f64)>,
    right: Vec<(f64, f64)>,
}

/// Computes offset vertices at a single interior join.
///
/// For a CCW polyline, `left_normal` points inward, so:
///   - left side = inner boundary (shrinks at convex corners, spikes at concave)
///   - right side = outer boundary (spikes at convex corners, shrinks at concave)
///
/// `cross = dir_in × dir_out`:
///   - `> 0`: left turn → convex for CCW → outer side is `right`, inner is `left`
///   - `< 0`: right turn → concave for CCW → outer side is `left`, inner is `right`
///   - `≈ 0`: collinear → single offset point, no miter needed
///
/// Both sides always take the miter intersection — sharp corners stay sharp.
/// At very acute angles where the miter would otherwise be cut off (bevel),
/// the spike is kept; the polygon-union arrangement downstream cleans up any
/// self-intersecting region the spike introduces.
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
        left: miter_chain(left_miter, lp_out),
        right: miter_chain(right_miter, rp_out),
    }
}

/// Always take the miter intersection. Falls back to the outgoing edge's
/// offset endpoint only when the two offset directions are numerically
/// parallel (`line_intersect` returned `None`); the single-point fallback
/// preserves continuity.
fn miter_chain(miter: Option<(f64, f64)>, fallback: (f64, f64)) -> Vec<(f64, f64)> {
    match miter {
        Some(m) => vec![m],
        None => vec![fallback],
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
        // Sharp acute turn — the outer side keeps the (potentially long)
        // miter point instead of beveling. Result must be non-degenerate.
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
    fn join_result_acute_keeps_sharp_miter_both_sides() {
        // Even at a very acute U-turn, both sides keep a single miter point
        // — no bevel fallback. Sharp acute corners are preserved.
        let dir_out = {
            let (dx, dy): (f64, f64) = (-1.0, 0.05);
            let l = (dx * dx + dy * dy).sqrt();
            (dx / l, dy / l)
        };
        let join = compute_join((0.0, 0.0), (1.0, 0.0), dir_out, 0.3, 0.3);
        assert_eq!(join.left.len(), 1, "left side keeps single miter point");
        assert_eq!(join.right.len(), 1, "right side keeps single miter point");
    }

    #[test]
    fn join_result_collinear_single_point() {
        let join = compute_join((1.0, 0.0), (1.0, 0.0), (1.0, 0.0), 0.3, 0.3);
        assert_eq!(join.left.len(), 1);
        assert_eq!(join.right.len(), 1);
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
