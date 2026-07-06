//! Per-segment provenance for wall footprints.
//!
//! [`super::WallOutline2D::execute_faces_with_provenance`] reports, for
//! every boundary segment of every output [`super::WallFootprint2D`],
//! **where that segment came from in the input centerlines** — which
//! polyline, which centerline edge, which side (or which end cap), and
//! which surviving fragment of that source. Consumers derive stable
//! face names (e.g. `SegmentTag`s for
//! `MakeSegmentedPrism::with_segment_tags`) from this structured data;
//! geolis itself stays identity-dumb and never accepts caller strings
//! here.
//!
//! # How provenance is obtained
//!
//! Labels are **threaded** through the 2D boolean pipeline: every raw
//! stroke-polygon edge enters the arrangement engine carrying its
//! source site, and the engine's split / snap / dedup / classify /
//! face-walk stages preserve the label per surviving sub-edge. No
//! geometric matching is ever performed, so provenance is exact even
//! under the engine's `WALL_EPS` vertex snapping.
//!
//! # Determinism and stability
//!
//! - Same input ⇒ same provenance, including fragment ordinals.
//! - Fragments of one source segment are ordered by their position
//!   **along that source segment** (ascending parameter; for tessellated
//!   arc edges, by tessellation chord first). An edit that does not
//!   change which pieces of a source segment survive therefore does not
//!   renumber that segment's fragments.
//! - Fragment ordinals of a source segment MAY shift when that same
//!   segment gains or loses surviving pieces (e.g. a junction is added
//!   or removed on it). This is inherent: ordinals are dense.
//! - When two strokes contribute geometrically identical (coincident
//!   collinear) boundary pieces, the surviving piece is attributed to
//!   the lexicographically smallest source — earliest input polyline,
//!   then ring, then edge — deterministically.

use crate::operations::boolean_2d::{RingRef, SegmentSite, TracedFace};

/// Which side of the centerline an offset segment lies on, relative to
/// the centerline's own traversal direction (`Left` is +90° from the
/// segment direction). For closed centerlines the side always refers to
/// the direction the caller supplied the vertices in, regardless of
/// their winding.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum OffsetSide {
    Left,
    Right,
}

/// Which end of an open centerline a flat end cap closes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum CapEnd {
    Start,
    End,
}

/// Where a footprint boundary segment originated in its source
/// centerline polyline.
///
/// There is no join variant: wall joins are miters, so a join
/// contributes a single shared vertex to the stroke polygon, never a
/// segment of its own. Every boundary segment is either a side offset
/// or an end cap.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum SegmentOrigin {
    /// Offset of centerline edge `edge` (connecting centerline vertex
    /// `edge` to vertex `edge + 1`, wrapping for closed polylines) on
    /// the given side. An arc (bulge) centerline edge keeps a single
    /// `edge` index for all of its tessellated chords.
    Side { edge: usize, side: OffsetSide },
    /// Flat end cap of an open centerline.
    Cap { end: CapEnd },
}

/// Provenance of one output boundary segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SegmentProvenance {
    /// Index of the source polyline in the `Vec<Pline>` passed to
    /// [`super::WallOutline2D`] (original input position, including any
    /// entries that were skipped as too short).
    pub pline: usize,
    /// Structural origin within that polyline.
    pub origin: SegmentOrigin,
    /// Deterministic ordinal of the surviving piece when union trimming
    /// split this source into several pieces, ordered along the source
    /// segment (see the module docs for the exact rule). `0` when the
    /// source survived in one piece.
    pub fragment: u32,
}

/// Per-ring provenance aligned 1:1 with a [`super::WallFootprint2D`]:
/// `outer()[k]` describes the outer-ring segment from vertex `k` to
/// vertex `(k + 1) % n`, and `holes()[h][k]` likewise for hole `h`.
/// Hole rings are union outputs of the same labelled arrangement, so
/// they carry full provenance exactly like the outer ring.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FootprintProvenance {
    pub(super) outer: Vec<SegmentProvenance>,
    pub(super) holes: Vec<Vec<SegmentProvenance>>,
}

impl FootprintProvenance {
    /// Per-segment provenance of the outer ring, aligned with
    /// `WallFootprint2D::outer()`'s segments.
    #[must_use]
    pub fn outer(&self) -> &[SegmentProvenance] {
        &self.outer
    }

    /// Per-segment provenance of each hole ring, aligned with
    /// `WallFootprint2D::holes()`.
    #[must_use]
    pub fn holes(&self) -> &[Vec<SegmentProvenance>] {
        &self.holes
    }
}

// === Crate-internal assembly ===

/// Source description of one stroke-polygon edge, built by
/// `WallOutline2D` before the union and resolved from the engine's
/// [`SegmentSite`]s afterwards.
pub(super) struct EdgeSource {
    pub pline: usize,
    pub origin: SegmentOrigin,
    /// Ordinal of the tessellated stroke segment this edge offsets
    /// (identity for line edges; chord index for arc edges). Used only
    /// to order fragments of one source along the source.
    pub tess_ord: usize,
    /// Supporting stroke-polygon edge geometry, for ordering fragments
    /// by their parameter along the source.
    pub a: (f64, f64),
    pub b: (f64, f64),
}

/// Per-input lookup table: `SegmentSite {input, ring, edge}` resolves to
/// `tables[input].get(ring, edge)`.
pub(super) struct InputEdgeSources {
    pub outer: Vec<EdgeSource>,
    pub holes: Vec<Vec<EdgeSource>>,
}

impl InputEdgeSources {
    fn get(&self, site: SegmentSite) -> &EdgeSource {
        match site.ring {
            RingRef::Outer => &self.outer[site.edge],
            RingRef::Hole(h) => &self.holes[h][site.edge],
        }
    }
}

/// Borrowed view of one output ring: `(ring points, per-edge sites)`.
type RingView<'a> = (&'a [(f64, f64)], &'a [SegmentSite]);

/// One maximal cyclic run of consecutive ring edges sharing the same
/// `(pline, origin)` key — i.e. one surviving fragment of one source.
struct Run {
    key: (usize, SegmentOrigin),
    /// Lexicographic-min `(tess_ord, parameter-along-source)` over the
    /// run's edges; orders fragments along their source segment.
    order: (usize, f64),
    face: usize,
    /// 0 = outer ring, `1 + h` = hole `h` (tie-break only).
    ring: usize,
    /// Edge indices within the ring, in ring order.
    edges: Vec<usize>,
}

/// Compute aligned [`FootprintProvenance`] for every traced face.
///
/// Every output edge is prefilled with its own `(pline, origin)` and
/// fragment `0`; a second pass groups edges into runs per source and
/// overwrites fragments with the deterministic run ordinal.
pub(super) fn footprint_provenances(
    faces: &[TracedFace],
    sources: &[InputEdgeSources],
) -> Vec<FootprintProvenance> {
    let source_of = |site: SegmentSite| -> &EdgeSource { sources[site.input].get(site) };

    // Prefill with fragment 0.
    let prefill = |sites: &[SegmentSite]| -> Vec<SegmentProvenance> {
        sites
            .iter()
            .map(|&s| {
                let src = source_of(s);
                SegmentProvenance {
                    pline: src.pline,
                    origin: src.origin,
                    fragment: 0,
                }
            })
            .collect()
    };
    let mut out: Vec<FootprintProvenance> = faces
        .iter()
        .map(|tf| FootprintProvenance {
            outer: prefill(&tf.outer_sites),
            holes: tf.hole_sites.iter().map(|s| prefill(s)).collect(),
        })
        .collect();

    // Collect runs over every ring of every face.
    let mut runs: Vec<Run> = Vec::new();
    for (fi, tf) in faces.iter().enumerate() {
        let mut rings: Vec<RingView<'_>> =
            vec![(tf.face.outer.as_slice(), tf.outer_sites.as_slice())];
        for (h, sites) in tf.hole_sites.iter().enumerate() {
            rings.push((tf.face.holes[h].as_slice(), sites.as_slice()));
        }
        for (ri, (pts, sites)) in rings.into_iter().enumerate() {
            collect_ring_runs(fi, ri, pts, sites, &source_of, &mut runs);
        }
    }

    // Deterministic ordering: group by key, order fragments along the
    // source, tie-break structurally.
    runs.sort_by(|x, y| {
        x.key
            .cmp(&y.key)
            .then_with(|| x.order.0.cmp(&y.order.0))
            .then_with(|| {
                x.order
                    .1
                    .partial_cmp(&y.order.1)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| (x.face, x.ring, x.edges[0]).cmp(&(y.face, y.ring, y.edges[0])))
    });

    let mut prev_key: Option<(usize, SegmentOrigin)> = None;
    let mut ordinal: u32 = 0;
    for run in &runs {
        if prev_key != Some(run.key) {
            prev_key = Some(run.key);
            ordinal = 0;
        }
        let prov = &mut out[run.face];
        for &e in &run.edges {
            let slot = if run.ring == 0 {
                &mut prov.outer[e]
            } else {
                &mut prov.holes[run.ring - 1][e]
            };
            slot.fragment = ordinal;
        }
        ordinal += 1;
    }

    out
}

/// Split one closed ring into maximal cyclic runs of equal
/// `(pline, origin)` keys and append them to `runs`.
fn collect_ring_runs<'a>(
    face: usize,
    ring: usize,
    pts: &[(f64, f64)],
    sites: &[SegmentSite],
    source_of: &impl Fn(SegmentSite) -> &'a EdgeSource,
    runs: &mut Vec<Run>,
) {
    let m = sites.len();
    if m == 0 {
        return;
    }
    debug_assert_eq!(pts.len(), m);

    let key_of = |e: usize| -> (usize, SegmentOrigin) {
        let src = source_of(sites[e]);
        (src.pline, src.origin)
    };
    // Parameter of point `p` along the source of edge `e` (unnormalised
    // — monotonic along the supporting line, which is all ordering
    // needs).
    let param = |e: usize, p: (f64, f64)| -> f64 {
        let src = source_of(sites[e]);
        (p.0 - src.a.0) * (src.b.0 - src.a.0) + (p.1 - src.a.1) * (src.b.1 - src.a.1)
    };
    let edge_order = |e: usize| -> (usize, f64) {
        let t0 = param(e, pts[e]);
        let t1 = param(e, pts[(e + 1) % m]);
        (source_of(sites[e]).tess_ord, t0.min(t1))
    };
    let min_order = |a: (usize, f64), b: (usize, f64)| -> (usize, f64) {
        if (b.0, b.1) < (a.0, a.1) {
            b
        } else {
            a
        }
    };

    // Rotate the scan start to a key boundary so no run is split by the
    // ring's arbitrary index origin.
    let start = (0..m).find(|&e| key_of(e) != key_of((e + m - 1) % m));
    let Some(start) = start else {
        // Whole ring is one source: single run.
        let mut order = edge_order(0);
        for e in 1..m {
            order = min_order(order, edge_order(e));
        }
        runs.push(Run {
            key: key_of(0),
            order,
            face,
            ring,
            edges: (0..m).collect(),
        });
        return;
    };

    let mut current_edges: Vec<usize> = Vec::new();
    let mut current_key = key_of(start);
    let mut current_order = (usize::MAX, f64::INFINITY);
    for k in 0..m {
        let e = (start + k) % m;
        let key = key_of(e);
        if key != current_key && !current_edges.is_empty() {
            runs.push(Run {
                key: current_key,
                order: current_order,
                face,
                ring,
                edges: std::mem::take(&mut current_edges),
            });
            current_order = (usize::MAX, f64::INFINITY);
        }
        current_key = key;
        current_order = min_order(current_order, edge_order(e));
        current_edges.push(e);
    }
    runs.push(Run {
        key: current_key,
        order: current_order,
        face,
        ring,
        edges: current_edges,
    });
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::super::{WallFootprint2D, WallOutline2D};
    use super::*;
    use crate::geometry::pline::{Pline, PlineVertex};
    use crate::math::Point3;

    fn open_pline(pts: &[(f64, f64)]) -> Pline {
        Pline::from_points(
            &pts.iter()
                .map(|&(x, y)| Point3::new(x, y, 0.0))
                .collect::<Vec<_>>(),
            false,
        )
    }

    fn closed_pline(pts: &[(f64, f64)]) -> Pline {
        Pline::from_points(
            &pts.iter()
                .map(|&(x, y)| Point3::new(x, y, 0.0))
                .collect::<Vec<_>>(),
            true,
        )
    }

    fn run(plines: Vec<Pline>, hw: f64) -> Vec<(WallFootprint2D, FootprintProvenance)> {
        WallOutline2D::new(plines, hw)
            .execute_faces_with_provenance()
            .expect("execute_faces_with_provenance must succeed")
    }

    /// One output ring as `(ring points, ring provenance)`.
    type OwnedRing = (Vec<(f64, f64)>, Vec<SegmentProvenance>);

    /// All rings of all faces as `(ring points, ring provenance)` pairs.
    fn all_rings(result: &[(WallFootprint2D, FootprintProvenance)]) -> Vec<OwnedRing> {
        let mut out = Vec::new();
        for (f, p) in result {
            let ring_pts = |pl: &Pline| -> Vec<(f64, f64)> {
                pl.vertices.iter().map(|v| (v.x, v.y)).collect()
            };
            out.push((ring_pts(f.outer()), p.outer().to_vec()));
            for (h, hp) in f.holes().iter().zip(p.holes()) {
                out.push((ring_pts(h), hp.clone()));
            }
        }
        out
    }

    /// Asserts every ring's provenance array aligns 1:1 with its segment
    /// count, and that every line-centerline edge lies on the supporting
    /// line its provenance names (within snap tolerance).
    fn assert_aligned_and_on_source_lines(
        result: &[(WallFootprint2D, FootprintProvenance)],
        plines: &[Pline],
        left_w: f64,
        right_w: f64,
    ) {
        for (pts, prov) in all_rings(result) {
            assert_eq!(pts.len(), prov.len(), "provenance must align with ring");
            for (e, sp) in prov.iter().enumerate() {
                let pl = &plines[sp.pline];
                let vert_count = pl.vertices.len();
                let (base, dir) = match sp.origin {
                    SegmentOrigin::Side { edge, side } => {
                        let seg_a = &pl.vertices[edge];
                        let seg_b = &pl.vertices[(edge + 1) % vert_count];
                        if seg_a.bulge.abs() > 1e-12 {
                            continue; // arc edges checked separately
                        }
                        let (dx, dy) = (seg_b.x - seg_a.x, seg_b.y - seg_a.y);
                        let len = (dx * dx + dy * dy).sqrt();
                        let seg_dir = (dx / len, dy / len);
                        let nn = (-seg_dir.1, seg_dir.0);
                        let width = match side {
                            OffsetSide::Left => left_w,
                            OffsetSide::Right => -right_w,
                        };
                        ((seg_a.x + width * nn.0, seg_a.y + width * nn.1), seg_dir)
                    }
                    SegmentOrigin::Cap { end } => {
                        let (cap_v, prev, next) = match end {
                            CapEnd::Start => (&pl.vertices[0], &pl.vertices[0], &pl.vertices[1]),
                            CapEnd::End => (
                                &pl.vertices[vert_count - 1],
                                &pl.vertices[vert_count - 2],
                                &pl.vertices[vert_count - 1],
                            ),
                        };
                        let (dx, dy) = (next.x - prev.x, next.y - prev.y);
                        let len = (dx * dx + dy * dy).sqrt();
                        // Cap line runs along the segment normal through cap_v.
                        ((cap_v.x, cap_v.y), (-dy / len, dx / len))
                    }
                };
                for p in [pts[e], pts[(e + 1) % pts.len()]] {
                    let perp = (p.0 - base.0) * dir.1 - (p.1 - base.1) * dir.0;
                    assert!(
                        perp.abs() < 1e-5,
                        "ring edge {e} with provenance {sp:?} not on its \
                         source supporting line: perp={perp}"
                    );
                }
            }
        }
    }

    // ===== W2a acceptance tests =====

    #[test]
    fn straight_two_edge_open_centerline_full_provenance() {
        let pline = open_pline(&[(0.0, 0.0), (4.0, 0.0), (4.0, 3.0)]);
        let result = run(vec![pline.clone()], 0.3);
        assert_eq!(result.len(), 1);
        assert!(result[0].0.holes().is_empty());
        assert_aligned_and_on_source_lines(&result, &[pline], 0.3, 0.3);

        let prov = result[0].1.outer();
        // No trimming: every fragment is 0, every segment is pline 0.
        for sp in prov {
            assert_eq!(sp.pline, 0);
            assert_eq!(sp.fragment, 0);
        }
        // All six origins of a 2-edge open stroke must appear.
        let mut origins: Vec<SegmentOrigin> = prov.iter().map(|sp| sp.origin).collect();
        origins.sort();
        origins.dedup();
        let mut expected = vec![
            SegmentOrigin::Side {
                edge: 0,
                side: OffsetSide::Left,
            },
            SegmentOrigin::Side {
                edge: 0,
                side: OffsetSide::Right,
            },
            SegmentOrigin::Side {
                edge: 1,
                side: OffsetSide::Left,
            },
            SegmentOrigin::Side {
                edge: 1,
                side: OffsetSide::Right,
            },
            SegmentOrigin::Cap { end: CapEnd::Start },
            SegmentOrigin::Cap { end: CapEnd::End },
        ];
        expected.sort();
        assert_eq!(origins, expected);
    }

    #[test]
    fn t_junction_trims_bar_left_side_into_two_ordered_fragments() {
        let bar = open_pline(&[(0.0, 0.0), (4.0, 0.0)]);
        let stem = open_pline(&[(2.0, 0.0), (2.0, 3.0)]);
        let result = run(vec![bar.clone(), stem.clone()], 0.15);
        assert_eq!(result.len(), 1, "T junction must union into one face");
        assert_aligned_and_on_source_lines(&result, &[bar, stem], 0.15, 0.15);

        let rings = all_rings(&result);
        let bar_left = SegmentOrigin::Side {
            edge: 0,
            side: OffsetSide::Left,
        };
        let bar_right = SegmentOrigin::Side {
            edge: 0,
            side: OffsetSide::Right,
        };

        let mut bar_left_fragments: Vec<(u32, f64)> = Vec::new(); // (fragment, mid x)
        let mut bar_right_fragments: Vec<u32> = Vec::new();
        let mut stem_origins: Vec<SegmentOrigin> = Vec::new();
        for (pts, prov) in &rings {
            for (e, sp) in prov.iter().enumerate() {
                let mid_x = (pts[e].0 + pts[(e + 1) % pts.len()].0) * 0.5;
                if sp.pline == 0 && sp.origin == bar_left {
                    bar_left_fragments.push((sp.fragment, mid_x));
                }
                if sp.pline == 0 && sp.origin == bar_right {
                    bar_right_fragments.push(sp.fragment);
                }
                if sp.pline == 1 {
                    stem_origins.push(sp.origin);
                }
            }
        }

        // Bar left side (y = +0.15) is trimmed by the stem into exactly
        // two fragments, numbered along the bar's direction (+x).
        let frags: Vec<u32> = {
            let mut f: Vec<u32> = bar_left_fragments.iter().map(|(f, _)| *f).collect();
            f.sort_unstable();
            f.dedup();
            f
        };
        assert_eq!(frags, vec![0, 1], "bar left side must split in two");
        for &(frag, mid_x) in &bar_left_fragments {
            if frag == 0 {
                assert!(
                    mid_x < 1.9,
                    "fragment 0 must be the -x piece; mid_x={mid_x}"
                );
            } else {
                assert!(
                    mid_x > 2.1,
                    "fragment 1 must be the +x piece; mid_x={mid_x}"
                );
            }
        }
        // Bar right side (y = -0.15) is untrimmed: single fragment 0.
        bar_right_fragments.sort_unstable();
        bar_right_fragments.dedup();
        assert_eq!(bar_right_fragments, vec![0]);

        // The stem's start cap is swallowed by the bar material; its end
        // cap and both sides survive with fragment 0.
        stem_origins.sort();
        stem_origins.dedup();
        assert_eq!(
            stem_origins,
            vec![
                SegmentOrigin::Side {
                    edge: 0,
                    side: OffsetSide::Left,
                },
                SegmentOrigin::Side {
                    edge: 0,
                    side: OffsetSide::Right,
                },
                SegmentOrigin::Cap { end: CapEnd::End },
            ]
        );
    }

    #[test]
    fn provenance_is_deterministic_across_runs() {
        let make = || {
            vec![
                open_pline(&[(0.0, 0.0), (4.0, 0.0)]),
                open_pline(&[(2.0, 0.0), (2.0, 3.0)]),
                closed_pline(&[(6.0, 0.0), (9.0, 0.0), (9.0, 3.0), (6.0, 3.0)]),
            ]
        };
        let a = run(make(), 0.15);
        let b = run(make(), 0.15);
        assert_eq!(a.len(), b.len());
        for ((fa, pa), (fb, pb)) in a.iter().zip(&b) {
            assert_eq!(fa.outer().vertices, fb.outer().vertices);
            assert_eq!(pa, pb, "provenance must be bit-identical across runs");
        }
    }

    /// Editing one wall must not change the provenance (origins AND
    /// fragment ordinals) reported for segments of unrelated walls whose
    /// surviving pieces did not change.
    #[test]
    fn unrelated_edit_preserves_other_walls_provenance() {
        let junction = |far_len: f64| {
            vec![
                open_pline(&[(0.0, 0.0), (4.0, 0.0)]),
                open_pline(&[(2.0, 0.0), (2.0, 3.0)]),
                // Unrelated far wall, edited between the two runs.
                open_pline(&[(0.0, 6.0), (far_len, 6.0)]),
            ]
        };
        let a = run(junction(4.0), 0.15);
        let b = run(junction(5.0), 0.15);

        // Collect (pline, origin, fragment, quantised edge midpoint) for
        // the untouched walls 0 and 1.
        let collect = |result: &[(WallFootprint2D, FootprintProvenance)]| {
            let mut items: Vec<(usize, SegmentOrigin, u32, (i64, i64))> = Vec::new();
            for (pts, prov) in all_rings(result) {
                for (e, sp) in prov.iter().enumerate() {
                    if sp.pline > 1 {
                        continue;
                    }
                    let m = (
                        (pts[e].0 + pts[(e + 1) % pts.len()].0) * 0.5,
                        (pts[e].1 + pts[(e + 1) % pts.len()].1) * 0.5,
                    );
                    #[allow(clippy::cast_possible_truncation)]
                    let q = ((m.0 * 1e6).round() as i64, (m.1 * 1e6).round() as i64);
                    items.push((sp.pline, sp.origin, sp.fragment, q));
                }
            }
            items.sort();
            items
        };
        assert_eq!(
            collect(&a),
            collect(&b),
            "editing wall 2 must not renumber walls 0/1 provenance"
        );
    }

    #[test]
    fn arc_bulge_centerline_provenance_on_curved_segments() {
        // Semicircle-ish arc from (0,0) to (4,0), bulge 0.5.
        let pline = Pline {
            vertices: vec![PlineVertex::new(0.0, 0.0, 0.5), PlineVertex::line(4.0, 0.0)],
            closed: false,
        };
        let result = run(vec![pline], 0.3);
        assert_eq!(result.len(), 1);
        let prov = result[0].1.outer();
        assert_eq!(prov.len(), result[0].0.outer().vertices.len());

        let mut left = 0usize;
        let mut right = 0usize;
        let mut caps = 0usize;
        for sp in prov {
            assert_eq!(sp.pline, 0);
            assert_eq!(sp.fragment, 0, "single arc stroke must not fragment");
            match sp.origin {
                SegmentOrigin::Side { edge, side } => {
                    assert_eq!(edge, 0, "all chords must map to the arc edge");
                    match side {
                        OffsetSide::Left => left += 1,
                        OffsetSide::Right => right += 1,
                    }
                }
                SegmentOrigin::Cap { .. } => caps += 1,
            }
        }
        assert!(left >= 2, "curved left offset must span several chords");
        assert!(right >= 2, "curved right offset must span several chords");
        assert_eq!(caps, 2, "both flat end caps must survive");
    }

    #[test]
    fn closed_square_sides_follow_caller_direction() {
        // CCW input: outer boundary is the RIGHT side of the traversal.
        let ccw = closed_pline(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]);
        let result = run(vec![ccw.clone()], 0.3);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0.holes().len(), 1);
        assert_aligned_and_on_source_lines(&result, &[ccw], 0.3, 0.3);
        let (_, prov) = &result[0];
        for sp in prov.outer() {
            assert!(
                matches!(
                    sp.origin,
                    SegmentOrigin::Side {
                        side: OffsetSide::Right,
                        ..
                    }
                ),
                "CCW closed input: outer must be Right sides; got {sp:?}"
            );
        }
        for sp in &prov.holes()[0] {
            assert!(
                matches!(
                    sp.origin,
                    SegmentOrigin::Side {
                        side: OffsetSide::Left,
                        ..
                    }
                ),
                "CCW closed input: hole must be Left sides; got {sp:?}"
            );
        }

        // CW input (same square reversed): sides flip because they refer
        // to the caller's traversal direction.
        let cw = closed_pline(&[(0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0)]);
        let result = run(vec![cw.clone()], 0.3);
        assert_aligned_and_on_source_lines(&result, &[cw], 0.3, 0.3);
        let (_, prov) = &result[0];
        for sp in prov.outer() {
            assert!(
                matches!(
                    sp.origin,
                    SegmentOrigin::Side {
                        side: OffsetSide::Left,
                        ..
                    }
                ),
                "CW closed input: outer must be Left sides; got {sp:?}"
            );
        }
    }

    #[test]
    fn execute_faces_matches_provenance_variant_geometry() {
        let plines = vec![
            open_pline(&[(0.0, 0.0), (4.0, 0.0)]),
            open_pline(&[(2.0, 0.0), (2.0, 3.0)]),
        ];
        let plain = WallOutline2D::new(plines.clone(), 0.15)
            .execute_faces()
            .unwrap();
        let traced = run(plines, 0.15);
        assert_eq!(plain.len(), traced.len());
        for (a, (b, _)) in plain.iter().zip(&traced) {
            assert_eq!(a.outer().vertices, b.outer().vertices);
            assert_eq!(a.holes().len(), b.holes().len());
        }
    }
}
