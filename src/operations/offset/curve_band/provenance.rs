//! Per-segment provenance for wall footprints.
//!
//! [`super::CurveBand2D::execute_faces_with_provenance`] reports, for
//! every boundary segment of every output [`super::BandFootprint2D`],
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

use crate::error::{OperationError, Result};
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
/// There is no join variant. A join that miters contributes a single
/// shared vertex to the stroke polygon and no segment of its own; a join
/// sharp enough to run past the miter limit is truncated by a chamfer,
/// and that chamfer is reported as `Side` on the centerline edge LEAVING
/// the join (see `build_edge_sources`). So every boundary segment is a
/// side offset or an end cap here, at the cost of a chamfer reading as
/// its segment's side — widening this enum is a breaking change for
/// every consumer that matches it, and is warranted only once one of
/// them has to tell the two apart.
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
    /// [`super::CurveBand2D`] (original input position, including any
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

/// Per-ring provenance aligned 1:1 with a [`super::BandFootprint2D`]:
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
    /// `BandFootprint2D::outer()`'s segments.
    #[must_use]
    pub fn outer(&self) -> &[SegmentProvenance] {
        &self.outer
    }

    /// Per-segment provenance of each hole ring, aligned with
    /// `BandFootprint2D::holes()`.
    #[must_use]
    pub fn holes(&self) -> &[Vec<SegmentProvenance>] {
        &self.holes
    }
}

// === Crate-internal assembly ===

/// Source description of one stroke-polygon edge, built by
/// `CurveBand2D` before the union and resolved from the engine's
/// [`SegmentSite`]s afterwards.
pub(super) struct EdgeSource {
    pub pline: usize,
    pub origin: SegmentOrigin,
    /// Ordinal of the tessellated stroke segment this edge belongs to
    /// (identity for line edges; chord index for arc edges). A `Side`
    /// edge offsets that segment; a bevelled `Join`'s chamfer cuts
    /// ACROSS the corner instead (see `stroke::StrokeOrigin::Join`) and
    /// takes the ordinal of the segment LEAVING the join, which keeps
    /// the two adjacent in the ordering. Used only to order fragments
    /// of one source along the source.
    pub tess_ord: usize,
    /// Supporting stroke-polygon edge geometry, for ordering fragments
    /// by their parameter along the source.
    pub a: (f64, f64),
    pub b: (f64, f64),
}

/// Per-input lookup table: `SegmentSite {input, ring, edge}` resolves to
/// `tables[input].get(site)`.
///
/// A stroke-expanded band is always a single ring (see
/// [`super::stroke::stroke_expand_labeled`]), so only `RingRef::Outer`
/// sites resolve. `None` marks an edge that carries no centerline
/// provenance — the closed band's internal slit — or a site the stroke
/// never emitted; both mean the same thing to the caller: an edge that
/// must not have survived the union.
pub(super) struct InputEdgeSources {
    pub edges: Vec<Option<EdgeSource>>,
}

impl InputEdgeSources {
    fn get(&self, site: SegmentSite) -> Option<&EdgeSource> {
        match site.ring {
            RingRef::Outer => self.edges.get(site.edge)?.as_ref(),
            RingRef::Hole(_) => None,
        }
    }
}

/// Source identity of one band-union output edge: which input polyline,
/// and which structural part of it.
type SourceKey = (usize, SegmentOrigin);

/// Where one output ring edge sits along the source it came from.
///
/// Compared lexicographically, and only ever against positions of the
/// *same* source key:
///
/// 1. **coarse ordinal** — which sub-source the edge lies on
///    (tessellation chord ordinal for the band union; band fragment
///    ordinal for the carve).
/// 2. **fine ordinal** — position of that sub-source within the source's
///    own run, for sources whose sub-sources are not collinear and whose
///    parameters therefore cannot be compared (`0` when the coarse
///    ordinal already resolves the sub-source).
/// 3. **parameter** — unnormalised (but monotonic) projection of the
///    edge's earliest endpoint onto its sub-source, which orders pieces
///    that were split out of one sub-source.
pub(super) type SourcePosition = (usize, usize, f64);

/// Per-ring `(source key, position)` pairs of one output face, aligned
/// with the face's rings exactly like [`FootprintProvenance`]: one entry
/// per ring edge, `outer` first, then each hole.
pub(super) struct FaceEdgeKeys<K> {
    pub outer: Vec<(K, SourcePosition)>,
    pub holes: Vec<Vec<(K, SourcePosition)>>,
}

/// Fragment ordinals aligned 1:1 with the rings of a [`FaceEdgeKeys`].
pub(super) struct FaceFragments {
    pub outer: Vec<u32>,
    pub holes: Vec<Vec<u32>>,
}

/// One maximal cyclic run of consecutive ring edges sharing the same
/// source key — i.e. one surviving fragment of one source.
struct Run<K> {
    key: K,
    /// Lexicographic-min [`SourcePosition`] over the run's edges; orders
    /// the fragments of one source along that source.
    position: SourcePosition,
    face: usize,
    /// 0 = outer ring, `1 + h` = hole `h` (tie-break only).
    ring: usize,
    /// Edge indices within the ring, in ring order.
    edges: Vec<usize>,
}

/// Lexicographic [`SourcePosition`] comparison. The parameter component
/// is a plain `f64`; the arrangement never produces NaN coordinates, and
/// a hypothetical NaN degrades to "equal" rather than poisoning the sort.
fn cmp_position(a: SourcePosition, b: SourcePosition) -> std::cmp::Ordering {
    a.0.cmp(&b.0)
        .then_with(|| a.1.cmp(&b.1))
        .then_with(|| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
}

/// Number the surviving fragments of every source deterministically.
///
/// Consecutive ring edges sharing a key form one run — one surviving
/// fragment. Runs are grouped by key, ordered along their source by
/// [`SourcePosition`], and numbered densely from `0`; the structural
/// `(face, ring, first edge)` triple breaks any remaining tie so the
/// result never depends on iteration order.
///
/// Shared by the band union ([`footprint_provenances`], keyed by
/// `(pline, origin)`) and the band carve ([`super::carve`], keyed by
/// base source or cutter index), so both obey one numbering contract.
pub(super) fn number_fragments<K: Copy + Ord>(faces: &[FaceEdgeKeys<K>]) -> Vec<FaceFragments> {
    let mut out: Vec<FaceFragments> = faces
        .iter()
        .map(|f| FaceFragments {
            outer: vec![0; f.outer.len()],
            holes: f.holes.iter().map(|h| vec![0; h.len()]).collect(),
        })
        .collect();

    // Collect runs over every ring of every face.
    let mut runs: Vec<Run<K>> = Vec::new();
    for (fi, f) in faces.iter().enumerate() {
        collect_ring_runs(fi, 0, &f.outer, &mut runs);
        for (h, ring) in f.holes.iter().enumerate() {
            collect_ring_runs(fi, 1 + h, ring, &mut runs);
        }
    }

    // Deterministic ordering: group by key, order fragments along the
    // source, tie-break structurally.
    runs.sort_by(|x, y| {
        x.key
            .cmp(&y.key)
            .then_with(|| cmp_position(x.position, y.position))
            .then_with(|| (x.face, x.ring, x.edges[0]).cmp(&(y.face, y.ring, y.edges[0])))
    });

    let mut prev_key: Option<K> = None;
    let mut ordinal: u32 = 0;
    for run in &runs {
        if prev_key != Some(run.key) {
            prev_key = Some(run.key);
            ordinal = 0;
        }
        let frag = &mut out[run.face];
        for &e in &run.edges {
            let slot = if run.ring == 0 {
                &mut frag.outer[e]
            } else {
                &mut frag.holes[run.ring - 1][e]
            };
            *slot = ordinal;
        }
        ordinal += 1;
    }

    out
}

/// Split one closed ring into maximal cyclic runs of equal key and
/// append them to `runs`.
fn collect_ring_runs<K: Copy + Ord>(
    face: usize,
    ring: usize,
    edges: &[(K, SourcePosition)],
    runs: &mut Vec<Run<K>>,
) {
    let m = edges.len();
    if m == 0 {
        return;
    }
    let min_position = |a: SourcePosition, b: SourcePosition| -> SourcePosition {
        if cmp_position(b, a) == std::cmp::Ordering::Less {
            b
        } else {
            a
        }
    };

    // Rotate the scan start to a key boundary so no run is split by the
    // ring's arbitrary index origin.
    let start = (0..m).find(|&e| edges[e].0 != edges[(e + m - 1) % m].0);
    let Some(start) = start else {
        // Whole ring is one source: single run.
        let mut position = edges[0].1;
        for edge in edges.iter().skip(1) {
            position = min_position(position, edge.1);
        }
        runs.push(Run {
            key: edges[0].0,
            position,
            face,
            ring,
            edges: (0..m).collect(),
        });
        return;
    };

    let mut current_edges: Vec<usize> = Vec::new();
    let mut current_key = edges[start].0;
    let mut current_position = (usize::MAX, usize::MAX, f64::INFINITY);
    for k in 0..m {
        let e = (start + k) % m;
        let key = edges[e].0;
        if key != current_key && !current_edges.is_empty() {
            runs.push(Run {
                key: current_key,
                position: current_position,
                face,
                ring,
                edges: std::mem::take(&mut current_edges),
            });
            current_position = (usize::MAX, usize::MAX, f64::INFINITY);
        }
        current_key = key;
        current_position = min_position(current_position, edges[e].1);
        current_edges.push(e);
    }
    runs.push(Run {
        key: current_key,
        position: current_position,
        face,
        ring,
        edges: current_edges,
    });
}

/// Compute aligned [`FootprintProvenance`] for every traced face.
///
/// Each output edge resolves its [`SegmentSite`] to the stroke edge it
/// came from, yielding a `(pline, origin)` key and a position along that
/// source; [`number_fragments`] then turns runs of equal key into dense
/// fragment ordinals.
///
/// # Errors
///
/// [`OperationError::Failed`] when an output edge resolves to a stroke
/// edge that carries no centerline provenance — the closed band's
/// internal slit, which has material on both sides and must always be
/// dropped by the arrangement. Reaching an output boundary means the
/// half-edge classification broke, and provenance is reported as broken
/// rather than guessed.
pub(super) fn footprint_provenances(
    faces: &[TracedFace],
    sources: &[InputEdgeSources],
) -> Result<Vec<FootprintProvenance>> {
    let source_of = |site: SegmentSite| -> Result<&EdgeSource> {
        sources[site.input].get(site).ok_or_else(|| {
            OperationError::Failed(format!(
                "curve_band provenance: output edge resolves to stroke edge \
                 {site:?}, which bounds no band material (internal slit) — \
                 the arrangement must drop it"
            ))
            .into()
        })
    };

    let ring_keys = |pts: &[(f64, f64)],
                     sites: &[SegmentSite]|
     -> Result<Vec<(SourceKey, SourcePosition)>> {
        debug_assert_eq!(pts.len(), sites.len());
        let m = sites.len();
        (0..m)
            .map(|e| {
                let src = source_of(sites[e])?;
                // Parameter along the stroke edge (unnormalised —
                // monotonic along the supporting line, which is all
                // ordering needs).
                let param = |p: (f64, f64)| -> f64 {
                    (p.0 - src.a.0) * (src.b.0 - src.a.0) + (p.1 - src.a.1) * (src.b.1 - src.a.1)
                };
                let t0 = param(pts[e]);
                let t1 = param(pts[(e + 1) % m]);
                // The tessellation ordinal already names the sub-source,
                // so the fine ordinal is unused here.
                Ok(((src.pline, src.origin), (src.tess_ord, 0, t0.min(t1))))
            })
            .collect()
    };

    let keys: Vec<FaceEdgeKeys<SourceKey>> = faces
        .iter()
        .map(|tf| -> Result<FaceEdgeKeys<SourceKey>> {
            Ok(FaceEdgeKeys {
                outer: ring_keys(&tf.face.outer, &tf.outer_sites)?,
                holes: tf
                    .face
                    .holes
                    .iter()
                    .zip(&tf.hole_sites)
                    .map(|(h, s)| ring_keys(h, s))
                    .collect::<Result<Vec<_>>>()?,
            })
        })
        .collect::<Result<Vec<_>>>()?;

    let fragments = number_fragments(&keys);

    let ring_provenance =
        |ks: &[(SourceKey, SourcePosition)], fs: &[u32]| -> Vec<SegmentProvenance> {
            ks.iter()
                .zip(fs)
                .map(|(&((pline, origin), _), &fragment)| SegmentProvenance {
                    pline,
                    origin,
                    fragment,
                })
                .collect()
        };
    Ok(keys
        .iter()
        .zip(fragments)
        .map(|(k, frag)| FootprintProvenance {
            outer: ring_provenance(&k.outer, &frag.outer),
            holes: k
                .holes
                .iter()
                .zip(&frag.holes)
                .map(|(ks, fs)| ring_provenance(ks, fs))
                .collect(),
        })
        .collect())
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::super::{BandFootprint2D, CurveBand2D};
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

    fn run(plines: Vec<Pline>, hw: f64) -> Vec<(BandFootprint2D, FootprintProvenance)> {
        CurveBand2D::new(plines, hw)
            .execute_faces_with_provenance()
            .expect("execute_faces_with_provenance must succeed")
    }

    /// One output ring as `(ring points, ring provenance)`.
    type OwnedRing = (Vec<(f64, f64)>, Vec<SegmentProvenance>);

    /// All rings of all faces as `(ring points, ring provenance)` pairs.
    fn all_rings(result: &[(BandFootprint2D, FootprintProvenance)]) -> Vec<OwnedRing> {
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
        result: &[(BandFootprint2D, FootprintProvenance)],
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
        let collect = |result: &[(BandFootprint2D, FootprintProvenance)]| {
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

    // ===== Closed centerlines — including self-crossing ones =====

    /// Material area of the whole result: outer rings minus their holes.
    fn band_area(result: &[(BandFootprint2D, FootprintProvenance)]) -> f64 {
        let ring_area = |pl: &Pline| -> f64 {
            let n = pl.vertices.len();
            (0..n)
                .map(|i| {
                    let (a, b) = (&pl.vertices[i], &pl.vertices[(i + 1) % n]);
                    a.x * b.y - b.x * a.y
                })
                .sum::<f64>()
                * 0.5
        };
        result
            .iter()
            .map(|(f, _)| ring_area(f.outer()) + f.holes().iter().map(ring_area).sum::<f64>())
            .sum()
    }

    /// Whether the band lays material over `p` (inside an outer ring and
    /// not inside one of its holes).
    fn band_covers(result: &[(BandFootprint2D, FootprintProvenance)], p: (f64, f64)) -> bool {
        use crate::operations::boolean_2d::{point_in_polygon_class, PointClass};
        let ring_pts =
            |pl: &Pline| -> Vec<(f64, f64)> { pl.vertices.iter().map(|v| (v.x, v.y)).collect() };
        result.iter().any(|(f, _)| {
            matches!(
                point_in_polygon_class(p, &ring_pts(f.outer())),
                PointClass::Inside | PointClass::Boundary
            ) && !f
                .holes()
                .iter()
                .any(|h| point_in_polygon_class(p, &ring_pts(h)) == PointClass::Inside)
        })
    }

    /// Every segment of a closed centerline must (a) carry material at
    /// its midpoint and (b) appear on the output boundary with BOTH of
    /// its offset sides — the two properties the annulus band model
    /// silently broke on a self-crossing ring.
    fn assert_every_segment_carries_material(
        result: &[(BandFootprint2D, FootprintProvenance)],
        ring: &[(f64, f64)],
    ) {
        let mut origins: Vec<SegmentOrigin> = all_rings(result)
            .iter()
            .flat_map(|(_, prov)| prov.iter().map(|sp| sp.origin))
            .collect();
        origins.sort();
        origins.dedup();
        for seg in 0..ring.len() {
            let a = ring[seg];
            let b = ring[(seg + 1) % ring.len()];
            let mid = ((a.0 + b.0) * 0.5, (a.1 + b.1) * 0.5);
            assert!(
                band_covers(result, mid),
                "segment {seg} (midpoint {mid:?}) carries no material"
            );
            for side in [OffsetSide::Left, OffsetSide::Right] {
                assert!(
                    origins.contains(&SegmentOrigin::Side { edge: seg, side }),
                    "segment {seg}'s {side:?} offset is absent from the \
                     output provenance; present origins: {origins:?}"
                );
            }
        }
        // A closed centerline has no ends, so the closed band's internal
        // slit must never reach an output boundary.
        assert!(
            !origins
                .iter()
                .any(|o| matches!(o, SegmentOrigin::Cap { .. })),
            "closed centerline reported an end cap: {origins:?}"
        );
    }

    /// The user's crossing loop: leg 2 crosses leg 0 at `(2.667, 0)`.
    /// The annulus band model dropped segments 0 and 3 entirely (area
    /// 1.68 of an expected 3.05).
    #[test]
    fn closed_crossing_ring_keeps_every_segment() {
        let ring = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (2.0, -2.0)];
        let result = run(vec![closed_pline(&ring)], 0.09);
        assert_every_segment_carries_material(&result, &ring);
        // Swept area = perimeter · thickness, less the crossing overlap,
        // less the corner at (4, 4). Its legs meet at 18.4°, so its miter
        // reached 0.5625 — past `MITER_LIMIT * 0.09 = 0.36` — and the
        // corner bevels: the 0.0486 spike triangle it used to add is gone.
        let area = band_area(&result);
        assert!(
            (area - 3.004_751).abs() < 1e-4,
            "crossing loop band area={area}"
        );
        // Both lobes of the crossing are enclosed voids, not material.
        assert!(!band_covers(&result, (3.5, 1.0)), "upper lobe must be void");
        assert!(
            !band_covers(&result, (1.5, -0.5)),
            "lower lobe must be void"
        );
    }

    /// A figure-8 ring's lobes wind in opposite directions, so its signed
    /// area is exactly zero — the annulus model's larger-|area| outer /
    /// hole pick was a coin flip and it lost three of four segments.
    #[test]
    fn closed_figure_8_keeps_every_segment() {
        let ring = [(0.0, 0.0), (2.0, 2.0), (0.0, 2.0), (2.0, 0.0)];
        let result = run(vec![closed_pline(&ring)], 0.09);
        assert_every_segment_carries_material(&result, &ring);
        let area = band_area(&result);
        assert!((area - 1.705_834).abs() < 1e-4, "figure-8 band area={area}");
        // The crossing carries material; both lobes stay void.
        assert!(band_covers(&result, (1.0, 1.0)));
        assert!(!band_covers(&result, (1.0, 1.6)));
        assert!(!band_covers(&result, (1.0, 0.4)));
    }

    /// An asymmetric crossing loop (five legs, the last one cutting back
    /// across the first). The annulus model lost the returning leg.
    #[test]
    fn closed_asymmetric_crossing_keeps_every_segment() {
        let ring = [(0.0, 0.0), (8.0, 0.0), (8.0, 5.0), (3.0, 5.0), (3.0, -3.0)];
        let result = run(vec![closed_pline(&ring)], 0.09);
        assert_every_segment_carries_material(&result, &ring);
        let area = band_area(&result);
        assert!(
            (area - 5.411_275).abs() < 1e-4,
            "asymmetric crossing band area={area}"
        );
    }

    /// A ring that crosses itself repeatedly keeps every segment too.
    #[test]
    fn closed_multi_crossing_ring_keeps_every_segment() {
        let ring = [
            (0.0, 0.0),
            (4.0, 4.0),
            (1.0, 4.0),
            (4.0, 0.0),
            (3.0, 4.0),
            (0.0, 2.0),
        ];
        let result = run(vec![closed_pline(&ring)], 0.1);
        assert_every_segment_carries_material(&result, &ring);
    }

    /// CONTROL for the crossing fixtures above: a simple closed ring
    /// still bands to the same annulus, vertex for vertex and provenance
    /// for provenance. The keyhole assembly reproduces it because its
    /// slit is dropped and its two offset rings enter the arrangement
    /// unchanged. (The two miter points the slit touches gather two
    /// extra members in the vertex-snap cluster average, so they can
    /// land one ulp off — twelve orders of magnitude below `WALL_EPS`.)
    #[test]
    fn simple_closed_square_band_is_unchanged() {
        let square = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)];
        let result = run(vec![closed_pline(&square)], 0.3);
        assert_eq!(result.len(), 1);
        let (face, prov) = &result[0];
        assert_eq!(face.holes().len(), 1);
        let assert_pts = |pl: &Pline, want: &[(f64, f64)]| {
            let got: Vec<(f64, f64)> = pl.vertices.iter().map(|v| (v.x, v.y)).collect();
            assert_eq!(got.len(), want.len(), "got {got:?}");
            for (g, w) in got.iter().zip(want) {
                assert!(
                    (g.0 - w.0).abs() < 1e-12 && (g.1 - w.1).abs() < 1e-12,
                    "got {got:?}, want {want:?}"
                );
            }
        };
        assert_pts(
            face.outer(),
            &[(-0.3, 10.3), (-0.3, -0.3), (10.3, -0.3), (10.3, 10.3)],
        );
        assert_pts(
            &face.holes()[0],
            &[(0.3, 0.3), (0.3, 9.7), (9.7, 9.7), (9.7, 0.3)],
        );
        let side = |edge, side| SegmentProvenance {
            pline: 0,
            origin: SegmentOrigin::Side { edge, side },
            fragment: 0,
        };
        assert_eq!(
            prov.outer(),
            [
                side(3, OffsetSide::Right),
                side(0, OffsetSide::Right),
                side(1, OffsetSide::Right),
                side(2, OffsetSide::Right),
            ]
        );
        assert_eq!(
            prov.holes()[0],
            [
                side(3, OffsetSide::Left),
                side(2, OffsetSide::Left),
                side(1, OffsetSide::Left),
                side(0, OffsetSide::Left),
            ]
        );
        assert!((band_area(&result) - 24.0).abs() < 1e-9);
    }

    #[test]
    fn execute_faces_matches_provenance_variant_geometry() {
        let plines = vec![
            open_pline(&[(0.0, 0.0), (4.0, 0.0)]),
            open_pline(&[(2.0, 0.0), (2.0, 3.0)]),
        ];
        let plain = CurveBand2D::new(plines.clone(), 0.15)
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
