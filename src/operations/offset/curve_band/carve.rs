//! Cross-band carving: remove other bands' regions from a band's own
//! footprint faces, in 2D, before extrusion.
//!
//! [`carve_band_faces`] subtracts a set of cutter footprints from the
//! faces [`super::CurveBand2D::execute_faces_with_provenance`] produced,
//! and reports for **every** boundary segment of every surviving face
//! whether that segment is still a piece of the original band boundary
//! or a freshly cut seam — see [`CarvedSegmentProvenance`].
//!
//! The typical caller is a BIM modeller giving walls a priority order: a
//! lower-priority wall's footprint has the footprints of the
//! higher-priority walls it crosses removed, so the two solids meet
//! without overlapping. geolis stays identity-dumb — a cut segment
//! names its cutter only by **index into the `cutters` slice**, never by
//! a caller-supplied name or id.
//!
//! # Result semantics
//!
//! `result = ⋃ base_faces ∩ ¬⋃ cutters`, re-split into faces.
//!
//! - One base face may split into several output faces (a cutter cutting
//!   clean through it), or into none (a cutter covering it).
//! - Material of a base face lying under a **hole** of a cutter survives
//!   — a cutter's holes are empty space, exactly as they are for
//!   [`crate::operations::boolean_2d::subtract_all_with_holes`].
//! - An empty result is `Ok(vec![])`, not an error.
//!
//! Base faces are carved independently. Band footprint faces are the
//! connected components of one union and are therefore disjoint, so this
//! is equivalent to one global subtraction while letting faces no cutter
//! can reach skip the arrangement entirely.
//!
//! # How provenance is obtained
//!
//! Labels are **threaded** through the same arrangement engine the band
//! union uses: [`crate::operations::boolean_2d::subtract_all_with_holes_traced`]
//! reports the input ring edge behind every output edge, and its input
//! layout (`0` = base, `1 + c` = cutter `c`) makes the base/cut split
//! exact. No geometric matching is performed anywhere, so provenance
//! survives the engine's `WALL_EPS` vertex snapping unharmed.
//!
//! # Determinism and stability
//!
//! - Same input ⇒ same output, including face order and every fragment
//!   ordinal.
//! - Fragments are numbered by the shared rule in
//!   [`super::provenance::number_fragments`]: runs of consecutive edges
//!   from one source are ordered along that source and numbered densely
//!   from `0`.
//!   - **Inherited** segments order by their band fragment first, then
//!     by position within that fragment's run, then along the base edge
//!     — so a carve that splits one band fragment in two numbers the
//!     pieces in band order.
//!   - **Cut** segments of one cutter order by the cutter's own ring
//!     (outer, then holes), edge index, and parameter along that edge.
//! - Ordinals of a source MAY shift when the carve changes which pieces
//!   of it survive. This is inherent: ordinals are dense. The identity
//!   fast path below is exempt — it returns the input provenance
//!   verbatim.
//!
//! # Identity fast path
//!
//! Cutters whose bounding box cannot reach the base (with `WALL_EPS`
//! slack) are dropped up front. If none remains, the input faces and
//! their provenance are returned unchanged — wrapped in
//! [`CarvedSegmentProvenance::Base`] — without running the arrangement,
//! so neither the vertices nor the fragment ordinals move.

use crate::error::{OperationError, Result};
use crate::operations::boolean_2d::{subtract_all_with_holes_traced, RingRef, SegmentSite};

use super::polygon_union::{PolygonWithHoles, WALL_EPS};
use super::provenance::{number_fragments, FaceEdgeKeys, SourcePosition};
use super::{pline_xy, BandFootprint2D, FootprintProvenance, SegmentOrigin, SegmentProvenance};

/// Where one boundary segment of a carved footprint came from.
///
/// Wraps — rather than extends — [`SegmentProvenance`], so the band
/// union's provenance vocabulary is unchanged for its existing
/// consumers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CarvedSegmentProvenance {
    /// A surviving piece of the base band's own boundary. Carries the
    /// base segment's [`SegmentProvenance`] with `pline` and `origin`
    /// untouched; only `fragment` is renumbered, because the carve may
    /// have split the base segment into several pieces.
    Base(SegmentProvenance),
    /// A seam the carve created: the segment lies on the boundary of a
    /// cutter, and is interior to the un-carved band. Consumers use this
    /// to suppress seam outlines.
    Carve {
        /// Index of the cutter in the `cutters` slice passed to
        /// [`carve_band_faces`].
        cutter: usize,
        /// Deterministic ordinal of this seam piece among the seams that
        /// cutter contributed (see the module docs for the ordering).
        fragment: u32,
    },
}

/// Per-ring provenance aligned 1:1 with a carved [`BandFootprint2D`]:
/// `outer()[k]` describes the outer-ring segment from vertex `k` to
/// vertex `(k + 1) % n`, and `holes()[h][k]` likewise for hole `h`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CarvedFootprintProvenance {
    outer: Vec<CarvedSegmentProvenance>,
    holes: Vec<Vec<CarvedSegmentProvenance>>,
}

impl CarvedFootprintProvenance {
    /// Per-segment provenance of the outer ring, aligned with
    /// [`BandFootprint2D::outer`]'s segments.
    #[must_use]
    pub fn outer(&self) -> &[CarvedSegmentProvenance] {
        &self.outer
    }

    /// Per-segment provenance of each hole ring, aligned with
    /// [`BandFootprint2D::holes`].
    #[must_use]
    pub fn holes(&self) -> &[Vec<CarvedSegmentProvenance>] {
        &self.holes
    }

    /// Wrap an un-carved footprint's provenance verbatim (identity fast
    /// path): every segment is still a base segment, with its band
    /// fragment ordinal preserved.
    fn inherited(p: &FootprintProvenance) -> Self {
        let ring = |r: &[SegmentProvenance]| -> Vec<CarvedSegmentProvenance> {
            r.iter()
                .copied()
                .map(CarvedSegmentProvenance::Base)
                .collect()
        };
        Self {
            outer: ring(p.outer()),
            holes: p.holes().iter().map(|h| ring(h)).collect(),
        }
    }
}

/// Source identity of one carved output edge — the grouping key the
/// fragment numbering runs over. The `fragment` field of the base
/// segments is deliberately **not** part of the key: it is a position
/// along the source, not an identity, so pieces split out of one band
/// fragment stay in the same numbering group as their siblings.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum CarveKey {
    Base { pline: usize, origin: SegmentOrigin },
    Cut { cutter: usize },
}

/// Axis-aligned 2D bounds as `(min_x, min_y, max_x, max_y)`.
type Bounds2 = (f64, f64, f64, f64);

/// Subtract `cutters` from the band faces in `base`, preserving
/// per-segment provenance.
///
/// `base` is the pairing
/// [`super::CurveBand2D::execute_faces_with_provenance`] returns.
/// `cutters` are other bands' footprints — the type they are already
/// held in by callers that carve one wall against its neighbours.
///
/// Returns the surviving faces paired with provenance aligned to their
/// rings. The result may be shorter **or longer** than `base` (a cutter
/// may split a face or consume it entirely); an empty result is `Ok`.
/// Every returned face satisfies the [`BandFootprint2D`] winding and
/// simplicity contract, because it is assembled by the same arrangement
/// engine that guarantees it for [`super::CurveBand2D::execute_faces`].
///
/// # Validation
///
/// Both footprint types uphold their invariants by construction
/// ([`BandFootprint2D::try_from_parts`] or the union pipeline), so no
/// geometry is re-validated here. What *is* checked, in every build, is
/// the one thing the types cannot guarantee: that each supplied
/// [`FootprintProvenance`] is aligned with its footprint's rings.
///
/// # Errors
///
/// - [`OperationError::InvalidInput`] — a `base` entry's provenance does
///   not align 1:1 with its footprint's rings (wrong hole count, or a
///   ring whose segment count differs from its provenance length).
/// - [`OperationError::Failed`] — propagated from the arrangement engine
///   on degenerate input (bilateral classification still ambiguous after
///   ε exhaustion, broken parent topology, orientation/depth parity
///   violation), exactly as for
///   [`crate::operations::boolean_2d::subtract_all_with_holes`].
pub fn carve_band_faces(
    base: &[(BandFootprint2D, FootprintProvenance)],
    cutters: &[BandFootprint2D],
) -> Result<Vec<(BandFootprint2D, CarvedFootprintProvenance)>> {
    validate_provenance_alignment(base)?;

    let Some(base_bounds) = base
        .iter()
        .map(|(f, _)| footprint_bounds(f))
        .reduce(merge_bounds)
    else {
        // Nothing to carve.
        return Ok(Vec::new());
    };

    // bbox prefilter: a cutter that cannot reach any base face can only
    // contribute edges the fill oracle would discard. Dropping it here
    // keeps its index out of the arrangement without changing the result.
    let live: Vec<(usize, PolygonWithHoles)> = cutters
        .iter()
        .enumerate()
        .filter(|(_, c)| bounds_overlap(base_bounds, footprint_bounds(c)))
        .map(|(i, c)| (i, footprint_to_pwh(c)))
        .collect();

    if live.is_empty() {
        // Identity fast path — no arrangement, nothing renumbered.
        return Ok(base
            .iter()
            .map(|(f, p)| (f.clone(), CarvedFootprintProvenance::inherited(p)))
            .collect());
    }

    let mut out_faces: Vec<PolygonWithHoles> = Vec::new();
    let mut out_keys: Vec<FaceEdgeKeys<CarveKey>> = Vec::new();

    for (footprint, prov) in base {
        let base_pwh = footprint_to_pwh(footprint);
        let face_bounds = pwh_bounds(&base_pwh);

        // Narrow the cutter set once more against this single face, so a
        // face on the far side of the band skips the engine.
        let mut cut_index: Vec<usize> = Vec::new();
        let mut cut_pwhs: Vec<PolygonWithHoles> = Vec::new();
        for (global, pwh) in &live {
            if bounds_overlap(face_bounds, pwh_bounds(pwh)) {
                cut_index.push(*global);
                cut_pwhs.push(pwh.clone());
            }
        }

        let traced = subtract_all_with_holes_traced(&base_pwh, &cut_pwhs)?;

        // Position of each base ring edge inside its own band fragment,
        // which orders the pieces the carve splits that fragment into.
        let base_runs = base_run_positions(prov);

        for tf in traced {
            let resolve = |site: SegmentSite, p0: (f64, f64), p1: (f64, f64)| {
                resolve_edge(
                    site, p0, p1, &base_pwh, prov, &base_runs, &cut_index, &cut_pwhs,
                )
            };
            out_keys.push(FaceEdgeKeys {
                outer: ring_keys(&tf.face.outer, &tf.outer_sites, &resolve),
                holes: tf
                    .face
                    .holes
                    .iter()
                    .zip(&tf.hole_sites)
                    .map(|(pts, sites)| ring_keys(pts, sites, &resolve))
                    .collect(),
            });
            out_faces.push(tf.face);
        }
    }

    let fragments = number_fragments(&out_keys);

    Ok(out_faces
        .into_iter()
        .zip(&out_keys)
        .zip(fragments)
        .map(|((face, keys), frag)| {
            (
                BandFootprint2D::from_polygon_with_holes_unchecked(face),
                CarvedFootprintProvenance {
                    outer: ring_provenance(&keys.outer, &frag.outer),
                    holes: keys
                        .holes
                        .iter()
                        .zip(&frag.holes)
                        .map(|(k, f)| ring_provenance(k, f))
                        .collect(),
                },
            )
        })
        .collect())
}

/// `(key, position)` for every edge of one output ring.
fn ring_keys(
    pts: &[(f64, f64)],
    sites: &[SegmentSite],
    resolve: &impl Fn(SegmentSite, (f64, f64), (f64, f64)) -> (CarveKey, SourcePosition),
) -> Vec<(CarveKey, SourcePosition)> {
    debug_assert_eq!(pts.len(), sites.len());
    let m = sites.len();
    (0..m)
        .map(|e| resolve(sites[e], pts[e], pts[(e + 1) % m]))
        .collect()
}

/// Turn the numbered keys of one ring into its public provenance.
fn ring_provenance(
    keys: &[(CarveKey, SourcePosition)],
    fragments: &[u32],
) -> Vec<CarvedSegmentProvenance> {
    keys.iter()
        .zip(fragments)
        .map(|(&(key, _), &fragment)| match key {
            CarveKey::Base { pline, origin } => CarvedSegmentProvenance::Base(SegmentProvenance {
                pline,
                origin,
                fragment,
            }),
            CarveKey::Cut { cutter } => CarvedSegmentProvenance::Carve { cutter, fragment },
        })
        .collect()
}

/// Resolve one output edge's [`SegmentSite`] into its carve key and its
/// position along that source.
///
/// `site.input == 0` is a base ring edge (the traced subtract's input
/// layout); anything else indexes `cut_index` / `cut_pwhs`.
#[allow(
    clippy::too_many_arguments,
    reason = "one call site; the alternative is a context struct that only \
              re-groups the same borrows"
)]
fn resolve_edge(
    site: SegmentSite,
    p0: (f64, f64),
    p1: (f64, f64),
    base_pwh: &PolygonWithHoles,
    prov: &FootprintProvenance,
    base_runs: &[Vec<usize>],
    cut_index: &[usize],
    cut_pwhs: &[PolygonWithHoles],
) -> (CarveKey, SourcePosition) {
    if site.input == 0 {
        let ring = ring_ordinal(site.ring);
        let (a, b) = ring_edge(base_pwh, site.ring, site.edge);
        let sp = match site.ring {
            RingRef::Outer => prov.outer()[site.edge],
            RingRef::Hole(h) => prov.holes()[h][site.edge],
        };
        (
            CarveKey::Base {
                pline: sp.pline,
                origin: sp.origin,
            },
            (
                sp.fragment as usize,
                base_runs[ring][site.edge],
                min_param(a, b, p0, p1),
            ),
        )
    } else {
        let local = site.input - 1;
        let (a, b) = ring_edge(&cut_pwhs[local], site.ring, site.edge);
        (
            CarveKey::Cut {
                cutter: cut_index[local],
            },
            (ring_ordinal(site.ring), site.edge, min_param(a, b, p0, p1)),
        )
    }
}

/// `0` for the outer ring, `1 + h` for hole `h` — the ring order every
/// per-ring table in this module uses.
fn ring_ordinal(ring: RingRef) -> usize {
    match ring {
        RingRef::Outer => 0,
        RingRef::Hole(h) => 1 + h,
    }
}

/// Endpoints of edge `edge` of the named ring of `pwh`.
fn ring_edge(pwh: &PolygonWithHoles, ring: RingRef, edge: usize) -> ((f64, f64), (f64, f64)) {
    let pts = match ring {
        RingRef::Outer => &pwh.outer,
        RingRef::Hole(h) => &pwh.holes[h],
    };
    (pts[edge], pts[(edge + 1) % pts.len()])
}

/// Smaller of the two endpoints' projections onto `a → b`. Unnormalised
/// — monotonic along the supporting line, which is all ordering needs.
fn min_param(a: (f64, f64), b: (f64, f64), p0: (f64, f64), p1: (f64, f64)) -> f64 {
    let t = |p: (f64, f64)| (p.0 - a.0) * (b.0 - a.0) + (p.1 - a.1) * (b.1 - a.1);
    t(p0).min(t(p1))
}

/// Position of each ring edge within its maximal cyclic run of edges
/// sharing the same `(pline, origin, fragment)` — i.e. within the band
/// fragment it belongs to. Indexed `[ring_ordinal][edge]`.
///
/// The band union numbers one fragment per maximal run, so all edges of
/// one `(pline, origin, fragment)` triple in one ring are contiguous and
/// the walk below recovers exactly that run.
fn base_run_positions(prov: &FootprintProvenance) -> Vec<Vec<usize>> {
    let mut out = vec![run_positions(prov.outer())];
    out.extend(prov.holes().iter().map(|h| run_positions(h)));
    out
}

fn run_positions(ring: &[SegmentProvenance]) -> Vec<usize> {
    let m = ring.len();
    let mut pos = vec![0usize; m];
    if m == 0 {
        return pos;
    }
    let key = |e: usize| (ring[e].pline, ring[e].origin, ring[e].fragment);

    // Rotate the scan start to a fragment boundary so no run is split by
    // the ring's arbitrary index origin.
    let Some(start) = (0..m).find(|&e| key(e) != key((e + m - 1) % m)) else {
        // Whole ring is one fragment: ring order is run order.
        for (e, slot) in pos.iter_mut().enumerate() {
            *slot = e;
        }
        return pos;
    };

    let mut current = key(start);
    let mut counter = 0usize;
    for k in 0..m {
        let e = (start + k) % m;
        if key(e) != current {
            current = key(e);
            counter = 0;
        }
        pos[e] = counter;
        counter += 1;
    }
    pos
}

/// Reject a `base` entry whose provenance is not aligned with its
/// footprint's rings. Everything else about both inputs is guaranteed by
/// the footprint type's own constructors.
fn validate_provenance_alignment(base: &[(BandFootprint2D, FootprintProvenance)]) -> Result<()> {
    for (i, (footprint, prov)) in base.iter().enumerate() {
        let outer_segments = footprint.outer().vertices.len();
        if outer_segments != prov.outer().len() {
            return Err(OperationError::InvalidInput(format!(
                "carve_band_faces: base[{i}] outer has {outer_segments} segments \
                 but its provenance has {}",
                prov.outer().len()
            ))
            .into());
        }
        if footprint.holes().len() != prov.holes().len() {
            return Err(OperationError::InvalidInput(format!(
                "carve_band_faces: base[{i}] has {} holes but its provenance \
                 has {}",
                footprint.holes().len(),
                prov.holes().len()
            ))
            .into());
        }
        for (h, (hole, hole_prov)) in footprint.holes().iter().zip(prov.holes()).enumerate() {
            if hole.vertices.len() != hole_prov.len() {
                return Err(OperationError::InvalidInput(format!(
                    "carve_band_faces: base[{i}] hole[{h}] has {} segments but \
                     its provenance has {}",
                    hole.vertices.len(),
                    hole_prov.len()
                ))
                .into());
            }
        }
    }
    Ok(())
}

fn footprint_to_pwh(f: &BandFootprint2D) -> PolygonWithHoles {
    PolygonWithHoles {
        outer: pline_xy(f.outer()),
        holes: f.holes().iter().map(pline_xy).collect(),
    }
}

/// Bounds of a footprint. The outer ring encloses every hole, so it
/// alone bounds the face.
fn footprint_bounds(f: &BandFootprint2D) -> Bounds2 {
    f.outer()
        .vertices
        .iter()
        .fold(EMPTY_BOUNDS, |b, v| grow_bounds(b, (v.x, v.y)))
}

fn pwh_bounds(p: &PolygonWithHoles) -> Bounds2 {
    p.outer.iter().fold(EMPTY_BOUNDS, |b, &v| grow_bounds(b, v))
}

/// Inverted seed bounds: growing it by any point yields that point's
/// degenerate box, and it never overlaps anything.
const EMPTY_BOUNDS: Bounds2 = (
    f64::INFINITY,
    f64::INFINITY,
    f64::NEG_INFINITY,
    f64::NEG_INFINITY,
);

fn grow_bounds(b: Bounds2, p: (f64, f64)) -> Bounds2 {
    (b.0.min(p.0), b.1.min(p.1), b.2.max(p.0), b.3.max(p.1))
}

fn merge_bounds(a: Bounds2, b: Bounds2) -> Bounds2 {
    (a.0.min(b.0), a.1.min(b.1), a.2.max(b.2), a.3.max(b.3))
}

/// Bounds overlap with `WALL_EPS` slack, so a cutter merely touching the
/// base within snap tolerance is still fed to the arrangement.
fn bounds_overlap(a: Bounds2, b: Bounds2) -> bool {
    a.0 <= b.2 + WALL_EPS && b.0 <= a.2 + WALL_EPS && a.1 <= b.3 + WALL_EPS && b.1 <= a.3 + WALL_EPS
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::super::{CapEnd, CurveBand2D, OffsetSide};
    use super::*;
    use crate::geometry::pline::{Pline, PlineVertex};
    use crate::math::Point3;
    use crate::operations::boolean_2d::{intersect_all_with_holes, signed_area};

    /// Snap tolerance budget: the arrangement may move a vertex by up to
    /// `WALL_EPS`, so an area assertion over a footprint of this scale
    /// allows a few orders of magnitude more than that.
    const AREA_EPS: f64 = 1e-4;

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

    /// Base band (faces + provenance) from one open centerline.
    fn band(pts: &[(f64, f64)], hw: f64) -> Vec<(BandFootprint2D, FootprintProvenance)> {
        CurveBand2D::new(vec![open_pline(pts)], hw)
            .execute_faces_with_provenance()
            .expect("band must build")
    }

    /// Cutter faces from one open centerline.
    fn cutter_band(pts: &[(f64, f64)], hw: f64) -> Vec<BandFootprint2D> {
        CurveBand2D::new(vec![open_pline(pts)], hw)
            .execute_faces()
            .expect("band must build")
    }

    /// Cutter faces from one closed centerline — a ring, i.e. a cutter
    /// with a hole.
    fn cutter_ring(pts: &[(f64, f64)], hw: f64) -> Vec<BandFootprint2D> {
        CurveBand2D::new(vec![closed_pline(pts)], hw)
            .execute_faces()
            .expect("ring band must build")
    }

    fn pwh_area(p: &PolygonWithHoles) -> f64 {
        signed_area(&p.outer).abs() - p.holes.iter().map(|h| signed_area(h).abs()).sum::<f64>()
    }

    fn face_area(f: &BandFootprint2D) -> f64 {
        pwh_area(&footprint_to_pwh(f))
    }

    fn total_area(out: &[(BandFootprint2D, CarvedFootprintProvenance)]) -> f64 {
        out.iter().map(|(f, _)| face_area(f)).sum()
    }

    fn ring_pts(p: &Pline) -> Vec<(f64, f64)> {
        p.vertices.iter().map(|v| (v.x, v.y)).collect()
    }

    /// Every output face must independently satisfy the full twelve-check
    /// [`BandFootprint2D`] contract, and its provenance must align with
    /// its rings. Asserted by every fixture below.
    fn assert_faces_valid(out: &[(BandFootprint2D, CarvedFootprintProvenance)]) {
        for (i, (f, p)) in out.iter().enumerate() {
            BandFootprint2D::try_from_parts(f.outer().clone(), f.holes().to_vec())
                .unwrap_or_else(|e| panic!("output face {i} violates the footprint contract: {e}"));
            assert_eq!(
                f.outer().vertices.len(),
                p.outer().len(),
                "face {i}: outer provenance must align with the ring"
            );
            assert_eq!(f.holes().len(), p.holes().len(), "face {i}: hole count");
            for (h, hole) in f.holes().iter().enumerate() {
                assert_eq!(
                    hole.vertices.len(),
                    p.holes()[h].len(),
                    "face {i} hole {h}: provenance must align with the ring"
                );
            }
        }
    }

    /// One output boundary segment as `(midpoint, provenance)`.
    type LabelledSegment = ((f64, f64), CarvedSegmentProvenance);

    fn labelled_segments(
        out: &[(BandFootprint2D, CarvedFootprintProvenance)],
    ) -> Vec<LabelledSegment> {
        let mut segs = Vec::new();
        for (f, p) in out {
            let mut rings: Vec<(&Pline, &[CarvedSegmentProvenance])> = vec![(f.outer(), p.outer())];
            for (h, hole) in f.holes().iter().enumerate() {
                rings.push((hole, &p.holes()[h]));
            }
            for (ring, prov) in rings {
                let n = ring.vertices.len();
                for (e, sp) in prov.iter().enumerate() {
                    let a = &ring.vertices[e];
                    let b = &ring.vertices[(e + 1) % n];
                    segs.push((((a.x + b.x) * 0.5, (a.y + b.y) * 0.5), *sp));
                }
            }
        }
        segs
    }

    fn near(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-5
    }

    // ===== Cross-wall carving fixtures =====

    /// X crossing: a horizontal band cut clean through by a vertical one.
    /// Pins the split count, the seam labels, and the re-fragmentation of
    /// the surviving side offsets.
    #[test]
    fn x_crossing_splits_the_band_into_two_stubs() {
        // Base spans [0, 4] x [0.85, 1.15]; cutter spans
        // [1.85, 2.15] x [-0.15, 3.15], so it cuts the base's full depth.
        let base = band(&[(0.0, 1.0), (4.0, 1.0)], 0.15);
        assert_eq!(base.len(), 1, "one straight wall is one face");
        let cutters = cutter_band(&[(2.0, 0.0), (2.0, 3.0)], 0.15);
        assert_eq!(cutters.len(), 1);

        let out = carve_band_faces(&base, &cutters).expect("carve must succeed");
        assert_eq!(out.len(), 2, "a through-cut leaves two stubs");
        assert_faces_valid(&out);

        for (f, _) in &out {
            assert!(f.holes().is_empty(), "a stub is a plain rectangle");
            assert_eq!(f.outer().vertices.len(), 4);
        }
        let expected = 4.0 * 0.3 - 0.3 * 0.3;
        let area = total_area(&out);
        assert!(
            (area - expected).abs() < AREA_EPS,
            "area={area}, expected {expected}"
        );

        // Seams: the cutter's two side walls, one per stub, densely
        // numbered because each is its own run.
        let segments = labelled_segments(&out);
        let mut seam_fragments: Vec<u32> = Vec::new();
        // (fragment, midpoint x) of the surviving pieces of each side.
        let mut right_side: Vec<(u32, f64)> = Vec::new();
        let mut left_side: Vec<(u32, f64)> = Vec::new();
        let mut caps: Vec<(CapEnd, f64)> = Vec::new();
        for (mid, prov) in &segments {
            match *prov {
                CarvedSegmentProvenance::Carve { cutter, fragment } => {
                    assert_eq!(cutter, 0, "only cutter 0 exists");
                    assert!(
                        near(mid.0, 1.85) || near(mid.0, 2.15),
                        "seam must lie on a cutter wall; got {mid:?}"
                    );
                    seam_fragments.push(fragment);
                }
                CarvedSegmentProvenance::Base(sp) => {
                    assert_eq!(sp.pline, 0);
                    // Surviving edges keep the exact Side / Cap origin the
                    // band union gave them.
                    if near(mid.1, 0.85) {
                        assert_eq!(
                            sp.origin,
                            SegmentOrigin::Side {
                                edge: 0,
                                side: OffsetSide::Right
                            }
                        );
                        right_side.push((sp.fragment, mid.0));
                    } else if near(mid.1, 1.15) {
                        assert_eq!(
                            sp.origin,
                            SegmentOrigin::Side {
                                edge: 0,
                                side: OffsetSide::Left
                            }
                        );
                        left_side.push((sp.fragment, mid.0));
                    } else if near(mid.0, 0.0) {
                        assert_eq!(sp.origin, SegmentOrigin::Cap { end: CapEnd::Start });
                        caps.push((CapEnd::Start, mid.0));
                    } else if near(mid.0, 4.0) {
                        assert_eq!(sp.origin, SegmentOrigin::Cap { end: CapEnd::End });
                        caps.push((CapEnd::End, mid.0));
                    } else {
                        panic!("unexpected surviving base segment at {mid:?}: {sp:?}");
                    }
                }
            }
        }
        seam_fragments.sort_unstable();
        assert_eq!(seam_fragments, vec![0, 1], "two seams, densely numbered");
        assert_eq!(caps.len(), 2, "both end caps survive");

        // The Right offset runs +x along the ring, so its fragments are
        // numbered left stub first; the Left offset runs -x, so the right
        // stub takes fragment 0.
        right_side.sort_by(|a, b| a.0.cmp(&b.0));
        left_side.sort_by(|a, b| a.0.cmp(&b.0));
        assert_eq!(right_side.len(), 2);
        assert_eq!(left_side.len(), 2);
        assert_eq!(right_side[0].0, 0);
        assert_eq!(right_side[1].0, 1);
        assert!(
            right_side[0].1 < right_side[1].1,
            "Right fragments ascend with x: {right_side:?}"
        );
        assert!(
            left_side[0].1 > left_side[1].1,
            "Left fragments descend with x (ring direction): {left_side:?}"
        );
    }

    /// T junction: a cutter that stops inside the base leaves one face
    /// with a notch. The notch's three edges are contiguous, so they are
    /// one seam fragment.
    #[test]
    fn t_junction_cuts_a_notch_into_a_single_face() {
        // The flat end cap sits flush at the centerline endpoint, so the
        // cutter spans [1.85, 2.15] x [1.05, 3.0]. The base's top edge is
        // at y = 1.15, so the notch is 0.10 deep and 0.20 of material
        // still bridges it.
        let base = band(&[(0.0, 1.0), (4.0, 1.0)], 0.15);
        let cutters = cutter_band(&[(2.0, 1.05), (2.0, 3.0)], 0.15);

        let out = carve_band_faces(&base, &cutters).expect("carve must succeed");
        assert_eq!(out.len(), 1, "a partial cut leaves one connected face");
        assert!(out[0].0.holes().is_empty(), "a notch is not a hole");
        assert_faces_valid(&out);
        assert_eq!(
            out[0].0.outer().vertices.len(),
            8,
            "four original corners plus the notch's four"
        );

        let expected = 4.0 * 0.3 - 0.3 * 0.10;
        let area = total_area(&out);
        assert!(
            (area - expected).abs() < AREA_EPS,
            "area={area}, expected {expected}"
        );

        let segments = labelled_segments(&out);
        let seams: Vec<_> = segments
            .iter()
            .filter_map(|(mid, p)| match *p {
                CarvedSegmentProvenance::Carve { cutter, fragment } => {
                    Some((*mid, cutter, fragment))
                }
                CarvedSegmentProvenance::Base(_) => None,
            })
            .collect();
        assert_eq!(seams.len(), 3, "the notch has three cut edges: {seams:?}");
        for (_, cutter, fragment) in &seams {
            assert_eq!(*cutter, 0);
            assert_eq!(
                *fragment, 0,
                "contiguous cut edges are one seam fragment: {seams:?}"
            );
        }

        // The notch splits the Left offset into two pieces, numbered along
        // the ring direction (-x).
        let mut left: Vec<(u32, f64)> = segments
            .iter()
            .filter_map(|(mid, p)| match *p {
                CarvedSegmentProvenance::Base(sp)
                    if sp.origin
                        == (SegmentOrigin::Side {
                            edge: 0,
                            side: OffsetSide::Left,
                        }) =>
                {
                    Some((sp.fragment, mid.0))
                }
                _ => None,
            })
            .collect();
        left.sort_by(|a, b| a.0.cmp(&b.0));
        assert_eq!(left.len(), 2, "the notch splits the left offset: {left:?}");
        assert_eq!(left[0].0, 0);
        assert_eq!(left[1].0, 1);
        assert!(left[0].1 > left[1].1, "left fragments descend with x");

        // The Right offset is untouched and stays a single fragment.
        let right: Vec<u32> = segments
            .iter()
            .filter_map(|(_, p)| match *p {
                CarvedSegmentProvenance::Base(sp)
                    if sp.origin
                        == (SegmentOrigin::Side {
                            edge: 0,
                            side: OffsetSide::Right,
                        }) =>
                {
                    Some(sp.fragment)
                }
                _ => None,
            })
            .collect();
        assert_eq!(right, vec![0], "the uncut right offset keeps fragment 0");
    }

    /// A cutter that cannot reach the base is dropped by the bbox
    /// prefilter: the faces and their provenance come back untouched,
    /// vertices included (no arrangement, so no snapping).
    #[test]
    fn non_overlapping_cutter_takes_the_identity_fast_path() {
        let base = band(&[(0.0, 0.0), (4.0, 0.0)], 0.15);
        let cutters = cutter_band(&[(20.0, 20.0), (24.0, 20.0)], 0.15);

        let out = carve_band_faces(&base, &cutters).expect("carve must succeed");
        assert_eq!(out.len(), base.len());
        assert_faces_valid(&out);

        for ((bf, bp), (of, op)) in base.iter().zip(&out) {
            assert_eq!(
                ring_pts(bf.outer()),
                ring_pts(of.outer()),
                "identity path must not perturb vertices"
            );
            assert_eq!(bf.holes().len(), of.holes().len());
            for (bh, oh) in bf.holes().iter().zip(of.holes()) {
                assert_eq!(ring_pts(bh), ring_pts(oh));
            }
            let inherited: Vec<CarvedSegmentProvenance> = bp
                .outer()
                .iter()
                .copied()
                .map(CarvedSegmentProvenance::Base)
                .collect();
            assert_eq!(
                op.outer(),
                inherited.as_slice(),
                "identity path must not renumber fragments"
            );
        }
    }

    #[test]
    fn empty_cutter_slice_is_the_identity_fast_path() {
        let base = band(&[(0.0, 0.0), (4.0, 0.0)], 0.15);
        let out = carve_band_faces(&base, &[]).expect("carve must succeed");
        assert_eq!(out.len(), 1);
        assert_eq!(ring_pts(out[0].0.outer()), ring_pts(base[0].0.outer()));
    }

    #[test]
    fn empty_base_yields_no_faces() {
        let cutters = cutter_band(&[(0.0, 0.0), (4.0, 0.0)], 0.15);
        let out = carve_band_faces(&[], &cutters).expect("carve must succeed");
        assert!(out.is_empty());
    }

    /// A cutter swallowing the base is `Ok` with zero faces, not an error.
    #[test]
    fn cutter_covering_the_base_leaves_no_faces() {
        let base = band(&[(0.0, 0.0), (4.0, 0.0)], 0.15);
        let cutters = cutter_band(&[(-2.0, 0.0), (6.0, 0.0)], 2.0);
        let out = carve_band_faces(&base, &cutters).expect("full coverage is Ok, not an error");
        assert!(out.is_empty(), "nothing survives, got {} faces", out.len());
    }

    /// Carving a band with its own footprint removes everything.
    #[test]
    fn identical_footprint_fully_carves_the_base() {
        let base = band(&[(0.0, 0.0), (4.0, 1.0), (7.0, 1.0)], 0.2);
        let cutters: Vec<BandFootprint2D> = base.iter().map(|(f, _)| f.clone()).collect();
        let out = carve_band_faces(&base, &cutters).expect("carve must succeed");
        assert!(
            out.is_empty(),
            "an identical cutter leaves nothing, got {} faces",
            out.len()
        );
    }

    /// A cutter with a hole is material only where it is filled: the base
    /// under the hole survives as a disconnected island.
    #[test]
    fn base_under_a_cutter_hole_survives_as_an_island() {
        // Base [0, 10] x [-3, 3]. Cutter ring: outer [1.7, 8.3] x
        // [-2.3, 2.3], hole [2.3, 7.7] x [-1.7, 1.7].
        let base = band(&[(0.0, 0.0), (10.0, 0.0)], 3.0);
        assert_eq!(base.len(), 1);
        let cutters = cutter_ring(&[(2.0, -2.0), (8.0, -2.0), (8.0, 2.0), (2.0, 2.0)], 0.3);
        assert_eq!(cutters.len(), 1);
        assert_eq!(cutters[0].holes().len(), 1, "the ring cutter has a hole");

        let out = carve_band_faces(&base, &cutters).expect("carve must succeed");
        assert_eq!(out.len(), 2, "surround plus island");
        assert_faces_valid(&out);

        let island = out
            .iter()
            .find(|(f, _)| f.holes().is_empty())
            .expect("the island has no holes");
        let surround = out
            .iter()
            .find(|(f, _)| f.holes().len() == 1)
            .expect("the surround has the cutter's outline as a hole");

        let island_area = face_area(&island.0);
        assert!(
            (island_area - 5.4 * 3.4).abs() < AREA_EPS,
            "island area={island_area}, expected {}",
            5.4 * 3.4
        );
        let surround_area = face_area(&surround.0);
        assert!(
            (surround_area - (60.0 - 6.6 * 4.6)).abs() < AREA_EPS,
            "surround area={surround_area}"
        );

        // The island's whole boundary is the cutter's hole ring: one run,
        // hence one fragment. The surround's hole is the cutter's outer
        // ring — also one run, and ordered first because cut fragments
        // order by the cutter's own ring (outer before holes).
        let island_labels: Vec<CarvedSegmentProvenance> = island.1.outer().to_vec();
        for label in &island_labels {
            assert_eq!(
                *label,
                CarvedSegmentProvenance::Carve {
                    cutter: 0,
                    fragment: 1
                },
                "island boundary must be the cutter's hole seam"
            );
        }
        for label in &surround.1.holes()[0] {
            assert_eq!(
                *label,
                CarvedSegmentProvenance::Carve {
                    cutter: 0,
                    fragment: 0
                },
                "the surround's hole is the cutter's outer seam"
            );
        }
        // The surround's outer ring is untouched base boundary.
        for label in surround.1.outer() {
            assert!(
                matches!(label, CarvedSegmentProvenance::Base(_)),
                "the surround's outer ring is base boundary; got {label:?}"
            );
        }
    }

    /// Curved base: the footprint of a bulge centerline is a tessellated
    /// polygon, and the carve operates on it like any other footprint.
    /// The area is anchored analytically by the subtract/intersect
    /// partition invariant.
    #[test]
    fn curved_base_band_carves_on_its_tessellated_footprint() {
        // Arc from (0,0) to (4,0), bulge 0.5 — the same fixture the
        // provenance suite uses for curved centerlines.
        let arc = Pline {
            vertices: vec![PlineVertex::new(0.0, 0.0, 0.5), PlineVertex::line(4.0, 0.0)],
            closed: false,
        };
        let base = CurveBand2D::new(vec![arc], 0.3)
            .execute_faces_with_provenance()
            .expect("curved band must build");
        assert_eq!(base.len(), 1);
        assert!(
            base[0].0.outer().vertices.len() > 8,
            "the arc must tessellate into many chords"
        );
        let base_area = face_area(&base[0].0);

        // Vertical through-cut across the arc's crown.
        let cutters = cutter_band(&[(2.0, -3.0), (2.0, 3.0)], 0.2);
        let out = carve_band_faces(&base, &cutters).expect("carve must succeed");
        assert_eq!(out.len(), 2, "the cut splits the arc band in two");
        assert_faces_valid(&out);

        // Partition invariant: what the carve kept plus what a plain
        // intersect selects must be the whole base.
        let removed = intersect_all_with_holes(
            &footprint_to_pwh(&base[0].0),
            &[footprint_to_pwh(&cutters[0])],
        )
        .expect("intersect must succeed");
        let removed_area: f64 = removed.iter().map(pwh_area).sum();
        assert!(removed_area > 0.1, "the cut must remove real material");
        let kept = total_area(&out);
        assert!(
            (kept + removed_area - base_area).abs() < AREA_EPS,
            "partition invariant violated: kept={kept}, removed={removed_area}, base={base_area}"
        );

        // Every chord of the arc keeps centerline edge 0; every seam names
        // cutter 0 and lies on one of its walls.
        for (mid, prov) in labelled_segments(&out) {
            match prov {
                CarvedSegmentProvenance::Base(sp) => {
                    assert_eq!(sp.pline, 0);
                    match sp.origin {
                        SegmentOrigin::Side { edge, .. } => {
                            assert_eq!(edge, 0, "all chords map back to the arc edge");
                        }
                        SegmentOrigin::Cap { .. } => {}
                    }
                }
                CarvedSegmentProvenance::Carve { cutter, .. } => {
                    assert_eq!(cutter, 0);
                    assert!(
                        near(mid.0, 1.8) || near(mid.0, 2.2),
                        "seam must lie on a cutter wall; got {mid:?}"
                    );
                }
            }
        }
    }

    /// Two cutters are reported by their index in the caller's slice —
    /// including after the bbox prefilter drops one in between.
    #[test]
    fn cut_seams_name_their_index_in_the_caller_slice() {
        let base = band(&[(0.0, 1.0), (12.0, 1.0)], 0.15);
        let mut cutters = cutter_band(&[(3.0, 0.0), (3.0, 3.0)], 0.15);
        // Cutter 1 is far away and is dropped by the prefilter; cutter 2
        // must still report index 2.
        cutters.extend(cutter_band(&[(50.0, 50.0), (50.0, 53.0)], 0.15));
        cutters.extend(cutter_band(&[(9.0, 0.0), (9.0, 3.0)], 0.15));
        assert_eq!(cutters.len(), 3);

        let out = carve_band_faces(&base, &cutters).expect("carve must succeed");
        assert_eq!(out.len(), 3, "two through-cuts leave three pieces");
        assert_faces_valid(&out);

        let mut seen: Vec<usize> = labelled_segments(&out)
            .iter()
            .filter_map(|(_, p)| match *p {
                CarvedSegmentProvenance::Carve { cutter, .. } => Some(cutter),
                CarvedSegmentProvenance::Base(_) => None,
            })
            .collect();
        seen.sort_unstable();
        seen.dedup();
        assert_eq!(
            seen,
            vec![0, 2],
            "the dropped cutter must not shift indices"
        );
    }

    #[test]
    fn carve_is_deterministic() {
        let base = band(&[(0.0, 1.0), (4.0, 1.0), (4.0, 5.0)], 0.15);
        let cutters = cutter_band(&[(2.0, 0.0), (2.0, 3.0)], 0.2);
        let first = carve_band_faces(&base, &cutters).expect("carve must succeed");
        let second = carve_band_faces(&base, &cutters).expect("carve must succeed");

        assert_eq!(first.len(), second.len());
        for ((fa, pa), (fb, pb)) in first.iter().zip(&second) {
            assert_eq!(ring_pts(fa.outer()), ring_pts(fb.outer()));
            assert_eq!(fa.holes().len(), fb.holes().len());
            for (ha, hb) in fa.holes().iter().zip(fb.holes()) {
                assert_eq!(ring_pts(ha), ring_pts(hb));
            }
            assert_eq!(pa, pb, "provenance including ordinals must be stable");
        }
    }

    // ===== Public-API validation =====

    #[test]
    fn misaligned_outer_provenance_is_rejected() {
        let mut base = band(&[(0.0, 1.0), (4.0, 1.0)], 0.15);
        base[0].1.outer.pop();
        let cutters = cutter_band(&[(2.0, 0.0), (2.0, 3.0)], 0.15);
        let err = carve_band_faces(&base, &cutters).expect_err("misalignment must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("provenance"), "{msg}");
    }

    #[test]
    fn misaligned_hole_count_is_rejected() {
        let mut base = CurveBand2D::new(
            vec![closed_pline(&[
                (0.0, 0.0),
                (10.0, 0.0),
                (10.0, 10.0),
                (0.0, 10.0),
            ])],
            0.3,
        )
        .execute_faces_with_provenance()
        .expect("ring band must build");
        assert_eq!(base[0].0.holes().len(), 1);
        base[0].1.holes.clear();
        let cutters = cutter_band(&[(5.0, -1.0), (5.0, 11.0)], 0.2);
        let err = carve_band_faces(&base, &cutters).expect_err("misalignment must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("holes"), "{msg}");
    }
}
