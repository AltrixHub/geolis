//! Bounds index over a fill oracle's inputs.
//!
//! Half-edge classification is the arrangement engine's hot loop: every
//! sub-edge probes the oracle four times (two directions × bilateral
//! sampling), and every probe asked *each* input for a full
//! [`point_in_polygon_class`](super::types::point_in_polygon_class) walk
//! — two passes over the ring, one of them with a `sqrt` per edge. That
//! made classification `O(sub-edges × total input vertices)`, which is
//! why a floor plan's fused outline got quadratically slower as elements
//! accumulated on it.
//!
//! Almost all of that work is provably wasted: a probe near one wall band
//! is nowhere near the other thirty. This index answers "can input `i`
//! possibly matter at `p`?" in four comparisons.
//!
//! # Why the rejection is exact
//!
//! [`Bounds::of`] takes the outer ring's axis-aligned extent and grows it
//! by at least [`WALL_EPS`] (see [`grown_extent`]). If `p` falls outside
//! that box then, for every ring of the input (outer *and* holes, since a
//! hole is contained in its outer):
//!
//! - **No edge is within `WALL_EPS`.** Every ring vertex lies inside the
//!   un-grown box, so every point of every edge does too; `p` is further
//!   than `WALL_EPS` from the box, hence further than `WALL_EPS` from any
//!   edge. So no ring can report [`PointClass::Boundary`](super::types::PointClass::Boundary).
//! - **The winding number is zero.** `p` lies outside the ring's convex
//!   hull, so a ray from it crosses the closed ring an even number of
//!   times. So the outer reports [`PointClass::Outside`](super::types::PointClass::Outside).
//!
//! An input whose outer is `Outside` is not filled and touches no ring —
//! exactly the contribution a skipped input would have made. The
//! rejection therefore changes no classification, only how long it takes
//! to reach it.

use super::types::{Polygon, WALL_EPS};

/// How far a locality prefilter must grow an extent at coordinate `v` so
/// that "outside the grown extent" still implies "further than `reach`
/// away".
///
/// The requested `reach` alone is not enough. The arrangement accepts coordinates up
/// to `MAX_ARRANGEMENT_COORD` (1e12), where one ulp is ~1.2e-4 — six
/// orders of magnitude COARSER than `WALL_EPS`. There `v - WALL_EPS`
/// rounds back to `v`, the extent is not grown at all, and a pair that is
/// genuinely within tolerance gets rejected: a dropped split point, a
/// corrupted arrangement, and a node that silently publishes nothing.
/// (That is exactly what the extreme-bulge wall chain hit — an arc whose
/// radius runs off to the coordinate ceiling.)
///
/// Padding by several ulps of the coordinate as well keeps the growth
/// real at every magnitude. Over-growing is always safe — it only makes
/// the filter reject fewer candidates — so taking the max of the three is
/// both sound and tight where it matters (building scale, where the
/// caller's own tolerance dominates).
fn eps_pad(v: f64, reach: f64) -> f64 {
    reach.max(WALL_EPS).max(v.abs() * (8.0 * f64::EPSILON))
}

/// `(min, max)` of `lo..=hi`, grown so that anything outside the result is
/// further than `reach` — and never less than [`WALL_EPS`] or a few ulps
/// of the coordinate, whichever is coarser.
pub(crate) fn grown_extent_reaching(lo: f64, hi: f64, reach: f64) -> (f64, f64) {
    (lo - eps_pad(lo, reach), hi + eps_pad(hi, reach))
}

/// [`grown_extent_reaching`] for a consumer whose tolerance is plain
/// [`WALL_EPS`] in world units — the fill oracles' case.
pub(crate) fn grown_extent(lo: f64, hi: f64) -> (f64, f64) {
    grown_extent_reaching(lo, hi, WALL_EPS)
}

/// How far outside its own extent a segment of length `len` can still
/// take part in a split, in world units.
///
/// The split helpers do NOT all work in world units, which is what makes
/// this more than `WALL_EPS`:
///
/// | Helper | Its tolerance | World reach |
/// |---|---|---|
/// | `seg_seg_intersect` | `t, u ∈ (-WALL_EPS, 1 + WALL_EPS)` — **parameter** space | `WALL_EPS · len` |
/// | `collinear_overlap_params` | `\|d × (b0 - a0)\| < WALL_EPS` — a **cross product**, i.e. `len ·` distance | `WALL_EPS / len` |
/// | `project_endpoint_on_interior` | perpendicular distance `< WALL_EPS` | `WALL_EPS` |
///
/// Both of the first two blow up away from unit length, in opposite
/// directions, so the reach is `WALL_EPS · max(1, len, 1/len)`. On
/// building geometry (`len` within a couple of decades of 1 m) that is
/// still micrometres and the prefilter loses nothing; on the degenerate
/// end — a wall arc whose radius runs to the coordinate ceiling, where
/// chords are kilometres long — it correctly stops rejecting.
///
/// A segment shorter than `WALL_EPS` is exempt from the `1/len` term:
/// both helpers that need it bail on `len_sq < WALL_EPS_SQ` before
/// computing anything, so only the parameter-space reach applies.
pub(crate) fn split_reach(len: f64) -> f64 {
    if !len.is_finite() || len < WALL_EPS {
        return WALL_EPS;
    }
    WALL_EPS * len.max(1.0 / len).max(1.0)
}

/// Axis-aligned extent of one input, already grown by at least
/// [`WALL_EPS`] (see [`grown_extent`]).
#[derive(Debug, Clone, Copy)]
pub(crate) struct Bounds {
    min_x: f64,
    min_y: f64,
    max_x: f64,
    max_y: f64,
}

impl Bounds {
    /// Bounds of `poly`, grown by [`grown_extent`]. `None` when the ring
    /// has no vertices or carries a non-finite coordinate — in both cases
    /// the caller must fall back to the exact classifier rather than
    /// reject.
    fn of(poly: &Polygon) -> Option<Self> {
        let mut b = Self {
            min_x: f64::INFINITY,
            min_y: f64::INFINITY,
            max_x: f64::NEG_INFINITY,
            max_y: f64::NEG_INFINITY,
        };
        for &(x, y) in poly {
            if !x.is_finite() || !y.is_finite() {
                return None;
            }
            b.min_x = b.min_x.min(x);
            b.min_y = b.min_y.min(y);
            b.max_x = b.max_x.max(x);
            b.max_y = b.max_y.max(y);
        }
        if b.min_x > b.max_x {
            return None;
        }
        (b.min_x, b.max_x) = grown_extent(b.min_x, b.max_x);
        (b.min_y, b.max_y) = grown_extent(b.min_y, b.max_y);
        Some(b)
    }

    fn contains(&self, p: (f64, f64)) -> bool {
        p.0 >= self.min_x && p.0 <= self.max_x && p.1 >= self.min_y && p.1 <= self.max_y
    }
}

/// Per-input bounds for a fill oracle.
///
/// `None` at an index means "never reject this input" — the ring gave no
/// usable extent, so only the exact classifier may speak for it.
#[derive(Debug)]
pub(crate) struct BoundsIndex {
    bounds: Vec<Option<Bounds>>,
}

impl BoundsIndex {
    pub(crate) fn new<'a>(outers: impl Iterator<Item = &'a Polygon>) -> Self {
        Self {
            bounds: outers.map(Bounds::of).collect(),
        }
    }

    /// Whether input `index` can possibly classify `p` as anything other
    /// than [`Outside`](super::types::PointClass::Outside). Conservative:
    /// an unknown index or an
    /// input with no usable bounds always answers `true`.
    pub(crate) fn may_contain(&self, index: usize, p: (f64, f64)) -> bool {
        match self.bounds.get(index) {
            Some(Some(b)) => b.contains(p),
            _ => true,
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::super::types::{point_in_polygon_class, PointClass};
    use super::*;

    fn square(cx: f64, cy: f64, r: f64) -> Polygon {
        vec![
            (cx - r, cy - r),
            (cx + r, cy - r),
            (cx + r, cy + r),
            (cx - r, cy + r),
        ]
    }

    /// The rejection must never disagree with the exact classifier: every
    /// point the index rejects must genuinely be `Outside`.
    #[test]
    fn rejection_agrees_with_the_exact_classifier() {
        let poly = square(0.0, 0.0, 1.0);
        let index = BoundsIndex::new(std::iter::once(&poly));
        let mut rejected = 0;
        for i in -40i32..=40 {
            for j in -40i32..=40 {
                let p = (f64::from(i) * 0.1, f64::from(j) * 0.1);
                if index.may_contain(0, p) {
                    continue;
                }
                rejected += 1;
                assert_eq!(
                    point_in_polygon_class(p, &poly),
                    PointClass::Outside,
                    "index rejected {p:?} but the classifier disagrees",
                );
            }
        }
        assert!(rejected > 0, "fixture must exercise the rejection path");
    }

    /// A point inside the boundary band is never rejected — that is the
    /// case the `WALL_EPS` growth exists for.
    #[test]
    fn the_boundary_band_is_never_rejected() {
        let poly = square(0.0, 0.0, 1.0);
        let index = BoundsIndex::new(std::iter::once(&poly));
        for d in [0.0, WALL_EPS * 0.5, WALL_EPS * 0.99] {
            assert!(index.may_contain(0, (1.0 + d, 0.0)));
            assert!(index.may_contain(0, (0.0, -1.0 - d)));
        }
    }

    /// Deterministic xorshift — the randomised fixture below must run
    /// the same sequence on every machine and every CI run.
    struct Rng(u64);

    impl Rng {
        fn next_u64(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }

        /// A coordinate quantised to 1/32, so probes land exactly on ring
        /// vertices and edges often enough to exercise the boundary band.
        fn coord(&mut self) -> f64 {
            let q = u32::try_from(self.next_u64() % 1024).unwrap_or(0);
            f64::from(q) / 32.0 - 16.0
        }
    }

    /// The soundness premise, over random rings and random probes: an
    /// input the index rejects is genuinely `PointClass::Outside` —
    /// never `Inside`, never `Boundary`. If this fails, the oracles are
    /// skipping an input that would have changed the fill verdict, and
    /// every boolean built on the engine is silently wrong.
    #[test]
    fn a_rejected_input_is_always_outside() {
        let mut rng = Rng(0x1234_5678_9abc_def0);
        let mut rejected = 0usize;
        for _ in 0..20_000 {
            // A random quadrilateral — convex or not, the classifier and
            // the bounds must agree either way.
            let poly: Polygon = (0..4).map(|_| (rng.coord(), rng.coord())).collect();
            let index = BoundsIndex::new(std::iter::once(&poly));
            for _ in 0..8 {
                let p = (rng.coord(), rng.coord());
                if index.may_contain(0, p) {
                    continue;
                }
                rejected += 1;
                assert_eq!(
                    point_in_polygon_class(p, &poly),
                    PointClass::Outside,
                    "index rejected {p:?} against {poly:?}",
                );
            }
        }
        assert!(
            rejected > 10_000,
            "fixture must exercise the rejection path; got {rejected}",
        );
    }

    /// The growth must be REAL at every magnitude the arrangement
    /// accepts: at the 1e12 coordinate ceiling one ulp is ~1.2e-4, so a
    /// bare `WALL_EPS` subtraction rounds away and the extent silently
    /// stops covering the tolerance band.
    #[test]
    fn the_extent_grows_at_every_accepted_magnitude() {
        for v in [0.0, 1.0, 1.0e3, 1.0e6, 1.0e9, 1.0e12, -1.0e12] {
            let (lo, hi) = grown_extent(v, v);
            assert!(lo < v, "extent low end did not grow below {v:e}");
            assert!(hi > v, "extent high end did not grow above {v:e}");
        }
    }

    /// And a ring at the coordinate ceiling still admits probes inside
    /// its band — the property the extreme-bulge wall chain depends on.
    #[test]
    fn a_ring_at_the_coordinate_ceiling_keeps_its_band() {
        let poly = square(1.0e12, 1.0e12, 1.0);
        let index = BoundsIndex::new(std::iter::once(&poly));
        assert!(index.may_contain(0, (1.0e12 + 1.0, 1.0e12)));
        assert!(index.may_contain(0, (1.0e12, 1.0e12 - 1.0)));
    }

    /// A degenerate or non-finite ring yields no bounds, so it is always
    /// handed to the exact classifier.
    #[test]
    fn an_unusable_ring_is_never_rejected() {
        let empty: Polygon = Vec::new();
        let nan: Polygon = vec![(0.0, 0.0), (f64::NAN, 1.0), (1.0, 1.0)];
        let index = BoundsIndex::new([&empty, &nan].into_iter());
        assert!(index.may_contain(0, (1e9, 1e9)));
        assert!(index.may_contain(1, (1e9, 1e9)));
        // Out-of-range indices are conservative too.
        assert!(index.may_contain(7, (1e9, 1e9)));
    }
}
