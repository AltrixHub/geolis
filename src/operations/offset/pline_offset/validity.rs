//! Result validity for an offset — closed ring or open path.
//!
//! The offset of a source by `d` is, by definition, carried by the points
//! at distance `|d|` from it. A result is therefore checked against that
//! distance from BOTH sides, and each side names a different artifact:
//!
//! | A result point | Is | Because |
//! |---|---|---|
//! | CLOSER than `\|d\|` (a closed RING only) | a **phantom** — the artifact a mitred join produces once it crosses the ring's medial axis | nothing at distance `< \|d\|` is inside the region a ring's offset bounds |
//! | FURTHER than [`RawOffset::MITER_LIMIT`] · `\|d\|` | a **runaway** — an offset segment that came back with the complementary sweep, so its arc takes the long way round and swings clear of the source | a join miters only while it stays within that of its own source vertex, and everything else sits at exactly `\|d\|`; nothing on a well-formed offset is further |
//!
//! [`super::filter`] already discards phantoms, but only where the RAW
//! offset self-intersects: that is the only path that slices and filters.
//! A source whose raw offset does NOT self-intersect skips every check —
//! an axis-aligned rectangle inset past its own half-extent came back
//! inside-out that way, and a four-legged arc chain offset past half its
//! leg length came back with an arc bulged the long way round, putting a
//! guard's rail line 3.8 m off a run it is meant to sit 0.73 m inside.
//! Applying the criterion to the finished result closes both holes for
//! every caller, on every path.
//!
//! Distances are measured against the true source geometry (arcs
//! included, via [`super::filter::min_dist_to_pline`]), so an arc-carrying
//! source needs no sampling slack.

use crate::geometry::pline::Pline;
use crate::math::arc_2d::arc_from_bulge;

use super::filter::min_dist_to_pline;
use super::RawOffset;

/// Relative slack on the distance criterion, so a result sitting at
/// exactly `|d|` (every corner of a well-formed mitred offset does)
/// survives floating-point noise.
const REL_TOLERANCE: f64 = 1e-6;

/// Keeps the LOOPS that are genuine offsets of `original` at `distance`
/// — no phantom, no runaway (see the module docs); a loop with fewer
/// than 3 vertices encloses no area and is dropped with them.
#[must_use]
pub fn keep_valid(loops: Vec<Pline>, original: &Pline, distance: f64) -> Vec<Pline> {
    let abs_d = distance.abs();
    let nearest = abs_d - tolerance(abs_d);
    keep_within(loops, original, distance, 3, nearest)
}

/// …and the RUNAWAY half alone for the open PATHS an open source
/// offsets to, which bound no area and so need only their two ends.
///
/// The phantom half cannot be asked of an open source: it is a claim
/// about a boundary, and an open path is not one. A chain that wraps
/// back past itself — two arc legs turning ninety degrees into each
/// other, the ordinary switchback — legitimately offsets to a line that
/// runs within centimetres of an EARLIER leg of the same source, and
/// there is nothing wrong with it: the offset of leg two is measured
/// from leg two. Only the upper bound survives the generalisation, and
/// it is the one that catches the runaway.
#[must_use]
pub fn keep_valid_open(paths: Vec<Pline>, original: &Pline, distance: f64) -> Vec<Pline> {
    keep_within(paths, original, distance, 2, 0.0)
}

/// Relative slack for a given distance.
fn tolerance(abs_d: f64) -> f64 {
    REL_TOLERANCE * abs_d.max(1.0)
}

/// The shared criterion: every probe point of a candidate lies between
/// `nearest` and ONE bound shared by the whole candidate — the furthest
/// any join is allowed to miter, [`RawOffset::MITER_LIMIT`] · `|distance|`
/// plus slack. Both are measured from the source.
fn keep_within(
    candidates: Vec<Pline>,
    original: &Pline,
    distance: f64,
    least_vertices: usize,
    nearest: f64,
) -> Vec<Pline> {
    let abs_d = distance.abs();
    let furthest = RawOffset::MITER_LIMIT * abs_d + tolerance(abs_d);
    candidates
        .into_iter()
        .filter(|candidate| {
            candidate.vertices.len() >= least_vertices
                && probe_points(candidate).all(|(px, py)| {
                    let d = min_dist_to_pline(px, py, original);
                    d >= nearest && d <= furthest
                })
        })
        .collect()
}

/// The points a candidate is measured at — a closed loop or an open
/// path alike: every vertex plus every segment's midpoint (the arc
/// midpoint for a bulge-carrying segment, which bows away from its
/// chord and is where an arc join lands closest to the source).
///
/// An open path has one segment fewer than it has vertices, so the
/// wrap-around below only ever fires for a loop.
fn probe_points(path: &Pline) -> impl Iterator<Item = (f64, f64)> + '_ {
    let n = path.vertices.len();
    let vertices = path.vertices.iter().map(|v| (v.x, v.y));
    let midpoints = (0..path.segment_count()).map(move |i| {
        let v0 = &path.vertices[i];
        let v1 = &path.vertices[(i + 1) % n];
        if v0.bulge.abs() < 1e-12 {
            ((v0.x + v1.x) * 0.5, (v0.y + v1.y) * 0.5)
        } else {
            let (cx, cy, r, start, sweep) = arc_from_bulge(v0.x, v0.y, v1.x, v1.y, v0.bulge);
            let mid = start + sweep * 0.5;
            (cx + r * mid.cos(), cy + r * mid.sin())
        }
    });
    vertices.chain(midpoints)
}
