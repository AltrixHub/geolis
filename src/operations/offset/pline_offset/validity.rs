//! Result validity for a closed-ring offset.
//!
//! The offset of a closed ring by `d` is, by definition, the boundary of
//! the set of points at distance `≥ |d|` from the source ring (inside it
//! for an inward offset, outside for an outward one). A result loop
//! carrying a point CLOSER than `|d|` to the source is therefore not an
//! offset of anything: it is a **phantom loop** — the artifact a mitred
//! join produces once it crosses the ring's medial axis.
//!
//! [`super::filter`] already discards such loops, but only for rings whose
//! RAW offset self-intersects: that is the only path that slices and
//! filters. A convex ring — an axis-aligned rectangle inset past its own
//! half-extent being the canonical case — produces a raw offset that never
//! self-intersects, so it skipped every check and the caller received an
//! inside-out rectangle. Applying the same distance criterion to the
//! finished loops closes that hole for every caller, on both paths.
//!
//! Distances are measured against the true source geometry (arcs
//! included, via [`super::filter::min_dist_to_pline`]), so an arc-carrying
//! ring needs no sampling slack.

use crate::geometry::pline::Pline;
use crate::math::arc_2d::arc_from_bulge;

use super::filter::min_dist_to_pline;

/// Relative slack on the distance criterion, so a result sitting at
/// exactly `|d|` (every corner of a well-formed mitred offset does)
/// survives floating-point noise.
const REL_TOLERANCE: f64 = 1e-6;

/// Keeps the loops that are genuine offsets of `original` at `distance`
/// (see the module docs); a loop with fewer than 3 vertices encloses no
/// area and is dropped with them.
#[must_use]
pub fn keep_valid(loops: Vec<Pline>, original: &Pline, distance: f64) -> Vec<Pline> {
    let abs_d = distance.abs();
    let tolerance = REL_TOLERANCE * abs_d.max(1.0);
    loops
        .into_iter()
        .filter(|candidate| {
            candidate.vertices.len() >= 3
                && probe_points(candidate)
                    .all(|(px, py)| min_dist_to_pline(px, py, original) + tolerance >= abs_d)
        })
        .collect()
}

/// The points a loop is measured at: every vertex plus every segment's
/// midpoint (the arc midpoint for a bulge-carrying segment, which bows
/// away from its chord and is where an arc join lands closest to the
/// source).
fn probe_points(ring: &Pline) -> impl Iterator<Item = (f64, f64)> + '_ {
    let n = ring.vertices.len();
    let vertices = ring.vertices.iter().map(|v| (v.x, v.y));
    let midpoints = (0..ring.segment_count()).map(move |i| {
        let v0 = &ring.vertices[i];
        let v1 = &ring.vertices[(i + 1) % n];
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
