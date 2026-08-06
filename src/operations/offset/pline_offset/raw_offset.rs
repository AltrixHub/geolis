//! The raw offset stage: every segment of a polyline offset on its own,
//! then joined at the corners.
//!
//! This is the first step of [`PlineOffset2D`](super::PlineOffset2D) and a
//! public kernel operation in its own right — see [`RawOffset`] for the
//! contract, and for what "raw" costs the caller.

use std::f64::consts::{PI, TAU};

use crate::error::{OperationError, Result};
use crate::geometry::pline::{Pline, PlineVertex};
use crate::math::arc_2d::{arc_from_bulge, arc_tangent_at, bulge_from_arc, offset_arc_segment};
use crate::math::intersect_2d::line_line_intersect_2d;
use crate::math::polygon_2d::{left_normal, segment_direction};
use crate::math::{Point3, TOLERANCE};

/// Threshold for flat cap: `cos(angle) < this` → near-180° reversal.
const FLAT_CAP_COS: f64 = -0.98;

/// The carrier curve an offset segment lies on — the infinite line /
/// full circle used for EXACT corner joins. Intersecting carriers (not
/// tangent-line approximations) keeps every joined arc endpoint ON its
/// offset circle, so the re-derived bulge encodes the exact concentric
/// offset arc.
#[derive(Clone, Copy)]
enum Carrier {
    Line,
    Circle {
        cx: f64,
        cy: f64,
        r: f64,
        ccw: bool,
        /// Signed sweep of the RAW (un-joined) offset arc, inherited from
        /// the source arc. Kept so a joined arc can be measured against
        /// where it started rather than guessed at from its endpoints,
        /// which cannot distinguish a sweep from its explement.
        sweep: f64,
    },
}

/// An offset segment with endpoints, carrier, and tangent directions.
struct OffsetSeg {
    start: (f64, f64),
    end: (f64, f64),
    carrier: Carrier,
    /// Unit tangent direction at the start of the segment.
    start_dir: (f64, f64),
    /// Unit tangent direction at the end of the segment.
    end_dir: (f64, f64),
}

/// A resolved corner between two consecutive offset segments.
enum Join {
    /// Single exact corner point shared by both segments.
    Miter((f64, f64)),
    /// Two points — the previous segment's own end, then the next
    /// segment's own start — connected by a straight bevel span
    /// (flat cap, miter-limit bevel, or disjoint-carrier fallback).
    Bevel((f64, f64), (f64, f64)),
}

impl Join {
    /// The point the PREVIOUS segment ends at.
    fn prev_end(&self) -> (f64, f64) {
        match self {
            Self::Miter(p) | Self::Bevel(p, _) => *p,
        }
    }

    /// The point the NEXT segment starts at.
    fn next_start(&self) -> (f64, f64) {
        match self {
            Self::Miter(p) | Self::Bevel(_, p) => *p,
        }
    }
}

/// One source segment's image under a raw offset.
///
/// Every segment of the input has exactly one of these, in input order, so
/// "how long did segment 3 end up" is answerable — the question a polyline
/// offset cannot answer once its slice-and-filter pass has merged and split
/// the output.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RawOffsetSegment {
    /// Index of the source segment: the one running from
    /// `pline.vertices[source]` to the next vertex.
    ///
    /// The raw stage prunes nothing, so today `segments[i].source == i`
    /// unconditionally. It is carried explicitly so a caller stays
    /// source-addressable if a later stage ever drops or splits a
    /// segment — the lineage question this type exists to answer.
    pub source: usize,
    /// Offset start, after the corner join with the previous segment.
    pub start: (f64, f64),
    /// Offset end, after the corner join with the next segment.
    pub end: (f64, f64),
    /// Bulge re-derived from `start` and `end` about the exact offset
    /// circle, so a joined arc stays concentric with its source. `0` for a
    /// straight segment.
    ///
    /// Meaningless when `signed_length` is negative: the endpoints then
    /// run backwards, so the bulge names the EXPLEMENT arc rather than the
    /// segment. Check `signed_length` first.
    pub bulge: f64,
    /// Length along the segment's OWN curve — the arc length for a bulged
    /// segment, which is strictly greater than the chord — signed by the
    /// source traversal direction.
    ///
    /// **Negative means the segment was consumed**: at a corner sharp
    /// enough, or on a run short enough, the two joins cross each other and
    /// the segment's image runs backwards. The raw stage reports that
    /// faithfully rather than hiding it, because "this face no longer
    /// exists" is exactly what a caller measuring faces needs to know, and
    /// exactly what the self-intersection pass downstream is there to
    /// resolve for a caller that wants a clean polyline instead.
    pub signed_length: f64,
    /// Whether the corner joining this segment to the previous one
    /// **bevelled** instead of mitering — the two carriers did not meet, or
    /// met further out than the miter limit allows, so the segments do not
    /// share an endpoint and a straight run bridges the previous segment's
    /// `end` to this one's `start`.
    ///
    /// Always `false` for the first segment of an open polyline: its start
    /// is a free end, not a corner.
    pub bevelled_entry: bool,
}

/// The raw offset of a polyline: every segment offset on its own, corners
/// joined, and **nothing cleaned up**.
///
/// This is the first stage of [`PlineOffset2D`](super::PlineOffset2D),
/// exposed because the stage is useful on its own to a caller that needs
/// per-source-segment lineage — which the finished offset destroys.
///
/// # What the stage guarantees
///
/// - A line segment is shifted along its left normal by `distance`
///   (positive is the left of travel).
/// - An arc segment stays **concentric** with its source, at radius
///   `r ± distance`; it is never approximated by a polyline or a
///   tangent-shifted arc.
/// - Corners are joined at the exact intersection of the segments' CARRIER
///   curves — the infinite line and the full offset circle, not the bounded
///   segments — so a joined arc endpoint lies on its offset circle to full
///   precision, and each segment's [`bulge`](RawOffsetSegment::bulge) is
///   re-derived from its FINAL endpoints. A corner whose carriers miss, or
///   whose meeting point runs past the miter limit, bevels instead of
///   spiking.
/// - Segments come back in source order, one per source segment, with
///   [`source`](RawOffsetSegment::source) naming the segment each came
///   from.
///
/// # What "raw" means — the caller owns these
///
/// - **No self-intersection cleanup.** Where the offset crosses itself the
///   overlapping segments are all still here, and a segment the corners ate
///   from both ends comes back with a negative
///   [`signed_length`](RawOffsetSegment::signed_length). Slicing the
///   crossings out and discarding what is closer to the source than
///   `|distance|` is the job of the stages after this one.
/// - **No validity check.** An offset that swallowed its own ring is
///   returned, not rejected.
///
/// A caller that wants a *clean* offset polyline wants
/// [`PlineOffset2D`](super::PlineOffset2D). A caller that wants to know
/// what became of each source segment wants this.
#[derive(Debug, Clone, PartialEq)]
pub struct RawOffset {
    /// One entry per source segment, in source order.
    pub segments: Vec<RawOffsetSegment>,
    /// Whether the source polyline was closed. A closed raw offset joins
    /// its last segment back to its first.
    pub closed: bool,
}

impl RawOffset {
    /// How far a corner's miter may reach from the source vertex, as a
    /// multiple of `|distance|`, before the corner bevels instead.
    ///
    /// Public because it is part of the contract, not a tuning knob: a
    /// caller that measures the resulting faces has to know that a corner
    /// sharper than this reads as two blunt ends rather than one spike, and
    /// a caller that reproduces the geometry alongside a band built from
    /// the same stage must agree with it exactly.
    pub const MITER_LIMIT: f64 = 4.0;

    /// Offsets every segment of `pline` by `distance` and joins the
    /// corners. See the type docs for the contract.
    ///
    /// # Errors
    ///
    /// Returns `OperationError::Failed` when the polyline has no segments
    /// or an arc collapses under the offset (its offset radius reaches its
    /// own centre), and `OperationError::InvalidInput` for a zero-length
    /// segment, which has no direction to offset along.
    pub fn build(pline: &Pline, distance: f64) -> Result<Self> {
        let seg_count = pline.segment_count();
        if seg_count == 0 {
            return Err(OperationError::Failed("no segments to offset".to_owned()).into());
        }

        // Phase A: offset each segment on its own.
        let offset_segs = offset_segments(pline, distance)?;

        // Phase B: resolve every corner, then measure each segment between
        // its FINAL endpoints (a joined arc endpoint moves along the offset
        // circle, so both the bulge and the length must encode the trimmed
        // or extended sweep, not the source one).
        //
        // A closed polyline has a corner at every vertex; an open one has
        // them at its interior vertices only, and its two free ends keep
        // their segments' own endpoints.
        let interior = usize::from(!pline.closed);
        let joins: Vec<Option<Join>> = (0..seg_count)
            .map(|i| {
                (i >= interior).then(|| {
                    let prev = (i + seg_count - 1) % seg_count;
                    corner_join(
                        &offset_segs[prev],
                        &offset_segs[i],
                        pline.vertices[i].x,
                        pline.vertices[i].y,
                        distance,
                    )
                })
            })
            .collect();

        let segments = (0..seg_count)
            .map(|i| {
                let seg = &offset_segs[i];
                let entry = joins[i].as_ref();
                let start = entry.map_or(seg.start, Join::next_start);
                // The corner CLOSING this segment is the one OPENING the
                // next. On an open polyline the wrap lands on the absent
                // first corner, which is the free end this segment keeps.
                let end = joins[(i + 1) % seg_count]
                    .as_ref()
                    .map_or(seg.end, Join::prev_end);
                RawOffsetSegment {
                    source: i,
                    start,
                    end,
                    bulge: seg_bulge(seg, start, end),
                    signed_length: seg_signed_length(seg, start, end),
                    bevelled_entry: matches!(entry, Some(Join::Bevel(..))),
                }
            })
            .collect();

        Ok(Self {
            segments,
            closed: pline.closed,
        })
    }

    /// Assembles the raw offset back into a polyline.
    ///
    /// A bevelled corner contributes the extra straight vertex that bridges
    /// the two segments it failed to join; an open polyline is terminated
    /// with its last segment's end.
    #[must_use]
    pub fn into_pline(self) -> Pline {
        let mut vertices = Vec::with_capacity(self.segments.len() * 2 + 1);
        for (i, seg) in self.segments.iter().enumerate() {
            if seg.bevelled_entry {
                let previous = &self.segments[(i + self.segments.len() - 1) % self.segments.len()];
                vertices.push(PlineVertex::line(previous.end.0, previous.end.1));
            }
            vertices.push(PlineVertex::new(seg.start.0, seg.start.1, seg.bulge));
        }
        if !self.closed {
            if let Some(last) = self.segments.last() {
                vertices.push(PlineVertex::line(last.end.0, last.end.1));
            }
        }
        Pline {
            vertices,
            closed: self.closed,
        }
    }
}

/// Offsets every segment of `pline` individually: lines shift along
/// their left normal, arcs stay concentric with the source arc at the
/// left-offset radius (`offset_arc_segment`).
///
/// # Errors
///
/// Returns `OperationError::InvalidInput` for zero-length line segments
/// or `OperationError::Failed` when an arc collapses under the offset.
fn offset_segments(pline: &Pline, distance: f64) -> Result<Vec<OffsetSeg>> {
    let n = pline.vertices.len();
    let seg_count = pline.segment_count();
    let mut offset_segs: Vec<OffsetSeg> = Vec::with_capacity(seg_count);

    for i in 0..seg_count {
        let v0 = &pline.vertices[i];
        let v1 = &pline.vertices[(i + 1) % n];

        if v0.bulge.abs() < 1e-12 {
            // Line segment: parallel offset.
            let p0 = Point3::new(v0.x, v0.y, 0.0);
            let p1 = Point3::new(v1.x, v1.y, 0.0);
            let dir = segment_direction(&p0, &p1)?;
            let normal = left_normal(dir);

            let start = (v0.x + normal.x * distance, v0.y + normal.y * distance);
            let end = (v1.x + normal.x * distance, v1.y + normal.y * distance);
            let d = (dir.x, dir.y);

            offset_segs.push(OffsetSeg {
                start,
                end,
                carrier: Carrier::Line,
                start_dir: d,
                end_dir: d,
            });
        } else {
            // Arc segment: change radius, preserve sweep — concentric
            // with the source arc.
            let seg = offset_arc_segment(v0.x, v0.y, v1.x, v1.y, v0.bulge, distance).ok_or_else(
                || OperationError::Failed("arc segment collapsed during offset".to_owned()),
            )?;

            let (ox0, oy0, ox1, oy1, ob) = seg;
            let (cx, cy, _, sa, sw) = arc_from_bulge(ox0, oy0, ox1, oy1, ob);
            let r = ((ox0 - cx).powi(2) + (oy0 - cy).powi(2)).sqrt();
            let sd = arc_tangent_at(sa, sw, 0.0);
            let ed = arc_tangent_at(sa, sw, 1.0);

            offset_segs.push(OffsetSeg {
                start: (ox0, oy0),
                end: (ox1, oy1),
                carrier: Carrier::Circle {
                    cx,
                    cy,
                    r,
                    ccw: ob > 0.0,
                    sweep: sw,
                },
                start_dir: sd,
                end_dir: ed,
            });
        }
    }
    Ok(offset_segs)
}

/// Bulge of one offset segment between its FINAL endpoints: `0` on a
/// line carrier; on a circle carrier, re-derived about the exact offset
/// circle in the source arc's winding.
fn seg_bulge(seg: &OffsetSeg, start: (f64, f64), end: (f64, f64)) -> f64 {
    match seg.carrier {
        Carrier::Line => 0.0,
        Carrier::Circle { cx, cy, ccw, .. } => {
            bulge_from_arc(start.0, start.1, end.0, end.1, cx, cy, ccw)
        }
    }
}

/// Length of one offset segment along its own carrier, between its FINAL
/// endpoints, signed by the source traversal direction.
///
/// Measuring along the carrier rather than as a chord is what makes a
/// consumed segment readable: when the joins at both ends crossed each
/// other the segment's image runs backwards, which a distance cannot
/// express and a signed extent can.
fn seg_signed_length(seg: &OffsetSeg, start: (f64, f64), end: (f64, f64)) -> f64 {
    match seg.carrier {
        Carrier::Line => (end.0 - start.0) * seg.start_dir.0 + (end.1 - start.1) * seg.start_dir.1,
        Carrier::Circle {
            cx,
            cy,
            r,
            ccw,
            sweep,
        } => {
            // Each join moved an endpoint ALONG the offset circle. Measure
            // that movement against where the raw arc started and ended,
            // rather than re-deriving a sweep from the final endpoints —
            // two points on a circle name two arcs, and only the raw sweep
            // says which of them this segment is.
            let angle = |p: (f64, f64)| (p.1 - cy).atan2(p.0 - cx);
            let raw_start = angle(seg.start);
            let at_start = traversal_delta(raw_start, angle(start), ccw);
            let at_end = traversal_delta(raw_start + sweep, angle(end), ccw);
            r * (sweep.abs() - at_start + at_end)
        }
    }
}

/// Angular displacement from `from` to `to` measured in the arc's traversal
/// direction and wrapped into `(−π, π]`.
///
/// A corner join moves an endpoint by the miter, never by half a turn, so
/// the short way round is always the intended one — which is what stops a
/// slightly extended tiny arc from reading as a nearly complete circle.
fn traversal_delta(from: f64, to: f64, ccw: bool) -> f64 {
    let raw = if ccw { to - from } else { from - to };
    let wrapped = raw.rem_euclid(TAU);
    if wrapped > PI {
        wrapped - TAU
    } else {
        wrapped
    }
}

/// Resolves the corner between two consecutive offset segments at the
/// original vertex `(orig_x, orig_y)`.
///
/// Handles three cases:
/// 1. Near-antiparallel (>~169°): flat cap (bevel).
/// 2. Miter too long: bevel.
/// 3. Normal corner: single exact corner point — line × line miter for
///    two straight segments (unchanged legacy math), carrier
///    intersection (line × circle / circle × circle) when an arc is
///    involved so the joined point lies exactly on the offset circle.
fn corner_join(
    seg_prev: &OffsetSeg,
    seg_next: &OffsetSeg,
    orig_x: f64,
    orig_y: f64,
    distance: f64,
) -> Join {
    let dir_prev = &seg_prev.end_dir;
    let dir_next = &seg_next.start_dir;
    let cos_angle = dir_prev.0 * dir_next.0 + dir_prev.1 * dir_next.1;

    if cos_angle < FLAT_CAP_COS {
        // Near-antiparallel: flat cap.
        return Join::Bevel(seg_prev.end, seg_next.start);
    }

    let corner = match (&seg_prev.carrier, &seg_next.carrier) {
        (Carrier::Line, Carrier::Line) => {
            // Legacy line-line miter via tangent intersection (exact for
            // straight carriers).
            let p_prev = Point3::new(seg_prev.end.0, seg_prev.end.1, 0.0);
            let d_prev = crate::math::Vector3::new(dir_prev.0, dir_prev.1, 0.0);
            let p_next = Point3::new(seg_next.start.0, seg_next.start.1, 0.0);
            let d_next = crate::math::Vector3::new(dir_next.0, dir_next.1, 0.0);
            let Some((t, _)) = line_line_intersect_2d(&p_prev, &d_prev, &p_next, &d_next) else {
                // Parallel: use offset of the original corner point
                // (no miter-limit check — matches the legacy path).
                let fallback_normal = left_normal(
                    crate::math::Vector3::new(d_prev.x, d_prev.y, 0.0)
                        .try_normalize(TOLERANCE)
                        .unwrap_or(crate::math::Vector3::new(1.0, 0.0, 0.0)),
                );
                return Join::Miter((
                    orig_x + fallback_normal.x * distance,
                    orig_y + fallback_normal.y * distance,
                ));
            };
            (p_prev.x + d_prev.x * t, p_prev.y + d_prev.y * t)
        }
        pair => {
            // At least one arc: intersect the exact carriers and keep
            // the root nearest the two segment endpoints. Disjoint
            // carriers (arc curving away) get a straight bevel.
            let candidates = match pair {
                (Carrier::Line, Carrier::Circle { cx, cy, r, .. }) => {
                    line_circle_intersections(seg_prev.start, *dir_prev, (*cx, *cy), *r)
                }
                (Carrier::Circle { cx, cy, r, .. }, Carrier::Line) => {
                    line_circle_intersections(seg_next.start, *dir_next, (*cx, *cy), *r)
                }
                (
                    Carrier::Circle { cx, cy, r, .. },
                    Carrier::Circle {
                        cx: bx,
                        cy: by,
                        r: br,
                        ..
                    },
                ) => {
                    // Two arcs split from the SAME carrier circle (a
                    // tangent-continuous junction): the carriers are
                    // coincident, so intersecting them yields nothing
                    // (concentric ⇒ empty) and would force a zero-length
                    // bevel between the coincident endpoints. The join is
                    // the shared offset endpoint itself — the legacy
                    // single-point join.
                    if (cx - bx).abs() < TOLERANCE
                        && (cy - by).abs() < TOLERANCE
                        && (r - br).abs() < TOLERANCE
                    {
                        vec![(
                            0.5 * (seg_prev.end.0 + seg_next.start.0),
                            0.5 * (seg_prev.end.1 + seg_next.start.1),
                        )]
                    } else {
                        circle_circle_intersections((*cx, *cy), *r, (*bx, *by), *br)
                    }
                }
                (Carrier::Line, Carrier::Line) => unreachable!("handled above"),
            };
            let nearest = candidates.into_iter().min_by(|a, b| {
                let score = |p: &(f64, f64)| {
                    (p.0 - seg_prev.end.0).powi(2)
                        + (p.1 - seg_prev.end.1).powi(2)
                        + (p.0 - seg_next.start.0).powi(2)
                        + (p.1 - seg_next.start.1).powi(2)
                };
                score(a).total_cmp(&score(b))
            });
            match nearest {
                Some(p) => p,
                None => return Join::Bevel(seg_prev.end, seg_next.start),
            }
        }
    };

    let dx = corner.0 - orig_x;
    let dy = corner.1 - orig_y;
    let miter_dist_sq = dx * dx + dy * dy;
    let limit = RawOffset::MITER_LIMIT * distance.abs();
    if miter_dist_sq > limit * limit {
        // Miter too long: bevel.
        Join::Bevel(seg_prev.end, seg_next.start)
    } else {
        Join::Miter(corner)
    }
}

/// Intersections of an INFINITE line (point + unit direction) with a
/// FULL circle. Tangential contact yields the single tangent point.
fn line_circle_intersections(
    p: (f64, f64),
    dir: (f64, f64),
    center: (f64, f64),
    radius: f64,
) -> Vec<(f64, f64)> {
    let fx = p.0 - center.0;
    let fy = p.1 - center.1;
    // dir is unit length: t² + 2bt + c = 0.
    let b = fx * dir.0 + fy * dir.1;
    let c = fx * fx + fy * fy - radius * radius;
    let disc = b * b - c;
    if disc < -TOLERANCE {
        return Vec::new();
    }
    if disc <= TOLERANCE {
        let t = -b;
        return vec![(p.0 + dir.0 * t, p.1 + dir.1 * t)];
    }
    let sqrt_disc = disc.sqrt();
    [-b - sqrt_disc, -b + sqrt_disc]
        .iter()
        .map(|t| (p.0 + dir.0 * t, p.1 + dir.1 * t))
        .collect()
}

/// Intersections of two FULL circles. Concentric or disjoint circles
/// yield none; tangential contact yields the single touch point.
fn circle_circle_intersections(a: (f64, f64), ra: f64, b: (f64, f64), rb: f64) -> Vec<(f64, f64)> {
    let dx = b.0 - a.0;
    let dy = b.1 - a.1;
    let dist = (dx * dx + dy * dy).sqrt();
    if dist < TOLERANCE || dist > ra + rb + TOLERANCE || dist < (ra - rb).abs() - TOLERANCE {
        return Vec::new();
    }
    let along = (ra * ra - rb * rb + dist * dist) / (2.0 * dist);
    let h_sq = ra * ra - along * along;
    let h = h_sq.max(0.0).sqrt();
    let mx = a.0 + along * dx / dist;
    let my = a.1 + along * dy / dist;
    if h < TOLERANCE {
        return vec![(mx, my)];
    }
    let px = -dy / dist;
    let py = dx / dist;
    vec![(mx + h * px, my + h * py), (mx - h * px, my - h * py)]
}

#[cfg(test)]
#[allow(clippy::expect_used, clippy::unwrap_used)]
mod tests {
    use std::f64::consts::FRAC_PI_2;

    use super::*;

    /// Tighter than a micrometre on metre-scale input, so a wrong join
    /// cannot hide behind the tolerance.
    const EPS: f64 = 1e-9;

    fn open(vertices: Vec<PlineVertex>) -> Pline {
        Pline {
            vertices,
            closed: false,
        }
    }

    /// Bulge of a CCW quarter turn.
    fn quarter() -> f64 {
        (PI / 8.0).tan()
    }

    /// The offset segment that came from source segment `source`.
    fn from_source(raw: &RawOffset, source: usize) -> RawOffsetSegment {
        *raw.segments
            .iter()
            .find(|s| s.source == source)
            .expect("every source segment has an image")
    }

    // ── lineage ──────────────────────────────────────────────────────

    /// The property the whole API exists for: one image per source
    /// segment, in source order, whatever the corners did.
    #[test]
    fn every_source_segment_gets_exactly_one_image_in_order() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(6.0, 0.0),
                PlineVertex::line(6.0, 4.0),
                PlineVertex::line(0.0, 4.0),
            ],
            closed: true,
        };
        let raw = RawOffset::build(&pline, 0.5).unwrap();
        assert!(raw.closed);
        assert_eq!(
            raw.segments.iter().map(|s| s.source).collect::<Vec<_>>(),
            vec![0, 1, 2, 3],
        );
    }

    // ── a straight segment ───────────────────────────────────────────

    /// A lone straight segment shifts along its left normal and keeps its
    /// length: an open polyline's free ends are not corners.
    #[test]
    fn a_straight_segment_shifts_along_its_left_normal() {
        let pline = open(vec![
            PlineVertex::line(1.0, 2.0),
            PlineVertex::line(5.0, 2.0),
        ]);
        for (distance, expected_y) in [(0.3, 2.3), (-0.3, 1.7)] {
            let seg = from_source(&RawOffset::build(&pline, distance).unwrap(), 0);
            assert!((seg.start.0 - 1.0).abs() < EPS && (seg.start.1 - expected_y).abs() < EPS);
            assert!((seg.end.0 - 5.0).abs() < EPS && (seg.end.1 - expected_y).abs() < EPS);
            assert!(seg.bulge.abs() < EPS, "a straight segment stays straight");
            assert!((seg.signed_length - 4.0).abs() < EPS);
            assert!(!seg.bevelled_entry);
        }
    }

    /// A right-angle corner takes a half-distance off the inner face and
    /// adds one to the outer — the miter, and the reason a source segment's
    /// image is never as long as the segment itself.
    #[test]
    fn a_corner_miters_the_two_faces_by_the_offset_distance() {
        let pline = open(vec![
            PlineVertex::line(0.0, 0.0),
            PlineVertex::line(4.0, 0.0),
            PlineVertex::line(4.0, 3.0),
        ]);
        let inner = from_source(&RawOffset::build(&pline, 0.25).unwrap(), 0);
        let outer = from_source(&RawOffset::build(&pline, -0.25).unwrap(), 0);
        assert!((inner.signed_length - 3.75).abs() < EPS);
        assert!((outer.signed_length - 4.25).abs() < EPS);
    }

    // ── an arc stays concentric ──────────────────────────────────────

    /// An arc's image is an arc on the exact concentric circle at
    /// `r ∓ distance` — never a tangent shift, never a polyline.
    #[test]
    fn an_arc_offsets_to_its_exact_concentric_circle() {
        let radius = 2.0;
        let pline = open(vec![
            PlineVertex::new(radius, 0.0, quarter()),
            PlineVertex::line(0.0, radius),
        ]);
        // The centre is on the left of a CCW arc's travel, so a positive
        // (leftward) offset shrinks the radius.
        for (distance, offset_radius) in [(0.4, radius - 0.4), (-0.4, radius + 0.4)] {
            let seg = from_source(&RawOffset::build(&pline, distance).unwrap(), 0);
            for end in [seg.start, seg.end] {
                assert!(
                    (end.0.hypot(end.1) - offset_radius).abs() < EPS,
                    "endpoint must sit on the offset circle of radius {offset_radius}"
                );
            }
            assert!(
                (seg.signed_length - FRAC_PI_2 * offset_radius).abs() < EPS,
                "an arc measures its arc length, not its chord"
            );
            let chord = (seg.end.0 - seg.start.0).hypot(seg.end.1 - seg.start.1);
            assert!(chord < seg.signed_length);
        }
    }

    /// The bulge is re-derived from the FINAL endpoints about the offset
    /// circle. A corner that MOVED the arc's start therefore leaves a bulge
    /// encoding the new sweep, and `arc_from_bulge` on it recovers the
    /// exact source centre and the offset radius — the property a
    /// bisector-miter approximation misses by `O(distance² / r)`.
    ///
    /// The corner is deliberately not tangent-continuous: a straight run
    /// arriving along the arc's own tangent joins at the arc's raw start
    /// and proves nothing.
    #[test]
    fn a_joined_arcs_bulge_is_re_derived_about_the_offset_circle() {
        let radius = 2.0;
        // A straight run into (2, 0) at 45° to the arc's tangent there,
        // then the CCW quarter arc to (0, 2).
        let pline = open(vec![
            PlineVertex::line(4.0, -2.0),
            PlineVertex::new(2.0, 0.0, quarter()),
            PlineVertex::line(0.0, radius),
        ]);
        let distance = 0.3;
        let raw = RawOffset::build(&pline, distance).unwrap();
        let straight = from_source(&raw, 0);
        let arc = from_source(&raw, 1);
        let offset_radius = radius - distance;

        // The join is ONE point, shared exactly, and it is on the circle.
        assert!((straight.end.0 - arc.start.0).abs() < 1e-12);
        assert!((straight.end.1 - arc.start.1).abs() < 1e-12);
        assert!(!arc.bevelled_entry);
        assert!(
            (arc.start.0.hypot(arc.start.1) - offset_radius).abs() < EPS,
            "the joint sits on the exact offset circle"
        );
        // ...and it is NOT where the un-joined offset arc started, so the
        // re-derivation is doing real work.
        assert!(
            (arc.start.1).abs() > 1e-6,
            "the corner must have moved the arc's start off (r, 0)"
        );

        let (cx, cy, r, _, sweep) =
            arc_from_bulge(arc.start.0, arc.start.1, arc.end.0, arc.end.1, arc.bulge);
        assert!(cx.abs() < EPS && cy.abs() < EPS, "centre is preserved");
        assert!((r - offset_radius).abs() < EPS, "radius is r − distance");
        assert!(
            (sweep - FRAC_PI_2).abs() > 1e-6,
            "the join re-swept the arc: {sweep}"
        );
        assert!(
            (arc.signed_length - r * sweep).abs() < EPS,
            "the measured length is the re-derived sweep's arc length"
        );
    }

    // ── what "raw" means ─────────────────────────────────────────────

    /// A segment the joins ate from both ends reports a NEGATIVE length
    /// rather than a mirrored positive one. Nothing prunes it — that is
    /// the caller's business.
    #[test]
    fn a_segment_consumed_by_its_joins_reports_a_negative_length() {
        // A hairpin: 0.2 of crossbar between two right turns, so both
        // corners eat 0.5 off the crossbar's right face.
        let pline = open(vec![
            PlineVertex::line(0.0, 0.0),
            PlineVertex::line(0.0, 3.0),
            PlineVertex::line(0.2, 3.0),
            PlineVertex::line(0.2, 0.0),
        ]);
        let raw = RawOffset::build(&pline, -0.5).unwrap();
        assert_eq!(raw.segments.len(), 3, "nothing is pruned");
        assert!(
            from_source(&raw, 1).signed_length < 0.0,
            "the crossbar's short face ran backwards"
        );
        // ...and the same crossbar's other face gains both miters.
        let wide = RawOffset::build(&pline, 0.5).unwrap();
        assert!((from_source(&wide, 1).signed_length - 1.2).abs() < EPS);
    }

    /// A near-180° reversal has no usable corner: it bevels, each segment
    /// keeps its own raw end, and the bevel is reported.
    #[test]
    fn a_reversal_bevels_and_says_so() {
        let pline = open(vec![
            PlineVertex::line(0.0, 0.0),
            PlineVertex::line(4.0, 0.0),
            PlineVertex::line(0.5, 0.0),
        ]);
        let raw = RawOffset::build(&pline, 0.15).unwrap();
        assert!(
            !from_source(&raw, 0).bevelled_entry,
            "a free end is not a corner"
        );
        assert!(from_source(&raw, 1).bevelled_entry);
        assert!((from_source(&raw, 0).signed_length - 4.0).abs() < EPS);
        assert!((from_source(&raw, 1).signed_length - 3.5).abs() < EPS);
    }

    /// An arc whose offset would reach its own centre has no image at all,
    /// and the whole stage says so rather than inventing a mirrored arc.
    #[test]
    fn an_arc_collapsing_under_the_offset_is_an_error() {
        let pline = open(vec![
            PlineVertex::new(0.5, 0.0, quarter()),
            PlineVertex::line(0.0, 0.5),
        ]);
        assert!(RawOffset::build(&pline, 0.8).is_err());
        assert!(RawOffset::build(&pline, -0.8).is_ok(), "outward is fine");
    }

    #[test]
    fn a_polyline_with_no_segments_is_an_error() {
        assert!(RawOffset::build(&open(vec![PlineVertex::line(0.0, 0.0)]), 1.0).is_err());
    }

    // ── the assembled polyline ───────────────────────────────────────

    /// `into_pline` is the shape the offset pipeline consumes: one vertex
    /// per mitered segment, plus the bridging vertex a bevel needs, plus
    /// the terminating vertex of an open run.
    #[test]
    fn into_pline_emits_one_vertex_per_mitered_segment() {
        let square = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(10.0, 0.0),
                PlineVertex::line(10.0, 10.0),
                PlineVertex::line(0.0, 10.0),
            ],
            closed: true,
        };
        let pline = RawOffset::build(&square, 1.0).unwrap().into_pline();
        assert!(pline.closed);
        assert_eq!(pline.vertices.len(), 4);
        for v in &pline.vertices {
            assert!(v.bulge.abs() < EPS);
            assert!((v.x - 1.0).abs() < EPS || (v.x - 9.0).abs() < EPS);
            assert!((v.y - 1.0).abs() < EPS || (v.y - 9.0).abs() < EPS);
        }
    }

    #[test]
    fn into_pline_bridges_a_bevel_and_terminates_an_open_run() {
        let pline = open(vec![
            PlineVertex::line(0.0, 0.0),
            PlineVertex::line(4.0, 0.0),
            PlineVertex::line(0.5, 0.0),
        ]);
        let raw = RawOffset::build(&pline, 0.15).unwrap();
        let bevels = raw.segments.iter().filter(|s| s.bevelled_entry).count();
        let assembled = raw.into_pline();
        // One vertex per segment, one more for the bevel bridge, one more
        // to terminate the open run.
        assert_eq!(assembled.vertices.len(), 2 + bevels + 1);
        assert!(!assembled.closed);
    }
}
