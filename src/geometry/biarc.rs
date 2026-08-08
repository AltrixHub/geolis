//! G1 Hermite interpolation with circular arcs (biarcs).
//!
//! A biarc joins two points with prescribed tangent *directions* using
//! at most two circular arcs that meet with a common tangent. Results
//! are expressed in this crate's bulge form (see
//! [`crate::math::arc_2d`]), so they drop straight into a
//! [`Pline`](crate::geometry::Pline).
//!
//! # Construction
//!
//! The construction follows the classical joint-locus property (Bolton,
//! "Biarc curves", *Computer-Aided Design* 7(2), 1975): for fixed
//! endpoints and tangents the admissible joint points form a circle
//! through both endpoints.
//!
//! With this crate's sign convention an arc's start tangent lies
//! `sweep / 2` clockwise of its chord and its end tangent `sweep / 2`
//! counter-clockwise of it (see
//! [`bulge_from_chord_tangent`](crate::math::arc_2d::bulge_from_chord_tangent)).
//! Writing `c1`, `c2` for the directions of the two sub-chords
//! `p0 -> joint` and `joint -> p1`, tangent continuity at the joint
//! reads `2*c1 - t0 = 2*c2 - t1`, i.e.
//!
//! ```text
//! c1 - c2 = (t0 - t1) / 2   (mod pi)
//! ```
//!
//! The angle subtended at the joint is therefore constant, so by the
//! inscribed-angle theorem the joint locus is a circle through `p0` and
//! `p1` — equivalently, the arc from `p0` to `p1` whose sweep is the
//! total turn `t1 - t0`. That arc makes the same angle with `t0` at
//! `p0` as it does with `t1` at `p1` (the tangent-bisector property),
//! and `joint_ratio` selects a point along it, which is what makes the
//! biarc family one-parameter.
//!
//! The joint is evaluated in closed form from the chord rather than
//! from a center/radius, because the locus circle degenerates to the
//! chord line as the total turn goes to zero — the most common input
//! for a pen tool — and a center-based evaluation loses all precision
//! there.

use std::f64::consts::{PI, TAU};

use crate::math::arc_2d::bulge_from_chord_tangent;

/// Degenerate-length threshold, matching [`crate::math::arc_2d`].
const EPS: f64 = 1e-12;

/// Angular tolerance for direction comparisons, in radians.
const ANGULAR_EPS: f64 = 1e-9;

/// Keeps the joint clear of both endpoints so that neither sub-chord
/// collapses, whatever `joint_ratio` the caller passes.
///
/// [`biarc_from_hermite`] clamps into `[JOINT_RATIO_MARGIN, 1 -
/// JOINT_RATIO_MARGIN]`, so that closed interval is the family a caller
/// can actually address — and the range
/// [`biarc_joint_ratio`] answers within.
pub const JOINT_RATIO_MARGIN: f64 = 1e-6;

/// Relative distance (against the chord length) within which a
/// re-solved joint counts as the requested one.
const JOINT_MATCH_TOL: f64 = 1e-9;

/// The chord frame both the forward and the inverse solve read: the
/// chord `p0 -> p1` plus the offset of each tangent from it.
#[derive(Debug, Clone, Copy)]
struct ChordFrame {
    chord_len: f64,
    chord_angle: f64,
    start_offset: f64,
    end_offset: f64,
}

impl ChordFrame {
    /// `None` for a degenerate element: a zero-length chord, a
    /// zero-length tangent, or any non-finite input.
    #[allow(clippy::too_many_arguments)]
    fn new(
        x0: f64,
        y0: f64,
        tx0: f64,
        ty0: f64,
        x1: f64,
        y1: f64,
        tx1: f64,
        ty1: f64,
    ) -> Option<Self> {
        if ![x0, y0, tx0, ty0, x1, y1, tx1, ty1]
            .iter()
            .all(|v| v.is_finite())
        {
            return None;
        }
        let dx = x1 - x0;
        let dy = y1 - y0;
        let chord_len = (dx * dx + dy * dy).sqrt();
        if chord_len < EPS
            || (tx0 * tx0 + ty0 * ty0).sqrt() < EPS
            || (tx1 * tx1 + ty1 * ty1).sqrt() < EPS
        {
            return None;
        }
        let chord_angle = dy.atan2(dx);
        Some(Self {
            chord_len,
            chord_angle,
            start_offset: wrap_signed(ty0.atan2(tx0) - chord_angle),
            end_offset: wrap_signed(ty1.atan2(tx1) - chord_angle),
        })
    }

    /// Both tangents already point along the chord, so the element is a
    /// straight segment.
    fn is_straight(&self) -> bool {
        self.start_offset.abs() < ANGULAR_EPS && self.end_offset.abs() < ANGULAR_EPS
    }

    /// Half the total turn from the start tangent to the end tangent —
    /// half the sweep of the locus arc the joint lives on.
    fn half_turn(&self) -> f64 {
        wrap_signed(self.end_offset - self.start_offset) * 0.5
    }

    /// The joint at `ratio` on the locus, as an offset from `p0`.
    ///
    /// Its sub-chord from `p0` has direction `chord_angle + half *
    /// (ratio - 1)` and length `chord_len * sin(half * ratio) /
    /// sin(half)`, both read off the tangent-chord relation on the locus
    /// arc; the length ratio tends to `ratio` as that arc flattens.
    fn joint_offset(&self, ratio: f64) -> (f64, f64) {
        let half = self.half_turn();
        let sin_half = half.sin();
        let scale = if sin_half.abs() < EPS {
            ratio
        } else {
            (half * ratio).sin() / sin_half
        };
        let angle = self.chord_angle + half * (ratio - 1.0);
        (
            self.chord_len * scale * angle.cos(),
            self.chord_len * scale * angle.sin(),
        )
    }
}

/// The curve produced by [`biarc_from_hermite`].
///
/// Every variant describes the run from `p0` to `p1`; the endpoints are
/// the caller's and are not repeated here.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BiarcShape {
    /// A straight segment from `p0` to `p1` (bulge 0).
    Straight,
    /// A single arc from `p0` to `p1` — no joint is needed because one
    /// arc already meets both tangents.
    SingleArc {
        /// Bulge of the `p0` -> `p1` segment.
        bulge: f64,
    },
    /// Two arcs, `p0` -> joint -> `p1`, meeting with a common tangent.
    Biarc {
        /// X coordinate of the joint between the two arcs.
        joint_x: f64,
        /// Y coordinate of the joint between the two arcs.
        joint_y: f64,
        /// Bulge of the `p0` -> joint segment.
        bulge_first: f64,
        /// Bulge of the joint -> `p1` segment.
        bulge_second: f64,
    },
}

/// Interpolates a G1 Hermite element — two points plus a tangent
/// direction at each — with at most two circular arcs.
///
/// Both tangent directions are matched exactly and the two arcs share a
/// tangent at their joint. Tangents need not be normalized; only their
/// directions are used.
///
/// `joint_ratio` selects the member of the one-parameter biarc family:
/// it is the parameter of the joint along the locus arc (`0.5` =
/// balanced), and is clamped into the open interval so the joint never
/// lands on an endpoint. A non-finite ratio falls back to `0.5`.
///
/// # Shapes
///
/// - [`BiarcShape::Straight`] — the chord is degenerate, a tangent is
///   degenerate, an input is non-finite, or both tangents already point
///   along the chord.
/// - [`BiarcShape::SingleArc`] — the arc from `p0` with tangent `t0`
///   that reaches `p1` already ends with tangent `t1` (within
///   [`ANGULAR_EPS`]). The biarc family degenerates continuously into
///   this case: as the configuration approaches it, the two arcs
///   approach sub-arcs of the single one.
/// - [`BiarcShape::Biarc`] — otherwise. S-shaped configurations come
///   out with opposite-sign bulges, C-shaped ones with same-sign
///   bulges.
///
/// Each arc's sweep is bounded by the
/// [`MAX_TANGENT_CHORD_ANGLE`](crate::math::arc_2d::MAX_TANGENT_CHORD_ANGLE)
/// guard. Configurations that need more than that — notably tangents
/// that point back along the chord at *both* ends, where no arc pair
/// can honour them — return a clamped, finite pair rather than an
/// infinite bulge; the tangent match is given up in exchange.
///
/// # Panics
///
/// Does not panic, and never returns a non-finite bulge or joint.
#[allow(clippy::too_many_arguments)]
#[must_use]
pub fn biarc_from_hermite(
    x0: f64,
    y0: f64,
    tx0: f64,
    ty0: f64,
    x1: f64,
    y1: f64,
    tx1: f64,
    ty1: f64,
    joint_ratio: f64,
) -> BiarcShape {
    let Some(frame) = ChordFrame::new(x0, y0, tx0, ty0, x1, y1, tx1, ty1) else {
        return BiarcShape::Straight;
    };
    if frame.is_straight() {
        return BiarcShape::Straight;
    }

    // One arc may already do the job; its sweep then reaches t1 exactly
    // (possibly after a full turn, e.g. a 270-degree arc).
    let single = bulge_from_chord_tangent(x0, y0, x1, y1, tx0, ty0);
    if wrap_signed(4.0 * single.atan() - 2.0 * frame.half_turn()).abs() < ANGULAR_EPS {
        return BiarcShape::SingleArc { bulge: single };
    }

    let ratio = if joint_ratio.is_finite() {
        joint_ratio.clamp(JOINT_RATIO_MARGIN, 1.0 - JOINT_RATIO_MARGIN)
    } else {
        0.5
    };

    let (offset_x, offset_y) = frame.joint_offset(ratio);
    let joint_x = x0 + offset_x;
    let joint_y = y0 + offset_y;

    // The first arc is built forward from p0; the second is built
    // backwards from p1 (and negated) so that both prescribed tangents
    // are honoured exactly. G1 at the joint follows from the locus.
    let bulge_first = bulge_from_chord_tangent(x0, y0, joint_x, joint_y, tx0, ty0);
    let bulge_second = -bulge_from_chord_tangent(x1, y1, joint_x, joint_y, -tx1, -ty1);

    BiarcShape::Biarc {
        joint_x,
        joint_y,
        bulge_first,
        bulge_second,
    }
}

/// Recovers the `joint_ratio` that puts a biarc's joint at `(joint_x,
/// joint_y)` — the inverse of [`biarc_from_hermite`].
///
/// This is what makes an authored biarc closed under splitting: cut one
/// of its arcs, and the surviving joint is reproduced exactly by the
/// ratio this returns for the remaining Hermite element, so the curve
/// does not move when a vertex is inserted on it.
///
/// Returns `None` when no admissible ratio produces that joint:
///
/// - the Hermite data is degenerate, or does not need two arcs at all
///   (see [`BiarcShape`]) — there is no joint to place;
/// - the requested joint is not on the biarc family's joint locus for
///   that data (checked by re-solving, within [`JOINT_MATCH_TOL`] of the
///   chord length);
/// - the ratio it would need lies outside
///   `[JOINT_RATIO_MARGIN, 1 - JOINT_RATIO_MARGIN]`, the interval
///   [`biarc_from_hermite`] clamps to — i.e. the joint sits (all but) on
///   an endpoint.
///
/// # Panics
///
/// Does not panic.
#[allow(clippy::too_many_arguments)]
#[must_use]
pub fn biarc_joint_ratio(
    x0: f64,
    y0: f64,
    tx0: f64,
    ty0: f64,
    x1: f64,
    y1: f64,
    tx1: f64,
    ty1: f64,
    joint_x: f64,
    joint_y: f64,
) -> Option<f64> {
    if !joint_x.is_finite() || !joint_y.is_finite() {
        return None;
    }
    let frame = ChordFrame::new(x0, y0, tx0, ty0, x1, y1, tx1, ty1)?;

    let jdx = joint_x - x0;
    let jdy = joint_y - y0;
    if (jdx * jdx + jdy * jdy).sqrt() < EPS {
        return None;
    }

    let half = frame.half_turn();
    let ratio = if half.sin().abs() < EPS {
        // The locus has flattened onto the chord, where the joint's
        // parameter is simply how far along the chord it sits.
        (jdx * frame.chord_angle.cos() + jdy * frame.chord_angle.sin()) / frame.chord_len
    } else {
        // `ChordFrame::joint_offset` aims the joint `half * (ratio - 1)`
        // off the chord direction. That offset is under a quarter turn
        // in magnitude (`|half| <= pi / 2`, `ratio` in the unit
        // interval), so measuring it back inverts without ambiguity.
        1.0 + wrap_signed(jdy.atan2(jdx) - frame.chord_angle) / half
    };
    if !(JOINT_RATIO_MARGIN..=1.0 - JOINT_RATIO_MARGIN).contains(&ratio) {
        return None;
    }

    // Direction alone cannot tell a joint that is ON the locus from one
    // merely aimed at it, so confirm by re-solving.
    match biarc_from_hermite(x0, y0, tx0, ty0, x1, y1, tx1, ty1, ratio) {
        BiarcShape::Biarc {
            joint_x: solved_x,
            joint_y: solved_y,
            ..
        } if (solved_x - joint_x).hypot(solved_y - joint_y)
            <= JOINT_MATCH_TOL * frame.chord_len =>
        {
            Some(ratio)
        }
        _ => None,
    }
}

/// Wraps an angle into `(-pi, pi]`.
fn wrap_signed(angle: f64) -> f64 {
    let wrapped = angle % TAU;
    if wrapped > PI {
        wrapped - TAU
    } else if wrapped <= -PI {
        wrapped + TAU
    } else {
        wrapped
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::arc_2d::{arc_from_bulge, arc_point_at, arc_tangent_at};

    const TOL: f64 = 1e-9;

    /// Start and end unit tangents of a bulge segment, taken from the
    /// existing arc evaluator.
    fn segment_tangents(
        x0: f64,
        y0: f64,
        x1: f64,
        y1: f64,
        bulge: f64,
    ) -> ((f64, f64), (f64, f64)) {
        if bulge.abs() < EPS {
            let len = ((x1 - x0).powi(2) + (y1 - y0).powi(2)).sqrt();
            let dir = ((x1 - x0) / len, (y1 - y0) / len);
            return (dir, dir);
        }
        let (_, _, _, start_angle, sweep) = arc_from_bulge(x0, y0, x1, y1, bulge);
        (
            arc_tangent_at(start_angle, sweep, 0.0),
            arc_tangent_at(start_angle, sweep, 1.0),
        )
    }

    /// End point of a bulge segment, taken from the existing arc evaluator.
    fn segment_end(x0: f64, y0: f64, x1: f64, y1: f64, bulge: f64) -> (f64, f64) {
        if bulge.abs() < EPS {
            return (x1, y1);
        }
        let (cx, cy, radius, start_angle, sweep) = arc_from_bulge(x0, y0, x1, y1, bulge);
        arc_point_at(cx, cy, radius, start_angle, sweep, 1.0)
    }

    fn unit(v: (f64, f64)) -> (f64, f64) {
        let len = (v.0 * v.0 + v.1 * v.1).sqrt();
        (v.0 / len, v.1 / len)
    }

    #[track_caller]
    fn assert_direction(actual: (f64, f64), expected: (f64, f64), what: &str) {
        let a = unit(actual);
        let e = unit(expected);
        assert!(
            (a.0 - e.0).abs() < TOL && (a.1 - e.1).abs() < TOL,
            "{what}: got {a:?}, expected {e:?}"
        );
    }

    #[track_caller]
    fn assert_point(actual: (f64, f64), expected: (f64, f64), what: &str) {
        assert!(
            (actual.0 - expected.0).abs() < TOL && (actual.1 - expected.1).abs() < TOL,
            "{what}: got {actual:?}, expected {expected:?}"
        );
    }

    /// Deconstructs a `Biarc`, failing on any other shape.
    #[track_caller]
    fn expect_biarc(shape: BiarcShape) -> (f64, f64, f64, f64) {
        match shape {
            BiarcShape::Biarc {
                joint_x,
                joint_y,
                bulge_first,
                bulge_second,
            } => (joint_x, joint_y, bulge_first, bulge_second),
            other => panic!("expected a Biarc, got {other:?}"),
        }
    }

    /// An S-shaped configuration: same tangent at both ends, offset chord.
    const S_CURVE: (f64, f64, f64, f64, f64, f64, f64, f64) =
        (0.0, 0.0, 1.0, 0.0, 2.0, 2.0, 1.0, 0.0);

    /// A C-shaped configuration: the tangent turns by -70 degrees, which
    /// no single arc from `t0` can absorb.
    fn c_curve() -> (f64, f64, f64, f64, f64, f64, f64, f64) {
        let t0 = 45.0_f64.to_radians();
        let t1 = -25.0_f64.to_radians();
        (0.0, 0.0, t0.cos(), t0.sin(), 2.0, 0.0, t1.cos(), t1.sin())
    }

    fn solve(config: (f64, f64, f64, f64, f64, f64, f64, f64), ratio: f64) -> BiarcShape {
        let (x0, y0, tx0, ty0, x1, y1, tx1, ty1) = config;
        biarc_from_hermite(x0, y0, tx0, ty0, x1, y1, tx1, ty1, ratio)
    }

    #[test]
    fn straight_when_both_tangents_follow_chord() {
        let shape = biarc_from_hermite(0.0, 0.0, 3.0, 4.0, 3.0, 4.0, 0.6, 0.8, 0.5);
        assert_eq!(shape, BiarcShape::Straight);
    }

    #[test]
    fn single_arc_when_tangents_come_from_one_arc() {
        // Tangents lifted off an existing arc must come back as that arc.
        for bulge in [0.3, -0.3, 1.0, -1.0, 0.0001, 2.414_213_562_373_095] {
            let (x0, y0, x1, y1) = (-1.3, 0.7, 2.4, -0.9);
            let (t0, t1) = segment_tangents(x0, y0, x1, y1, bulge);
            let shape = biarc_from_hermite(x0, y0, t0.0, t0.1, x1, y1, t1.0, t1.1, 0.5);
            match shape {
                BiarcShape::SingleArc { bulge: recovered } => assert!(
                    (recovered - bulge).abs() < TOL,
                    "bulge={bulge} recovered={recovered}"
                ),
                other => panic!("expected SingleArc for bulge={bulge}, got {other:?}"),
            }
        }
    }

    #[test]
    fn biarc_matches_endpoint_tangents() {
        for config in [S_CURVE, c_curve()] {
            let (x0, y0, tx0, ty0, x1, y1, tx1, ty1) = config;
            for ratio in [0.25, 0.5, 0.8] {
                let (jx, jy, b1, b2) = expect_biarc(solve(config, ratio));
                let (start, _) = segment_tangents(x0, y0, jx, jy, b1);
                let (_, end) = segment_tangents(jx, jy, x1, y1, b2);
                assert_direction(start, (tx0, ty0), "start tangent");
                assert_direction(end, (tx1, ty1), "end tangent");
            }
        }
    }

    #[test]
    fn biarc_is_g1_at_joint() {
        for config in [S_CURVE, c_curve()] {
            let (x0, y0, _, _, x1, y1, _, _) = config;
            for ratio in [0.25, 0.5, 0.8] {
                let (jx, jy, b1, b2) = expect_biarc(solve(config, ratio));
                let (_, first_end) = segment_tangents(x0, y0, jx, jy, b1);
                let (second_start, _) = segment_tangents(jx, jy, x1, y1, b2);
                assert_direction(first_end, second_start, "joint tangent");
            }
        }
    }

    #[test]
    fn biarc_segments_connect_the_endpoints() {
        for config in [S_CURVE, c_curve()] {
            let (x0, y0, _, _, x1, y1, _, _) = config;
            let (jx, jy, b1, b2) = expect_biarc(solve(config, 0.35));
            assert_point(segment_end(x0, y0, jx, jy, b1), (jx, jy), "first arc end");
            assert_point(segment_end(jx, jy, x1, y1, b2), (x1, y1), "second arc end");
        }
    }

    #[test]
    fn s_curve_bulges_have_opposite_signs() {
        let (_, _, b1, b2) = expect_biarc(solve(S_CURVE, 0.5));
        assert!(b1 > 0.0 && b2 < 0.0, "b1={b1} b2={b2}");
    }

    #[test]
    fn c_curve_bulges_have_the_same_sign() {
        let (_, _, b1, b2) = expect_biarc(solve(c_curve(), 0.5));
        assert!(b1 < 0.0 && b2 < 0.0, "b1={b1} b2={b2}");
    }

    #[test]
    fn joint_ratio_moves_the_joint_along_the_chord() {
        // The chord of the C configuration runs along +x.
        let mut previous = f64::NEG_INFINITY;
        for ratio in [0.05, 0.2, 0.4, 0.6, 0.8, 0.95] {
            let (jx, _, _, _) = expect_biarc(solve(c_curve(), ratio));
            assert!(jx > previous, "ratio={ratio} jx={jx} previous={previous}");
            previous = jx;
        }
    }

    #[test]
    fn balanced_joint_ratio_is_symmetric() {
        // The S configuration is invariant under a half turn about the
        // chord midpoint, so the balanced joint sits on that midpoint
        // and the two arcs mirror each other.
        let (jx, jy, b1, b2) = expect_biarc(solve(S_CURVE, 0.5));
        assert_point((jx, jy), (1.0, 1.0), "balanced joint");
        assert!((b1 + b2).abs() < TOL, "b1={b1} b2={b2}");
    }

    #[test]
    fn reversal_mirrors_the_solution() {
        let (x0, y0, tx0, ty0, x1, y1, tx1, ty1) = c_curve();
        let (jx, jy, b1, b2) =
            expect_biarc(biarc_from_hermite(x0, y0, tx0, ty0, x1, y1, tx1, ty1, 0.3));
        // Walk the same element backwards: swapped endpoints, negated
        // tangents, mirrored ratio.
        let (rjx, rjy, rb1, rb2) = expect_biarc(biarc_from_hermite(
            x1, y1, -tx1, -ty1, x0, y0, -tx0, -ty0, 0.7,
        ));
        assert_point((rjx, rjy), (jx, jy), "reversed joint");
        assert!((rb1 + b2).abs() < TOL, "rb1={rb1} b2={b2}");
        assert!((rb2 + b1).abs() < TOL, "rb2={rb2} b1={b1}");
    }

    #[test]
    fn near_single_arc_degrades_continuously() {
        // Perturb an exact single-arc configuration: the pair must sit
        // next to the sub-arcs of the arc it came from.
        let (x0, y0, x1, y1) = (0.0, 0.0, 2.0, 0.5);
        let bulge = 0.4;
        let (t0, t1) = segment_tangents(x0, y0, x1, y1, bulge);
        let nudge = 1e-6_f64;
        let end_angle = t1.1.atan2(t1.0) + nudge;
        let (jx, jy, b1, b2) = expect_biarc(biarc_from_hermite(
            x0,
            y0,
            t0.0,
            t0.1,
            x1,
            y1,
            end_angle.cos(),
            end_angle.sin(),
            0.5,
        ));
        let sweep = 4.0 * bulge.atan();
        let half_bulge = (sweep * 0.5 / 4.0).tan();
        assert!((b1 - half_bulge).abs() < 1e-4, "b1={b1} vs {half_bulge}");
        assert!((b2 - half_bulge).abs() < 1e-4, "b2={b2} vs {half_bulge}");
        let (cx, cy, radius, start_angle, sw) = arc_from_bulge(x0, y0, x1, y1, bulge);
        let mid = arc_point_at(cx, cy, radius, start_angle, sw, 0.5);
        assert!(
            (jx - mid.0).abs() < 1e-4 && (jy - mid.1).abs() < 1e-4,
            "joint=({jx},{jy}) vs {mid:?}"
        );
    }

    #[test]
    fn antiparallel_tangents_produce_a_finite_biarc() {
        // Out and back along the same line: no single arc works, but a
        // pair does — and it is still G1.
        let (jx, jy, b1, b2) = expect_biarc(biarc_from_hermite(
            0.0, 0.0, 1.0, 0.0, 2.0, 0.0, -1.0, 0.0, 0.5,
        ));
        assert!(b1.is_finite() && b2.is_finite(), "b1={b1} b2={b2}");
        assert_point((jx, jy), (1.0, -1.0), "antiparallel joint");
        assert!(b1 < 0.0 && b2 > 0.0, "b1={b1} b2={b2}");
        let (start, first_end) = segment_tangents(0.0, 0.0, jx, jy, b1);
        let (second_start, end) = segment_tangents(jx, jy, 2.0, 0.0, b2);
        assert_direction(start, (1.0, 0.0), "start tangent");
        assert_direction(end, (-1.0, 0.0), "end tangent");
        assert_direction(first_end, second_start, "joint tangent");
    }

    #[test]
    fn tangents_pointing_back_along_the_chord_stay_finite() {
        // No arc pair can honour tangents that both point away from the
        // chord along it; the guard keeps the result finite.
        let (jx, jy, b1, b2) = expect_biarc(biarc_from_hermite(
            0.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.5,
        ));
        assert!(
            jx.is_finite() && jy.is_finite() && b1.is_finite() && b2.is_finite(),
            "joint=({jx},{jy}) b1={b1} b2={b2}"
        );
    }

    #[test]
    fn degenerate_inputs_are_straight() {
        // Zero-length chord.
        assert_eq!(
            biarc_from_hermite(1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.5),
            BiarcShape::Straight
        );
        // Zero-length tangents.
        assert_eq!(
            biarc_from_hermite(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.5),
            BiarcShape::Straight
        );
        assert_eq!(
            biarc_from_hermite(0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.5),
            BiarcShape::Straight
        );
        // Non-finite inputs.
        assert_eq!(
            biarc_from_hermite(f64::NAN, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.5),
            BiarcShape::Straight
        );
        assert_eq!(
            biarc_from_hermite(0.0, 0.0, 1.0, 0.0, f64::INFINITY, 0.0, 1.0, 0.0, 0.5),
            BiarcShape::Straight
        );
        assert_eq!(
            biarc_from_hermite(0.0, 0.0, 1.0, f64::NAN, 1.0, 0.0, 1.0, 0.0, 0.5),
            BiarcShape::Straight
        );
    }

    #[test]
    fn joint_ratio_out_of_range_is_clamped() {
        for ratio in [0.0, 1.0, -5.0, 7.5] {
            let (jx, jy, b1, b2) = expect_biarc(solve(c_curve(), ratio));
            assert!(
                jx.is_finite() && jy.is_finite() && b1.is_finite() && b2.is_finite(),
                "ratio={ratio} joint=({jx},{jy}) b1={b1} b2={b2}"
            );
            assert!(jx > 0.0 && jx < 2.0, "ratio={ratio} jx={jx}");
        }
        // A non-finite ratio falls back to the balanced joint.
        assert_eq!(
            solve(c_curve(), f64::NAN),
            solve(c_curve(), 0.5),
            "NaN ratio must behave like 0.5"
        );
    }

    /// The inverse solve, applied to a configuration's own joint.
    fn recover(config: (f64, f64, f64, f64, f64, f64, f64, f64), jx: f64, jy: f64) -> Option<f64> {
        let (x0, y0, tx0, ty0, x1, y1, tx1, ty1) = config;
        biarc_joint_ratio(x0, y0, tx0, ty0, x1, y1, tx1, ty1, jx, jy)
    }

    #[test]
    fn joint_ratio_round_trips_through_the_solver() {
        // S_CURVE turns by zero in total, so it exercises the flattened
        // locus; c_curve exercises the curved one.
        for config in [S_CURVE, c_curve()] {
            for ratio in [0.05, 0.25, 0.5, 0.8, 0.95] {
                let (jx, jy, _, _) = expect_biarc(solve(config, ratio));
                let recovered = recover(config, jx, jy);
                match recovered {
                    Some(value) => assert!(
                        (value - ratio).abs() < TOL,
                        "ratio={ratio} recovered={value}"
                    ),
                    None => panic!("no ratio recovered for {ratio}"),
                }
            }
        }
    }

    #[test]
    fn joint_ratio_rejects_shapes_without_a_joint() {
        // A straight element: nothing to place.
        assert_eq!(
            biarc_joint_ratio(0.0, 0.0, 3.0, 4.0, 3.0, 4.0, 0.6, 0.8, 1.5, 2.0),
            None
        );
        // A single arc already meets both tangents.
        let (x0, y0, x1, y1) = (-1.3, 0.7, 2.4, -0.9);
        let (t0, t1) = segment_tangents(x0, y0, x1, y1, 0.4);
        assert_eq!(
            biarc_joint_ratio(x0, y0, t0.0, t0.1, x1, y1, t1.0, t1.1, 0.5, 0.0),
            None
        );
    }

    #[test]
    fn joint_ratio_rejects_a_joint_off_the_locus() {
        let (jx, jy, _, _) = expect_biarc(solve(c_curve(), 0.4));
        assert!(recover(c_curve(), jx, jy).is_some(), "the locus joint");
        assert_eq!(recover(c_curve(), jx, jy + 0.2), None, "displaced joint");
        // An endpoint is on the locus but outside the addressable range.
        assert_eq!(recover(c_curve(), 2.0, 0.0), None, "end point");
    }

    #[test]
    fn joint_ratio_rejects_degenerate_input() {
        let (jx, jy, _, _) = expect_biarc(solve(c_curve(), 0.4));
        // Non-finite joint.
        assert_eq!(recover(c_curve(), f64::NAN, jy), None);
        assert_eq!(recover(c_curve(), jx, f64::INFINITY), None);
        // Zero-length chord and zero-length tangent.
        assert_eq!(
            biarc_joint_ratio(1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.5),
            None
        );
        assert_eq!(
            biarc_joint_ratio(0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0, 1.0, 0.5),
            None
        );
        // A joint sitting on p0 has no direction to read.
        assert_eq!(recover(c_curve(), 0.0, 0.0), None);
    }

    #[test]
    fn wrap_signed_covers_the_half_open_range() {
        assert!((wrap_signed(0.0)).abs() < TOL);
        assert!((wrap_signed(PI) - PI).abs() < TOL);
        assert!((wrap_signed(-PI) - PI).abs() < TOL);
        assert!((wrap_signed(1.5 * PI) + 0.5 * PI).abs() < TOL);
        assert!((wrap_signed(-3.0 * PI) - PI).abs() < TOL);
    }
}
