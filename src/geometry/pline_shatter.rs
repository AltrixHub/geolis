//! Splitting a [`Pline`] at arc-length stations, plus the turn angle at
//! each of its joints.
//!
//! Both are polyline-structure primitives rather than stair / wall / road
//! specifics: `shatter_at_lengths` is the kernel behind "cut this path
//! into the pieces between these positions" (paving runs, panel splits,
//! chainage stations), and `turn_angles` is what every mitre, chamfer,
//! fillet-radius or corner-trim rule measures before it decides how much
//! material a joint needs.
//!
//! Arc segments survive both operations exactly: a split arc keeps its
//! centre and radius because the sub-segment's bulge is recomputed from
//! the swept fraction (`bulge = tan(sweep / 4)`), and the turn angle is
//! read off the true end / start tangents rather than the chord.

use crate::error::{GeometryError, Result};
use crate::math::arc_2d::sub_arc_bulge;
use crate::math::Vector3;

use super::pline::{Pline, PlineVertex};

/// Arc-length tolerance for treating two stations (or a station and a
/// polyline end) as the same position.
const STATION_EPS: f64 = 1e-9;

/// Fraction tolerance along one edge, used at both of its ends: a
/// sub-segment shorter than it is empty, and a station within it of an
/// edge's far end sits ON the joint rather than inside the edge (which
/// is where a piece opening there is re-seated onto the next edge).
const FRACTION_EPS: f64 = 1e-12;

/// Signed angle (radians, CCW positive) from `from` to `to` in the XY
/// plane. Returns `0.0` when either vector is degenerate.
fn signed_turn(from: Vector3, to: Vector3) -> f64 {
    let cross = from.x * to.y - from.y * to.x;
    let dot = from.x * to.x + from.y * to.y;
    if cross.abs() < f64::EPSILON && dot.abs() < f64::EPSILON {
        return 0.0;
    }
    cross.atan2(dot)
}

/// One piece of a shattered polyline, positioned on its source.
///
/// The arc-length bounds are what make a piece usable beyond its own
/// XY: a caller carrying per-vertex attributes (elevation, chainage,
/// material runs) re-samples the source at `start_along + offset` to
/// recover them.
#[derive(Clone, Debug)]
pub struct PlineSpan {
    /// The piece itself, always an open polyline.
    pub pline: Pline,
    /// Arc length of the piece's first vertex along the source.
    pub start_along: f64,
    /// Arc length of the piece's last vertex along the source. A closed
    /// polyline's seam-crossing piece reports a value beyond the
    /// source's `arc_length()` (it wrapped once).
    pub end_along: f64,
}

impl Pline {
    /// Signed turn angle (radians, CCW positive) at each **joint** — a
    /// vertex where two segments meet.
    ///
    /// An open polyline of `n` vertices has `n - 2` joints (its two free
    /// ends are not joints and are therefore not reported); a closed
    /// polyline has one per vertex, the seam included. The angle is
    /// measured from the incoming segment's end tangent to the outgoing
    /// segment's start tangent, so an arc's true tangent is used rather
    /// than its chord, and a straight continuation reports `0.0`.
    #[must_use]
    pub fn turn_angles(&self) -> Vec<f64> {
        let segments = self.segment_count();
        if segments < 2 {
            return Vec::new();
        }
        // Open: joint j sits between segment j and segment j + 1, so
        // there are `segments - 1` of them. Closed: every vertex is a
        // joint, and joint 0 (the seam) joins the last segment to the
        // first.
        let joints = if self.closed { segments } else { segments - 1 };
        (0..joints)
            .map(|j| {
                let (incoming, outgoing) = if self.closed {
                    ((j + segments - 1) % segments, j)
                } else {
                    (j, j + 1)
                };
                signed_turn(
                    self.sample_segment(incoming, 1.0, 0.0).tangent,
                    self.sample_segment(outgoing, 0.0, 0.0).tangent,
                )
            })
            .collect()
    }

    /// Split the polyline into the sub-polylines lying between the given
    /// arc-length `stations`.
    ///
    /// Stations are clamped into `0..=arc_length()`, sorted and deduped
    /// before use, so a caller may hand over an unordered or
    /// out-of-range list. Every emitted piece is an **open** polyline
    /// carrying its arc-length bounds on the source (see [`PlineSpan`]),
    /// and arcs keep their exact centre and radius across a cut.
    ///
    /// - **Open** polyline: the ends are implicit boundaries, so
    ///   `k` interior stations produce `k + 1` pieces. With no interior
    ///   station the whole polyline is returned as the single piece.
    /// - **Closed** polyline: the stations are the only boundaries, so
    ///   `k >= 1` stations produce `k` pieces (the last wraps through the
    ///   seam). With no station the loop is returned opened at its seam.
    ///
    /// # Errors
    ///
    /// Returns an error when the polyline has no length or a station is
    /// not finite.
    pub fn shatter_at_lengths(&self, stations: &[f64]) -> Result<Vec<PlineSpan>> {
        let total = self.arc_length();
        if total <= STATION_EPS {
            return Err(GeometryError::Degenerate("polyline has no length".into()).into());
        }
        if let Some(bad) = stations.iter().copied().find(|s| !s.is_finite()) {
            return Err(GeometryError::Degenerate(format!(
                "shatter station must be finite, got {bad}"
            ))
            .into());
        }

        let mut cuts: Vec<f64> = stations.iter().map(|s| s.clamp(0.0, total)).collect();
        cuts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        cuts.dedup_by(|a, b| (*a - *b).abs() <= STATION_EPS);

        if self.closed {
            // The seam is not a boundary on its own: with no station the
            // loop opens at its seam, otherwise each station starts a
            // piece and the last one wraps back to the first.
            if cuts.is_empty() {
                return Ok(vec![self.piece(0.0, total, total)?]);
            }
            return (0..cuts.len())
                .map(|k| {
                    let start = cuts[k];
                    let end = if k + 1 < cuts.len() {
                        cuts[k + 1]
                    } else {
                        cuts[0] + total
                    };
                    self.piece(start, end, total)
                })
                .collect();
        }

        // Open: drop the stations that land on the free ends — they add
        // no cut, only empty pieces.
        cuts.retain(|s| *s > STATION_EPS && *s < total - STATION_EPS);
        let mut bounds = Vec::with_capacity(cuts.len() + 2);
        bounds.push(0.0);
        bounds.extend(cuts);
        bounds.push(total);
        bounds
            .windows(2)
            .map(|pair| self.piece(pair[0], pair[1], total))
            .collect()
    }

    /// The open sub-polyline from arc length `start` to `end`.
    ///
    /// `end` may exceed `total` by up to one lap: that is how a closed
    /// polyline's seam-crossing piece is requested.
    fn piece(&self, start: f64, end: f64, total: f64) -> Result<PlineSpan> {
        let wrapped = end > total + STATION_EPS;
        let end_along = if wrapped { end - total } else { end };
        let mut from = self.sample_at_length(start)?;
        let to = self.sample_at_length(end_along)?;
        let segments = self.segment_count();

        // A station that lands exactly ON a joint is reported as the END
        // of the edge arriving at it. A piece OPENING there walks the
        // edge that LEAVES it, so re-seat the sample on that edge: the
        // point is the same one, and the piece would otherwise open with
        // a zero-length first segment — a vertex duplicated at the joint,
        // which every consumer that measures per-edge geometry then reads
        // as a degenerate edge.
        if from.edge_fraction > 1.0 - FRACTION_EPS {
            from.edge_index = (from.edge_index + 1) % segments;
            from.edge_fraction = 0.0;
        }

        let mut steps = (to.edge_index + segments - from.edge_index) % segments;
        if steps == 0 && (wrapped || to.edge_fraction + FRACTION_EPS < from.edge_fraction) {
            // Same segment, but the piece walks the whole loop first.
            steps = segments;
        }

        let mut vertices: Vec<PlineVertex> = Vec::with_capacity(steps + 2);
        let mut edge = from.edge_index;
        let mut fraction = from.edge_fraction;
        let mut point = (from.point.x, from.point.y);
        for _ in 0..steps {
            vertices.push(PlineVertex::new(
                point.0,
                point.1,
                sub_arc_bulge(self.vertices[edge].bulge, fraction, 1.0),
            ));
            edge = (edge + 1) % segments;
            fraction = 0.0;
            let next = self.vertices[edge];
            point = (next.x, next.y);
        }
        if to.edge_fraction > fraction + FRACTION_EPS {
            vertices.push(PlineVertex::new(
                point.0,
                point.1,
                sub_arc_bulge(self.vertices[edge].bulge, fraction, to.edge_fraction),
            ));
        }
        vertices.push(PlineVertex::line(to.point.x, to.point.y));
        Ok(PlineSpan {
            pline: Pline {
                vertices,
                closed: false,
            },
            start_along: start,
            end_along: end,
        })
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use std::f64::consts::{FRAC_PI_2, PI};

    const TOL: f64 = 1e-9;

    /// (0,0) → (4,0) → (4,3): one right-angle joint, total length 7.
    fn l_path(closed: bool) -> Pline {
        Pline::from_points(
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(4.0, 0.0, 0.0),
                Point3::new(4.0, 3.0, 0.0),
            ],
            closed,
        )
    }

    /// CCW semicircle (0,0) → (2,0) through (1,-1); arc length π.
    fn semicircle() -> Pline {
        Pline {
            vertices: vec![PlineVertex::new(0.0, 0.0, 1.0), PlineVertex::line(2.0, 0.0)],
            closed: false,
        }
    }

    // ── turn angles ─────────────────────────────────────────────────

    #[test]
    fn open_polyline_reports_only_its_interior_joints() {
        let angles = l_path(false).turn_angles();
        assert_eq!(angles.len(), 1, "3 vertices → 1 joint");
        assert!(
            (angles[0] - FRAC_PI_2).abs() < TOL,
            "a left turn is +90°, got {}",
            angles[0]
        );
    }

    #[test]
    fn right_turns_are_negative() {
        let path = Pline::from_points(
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(4.0, 0.0, 0.0),
                Point3::new(4.0, -3.0, 0.0),
            ],
            false,
        );
        assert!((path.turn_angles()[0] + FRAC_PI_2).abs() < TOL);
    }

    #[test]
    fn collinear_continuation_does_not_turn() {
        let path = Pline::from_points(
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(3.0, 0.0, 0.0),
            ],
            false,
        );
        assert!(path.turn_angles()[0].abs() < TOL);
    }

    #[test]
    fn closed_polyline_reports_every_vertex_including_the_seam() {
        let square = Pline::from_points(
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(4.0, 0.0, 0.0),
                Point3::new(4.0, 4.0, 0.0),
                Point3::new(0.0, 4.0, 0.0),
            ],
            true,
        );
        let angles = square.turn_angles();
        assert_eq!(angles.len(), 4, "the seam is a joint too");
        for angle in angles {
            assert!((angle - FRAC_PI_2).abs() < TOL, "CCW square turns +90°");
        }
    }

    #[test]
    fn arc_joints_measure_the_true_tangent_not_the_chord() {
        // Semicircle (tangent at its end is +Y) into a straight +Y leg:
        // the tangents agree, so the joint does not turn — a chord-based
        // measurement would report a large angle instead.
        let path = Pline {
            vertices: vec![
                PlineVertex::new(0.0, 0.0, 1.0),
                PlineVertex::line(2.0, 0.0),
                PlineVertex::line(2.0, 3.0),
            ],
            closed: false,
        };
        assert!(path.turn_angles()[0].abs() < TOL);
    }

    #[test]
    fn a_single_segment_has_no_joints() {
        let line = Pline::from_points(
            &[Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            false,
        );
        assert!(line.turn_angles().is_empty());
    }

    // ── shatter ─────────────────────────────────────────────────────

    fn ends(span: &PlineSpan) -> ((f64, f64), (f64, f64)) {
        let first = span.pline.vertices.first().unwrap();
        let last = span.pline.vertices.last().unwrap();
        ((first.x, first.y), (last.x, last.y))
    }

    fn lengths(spans: &[PlineSpan]) -> Vec<f64> {
        spans.iter().map(|s| s.pline.arc_length()).collect()
    }

    #[test]
    fn open_polyline_splits_into_one_more_piece_than_stations() {
        let pieces = l_path(false).shatter_at_lengths(&[2.0, 5.0]).unwrap();
        assert_eq!(pieces.len(), 3);
        let lengths = lengths(&pieces);
        assert!((lengths[0] - 2.0).abs() < TOL, "{lengths:?}");
        assert!((lengths[1] - 3.0).abs() < TOL, "{lengths:?}");
        assert!((lengths[2] - 2.0).abs() < TOL, "{lengths:?}");
        // Every piece knows where it sits on the source.
        let bounds: Vec<(f64, f64)> = pieces
            .iter()
            .map(|s| (s.start_along, s.end_along))
            .collect();
        assert_eq!(bounds, vec![(0.0, 2.0), (2.0, 5.0), (5.0, 7.0)]);
    }

    #[test]
    fn pieces_tile_the_source_end_to_end() {
        let pieces = l_path(false).shatter_at_lengths(&[2.0, 5.0]).unwrap();
        let (start, _) = ends(&pieces[0]);
        assert!(start.0.abs() < TOL && start.1.abs() < TOL);
        for pair in pieces.windows(2) {
            let (_, left_end) = ends(&pair[0]);
            let (right_start, _) = ends(&pair[1]);
            assert!((left_end.0 - right_start.0).abs() < TOL);
            assert!((left_end.1 - right_start.1).abs() < TOL);
        }
        let (_, end) = ends(pieces.last().unwrap());
        assert!((end.0 - 4.0).abs() < TOL && (end.1 - 3.0).abs() < TOL);
    }

    #[test]
    fn a_piece_crossing_a_joint_keeps_the_joint_vertex() {
        // 3.0 .. 5.0 straddles the corner at arc length 4.
        let pieces = l_path(false).shatter_at_lengths(&[3.0, 5.0]).unwrap();
        assert_eq!(pieces[1].pline.vertices.len(), 3, "start, corner, end");
        let corner = pieces[1].pline.vertices[1];
        assert!((corner.x - 4.0).abs() < TOL && corner.y.abs() < TOL);
    }

    #[test]
    fn stations_are_sorted_deduped_and_clamped() {
        let pieces = l_path(false)
            .shatter_at_lengths(&[5.0, -3.0, 2.0, 2.0, 99.0])
            .unwrap();
        // -3 and 99 clamp onto the ends (dropped), 2.0 dedups.
        assert_eq!(pieces.len(), 3);
    }

    #[test]
    fn no_station_returns_the_whole_open_polyline() {
        let pieces = l_path(false).shatter_at_lengths(&[]).unwrap();
        assert_eq!(pieces.len(), 1);
        assert!((pieces[0].pline.arc_length() - 7.0).abs() < TOL);
    }

    #[test]
    fn splitting_an_arc_preserves_its_radius() {
        let pieces = semicircle().shatter_at_lengths(&[PI / 2.0]).unwrap();
        assert_eq!(pieces.len(), 2);
        for piece in &pieces {
            assert!(
                (piece.pline.arc_length() - PI / 2.0).abs() < TOL,
                "each half sweeps a quarter turn: {}",
                piece.pline.arc_length()
            );
            // A quarter turn of a unit circle: bulge = tan(90° / 4).
            let bulge = piece.pline.vertices[0].bulge;
            assert!(
                (bulge - (FRAC_PI_2 / 4.0).tan()).abs() < TOL,
                "split arc must stay an arc, got bulge {bulge}"
            );
        }
        // The cut lands on the bottom of the circle.
        let (_, mid) = ends(&pieces[0]);
        assert!((mid.0 - 1.0).abs() < TOL && (mid.1 + 1.0).abs() < TOL);
    }

    #[test]
    fn closed_polyline_yields_one_piece_per_station() {
        let square = Pline::from_points(
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(4.0, 0.0, 0.0),
                Point3::new(4.0, 4.0, 0.0),
                Point3::new(0.0, 4.0, 0.0),
            ],
            true,
        );
        let pieces = square.shatter_at_lengths(&[2.0, 8.0]).unwrap();
        assert_eq!(pieces.len(), 2, "2 stations cut a loop into 2 pieces");
        let lengths = lengths(&pieces);
        assert!((lengths[0] - 6.0).abs() < TOL, "{lengths:?}");
        assert!((lengths[1] - 10.0).abs() < TOL, "seam-crossing piece");
        // The wrap piece reports bounds past the source length.
        assert!((pieces[1].start_along - 8.0).abs() < TOL);
        assert!((pieces[1].end_along - 18.0).abs() < TOL);
        // The wrapping piece starts at station 8 and ends back at 2.
        let (start, end) = ends(&pieces[1]);
        assert!((start.0 - 4.0).abs() < TOL && (start.1 - 4.0).abs() < TOL);
        assert!((end.0 - 2.0).abs() < TOL && end.1.abs() < TOL);
    }

    /// A station that lands exactly ON a joint cuts there and NOTHING
    /// else: the piece that opens at the joint starts on the edge that
    /// leaves it, not with a duplicate of the joint vertex.
    ///
    /// The duplicate is invisible in XY — the piece traces the same path
    /// either way — and lethal to every consumer that measures per edge:
    /// a leading zero-length edge is a segment whose direction, length
    /// and parameterisation are all undefined.
    #[test]
    fn a_station_on_a_joint_opens_the_next_piece_on_the_edge_that_leaves_it() {
        // The L's joint is at station 4; cut a hair either side of it too,
        // so the same walk is read with the joint interior to a piece.
        for stations in [vec![4.0], vec![3.5, 4.0, 4.5]] {
            let pieces = l_path(false).shatter_at_lengths(&stations).unwrap();
            assert_eq!(pieces.len(), stations.len() + 1, "{stations:?}");
            for piece in &pieces {
                for (index, pair) in piece.pline.vertices.windows(2).enumerate() {
                    let (dx, dy) = (pair[1].x - pair[0].x, pair[1].y - pair[0].y);
                    assert!(
                        dx.hypot(dy) > TOL,
                        "{stations:?}: edge {index} of {:?} has no length",
                        piece.pline.vertices,
                    );
                }
            }
        }
        // …and the cut is still exactly at the joint, on both sides.
        let pieces = l_path(false).shatter_at_lengths(&[4.0]).unwrap();
        let lengths = lengths(&pieces);
        assert!((lengths[0] - 4.0).abs() < TOL && (lengths[1] - 3.0).abs() < TOL);
        let (_, first_end) = ends(&pieces[0]);
        let (second_start, _) = ends(&pieces[1]);
        assert!((first_end.0 - 4.0).abs() < TOL && first_end.1.abs() < TOL);
        assert!((second_start.0 - 4.0).abs() < TOL && second_start.1.abs() < TOL);
    }

    #[test]
    fn degenerate_and_non_finite_inputs_are_rejected() {
        assert!(Pline::from_points(&[], false)
            .shatter_at_lengths(&[1.0])
            .is_err());
        assert!(l_path(false).shatter_at_lengths(&[f64::NAN]).is_err());
    }
}
