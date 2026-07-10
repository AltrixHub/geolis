//! Arc-length sampling and division of [`Pline`]s.
//!
//! Exact per-segment arc lengths (lines: chord length, arcs:
//! `radius * |sweep|`) drive parameterization by distance along the
//! polyline — the kernel primitive behind curve evaluation, division,
//! and frame generation in downstream node graphs.

use crate::error::{GeometryError, Result};
use crate::math::arc_2d::{arc_from_bulge, arc_point_at};
use crate::math::{Point3, Vector3};

use super::pline::Pline;

const EPS: f64 = 1e-12;

/// Guard against pathological `divide_by_length` requests
/// (`total / segment_length` explosions) before allocating.
const MAX_DIVISIONS: usize = 100_000;

/// One sample on a polyline, positioned by arc length.
///
/// `edge_index` / `edge_fraction` identify the source segment and the
/// arc-length fraction within it, so callers that carry per-vertex
/// attributes (e.g. elevation) can interpolate them exactly.
#[derive(Clone, Debug)]
pub struct PlineSample {
    /// Position (z = 0; `Pline` is planar).
    pub point: Point3,
    /// Unit tangent in the XY plane (z = 0).
    pub tangent: Vector3,
    /// Index of the source segment.
    pub edge_index: usize,
    /// Arc-length fraction within the source segment, `0..=1`.
    pub edge_fraction: f64,
    /// Arc length from the polyline start.
    pub length_along: f64,
}

impl Pline {
    /// Total arc length (exact for arcs: `radius * |sweep|`).
    #[must_use]
    pub fn arc_length(&self) -> f64 {
        self.segment_arc_lengths().iter().sum()
    }

    /// Sample the polyline at arc length `s` from its start.
    ///
    /// # Errors
    ///
    /// Returns an error when the polyline has no length or `s` lies
    /// outside `0..=arc_length()`.
    pub fn sample_at_length(&self, s: f64) -> Result<PlineSample> {
        let lengths = self.segment_arc_lengths();
        let total: f64 = lengths.iter().sum();
        if total <= EPS {
            return Err(GeometryError::Degenerate("polyline has no length".into()).into());
        }
        if !s.is_finite() || s < -EPS || s > total + EPS {
            return Err(GeometryError::ParameterOutOfRange {
                parameter: "arc length",
                value: s,
                min: 0.0,
                max: total,
            }
            .into());
        }
        let s = s.clamp(0.0, total);

        let mut cum = 0.0;
        let mut last_real = 0;
        for (i, &len) in lengths.iter().enumerate() {
            if len <= EPS {
                continue;
            }
            last_real = i;
            if s <= cum + len + EPS {
                let fraction = ((s - cum) / len).clamp(0.0, 1.0);
                return Ok(self.sample_segment(i, fraction, s));
            }
            cum += len;
        }
        // Numeric tail: land on the end of the last real segment.
        Ok(self.sample_segment(last_real, 1.0, total))
    }

    /// Divide into `count` equal-arc-length segments. Open polylines
    /// yield `count + 1` samples (both ends included); closed polylines
    /// yield `count` samples (seam emitted once).
    ///
    /// # Errors
    ///
    /// Returns an error when `count == 0` or the polyline has no length.
    pub fn divide_by_count(&self, count: usize) -> Result<Vec<PlineSample>> {
        if count == 0 {
            return Err(GeometryError::Degenerate("division count must be >= 1".into()).into());
        }
        let total = self.arc_length();
        if total <= EPS {
            return Err(GeometryError::Degenerate("polyline has no length".into()).into());
        }
        let sample_count = if self.closed { count } else { count + 1 };
        #[allow(clippy::cast_precision_loss)]
        let samples = (0..sample_count)
            .map(|i| self.sample_at_length(total * (i as f64) / (count as f64)))
            .collect::<Result<Vec<_>>>()?;
        Ok(samples)
    }

    /// Samples every `segment_length` of arc length starting at the
    /// polyline start (the end is included only when it falls on an
    /// exact multiple, matching the GH `Divide Length` contract).
    ///
    /// # Errors
    ///
    /// Returns an error when `segment_length` is not strictly positive
    /// and finite, the polyline has no length, or the request would
    /// produce more than [`MAX_DIVISIONS`] samples.
    pub fn divide_by_length(&self, segment_length: f64) -> Result<Vec<PlineSample>> {
        if !segment_length.is_finite() || segment_length <= 0.0 {
            return Err(GeometryError::Degenerate(format!(
                "segment length must be strictly positive, got {segment_length}"
            ))
            .into());
        }
        let total = self.arc_length();
        if total <= EPS {
            return Err(GeometryError::Degenerate("polyline has no length".into()).into());
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let steps = (total / segment_length + EPS).floor() as usize;
        if steps + 1 > MAX_DIVISIONS {
            return Err(GeometryError::Degenerate(format!(
                "division would produce {} samples (max {MAX_DIVISIONS})",
                steps + 1
            ))
            .into());
        }
        // On a closed polyline the seam sits at both 0 and `total`; skip
        // the duplicate when the last step lands exactly on the end.
        #[allow(clippy::cast_precision_loss)]
        let last = if self.closed && ((steps as f64) * segment_length - total).abs() <= EPS {
            steps.saturating_sub(1)
        } else {
            steps
        };
        #[allow(clippy::cast_precision_loss)]
        let samples = (0..=last)
            .map(|k| self.sample_at_length(((k as f64) * segment_length).min(total)))
            .collect::<Result<Vec<_>>>()?;
        Ok(samples)
    }

    /// Exact arc length per segment (lines: chord, arcs:
    /// `radius * |sweep|`, degenerate arcs fall back to the chord).
    fn segment_arc_lengths(&self) -> Vec<f64> {
        let n = self.vertices.len();
        (0..self.segment_count())
            .map(|i| {
                let v0 = &self.vertices[i];
                let v1 = &self.vertices[(i + 1) % n];
                let chord = ((v1.x - v0.x).powi(2) + (v1.y - v0.y).powi(2)).sqrt();
                if v0.bulge.abs() < EPS {
                    chord
                } else {
                    let (_, _, radius, _, sweep) =
                        arc_from_bulge(v0.x, v0.y, v1.x, v1.y, v0.bulge);
                    if radius < EPS {
                        chord
                    } else {
                        radius * sweep.abs()
                    }
                }
            })
            .collect()
    }

    /// Point + unit tangent at arc-length fraction `t` of segment `i`.
    fn sample_segment(&self, edge_index: usize, fraction: f64, length_along: f64) -> PlineSample {
        let n = self.vertices.len();
        let v0 = &self.vertices[edge_index];
        let v1 = &self.vertices[(edge_index + 1) % n];
        let (point, tangent) = if v0.bulge.abs() < EPS {
            let dx = v1.x - v0.x;
            let dy = v1.y - v0.y;
            let len = (dx * dx + dy * dy).sqrt();
            let (tx, ty) = if len > EPS {
                (dx / len, dy / len)
            } else {
                (1.0, 0.0)
            };
            (
                Point3::new(v0.x + dx * fraction, v0.y + dy * fraction, 0.0),
                Vector3::new(tx, ty, 0.0),
            )
        } else {
            let (cx, cy, radius, start_angle, sweep) =
                arc_from_bulge(v0.x, v0.y, v1.x, v1.y, v0.bulge);
            if radius < EPS {
                (Point3::new(v0.x, v0.y, 0.0), Vector3::new(1.0, 0.0, 0.0))
            } else {
                let (px, py) = arc_point_at(cx, cy, radius, start_angle, sweep, fraction);
                let angle = start_angle + sweep * fraction;
                let sign = sweep.signum();
                (
                    Point3::new(px, py, 0.0),
                    Vector3::new(-angle.sin() * sign, angle.cos() * sign, 0.0),
                )
            }
        };
        PlineSample {
            point,
            tangent,
            edge_index,
            edge_fraction: fraction,
            length_along,
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::pline::PlineVertex;
    use std::f64::consts::PI;

    const TOL: f64 = 1e-9;

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

    fn semicircle() -> Pline {
        Pline {
            vertices: vec![PlineVertex::new(0.0, 0.0, 1.0), PlineVertex::line(2.0, 0.0)],
            closed: false,
        }
    }

    #[test]
    fn arc_length_sums_lines_and_arcs_exactly() {
        assert!((l_path(false).arc_length() - 7.0).abs() < TOL);
        assert!((l_path(true).arc_length() - 12.0).abs() < TOL);
        assert!((semicircle().arc_length() - PI).abs() < TOL);
    }

    #[test]
    fn sample_walks_across_segments() {
        let sample = l_path(false).sample_at_length(5.0).unwrap();
        assert!((sample.point.x - 4.0).abs() < TOL);
        assert!((sample.point.y - 1.0).abs() < TOL);
        assert!(sample.tangent.x.abs() < TOL);
        assert!((sample.tangent.y - 1.0).abs() < TOL);
        assert_eq!(sample.edge_index, 1);
        assert!((sample.edge_fraction - 1.0 / 3.0).abs() < TOL);
        assert!((sample.length_along - 5.0).abs() < TOL);
    }

    #[test]
    fn sample_at_start_and_end() {
        let start = l_path(false).sample_at_length(0.0).unwrap();
        assert!((start.tangent.x - 1.0).abs() < TOL);
        assert_eq!(start.edge_index, 0);
        let end = l_path(false).sample_at_length(7.0).unwrap();
        assert!((end.point.x - 4.0).abs() < TOL);
        assert!((end.point.y - 3.0).abs() < TOL);
    }

    #[test]
    fn sample_rejects_out_of_range() {
        assert!(l_path(false).sample_at_length(-0.1).is_err());
        assert!(l_path(false).sample_at_length(7.1).is_err());
        assert!(Pline::from_points(&[], false).sample_at_length(0.0).is_err());
    }

    #[test]
    fn arc_tangents_follow_sweep_direction() {
        // CCW semicircle (0,0)→(2,0) through the bottom (1,-1).
        let arc = semicircle();
        let start = arc.sample_at_length(0.0).unwrap();
        assert!(start.tangent.x.abs() < TOL, "tx={}", start.tangent.x);
        assert!((start.tangent.y + 1.0).abs() < TOL, "ty={}", start.tangent.y);
        let bottom = arc.sample_at_length(PI / 2.0).unwrap();
        assert!((bottom.point.x - 1.0).abs() < TOL);
        assert!((bottom.point.y + 1.0).abs() < TOL);
        assert!((bottom.tangent.x - 1.0).abs() < TOL, "tx={}", bottom.tangent.x);
        assert!(bottom.tangent.y.abs() < TOL, "ty={}", bottom.tangent.y);
    }

    #[test]
    fn divide_by_count_open_includes_both_ends() {
        let samples = l_path(false).divide_by_count(4).unwrap();
        assert_eq!(samples.len(), 5);
        assert!((samples[0].point.x).abs() < TOL);
        assert!((samples[2].length_along - 3.5).abs() < TOL);
        let last = samples.last().unwrap();
        assert!((last.point.x - 4.0).abs() < TOL);
        assert!((last.point.y - 3.0).abs() < TOL);
    }

    #[test]
    fn divide_by_count_closed_emits_seam_once() {
        let square = Pline::from_points(
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(4.0, 0.0, 0.0),
                Point3::new(4.0, 4.0, 0.0),
                Point3::new(0.0, 4.0, 0.0),
            ],
            true,
        );
        let samples = square.divide_by_count(4).unwrap();
        assert_eq!(samples.len(), 4);
        assert!((samples[1].point.x - 4.0).abs() < TOL);
        assert!(samples[1].point.y.abs() < TOL);
        assert!((samples[3].point.y - 4.0).abs() < TOL);
        assert!(samples[3].point.x.abs() < TOL);
    }

    #[test]
    fn divide_by_count_rejects_zero_and_degenerate() {
        assert!(l_path(false).divide_by_count(0).is_err());
        assert!(Pline::from_points(&[], false).divide_by_count(2).is_err());
    }

    #[test]
    fn divide_by_length_excludes_inexact_end() {
        let line = Pline::from_points(
            &[Point3::new(0.0, 0.0, 0.0), Point3::new(5.0, 0.0, 0.0)],
            false,
        );
        let samples = line.divide_by_length(2.0).unwrap();
        assert_eq!(samples.len(), 3);
        assert!((samples[2].point.x - 4.0).abs() < TOL);
    }

    #[test]
    fn divide_by_length_rejects_bad_input() {
        let line = Pline::from_points(
            &[Point3::new(0.0, 0.0, 0.0), Point3::new(5.0, 0.0, 0.0)],
            false,
        );
        assert!(line.divide_by_length(0.0).is_err());
        assert!(line.divide_by_length(-1.0).is_err());
        assert!(line.divide_by_length(f64::NAN).is_err());
        assert!(line.divide_by_length(1e-9).is_err(), "MAX_DIVISIONS guard");
    }
}
