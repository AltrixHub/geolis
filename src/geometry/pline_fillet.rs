//! Corner filleting of [`Pline`]s.
//!
//! Replaces each interior corner between two straight segments with a
//! tangent circular arc (bulge edge) of the requested radius.

use crate::error::{GeometryError, Result};

use super::pline::{Pline, PlineVertex};

const EPS: f64 = 1e-12;

impl Pline {
    /// Fillet every corner between two straight segments with a
    /// tangent arc of `radius`. Corners adjacent to an existing arc
    /// segment are left untouched; near-collinear corners are skipped.
    ///
    /// # Errors
    ///
    /// Returns an error when `radius` is not strictly positive and
    /// finite, or when a fillet would consume more than half of an
    /// adjacent segment (the tangent points of neighbouring fillets
    /// would overlap).
    pub fn fillet(&self, radius: f64) -> Result<Pline> {
        if !radius.is_finite() || radius <= 0.0 {
            return Err(GeometryError::Degenerate(format!(
                "fillet radius must be strictly positive, got {radius}"
            ))
            .into());
        }
        let n = self.vertices.len();
        if n < 3 {
            return Ok(self.clone());
        }
        let seg_count = self.segment_count();
        // Interior corners: vertex i joins segment (i-1) and segment i.
        // Open plines exclude the endpoints; closed plines fillet all.
        let corner_range: Vec<usize> = if self.closed {
            (0..n).collect()
        } else {
            (1..n - 1).collect()
        };

        let mut out: Vec<PlineVertex> = Vec::with_capacity(n * 2);
        // Track trims per segment end so we can validate overlap.
        let mut result = self.vertices.clone();

        // We build a fresh vertex list corner by corner.
        let vertex = |i: usize| &self.vertices[i % n];
        let seg_vec = |i: usize| {
            let a = vertex(i);
            let b = vertex(i + 1);
            (b.x - a.x, b.y - a.y)
        };
        let seg_len = |i: usize| {
            let (dx, dy) = seg_vec(i);
            (dx * dx + dy * dy).sqrt()
        };

        // For each vertex, compute the trim distance (0 for untouched).
        let mut trims = vec![0.0_f64; n];
        let mut sweeps = vec![0.0_f64; n];
        for &i in &corner_range {
            let prev_seg = (i + n - 1) % n;
            let next_seg = i % n;
            if !self.closed && (i == 0 || i == n - 1) {
                continue;
            }
            if prev_seg >= seg_count || next_seg >= seg_count {
                continue;
            }
            // Only fillet corners between straight segments.
            if vertex(prev_seg).bulge.abs() > EPS || vertex(next_seg).bulge.abs() > EPS {
                continue;
            }
            let (ax, ay) = seg_vec(prev_seg);
            let (bx, by) = seg_vec(next_seg);
            let (la, lb) = (seg_len(prev_seg), seg_len(next_seg));
            if la < EPS || lb < EPS {
                continue;
            }
            let (uax, uay) = (ax / la, ay / la);
            let (ubx, uby) = (bx / lb, by / lb);
            let cross = uax * uby - uay * ubx;
            let dot = uax * ubx + uay * uby;
            let turn = cross.atan2(dot); // signed turn angle at the corner
            if turn.abs() < 1e-9 {
                continue; // collinear
            }
            let half = (std::f64::consts::PI - turn.abs()) / 2.0;
            let trim = radius / half.tan();
            trims[i] = trim;
            sweeps[i] = if cross > 0.0 {
                std::f64::consts::PI - 2.0 * half
            } else {
                -(std::f64::consts::PI - 2.0 * half)
            };
        }

        // Validate segment budgets: trims at both ends must fit.
        for seg in 0..seg_count {
            let start_trim = trims[(seg + 1) % n]; // corner at segment end
            let end_trim = trims[seg]; // corner at segment start
            if end_trim + start_trim > seg_len(seg) + EPS {
                return Err(GeometryError::Degenerate(format!(
                    "fillet radius {radius} does not fit on segment {seg}"
                ))
                .into());
            }
        }

        // Emit vertices: for each original vertex, either passthrough or
        // (arc start with bulge) + (arc end).
        result.clear();
        for i in 0..n {
            let v = vertex(i);
            if trims[i] <= EPS {
                result.push(*v);
                continue;
            }
            let prev_seg = (i + n - 1) % n;
            let next_seg = i % n;
            let (ax, ay) = seg_vec(prev_seg);
            let la = seg_len(prev_seg);
            let (bx, by) = seg_vec(next_seg);
            let lb = seg_len(next_seg);
            let t = trims[i];
            // Arc start: back along the incoming segment.
            let start = PlineVertex::new(
                v.x - ax / la * t,
                v.y - ay / la * t,
                (sweeps[i] / 4.0).tan(),
            );
            // Arc end: forward along the outgoing segment (straight from here).
            let end = PlineVertex::line(v.x + bx / lb * t, v.y + by / lb * t);
            result.push(start);
            result.push(end);
        }
        out.extend(result);
        Ok(Pline {
            vertices: out,
            closed: self.closed,
        })
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use std::f64::consts::PI;

    #[test]
    fn fillets_a_right_angle_with_a_quarter_arc() {
        let pline = Pline::from_points(
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(4.0, 0.0, 0.0),
                Point3::new(4.0, 4.0, 0.0),
            ],
            false,
        );
        let filleted = pline.fillet(1.0).unwrap();
        // 4 vertices now: start, arc start, arc end, end.
        assert_eq!(filleted.vertices.len(), 4);
        // Arc start sits 1.0 before the corner on the incoming segment.
        assert!((filleted.vertices[1].x - 3.0).abs() < 1e-9);
        assert!(filleted.vertices[1].y.abs() < 1e-9);
        // Quarter arc: |bulge| = tan(90deg / 4).
        assert!((filleted.vertices[1].bulge.abs() - (PI / 8.0).tan()).abs() < 1e-9);
        // Total length: two trimmed legs + quarter arc of radius 1.
        let expected = 3.0 + 3.0 + PI / 2.0;
        assert!((filleted.arc_length() - expected).abs() < 1e-9);
    }

    #[test]
    fn closed_square_fillets_every_corner() {
        let square = Pline::from_points(
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(4.0, 0.0, 0.0),
                Point3::new(4.0, 4.0, 0.0),
                Point3::new(0.0, 4.0, 0.0),
            ],
            true,
        );
        let filleted = square.fillet(1.0).unwrap();
        assert_eq!(filleted.vertices.len(), 8);
        // Perimeter: 4 * (4 - 2) straight + full circle of radius 1.
        let expected = 8.0 + 2.0 * PI;
        assert!((filleted.arc_length() - expected).abs() < 1e-9);
    }

    #[test]
    fn oversized_radius_is_rejected() {
        let pline = Pline::from_points(
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
                Point3::new(2.0, 2.0, 0.0),
            ],
            false,
        );
        assert!(pline.fillet(10.0).is_err());
        assert!(pline.fillet(0.0).is_err());
        assert!(pline.fillet(f64::NAN).is_err());
    }

    #[test]
    fn collinear_corners_pass_through() {
        let pline = Pline::from_points(
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
                Point3::new(4.0, 0.0, 0.0),
            ],
            false,
        );
        let filleted = pline.fillet(0.5).unwrap();
        assert_eq!(filleted.vertices.len(), 3);
        assert!((filleted.arc_length() - 4.0).abs() < 1e-9);
    }
}
