use crate::math::arc_2d::{arc_from_bulge, arc_point_at};
use crate::math::Point3;

/// Bulge-encoded polyline vertex for mixed line/arc segments.
///
/// `bulge = tan(sweep_angle / 4)`:
/// - `0` = straight line to next vertex
/// - `> 0` = counter-clockwise arc to next vertex
/// - `< 0` = clockwise arc to next vertex
/// - `|bulge| = 1` = semicircle
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlineVertex {
    pub x: f64,
    pub y: f64,
    pub bulge: f64,
}

impl PlineVertex {
    /// Creates a new vertex with the given coordinates and bulge.
    #[must_use]
    pub fn new(x: f64, y: f64, bulge: f64) -> Self {
        Self { x, y, bulge }
    }

    /// Creates a line vertex (bulge = 0).
    #[must_use]
    pub fn line(x: f64, y: f64) -> Self {
        Self { x, y, bulge: 0.0 }
    }
}

/// A polyline with mixed straight-line and circular-arc segments.
///
/// Each segment between consecutive vertices is either a line (bulge=0)
/// or a circular arc (bulge≠0). For closed polylines, the last vertex
/// connects back to the first.
#[derive(Debug, Clone)]
pub struct Pline {
    pub vertices: Vec<PlineVertex>,
    pub closed: bool,
}

impl Pline {
    /// Creates a `Pline` from `Point3` vertices with all-zero bulges (line segments only).
    #[must_use]
    pub fn from_points(points: &[Point3], closed: bool) -> Self {
        let vertices = points
            .iter()
            .map(|p| PlineVertex::line(p.x, p.y))
            .collect();
        Self { vertices, closed }
    }

    /// Converts this polyline to a list of `Point3` by tessellating arcs into line segments.
    ///
    /// `tolerance` controls the maximum deviation between the arc and its chord approximation.
    #[must_use]
    pub fn to_points(&self, tolerance: f64) -> Vec<Point3> {
        let n = self.vertices.len();
        if n == 0 {
            return Vec::new();
        }

        let seg_count = self.segment_count();
        let mut points = Vec::with_capacity(n * 2);

        for i in 0..seg_count {
            let v0 = &self.vertices[i];
            let v1 = &self.vertices[(i + 1) % n];

            // Always add start point of segment.
            if i == 0 {
                points.push(Point3::new(v0.x, v0.y, 0.0));
            }

            if v0.bulge.abs() < 1e-12 {
                // Straight line: just add endpoint.
                points.push(Point3::new(v1.x, v1.y, 0.0));
            } else {
                // Arc: tessellate into line segments.
                let (cx, cy, radius, start_angle, sweep) =
                    arc_from_bulge(v0.x, v0.y, v1.x, v1.y, v0.bulge);

                if radius < 1e-12 {
                    points.push(Point3::new(v1.x, v1.y, 0.0));
                    continue;
                }

                // Number of subdivisions based on tolerance.
                let n_sub = arc_subdivision_count(radius, sweep.abs(), tolerance);

                for j in 1..n_sub {
                    let t = f64::from(j) / f64::from(n_sub);
                    let (px, py) = arc_point_at(cx, cy, radius, start_angle, sweep, t);
                    points.push(Point3::new(px, py, 0.0));
                }
                points.push(Point3::new(v1.x, v1.y, 0.0));
            }
        }

        points
    }

    /// Returns a new polyline with vertices in reverse order and negated bulges.
    ///
    /// For a segment `v[i] → v[i+1]` with bulge `b`, the reversed segment
    /// `v[i+1] → v[i]` has bulge `-b` (arc direction flips).
    #[must_use]
    pub fn reversed(&self) -> Self {
        let m = self.vertices.len();
        if m == 0 {
            return self.clone();
        }
        let mut new_verts = Vec::with_capacity(m);
        for j in 0..m {
            let orig_idx = m - 1 - j;
            // In the reversed polyline, vertex j connects to vertex j+1,
            // which corresponds to the reverse of original segment (m-2-j).
            let bulge = if j < m - 1 {
                -self.vertices[m - 2 - j].bulge
            } else {
                0.0
            };
            new_verts.push(PlineVertex::new(
                self.vertices[orig_idx].x,
                self.vertices[orig_idx].y,
                bulge,
            ));
        }
        Self {
            vertices: new_verts,
            closed: self.closed,
        }
    }

    /// Returns the number of segments in this polyline.
    #[must_use]
    pub fn segment_count(&self) -> usize {
        let n = self.vertices.len();
        if n < 2 {
            return 0;
        }
        if self.closed {
            n
        } else {
            n - 1
        }
    }
}

/// Computes the number of line segments needed to approximate an arc
/// within the given tolerance.
fn arc_subdivision_count(radius: f64, abs_sweep: f64, tolerance: f64) -> u32 {
    if radius < 1e-12 || abs_sweep < 1e-12 || tolerance <= 0.0 {
        return 1;
    }
    // From the sagitta formula: sagitta = r * (1 - cos(θ/2))
    // For a given tolerance: θ = 2 * acos(1 - tolerance/r)
    let max_angle = if tolerance >= radius {
        std::f64::consts::PI
    } else {
        2.0 * (1.0 - tolerance / radius).acos()
    };
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let n = (abs_sweep / max_angle).ceil() as u32;
    n.max(1)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn from_points_creates_line_only() {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ];
        let pline = Pline::from_points(&pts, false);
        assert_eq!(pline.vertices.len(), 3);
        assert_eq!(pline.segment_count(), 2);
        for v in &pline.vertices {
            assert!((v.bulge).abs() < 1e-12);
        }
    }

    #[test]
    fn from_points_closed() {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ];
        let pline = Pline::from_points(&pts, true);
        assert_eq!(pline.segment_count(), 3); // 3 sides of triangle
    }

    #[test]
    fn to_points_line_only_roundtrip() {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ];
        let pline = Pline::from_points(&pts, false);
        let result = pline.to_points(0.01);
        assert_eq!(result.len(), 3);
        for (a, b) in result.iter().zip(pts.iter()) {
            assert!((a.x - b.x).abs() < 1e-10);
            assert!((a.y - b.y).abs() < 1e-10);
        }
    }

    #[test]
    fn to_points_semicircle_arc() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::new(0.0, 0.0, 1.0), // semicircle
                PlineVertex::new(2.0, 0.0, 0.0),
            ],
            closed: false,
        };
        let pts = pline.to_points(0.01);
        // Should have start, some intermediate points, and end.
        assert!(pts.len() > 2, "expected more than 2 points, got {}", pts.len());
        // First and last points should match vertices.
        assert!((pts[0].x).abs() < 1e-10);
        assert!((pts[0].y).abs() < 1e-10);
        assert!((pts.last().unwrap().x - 2.0).abs() < 1e-10);
        assert!((pts.last().unwrap().y).abs() < 1e-10);
    }

    #[test]
    fn segment_count_empty() {
        let pline = Pline {
            vertices: vec![],
            closed: false,
        };
        assert_eq!(pline.segment_count(), 0);
    }

    #[test]
    fn segment_count_single_vertex() {
        let pline = Pline {
            vertices: vec![PlineVertex::line(0.0, 0.0)],
            closed: false,
        };
        assert_eq!(pline.segment_count(), 0);
    }

    #[test]
    fn reversed_line_only() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(1.0, 0.0),
                PlineVertex::line(1.0, 1.0),
            ],
            closed: false,
        };
        let rev = pline.reversed();
        assert_eq!(rev.vertices.len(), 3);
        assert!((rev.vertices[0].x - 1.0).abs() < 1e-12);
        assert!((rev.vertices[0].y - 1.0).abs() < 1e-12);
        assert!((rev.vertices[2].x).abs() < 1e-12);
        assert!((rev.vertices[2].y).abs() < 1e-12);
        // All bulges should be 0 for line-only.
        for v in &rev.vertices {
            assert!(v.bulge.abs() < 1e-12);
        }
    }

    #[test]
    fn reversed_with_arc() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::new(2.0, 0.0, 1.0), // semicircle CCW
                PlineVertex::line(4.0, 0.0),
            ],
            closed: false,
        };
        let rev = pline.reversed();
        assert_eq!(rev.vertices.len(), 3);
        // Reversed: (4,0) → (2,0) → (0,0)
        // Seg 0 (4→2): reverse of original seg 1 (v[1].bulge=1.0) → bulge = -1
        // Seg 1 (2→0): reverse of original seg 0 (v[0].bulge=0) → bulge = 0
        assert!((rev.vertices[0].bulge - (-1.0)).abs() < 1e-12); // (4,0), CW semicircle to (2,0)
        assert!(rev.vertices[1].bulge.abs() < 1e-12); // (2,0), line to (0,0)
    }

    #[test]
    fn arc_subdivision_count_large_tolerance() {
        // Large tolerance → fewer subdivisions.
        let n = arc_subdivision_count(1.0, std::f64::consts::PI, 10.0);
        assert_eq!(n, 1);
    }

    #[test]
    fn arc_subdivision_count_small_tolerance() {
        // Small tolerance → more subdivisions.
        let n = arc_subdivision_count(1.0, std::f64::consts::PI, 0.001);
        assert!(n > 10, "expected many subdivisions, got {n}");
    }
}
