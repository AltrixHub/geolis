use crate::math::arc_2d::{arc_from_bulge, arc_point_at};
use crate::math::Point3;

/// Self-intersection detection primitives. `find_self_intersection` is
/// reused by the `WallOutline2D` test oracle (P3.1 S2) and by the
/// figure-8 / multi-self-crossing fixture assertions; consequently the
/// module is test-only — `polygon_union` no longer relies on it for
/// production output.
#[cfg(test)]
pub(crate) mod self_intersection;

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
        let vertices = points.iter().map(|p| PlineVertex::line(p.x, p.y)).collect();
        Self { vertices, closed }
    }

    /// Converts this polyline to a list of `Point3` by tessellating arcs into line segments.
    ///
    /// `tolerance` controls the maximum deviation between the arc and its chord approximation.
    #[must_use]
    pub fn to_points(&self, tolerance: f64) -> Vec<Point3> {
        self.to_points_with_sources(tolerance).0
    }

    /// [`Self::to_points`] variant that additionally reports, for each
    /// tessellated edge, the index of the source polyline segment it
    /// approximates.
    ///
    /// Returns `(points, sources)` where `sources[j]` is the source
    /// segment index of the edge `points[j] → points[j + 1]`
    /// (`sources.len() == points.len() - 1` whenever at least one
    /// segment exists). A line segment contributes one edge; an arc
    /// segment contributes one edge per tessellation chord, all mapped
    /// to the same source segment index.
    #[must_use]
    pub fn to_points_with_sources(&self, tolerance: f64) -> (Vec<Point3>, Vec<usize>) {
        let n = self.vertices.len();
        if n == 0 {
            return (Vec::new(), Vec::new());
        }

        let seg_count = self.segment_count();
        let mut points = Vec::with_capacity(n * 2);
        let mut sources = Vec::with_capacity(n * 2);

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
                sources.push(i);
            } else {
                // Arc: tessellate into line segments.
                let (cx, cy, radius, start_angle, sweep) =
                    arc_from_bulge(v0.x, v0.y, v1.x, v1.y, v0.bulge);

                if radius < 1e-12 {
                    points.push(Point3::new(v1.x, v1.y, 0.0));
                    sources.push(i);
                    continue;
                }

                // Number of subdivisions based on tolerance.
                let n_sub = arc_subdivision_count(radius, sweep.abs(), tolerance);

                for j in 1..n_sub {
                    let t = f64::from(j) / f64::from(n_sub);
                    let (px, py) = arc_point_at(cx, cy, radius, start_angle, sweep, t);
                    points.push(Point3::new(px, py, 0.0));
                    sources.push(i);
                }
                points.push(Point3::new(v1.x, v1.y, 0.0));
                sources.push(i);
            }
        }

        (points, sources)
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

    /// Returns the signed area enclosed by this polyline.
    ///
    /// Counter-clockwise orientation yields a positive area, clockwise a
    /// negative one. Circular-arc segments are accounted for exactly:
    /// each bulged segment contributes its chord to the shoelace sum plus
    /// the signed circular-segment area between chord and arc,
    /// `sign(bulge) · r²/2 · (θ − sin θ)` with `θ = 4·atan(|bulge|)`.
    ///
    /// An open polyline is treated as implicitly closed by a straight
    /// chord from the last vertex back to the first (shoelace
    /// convention); the last vertex's bulge is ignored because it has no
    /// segment. Fewer than two vertices → `0.0`.
    #[must_use]
    pub fn signed_area(&self) -> f64 {
        let n = self.vertices.len();
        if n < 2 {
            return 0.0;
        }

        // Chord shoelace over the (implicitly) closed ring.
        let mut twice_chord_area = 0.0;
        for i in 0..n {
            let v0 = &self.vertices[i];
            let v1 = &self.vertices[(i + 1) % n];
            twice_chord_area += v0.x * v1.y - v1.x * v0.y;
        }
        let mut area = 0.5 * twice_chord_area;

        // Circular-segment corrections on real segments only: for an
        // open polyline the implicit closing chord stays straight.
        for i in 0..self.segment_count() {
            let v0 = &self.vertices[i];
            let v1 = &self.vertices[(i + 1) % n];
            if v0.bulge.abs() < 1e-12 {
                continue;
            }
            let (_, _, radius, _, sweep) = arc_from_bulge(v0.x, v0.y, v1.x, v1.y, v0.bulge);
            if radius < 1e-12 {
                // Degenerate chord: no well-defined arc.
                continue;
            }
            let theta = sweep.abs();
            area += v0.bulge.signum() * 0.5 * radius * radius * (theta - theta.sin());
        }

        area
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
        assert!(
            pts.len() > 2,
            "expected more than 2 points, got {}",
            pts.len()
        );
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
    fn signed_area_unit_square_ccw() {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let pline = Pline::from_points(&pts, true);
        let area = pline.signed_area();
        assert!((area - 1.0).abs() < 1e-12, "area={area}");
    }

    #[test]
    fn signed_area_cw_rect_negative() {
        // 2x3 rectangle traversed clockwise → -6.
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 3.0, 0.0),
            Point3::new(2.0, 3.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ];
        let pline = Pline::from_points(&pts, true);
        let area = pline.signed_area();
        assert!((area + 6.0).abs() < 1e-12, "area={area}");
    }

    #[test]
    fn signed_area_bulged_edge_adds_circular_segment() {
        // 4x3 CCW rectangle with bulge=0.5 on the bottom edge (0,0)→(4,0).
        // Positive bulge bows outward (below the chord), so the circular
        // segment adds to the chord-polygon area.
        let bulge = 0.5;
        let pline = Pline {
            vertices: vec![
                PlineVertex::new(0.0, 0.0, bulge),
                PlineVertex::line(4.0, 0.0),
                PlineVertex::line(4.0, 3.0),
                PlineVertex::line(0.0, 3.0),
            ],
            closed: true,
        };
        // Closed form from the same arc conversion: r²/2 · (θ − sin θ).
        let (_, _, r, _, sweep) = arc_from_bulge(0.0, 0.0, 4.0, 0.0, bulge);
        let theta = sweep.abs();
        let expected = 12.0 + 0.5 * r * r * (theta - theta.sin());
        let area = pline.signed_area();
        assert!(
            (area - expected).abs() < 1e-10,
            "area={area} expected={expected}"
        );
    }

    #[test]
    fn signed_area_negative_bulge_subtracts() {
        // Same rectangle but the bottom edge bows inward (bulge=-0.5):
        // the circular segment is subtracted from the chord-polygon area.
        let bulge = -0.5;
        let pline = Pline {
            vertices: vec![
                PlineVertex::new(0.0, 0.0, bulge),
                PlineVertex::line(4.0, 0.0),
                PlineVertex::line(4.0, 3.0),
                PlineVertex::line(0.0, 3.0),
            ],
            closed: true,
        };
        let (_, _, r, _, sweep) = arc_from_bulge(0.0, 0.0, 4.0, 0.0, bulge);
        let theta = sweep.abs();
        let expected = 12.0 - 0.5 * r * r * (theta - theta.sin());
        let area = pline.signed_area();
        assert!(
            area < 12.0,
            "area={area} should be less than the straight rect"
        );
        assert!(
            (area - expected).abs() < 1e-10,
            "area={area} expected={expected}"
        );
    }

    #[test]
    fn signed_area_two_vertex_circle() {
        // Two semicircles (bulge=1 each) form a full CCW circle of
        // radius 1 (half the chord distance) → area = π·r².
        let pline = Pline {
            vertices: vec![
                PlineVertex::new(0.0, 0.0, 1.0),
                PlineVertex::new(2.0, 0.0, 1.0),
            ],
            closed: true,
        };
        let area = pline.signed_area();
        assert!(
            (area - std::f64::consts::PI).abs() < 1e-10,
            "area={area} expected=π"
        );
    }

    #[test]
    fn signed_area_reversal_negates() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::new(0.0, 0.0, 0.5),
                PlineVertex::line(4.0, 0.0),
                PlineVertex::line(4.0, 3.0),
                PlineVertex::line(0.0, 3.0),
            ],
            closed: true,
        };
        let area = pline.signed_area();
        let rev_area = pline.reversed().signed_area();
        assert!(
            (rev_area + area).abs() < 1e-10,
            "area={area} rev_area={rev_area}"
        );
    }

    #[test]
    fn signed_area_open_pline_uses_implicit_closing_chord() {
        // Open polylines are treated as implicitly closed by a straight
        // chord from the last vertex back to the first.
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let pline = Pline::from_points(&pts, false);
        let area = pline.signed_area();
        assert!((area - 1.0).abs() < 1e-12, "area={area}");
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
