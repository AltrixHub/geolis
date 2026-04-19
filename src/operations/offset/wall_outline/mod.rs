pub(crate) mod polygon_union;
mod stroke;

use crate::error::{OperationError, Result};
use crate::geometry::pline::{Pline, PlineVertex};

/// Generates wall outlines from one or more centerline polylines.
///
/// Given a collection of `Pline`s representing wall centerlines (potentially with
/// self-intersecting paths), produces closed outline polygons at the
/// specified distances. When multiple polylines are provided,
/// their segments are merged into a single network so that intersections
/// between separate walls are properly trimmed.
///
/// The wall material spans from `left_width` to the left of each segment to
/// `right_width` to the right (using the segment's forward direction).
/// Use [`WallOutline2D::new`] for a centred wall (`left == right == half_thickness`)
/// or [`WallOutline2D::new_asymmetric`] for a wall aligned to one side of the baseline.
#[derive(Debug)]
pub struct WallOutline2D {
    plines: Vec<Pline>,
    left_width: f64,
    right_width: f64,
}

impl WallOutline2D {
    /// Creates a centred wall outline (equal offset on both sides of the baseline).
    #[must_use]
    pub fn new(plines: Vec<Pline>, half_width: f64) -> Self {
        Self {
            plines,
            left_width: half_width,
            right_width: half_width,
        }
    }

    /// Creates a wall outline with independent left and right offsets.
    ///
    /// - `left_width = 0, right_width = thickness`: baseline is the left (inner) boundary;
    ///   wall material extends entirely to the right.
    /// - `left_width = thickness, right_width = 0`: baseline is the right (outer) boundary;
    ///   wall material extends entirely to the left.
    #[must_use]
    pub fn new_asymmetric(plines: Vec<Pline>, left_width: f64, right_width: f64) -> Self {
        Self {
            plines,
            left_width,
            right_width,
        }
    }

    /// Executes the wall outline generation.
    ///
    /// # Errors
    ///
    /// Returns `OperationError::InvalidInput` if no polyline has at least
    /// 2 vertices, or `OperationError::Failed` if no outline can be generated.
    pub fn execute(&self) -> Result<Vec<Pline>> {
        let valid: Vec<&Pline> = self
            .plines
            .iter()
            .filter(|p| p.vertices.len() >= 2)
            .collect();

        if valid.is_empty() {
            return Err(OperationError::InvalidInput(
                "at least 2 vertices required for wall outline".to_owned(),
            )
            .into());
        }

        if self.left_width.abs() < crate::math::TOLERANCE
            && self.right_width.abs() < crate::math::TOLERANCE
        {
            return Ok(self.plines.clone());
        }

        // Step 1: Stroke-expand each polyline into a wall polygon.
        let mut wall_polys: Vec<polygon_union::PolygonWithHoles> = Vec::new();

        for pline in &valid {
            // Tessellate arc segments into line segments.
            // Tolerance scales with wall width for consistent arc resolution.
            let has_arcs = pline.vertices.iter().any(|v| v.bulge.abs() > 1e-12);
            let arc_tolerance = self.left_width.max(self.right_width) * 0.1;
            let mut verts: Vec<(f64, f64)> = if has_arcs {
                let pts = pline.to_points(arc_tolerance.max(polygon_union::WALL_EPS));
                pts.iter().map(|p| (p.x, p.y)).collect()
            } else {
                pline.vertices.iter().map(|v| (v.x, v.y)).collect()
            };
            // For closed polylines, to_points() may duplicate the start point at
            // the end. Strip trailing duplicate to avoid a zero-length segment.
            if pline.closed && verts.len() >= 2 {
                let first = verts[0];
                let last = verts[verts.len() - 1];
                if (first.0 - last.0).powi(2) + (first.1 - last.1).powi(2)
                    < polygon_union::WALL_EPS * polygon_union::WALL_EPS
                {
                    verts.pop();
                }
            }
            let pwh = stroke::stroke_expand(&verts, pline.closed, self.left_width, self.right_width);
            if pwh.outer.len() >= 3 {
                wall_polys.push(pwh);
            }
        }

        if wall_polys.is_empty() {
            return Err(OperationError::Failed("no valid wall polygons".to_owned()).into());
        }

        // Step 2: Union all wall polygons.
        let union_result = polygon_union::union_all_with_holes(&wall_polys);

        if union_result.boundaries.is_empty() {
            return Err(OperationError::Failed(
                "wall outline union produced no results".to_owned(),
            )
            .into());
        }

        // Step 3: Convert to Pline boundaries.
        let outlines: Vec<Pline> = union_result
            .boundaries
            .into_iter()
            .filter(|b| b.len() >= 3)
            .map(|b| {
                let vertices = b
                    .into_iter()
                    .map(|(x, y)| PlineVertex::line(x, y))
                    .collect();
                Pline {
                    vertices,
                    closed: true,
                }
            })
            .collect();

        if outlines.is_empty() {
            return Err(OperationError::Failed(
                "wall outline union produced no valid boundaries".to_owned(),
            )
            .into());
        }

        Ok(outlines)
    }
}


#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::distance_2d::point_to_segment_dist;
    use crate::math::Point3;

    fn run_outline(plines: Vec<Pline>, d: f64) -> Vec<Pline> {
        WallOutline2D::new(plines, d).execute().unwrap()
    }

    fn total_area(boundaries: &[Pline]) -> f64 {
        boundaries
            .iter()
            .map(|b| {
                let n = b.vertices.len();
                let mut a = 0.0;
                for i in 0..n {
                    let j = (i + 1) % n;
                    a += b.vertices[i].x * b.vertices[j].y;
                    a -= b.vertices[j].x * b.vertices[i].y;
                }
                a * 0.5
            })
            .sum()
    }

    fn max_dist_to_centerlines(boundaries: &[Pline], centerlines: &[((f64,f64),(f64,f64))]) -> f64 {
        let mut max_d = 0.0_f64;
        for b in boundaries {
            for v in &b.vertices {
                let d = centerlines
                    .iter()
                    .map(|&(a, b)| point_to_segment_dist(v.x, v.y, a.0, a.1, b.0, b.1))
                    .fold(f64::MAX, f64::min);
                max_d = max_d.max(d);
            }
        }
        max_d
    }

    fn pline_to_centerlines(pline: &Pline) -> Vec<((f64,f64),(f64,f64))> {
        let n = pline.vertices.len();
        let seg_count = if pline.closed { n } else { n.saturating_sub(1) };
        (0..seg_count)
            .map(|i| {
                let a = &pline.vertices[i];
                let b = &pline.vertices[(i + 1) % n];
                ((a.x, a.y), (b.x, b.y))
            })
            .collect()
    }

    #[test]
    fn single_segment() {
        let pline = Pline::from_points(
            &[Point3::new(0.0, 0.0, 0.0), Point3::new(5.0, 0.0, 0.0)],
            false,
        );
        let d = 0.3;
        let result = run_outline(vec![pline.clone()], d);
        assert!(!result.is_empty());
        let area = total_area(&result);
        let expected = 5.0 * 0.6; // length * thickness
        assert!((area.abs() - expected).abs() < 0.5, "area={area}, expected≈{expected}");
    }

    #[test]
    fn l_shape() {
        let pline = Pline::from_points(
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(3.0, 0.0, 0.0),
                Point3::new(3.0, 3.0, 0.0),
            ],
            false,
        );
        let d = 0.3;
        let result = run_outline(vec![pline.clone()], d);
        assert!(!result.is_empty());
        let cls = pline_to_centerlines(&pline);
        let max_d = max_dist_to_centerlines(&result, &cls);
        assert!(max_d < d * 3.0, "max_d={max_d}");
    }

    #[test]
    fn closed_square() {
        let pline = Pline::from_points(
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(10.0, 0.0, 0.0),
                Point3::new(10.0, 10.0, 0.0),
                Point3::new(0.0, 10.0, 0.0),
            ],
            true,
        );
        let d = 0.3;
        let result = run_outline(vec![pline], d);
        assert!(result.len() >= 2, "outer + hole, got {}", result.len());
        let area = total_area(&result).abs();
        // Wall ring area ≈ perimeter * thickness = 40 * 0.6 = 24
        assert!(area > 15.0 && area < 30.0, "area={area}");
    }

    #[test]
    fn closed_l_room() {
        let pline = Pline::from_points(
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(5.0, 0.0, 0.0),
                Point3::new(5.0, 3.0, 0.0),
                Point3::new(3.0, 3.0, 0.0),
                Point3::new(3.0, 5.0, 0.0),
                Point3::new(0.0, 5.0, 0.0),
            ],
            true,
        );
        let d = 0.3;
        let result = run_outline(vec![pline], d);
        assert!(result.len() >= 2, "outer + hole(s)");
    }

    #[test]
    fn t_junction_single_pline() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 3.0),
                PlineVertex::line(8.0, 3.0),
                PlineVertex::line(4.0, 3.0),
                PlineVertex::line(4.0, 5.0),
            ],
            closed: false,
        };
        let d = 1.0;
        let result = run_outline(vec![pline.clone()], d);
        assert!(!result.is_empty());
        let cls = pline_to_centerlines(&pline);
        let max_d = max_dist_to_centerlines(&result, &cls);
        assert!(max_d < d * 3.0, "max_d={max_d}");
    }

    #[test]
    fn t_junction_independent_plines() {
        let d = 0.3;
        let plines = vec![
            Pline::from_points(
                &[Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 0.0, 0.0)],
                false,
            ),
            Pline::from_points(
                &[Point3::new(4.0, 0.0, 0.0), Point3::new(4.0, 3.0, 0.0)],
                false,
            ),
            Pline::from_points(
                &[Point3::new(4.0, 0.0, 0.0), Point3::new(8.0, 0.0, 0.0)],
                false,
            ),
        ];
        let result = run_outline(plines, d);
        assert!(!result.is_empty());
    }

    /// Two adjacent **closed** zone footprints (same Y extent). Each
    /// footprint strokes to an annulus; the combined output must be:
    ///   - 1 outer rectangle (combined perimeter)
    ///   - 2 holes (one per room, separated by the shared wall material)
    ///
    /// This mirrors `BoundarySolver` emitting one closed-ring `WallBaseline`
    /// per zone into WallLayer's Rings slot.
    #[test]
    fn two_adjacent_zones_one_outer_two_holes() {
        let d = 0.15;
        let plines = vec![
            // Zone A footprint: (0,0) to (5,3)
            Pline::from_points(
                &[
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(5.0, 0.0, 0.0),
                    Point3::new(5.0, 3.0, 0.0),
                    Point3::new(0.0, 3.0, 0.0),
                ],
                true,
            ),
            // Zone B footprint: (5,0) to (8,3)
            Pline::from_points(
                &[
                    Point3::new(5.0, 0.0, 0.0),
                    Point3::new(8.0, 0.0, 0.0),
                    Point3::new(8.0, 3.0, 0.0),
                    Point3::new(5.0, 3.0, 0.0),
                ],
                true,
            ),
        ];
        let result = run_outline(plines, d);
        let outer_count = result
            .iter()
            .filter(|p| {
                let pts: Vec<Point3> = p
                    .vertices
                    .iter()
                    .map(|v| Point3::new(v.x, v.y, 0.0))
                    .collect();
                let mut area = 0.0;
                let n = pts.len();
                for i in 0..n {
                    let j = (i + 1) % n;
                    area += pts[i].x * pts[j].y - pts[j].x * pts[i].y;
                }
                area > 0.0
            })
            .count();
        let hole_count = result.len() - outer_count;
        assert_eq!(outer_count, 1, "two adjacent zones: one combined outer");
        assert_eq!(hole_count, 2, "two adjacent zones: two separate rooms");

        // Dump the outer boundary's vertices to stderr for diagnosis. The
        // combined perimeter is geometrically a 4-corner rectangle —
        // polygon_union may leave extra colinear split vertices, but the
        // crease filter in WallLayer must drop those from the 3D wireframe.
        for (i, b) in result.iter().enumerate() {
            eprintln!(
                "boundary[{i}] verts={} area_sign={:+}",
                b.vertices.len(),
                {
                    let n = b.vertices.len();
                    let mut a = 0.0;
                    for k in 0..n {
                        let j = (k + 1) % n;
                        a += b.vertices[k].x * b.vertices[j].y
                            - b.vertices[j].x * b.vertices[k].y;
                    }
                    if a > 0.0 { 1 } else { -1 }
                }
            );
            for (k, v) in b.vertices.iter().enumerate() {
                eprintln!("  v[{k}] = ({:.3}, {:.3})", v.x, v.y);
            }
        }
    }

    /// Two open 2-vertex walls: one horizontal through (0,0)-(4,0), one
    /// vertical stem at (2,0)-(2,3). They form a T; the full pipeline
    /// must return exactly one outline boundary.
    #[test]
    fn two_open_walls_forming_t_merge_into_one_boundary() {
        let d = 0.15;
        let plines = vec![
            Pline::from_points(
                &[Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 0.0, 0.0)],
                false,
            ),
            Pline::from_points(
                &[Point3::new(2.0, 0.0, 0.0), Point3::new(2.0, 3.0, 0.0)],
                false,
            ),
        ];
        let result = run_outline(plines, d);
        assert_eq!(
            result.len(),
            1,
            "two overlapping stroke rectangles must merge into a single T boundary, got {}",
            result.len()
        );
    }

    #[test]
    fn two_adjacent_rectangles() {
        let d = 0.15;
        let plines = vec![
            Pline::from_points(
                &[
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(4.0, 0.0, 0.0),
                    Point3::new(4.0, 3.0, 0.0),
                    Point3::new(0.0, 3.0, 0.0),
                ],
                true,
            ),
            Pline::from_points(
                &[
                    Point3::new(4.0, 0.0, 0.0),
                    Point3::new(8.0, 0.0, 0.0),
                    Point3::new(8.0, 3.0, 0.0),
                    Point3::new(4.0, 3.0, 0.0),
                ],
                true,
            ),
        ];
        let result = run_outline(plines, d);
        assert!(!result.is_empty());
    }

    #[test]
    fn angled_per_segment_walls() {
        let d = 0.15;
        let plines = vec![
            Pline::from_points(&[Point3::new(-3.217, -4.144, 0.0), Point3::new(-2.635, 2.085, 0.0)], false),
            Pline::from_points(&[Point3::new(-3.217, -4.144, 0.0), Point3::new(2.002, -4.631, 0.0)], false),
            Pline::from_points(&[Point3::new(-2.635, 2.085, 0.0), Point3::new(2.578, 1.534, 0.0)], false),
            Pline::from_points(&[Point3::new(2.002, -4.631, 0.0), Point3::new(2.578, 1.534, 0.0)], false),
            Pline::from_points(&[Point3::new(2.002, -4.631, 0.0), Point3::new(6.473, -5.049, 0.0)], false),
            Pline::from_points(&[Point3::new(2.578, 1.534, 0.0), Point3::new(6.861, -0.896, 0.0)], false),
            Pline::from_points(&[Point3::new(6.473, -5.049, 0.0), Point3::new(6.861, -0.896, 0.0)], false),
        ];
        let result = run_outline(plines, d);
        assert!(!result.is_empty());
        for b in &result {
            for v in &b.vertices {
                assert!(
                    v.x >= -4.0 && v.x <= 8.0 && v.y >= -6.0 && v.y <= 3.0,
                    "vertex ({:.3}, {:.3}) out of range",
                    v.x, v.y,
                );
            }
        }
    }

    #[test]
    fn double_cross() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(3.0, 0.0),
                PlineVertex::line(3.0, 10.0),
                PlineVertex::line(3.0, 7.0),
                PlineVertex::line(0.0, 7.0),
                PlineVertex::line(10.0, 7.0),
                PlineVertex::line(7.0, 7.0),
                PlineVertex::line(7.0, 10.0),
                PlineVertex::line(7.0, 0.0),
                PlineVertex::line(7.0, 3.0),
                PlineVertex::line(10.0, 3.0),
                PlineVertex::line(0.0, 3.0),
            ],
            closed: false,
        };
        let d = 0.3;
        let result = run_outline(vec![pline], d);
        assert!(!result.is_empty());
    }

    #[test]
    fn cross_shape() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 3.0),
                PlineVertex::line(10.0, 3.0),
                PlineVertex::line(5.0, 3.0),
                PlineVertex::line(5.0, 0.0),
                PlineVertex::line(5.0, 10.0),
            ],
            closed: false,
        };
        let d = 0.5;
        let result = run_outline(vec![pline], d);
        assert!(!result.is_empty());
    }
}
