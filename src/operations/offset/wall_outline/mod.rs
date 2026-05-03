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
    /// # Output guarantee
    ///
    /// Each returned [`Pline`] is **closed** and the boundary set is the
    /// boolean-union outline of the stroke-expanded inputs (every output
    /// edge separates filled material from empty). The output is safe
    /// for direct ingestion by `spade::cdt`-based tessellation; this is
    /// guaranteed at the [`polygon_union::union_all_with_holes`] layer
    /// by an arrangement-based construction (split → snap → bilateral
    /// half-edge classification → polar-angle face-walk) and verified
    /// by a `#[cfg(debug_assertions)]` CDT-safe post-condition.
    ///
    /// Self-intersecting offset boundaries — which can arise from sharp
    /// zigzag or self-crossing centerlines — are flattened by the same
    /// arrangement: any internal seam edges (filled-on-both sides) are
    /// dropped during half-edge classification.
    ///
    /// # Winding
    ///
    /// Outer boundaries are returned CCW (signed area > 0); hole
    /// boundaries (when nested inside an outer) are returned CW
    /// (signed area < 0). For depth ≥ 2 nested-island outputs, parity
    /// alternates: even depth = CCW (outer), odd depth = CW (hole).
    ///
    /// # Errors
    ///
    /// - `OperationError::InvalidInput` — no polyline has at least 2 vertices.
    /// - `OperationError::Failed` — no outline can be generated, or
    ///   `polygon_union` returned ambiguous half-edge classifications
    ///   that survived the ε-shrink retries (typically degenerate input
    ///   where multiple inputs share a tangent boundary at a sampled
    ///   edge midpoint).
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
            let pwh =
                stroke::stroke_expand(&verts, pline.closed, self.left_width, self.right_width);
            if pwh.outer.len() >= 3 {
                wall_polys.push(pwh);
            }
        }

        if wall_polys.is_empty() {
            return Err(OperationError::Failed("no valid wall polygons".to_owned()).into());
        }

        // Step 2: Union all wall polygons.
        let union_result = polygon_union::union_all_with_holes(&wall_polys)?;

        if union_result.boundaries.is_empty() {
            return Err(OperationError::Failed(
                "wall outline union produced no results".to_owned(),
            )
            .into());
        }

        // Step 3: Convert union output to closed Plines. The
        // arrangement-based union already guarantees CDT safety; a
        // `#[cfg(debug_assertions)]` post-condition inside
        // `polygon_union::union_all_with_holes` re-verifies the spade
        // dry-run for defense in depth.
        let outlines: Vec<Pline> = union_result
            .boundaries
            .into_iter()
            .filter(|b| b.len() >= 3)
            .map(|b| Pline {
                vertices: b
                    .into_iter()
                    .map(|(x, y)| PlineVertex::line(x, y))
                    .collect(),
                closed: true,
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

    fn max_dist_to_centerlines(
        boundaries: &[Pline],
        centerlines: &[((f64, f64), (f64, f64))],
    ) -> f64 {
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

    fn pline_to_centerlines(pline: &Pline) -> Vec<((f64, f64), (f64, f64))> {
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
        assert!(
            (area.abs() - expected).abs() < 0.5,
            "area={area}, expected≈{expected}"
        );
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
            eprintln!("boundary[{i}] verts={} area_sign={:+}", b.vertices.len(), {
                let n = b.vertices.len();
                let mut a = 0.0;
                for k in 0..n {
                    let j = (k + 1) % n;
                    a += b.vertices[k].x * b.vertices[j].y - b.vertices[j].x * b.vertices[k].y;
                }
                if a > 0.0 {
                    1
                } else {
                    -1
                }
            });
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
            Pline::from_points(
                &[
                    Point3::new(-3.217, -4.144, 0.0),
                    Point3::new(-2.635, 2.085, 0.0),
                ],
                false,
            ),
            Pline::from_points(
                &[
                    Point3::new(-3.217, -4.144, 0.0),
                    Point3::new(2.002, -4.631, 0.0),
                ],
                false,
            ),
            Pline::from_points(
                &[
                    Point3::new(-2.635, 2.085, 0.0),
                    Point3::new(2.578, 1.534, 0.0),
                ],
                false,
            ),
            Pline::from_points(
                &[
                    Point3::new(2.002, -4.631, 0.0),
                    Point3::new(2.578, 1.534, 0.0),
                ],
                false,
            ),
            Pline::from_points(
                &[
                    Point3::new(2.002, -4.631, 0.0),
                    Point3::new(6.473, -5.049, 0.0),
                ],
                false,
            ),
            Pline::from_points(
                &[
                    Point3::new(2.578, 1.534, 0.0),
                    Point3::new(6.861, -0.896, 0.0),
                ],
                false,
            ),
            Pline::from_points(
                &[
                    Point3::new(6.473, -5.049, 0.0),
                    Point3::new(6.861, -0.896, 0.0),
                ],
                false,
            ),
        ];
        let result = run_outline(plines, d);
        assert!(!result.is_empty());
        for b in &result {
            for v in &b.vertices {
                assert!(
                    v.x >= -4.0 && v.x <= 8.0 && v.y >= -6.0 && v.y <= 3.0,
                    "vertex ({:.3}, {:.3}) out of range",
                    v.x,
                    v.y,
                );
            }
        }
    }

    /// 11-vertex comb pattern with many self-intersections. Previously
    /// failed before the `union_all_with_holes` fix that now detects
    /// same-ring self-crossings; passes after that fix.
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

    /// Integration: a centerline that crosses itself must yield only
    /// simple boundaries. The arrangement-based `polygon_union` drops
    /// internal seam edges (filled-on-both sides) during half-edge
    /// classification, so output boundaries are simple by construction
    /// and safe for downstream `spade::cdt`-based tessellation.
    #[test]
    fn closed_self_intersecting_centerline_returns_simple_boundaries() {
        // Figure-8 closed centerline.
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(2.0, 2.0),
                PlineVertex::line(0.0, 2.0),
                PlineVertex::line(2.0, 0.0),
            ],
            closed: true,
        };
        let result = run_outline(vec![pline], 0.1);

        assert!(
            !result.is_empty(),
            "self-intersecting centerline should still produce boundaries"
        );

        // Use the crate-rooted path — `wall_outline/tests` is at
        // `crate::operations::offset::wall_outline::tests`; its `super`
        // is `wall_outline`, not `pline`. Only `crate::geometry::pline::*`
        // resolves correctly.
        for (idx, b) in result.iter().enumerate() {
            assert!(
                b.closed,
                "boundary[{idx}] should be closed; got open with {} vertices",
                b.vertices.len()
            );
            assert!(
                b.vertices.len() >= 3,
                "boundary[{idx}] should have >=3 vertices; got {}",
                b.vertices.len()
            );
            assert!(
                crate::geometry::pline::self_intersection::find_self_intersection(b).is_none(),
                "boundary[{idx}] should be simple after split; \
                 still contains a self-intersection"
            );
        }
    }

    /// Regression for plan-13k: a continuous Wall centerline that crosses
    /// itself MULTIPLE times must not produce a CDT-unsafe output set,
    /// even when polygon_union returns nested-island depth-2 structure.
    /// Before T6-T10, this case panicked spade::cdt with a 2nd-crossing
    /// input.
    ///
    /// Success criterion: `WallOutline2D::execute` returns Ok(non-empty)
    /// AND every output boundary is intra-simple. The set-level CDT-safe
    /// guarantee is enforced inside `polygon_union::union_all_with_holes`
    /// by a `#[cfg(debug_assertions)]` post-condition that re-runs spade's
    /// constraint dry-run on the output — a successful `run_outline`
    /// result in debug builds is therefore evidence that the CDT
    /// constraints pass.
    #[test]
    fn multi_self_intersecting_centerline_returns_cdt_safe_set() {
        // 6-vertex closed centerline with 2 crossings (zigzag-cross shape).
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(4.0, 4.0),
                PlineVertex::line(1.0, 4.0),
                PlineVertex::line(4.0, 0.0),
                PlineVertex::line(3.0, 4.0),
                PlineVertex::line(0.0, 2.0),
            ],
            closed: true,
        };
        let result = run_outline(vec![pline], 0.1);
        assert!(!result.is_empty());
        for b in &result {
            assert!(
                crate::geometry::pline::self_intersection::find_self_intersection(b).is_none(),
                "boundary with {} vertices still has a transverse \
                 self-intersection after CDT-safe pass",
                b.vertices.len()
            );
        }
    }

    /// Semantic union check: a self-crossing OPEN centerline (the user's
    /// continuous-Wall click sequence in revion) must produce a merged
    /// outline whose interior contains the crossing region. For an open
    /// polyline `[(0,0)→(2,2)→(0,2)→(2,0)]` stroked at width 0.3 (each
    /// side), segment 0 crosses segment 2 at (1, 1). The merged outline
    /// must:
    /// - contain (1, 1) as Inside (the crossing point is wall material)
    /// - contain (0.5, 0.5) as Inside (centerline of segment 0)
    /// - have NO self-intersection in any output boundary
    ///
    /// Superseded by P3.1 `figure_8_open_centerline_no_internal_seam`:
    /// the centroid-Inside check below is satisfied even when an internal
    /// seam is emitted at the crossing, so this assertion does not
    /// distinguish a true boolean-union outline from the buggy state.
    #[test]
    #[ignore = "P3.1 supersedes: only checks centroid Inside, not no-internal-seam"]
    fn self_crossing_open_centerline_unions_to_merged_filled_region() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(2.0, 2.0),
                PlineVertex::line(0.0, 2.0),
                PlineVertex::line(2.0, 0.0),
            ],
            closed: false,
        };
        let result = run_outline(vec![pline], 0.3);
        assert!(!result.is_empty());

        // The crossing point (1, 1) and a band-arm interior (0.5, 0.5)
        // must both be classified as Inside the union region (i.e.
        // covered by an outer boundary and not punched out by a hole).
        let xs = [(1.0, 1.0), (0.5, 0.5)];
        for (x, y) in xs {
            let mut inside = false;
            for b in &result {
                let pts: Vec<(f64, f64)> = b.vertices.iter().map(|v| (v.x, v.y)).collect();
                let cls = polygon_union::point_in_polygon_class((x, y), &pts);
                let area = polygon_union::signed_area(&pts);
                if cls == polygon_union::PointClass::Inside {
                    if area > 0.0 {
                        inside = true; // outer boundary covers this point
                    } else {
                        inside = false; // hole punches it out → seam still present
                        break;
                    }
                }
            }
            assert!(
                inside,
                "({x}, {y}) should be Inside the merged wall material; \
                 internal seam still present?"
            );
        }
    }

    /// Two crossing wall polylines (X shape) must merge into a single
    /// connected wall region. The crossing-center diamond should be
    /// covered by wall material — not represented as a hole that punches
    /// it out.
    ///
    /// Superseded by P3.1 `two_crossing_separate_walls_no_internal_seam`:
    /// the centroid-Inside check is satisfied even when internal seams
    /// are emitted, so this assertion does not catch the actual bug.
    #[test]
    #[ignore = "P3.1 supersedes: only checks centroid Inside, not no-internal-seam"]
    fn two_crossing_walls_merge_with_filled_center() {
        let h = Pline {
            vertices: vec![PlineVertex::line(0.0, 1.0), PlineVertex::line(2.0, 1.0)],
            closed: false,
        };
        let v = Pline {
            vertices: vec![PlineVertex::line(1.0, 0.0), PlineVertex::line(1.0, 2.0)],
            closed: false,
        };
        let result = run_outline(vec![h, v], 0.2);
        assert!(!result.is_empty());

        // The crossing center (1, 1) must be Inside the wall material.
        let center = (1.0, 1.0);
        let mut covered_by_outer = false;
        for b in &result {
            let pts: Vec<(f64, f64)> = b.vertices.iter().map(|v| (v.x, v.y)).collect();
            let cls = polygon_union::point_in_polygon_class(center, &pts);
            let area = polygon_union::signed_area(&pts);
            if cls == polygon_union::PointClass::Inside {
                if area > 0.0 {
                    covered_by_outer = true;
                } else {
                    panic!("(1, 1) is Inside a CW boundary (hole) — center is punched out");
                }
            }
        }
        assert!(
            covered_by_outer,
            "X-crossing center (1, 1) should be wall material, but no outer boundary covers it"
        );
    }

    // -----------------------------------------------------------------
    // P3.1 — Sample-based "outline only" oracle + fixtures
    //
    // Phase 3 contract: WallOutline2D::execute output edges must be
    // exactly the boolean-union outline of the stroke-expanded inputs.
    // For every directed edge of every output boundary, "filled" material
    // must be on EXACTLY one side at a small perpendicular ε. The S1
    // bilateral-sample oracle below verifies that property; S2 verifies
    // no transverse crossing exists between or within boundaries; S3
    // re-runs spade's CDT dry-run to verify the output is constraint-safe.
    //
    // The current `union_all_with_holes` implementation fails S1 on the
    // self-crossing fixtures (figure-8, zigzag, X-cross, T-junction)
    // because it (a) drops source/direction info before tracing,
    // (b) keeps the midpoint-Inside filter at polygon_union.rs:120 which
    // treats `Boundary` as not-Inside, and (c) uses a "first matching"
    // tiebreak in trace_one_loop that picks wrong loops at degree-4
    // vertices. P3.3 replaces all three with a directed half-edge
    // arrangement; once that's in place these tests turn green.
    // -----------------------------------------------------------------

    #[derive(Debug, Clone, PartialEq, Eq)]
    enum BoundaryRef {
        Outer,
        Hole(usize),
    }

    #[derive(Debug, Clone)]
    enum FilledClass {
        Filled,
        Empty,
        AmbiguousOnBoundary {
            #[allow(
                dead_code,
                reason = "consumed via Debug formatting in panic diagnostics"
            )]
            touched: Vec<(usize, BoundaryRef)>,
        },
    }

    fn is_filled_in_input_set(
        p: (f64, f64),
        inputs: &[polygon_union::PolygonWithHoles],
    ) -> FilledClass {
        let mut touched: Vec<(usize, BoundaryRef)> = Vec::new();
        let mut any_filled = false;
        for (idx, pwh) in inputs.iter().enumerate() {
            let outer_class = polygon_union::point_in_polygon_class(p, &pwh.outer);
            if outer_class == polygon_union::PointClass::Boundary {
                touched.push((idx, BoundaryRef::Outer));
            }
            let hole_classes: Vec<polygon_union::PointClass> = pwh
                .holes
                .iter()
                .map(|h| polygon_union::point_in_polygon_class(p, h))
                .collect();
            for (hi, hc) in hole_classes.iter().enumerate() {
                if *hc == polygon_union::PointClass::Boundary {
                    touched.push((idx, BoundaryRef::Hole(hi)));
                }
            }
            if outer_class == polygon_union::PointClass::Inside
                && hole_classes
                    .iter()
                    .all(|c| *c == polygon_union::PointClass::Outside)
            {
                any_filled = true;
            }
        }
        if !touched.is_empty() {
            FilledClass::AmbiguousOnBoundary { touched }
        } else if any_filled {
            FilledClass::Filled
        } else {
            FilledClass::Empty
        }
    }

    fn build_oracle_inputs(
        plines: &[Pline],
        half_width: f64,
    ) -> Vec<polygon_union::PolygonWithHoles> {
        plines
            .iter()
            .filter_map(|p| {
                let verts: Vec<(f64, f64)> = p.vertices.iter().map(|v| (v.x, v.y)).collect();
                let pwh = super::stroke::stroke_expand(&verts, p.closed, half_width, half_width);
                if pwh.outer.len() >= 3 {
                    Some(pwh)
                } else {
                    None
                }
            })
            .collect()
    }

    /// S1: bilateral perpendicular sampling on every directed output edge.
    ///
    /// For each edge, sample at multiple positions along the edge
    /// (t = 0.5, 0.25, 0.75, 0.1, 0.9) with adaptive ε halving (3
    /// retries each). The first position that produces an unambiguous
    /// bilateral verdict is used. Multiple sample positions handle the
    /// case where a polygon edge's midpoint happens to lie on a tangent
    /// of an unrelated input boundary (e.g. when two adjacent input
    /// rectangles share a vertical boundary that the polygon edge
    /// crosses perpendicularly at its midpoint).
    ///
    /// Asserts EXACTLY one side is `Filled` (the other `Empty`).
    /// `Filled on both` = internal seam; `Empty on both` = spurious loop.
    fn assert_s1_bilateral_outline_only(
        outputs: &[Pline],
        inputs: &[polygon_union::PolygonWithHoles],
        half_width: f64,
        fixture: &str,
    ) {
        let try_sample_at = |t: f64,
                             v0: (f64, f64),
                             v1: (f64, f64),
                             nx: f64,
                             ny: f64,
                             initial_eps: f64|
         -> Option<(FilledClass, FilledClass)> {
            let sx = v0.0 + t * (v1.0 - v0.0);
            let sy = v0.1 + t * (v1.1 - v0.1);
            let mut eps = initial_eps;
            for _ in 0..4 {
                let l = is_filled_in_input_set((sx + eps * nx, sy + eps * ny), inputs);
                let r = is_filled_in_input_set((sx - eps * nx, sy - eps * ny), inputs);
                let l_amb = matches!(l, FilledClass::AmbiguousOnBoundary { .. });
                let r_amb = matches!(r, FilledClass::AmbiguousOnBoundary { .. });
                if !l_amb && !r_amb {
                    return Some((l, r));
                }
                eps *= 0.5;
                if eps < polygon_union::WALL_EPS * 0.5 {
                    break;
                }
            }
            None
        };

        for (bi, b) in outputs.iter().enumerate() {
            let n = b.vertices.len();
            if n < 2 {
                continue;
            }
            let seg_count = if b.closed { n } else { n.saturating_sub(1) };
            for i in 0..seg_count {
                let v0 = (b.vertices[i].x, b.vertices[i].y);
                let v1 = (b.vertices[(i + 1) % n].x, b.vertices[(i + 1) % n].y);
                let dx = v1.0 - v0.0;
                let dy = v1.1 - v0.1;
                let edge_len = (dx * dx + dy * dy).sqrt();
                if edge_len < polygon_union::WALL_EPS {
                    continue;
                }
                let nx = -dy / edge_len;
                let ny = dx / edge_len;
                let initial_eps = (polygon_union::WALL_EPS * 10.0)
                    .min(edge_len * 0.1)
                    .min(half_width * 0.1);

                let mut resolved: Option<(FilledClass, FilledClass)> = None;
                for &t in &[0.5_f64, 0.25, 0.75, 0.1, 0.9] {
                    if let Some(pair) = try_sample_at(t, v0, v1, nx, ny, initial_eps) {
                        resolved = Some(pair);
                        break;
                    }
                }
                let (left, right) = resolved.unwrap_or_else(|| {
                    let mid = ((v0.0 + v1.0) * 0.5, (v0.1 + v1.1) * 0.5);
                    panic!(
                        "[{fixture}] P3.1 S1: ambiguous bilateral sample at \
                         boundary={bi} edge={i} mid=({:.6},{:.6}); \
                         all sample positions exhausted (t=0.5/0.25/0.75/0.1/0.9)",
                        mid.0, mid.1
                    )
                });
                let l_filled = matches!(left, FilledClass::Filled);
                let r_filled = matches!(right, FilledClass::Filled);
                assert!(
                    l_filled != r_filled,
                    "[{fixture}] P3.1 S1: edge at boundary={bi} edge={i} \
                     v0=({:.6},{:.6}) v1=({:.6},{:.6}) \
                     has filled=left:{} right:{} (must be exactly one filled)",
                    v0.0,
                    v0.1,
                    v1.0,
                    v1.1,
                    l_filled,
                    r_filled
                );
            }
        }
    }

    /// S2: no transverse crossing within or between output boundaries.
    fn assert_s2_no_transverse_crossings(outputs: &[Pline], fixture: &str) {
        use crate::geometry::pline::self_intersection::{
            find_self_intersection, segment_segment_intersection_2d,
        };
        for (bi, b) in outputs.iter().enumerate() {
            if let Some((i, j, x, y)) = find_self_intersection(b) {
                panic!(
                    "[{fixture}] P3.1 S2: boundary {bi} self-intersects at \
                     edges {i}-{j} ({x:.4}, {y:.4})"
                );
            }
        }
        for (ai, a) in outputs.iter().enumerate() {
            let an = a.vertices.len();
            let aseg = if a.closed { an } else { an.saturating_sub(1) };
            for (bi, b) in outputs.iter().enumerate().skip(ai + 1) {
                let bn = b.vertices.len();
                let bseg = if b.closed { bn } else { bn.saturating_sub(1) };
                for ai_e in 0..aseg {
                    let p0 = (a.vertices[ai_e].x, a.vertices[ai_e].y);
                    let p1 = (a.vertices[(ai_e + 1) % an].x, a.vertices[(ai_e + 1) % an].y);
                    for bj_e in 0..bseg {
                        let q0 = (b.vertices[bj_e].x, b.vertices[bj_e].y);
                        let q1 = (b.vertices[(bj_e + 1) % bn].x, b.vertices[(bj_e + 1) % bn].y);
                        if let Some((t, _u)) = segment_segment_intersection_2d(p0, p1, q0, q1) {
                            let x = p0.0 + t * (p1.0 - p0.0);
                            let y = p0.1 + t * (p1.1 - p0.1);
                            panic!(
                                "[{fixture}] P3.1 S2: boundary {ai} edge {ai_e} \
                                 crosses boundary {bi} edge {bj_e} at ({x:.4}, {y:.4})"
                            );
                        }
                    }
                }
            }
        }
    }

    /// S3: spade CDT dry-run accepts every output edge as a constraint.
    fn assert_s3_cdt_dry_run(outputs: &[Pline], fixture: &str) {
        use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};
        let mut cdt: ConstrainedDelaunayTriangulation<Point2<f64>> =
            ConstrainedDelaunayTriangulation::new();
        for (bi, b) in outputs.iter().enumerate() {
            let n = b.vertices.len();
            if n < 3 {
                continue;
            }
            let mut handles = Vec::with_capacity(n);
            for (vi, v) in b.vertices.iter().enumerate() {
                match cdt.insert(Point2::new(v.x, v.y)) {
                    Ok(h) => handles.push(h),
                    Err(e) => panic!(
                        "[{fixture}] P3.1 S3: CDT vertex insert rejected (b={bi}, v={vi}): {e:?}"
                    ),
                }
            }
            for k in 0..n {
                let from = handles[k];
                let to = handles[(k + 1) % n];
                if from == to {
                    continue;
                }
                let added = cdt.try_add_constraint(from, to);
                assert!(
                    !added.is_empty(),
                    "[{fixture}] P3.1 S3: CDT constraint rejected (b={bi}, edge={k})"
                );
            }
        }
    }

    fn run_p3_oracle(plines: Vec<Pline>, half_width: f64, fixture: &str) -> Vec<Pline> {
        let inputs = build_oracle_inputs(&plines, half_width);
        let outputs = WallOutline2D::new(plines, half_width)
            .execute()
            .unwrap_or_else(|e| panic!("[{fixture}] WallOutline2D::execute failed: {e}"));
        assert!(
            !outputs.is_empty(),
            "[{fixture}] WallOutline2D produced no boundaries"
        );
        assert_s1_bilateral_outline_only(&outputs, &inputs, half_width, fixture);
        assert_s2_no_transverse_crossings(&outputs, fixture);
        assert_s3_cdt_dry_run(&outputs, fixture);
        outputs
    }

    // ---- NEW fixtures (target failures of the current bug) ----

    #[test]
    fn figure_8_open_centerline_no_internal_seam() {
        // S1' independent oracle is checked separately at the bottom.
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(2.0, 2.0),
                PlineVertex::line(0.0, 2.0),
                PlineVertex::line(2.0, 0.0),
            ],
            closed: false,
        };
        let half_width = 0.3;
        let outputs = run_p3_oracle(vec![pline.clone()], half_width, "figure_8_open");

        // S1' (analytic centerline-band oracle) — independent of stroke_expand.
        //
        // For figure-8 the centerline crossing at (1, 1) creates band overlaps
        // whose perimeter passes through 4 inner corners arranged as a
        // diamond around (1, 1). Sample those crossing-region edges only.
        //
        // The band oracle is approximate near miter overshoots at convex
        // corners and small triangular hole pockets between band arms. So
        // we only assert agreement where BOTH oracles produce a clear
        // exactly-one-side-filled verdict — disagreements where the band
        // oracle says "both sides outside" mean we're in a miter-artifact
        // region the analytic model doesn't cover (skip silently).
        let centerline: Vec<(f64, f64)> = pline.vertices.iter().map(|v| (v.x, v.y)).collect();
        let inputs = build_oracle_inputs(&[pline], half_width);
        let mut sampled = 0usize;
        for b in &outputs {
            let n = b.vertices.len();
            let seg_count = if b.closed { n } else { n.saturating_sub(1) };
            for i in 0..seg_count {
                let v0 = (b.vertices[i].x, b.vertices[i].y);
                let v1 = (b.vertices[(i + 1) % n].x, b.vertices[(i + 1) % n].y);
                let mid = ((v0.0 + v1.0) * 0.5, (v0.1 + v1.1) * 0.5);
                let dx = v1.0 - v0.0;
                let dy = v1.1 - v0.1;
                let edge_len = (dx * dx + dy * dy).sqrt();
                if edge_len < polygon_union::WALL_EPS {
                    continue;
                }
                let nx = -dy / edge_len;
                let ny = dx / edge_len;
                let eps = (polygon_union::WALL_EPS * 10.0).min(half_width * 0.1);
                let lp = (mid.0 + eps * nx, mid.1 + eps * ny);
                let rp = (mid.0 - eps * nx, mid.1 - eps * ny);
                let stroke_l = matches!(is_filled_in_input_set(lp, &inputs), FilledClass::Filled);
                let stroke_r = matches!(is_filled_in_input_set(rp, &inputs), FilledClass::Filled);
                let band_l = is_within_centerline_band(lp, &centerline, half_width);
                let band_r = is_within_centerline_band(rp, &centerline, half_width);
                let stroke_unambiguous = stroke_l != stroke_r;
                let band_unambiguous = band_l != band_r;
                if stroke_unambiguous && band_unambiguous {
                    assert_eq!(
                        (stroke_l, stroke_r),
                        (band_l, band_r),
                        "S1': stroke_expand and analytic centerline-band disagree on \
                         which side is filled at edge mid=({:.4},{:.4})",
                        mid.0,
                        mid.1
                    );
                    sampled += 1;
                }
            }
        }
        assert!(
            sampled > 0,
            "S1' must agree with stroke_expand on at least one polygon edge"
        );
    }

    #[allow(
        clippy::many_single_char_names,
        reason = "p/a/b/d/n are domain-standard names for point/endpoints/distance/count \
                  in 2D segment-distance geometry"
    )]
    fn is_within_centerline_band(p: (f64, f64), vertices: &[(f64, f64)], half_width: f64) -> bool {
        let n = vertices.len();
        if n < 2 {
            return false;
        }
        for i in 0..(n - 1) {
            let a = vertices[i];
            let b = vertices[i + 1];
            let d = point_to_segment_dist(p.0, p.1, a.0, a.1, b.0, b.1);
            if d < half_width - polygon_union::WALL_EPS * 10.0 {
                return true;
            }
        }
        false
    }

    #[test]
    fn zigzag_with_inner_crossings_no_internal_seam() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(4.0, 4.0),
                PlineVertex::line(1.0, 4.0),
                PlineVertex::line(4.0, 0.0),
                PlineVertex::line(3.0, 4.0),
                PlineVertex::line(0.0, 2.0),
            ],
            closed: false,
        };
        run_p3_oracle(vec![pline], 0.1, "zigzag_with_inner_crossings");
    }

    #[test]
    fn two_crossing_separate_walls_no_internal_seam() {
        let h = Pline {
            vertices: vec![PlineVertex::line(0.0, 1.0), PlineVertex::line(2.0, 1.0)],
            closed: false,
        };
        let v = Pline {
            vertices: vec![PlineVertex::line(1.0, 0.0), PlineVertex::line(1.0, 2.0)],
            closed: false,
        };
        run_p3_oracle(vec![h, v], 0.2, "two_crossing_separate_walls");
    }

    #[test]
    fn t_junction_two_open_walls_no_internal_seam() {
        let horiz = Pline {
            vertices: vec![PlineVertex::line(0.0, 0.0), PlineVertex::line(4.0, 0.0)],
            closed: false,
        };
        let stem = Pline {
            vertices: vec![PlineVertex::line(2.0, 0.0), PlineVertex::line(2.0, 3.0)],
            closed: false,
        };
        run_p3_oracle(vec![horiz, stem], 0.15, "t_junction_two_open_walls");
    }

    // ---- EXISTING-regression fixtures (no-hole simple cases) ----

    #[test]
    fn two_overlapping_disjoint_walls_outline_only() {
        let a = Pline {
            vertices: vec![PlineVertex::line(0.0, 0.0), PlineVertex::line(3.0, 0.0)],
            closed: false,
        };
        let b = Pline {
            vertices: vec![PlineVertex::line(2.0, 0.0), PlineVertex::line(5.0, 0.0)],
            closed: false,
        };
        run_p3_oracle(vec![a, b], 0.15, "two_overlapping_disjoint_walls");
    }

    #[test]
    fn two_walls_sharing_one_endpoint_outline_only() {
        let a = Pline {
            vertices: vec![PlineVertex::line(0.0, 0.0), PlineVertex::line(4.0, 0.0)],
            closed: false,
        };
        let b = Pline {
            vertices: vec![PlineVertex::line(4.0, 0.0), PlineVertex::line(8.0, 0.0)],
            closed: false,
        };
        run_p3_oracle(vec![a, b], 0.15, "two_walls_sharing_one_endpoint");
    }

    #[test]
    fn three_walls_t_configuration_outline_only() {
        let a = Pline {
            vertices: vec![PlineVertex::line(0.0, 0.0), PlineVertex::line(4.0, 0.0)],
            closed: false,
        };
        let b = Pline {
            vertices: vec![PlineVertex::line(4.0, 0.0), PlineVertex::line(4.0, 3.0)],
            closed: false,
        };
        let c = Pline {
            vertices: vec![PlineVertex::line(4.0, 0.0), PlineVertex::line(8.0, 0.0)],
            closed: false,
        };
        run_p3_oracle(vec![a, b, c], 0.15, "three_walls_t_configuration");
    }

    #[test]
    fn closed_rectangle_ring_annular_outline_only() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(10.0, 0.0),
                PlineVertex::line(10.0, 10.0),
                PlineVertex::line(0.0, 10.0),
            ],
            closed: true,
        };
        run_p3_oracle(vec![pline], 0.3, "closed_rectangle_ring");
    }
}
