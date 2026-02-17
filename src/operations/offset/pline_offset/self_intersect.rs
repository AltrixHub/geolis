use crate::geometry::pline::Pline;
use crate::math::arc_2d::arc_from_bulge;
use crate::math::intersect_2d::{
    arc_arc_intersect_2d, line_arc_intersect_2d, segment_segment_intersect_2d,
};
use crate::math::{Point3, TOLERANCE};

/// A self-intersection between two segments of a polyline.
#[derive(Debug, Clone)]
pub struct Intersection {
    /// Index of the first segment.
    pub seg_i: usize,
    /// Index of the second segment (always > `seg_i`).
    pub seg_j: usize,
    /// Parameter on segment i (0..1).
    pub t_i: f64,
    /// Parameter on segment j (0..1).
    pub t_j: f64,
    /// Intersection point.
    pub point: (f64, f64),
}

/// Finds all self-intersections between non-adjacent segments of a polyline.
///
/// Handles line-line, line-arc, arc-line, and arc-arc intersections.
/// Skips endpoint-to-endpoint touches (both parameters near 0 or 1).
#[must_use]
pub fn find_all(pline: &Pline) -> Vec<Intersection> {
    let n = pline.vertices.len();
    let seg_count = pline.segment_count();
    if seg_count < 3 {
        return Vec::new();
    }

    let eps = TOLERANCE * 100.0;
    let mut results = Vec::new();

    for i in 0..seg_count {
        let i_next = (i + 1) % n;

        for j in (i + 2)..seg_count {
            // Skip adjacent segments.
            if pline.closed && i == 0 && j == seg_count - 1 {
                continue;
            }

            let j_next = (j + 1) % n;
            let vi0 = &pline.vertices[i];
            let vi1 = &pline.vertices[i_next];
            let vj0 = &pline.vertices[j];
            let vj1 = &pline.vertices[j_next];

            let i_is_arc = vi0.bulge.abs() >= 1e-12;
            let j_is_arc = vj0.bulge.abs() >= 1e-12;

            let hits: Vec<((f64, f64), f64, f64)> = match (i_is_arc, j_is_arc) {
                (false, false) => {
                    // Line-line.
                    let a0 = Point3::new(vi0.x, vi0.y, 0.0);
                    let a1 = Point3::new(vi1.x, vi1.y, 0.0);
                    let b0 = Point3::new(vj0.x, vj0.y, 0.0);
                    let b1 = Point3::new(vj1.x, vj1.y, 0.0);
                    segment_segment_intersect_2d(&a0, &a1, &b0, &b1)
                        .map(|(pt, t, u)| vec![((pt.x, pt.y), t, u)])
                        .unwrap_or_default()
                }
                (false, true) => {
                    // Line-arc: segment i is line, segment j is arc.
                    let (cx, cy, r, sa, sw) =
                        arc_from_bulge(vj0.x, vj0.y, vj1.x, vj1.y, vj0.bulge);
                    line_arc_intersect_2d(vi0.x, vi0.y, vi1.x, vi1.y, cx, cy, r, sa, sw)
                }
                (true, false) => {
                    // Arc-line: segment i is arc, segment j is line.
                    let (cx, cy, r, sa, sw) =
                        arc_from_bulge(vi0.x, vi0.y, vi1.x, vi1.y, vi0.bulge);
                    // line_arc returns (point, t_line, t_arc); we need (point, t_arc, t_line).
                    line_arc_intersect_2d(vj0.x, vj0.y, vj1.x, vj1.y, cx, cy, r, sa, sw)
                        .into_iter()
                        .map(|(pt, t_line, t_arc)| (pt, t_arc, t_line))
                        .collect()
                }
                (true, true) => {
                    // Arc-arc.
                    let (c1x, c1y, r1, s1, sw1) =
                        arc_from_bulge(vi0.x, vi0.y, vi1.x, vi1.y, vi0.bulge);
                    let (c2x, c2y, r2, s2, sw2) =
                        arc_from_bulge(vj0.x, vj0.y, vj1.x, vj1.y, vj0.bulge);
                    arc_arc_intersect_2d(c1x, c1y, r1, s1, sw1, c2x, c2y, r2, s2, sw2)
                }
            };

            for (pt, t, u) in hits {
                // Skip vertex touches: any intersection where either parameter
                // is at a segment endpoint is a vertex-on-segment touch, not a
                // genuine crossing.
                let t_at_end = t < eps || t > 1.0 - eps;
                let u_at_end = u < eps || u > 1.0 - eps;
                if t_at_end || u_at_end {
                    continue;
                }

                results.push(Intersection {
                    seg_i: i,
                    seg_j: j,
                    t_i: t,
                    t_j: u,
                    point: pt,
                });
            }
        }
    }

    // Sort by segment index, then by parameter.
    results.sort_by(|a, b| {
        a.seg_i
            .cmp(&b.seg_i)
            .then(a.t_i.partial_cmp(&b.t_i).unwrap_or(std::cmp::Ordering::Equal))
    });

    results
}
