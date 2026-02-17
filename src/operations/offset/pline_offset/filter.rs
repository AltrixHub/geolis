use crate::geometry::pline::Pline;
use crate::math::arc_2d::arc_from_bulge;
use crate::math::distance_2d::{point_to_arc_dist, point_to_segment_dist};

use super::slice::PlineSlice;

/// Filters slices, keeping only those whose midpoints are at least `|distance| - eps`
/// from the original polyline.
///
/// This removes self-intersection loops that are "too close" to the original,
/// which are artifacts of the offset rather than valid geometry.
#[must_use]
pub fn apply<'a>(
    slices: &'a [PlineSlice],
    original: &Pline,
    distance: f64,
) -> Vec<&'a PlineSlice> {
    let abs_d = distance.abs();
    let threshold = abs_d * 0.5; // Accept slices at ≥ 50% of offset distance.

    slices
        .iter()
        .filter(|s| {
            if s.vertices.len() < 2 {
                return false;
            }
            // Check the midpoint of the slice.
            let mid_idx = s.vertices.len() / 2;
            let mid = &s.vertices[mid_idx];
            let dist = min_dist_to_pline(mid.x, mid.y, original);
            dist >= threshold
        })
        .collect()
}

/// Computes the minimum distance from a point to a polyline.
///
/// Handles both line segments (bulge=0) and arc segments (bulge≠0).
fn min_dist_to_pline(px: f64, py: f64, pline: &Pline) -> f64 {
    let n = pline.vertices.len();
    let seg_count = pline.segment_count();
    let mut min_d = f64::MAX;

    for i in 0..seg_count {
        let v0 = &pline.vertices[i];
        let v1 = &pline.vertices[(i + 1) % n];

        let d = if v0.bulge.abs() < 1e-12 {
            point_to_segment_dist(px, py, v0.x, v0.y, v1.x, v1.y)
        } else {
            let (cx, cy, r, sa, sw) = arc_from_bulge(v0.x, v0.y, v1.x, v1.y, v0.bulge);
            point_to_arc_dist(px, py, cx, cy, r, sa, sw)
        };

        if d < min_d {
            min_d = d;
        }
    }

    min_d
}
