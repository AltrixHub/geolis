use crate::error::{OperationError, Result};
use crate::geometry::pline::{Pline, PlineVertex};
use crate::math::arc_2d::{arc_from_bulge, arc_tangent_at, offset_arc_segment};
use crate::math::intersect_2d::line_line_intersect_2d;
use crate::math::polygon_2d::{left_normal, segment_direction};
use crate::math::{Point3, TOLERANCE};

/// Maximum miter distance as a multiple of `|distance|`.
const MITER_LIMIT: f64 = 4.0;

/// Threshold for flat cap: `cos(angle) < this` → near-180° reversal.
const FLAT_CAP_COS: f64 = -0.98;

/// An offset segment with endpoints, bulge, and tangent directions.
struct OffsetSeg {
    start: (f64, f64),
    end: (f64, f64),
    bulge: f64,
    /// Unit tangent direction at the start of the segment.
    start_dir: (f64, f64),
    /// Unit tangent direction at the end of the segment.
    end_dir: (f64, f64),
}

/// Builds the raw (untrimmed) offset polyline by offsetting each segment
/// and connecting them at corners.
///
/// Supports both line segments (bulge=0) and arc segments (bulge≠0).
///
/// # Errors
///
/// Returns `OperationError::InvalidInput` for zero-length segments or
/// `OperationError::Failed` if no valid segments exist.
pub fn build(pline: &Pline, distance: f64) -> Result<Pline> {
    let n = pline.vertices.len();
    let seg_count = pline.segment_count();
    if seg_count == 0 {
        return Err(OperationError::Failed("no segments to offset".to_owned()).into());
    }

    // Phase A: Compute offset segments.
    let mut offset_segs: Vec<OffsetSeg> = Vec::with_capacity(seg_count);

    for i in 0..seg_count {
        let v0 = &pline.vertices[i];
        let v1 = &pline.vertices[(i + 1) % n];

        if v0.bulge.abs() < 1e-12 {
            // Line segment: parallel offset.
            let p0 = Point3::new(v0.x, v0.y, 0.0);
            let p1 = Point3::new(v1.x, v1.y, 0.0);
            let dir = segment_direction(&p0, &p1)?;
            let normal = left_normal(dir);

            let start = (v0.x + normal.x * distance, v0.y + normal.y * distance);
            let end = (v1.x + normal.x * distance, v1.y + normal.y * distance);
            let d = (dir.x, dir.y);

            offset_segs.push(OffsetSeg {
                start,
                end,
                bulge: 0.0,
                start_dir: d,
                end_dir: d,
            });
        } else {
            // Arc segment: change radius, preserve sweep.
            let seg = offset_arc_segment(v0.x, v0.y, v1.x, v1.y, v0.bulge, distance)
                .ok_or_else(|| {
                    OperationError::Failed("arc segment collapsed during offset".to_owned())
                })?;

            let (ox0, oy0, ox1, oy1, ob) = seg;
            let (_, _, _, sa, sw) = arc_from_bulge(ox0, oy0, ox1, oy1, ob);
            let sd = arc_tangent_at(sa, sw, 0.0);
            let ed = arc_tangent_at(sa, sw, 1.0);

            offset_segs.push(OffsetSeg {
                start: (ox0, oy0),
                end: (ox1, oy1),
                bulge: ob,
                start_dir: sd,
                end_dir: ed,
            });
        }
    }

    // Phase B: Build raw offset by connecting consecutive offset segments at corners.
    let mut verts = Vec::with_capacity(n * 2);

    if pline.closed {
        for i in 0..seg_count {
            let prev = if i == 0 { seg_count - 1 } else { i - 1 };
            push_corner_and_seg_start(
                &mut verts,
                &offset_segs[prev],
                &offset_segs[i],
                pline.vertices[i].x,
                pline.vertices[i].y,
                distance,
            );
        }
    } else {
        // Open polyline: start with the first offset segment's start point.
        verts.push(PlineVertex::new(
            offset_segs[0].start.0,
            offset_segs[0].start.1,
            offset_segs[0].bulge,
        ));

        for i in 1..seg_count {
            push_corner_and_seg_start(
                &mut verts,
                &offset_segs[i - 1],
                &offset_segs[i],
                pline.vertices[i].x,
                pline.vertices[i].y,
                distance,
            );
        }

        // End with the last offset segment's end point.
        let last = &offset_segs[seg_count - 1];
        verts.push(PlineVertex::line(last.end.0, last.end.1));
    }

    Ok(Pline {
        vertices: verts,
        closed: pline.closed,
    })
}

/// Pushes corner vertex/vertices between two consecutive offset segments,
/// then sets the last pushed vertex's bulge to the next segment's bulge.
///
/// Handles three cases:
/// 1. Near-antiparallel (>~169°): flat cap (two vertices)
/// 2. Miter too long: bevel (two vertices)
/// 3. Normal corner: single miter intersection point
#[allow(clippy::too_many_arguments)]
fn push_corner_and_seg_start(
    verts: &mut Vec<PlineVertex>,
    seg_prev: &OffsetSeg,
    seg_next: &OffsetSeg,
    orig_x: f64,
    orig_y: f64,
    distance: f64,
) {
    let dir_prev = &seg_prev.end_dir;
    let dir_next = &seg_next.start_dir;
    let cos_angle = dir_prev.0 * dir_next.0 + dir_prev.1 * dir_next.1;

    if cos_angle < FLAT_CAP_COS {
        // Near-antiparallel: flat cap.
        verts.push(PlineVertex::line(seg_prev.end.0, seg_prev.end.1));
        verts.push(PlineVertex::new(
            seg_next.start.0,
            seg_next.start.1,
            seg_next.bulge,
        ));
        return;
    }

    // Try miter intersection using tangent directions at the join point.
    let p_prev = Point3::new(seg_prev.end.0, seg_prev.end.1, 0.0);
    let d_prev = crate::math::Vector3::new(dir_prev.0, dir_prev.1, 0.0);
    let p_next = Point3::new(seg_next.start.0, seg_next.start.1, 0.0);
    let d_next = crate::math::Vector3::new(dir_next.0, dir_next.1, 0.0);

    if let Some((t, _)) = line_line_intersect_2d(&p_prev, &d_prev, &p_next, &d_next) {
        let corner_x = p_prev.x + d_prev.x * t;
        let corner_y = p_prev.y + d_prev.y * t;

        let dx = corner_x - orig_x;
        let dy = corner_y - orig_y;
        let miter_dist_sq = dx * dx + dy * dy;
        let limit = MITER_LIMIT * distance.abs();

        if miter_dist_sq > limit * limit {
            // Miter too long: bevel.
            verts.push(PlineVertex::line(seg_prev.end.0, seg_prev.end.1));
            verts.push(PlineVertex::new(
                seg_next.start.0,
                seg_next.start.1,
                seg_next.bulge,
            ));
        } else {
            verts.push(PlineVertex::new(corner_x, corner_y, seg_next.bulge));
        }
    } else {
        // Parallel: use offset of the original corner point.
        let fallback_normal = left_normal(
            crate::math::Vector3::new(d_prev.x, d_prev.y, 0.0)
                .try_normalize(TOLERANCE)
                .unwrap_or(crate::math::Vector3::new(1.0, 0.0, 0.0)),
        );
        verts.push(PlineVertex::new(
            orig_x + fallback_normal.x * distance,
            orig_y + fallback_normal.y * distance,
            seg_next.bulge,
        ));
    }
}
