use crate::error::{OperationError, Result};
use crate::geometry::pline::{Pline, PlineVertex};
use crate::math::arc_2d::{arc_from_bulge, arc_tangent_at, bulge_from_arc, offset_arc_segment};
use crate::math::intersect_2d::line_line_intersect_2d;
use crate::math::polygon_2d::{left_normal, segment_direction};
use crate::math::{Point3, TOLERANCE};

/// Maximum miter distance as a multiple of `|distance|`.
const MITER_LIMIT: f64 = 4.0;

/// Threshold for flat cap: `cos(angle) < this` → near-180° reversal.
const FLAT_CAP_COS: f64 = -0.98;

/// The carrier curve an offset segment lies on — the infinite line /
/// full circle used for EXACT corner joins. Intersecting carriers (not
/// tangent-line approximations) keeps every joined arc endpoint ON its
/// offset circle, so the re-derived bulge encodes the exact concentric
/// offset arc.
#[derive(Clone, Copy)]
enum Carrier {
    Line,
    Circle { cx: f64, cy: f64, r: f64, ccw: bool },
}

/// An offset segment with endpoints, carrier, and tangent directions.
struct OffsetSeg {
    start: (f64, f64),
    end: (f64, f64),
    carrier: Carrier,
    /// Unit tangent direction at the start of the segment.
    start_dir: (f64, f64),
    /// Unit tangent direction at the end of the segment.
    end_dir: (f64, f64),
}

/// A resolved corner between two consecutive offset segments.
enum Join {
    /// Single exact corner point shared by both segments.
    Miter((f64, f64)),
    /// Two points — the previous segment's own end, then the next
    /// segment's own start — connected by a straight bevel span
    /// (flat cap, miter-limit bevel, or disjoint-carrier fallback).
    Bevel((f64, f64), (f64, f64)),
}

impl Join {
    /// The point the PREVIOUS segment ends at.
    fn prev_end(&self) -> (f64, f64) {
        match self {
            Self::Miter(p) | Self::Bevel(p, _) => *p,
        }
    }

    /// The point the NEXT segment starts at.
    fn next_start(&self) -> (f64, f64) {
        match self {
            Self::Miter(p) | Self::Bevel(_, p) => *p,
        }
    }
}

/// Builds the raw (untrimmed) offset polyline by offsetting each segment
/// and connecting them at corners.
///
/// Supports both line segments (bulge=0) and arc segments (bulge≠0).
/// Corners involving an arc are joined at the exact intersection of the
/// segments' carrier curves (line × circle / circle × circle), and every
/// arc's bulge is re-derived from its final endpoints about its exact
/// offset circle — the joined arc stays concentric with its source at
/// `r ± distance`.
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
    let offset_segs = offset_segments(pline, distance)?;

    // Phase B: resolve every corner, then assemble vertices with each
    // segment's bulge re-derived from its FINAL endpoints (a joined arc
    // endpoint moves along the offset circle, so the bulge must encode
    // the trimmed / extended sweep, not the source sweep).
    let mut verts = Vec::with_capacity(n * 2);

    if pline.closed {
        let joins: Vec<Join> = (0..seg_count)
            .map(|i| {
                let prev = if i == 0 { seg_count - 1 } else { i - 1 };
                corner_join(
                    &offset_segs[prev],
                    &offset_segs[i],
                    pline.vertices[i].x,
                    pline.vertices[i].y,
                    distance,
                )
            })
            .collect();
        for i in 0..seg_count {
            let start = joins[i].next_start();
            let end = joins[(i + 1) % seg_count].prev_end();
            if let Join::Bevel(a, _) = &joins[i] {
                verts.push(PlineVertex::line(a.0, a.1));
            }
            verts.push(PlineVertex::new(
                start.0,
                start.1,
                seg_bulge(&offset_segs[i], start, end),
            ));
        }
    } else {
        // Joins at interior vertices 1..seg_count; the open ends keep
        // their segments' own endpoints.
        let joins: Vec<Join> = (1..seg_count)
            .map(|i| {
                corner_join(
                    &offset_segs[i - 1],
                    &offset_segs[i],
                    pline.vertices[i].x,
                    pline.vertices[i].y,
                    distance,
                )
            })
            .collect();
        for i in 0..seg_count {
            let start = if i == 0 {
                offset_segs[0].start
            } else {
                joins[i - 1].next_start()
            };
            let end = if i + 1 < seg_count {
                joins[i].prev_end()
            } else {
                offset_segs[seg_count - 1].end
            };
            if i > 0 {
                if let Join::Bevel(a, _) = &joins[i - 1] {
                    verts.push(PlineVertex::line(a.0, a.1));
                }
            }
            verts.push(PlineVertex::new(
                start.0,
                start.1,
                seg_bulge(&offset_segs[i], start, end),
            ));
        }
        let last = &offset_segs[seg_count - 1];
        verts.push(PlineVertex::line(last.end.0, last.end.1));
    }

    Ok(Pline {
        vertices: verts,
        closed: pline.closed,
    })
}

/// Offsets every segment of `pline` individually: lines shift along
/// their left normal, arcs stay concentric with the source arc at the
/// left-offset radius (`offset_arc_segment`).
///
/// # Errors
///
/// Returns `OperationError::InvalidInput` for zero-length line segments
/// or `OperationError::Failed` when an arc collapses under the offset.
fn offset_segments(pline: &Pline, distance: f64) -> Result<Vec<OffsetSeg>> {
    let n = pline.vertices.len();
    let seg_count = pline.segment_count();
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
                carrier: Carrier::Line,
                start_dir: d,
                end_dir: d,
            });
        } else {
            // Arc segment: change radius, preserve sweep — concentric
            // with the source arc.
            let seg = offset_arc_segment(v0.x, v0.y, v1.x, v1.y, v0.bulge, distance).ok_or_else(
                || OperationError::Failed("arc segment collapsed during offset".to_owned()),
            )?;

            let (ox0, oy0, ox1, oy1, ob) = seg;
            let (cx, cy, _, sa, sw) = arc_from_bulge(ox0, oy0, ox1, oy1, ob);
            let r = ((ox0 - cx).powi(2) + (oy0 - cy).powi(2)).sqrt();
            let sd = arc_tangent_at(sa, sw, 0.0);
            let ed = arc_tangent_at(sa, sw, 1.0);

            offset_segs.push(OffsetSeg {
                start: (ox0, oy0),
                end: (ox1, oy1),
                carrier: Carrier::Circle {
                    cx,
                    cy,
                    r,
                    ccw: ob > 0.0,
                },
                start_dir: sd,
                end_dir: ed,
            });
        }
    }
    Ok(offset_segs)
}

/// Bulge of one offset segment between its FINAL endpoints: `0` on a
/// line carrier; on a circle carrier, re-derived about the exact offset
/// circle in the source arc's winding.
fn seg_bulge(seg: &OffsetSeg, start: (f64, f64), end: (f64, f64)) -> f64 {
    match seg.carrier {
        Carrier::Line => 0.0,
        Carrier::Circle { cx, cy, ccw, .. } => {
            bulge_from_arc(start.0, start.1, end.0, end.1, cx, cy, ccw)
        }
    }
}

/// Resolves the corner between two consecutive offset segments at the
/// original vertex `(orig_x, orig_y)`.
///
/// Handles three cases:
/// 1. Near-antiparallel (>~169°): flat cap (bevel).
/// 2. Miter too long: bevel.
/// 3. Normal corner: single exact corner point — line × line miter for
///    two straight segments (unchanged legacy math), carrier
///    intersection (line × circle / circle × circle) when an arc is
///    involved so the joined point lies exactly on the offset circle.
fn corner_join(
    seg_prev: &OffsetSeg,
    seg_next: &OffsetSeg,
    orig_x: f64,
    orig_y: f64,
    distance: f64,
) -> Join {
    let dir_prev = &seg_prev.end_dir;
    let dir_next = &seg_next.start_dir;
    let cos_angle = dir_prev.0 * dir_next.0 + dir_prev.1 * dir_next.1;

    if cos_angle < FLAT_CAP_COS {
        // Near-antiparallel: flat cap.
        return Join::Bevel(seg_prev.end, seg_next.start);
    }

    let corner = match (&seg_prev.carrier, &seg_next.carrier) {
        (Carrier::Line, Carrier::Line) => {
            // Legacy line-line miter via tangent intersection (exact for
            // straight carriers).
            let p_prev = Point3::new(seg_prev.end.0, seg_prev.end.1, 0.0);
            let d_prev = crate::math::Vector3::new(dir_prev.0, dir_prev.1, 0.0);
            let p_next = Point3::new(seg_next.start.0, seg_next.start.1, 0.0);
            let d_next = crate::math::Vector3::new(dir_next.0, dir_next.1, 0.0);
            let Some((t, _)) = line_line_intersect_2d(&p_prev, &d_prev, &p_next, &d_next) else {
                // Parallel: use offset of the original corner point
                // (no miter-limit check — matches the legacy path).
                let fallback_normal = left_normal(
                    crate::math::Vector3::new(d_prev.x, d_prev.y, 0.0)
                        .try_normalize(TOLERANCE)
                        .unwrap_or(crate::math::Vector3::new(1.0, 0.0, 0.0)),
                );
                return Join::Miter((
                    orig_x + fallback_normal.x * distance,
                    orig_y + fallback_normal.y * distance,
                ));
            };
            (p_prev.x + d_prev.x * t, p_prev.y + d_prev.y * t)
        }
        pair => {
            // At least one arc: intersect the exact carriers and keep
            // the root nearest the two segment endpoints. Disjoint
            // carriers (arc curving away) get a straight bevel.
            let candidates = match pair {
                (Carrier::Line, Carrier::Circle { cx, cy, r, .. }) => {
                    line_circle_intersections(seg_prev.start, *dir_prev, (*cx, *cy), *r)
                }
                (Carrier::Circle { cx, cy, r, .. }, Carrier::Line) => {
                    line_circle_intersections(seg_next.start, *dir_next, (*cx, *cy), *r)
                }
                (
                    Carrier::Circle { cx, cy, r, .. },
                    Carrier::Circle {
                        cx: bx,
                        cy: by,
                        r: br,
                        ..
                    },
                ) => {
                    // Two arcs split from the SAME carrier circle (a
                    // tangent-continuous junction): the carriers are
                    // coincident, so intersecting them yields nothing
                    // (concentric ⇒ empty) and would force a zero-length
                    // bevel between the coincident endpoints. The join is
                    // the shared offset endpoint itself — the legacy
                    // single-point join.
                    if (cx - bx).abs() < TOLERANCE
                        && (cy - by).abs() < TOLERANCE
                        && (r - br).abs() < TOLERANCE
                    {
                        vec![(
                            0.5 * (seg_prev.end.0 + seg_next.start.0),
                            0.5 * (seg_prev.end.1 + seg_next.start.1),
                        )]
                    } else {
                        circle_circle_intersections((*cx, *cy), *r, (*bx, *by), *br)
                    }
                }
                (Carrier::Line, Carrier::Line) => unreachable!("handled above"),
            };
            let nearest = candidates.into_iter().min_by(|a, b| {
                let score = |p: &(f64, f64)| {
                    (p.0 - seg_prev.end.0).powi(2)
                        + (p.1 - seg_prev.end.1).powi(2)
                        + (p.0 - seg_next.start.0).powi(2)
                        + (p.1 - seg_next.start.1).powi(2)
                };
                score(a).total_cmp(&score(b))
            });
            match nearest {
                Some(p) => p,
                None => return Join::Bevel(seg_prev.end, seg_next.start),
            }
        }
    };

    let dx = corner.0 - orig_x;
    let dy = corner.1 - orig_y;
    let miter_dist_sq = dx * dx + dy * dy;
    let limit = MITER_LIMIT * distance.abs();
    if miter_dist_sq > limit * limit {
        // Miter too long: bevel.
        Join::Bevel(seg_prev.end, seg_next.start)
    } else {
        Join::Miter(corner)
    }
}

/// Intersections of an INFINITE line (point + unit direction) with a
/// FULL circle. Tangential contact yields the single tangent point.
fn line_circle_intersections(
    p: (f64, f64),
    dir: (f64, f64),
    center: (f64, f64),
    radius: f64,
) -> Vec<(f64, f64)> {
    let fx = p.0 - center.0;
    let fy = p.1 - center.1;
    // dir is unit length: t² + 2bt + c = 0.
    let b = fx * dir.0 + fy * dir.1;
    let c = fx * fx + fy * fy - radius * radius;
    let disc = b * b - c;
    if disc < -TOLERANCE {
        return Vec::new();
    }
    if disc <= TOLERANCE {
        let t = -b;
        return vec![(p.0 + dir.0 * t, p.1 + dir.1 * t)];
    }
    let sqrt_disc = disc.sqrt();
    [-b - sqrt_disc, -b + sqrt_disc]
        .iter()
        .map(|t| (p.0 + dir.0 * t, p.1 + dir.1 * t))
        .collect()
}

/// Intersections of two FULL circles. Concentric or disjoint circles
/// yield none; tangential contact yields the single touch point.
fn circle_circle_intersections(a: (f64, f64), ra: f64, b: (f64, f64), rb: f64) -> Vec<(f64, f64)> {
    let dx = b.0 - a.0;
    let dy = b.1 - a.1;
    let dist = (dx * dx + dy * dy).sqrt();
    if dist < TOLERANCE || dist > ra + rb + TOLERANCE || dist < (ra - rb).abs() - TOLERANCE {
        return Vec::new();
    }
    let along = (ra * ra - rb * rb + dist * dist) / (2.0 * dist);
    let h_sq = ra * ra - along * along;
    let h = h_sq.max(0.0).sqrt();
    let mx = a.0 + along * dx / dist;
    let my = a.1 + along * dy / dist;
    if h < TOLERANCE {
        return vec![(mx, my)];
    }
    let px = -dy / dist;
    let py = dx / dist;
    vec![(mx + h * px, my + h * py), (mx - h * px, my - h * py)]
}
