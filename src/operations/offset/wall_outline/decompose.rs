use crate::geometry::pline::Pline;
use crate::math::TOLERANCE;

/// A unique line segment in the centerline network.
#[derive(Debug, Clone)]
pub struct UniqueSegment {
    pub start: (f64, f64),
    pub end: (f64, f64),
}

/// Collects raw segments from a single Pline.
fn collect_raw_segments(pline: &Pline) -> Vec<((f64, f64), (f64, f64))> {
    let verts = &pline.vertices;
    if verts.len() < 2 {
        return Vec::new();
    }

    let seg_count = pline.segment_count();
    let mut raw_segments = Vec::new();
    for i in 0..seg_count {
        let a = (verts[i].x, verts[i].y);
        let next_i = if pline.closed { (i + 1) % verts.len() } else { i + 1 };
        let b = (verts[next_i].x, verts[next_i].y);
        let dx = b.0 - a.0;
        let dy = b.1 - a.1;
        if dx * dx + dy * dy < TOLERANCE * TOLERANCE {
            continue;
        }
        raw_segments.push((a, b));
    }
    raw_segments
}

/// Groups raw segments by supporting line and merges overlapping extents.
fn merge_to_unique_segments(raw_segments: &[((f64, f64), (f64, f64))]) -> Vec<UniqueSegment> {
    if raw_segments.is_empty() {
        return Vec::new();
    }

    let mut groups: Vec<SupportingLine> = Vec::new();

    for &(a, b) in raw_segments {
        let key = supporting_line_key(a, b);
        let (t_a, t_b) = project_extent(&key, a, b);
        let t_min = t_a.min(t_b);
        let t_max = t_a.max(t_b);

        let mut merged = false;
        for g in &mut groups {
            if same_supporting_line(&g.key, &key) {
                g.intervals.push((t_min, t_max));
                merged = true;
                break;
            }
        }
        if !merged {
            groups.push(SupportingLine {
                key,
                intervals: vec![(t_min, t_max)],
            });
        }
    }

    let mut result = Vec::new();
    for g in &mut groups {
        let merged = merge_intervals(&mut g.intervals);
        for &(t_min, t_max) in &merged {
            let start = unproject(&g.key, t_min);
            let end = unproject(&g.key, t_max);
            result.push(UniqueSegment { start, end });
        }
    }

    result
}

/// Decomposes multiple Plines into unique non-overlapping line segments.
///
/// Walks each Pline, groups collinear segments by their supporting line,
/// and merges overlapping extents into a minimal set of unique segments.
pub fn decompose(plines: &[&Pline]) -> Vec<UniqueSegment> {
    let mut raw = Vec::new();
    for pline in plines {
        raw.extend(collect_raw_segments(pline));
    }
    merge_to_unique_segments(&raw)
}

/// Supporting line representation: a point on the line + normalized direction.
///
/// The direction is canonicalized so that `dx > 0`, or `dx == 0 && dy > 0`.
#[derive(Debug, Clone)]
struct LineKey {
    /// A reference point on the line.
    origin: (f64, f64),
    /// Normalized, canonicalized direction.
    dir: (f64, f64),
}

struct SupportingLine {
    key: LineKey,
    intervals: Vec<(f64, f64)>,
}

/// Computes a canonical supporting line key for a segment.
fn supporting_line_key(a: (f64, f64), b: (f64, f64)) -> LineKey {
    let dx = b.0 - a.0;
    let dy = b.1 - a.1;
    let len = (dx * dx + dy * dy).sqrt();
    let (mut nx, mut ny) = (dx / len, dy / len);

    // Canonicalize direction: prefer positive x, or positive y if x is ~0.
    if nx < -TOLERANCE || (nx.abs() < TOLERANCE && ny < 0.0) {
        nx = -nx;
        ny = -ny;
    }

    // Project origin to the foot of the perpendicular from (0,0).
    // foot = a - (a · dir) * dir
    let dot = a.0 * nx + a.1 * ny;
    let origin = (a.0 - dot * nx, a.1 - dot * ny);

    LineKey {
        origin,
        dir: (nx, ny),
    }
}

/// Projects a point onto the supporting line, returning the parameter t.
fn project_param(key: &LineKey, p: (f64, f64)) -> f64 {
    (p.0 - key.origin.0) * key.dir.0 + (p.1 - key.origin.1) * key.dir.1
}

/// Projects both endpoints of a segment onto the supporting line.
fn project_extent(key: &LineKey, a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (project_param(key, a), project_param(key, b))
}

/// Recovers 2D coordinates from a parameter on the supporting line.
fn unproject(key: &LineKey, t: f64) -> (f64, f64) {
    (key.origin.0 + t * key.dir.0, key.origin.1 + t * key.dir.1)
}

/// Checks if two line keys represent the same supporting line.
fn same_supporting_line(a: &LineKey, b: &LineKey) -> bool {
    // Directions must be parallel (same canonical direction).
    let cross = a.dir.0 * b.dir.1 - a.dir.1 * b.dir.0;
    if cross.abs() > TOLERANCE * 100.0 {
        return false;
    }
    let dot = a.dir.0 * b.dir.0 + a.dir.1 * b.dir.1;
    if dot < 1.0 - TOLERANCE * 100.0 {
        return false;
    }

    // Origins must be at the same perpendicular distance from the origin.
    let d = (a.origin.0 - b.origin.0).powi(2) + (a.origin.1 - b.origin.1).powi(2);
    d < TOLERANCE * 100.0
}

/// Merges overlapping intervals and returns the union.
fn merge_intervals(intervals: &mut [(f64, f64)]) -> Vec<(f64, f64)> {
    if intervals.is_empty() {
        return Vec::new();
    }
    intervals.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut merged = vec![intervals[0]];
    for &(lo, hi) in &intervals[1..] {
        let last = merged.last_mut().unwrap_or_else(|| unreachable!());
        if lo <= last.1 + TOLERANCE {
            last.1 = last.1.max(hi);
        } else {
            merged.push((lo, hi));
        }
    }
    merged
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::pline::PlineVertex;

    /// 井-shaped centerline: 4 crossing segments.
    fn double_cross_pline() -> Pline {
        Pline {
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
        }
    }

    #[test]
    fn double_cross_decompose_4_segments() {
        let pline = double_cross_pline();
        let result = decompose(&[&pline]);
        assert_eq!(result.len(), 4, "expected 4 unique segments, got {}", result.len());
    }

    #[test]
    fn single_line_decompose() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(5.0, 0.0),
            ],
            closed: false,
        };
        let result = decompose(&[&pline]);
        assert_eq!(result.len(), 1);
        assert!((result[0].start.0).abs() < 1e-6);
        assert!((result[0].end.0 - 5.0).abs() < 1e-6);
    }

    #[test]
    fn l_shape_decompose() {
        // Two perpendicular segments meeting at a point.
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(5.0, 0.0),
                PlineVertex::line(5.0, 5.0),
            ],
            closed: false,
        };
        let result = decompose(&[&pline]);
        assert_eq!(result.len(), 2, "expected 2 unique segments, got {}", result.len());
    }

    #[test]
    fn closed_square_decompose() {
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(10.0, 0.0),
                PlineVertex::line(10.0, 10.0),
                PlineVertex::line(0.0, 10.0),
            ],
            closed: true,
        };
        let result = decompose(&[&pline]);
        assert_eq!(result.len(), 4, "expected 4 unique segments, got {}", result.len());
    }

    #[test]
    fn collinear_overlap_merges() {
        // Segment goes back along itself: should merge into one.
        let pline = Pline {
            vertices: vec![
                PlineVertex::line(0.0, 0.0),
                PlineVertex::line(10.0, 0.0),
                PlineVertex::line(5.0, 0.0),
            ],
            closed: false,
        };
        let result = decompose(&[&pline]);
        assert_eq!(result.len(), 1, "overlapping collinear should merge");
        // Should cover (0,0) to (10,0).
        let s = &result[0];
        let min_x = s.start.0.min(s.end.0);
        let max_x = s.start.0.max(s.end.0);
        assert!((min_x).abs() < 1e-6);
        assert!((max_x - 10.0).abs() < 1e-6);
    }
}
