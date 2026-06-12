//! Bounding-box subdivision seeding for the intersection solvers.
//!
//! Each candidate box carries a parameter range plus a tight AABB obtained by
//! *direct sampling* of the curve/surface over that range. Sampling (rather than
//! the control-hull AABB exposed by `NurbsCurve::bounding_box` /
//! `NurbsSurface::bounding_box`) is essential: rational arcs have loose control
//! hulls that fail to shrink at axis-aligned tangents, which explodes the
//! subdivision frontier. A small `pad` at the overlap test plus a coarse leaf
//! size cover the (non-conservative) sampling gap; downstream Newton refinement
//! converges from the resulting coarse seeds. Recursive subdivision prunes
//! parameter space to small overlapping box pairs that seed the solvers.

use std::collections::VecDeque;

use crate::error::Result;
use crate::geometry::nurbs::{NurbsCurve, NurbsCurve3D, NurbsSurface};

/// A parameter interval on one curve, paired with its tight sampled 3D AABB.
#[derive(Debug, Clone, Copy)]
pub(super) struct CurveBox<const D: usize> {
    /// Parameter range `[t0, t1]`.
    pub t0: f64,
    pub t1: f64,
    /// Sampled AABB min/max as flat `D`-arrays.
    pub min: [f64; D],
    pub max: [f64; D],
}

/// A parameter rectangle on one surface, paired with its tight sampled 3D AABB.
#[derive(Debug, Clone, Copy)]
pub(super) struct SurfaceBox {
    pub u0: f64,
    pub u1: f64,
    pub v0: f64,
    pub v1: f64,
    pub min: [f64; 3],
    pub max: [f64; 3],
}

/// Tests whether two `D`-dimensional AABBs overlap, with a tolerance pad so
/// tangential contacts are not missed by floating-point gaps.
fn aabb_overlap<const D: usize>(
    a_min: &[f64; D],
    a_max: &[f64; D],
    b_min: &[f64; D],
    b_max: &[f64; D],
    pad: f64,
) -> bool {
    for k in 0..D {
        if a_max[k] + pad < b_min[k] || b_max[k] + pad < a_min[k] {
            return false;
        }
    }
    true
}

/// Diagonal extent of a `D`-dimensional AABB (used as the size threshold).
fn box_extent<const D: usize>(min: &[f64; D], max: &[f64; D]) -> f64 {
    let mut sum = 0.0;
    for k in 0..D {
        let d = max[k] - min[k];
        sum += d * d;
    }
    sum.sqrt()
}

/// Builds a `CurveBox` for the whole curve (tight sampled AABB).
fn curve_root<const D: usize>(c: &NurbsCurve<D>) -> Result<CurveBox<D>> {
    let (t0, t1) = c.parameter_domain();
    sampled_curve_box(c, t0, t1)
}

/// Hard cap on the number of candidate pairs any seeder may emit. Far above
/// any transversal-intersection count; a runaway (near-coincident input) hits
/// it and terminates instead of hanging.
const MAX_SEED_PAIRS: usize = 4096;

/// Hard cap on the number of subdivision nodes processed. Bounds the breadth of
/// the BFS frontier so axis-aligned tangents (where AABB pruning is weak) cannot
/// explode the work. When hit, remaining frontier boxes are flushed as coarse
/// seeds — Newton refinement downstream still converges from them.
const MAX_SEED_ITERS: usize = 200_000;

/// Recursively subdivides two curves (BFS), collecting overlapping
/// parameter-box pairs whose sampled AABBs overlap and whose extents fall below
/// `leaf_extent` (or once `depth` / the iteration guard is exhausted). `pad`
/// accounts for tangential contact and the sampling gap.
pub(super) fn seed_curve_curve<const D: usize>(
    a: &NurbsCurve<D>,
    b: &NurbsCurve<D>,
    leaf_extent: f64,
    pad: f64,
    max_depth: usize,
) -> Result<Vec<(CurveBox<D>, CurveBox<D>)>> {
    let mut out = Vec::new();
    let mut stack: VecDeque<_> = VecDeque::new();
    stack.push_back((curve_root(a)?, curve_root(b)?, 0usize));
    let mut iters = 0usize;
    while let Some((ba, bb, depth)) = stack.pop_front() {
        iters += 1;
        if out.len() >= MAX_SEED_PAIRS || iters > MAX_SEED_ITERS {
            out.push((ba, bb));
            out.extend(stack.into_iter().map(|(x, y, _)| (x, y)));
            break;
        }
        if !aabb_overlap(&ba.min, &ba.max, &bb.min, &bb.max, pad) {
            continue;
        }
        let small_a = box_extent(&ba.min, &ba.max) < leaf_extent || !can_split_curve(&ba);
        let small_b = box_extent(&bb.min, &bb.max) < leaf_extent || !can_split_curve(&bb);
        if (small_a && small_b) || depth >= max_depth {
            out.push((ba, bb));
            continue;
        }
        // Split BOTH splittable boxes and recurse on overlapping child pairs.
        // Splitting both (rather than one) keeps intermediate nodes from
        // duplicating an unsplit large box across many small children, which is
        // what explodes the frontier at axis-aligned tangents.
        let children_a = if small_a {
            vec![ba]
        } else {
            split_curve_box(a, &ba)?
        };
        let children_b = if small_b {
            vec![bb]
        } else {
            split_curve_box(b, &bb)?
        };
        for &ca in &children_a {
            for &cb in &children_b {
                if aabb_overlap(&ca.min, &ca.max, &cb.min, &cb.max, pad) {
                    stack.push_back((ca, cb, depth + 1));
                }
            }
        }
    }
    out.truncate(MAX_SEED_PAIRS);
    Ok(out)
}

/// Number of samples used to bound a curve sub-interval. A tight sampled AABB
/// (vs. the loose rational control-hull) is what keeps the subdivision frontier
/// from exploding at axis-aligned tangents.
const CURVE_SAMPLES: usize = 8;

/// Whether a curve box's parameter span is wide enough to split.
fn can_split_curve<const D: usize>(b: &CurveBox<D>) -> bool {
    b.t1 - b.t0 > 1e-7
}

/// Bisects the box's parameter range and re-bounds each half by direct sampling
/// of the original curve (no geometric split — cheaper and tighter).
fn split_curve_box<const D: usize>(c: &NurbsCurve<D>, b: &CurveBox<D>) -> Result<Vec<CurveBox<D>>> {
    let mid = 0.5 * (b.t0 + b.t1);
    Ok(vec![
        sampled_curve_box(c, b.t0, mid)?,
        sampled_curve_box(c, mid, b.t1)?,
    ])
}

/// Builds a tight AABB for `c` over `[t0, t1]` by sampling `CURVE_SAMPLES`
/// interior+endpoint points. The convex-hull property guarantees the true
/// sub-curve hull is contained in the control hull, but sampling gives a far
/// tighter (though not strictly conservative) box; a `pad` at the overlap test
/// covers the sampling gap.
// Sample-index to f64 conversions are exact for the small sample counts.
#[allow(clippy::cast_precision_loss)]
fn sampled_curve_box<const D: usize>(c: &NurbsCurve<D>, t0: f64, t1: f64) -> Result<CurveBox<D>> {
    let mut min = [f64::INFINITY; D];
    let mut max = [f64::NEG_INFINITY; D];
    for i in 0..=CURVE_SAMPLES {
        let t = t0 + (t1 - t0) * (i as f64) / (CURVE_SAMPLES as f64);
        let p = c.point_at(t)?;
        for k in 0..D {
            min[k] = min[k].min(p.coords[k]);
            max[k] = max[k].max(p.coords[k]);
        }
    }
    Ok(CurveBox { t0, t1, min, max })
}

/// Samples per parametric direction when bounding a surface sub-rectangle.
const SURFACE_SAMPLES: usize = 6;

/// Builds a `SurfaceBox` for the whole surface (tight sampled AABB).
fn surface_root(s: &NurbsSurface) -> Result<SurfaceBox> {
    let ((u0, u1), (v0, v1)) = s.parameter_domain();
    sampled_surface_box(s, u0, u1, v0, v1)
}

/// Builds a tight AABB for `s` over the parameter rectangle `[u0,u1]x[v0,v1]`
/// by direct sampling (no geometric split — cheaper and far tighter than the
/// rational control hull).
// Sample-index to f64 conversions are exact for the small sample counts.
#[allow(clippy::cast_precision_loss)]
fn sampled_surface_box(s: &NurbsSurface, u0: f64, u1: f64, v0: f64, v1: f64) -> Result<SurfaceBox> {
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for iu in 0..=SURFACE_SAMPLES {
        let u = u0 + (u1 - u0) * (iu as f64) / (SURFACE_SAMPLES as f64);
        for iv in 0..=SURFACE_SAMPLES {
            let v = v0 + (v1 - v0) * (iv as f64) / (SURFACE_SAMPLES as f64);
            let p = s.point_at(u, v)?;
            min[0] = min[0].min(p.x);
            min[1] = min[1].min(p.y);
            min[2] = min[2].min(p.z);
            max[0] = max[0].max(p.x);
            max[1] = max[1].max(p.y);
            max[2] = max[2].max(p.z);
        }
    }
    Ok(SurfaceBox {
        u0,
        u1,
        v0,
        v1,
        min,
        max,
    })
}

/// Whether a surface box still has a splittable parameter span on either axis.
fn can_split_surface(b: &SurfaceBox) -> bool {
    (b.u1 - b.u0) > 1e-7 || (b.v1 - b.v0) > 1e-7
}

/// Bisects the box along its longer parameter side and re-bounds each half by
/// direct sampling of the original surface.
fn split_surface_box(s: &NurbsSurface, b: &SurfaceBox) -> Result<Vec<SurfaceBox>> {
    let du = b.u1 - b.u0;
    let dv = b.v1 - b.v0;
    if du >= dv && du > 1e-7 {
        let mid = 0.5 * (b.u0 + b.u1);
        Ok(vec![
            sampled_surface_box(s, b.u0, mid, b.v0, b.v1)?,
            sampled_surface_box(s, mid, b.u1, b.v0, b.v1)?,
        ])
    } else if dv > 1e-7 {
        let mid = 0.5 * (b.v0 + b.v1);
        Ok(vec![
            sampled_surface_box(s, b.u0, b.u1, b.v0, mid)?,
            sampled_surface_box(s, b.u0, b.u1, mid, b.v1)?,
        ])
    } else {
        Ok(vec![*b])
    }
}

/// Recursively subdivides a curve against a surface, collecting overlapping
/// `(CurveBox, SurfaceBox)` pairs below `leaf_extent`.
pub(super) fn seed_curve_surface(
    c: &NurbsCurve3D,
    s: &NurbsSurface,
    leaf_extent: f64,
    pad: f64,
    max_depth: usize,
) -> Result<Vec<(CurveBox<3>, SurfaceBox)>> {
    let mut out = Vec::new();
    let mut stack: VecDeque<_> = VecDeque::new();
    stack.push_back((curve_root(c)?, surface_root(s)?, 0usize));
    let mut iters = 0usize;
    while let Some((cb, sb, depth)) = stack.pop_front() {
        iters += 1;
        if out.len() >= MAX_SEED_PAIRS || iters > MAX_SEED_ITERS {
            out.push((cb, sb));
            out.extend(stack.into_iter().map(|(x, y, _)| (x, y)));
            break;
        }
        if !aabb_overlap(&cb.min, &cb.max, &sb.min, &sb.max, pad) {
            continue;
        }
        let small_c = box_extent(&cb.min, &cb.max) < leaf_extent || !can_split_curve(&cb);
        let small_s = box_extent(&sb.min, &sb.max) < leaf_extent || !can_split_surface(&sb);
        if (small_c && small_s) || depth >= max_depth {
            out.push((cb, sb));
            continue;
        }
        let children_c = if small_c {
            vec![cb]
        } else {
            split_curve_box(c, &cb)?
        };
        let children_s = if small_s {
            vec![sb]
        } else {
            split_surface_box(s, &sb)?
        };
        for &cc in &children_c {
            for &cs in &children_s {
                if aabb_overlap(&cc.min, &cc.max, &cs.min, &cs.max, pad) {
                    stack.push_back((cc, cs, depth + 1));
                }
            }
        }
    }
    out.truncate(MAX_SEED_PAIRS);
    Ok(out)
}

/// Recursively subdivides two surfaces, collecting overlapping
/// `(SurfaceBox, SurfaceBox)` pairs below `leaf_extent`.
pub(super) fn seed_surface_surface(
    a: &NurbsSurface,
    b: &NurbsSurface,
    leaf_extent: f64,
    pad: f64,
    max_depth: usize,
) -> Result<Vec<(SurfaceBox, SurfaceBox)>> {
    let mut out = Vec::new();
    let mut stack: VecDeque<_> = VecDeque::new();
    stack.push_back((surface_root(a)?, surface_root(b)?, 0usize));
    let mut iters = 0usize;
    while let Some((ba, bb, depth)) = stack.pop_front() {
        iters += 1;
        if out.len() >= MAX_SEED_PAIRS || iters > MAX_SEED_ITERS {
            out.push((ba, bb));
            out.extend(stack.into_iter().map(|(x, y, _)| (x, y)));
            break;
        }
        if !aabb_overlap(&ba.min, &ba.max, &bb.min, &bb.max, pad) {
            continue;
        }
        let small_a = box_extent(&ba.min, &ba.max) < leaf_extent || !can_split_surface(&ba);
        let small_b = box_extent(&bb.min, &bb.max) < leaf_extent || !can_split_surface(&bb);
        if (small_a && small_b) || depth >= max_depth {
            out.push((ba, bb));
            continue;
        }
        let children_a = if small_a {
            vec![ba]
        } else {
            split_surface_box(a, &ba)?
        };
        let children_b = if small_b {
            vec![bb]
        } else {
            split_surface_box(b, &bb)?
        };
        for &ca in &children_a {
            for &cb in &children_b {
                if aabb_overlap(&ca.min, &ca.max, &cb.min, &cb.max, pad) {
                    stack.push_back((ca, cb, depth + 1));
                }
            }
        }
    }
    out.truncate(MAX_SEED_PAIRS);
    Ok(out)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::nurbs::KnotVector;
    use crate::math::{Point2, Point3, Vector3};
    use std::f64::consts::FRAC_PI_2;

    fn quarter_circle_2d(center: Point2, start: f64) -> NurbsCurve<2> {
        // Build a 2D quarter circle by projecting the exact 3D arc's XY data.
        let arc = NurbsCurve3D::arc(
            Point3::new(center.x, center.y, 0.0),
            1.0,
            Vector3::z(),
            Vector3::x(),
            start,
            start + FRAC_PI_2,
        )
        .unwrap();
        let pts: Vec<Point2> = arc
            .control_points()
            .iter()
            .map(|p| Point2::new(p.x, p.y))
            .collect();
        NurbsCurve::<2>::new(
            pts,
            arc.weights().to_vec(),
            arc.knots().clone(),
            arc.degree(),
        )
        .unwrap()
    }

    #[test]
    fn crossing_quarter_circles_seed_near_true_intersection() {
        // Circle A centered at origin, quarter in first quadrant.
        let a = quarter_circle_2d(Point2::new(0.0, 0.0), 0.0);
        // Circle B centered at (1,0): its left quarter sweeps through (0,0)..(1,1).
        // They cross near (0.5, sqrt(3)/2) ~ (0.5, 0.866).
        let b = quarter_circle_2d(Point2::new(1.0, 0.0), FRAC_PI_2);
        let seeds = seed_curve_curve(&a, &b, 0.05, 1e-9, 24).unwrap();
        assert!(!seeds.is_empty(), "expected candidate boxes");
        // The true crossing is at x^2+y^2=1 and (x-1)^2+y^2=1 -> x=0.5, y=sqrt(3)/2.
        let tx = 0.5;
        let ty = (3.0_f64).sqrt() / 2.0;
        let hit = seeds.iter().any(|(ba, bb)| {
            tx >= ba.min[0] - 0.06
                && tx <= ba.max[0] + 0.06
                && ty >= ba.min[1] - 0.06
                && ty <= ba.max[1] + 0.06
                && tx >= bb.min[0] - 0.06
                && tx <= bb.max[0] + 0.06
        });
        assert!(hit, "no candidate box brackets the true intersection");
    }

    #[test]
    fn disjoint_curves_seed_empty() {
        let a = quarter_circle_2d(Point2::new(0.0, 0.0), 0.0);
        let b = quarter_circle_2d(Point2::new(10.0, 10.0), 0.0);
        let seeds = seed_curve_curve(&a, &b, 0.05, 1e-9, 24).unwrap();
        assert!(seeds.is_empty(), "disjoint curves must not seed");
    }

    #[test]
    fn line_through_patch_seeds_surface_box() {
        let s = NurbsSurface::from_unweighted(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 2.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
                Point3::new(2.0, 2.0, 0.0),
            ],
            2,
            2,
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
            1,
        )
        .unwrap();
        let c = NurbsCurve3D::from_unweighted(
            vec![Point3::new(1.0, 1.0, -1.0), Point3::new(1.0, 1.0, 1.0)],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
        )
        .unwrap();
        let seeds = seed_curve_surface(&c, &s, 0.2, 1e-9, 20).unwrap();
        assert!(!seeds.is_empty());
    }

    fn planar_patch(z: f64, x_lo: f64, x_hi: f64) -> NurbsSurface {
        NurbsSurface::from_unweighted(
            vec![
                Point3::new(x_lo, 0.0, z),
                Point3::new(x_lo, 2.0, z),
                Point3::new(x_hi, 0.0, z),
                Point3::new(x_hi, 2.0, z),
            ],
            2,
            2,
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
            1,
        )
        .unwrap()
    }

    #[test]
    fn overlapping_patches_seed_pairs() {
        let a = planar_patch(0.0, 0.0, 2.0);
        // A vertical patch crossing through a -> their hulls overlap.
        let b = NurbsSurface::from_unweighted(
            vec![
                Point3::new(1.0, 0.0, -1.0),
                Point3::new(1.0, 2.0, -1.0),
                Point3::new(1.0, 0.0, 1.0),
                Point3::new(1.0, 2.0, 1.0),
            ],
            2,
            2,
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
            1,
        )
        .unwrap();
        let seeds = seed_surface_surface(&a, &b, 0.5, 1e-9, 16).unwrap();
        assert!(!seeds.is_empty());
    }

    #[test]
    fn disjoint_patches_seed_empty() {
        let a = planar_patch(0.0, 0.0, 2.0);
        let b = planar_patch(5.0, 0.0, 2.0);
        let seeds = seed_surface_surface(&a, &b, 0.5, 1e-9, 16).unwrap();
        assert!(seeds.is_empty());
    }
}
