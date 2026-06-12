//! Bounding-box subdivision seeding for the intersection solvers.
//!
//! The convex-hull property of NURBS gives a cheap conservative AABB from the
//! control net (see `NurbsCurve::bounding_box` / `NurbsSurface::bounding_box`).
//! Recursive subdivision prunes parameter space to small boxes whose control
//! hulls still overlap; the box centres seed the Newton solvers downstream.

use crate::error::Result;
use crate::geometry::nurbs::{NurbsCurve, NurbsCurve3D, NurbsSurface};

/// A parameter interval on one curve, paired with its 3D control-hull AABB.
#[derive(Debug, Clone, Copy)]
pub(super) struct CurveBox<const D: usize> {
    /// Parameter range `[t0, t1]`.
    pub t0: f64,
    pub t1: f64,
    /// Control-hull AABB min/max as flat `D`-arrays.
    pub min: [f64; D],
    pub max: [f64; D],
}

/// A parameter rectangle on one surface, paired with its 3D control-hull AABB.
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

/// Builds a `CurveBox` for the whole curve.
fn curve_root<const D: usize>(c: &NurbsCurve<D>) -> CurveBox<D> {
    let (t0, t1) = c.parameter_domain();
    let (lo, hi) = c.bounding_box();
    let mut min = [0.0; D];
    let mut max = [0.0; D];
    for k in 0..D {
        min[k] = lo.coords[k];
        max[k] = hi.coords[k];
    }
    CurveBox { t0, t1, min, max }
}

/// Hard cap on the number of candidate pairs any seeder may emit. Far above
/// any transversal-intersection count; a runaway (near-coincident input) hits
/// it and terminates instead of hanging.
const MAX_SEED_PAIRS: usize = 4096;

/// Recursively subdivides two curves, collecting overlapping parameter-box
/// pairs whose 3D hulls overlap and whose extents fall below `leaf_extent`
/// (or once `depth` is exhausted). `pad` accounts for tangential contact.
///
/// Only the *larger* of the two boxes is split at each step (interval
/// subdivision), so the work scales with the number of surviving leaf pairs
/// rather than exponentially in depth.
pub(super) fn seed_curve_curve<const D: usize>(
    a: &NurbsCurve<D>,
    b: &NurbsCurve<D>,
    leaf_extent: f64,
    pad: f64,
    max_depth: usize,
) -> Result<Vec<(CurveBox<D>, CurveBox<D>)>> {
    let mut out = Vec::new();
    let mut stack = vec![(curve_root(a), curve_root(b), 0usize)];
    while let Some((ba, bb, depth)) = stack.pop() {
        if out.len() >= MAX_SEED_PAIRS {
            break;
        }
        if !aabb_overlap(&ba.min, &ba.max, &bb.min, &bb.max, pad) {
            continue;
        }
        let ext_a = box_extent(&ba.min, &ba.max);
        let ext_b = box_extent(&bb.min, &bb.max);
        let small_a = ext_a < leaf_extent || !can_split_curve(&ba);
        let small_b = ext_b < leaf_extent || !can_split_curve(&bb);
        if (small_a && small_b) || depth >= max_depth {
            out.push((ba, bb));
            continue;
        }
        // Split whichever box is geometrically larger and still splittable.
        if !small_a && (small_b || ext_a >= ext_b) {
            for ca in split_curve_box(a, &ba)? {
                stack.push((ca, bb, depth + 1));
            }
        } else {
            for cb in split_curve_box(b, &bb)? {
                stack.push((ba, cb, depth + 1));
            }
        }
    }
    Ok(out)
}

/// Whether a curve box's parameter span is wide enough to split.
fn can_split_curve<const D: usize>(b: &CurveBox<D>) -> bool {
    b.t1 - b.t0 > 1e-7
}

/// Splits a curve box at its parameter midpoint, recomputing child hulls.
fn split_curve_box<const D: usize>(c: &NurbsCurve<D>, b: &CurveBox<D>) -> Result<Vec<CurveBox<D>>> {
    let mid = 0.5 * (b.t0 + b.t1);
    let (left, right) = c.split(mid)?;
    Ok(vec![
        sub_curve_box(&left, b.t0, mid),
        sub_curve_box(&right, mid, b.t1),
    ])
}

/// Wraps a sub-curve's hull with the parameter range it covers.
fn sub_curve_box<const D: usize>(sub: &NurbsCurve<D>, t0: f64, t1: f64) -> CurveBox<D> {
    let (lo, hi) = sub.bounding_box();
    let mut min = [0.0; D];
    let mut max = [0.0; D];
    for k in 0..D {
        min[k] = lo.coords[k];
        max[k] = hi.coords[k];
    }
    CurveBox { t0, t1, min, max }
}

/// Builds a `SurfaceBox` for the whole surface.
fn surface_root(s: &NurbsSurface) -> SurfaceBox {
    let ((u0, u1), (v0, v1)) = s.parameter_domain();
    let (lo, hi) = s.bounding_box();
    SurfaceBox {
        u0,
        u1,
        v0,
        v1,
        min: [lo.x, lo.y, lo.z],
        max: [hi.x, hi.y, hi.z],
    }
}

/// Wraps a sub-surface's hull with the parameter rectangle it covers.
fn sub_surface_box(sub: &NurbsSurface, u0: f64, u1: f64, v0: f64, v1: f64) -> SurfaceBox {
    let (lo, hi) = sub.bounding_box();
    SurfaceBox {
        u0,
        u1,
        v0,
        v1,
        min: [lo.x, lo.y, lo.z],
        max: [hi.x, hi.y, hi.z],
    }
}

/// Whether a surface box still has a splittable parameter span on either axis.
fn can_split_surface(b: &SurfaceBox) -> bool {
    (b.u1 - b.u0) > 1e-7 || (b.v1 - b.v0) > 1e-7
}

/// Splits the sub-patch `sub` (whose domain equals box `b`'s rectangle) along
/// its longer parameter side, recomputing child hulls.
fn split_surface_box(sub: &NurbsSurface, b: &SurfaceBox) -> Result<Vec<SurfaceBox>> {
    let du = b.u1 - b.u0;
    let dv = b.v1 - b.v0;
    if du >= dv && du > 1e-7 {
        let mid = 0.5 * (b.u0 + b.u1);
        let (left, right) = sub.split_u(mid)?;
        Ok(vec![
            sub_surface_box(&left, b.u0, mid, b.v0, b.v1),
            sub_surface_box(&right, mid, b.u1, b.v0, b.v1),
        ])
    } else if dv > 1e-7 {
        let mid = 0.5 * (b.v0 + b.v1);
        let (left, right) = sub.split_v(mid)?;
        Ok(vec![
            sub_surface_box(&left, b.u0, b.u1, b.v0, mid),
            sub_surface_box(&right, b.u0, b.u1, mid, b.v1),
        ])
    } else {
        Ok(vec![*b])
    }
}

/// The surface sub-patch covering box `b` (for re-hulling its children).
fn sub_surface(s: &NurbsSurface, b: &SurfaceBox) -> Result<NurbsSurface> {
    let ((u0, u1), (v0, v1)) = s.parameter_domain();
    let mut cur = s.clone();
    let (mut cu0, mut cu1, mut cv0, mut cv1) = (u0, u1, v0, v1);
    if b.u0 > cu0 + 1e-9 {
        let (_, right) = cur.split_u(b.u0)?;
        cur = right;
        cu0 = b.u0;
    }
    if b.u1 < cu1 - 1e-9 {
        // Re-derive the local domain after the first split.
        let ((nu0, nu1), _) = cur.parameter_domain();
        let _ = (nu0, nu1, cu1);
        let (left, _) = cur.split_u(b.u1)?;
        cur = left;
        cu1 = b.u1;
    }
    if b.v0 > cv0 + 1e-9 {
        let (_, right) = cur.split_v(b.v0)?;
        cur = right;
        cv0 = b.v0;
    }
    if b.v1 < cv1 - 1e-9 {
        let (left, _) = cur.split_v(b.v1)?;
        cur = left;
        cv1 = b.v1;
    }
    let _ = (cu0, cu1, cv0, cv1);
    Ok(cur)
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
    let mut stack = vec![(curve_root(c), surface_root(s), 0usize)];
    while let Some((cb, sb, depth)) = stack.pop() {
        if out.len() >= MAX_SEED_PAIRS {
            break;
        }
        if !aabb_overlap(&cb.min, &cb.max, &sb.min, &sb.max, pad) {
            continue;
        }
        let ext_c = box_extent(&cb.min, &cb.max);
        let ext_s = box_extent(&sb.min, &sb.max);
        let small_c = ext_c < leaf_extent || !can_split_curve(&cb);
        let small_s = ext_s < leaf_extent || !can_split_surface(&sb);
        if (small_c && small_s) || depth >= max_depth {
            out.push((cb, sb));
            continue;
        }
        if !small_c && (small_s || ext_c >= ext_s) {
            for cc in split_curve_box(c, &cb)? {
                stack.push((cc, sb, depth + 1));
            }
        } else {
            let sub = sub_surface(s, &sb)?;
            for cs in split_surface_box(&sub, &sb)? {
                stack.push((cb, cs, depth + 1));
            }
        }
    }
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
    let mut stack = vec![(surface_root(a), surface_root(b), 0usize)];
    while let Some((ba, bb, depth)) = stack.pop() {
        if out.len() >= MAX_SEED_PAIRS {
            break;
        }
        if !aabb_overlap(&ba.min, &ba.max, &bb.min, &bb.max, pad) {
            continue;
        }
        let ext_a = box_extent(&ba.min, &ba.max);
        let ext_b = box_extent(&bb.min, &bb.max);
        let small_a = ext_a < leaf_extent || !can_split_surface(&ba);
        let small_b = ext_b < leaf_extent || !can_split_surface(&bb);
        if (small_a && small_b) || depth >= max_depth {
            out.push((ba, bb));
            continue;
        }
        if !small_a && (small_b || ext_a >= ext_b) {
            let sub = sub_surface(a, &ba)?;
            for ca in split_surface_box(&sub, &ba)? {
                stack.push((ca, bb, depth + 1));
            }
        } else {
            let sub = sub_surface(b, &bb)?;
            for cb in split_surface_box(&sub, &bb)? {
                stack.push((ba, cb, depth + 1));
            }
        }
    }
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
