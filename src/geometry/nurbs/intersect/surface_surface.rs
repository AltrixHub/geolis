//! Surface×surface intersection (SSI) via subdivision seeding and marching.
//!
//! ## Pipeline
//!
//! 1. **Seed.** Bounding-box subdivision ([`super::bbox::seed_surface_surface`])
//!    prunes parameter space to candidate box pairs. From each pair's center a
//!    damped Gauss-Newton minimizes `|S_a(u_a,v_a) - S_b(u_b,v_b)|²` over the
//!    four parameters (3 equations, 4 unknowns; the null direction along the
//!    intersection tangent is handled by Levenberg-Marquardt damping). Seeds
//!    whose residual exceeds tolerance are discarded; the rest are deduped.
//! 2. **March.** The intersection direction is `d = n_a × n_b` (the cross product
//!    of the two surface normals). `|d| ≈ 0` ⇒ tangential contact: the branch
//!    stops and records what exists. The predictor steps `step` along `d`; the
//!    corrector is a 4×4 Newton solving `S_a = S_b` plus the marching-plane
//!    constraint `(S_a - anchor)·d = 0` that pins the point's position along the
//!    curve.
//! 3. **Terminate.** Domain boundary of either surface (clamp + emit), closure
//!    (return within tolerance of the start after >3 points), or the
//!    `max_points` runaway guard.
//! 4. **Dedupe branches.** A seed already lying on an existing branch is skipped,
//!    so the two surface families do not yield duplicate traces.
//!
//! ## Corrector formulation
//!
//! The marching corrector uses the **augmented 4×4 system**
//! `[S_a - S_b; (S_a - anchor)·d] = 0` rather than a least-squares relaxation:
//! with the marching-plane row the system is square and well-conditioned for
//! transversal intersections (`|d|` bounded away from zero), giving quadratic
//! Newton convergence. Tangencies (`|d| → 0`) are out of scope beyond clean
//! termination.

use nalgebra::{Matrix4, Vector4};

use crate::error::Result;
use crate::geometry::nurbs::NurbsSurface;
use crate::math::{Point2, Point3, Vector3};

use super::bbox::seed_surface_surface;
use super::types::{IntersectionOptions, SurfaceIntersectionCurve};

/// Computes the intersection branches of two NURBS surfaces.
///
/// Each branch is an ordered 3D polyline with synchronized UV traces on both
/// surfaces (`uv_a`, `uv_b`). Transversal intersections are the quality target;
/// tangential contacts and coincident (overlapping) regions are required only to
/// terminate cleanly and may yield empty, partial, or boundary-only output.
/// Disjoint surfaces yield an empty vector.
///
/// # Errors
///
/// Returns an error only if a sub-box construction during seeding or a surface
/// evaluation fails (e.g. a vanishing rational denominator).
pub fn intersect_surfaces(
    a: &NurbsSurface,
    b: &NurbsSurface,
    options: &IntersectionOptions,
) -> Result<Vec<SurfaceIntersectionCurve>> {
    let leaf = seed_leaf_extent(a, b);
    let pad = 1e-7;
    let box_pairs = seed_surface_surface(a, b, leaf, pad, 40)?;

    // Refine each candidate box-pair center to a point on the intersection.
    let mut seeds: Vec<Seed> = Vec::new();
    for (ba, bb) in box_pairs {
        let guess = Seed {
            ua: 0.5 * (ba.u0 + ba.u1),
            va: 0.5 * (ba.v0 + ba.v1),
            ub: 0.5 * (bb.u0 + bb.u1),
            vb: 0.5 * (bb.v0 + bb.v1),
        };
        if let Some(s) = refine_seed(a, b, guess, options)? {
            push_unique_seed(&mut seeds, a, s, leaf);
        }
    }

    let step = (leaf * options.step_factor).max(1e-6);

    let mut branches: Vec<SurfaceIntersectionCurve> = Vec::new();
    for seed in seeds {
        if seed_covered(&branches, a, seed, step)? {
            continue;
        }
        if let Some(branch) = march_branch(a, b, seed, step, options)? {
            branches.push(branch);
        }
    }
    Ok(branches)
}

/// A point on the intersection, parameterized on both surfaces.
#[derive(Debug, Clone, Copy)]
struct Seed {
    ua: f64,
    va: f64,
    ub: f64,
    vb: f64,
}

/// Leaf-extent heuristic: 1/12 of the larger control-hull diagonal.
fn seed_leaf_extent(a: &NurbsSurface, b: &NurbsSurface) -> f64 {
    let (a_lo, a_hi) = a.bounding_box();
    let (b_lo, b_hi) = b.bounding_box();
    let da = (a_hi - a_lo).norm();
    let db = (b_hi - b_lo).norm();
    (da.max(db) / 12.0).max(1e-6)
}

/// Damped Gauss-Newton refinement of a box-pair center onto `S_a = S_b`.
///
/// Minimizes `|S_a - S_b|²`: residual `r = S_a - S_b` (3-vector), Jacobian
/// `J = [S_a_u, S_a_v, -S_b_u, -S_b_v]` (3×4). The normal equations
/// `(JᵀJ + λI) Δ = -Jᵀ r` are solved with Levenberg-Marquardt damping `λ` to
/// regularize the rank-3 system (the null direction is the intersection
/// tangent). Returns `None` if the residual stays above tolerance.
#[allow(clippy::many_single_char_names, clippy::similar_names)]
fn refine_seed(
    a: &NurbsSurface,
    b: &NurbsSurface,
    guess: Seed,
    options: &IntersectionOptions,
) -> Result<Option<Seed>> {
    let ((au0, au1), (av0, av1)) = a.parameter_domain();
    let ((bu0, bu1), (bv0, bv1)) = b.parameter_domain();
    let mut s = guess;
    let tol = options.tolerance.max(1e-12);

    for _ in 0..options.max_iterations {
        let (pa, sau, sav) = a.partials(s.ua, s.va)?;
        let (pb, sbu, sbv) = b.partials(s.ub, s.vb)?;
        let r = pa - pb;
        if r.norm() < tol {
            return Ok(Some(s));
        }
        // J columns (3x4): d r / d(ua, va, ub, vb).
        let cols = [sau, sav, -sbu, -sbv];
        // Build JᵀJ (4x4) and Jᵀr (4).
        let mut jtj = [[0.0_f64; 4]; 4];
        let mut jtr = [0.0_f64; 4];
        for i in 0..4 {
            jtr[i] = cols[i].dot(&r);
            for j in 0..4 {
                jtj[i][j] = cols[i].dot(&cols[j]);
            }
        }
        // Levenberg-Marquardt damping.
        let lambda = 1e-9 * (jtj[0][0] + jtj[1][1] + jtj[2][2] + jtj[3][3]).max(1.0);
        let mut m = Matrix4::from_fn(|i, j| jtj[i][j]);
        for i in 0..4 {
            m[(i, i)] += lambda;
        }
        let rhs = Vector4::new(-jtr[0], -jtr[1], -jtr[2], -jtr[3]);
        let Some(minv) = m.try_inverse() else {
            break;
        };
        let delta = minv * rhs;
        s.ua = (s.ua + delta[0]).clamp(au0, au1);
        s.va = (s.va + delta[1]).clamp(av0, av1);
        s.ub = (s.ub + delta[2]).clamp(bu0, bu1);
        s.vb = (s.vb + delta[3]).clamp(bv0, bv1);
        if delta.norm() < tol {
            break;
        }
    }
    let pa = a.point_at(s.ua, s.va)?;
    let pb = b.point_at(s.ub, s.vb)?;
    if (pa - pb).norm() < tol.max(1e-7) {
        Ok(Some(s))
    } else {
        Ok(None)
    }
}

/// Inserts a seed unless an existing one is within `leaf` (3D) of it.
fn push_unique_seed(seeds: &mut Vec<Seed>, a: &NurbsSurface, candidate: Seed, leaf: f64) {
    if let Ok(p) = a.point_at(candidate.ua, candidate.va) {
        for s in seeds.iter() {
            if let Ok(q) = a.point_at(s.ua, s.va) {
                if (p - q).norm() < leaf * 0.5 {
                    return;
                }
            }
        }
        seeds.push(candidate);
    }
}

/// Whether a seed already lies on an existing branch (within one marching step).
///
/// Assumption: distinct intersection branches are separated by more than one
/// marching step; branches closer than that are merged into one. Tighten the
/// `step_factor` in [`IntersectionOptions`] when intersecting near-coincident
/// sheets.
fn seed_covered(
    branches: &[SurfaceIntersectionCurve],
    a: &NurbsSurface,
    seed: Seed,
    step: f64,
) -> Result<bool> {
    let p = a.point_at(seed.ua, seed.va)?;
    for branch in branches {
        for q in &branch.points {
            if (p - q).norm() < step {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// Marches a full branch from `seed`: forward, then (if open) backward, joined.
fn march_branch(
    a: &NurbsSurface,
    b: &NurbsSurface,
    seed: Seed,
    step: f64,
    options: &IntersectionOptions,
) -> Result<Option<SurfaceIntersectionCurve>> {
    let half = options.max_points / 2;
    let (fwd, closed) = march_direction(a, b, seed, step, 1.0, half, options)?;
    if closed {
        let mut branch = assemble(a, &fwd)?;
        branch.closed = true;
        return Ok(Some(branch));
    }
    let (mut bwd, _) = march_direction(a, b, seed, step, -1.0, half, options)?;
    bwd.reverse();
    let mut trace = bwd;
    trace.extend_from_slice(&fwd[1.min(fwd.len())..]);
    if trace.len() < 2 {
        return Ok(None);
    }
    Ok(Some(assemble(a, &trace)?))
}

/// Marches in one direction collecting intersection points. Returns the trace
/// (starting at the seed) and whether the branch closed back onto the seed.
///
/// `sign = +1` walks forward along `d = n_a × n_b`; `-1` walks backward.
#[allow(clippy::many_single_char_names, clippy::similar_names)]
fn march_direction(
    a: &NurbsSurface,
    b: &NurbsSurface,
    seed: Seed,
    step: f64,
    sign: f64,
    max_points: usize,
    options: &IntersectionOptions,
) -> Result<(Vec<Seed>, bool)> {
    let mut pts = vec![seed];
    let mut cur = seed;
    let mut prev_dir: Option<Vector3> = None;
    let seed_point = a.point_at(seed.ua, seed.va)?;

    for _ in 0..max_points {
        // Marching direction d = n_a × n_b at the current point.
        let Some((p, d)) = march_dir(a, b, cur)? else {
            // |d| ~ 0: tangential contact — stop, keeping what exists.
            break;
        };
        let mut dir = sign * d;
        if let Some(pd) = prev_dir {
            if dir.dot(&pd) < 0.0 {
                dir = -dir;
            }
        }
        let anchor = p + dir * step;

        // Corrector: 4x4 Newton onto S_a = S_b plus the marching-plane constraint.
        let Some(next) = correct(a, b, cur, anchor, dir, options)? else {
            break;
        };

        let np = a.point_at(next.ua, next.va)?;

        // Boundary termination: either surface's parameters landed on a domain
        // edge (the predictor walked out and the clamp pinned it there).
        //
        // Limitation: surfaces that are geometrically closed but parametrically
        // non-periodic (e.g. an extruded full circle with a UV seam) terminate
        // here at the seam instead of wrapping across it, so `closed` stays
        // false even though the branch geometrically returns to its start.
        // Consumers needing seam-aware closure must reclassify the open branch
        // (see `operations::boolean::nurbs::loops`). Future work: periodic-domain
        // wrapping in the marcher.
        let on_boundary = at_boundary(a, next.ua, next.va) || at_boundary(b, next.ub, next.vb);

        // Closure.
        if pts.len() > 3 && (np - seed_point).norm() < step * 0.75 {
            return Ok((pts, true));
        }

        prev_dir = Some(dir);
        pts.push(next);
        cur = next;

        if on_boundary {
            break;
        }
    }
    Ok((pts, false))
}

/// Current point on `a` and the unit marching direction `d = n_a × n_b`.
/// Returns `None` when the surfaces are tangent (`|d|` below tolerance).
fn march_dir(a: &NurbsSurface, b: &NurbsSurface, s: Seed) -> Result<Option<(Point3, Vector3)>> {
    let (pa, sau, sav) = a.partials(s.ua, s.va)?;
    let (_pb, sbu, sbv) = b.partials(s.ub, s.vb)?;
    let na = sau.cross(&sav);
    let nb = sbu.cross(&sbv);
    let d = na.cross(&nb);
    let dn = d.norm();
    if dn < 1e-9 {
        return Ok(None);
    }
    Ok(Some((pa, d / dn)))
}

/// 4×4 Newton corrector: solve `S_a = S_b` and `(S_a - anchor)·dir = 0`.
///
/// Unknowns `(ua, va, ub, vb)`; residual rows 1-3 = `S_a - S_b`, row 4 =
/// `(S_a - anchor)·dir`. Jacobian rows 1-3 = `[S_a_u, S_a_v, -S_b_u, -S_b_v]`,
/// row 4 = `[dir·S_a_u, dir·S_a_v, 0, 0]`. Clamped to both domains.
#[allow(clippy::many_single_char_names, clippy::similar_names)]
fn correct(
    a: &NurbsSurface,
    b: &NurbsSurface,
    start: Seed,
    anchor: Point3,
    dir: Vector3,
    options: &IntersectionOptions,
) -> Result<Option<Seed>> {
    let ((au0, au1), (av0, av1)) = a.parameter_domain();
    let ((bu0, bu1), (bv0, bv1)) = b.parameter_domain();
    let mut s = start;
    let tol = options.tolerance.max(1e-12);

    for _ in 0..options.max_iterations {
        let (pa, sau, sav) = a.partials(s.ua, s.va)?;
        let (pb, sbu, sbv) = b.partials(s.ub, s.vb)?;
        let r3 = pa - pb;
        let r4 = (pa - anchor).dot(&dir);
        if r3.norm() < tol && r4.abs() < tol {
            return Ok(Some(s));
        }
        let f = Vector4::new(r3.x, r3.y, r3.z, r4);
        // Jacobian columns.
        let j = Matrix4::new(
            sau.x,
            sav.x,
            -sbu.x,
            -sbv.x, //
            sau.y,
            sav.y,
            -sbu.y,
            -sbv.y, //
            sau.z,
            sav.z,
            -sbu.z,
            -sbv.z, //
            dir.dot(&sau),
            dir.dot(&sav),
            0.0,
            0.0,
        );
        let Some(jinv) = j.try_inverse() else {
            break;
        };
        let delta = jinv * f;
        s.ua = (s.ua - delta[0]).clamp(au0, au1);
        s.va = (s.va - delta[1]).clamp(av0, av1);
        s.ub = (s.ub - delta[2]).clamp(bu0, bu1);
        s.vb = (s.vb - delta[3]).clamp(bv0, bv1);
        if delta.norm() < tol {
            break;
        }
    }
    let pa = a.point_at(s.ua, s.va)?;
    let pb = b.point_at(s.ub, s.vb)?;
    if (pa - pb).norm() < tol.max(1e-7) {
        Ok(Some(s))
    } else {
        Ok(None)
    }
}

/// Boundary-proximity tolerance in parameter space.
const BOUNDARY_EPS: f64 = 1e-7;

/// Whether `(u, v)` sits on a parameter-domain boundary of `surface`.
fn at_boundary(surface: &NurbsSurface, u: f64, v: f64) -> bool {
    let ((u0, u1), (v0, v1)) = surface.parameter_domain();
    (u - u0).abs() < BOUNDARY_EPS
        || (u1 - u).abs() < BOUNDARY_EPS
        || (v - v0).abs() < BOUNDARY_EPS
        || (v1 - v).abs() < BOUNDARY_EPS
}

/// Assembles an (open) branch from a parameter trace. A closed branch is built
/// by the caller setting `closed = true` on the returned value.
fn assemble(a: &NurbsSurface, trace: &[Seed]) -> Result<SurfaceIntersectionCurve> {
    let mut points = Vec::with_capacity(trace.len());
    let mut uv_a = Vec::with_capacity(trace.len());
    let mut uv_b = Vec::with_capacity(trace.len());
    for s in trace {
        points.push(a.point_at(s.ua, s.va)?);
        uv_a.push(Point2::new(s.ua, s.va));
        uv_b.push(Point2::new(s.ub, s.vb));
    }
    Ok(SurfaceIntersectionCurve {
        points,
        uv_a,
        uv_b,
        closed: false,
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::nurbs::{KnotVector, NurbsSurface};
    use crate::math::Point3;

    /// 2x2 bilinear patch in the z = 0 plane, spanning `[x_lo,x_hi]x[y_lo,y_hi]`.
    fn z0_patch(x_lo: f64, x_hi: f64, y_lo: f64, y_hi: f64) -> NurbsSurface {
        NurbsSurface::from_unweighted(
            vec![
                Point3::new(x_lo, y_lo, 0.0),
                Point3::new(x_lo, y_hi, 0.0),
                Point3::new(x_hi, y_lo, 0.0),
                Point3::new(x_hi, y_hi, 0.0),
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

    /// A vertical bilinear patch in the plane `x = x0`, spanning y in
    /// `[y_lo,y_hi]` and z in `[z_lo,z_hi]`.
    fn x_const_patch(x0: f64, y_lo: f64, y_hi: f64, z_lo: f64, z_hi: f64) -> NurbsSurface {
        NurbsSurface::from_unweighted(
            vec![
                Point3::new(x0, y_lo, z_lo),
                Point3::new(x0, y_lo, z_hi),
                Point3::new(x0, y_hi, z_lo),
                Point3::new(x0, y_hi, z_hi),
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

    /// Quarter-cylinder shell: rational quadratic quarter circle in u (XY plane,
    /// radius 1) extruded along +Z by `h` in v.
    fn quarter_cylinder(h: f64) -> NurbsSurface {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        NurbsSurface::new(
            vec![
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, h),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(1.0, 1.0, h),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, h),
            ],
            vec![1.0, 1.0, w, w, 1.0, 1.0],
            3,
            2,
            KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]).unwrap(),
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            2,
            1,
        )
        .unwrap()
    }

    #[test]
    fn two_crossing_planar_patches_recover_the_line() {
        // z=0 patch over [-1,2]x[-1,2] crossed by the x=0.5 vertical patch
        // (y in [-1,2], z in [-1,1]). Intersection is the line x=0.5, z=0,
        // y in [-1, 2].
        let a = z0_patch(-1.0, 2.0, -1.0, 2.0);
        let b = x_const_patch(0.5, -1.0, 2.0, -1.0, 1.0);
        let branches = intersect_surfaces(&a, &b, &IntersectionOptions::default()).unwrap();
        assert_eq!(
            branches.len(),
            1,
            "expected one branch, got {}",
            branches.len()
        );
        let br = &branches[0];
        assert!(!br.closed);
        assert!(br.points.len() >= 2);
        for p in &br.points {
            assert!((p.x - 0.5).abs() < 1e-9, "x off line: {p:?}");
            assert!(p.z.abs() < 1e-9, "z off line: {p:?}");
            assert!(
                p.y >= -1.0 - 1e-7 && p.y <= 2.0 + 1e-7,
                "y out of range: {p:?}"
            );
        }
        // Endpoints reach the y-extent of the overlap region [-1, 2].
        let ys: Vec<f64> = br.points.iter().map(|p| p.y).collect();
        let ymin = ys.iter().copied().fold(f64::INFINITY, f64::min);
        let ymax = ys.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        assert!(ymin < -1.0 + 1e-3, "branch should reach y=-1, got {ymin}");
        assert!(ymax > 2.0 - 1e-3, "branch should reach y=2, got {ymax}");
    }

    #[test]
    fn quarter_cylinder_crossed_by_z0_patch_is_unit_arc() {
        // Quarter cylinder (x^2+y^2=1, z in [0,2]) crossed by the z=0 plane piece
        // over [-1,2]^2. Intersection is the quarter arc x^2+y^2=1 at z=0.
        let cyl = quarter_cylinder(2.0);
        let plane = z0_patch(-1.0, 2.0, -1.0, 2.0);
        let branches = intersect_surfaces(&cyl, &plane, &IntersectionOptions::default()).unwrap();
        assert!(!branches.is_empty(), "expected at least one branch");
        let mut total = 0;
        for br in &branches {
            for p in &br.points {
                let radial = (p.x * p.x + p.y * p.y).sqrt();
                assert!((radial - 1.0).abs() < 1e-6, "off unit circle: {p:?}");
                assert!(p.z.abs() < 1e-6, "off z=0: {p:?}");
                total += 1;
            }
        }
        assert!(total >= 4, "too few intersection points: {total}");
    }

    #[test]
    fn two_quarter_cylinders_crossing_satisfy_both_equations() {
        // Cylinder A: x^2 + y^2 = 1, axis +Z, z in [0,2].
        // Cylinder B: same shell rotated so its axis is +X — i.e. y^2 + z^2 = 1,
        // built by swapping the roles of x and z in the control net, offset so a
        // genuine transversal crossing exists in the first octant.
        let ca = quarter_cylinder(2.0);
        // Build B as a quarter cylinder about the X axis: y^2 + z^2 = 1, x in [0,2].
        let w = std::f64::consts::FRAC_1_SQRT_2;
        let cb = NurbsSurface::new(
            vec![
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(2.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 1.0),
                Point3::new(2.0, 1.0, 1.0),
                Point3::new(0.0, 0.0, 1.0),
                Point3::new(2.0, 0.0, 1.0),
            ],
            vec![1.0, 1.0, w, w, 1.0, 1.0],
            3,
            2,
            KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]).unwrap(),
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            2,
            1,
        )
        .unwrap();
        let branches = intersect_surfaces(&ca, &cb, &IntersectionOptions::default()).unwrap();
        assert!(!branches.is_empty(), "expected a crossing branch");
        let mut count = 0;
        for br in &branches {
            for p in &br.points {
                let ra = (p.x * p.x + p.y * p.y).sqrt();
                let rb = (p.y * p.y + p.z * p.z).sqrt();
                assert!((ra - 1.0).abs() < 1e-6, "off cylinder A: {p:?} (r={ra})");
                assert!((rb - 1.0).abs() < 1e-6, "off cylinder B: {p:?} (r={rb})");
                count += 1;
            }
        }
        assert!(count >= 2, "too few points on the crossing: {count}");
    }

    #[test]
    fn disjoint_patches_yield_empty() {
        // The z = 0 patch and a parallel patch far above (z = 5): no crossing.
        let a = z0_patch(0.0, 1.0, 0.0, 1.0);
        let b = NurbsSurface::from_unweighted(
            vec![
                Point3::new(0.0, 0.0, 5.0),
                Point3::new(0.0, 1.0, 5.0),
                Point3::new(1.0, 0.0, 5.0),
                Point3::new(1.0, 1.0, 5.0),
            ],
            2,
            2,
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
            1,
        )
        .unwrap();
        let branches = intersect_surfaces(&a, &b, &IntersectionOptions::default()).unwrap();
        assert!(
            branches.is_empty(),
            "expected no branches, got {}",
            branches.len()
        );
    }

    #[test]
    fn identical_patches_terminate_within_guard() {
        // Two identical patches: a degenerate (coincident) case. Out of P4 scope
        // beyond clean termination — assert it returns within the runaway guard
        // (bounded max_points) without hanging, and document the behavior.
        let a = z0_patch(0.0, 1.0, 0.0, 1.0);
        let b = z0_patch(0.0, 1.0, 0.0, 1.0);
        let opts = IntersectionOptions {
            max_points: 200,
            ..IntersectionOptions::default()
        };
        let branches = intersect_surfaces(&a, &b, &opts).unwrap();
        // Coincident surfaces have no transversal intersection direction
        // (n_a × n_b = 0 everywhere), so marching stops immediately. Any output
        // must be bounded by the runaway guard.
        for br in &branches {
            assert!(
                br.points.len() <= opts.max_points,
                "branch exceeded runaway guard: {}",
                br.points.len()
            );
        }
    }
}
