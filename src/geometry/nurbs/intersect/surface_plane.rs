//! Surface×plane intersection via marching on the implicit field.
//!
//! The intersection of a surface `S(u, v)` with the plane through `origin` with
//! unit `normal` is the zero set of the scalar field
//!
//! ```text
//! f(u, v) = (S(u, v) - origin) · normal.
//! ```
//!
//! ## Pipeline
//!
//! 1. **Seed.** Evaluate `f` on a knot-aware `(u, v)` grid. Each grid *edge*
//!    whose endpoints straddle zero brackets one crossing; a 1D Newton along the
//!    edge refines it to a starting point on `f = 0`.
//! 2. **March.** From a seed, the level-curve tangent in parameter space is
//!    `(f_v, -f_u)` (perpendicular to `∇f`). The predictor steps along it, scaled
//!    so the 3D arc step is `step_factor` × a local feature scale; the corrector
//!    is a 1D Newton back onto `f = 0` along `∇f`.
//! 3. **Terminate.** A branch ends on the parameter-domain boundary (clamped, the
//!    boundary point emitted), on closure (return within tolerance of the start
//!    after enough points), or at the `max_points` runaway guard.
//!
//! `uv_a` holds the surface parameters; `uv_b` holds the point's own 2D
//! coordinates in an orthonormal in-plane basis of the cutting plane.
//!
//! Tangential contact (the surface grazing the plane over a region) is an
//! unsupported degeneracy: the field has no transversal zero crossing there, so
//! no branch is produced — the solver terminates cleanly.

use crate::error::{GeometryError, Result};
use crate::geometry::nurbs::NurbsSurface;
use crate::math::{Point2, Point3, Vector3};

use super::types::{IntersectionOptions, SurfaceIntersectionCurve};

/// Computes the intersection branches of a NURBS surface with an infinite plane.
///
/// Each branch is an ordered polyline carrying synchronized UV traces: `uv_a` on
/// the surface, `uv_b` in the plane's own orthonormal 2D frame. Branches that
/// return to their start are marked `closed`.
///
/// # Errors
///
/// Returns an error if `normal` is the zero vector, or if a surface evaluation
/// fails (e.g. a vanishing rational denominator).
pub fn intersect_surface_plane(
    surface: &NurbsSurface,
    plane_origin: Point3,
    plane_normal: Vector3,
    options: &IntersectionOptions,
) -> Result<Vec<SurfaceIntersectionCurve>> {
    let nlen = plane_normal.norm();
    if nlen < crate::math::TOLERANCE {
        return Err(GeometryError::ZeroVector.into());
    }
    let normal = plane_normal / nlen;
    // Orthonormal in-plane basis (e1, e2) for projecting points to plane 2D.
    let (e1, e2) = in_plane_basis(normal);

    let field = Field {
        surface,
        origin: plane_origin,
        normal,
    };

    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_domain();
    let feature = ((u_max - u_min) + (v_max - v_min)) * 0.5;
    // Step in parameter space, scaled by step_factor. Marching corrects each
    // point back onto f = 0, so the step controls density, not accuracy.
    let step = (feature * options.step_factor * 0.05).max(1e-6);

    // Seed: collect edge crossings on a knot-aware grid.
    let seeds = seed_crossings(&field, u_min, u_max, v_min, v_max, options)?;

    let mut branches: Vec<SurfaceIntersectionCurve> = Vec::new();
    for seed in seeds {
        // Skip seeds already traversed by an existing branch. A seed within one
        // marching step of any branch point lies on that branch — marching it
        // again only re-walks the same curve.
        if seed_covered(&branches, surface, seed, step)? {
            continue;
        }
        if let Some(branch) =
            march_branch(&field, surface, seed, step, &e1, &e2, plane_origin, options)?
        {
            branches.push(branch);
        }
    }
    Ok(branches)
}

/// The scalar field `f(u, v) = (S(u, v) - origin)·normal` and its gradient.
struct Field<'a> {
    surface: &'a NurbsSurface,
    origin: Point3,
    normal: Vector3,
}

impl Field<'_> {
    /// Field value at `(u, v)`.
    fn value(&self, u: f64, v: f64) -> Result<f64> {
        let p = self.surface.point_at(u, v)?;
        Ok((p - self.origin).dot(&self.normal))
    }

    /// Field value and its parameter-space gradient `(f, f_u, f_v)`.
    fn value_grad(&self, u: f64, v: f64) -> Result<(f64, f64, f64)> {
        let skl = self.surface.derivatives(u, v, 1)?;
        let f = (Point3::from(skl[0][0]) - self.origin).dot(&self.normal);
        let fu = skl[1][0].dot(&self.normal);
        let fv = skl[0][1].dot(&self.normal);
        Ok((f, fu, fv))
    }
}

/// A refined starting point on `f = 0`.
#[derive(Debug, Clone, Copy)]
struct Seed {
    u: f64,
    v: f64,
}

/// Builds a knot-aware sample grid, finds grid-edge sign changes of `f`, and
/// refines each to a point on `f = 0`.
fn seed_crossings(
    field: &Field,
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
    options: &IntersectionOptions,
) -> Result<Vec<Seed>> {
    let us = sample_axis(field.surface.knots_u(), u_min, u_max);
    let vs = sample_axis(field.surface.knots_v(), v_min, v_max);

    // Field values on the grid.
    let mut grid = vec![0.0; us.len() * vs.len()];
    for (iu, &u) in us.iter().enumerate() {
        for (iv, &v) in vs.iter().enumerate() {
            grid[iu * vs.len() + iv] = field.value(u, v)?;
        }
    }

    let mut seeds = Vec::new();
    // Horizontal edges (varying v at fixed u).
    for iu in 0..us.len() {
        for iv in 0..vs.len() - 1 {
            let a = grid[iu * vs.len() + iv];
            let b = grid[iu * vs.len() + iv + 1];
            if straddles_zero(a, b) {
                if let Some(s) = refine_edge_v(field, us[iu], vs[iv], vs[iv + 1], a, b, options)? {
                    seeds.push(s);
                }
            }
        }
    }
    // Vertical edges (varying u at fixed v).
    for iv in 0..vs.len() {
        for iu in 0..us.len() - 1 {
            let a = grid[iu * vs.len() + iv];
            let b = grid[(iu + 1) * vs.len() + iv];
            if straddles_zero(a, b) {
                if let Some(s) = refine_edge_u(field, vs[iv], us[iu], us[iu + 1], a, b, options)? {
                    seeds.push(s);
                }
            }
        }
    }
    Ok(seeds)
}

/// Number of interior samples per knot span; resolves curvature inside a span.
const SUBDIV: usize = 8;

/// Distinct knot values plus interior samples between consecutive knots, so the
/// grid resolves curvature inside each knot span.
// Sample-index to f64 conversions are exact for the small sample counts.
#[allow(clippy::cast_precision_loss)]
fn sample_axis(knots: &crate::geometry::nurbs::KnotVector, lo: f64, hi: f64) -> Vec<f64> {
    let mut breaks: Vec<f64> = Vec::new();
    for &k in knots.as_slice() {
        if k >= lo - crate::math::TOLERANCE
            && k <= hi + crate::math::TOLERANCE
            && breaks.last().is_none_or(|&p| (p - k).abs() > 1e-9)
        {
            breaks.push(k.clamp(lo, hi));
        }
    }
    if breaks.len() < 2 {
        breaks = vec![lo, hi];
    }
    let mut out = Vec::new();
    for w in breaks.windows(2) {
        let (a, b) = (w[0], w[1]);
        for i in 0..SUBDIV {
            let t = i as f64 / SUBDIV as f64;
            out.push(a + (b - a) * t);
        }
    }
    out.push(hi);
    out
}

/// Whether two field values bracket a zero (strict sign change; an exact zero
/// endpoint counts).
fn straddles_zero(a: f64, b: f64) -> bool {
    (a <= 0.0 && b >= 0.0) || (a >= 0.0 && b <= 0.0)
}

/// Refines a crossing along a v-edge at fixed `u` by 1D Newton (falling back to
/// bisection when the gradient is flat).
fn refine_edge_v(
    field: &Field,
    u: f64,
    v0: f64,
    v1: f64,
    f0: f64,
    f1: f64,
    options: &IntersectionOptions,
) -> Result<Option<Seed>> {
    let mut lo = v0;
    let mut hi = v1;
    let mut flo = f0;
    let mut fhi = f1;
    // Secant/bisection hybrid for robustness on the 1D edge.
    let mut v = if (fhi - flo).abs() > 1e-300 {
        v0 - f0 * (v1 - v0) / (f1 - f0)
    } else {
        0.5 * (v0 + v1)
    };
    for _ in 0..options.max_iterations {
        v = v.clamp(lo.min(hi), lo.max(hi));
        let (f, _fu, fv) = field.value_grad(u, v)?;
        if f.abs() < options.tolerance {
            return Ok(Some(Seed { u, v }));
        }
        if (f <= 0.0) == (flo <= 0.0) {
            lo = v;
            flo = f;
        } else {
            hi = v;
            fhi = f;
        }
        let next = if fv.abs() > 1e-12 {
            v - f / fv
        } else {
            0.5 * (lo + hi)
        };
        if next < lo.min(hi) || next > lo.max(hi) {
            v = 0.5 * (lo + hi);
        } else {
            v = next;
        }
        let _ = fhi;
    }
    let (f, _, _) = field.value_grad(u, v)?;
    if f.abs() < options.tolerance.max(1e-7) {
        Ok(Some(Seed { u, v }))
    } else {
        Ok(None)
    }
}

/// Refines a crossing along a u-edge at fixed `v`.
fn refine_edge_u(
    field: &Field,
    v: f64,
    u0: f64,
    u1: f64,
    f0: f64,
    f1: f64,
    options: &IntersectionOptions,
) -> Result<Option<Seed>> {
    let mut lo = u0;
    let mut hi = u1;
    let mut flo = f0;
    let mut fhi = f1;
    let mut u = if (fhi - flo).abs() > 1e-300 {
        u0 - f0 * (u1 - u0) / (f1 - f0)
    } else {
        0.5 * (u0 + u1)
    };
    for _ in 0..options.max_iterations {
        u = u.clamp(lo.min(hi), lo.max(hi));
        let (f, fu, _fv) = field.value_grad(u, v)?;
        if f.abs() < options.tolerance {
            return Ok(Some(Seed { u, v }));
        }
        if (f <= 0.0) == (flo <= 0.0) {
            lo = u;
            flo = f;
        } else {
            hi = u;
            fhi = f;
        }
        let next = if fu.abs() > 1e-12 {
            u - f / fu
        } else {
            0.5 * (lo + hi)
        };
        if next < lo.min(hi) || next > lo.max(hi) {
            u = 0.5 * (lo + hi);
        } else {
            u = next;
        }
        let _ = fhi;
    }
    let (f, _, _) = field.value_grad(u, v)?;
    if f.abs() < options.tolerance.max(1e-7) {
        Ok(Some(Seed { u, v }))
    } else {
        Ok(None)
    }
}

/// Whether a seed already lies on an existing branch (within roughly one
/// marching step in 3D). Marched points are spaced `~step` apart, so a seed on
/// the same curve sits within `step` of the nearest branch point.
fn seed_covered(
    branches: &[SurfaceIntersectionCurve],
    surface: &NurbsSurface,
    seed: Seed,
    step: f64,
) -> Result<bool> {
    let p = surface.point_at(seed.u, seed.v)?;
    let cover = step;
    for branch in branches {
        for q in &branch.points {
            if (p - q).norm() < cover {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// Marches one full branch from `seed`, walking forward and (for open branches)
/// backward, joining the two half-traces.
#[allow(clippy::too_many_arguments)]
fn march_branch(
    field: &Field,
    surface: &NurbsSurface,
    seed: Seed,
    step: f64,
    e1: &Vector3,
    e2: &Vector3,
    origin: Point3,
    options: &IntersectionOptions,
) -> Result<Option<SurfaceIntersectionCurve>> {
    let half = options.max_points / 2;

    let (fwd, closed) = march_direction(field, seed, step, 1.0, half, options)?;
    if closed {
        return Ok(Some(assemble(surface, &fwd, true, e1, e2, origin)?));
    }
    // Open: also march backward from the seed and prepend.
    let (mut bwd, _) = march_direction(field, seed, step, -1.0, half, options)?;
    bwd.reverse();
    // bwd ends with the seed; fwd starts with the seed — drop the duplicate.
    let mut params = bwd;
    params.extend_from_slice(&fwd[1.min(fwd.len())..]);
    if params.len() < 2 {
        return Ok(None);
    }
    Ok(Some(assemble(surface, &params, false, e1, e2, origin)?))
}

/// Marches in one direction (`sign` = +1 forward, -1 backward) collecting
/// `(u, v)` points on `f = 0`. Returns the trace (starting at the seed) and
/// whether the branch closed back onto the seed.
// su/sv (surface partials) and the tangent bindings follow The NURBS Book.
#[allow(clippy::many_single_char_names, clippy::similar_names)]
fn march_direction(
    field: &Field,
    seed: Seed,
    step: f64,
    sign: f64,
    max_points: usize,
    options: &IntersectionOptions,
) -> Result<(Vec<(f64, f64)>, bool)> {
    let ((u_min, u_max), (v_min, v_max)) = field.surface.parameter_domain();
    let mut pts = vec![(seed.u, seed.v)];
    let mut u = seed.u;
    let mut v = seed.v;
    let mut prev_dir: Option<(f64, f64)> = None;

    for _ in 0..max_points {
        let (_, fu, fv) = field.value_grad(u, v)?;
        let gnorm = (fu * fu + fv * fv).sqrt();
        if gnorm < 1e-12 {
            // Flat field: cannot define a marching tangent — stop.
            break;
        }
        // Level-curve tangent in parameter space: perpendicular to ∇f.
        let mut tu = fv;
        let mut tv = -fu;
        let tnorm = (tu * tu + tv * tv).sqrt();
        tu /= tnorm;
        tv /= tnorm;
        // Keep a consistent walking orientation.
        let mut dir = (sign * tu, sign * tv);
        if let Some((pu, pv)) = prev_dir {
            if dir.0 * pu + dir.1 * pv < 0.0 {
                dir = (-dir.0, -dir.1);
            }
        }
        // Predictor: parameter step scaled so the 3D arc step ≈ `step`.
        let (_p, su, sv) = field.surface.partials(u, v)?;
        let speed = (su * dir.0 + sv * dir.1).norm().max(1e-9);
        let h = step / speed;
        let mut nu = u + dir.0 * h;
        let mut nv = v + dir.1 * h;

        let mut on_boundary = false;
        if nu < u_min || nu > u_max || nv < v_min || nv > v_max {
            on_boundary = true;
            nu = nu.clamp(u_min, u_max);
            nv = nv.clamp(v_min, v_max);
        }

        // Corrector: Newton along ∇f back onto f = 0.
        match correct(field, nu, nv, u_min, u_max, v_min, v_max, options)? {
            Some((cu, cv)) => {
                nu = cu;
                nv = cv;
            }
            None => break,
        }

        // Closure: returned near the seed after a few points.
        if pts.len() > 3 {
            let sp = field.surface.point_at(seed.u, seed.v)?;
            let np = field.surface.point_at(nu, nv)?;
            if (np - sp).norm() < step * 0.75 {
                return Ok((pts, true));
            }
        }

        prev_dir = Some(dir);
        pts.push((nu, nv));
        u = nu;
        v = nv;

        if on_boundary {
            break;
        }
    }
    Ok((pts, false))
}

/// Newton corrector: drives `f(u, v)` to zero by stepping along `∇f`, clamped to
/// the parameter domain. Returns `None` if it cannot reach `f = 0`.
#[allow(clippy::too_many_arguments)]
fn correct(
    field: &Field,
    u0: f64,
    v0: f64,
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
    options: &IntersectionOptions,
) -> Result<Option<(f64, f64)>> {
    let mut u = u0;
    let mut v = v0;
    let tol = options.tolerance.max(1e-12);
    for _ in 0..options.max_iterations {
        let (f, fu, fv) = field.value_grad(u, v)?;
        if f.abs() < tol {
            return Ok(Some((u, v)));
        }
        let g2 = fu * fu + fv * fv;
        if g2 < 1e-20 {
            return Ok(None);
        }
        // Gauss-Newton step along the gradient: minimal-norm solve of
        // f + ∇f·Δ = 0.
        let scale = -f / g2;
        u = (u + fu * scale).clamp(u_min, u_max);
        v = (v + fv * scale).clamp(v_min, v_max);
    }
    let (f, _, _) = field.value_grad(u, v)?;
    if f.abs() < tol.max(1e-7) {
        Ok(Some((u, v)))
    } else {
        Ok(None)
    }
}

/// Assembles a `SurfaceIntersectionCurve` from a `(u, v)` trace.
fn assemble(
    surface: &NurbsSurface,
    params: &[(f64, f64)],
    closed: bool,
    e1: &Vector3,
    e2: &Vector3,
    origin: Point3,
) -> Result<SurfaceIntersectionCurve> {
    let mut points = Vec::with_capacity(params.len());
    let mut uv_a = Vec::with_capacity(params.len());
    let mut uv_b = Vec::with_capacity(params.len());
    for &(u, v) in params {
        let p = surface.point_at(u, v)?;
        let d = p - origin;
        points.push(p);
        uv_a.push(Point2::new(u, v));
        uv_b.push(Point2::new(d.dot(e1), d.dot(e2)));
    }
    Ok(SurfaceIntersectionCurve {
        points,
        uv_a,
        uv_b,
        closed,
    })
}

/// Builds an orthonormal basis `(e1, e2)` spanning the plane with unit `normal`.
fn in_plane_basis(normal: Vector3) -> (Vector3, Vector3) {
    // Choose the world axis least aligned with the normal to seed e1.
    let a = if normal.x.abs() <= normal.y.abs() && normal.x.abs() <= normal.z.abs() {
        Vector3::x()
    } else if normal.y.abs() <= normal.z.abs() {
        Vector3::y()
    } else {
        Vector3::z()
    };
    let e1 = (a - normal * a.dot(&normal)).normalize();
    let e2 = normal.cross(&e1);
    (e1, e2)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::nurbs::{KnotVector, NurbsSurface};
    use crate::math::{Point3, Vector3};

    /// Quarter-cylinder shell: rational quadratic quarter circle in u (XY plane,
    /// radius 1) extruded along +Z by 2 in v. Surface points lie on x^2+y^2=1.
    fn quarter_cylinder_patch() -> NurbsSurface {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        NurbsSurface::new(
            vec![
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 2.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(1.0, 1.0, 2.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 2.0),
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

    /// 2x2 bilinear patch at z = 0 spanning [0,2]x[0,2] in XY.
    fn bilinear_patch() -> NurbsSurface {
        NurbsSurface::from_unweighted(
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
        .unwrap()
    }

    /// Full closed cylinder (radius 1, height 2) as a rational quadratic in u
    /// (nine-point full circle pattern) extruded along +Z in v. Points lie on
    /// x^2 + y^2 = 1 for all (u, v).
    fn full_cylinder() -> NurbsSurface {
        // Nine control points (degree-2 full circle), weights alternating
        // 1, 1/sqrt(2). Stored u-major with nv = 2 (bottom z=0, top z=2).
        let w = std::f64::consts::FRAC_1_SQRT_2;
        let ring = [
            (1.0, 0.0, 1.0),
            (1.0, 1.0, w),
            (0.0, 1.0, 1.0),
            (-1.0, 1.0, w),
            (-1.0, 0.0, 1.0),
            (-1.0, -1.0, w),
            (0.0, -1.0, 1.0),
            (1.0, -1.0, w),
            (1.0, 0.0, 1.0),
        ];
        let mut pts = Vec::new();
        let mut weights = Vec::new();
        for &(x, y, wt) in &ring {
            // v = 0 (z = 0) then v = 1 (z = 2).
            pts.push(Point3::new(x, y, 0.0));
            weights.push(wt);
            pts.push(Point3::new(x, y, 2.0));
            weights.push(wt);
        }
        NurbsSurface::new(
            pts,
            weights,
            9,
            2,
            KnotVector::new(vec![
                0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
            ])
            .unwrap(),
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            2,
            1,
        )
        .unwrap()
    }

    #[test]
    fn quarter_cylinder_cut_by_z_plane_is_unit_arc() {
        let s = quarter_cylinder_patch();
        // Plane z = 1: cut is the quarter arc x^2 + y^2 = 1 at z = 1.
        let branches = intersect_surface_plane(
            &s,
            Point3::new(0.0, 0.0, 1.0),
            Vector3::z(),
            &IntersectionOptions::default(),
        )
        .unwrap();
        assert_eq!(branches.len(), 1, "expected one open arc branch");
        let b = &branches[0];
        assert!(!b.closed, "quarter arc is an open branch");
        assert!(b.points.len() >= 4, "too few marched points");
        for (k, p) in b.points.iter().enumerate() {
            let radial = (p.x * p.x + p.y * p.y).sqrt();
            assert!(
                (radial - 1.0).abs() < 1e-6,
                "point {k} {p:?} off unit circle (r={radial})"
            );
            assert!((p.z - 1.0).abs() < 1e-6, "point {k} {p:?} off z=1");
        }
        // UV traces in-domain.
        for uv in &b.uv_a {
            assert!((0.0..=1.0).contains(&uv.x) && (0.0..=1.0).contains(&uv.y));
        }
        // Endpoints on the v-boundaries of the arc sweep (u = 0 and u = 1).
        let first = b.uv_a.first().unwrap();
        let last = b.uv_a.last().unwrap();
        let mut us = [first.x, last.x];
        us.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!(us[0] < 1e-4, "branch should start at u=0, got {}", us[0]);
        assert!(
            us[1] > 1.0 - 1e-4,
            "branch should end at u=1, got {}",
            us[1]
        );
    }

    #[test]
    fn bilinear_patch_parallel_plane_above_is_empty() {
        let s = bilinear_patch();
        // Plane z = 1 is parallel to and above the z = 0 patch: no crossing.
        let branches = intersect_surface_plane(
            &s,
            Point3::new(0.0, 0.0, 1.0),
            Vector3::z(),
            &IntersectionOptions::default(),
        )
        .unwrap();
        assert!(
            branches.is_empty(),
            "expected no branches, got {}",
            branches.len()
        );
    }

    #[test]
    fn full_cylinder_cut_by_z_plane_is_closed_loop() {
        let s = full_cylinder();
        // Plane z = 1: cut is the full unit circle at z = 1 -> a closed loop.
        let branches = intersect_surface_plane(
            &s,
            Point3::new(0.0, 0.0, 1.0),
            Vector3::z(),
            &IntersectionOptions::default(),
        )
        .unwrap();
        assert_eq!(branches.len(), 1, "expected one closed branch");
        let b = &branches[0];
        assert!(b.closed, "full-circle cut must be a closed loop");
        for p in &b.points {
            let radial = (p.x * p.x + p.y * p.y).sqrt();
            assert!((radial - 1.0).abs() < 1e-6, "{p:?} off unit circle");
            assert!((p.z - 1.0).abs() < 1e-6, "{p:?} off z=1");
        }
    }

    #[test]
    fn zero_normal_is_error() {
        let s = bilinear_patch();
        let r = intersect_surface_plane(
            &s,
            Point3::origin(),
            Vector3::zeros(),
            &IntersectionOptions::default(),
        );
        assert!(r.is_err());
    }
}
