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
//! 3. **Terminate.** Open-direction domain boundary of either surface (clamp +
//!    emit), closure (return within tolerance of the start after >3 points), or
//!    the `max_points` runaway guard. Geometrically closed directions (e.g. the
//!    u-seam of an extruded full circle) are periodic: iterates wrap modulo the
//!    period instead of clamping, the seam does not terminate the branch, and
//!    loops on closed surfaces come out genuinely `closed`.
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

use nalgebra::{Matrix3, Matrix4, Vector4};

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
    let pair = SurfacePair {
        a,
        b,
        dom_a: SurfaceDomain::of(a),
        dom_b: SurfaceDomain::of(b),
    };

    // Refine each candidate box-pair center to a point on the intersection.
    let mut seeds: Vec<Seed> = Vec::new();
    for (ba, bb) in box_pairs {
        let guess = Seed {
            ua: 0.5 * (ba.u0 + ba.u1),
            va: 0.5 * (ba.v0 + ba.v1),
            ub: 0.5 * (bb.u0 + bb.u1),
            vb: 0.5 * (bb.v0 + bb.v1),
        };
        if let Some(s) = refine_seed(&pair, guess, options)? {
            push_unique_seed(&mut seeds, a, s, leaf);
        }
    }

    let step = (leaf * options.step_factor).max(1e-6);

    let mut branches: Vec<SurfaceIntersectionCurve> = Vec::new();
    for seed in seeds {
        if seed_covered(&branches, a, seed, step)? {
            continue;
        }
        if let Some(branch) = march_branch(&pair, seed, step, options)? {
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

/// The two surfaces under intersection plus their periodic-domain descriptors,
/// threaded through the seeding / marching helpers as one unit.
struct SurfacePair<'s> {
    a: &'s NurbsSurface,
    b: &'s NurbsSurface,
    dom_a: SurfaceDomain,
    dom_b: SurfaceDomain,
}

/// Parameter domain of one surface plus which directions are geometrically
/// closed (the parametric boundaries map to the same 3D seam curve).
///
/// Closed directions are periodic for the marcher: iterates wrap modulo the
/// period (the same convention as [`NurbsSurface::closest_point`]) and the
/// parametric seam does not terminate a branch. Open directions clamp and
/// terminate as before.
#[derive(Debug, Clone, Copy)]
struct SurfaceDomain {
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    u_closed: bool,
    v_closed: bool,
}

impl SurfaceDomain {
    fn of(s: &NurbsSurface) -> Self {
        let ((u0, u1), (v0, v1)) = s.parameter_domain();
        Self {
            u0,
            u1,
            v0,
            v1,
            u_closed: s.is_closed_in_u(),
            v_closed: s.is_closed_in_v(),
        }
    }

    /// Wraps closed directions modulo their period; clamps open ones.
    fn apply(&self, u: f64, v: f64) -> (f64, f64) {
        let wu = if self.u_closed {
            self.u0 + (u - self.u0).rem_euclid(self.u1 - self.u0)
        } else {
            u.clamp(self.u0, self.u1)
        };
        let wv = if self.v_closed {
            self.v0 + (v - self.v0).rem_euclid(self.v1 - self.v0)
        } else {
            v.clamp(self.v0, self.v1)
        };
        (wu, wv)
    }

    /// Boundary proximity that terminates marching: OPEN directions only. A
    /// closed direction's parametric seam is not a geometric boundary.
    fn at_open_boundary(&self, u: f64, v: f64) -> bool {
        (!self.u_closed
            && ((u - self.u0).abs() < BOUNDARY_EPS || (self.u1 - u).abs() < BOUNDARY_EPS))
            || (!self.v_closed
                && ((v - self.v0).abs() < BOUNDARY_EPS || (self.v1 - v).abs() < BOUNDARY_EPS))
    }
}

/// Leaf-extent heuristic: 1/12 of the SMALLER control-hull diagonal.
///
/// The intersection curve lies on both surfaces, so its geometric extent (and
/// the feature scale the marcher must resolve) is bounded by the smaller
/// surface. Scaling by the larger diagonal made the marching step blow up when
/// a small tool cuts a large target (e.g. a window prism through a long curved
/// wall): the step became comparable to the tool's corner radii, the corrector
/// jumped between loop lobes, and branches fragmented into overlapping open
/// pieces.
fn seed_leaf_extent(a: &NurbsSurface, b: &NurbsSurface) -> f64 {
    let (a_lo, a_hi) = a.bounding_box();
    let (b_lo, b_hi) = b.bounding_box();
    let da = (a_hi - a_lo).norm();
    let db = (b_hi - b_lo).norm();
    (da.min(db) / 12.0).max(1e-6)
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
    pair: &SurfacePair,
    guess: Seed,
    options: &IntersectionOptions,
) -> Result<Option<Seed>> {
    let mut s = guess;
    let tol = options.tolerance.max(1e-12);

    for _ in 0..options.max_iterations {
        let (pa, sau, sav) = pair.a.partials(s.ua, s.va)?;
        let (pb, sbu, sbv) = pair.b.partials(s.ub, s.vb)?;
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
        (s.ua, s.va) = pair.dom_a.apply(s.ua + delta[0], s.va + delta[1]);
        (s.ub, s.vb) = pair.dom_b.apply(s.ub + delta[2], s.vb + delta[3]);
        if delta.norm() < tol {
            break;
        }
    }
    let pa = pair.a.point_at(s.ua, s.va)?;
    let pb = pair.b.point_at(s.ub, s.vb)?;
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
    pair: &SurfacePair,
    seed: Seed,
    step: f64,
    options: &IntersectionOptions,
) -> Result<Option<SurfaceIntersectionCurve>> {
    let half = options.max_points / 2;
    let (fwd, closed) = march_direction(pair, seed, step, 1.0, half, options)?;
    if closed {
        let mut branch = assemble(pair.a, &fwd)?;
        branch.closed = true;
        return Ok(Some(branch));
    }
    let (mut bwd, _) = march_direction(pair, seed, step, -1.0, half, options)?;
    bwd.reverse();
    let mut trace = bwd;
    trace.extend_from_slice(&fwd[1.min(fwd.len())..]);
    if trace.len() < 2 {
        return Ok(None);
    }
    Ok(Some(assemble(pair.a, &trace)?))
}

/// Marches in one direction collecting intersection points. Returns the trace
/// (starting at the seed) and whether the branch closed back onto the seed.
///
/// `sign = +1` walks forward along `d = n_a × n_b`; `-1` walks backward.
#[allow(clippy::many_single_char_names, clippy::similar_names)]
fn march_direction(
    pair: &SurfacePair,
    seed: Seed,
    step: f64,
    sign: f64,
    max_points: usize,
    options: &IntersectionOptions,
) -> Result<(Vec<Seed>, bool)> {
    let mut pts = vec![seed];
    let mut cur = seed;
    let mut prev_dir: Option<Vector3> = None;
    let seed_point = pair.a.point_at(seed.ua, seed.va)?;

    for _ in 0..max_points {
        // Marching direction d = n_a × n_b at the current point.
        let Some((p, d)) = march_dir(pair, cur)? else {
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
        let next = match correct(pair, cur, anchor, dir, options)? {
            Corrected::Point(s) => s,
            // The Newton iterates kept escaping an OPEN parametric direction:
            // the branch leaves that domain between the current point and the
            // clamped bound, so the marching-plane system is infeasible there.
            // Solve the exact boundary crossing instead (the same pinned 3×3
            // Newton the seam samples use) so the branch endpoint lands ON the
            // boundary rather than up to one marching step short of it.
            Corrected::ClampedOut { param, bound } => {
                if let Some(end) = pinned_correct(pair, param, bound, cur, options)? {
                    let end_p = pair.a.point_at(end.ua, end.va)?;
                    // Sanity bound: the crossing lies within the current
                    // marching neighborhood. A far-away pinned solution means
                    // Newton escaped toward an unrelated boundary; keep the
                    // honest sub-step gap instead of fabricating a long edge.
                    if (end_p - p).norm() <= 2.0 * step {
                        if let Some(crossing) = seam_crossing(pair, cur, end) {
                            if let Some((near, far)) = seam_samples(pair, crossing, cur, options)? {
                                pts.push(near);
                                pts.push(far);
                            }
                        }
                        pts.push(end);
                    }
                }
                break;
            }
            Corrected::Lost => break,
        };

        let np = pair.a.point_at(next.ua, next.va)?;

        // Boundary termination: an OPEN direction's parameters landed on a
        // domain edge (the predictor walked out and the clamp pinned it there).
        // Geometrically closed directions wrap across their parametric seam
        // instead (see `SurfaceDomain::apply`), so the march continues until
        // 3D closure.
        let on_boundary = pair.dom_a.at_open_boundary(next.ua, next.va)
            || pair.dom_b.at_open_boundary(next.ub, next.vb);

        // The corrector wrapped a closed direction across its seam: insert the
        // exact seam point in BOTH parametric representations (same 3D point at
        // the low and high boundary) so downstream consumers sorting by the
        // wrapped parameter span the full domain with no wedge gap. Inserted
        // before the closure check so a crossing on the closing step is kept.
        if let Some(crossing) = seam_crossing(pair, cur, next) {
            if let Some((near, far)) = seam_samples(pair, crossing, cur, options)? {
                pts.push(near);
                pts.push(far);
            }
        }

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

/// Which of the four marching parameters crossed a closed direction's seam.
///
/// The value is the parameter's index in the `(ua, va, ub, vb)` layout shared
/// with the corrector Jacobian columns.
#[derive(Debug, Clone, Copy)]
struct SeamCrossing {
    /// Parameter index: 0 = `ua`, 1 = `va`, 2 = `ub`, 3 = `vb`.
    param: usize,
    /// Domain low / high bounds of the crossed direction.
    lo: f64,
    hi: f64,
}

/// Seed parameters as an indexable `(ua, va, ub, vb)` array.
fn seed_params(s: Seed) -> [f64; 4] {
    [s.ua, s.va, s.ub, s.vb]
}

/// Rebuilds a `Seed` from the indexable parameter layout.
fn seed_from_params(p: [f64; 4]) -> Seed {
    Seed {
        ua: p[0],
        va: p[1],
        ub: p[2],
        vb: p[3],
    }
}

/// Detects a seam crossing between two consecutive marched points: a closed
/// direction whose wrapped step is longer than half its period. Returns the
/// first crossing found (a simultaneous crossing of two directions in one
/// step is a corner-degenerate rarity; its second direction simply keeps a
/// sub-step gap).
fn seam_crossing(pair: &SurfacePair, prev: Seed, next: Seed) -> Option<SeamCrossing> {
    let (dom_a, dom_b) = (&pair.dom_a, &pair.dom_b);
    let closed = [
        (dom_a.u_closed, dom_a.u0, dom_a.u1),
        (dom_a.v_closed, dom_a.v0, dom_a.v1),
        (dom_b.u_closed, dom_b.u0, dom_b.u1),
        (dom_b.v_closed, dom_b.v0, dom_b.v1),
    ];
    let p = seed_params(prev);
    let n = seed_params(next);
    for (param, &(is_closed, lo, hi)) in closed.iter().enumerate() {
        if is_closed && (n[param] - p[param]).abs() > 0.5 * (hi - lo) {
            return Some(SeamCrossing { param, lo, hi });
        }
    }
    None
}

/// Solves the exact seam point for a crossing and returns it in both
/// parametric representations: `(near, far)` where `near` pins the crossed
/// parameter at the boundary closer to `prev` (trace order) and `far` is the
/// identical 3D point re-expressed at the opposite boundary.
///
/// One pinned 3×3 Newton solves `S_a = S_b` over the three non-pinned
/// parameters; because the direction is geometrically closed, the solution at
/// one boundary IS the solution at the other, so no second solve is needed.
/// Returns `None` on non-convergence — the caller inserts nothing and the
/// trace keeps an honest sub-step gap at the seam.
fn seam_samples(
    pair: &SurfacePair,
    crossing: SeamCrossing,
    prev: Seed,
    options: &IntersectionOptions,
) -> Result<Option<(Seed, Seed)>> {
    let prev_p = seed_params(prev);
    let mid = 0.5 * (crossing.lo + crossing.hi);
    let (near_bound, far_bound) = if prev_p[crossing.param] >= mid {
        (crossing.hi, crossing.lo)
    } else {
        (crossing.lo, crossing.hi)
    };
    let Some(pinned) = pinned_correct(pair, crossing.param, near_bound, prev, options)? else {
        return Ok(None);
    };
    let mut far = seed_params(pinned);
    far[crossing.param] = far_bound;
    Ok(Some((pinned, seed_from_params(far))))
}

/// Pinned 3×3 Newton: solve `S_a(u_a, v_a) = S_b(u_b, v_b)` with the parameter
/// at `pin_param` held at `pin_value` and the remaining three free.
///
/// Residual `r = S_a - S_b` (3 components); Jacobian columns are the three
/// non-pinned columns of `[S_a_u, S_a_v, -S_b_u, -S_b_v]`. Free iterates wrap /
/// clamp via the domain descriptors; the pinned parameter is re-asserted after
/// every domain application (a closed direction's `apply` would wrap the high
/// boundary back to the low one). Returns `None` if the residual stays above
/// tolerance.
#[allow(clippy::similar_names)]
fn pinned_correct(
    pair: &SurfacePair,
    pin_param: usize,
    pin_value: f64,
    start: Seed,
    options: &IntersectionOptions,
) -> Result<Option<Seed>> {
    let tol = options.tolerance.max(1e-12);
    let mut p = seed_params(start);
    p[pin_param] = pin_value;
    let free: Vec<usize> = (0..4).filter(|&i| i != pin_param).collect();

    for _ in 0..options.max_iterations {
        let (pa, sau, sav) = pair.a.partials(p[0], p[1])?;
        let (pb, sbu, sbv) = pair.b.partials(p[2], p[3])?;
        let r = pa - pb;
        if r.norm() < tol {
            return Ok(Some(seed_from_params(p)));
        }
        let cols = [sau, sav, -sbu, -sbv];
        let j = Matrix3::from_columns(&[cols[free[0]], cols[free[1]], cols[free[2]]]);
        let Some(jinv) = j.try_inverse() else {
            break;
        };
        let delta = jinv * r;
        for (t, &i) in free.iter().enumerate() {
            p[i] -= delta[t];
        }
        (p[0], p[1]) = pair.dom_a.apply(p[0], p[1]);
        (p[2], p[3]) = pair.dom_b.apply(p[2], p[3]);
        p[pin_param] = pin_value;
        if delta.norm() < tol {
            break;
        }
    }
    let pa = pair.a.point_at(p[0], p[1])?;
    let pb = pair.b.point_at(p[2], p[3])?;
    if (pa - pb).norm() < tol.max(1e-7) {
        Ok(Some(seed_from_params(p)))
    } else {
        Ok(None)
    }
}

/// Current point on `a` and the unit marching direction `d = n_a × n_b`.
/// Returns `None` when the surfaces are tangent (`|d|` below tolerance).
fn march_dir(pair: &SurfacePair, s: Seed) -> Result<Option<(Point3, Vector3)>> {
    let (pa, sau, sav) = pair.a.partials(s.ua, s.va)?;
    let (_pb, sbu, sbv) = pair.b.partials(s.ub, s.vb)?;
    let na = sau.cross(&sav);
    let nb = sbu.cross(&sbv);
    let d = na.cross(&nb);
    let dn = d.norm();
    if dn < 1e-9 {
        return Ok(None);
    }
    Ok(Some((pa, d / dn)))
}

/// Outcome of the marching corrector.
#[derive(Debug, Clone, Copy)]
enum Corrected {
    /// Converged onto the intersection at the marching plane.
    Point(Seed),
    /// The Newton iterates repeatedly pushed an OPEN parametric direction past
    /// a domain bound (the clamp kept pinning it back): the intersection curve
    /// exits the domain within this step. `param` is the escaping parameter's
    /// index in the `(ua, va, ub, vb)` layout; `bound` the bound it hit.
    ClampedOut { param: usize, bound: f64 },
    /// No convergence and no clamping — tangency or runaway.
    Lost,
}

/// 4×4 Newton corrector: solve `S_a = S_b` and `(S_a - anchor)·dir = 0`.
///
/// Unknowns `(ua, va, ub, vb)`; residual rows 1-3 = `S_a - S_b`, row 4 =
/// `(S_a - anchor)·dir`. Jacobian rows 1-3 = `[S_a_u, S_a_v, -S_b_u, -S_b_v]`,
/// row 4 = `[dir·S_a_u, dir·S_a_v, 0, 0]`. Iterates wrap on closed directions
/// and clamp on open ones (see [`SurfaceDomain::apply`]); the last open-domain
/// clamp is recorded so a divergence caused by the curve leaving a domain is
/// reported as [`Corrected::ClampedOut`] instead of a silent failure.
#[allow(clippy::many_single_char_names, clippy::similar_names)]
fn correct(
    pair: &SurfacePair,
    start: Seed,
    anchor: Point3,
    dir: Vector3,
    options: &IntersectionOptions,
) -> Result<Corrected> {
    let mut s = start;
    let tol = options.tolerance.max(1e-12);
    let mut last_clamp: Option<(usize, f64)> = None;

    for _ in 0..options.max_iterations {
        let (pa, sau, sav) = pair.a.partials(s.ua, s.va)?;
        let (pb, sbu, sbv) = pair.b.partials(s.ub, s.vb)?;
        let r3 = pa - pb;
        let r4 = (pa - anchor).dot(&dir);
        if r3.norm() < tol && r4.abs() < tol {
            return Ok(Corrected::Point(s));
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
        let raw = [
            s.ua - delta[0],
            s.va - delta[1],
            s.ub - delta[2],
            s.vb - delta[3],
        ];
        (s.ua, s.va) = pair.dom_a.apply(raw[0], raw[1]);
        (s.ub, s.vb) = pair.dom_b.apply(raw[2], raw[3]);
        if let Some(clamp) = detect_open_clamp(pair, raw) {
            last_clamp = Some(clamp);
        }
        if delta.norm() < tol {
            break;
        }
    }
    let pa = pair.a.point_at(s.ua, s.va)?;
    let pb = pair.b.point_at(s.ub, s.vb)?;
    if (pa - pb).norm() < tol.max(1e-7) {
        Ok(Corrected::Point(s))
    } else if let Some((param, bound)) = last_clamp {
        Ok(Corrected::ClampedOut { param, bound })
    } else {
        Ok(Corrected::Lost)
    }
}

/// Detects which OPEN parametric direction the last Newton update escaped:
/// a raw (pre-domain) parameter value outside its open domain means the clamp
/// pinned it at the bound it crossed. Wrapping of closed directions is not a
/// clamp. Returns the parameter index in the `(ua, va, ub, vb)` layout and
/// the bound that pinned it.
fn detect_open_clamp(pair: &SurfacePair, raw: [f64; 4]) -> Option<(usize, f64)> {
    let domains = [
        (!pair.dom_a.u_closed, pair.dom_a.u0, pair.dom_a.u1),
        (!pair.dom_a.v_closed, pair.dom_a.v0, pair.dom_a.v1),
        (!pair.dom_b.u_closed, pair.dom_b.u0, pair.dom_b.u1),
        (!pair.dom_b.v_closed, pair.dom_b.v0, pair.dom_b.v1),
    ];
    for (param, &(open, lo, hi)) in domains.iter().enumerate() {
        if !open {
            continue;
        }
        if raw[param] < lo {
            return Some((param, lo));
        }
        if raw[param] > hi {
            return Some((param, hi));
        }
    }
    None
}

/// Boundary-proximity tolerance in parameter space. Shared with the boolean
/// loop extraction so "this open branch endpoint sits on a parametric
/// boundary" means exactly what the marcher's own termination meant.
pub(crate) const BOUNDARY_EPS: f64 = 1e-7;

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
    use crate::geometry::nurbs::{KnotVector, NurbsCurve3D, NurbsSurface};
    use crate::math::Point3;

    /// 2x2 bilinear patch in the z = `z` plane, spanning `[x_lo,x_hi]x[y_lo,y_hi]`.
    fn z_patch(x_lo: f64, x_hi: f64, y_lo: f64, y_hi: f64, z: f64) -> NurbsSurface {
        NurbsSurface::from_unweighted(
            vec![
                Point3::new(x_lo, y_lo, z),
                Point3::new(x_lo, y_hi, z),
                Point3::new(x_hi, y_lo, z),
                Point3::new(x_hi, y_hi, z),
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

    /// 2x2 bilinear patch in the z = 0 plane, spanning `[x_lo,x_hi]x[y_lo,y_hi]`.
    fn z0_patch(x_lo: f64, x_hi: f64, y_lo: f64, y_hi: f64) -> NurbsSurface {
        z_patch(x_lo, x_hi, y_lo, y_hi, 0.0)
    }

    /// Full closed cylinder shell: the unit circle in the XY plane extruded `h`
    /// along +Z. Geometrically closed in u (`is_closed_in_u()` is true) with a
    /// parametric seam at the +X azimuth.
    fn full_cylinder(h: f64) -> NurbsSurface {
        let circle =
            NurbsCurve3D::circle(Point3::origin(), 1.0, Vector3::z(), Vector3::x()).unwrap();
        NurbsSurface::extrude(&circle, Vector3::new(0.0, 0.0, h)).unwrap()
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
    fn closed_cylinder_crossed_by_plane_patch_yields_one_closed_branch() {
        // Cylinder x²+y²=1, z∈[0,2] crossed by the z=1 plane patch over [-2,2]².
        // The marcher must wrap the cylinder's parametric u-seam and close the
        // branch instead of terminating there with an open trace.
        let cyl = full_cylinder(2.0);
        let plane = z_patch(-2.0, 2.0, -2.0, 2.0, 1.0);
        let branches = intersect_surfaces(&cyl, &plane, &IntersectionOptions::default()).unwrap();
        assert_eq!(
            branches.len(),
            1,
            "expected one branch, got {}",
            branches.len()
        );
        let br = &branches[0];
        assert!(br.closed, "wrap must close the branch across the seam");
        assert!(br.points.len() > 8, "too few points: {}", br.points.len());
        let ((u0, u1), _) = cyl.parameter_domain();
        for (p, uv) in br.points.iter().zip(&br.uv_a) {
            assert!(
                ((p.x * p.x + p.y * p.y).sqrt() - 1.0).abs() < 1e-6,
                "off unit circle: {p:?}"
            );
            assert!((p.z - 1.0).abs() < 1e-6, "off z=1: {p:?}");
            assert!(
                uv.x >= u0 - 1e-9 && uv.x <= u1 + 1e-9,
                "uv_a u out of domain: {}",
                uv.x
            );
        }
    }

    #[test]
    fn wrapped_branch_contains_exact_seam_samples() {
        // Downstream trim/band construction sorts closed-direction traces by the
        // wrapped parameter and needs samples at EXACTLY both parametric
        // boundaries (same 3D point, two representations) to span the full
        // domain without a wedge gap at the seam.
        let cyl = full_cylinder(2.0);
        let plane = z_patch(-2.0, 2.0, -2.0, 2.0, 1.0);
        let branches = intersect_surfaces(&cyl, &plane, &IntersectionOptions::default()).unwrap();
        assert_eq!(branches.len(), 1);
        let br = &branches[0];
        assert!(br.closed);
        let ((u0, u1), _) = cyl.parameter_domain();
        let at_lo: Vec<usize> = (0..br.uv_a.len())
            .filter(|&i| (br.uv_a[i].x - u0).abs() < 1e-9)
            .collect();
        let at_hi: Vec<usize> = (0..br.uv_a.len())
            .filter(|&i| (u1 - br.uv_a[i].x).abs() < 1e-9)
            .collect();
        assert!(
            !at_lo.is_empty() && !at_hi.is_empty(),
            "trace must touch both u boundaries exactly (lo hits: {}, hi hits: {})",
            at_lo.len(),
            at_hi.len()
        );
        // The two representations of one crossing share their 3D point.
        assert!(
            (br.points[at_lo[0]] - br.points[at_hi[0]]).norm() < 1e-9,
            "seam sample pair must coincide in 3D"
        );
    }

    #[test]
    fn seam_samples_satisfy_both_surface_equations() {
        // Every point of the closed branch (seam samples included) lies on BOTH
        // surfaces within 1e-6.
        let cyl = full_cylinder(2.0);
        let plane = z_patch(-2.0, 2.0, -2.0, 2.0, 1.0);
        let branches = intersect_surfaces(&cyl, &plane, &IntersectionOptions::default()).unwrap();
        let br = &branches[0];
        for (i, p) in br.points.iter().enumerate() {
            let pa = cyl.point_at(br.uv_a[i].x, br.uv_a[i].y).unwrap();
            let pb = plane.point_at(br.uv_b[i].x, br.uv_b[i].y).unwrap();
            assert!(
                (pa - p).norm() < 1e-6,
                "point {i} off surface a: {}",
                (pa - p).norm()
            );
            assert!(
                (pb - p).norm() < 1e-6,
                "point {i} off surface b: {}",
                (pb - p).norm()
            );
        }
    }

    #[test]
    fn closed_cylinder_as_surface_b_also_closes() {
        // Same crossing with the argument order swapped: the closed surface sits
        // in the tool position (`b`), so the wrap acts on the uv_b trace.
        let cyl = full_cylinder(2.0);
        let plane = z_patch(-2.0, 2.0, -2.0, 2.0, 1.0);
        let branches = intersect_surfaces(&plane, &cyl, &IntersectionOptions::default()).unwrap();
        assert_eq!(
            branches.len(),
            1,
            "expected one branch, got {}",
            branches.len()
        );
        assert!(branches[0].closed, "closed tool branch must close");
    }

    #[test]
    fn open_branch_endpoints_land_exactly_on_the_limiting_boundary() {
        // z=0 patch over [-1,2]x[-1,2] crossed by the x=0.5 vertical patch with
        // y limited to [0.25, 1.25] (surface `b` truncates the intersection
        // line; `a` extends beyond on both sides). The open branch must
        // terminate ON b's u boundaries — the boundary-pinned endpoint solve —
        // not one marching step short of them. This is the contract the
        // multi-face (box) tool chaining relies on: branch endpoints land on
        // the tool's kink edges exactly.
        let a = z0_patch(-1.0, 2.0, -1.0, 2.0);
        let b = x_const_patch(0.5, 0.25, 1.25, -1.0, 1.0);
        let branches = intersect_surfaces(&a, &b, &IntersectionOptions::default()).unwrap();
        assert_eq!(
            branches.len(),
            1,
            "expected one branch, got {}",
            branches.len()
        );
        let br = &branches[0];
        assert!(!br.closed);
        // Endpoints reach y = 0.25 and y = 1.25 exactly (b's u boundaries map
        // to those y values).
        let y_first = br.points.first().unwrap().y;
        let y_last = br.points.last().unwrap().y;
        let (y_lo, y_hi) = if y_first < y_last {
            (y_first, y_last)
        } else {
            (y_last, y_first)
        };
        assert!(
            (y_lo - 0.25).abs() < BOUNDARY_EPS,
            "low endpoint must land on b's boundary y=0.25, got {y_lo}"
        );
        assert!(
            (y_hi - 1.25).abs() < BOUNDARY_EPS,
            "high endpoint must land on b's boundary y=1.25, got {y_hi}"
        );
        // The uv_b endpoints sit at b's u domain bounds.
        let ((u0, u1), _) = b.parameter_domain();
        let ub_first = br.uv_b.first().unwrap().x;
        let ub_last = br.uv_b.last().unwrap().x;
        for ub in [ub_first, ub_last] {
            assert!(
                (ub - u0).abs() < BOUNDARY_EPS || (u1 - ub).abs() < BOUNDARY_EPS,
                "uv_b endpoint {ub} must sit at a u domain bound [{u0}, {u1}]"
            );
        }
        assert!(
            (ub_first - ub_last).abs() > 0.5 * (u1 - u0),
            "the two endpoints must sit at OPPOSITE u bounds"
        );
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
