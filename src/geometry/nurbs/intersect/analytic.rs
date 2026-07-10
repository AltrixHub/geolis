//! Closed-form surface×surface intersection for classified pairs.
//!
//! The BIM boolean workload (wall prisms × opening cutter prisms) consists
//! entirely of *extruded* NURBS surfaces — vertical extrusions of degree-≤2
//! wall profiles and planar cutter patches (extrusions of a line). Every
//! such pair intersects in curves with closed forms:
//!
//! - **extrusion × plane, transversal** (`n·d ≠ 0`): the curve is the
//!   explicit graph `t(u) = n·(q − P(u)) / (n·d)` over the profile
//!   parameter, where `S(u, t) = P(u) + t·d`.
//! - **extrusion × plane, parallel** (`n·d = 0`): the curve set is the
//!   iso-lines `u = u*` at the roots of `n·(P(u) − q) = 0`.
//!
//! [`try_intersect`] recognizes these pairs via
//! [`NurbsSurface::extrusion_form`] / [`NurbsSurface::parallelogram_form`]
//! and produces branches with the same conventions the numerical marcher
//! guarantees (see `surface_surface.rs`): endpoints of open branches land
//! exactly on the violated domain bound, sampling density follows
//! `seed_leaf_extent × step_factor`, a closed-in-u seam crossing carries
//! the seam point in both parametric representations, and full-period
//! branches come out `closed`. Unclassified pairs (curved × curved,
//! revolved, free-form) and the ill-conditioned near-parallel band return
//! `None` and fall back to seeding + marching.

use crate::error::Result;
use crate::geometry::nurbs::classify::{ExtrusionForm, ParallelogramForm};
use crate::geometry::nurbs::{NurbsCurve3D, NurbsSurface};
use crate::math::{Point2, Point3, Vector3};

use super::surface_surface::seed_leaf_extent;
use super::types::{IntersectionOptions, SurfaceIntersectionCurve};

/// Transversality gate on `|n̂ · d̂|`: above this the graph form is
/// well-conditioned; below [`PARALLEL_EPS`] the pair is treated as exactly
/// parallel. The band in between (a plane nearly containing the extrusion
/// direction) is analytically valid for BOTH forms but the graph curve gets
/// arbitrarily steep — those pairs fall back to the marcher.
const TRANSVERSAL_MIN: f64 = 1e-2;
/// Below this `|n̂ · d̂|` the plane is parallel to the extrusion direction.
const PARALLEL_EPS: f64 = 1e-9;

/// Minimum samples per emitted branch (loops.rs drops branches with fewer
/// than 3 points; degenerate 2-point slivers would silently vanish).
const MIN_BRANCH_SAMPLES: usize = 4;

/// Bisection iterations for domain-boundary crossings (halves one sweep
/// step to well below `BOUNDARY_EPS`).
const BISECT_ITERS: usize = 60;

/// Attempts the closed-form intersection of `a` × `b`.
///
/// Returns `Ok(None)` when the pair is not classified or is in the
/// near-parallel band (callers run the numerical pipeline instead).
/// `Ok(Some(vec![]))` is a genuine "no intersection" result.
///
/// # Errors
///
/// Propagates profile-curve evaluation failures (rational denominators are
/// validated at surface construction; parameters stay in-domain by
/// construction).
pub(crate) fn try_intersect(
    a: &NurbsSurface,
    b: &NurbsSurface,
    options: &IntersectionOptions,
) -> Result<Option<Vec<SurfaceIntersectionCurve>>> {
    let Some(ext_a) = a.extrusion_form() else {
        return Ok(None);
    };
    let Some(ext_b) = b.extrusion_form() else {
        return Ok(None);
    };

    // Keep the (possibly curved) swept side as `a`-traced; the planar side
    // supplies the plane. When both are planar either choice is exact.
    if let Some(para_b) = b.parallelogram_form() {
        return intersect_extrusion_plane(a, &ext_a, b, &para_b, options);
    }
    if let Some(para_a) = a.parallelogram_form() {
        let Some(mut branches) = intersect_extrusion_plane(b, &ext_b, a, &para_a, options)? else {
            return Ok(None);
        };
        for branch in &mut branches {
            std::mem::swap(&mut branch.uv_a, &mut branch.uv_b);
        }
        return Ok(Some(branches));
    }
    // Curved × curved (e.g. two arc walls) — out of the analytic scope.
    Ok(None)
}

/// The extruded surface `E(u, v) = P(u) + t(v)·d` with its domains and the
/// shared sweep resolution.
struct Sweep<'s> {
    profile: &'s NurbsCurve3D,
    dir: Vector3,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    u_closed: bool,
    /// Sample count over the full u domain.
    samples: usize,
    /// Marcher-equivalent 3D step (drives per-branch densities).
    step: f64,
}

impl Sweep<'_> {
    #[allow(clippy::cast_precision_loss)]
    fn u_at(&self, i: usize) -> f64 {
        if i == self.samples {
            self.u1 // hit the final bound exactly
        } else {
            self.u0 + (self.u1 - self.u0) * (i as f64) / (self.samples as f64)
        }
    }

    fn v_of_t(&self, t: f64) -> f64 {
        self.v0 + t * (self.v1 - self.v0)
    }
}

/// The planar patch with its parameter rectangle.
struct PlanePatch<'p> {
    para: &'p ParallelogramForm,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
}

impl PlanePatch<'_> {
    /// Maps normalized fractions to the patch's actual parameters.
    fn uv_of(&self, fu: f64, fv: f64) -> Point2 {
        Point2::new(
            self.u0 + fu * (self.u1 - self.u0),
            self.v0 + fv * (self.v1 - self.v0),
        )
    }
}

/// Intersects the extruded surface `e_surf` with the planar patch
/// `p_surf`, tracing `uv_a` on the extrusion and `uv_b` on the plane.
/// Returns `None` for the ill-conditioned near-parallel band.
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn intersect_extrusion_plane(
    e_surf: &NurbsSurface,
    ext: &ExtrusionForm,
    p_surf: &NurbsSurface,
    para: &ParallelogramForm,
    options: &IntersectionOptions,
) -> Result<Option<Vec<SurfaceIntersectionCurve>>> {
    let ((eu0, eu1), (ev0, ev1)) = e_surf.parameter_domain();
    let ((pu0, pu1), (pv0, pv1)) = p_surf.parameter_domain();

    let step = (seed_leaf_extent(e_surf, p_surf) * options.step_factor).max(1e-6);
    let profile_len = chord_length(&ext.profile, eu0, eu1)?;
    let samples = ((profile_len / step).ceil() as usize).clamp(32, 100_000);

    let sweep = Sweep {
        profile: &ext.profile,
        dir: ext.direction,
        u0: eu0,
        u1: eu1,
        v0: ev0,
        v1: ev1,
        u_closed: e_surf.is_closed_in_u(),
        samples,
        step,
    };
    let plane = PlanePatch {
        para,
        u0: pu0,
        u1: pu1,
        v0: pv0,
        v1: pv1,
    };

    let n_dot_d = para.normal.dot(&ext.direction.normalize());
    if n_dot_d.abs() >= TRANSVERSAL_MIN {
        graph_branches(&sweep, &plane, options).map(Some)
    } else if n_dot_d.abs() < PARALLEL_EPS {
        parallel_branches(&sweep, &plane, options).map(Some)
    } else {
        Ok(None)
    }
}

/// Estimates the profile's length by a fixed chord sweep.
#[allow(clippy::cast_precision_loss)]
fn chord_length(profile: &NurbsCurve3D, u0: f64, u1: f64) -> Result<f64> {
    const CHORDS: usize = 64;
    let mut len = 0.0;
    let mut prev = profile.point_at(u0)?;
    for i in 1..=CHORDS {
        let u = u0 + (u1 - u0) * (i as f64) / (CHORDS as f64);
        let p = profile.point_at(u)?;
        len += (p - prev).norm();
        prev = p;
    }
    Ok(len)
}

// ---------------------------------------------------------------------------
// Transversal case: explicit graph curve t(u).
// ---------------------------------------------------------------------------

/// One evaluated sample of the transversal graph curve.
#[derive(Clone, Copy)]
struct GraphSample {
    u: f64,
    /// Extrusion fraction at the plane crossing.
    t: f64,
    /// Plane-patch fractions.
    fu: f64,
    fv: f64,
    point: Point3,
}

impl GraphSample {
    fn valid(&self) -> bool {
        (0.0..=1.0).contains(&self.t)
            && (0.0..=1.0).contains(&self.fu)
            && (0.0..=1.0).contains(&self.fv)
    }
}

/// Evaluates the graph curve at `u`.
#[allow(clippy::many_single_char_names)]
fn graph_sample(sweep: &Sweep, plane: &PlanePatch, u: f64) -> Result<GraphSample> {
    let p = sweep.profile.point_at(u)?;
    let n = plane.para.normal;
    let t = n.dot(&(plane.para.origin - p)) / n.dot(&sweep.dir);
    let point = p + sweep.dir * t;
    let (fu, fv) = plane.para.invert(point);
    Ok(GraphSample {
        u,
        t,
        fu,
        fv,
        point,
    })
}

/// Emits the maximal in-domain runs of the graph curve as branches.
fn graph_branches(
    sweep: &Sweep,
    plane: &PlanePatch,
    options: &IntersectionOptions,
) -> Result<Vec<SurfaceIntersectionCurve>> {
    let n = sweep.samples;
    let mut samples = Vec::with_capacity(n + 1);
    for i in 0..=n {
        samples.push(graph_sample(sweep, plane, sweep.u_at(i))?);
    }

    // Full-period closed curve on a closed profile: every sample valid.
    // The first (u0) and last (u1) samples are the same 3D point in both
    // seam representations — the marcher's seam-sample convention.
    if sweep.u_closed && samples.iter().all(GraphSample::valid) {
        return Ok(vec![emit_branch(sweep, plane, &samples, true)]);
    }

    // Maximal valid runs as inclusive index ranges [start, end].
    let mut runs: Vec<(usize, usize)> = Vec::new();
    let mut start: Option<usize> = None;
    for (i, s) in samples.iter().enumerate() {
        match (start, s.valid()) {
            (None, true) => start = Some(i),
            (Some(s0), false) => {
                runs.push((s0, i - 1));
                start = None;
            }
            _ => {}
        }
    }
    if let Some(s0) = start {
        runs.push((s0, n));
    }

    // A closed profile's sweep is circular: a run ending at the last
    // sample continues into the run starting at the first sample (the
    // u-seam; samples 0 and n are one 3D point). Merge them into one
    // branch carrying the seam point in BOTH parametric representations.
    let mut seam_merge: Option<((usize, usize), (usize, usize))> = None;
    if sweep.u_closed && runs.len() >= 2 {
        let first = runs[0];
        let last = runs[runs.len() - 1];
        if first.0 == 0 && last.1 == n {
            seam_merge = Some((last, first));
            runs.pop();
            runs.remove(0);
        }
    }

    let mut branches = Vec::new();
    for &(s0, s1) in &runs {
        if let Some(branch) = graph_run_branch(sweep, plane, &samples, s0, s1, options)? {
            branches.push(branch);
        }
    }
    if let Some(((t0, t1), (h0, h1))) = seam_merge {
        // Tail run [t0, n] then head run [0, h1]: the seam pair (u = u1
        // sample followed by u = u0 sample) sits at the junction.
        let tail = graph_run_branch(sweep, plane, &samples, t0, t1, options)?;
        let head = graph_run_branch(sweep, plane, &samples, h0, h1, options)?;
        match (tail, head) {
            (Some(mut tail), Some(head)) => {
                tail.points.extend(head.points);
                tail.uv_a.extend(head.uv_a);
                tail.uv_b.extend(head.uv_b);
                branches.push(tail);
            }
            (Some(only), None) | (None, Some(only)) => branches.push(only),
            (None, None) => {}
        }
    }
    Ok(branches)
}

/// Emits one valid run `[s0, s1]` as a branch, bisecting the boundary
/// crossings so open endpoints land exactly on the violated domain bound.
#[allow(clippy::cast_precision_loss)]
fn graph_run_branch(
    sweep: &Sweep,
    plane: &PlanePatch,
    samples: &[GraphSample],
    s0: usize,
    s1: usize,
    options: &IntersectionOptions,
) -> Result<Option<SurfaceIntersectionCurve>> {
    let mut run: Vec<GraphSample> = Vec::with_capacity(s1 - s0 + 3);

    // Entry crossing between the invalid neighbor and the first valid
    // sample; exit crossing symmetrically.
    if s0 > 0 {
        run.push(bisect_crossing(
            sweep,
            plane,
            samples[s0].u,
            samples[s0 - 1].u,
        )?);
    }
    run.extend_from_slice(&samples[s0..=s1]);
    if s1 + 1 < samples.len() {
        run.push(bisect_crossing(
            sweep,
            plane,
            samples[s1].u,
            samples[s1 + 1].u,
        )?);
    }

    // Degenerate corner touches: shorter than the coincidence scale.
    let len: f64 = run
        .windows(2)
        .map(|w| (w[1].point - w[0].point).norm())
        .sum();
    if len < options.tolerance * 100.0 {
        return Ok(None);
    }

    // Guarantee the minimum sample count by subdividing the run interval.
    if run.len() < MIN_BRANCH_SAMPLES {
        let (ua, ub) = (run[0].u, run[run.len() - 1].u);
        let mut dense = Vec::with_capacity(MIN_BRANCH_SAMPLES + 1);
        dense.push(run[0]);
        for k in 1..MIN_BRANCH_SAMPLES {
            let u = ua + (ub - ua) * (k as f64) / (MIN_BRANCH_SAMPLES as f64);
            dense.push(graph_sample(sweep, plane, u)?);
        }
        dense.push(run[run.len() - 1]);
        run = dense;
    }

    Ok(Some(emit_branch(sweep, plane, &run, false)))
}

/// Bisects the validity boundary between a valid parameter `u_in` and an
/// invalid `u_out`, snapping the crossing's near-zero constraints exactly
/// onto their bound.
fn bisect_crossing(
    sweep: &Sweep,
    plane: &PlanePatch,
    u_in: f64,
    u_out: f64,
) -> Result<GraphSample> {
    let mut lo = u_in;
    let mut hi = u_out;
    for _ in 0..BISECT_ITERS {
        let mid = 0.5 * (lo + hi);
        if graph_sample(sweep, plane, mid)?.valid() {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let mut s = graph_sample(sweep, plane, lo)?;
    s.t = snap_unit(s.t);
    s.fu = snap_unit(s.fu);
    s.fv = snap_unit(s.fv);
    Ok(s)
}

/// Snaps a normalized fraction onto 0 / 1 when within bisection precision,
/// so endpoint classification (`BOUNDARY_EPS`) sees an exact boundary hit.
fn snap_unit(f: f64) -> f64 {
    const SNAP: f64 = 1e-9;
    if f.abs() < SNAP {
        return 0.0;
    }
    if (f - 1.0).abs() < SNAP {
        return 1.0;
    }
    f
}

/// Converts evaluated samples into a `SurfaceIntersectionCurve`.
fn emit_branch(
    sweep: &Sweep,
    plane: &PlanePatch,
    run: &[GraphSample],
    closed: bool,
) -> SurfaceIntersectionCurve {
    SurfaceIntersectionCurve {
        points: run.iter().map(|s| s.point).collect(),
        uv_a: run
            .iter()
            .map(|s| Point2::new(s.u, sweep.v_of_t(s.t)))
            .collect(),
        uv_b: run.iter().map(|s| plane.uv_of(s.fu, s.fv)).collect(),
        closed,
    }
}

// ---------------------------------------------------------------------------
// Parallel case: iso-lines u = u* at the roots of n·(P(u) − q) = 0.
// ---------------------------------------------------------------------------

/// Appends a root unless an existing one lies within `min_sep` of it (an
/// exactly-zero sample would otherwise duplicate its adjacent bisected
/// crossing).
fn push_root(roots: &mut Vec<f64>, min_sep: f64, u_star: f64) {
    if roots.iter().all(|&r| (r - u_star).abs() > min_sep) {
        roots.push(u_star);
    }
}

/// Emits the iso-line branches of a plane parallel to the extrusion.
#[allow(clippy::many_single_char_names)]
fn parallel_branches(
    sweep: &Sweep,
    plane: &PlanePatch,
    options: &IntersectionOptions,
) -> Result<Vec<SurfaceIntersectionCurve>> {
    let n = plane.para.normal;
    let q = plane.para.origin;
    let f = |u: f64| -> Result<f64> { Ok(n.dot(&(sweep.profile.point_at(u)? - q))) };

    let mut values = Vec::with_capacity(sweep.samples + 1);
    for i in 0..=sweep.samples {
        let u = sweep.u_at(i);
        values.push((u, f(u)?));
    }

    // Coincidence guard: the whole profile lying inside the plane is a
    // tangential/overlap contact — same contract as the marcher: terminate
    // cleanly with no branches.
    let scale = values.iter().map(|(_, v)| v.abs()).fold(0.0, f64::max);
    if scale < options.tolerance * 100.0 {
        return Ok(vec![]);
    }

    // Roots: exactly-zero samples are taken verbatim; STRICT sign changes
    // are bisected. A window with a zero endpoint and no sign change (the
    // curve touching zero at a sample) must NOT be bisected — bisection
    // assumes a crossing and would converge onto the window's far end.
    let min_sep = (sweep.u1 - sweep.u0) * 1e-9;
    let mut roots: Vec<f64> = Vec::new();
    for &(u, fu) in &values {
        if fu == 0.0 {
            push_root(&mut roots, min_sep, u);
        }
    }
    for w in values.windows(2) {
        let ((ua, fa), (ub, fb)) = (w[0], w[1]);
        if fa * fb >= 0.0 {
            continue;
        }
        let mut lo = ua;
        let mut hi = ub;
        let mut flo = fa;
        for _ in 0..BISECT_ITERS {
            let mid = 0.5 * (lo + hi);
            let fm = f(mid)?;
            if (fm > 0.0) == (flo > 0.0) {
                lo = mid;
                flo = fm;
            } else {
                hi = mid;
            }
        }
        push_root(&mut roots, min_sep, 0.5 * (lo + hi));
    }
    roots.sort_by(f64::total_cmp);

    let mut branches = Vec::new();
    for u_star in roots {
        if let Some(branch) = iso_line_branch(sweep, plane, u_star, options)? {
            branches.push(branch);
        }
    }
    Ok(branches)
}

/// Builds the `u = u*` iso-line branch, clipped to the plane patch.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    clippy::many_single_char_names
)]
fn iso_line_branch(
    sweep: &Sweep,
    plane: &PlanePatch,
    u_star: f64,
    options: &IntersectionOptions,
) -> Result<Option<SurfaceIntersectionCurve>> {
    let p = sweep.profile.point_at(u_star)?;

    // fu(t) and fv(t) are affine in the extrusion fraction t: intersect
    // t ∈ [0, 1] with both patch constraints analytically.
    let (fu0, fv0) = plane.para.invert(p);
    let (fu1, fv1) = plane.para.invert(p + sweep.dir);
    let mut t_lo = 0.0_f64;
    let mut t_hi = 1.0_f64;
    for (c0, c1) in [(fu0, fu1), (fv0, fv1)] {
        let d = c1 - c0;
        if d.abs() < 1e-15 {
            if !(0.0..=1.0).contains(&c0) {
                return Ok(None);
            }
            continue;
        }
        let (t_a, t_b) = ((0.0 - c0) / d, (1.0 - c0) / d);
        let (lo, hi) = if t_a <= t_b { (t_a, t_b) } else { (t_b, t_a) };
        t_lo = t_lo.max(lo);
        t_hi = t_hi.min(hi);
    }
    if t_hi <= t_lo {
        return Ok(None);
    }

    let seg_len = sweep.dir.norm() * (t_hi - t_lo);
    if seg_len < options.tolerance * 100.0 {
        return Ok(None);
    }

    let count = ((seg_len / sweep.step).ceil() as usize).clamp(MIN_BRANCH_SAMPLES, 100_000);
    let mut points = Vec::with_capacity(count + 1);
    let mut uv_a = Vec::with_capacity(count + 1);
    let mut uv_b = Vec::with_capacity(count + 1);
    for k in 0..=count {
        let mut t = t_lo + (t_hi - t_lo) * (k as f64) / (count as f64);
        // Exact segment ends: snap onto the extrusion's own v bounds when
        // they coincide (endpoint classification needs exact hits).
        if t.abs() < 1e-12 {
            t = 0.0;
        }
        if (t - 1.0).abs() < 1e-12 {
            t = 1.0;
        }
        let x = p + sweep.dir * t;
        let (fu, fv) = plane.para.invert(x);
        points.push(x);
        uv_a.push(Point2::new(u_star, sweep.v_of_t(t)));
        uv_b.push(plane.uv_of(snap_unit(fu), snap_unit(fv)));
    }
    Ok(Some(SurfaceIntersectionCurve {
        points,
        uv_a,
        uv_b,
        closed: false,
    }))
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use crate::geometry::nurbs::intersect::surface_surface::intersect_surfaces_marching;
    use crate::geometry::nurbs::KnotVector;

    /// Cutter-style vertical parallelogram patch in the plane `x = x0`.
    fn x_patch(x0: f64, y_lo: f64, y_hi: f64, z_lo: f64, z_hi: f64) -> NurbsSurface {
        let line = NurbsCurve3D::from_unweighted(
            vec![Point3::new(x0, y_lo, z_lo), Point3::new(x0, y_hi, z_lo)],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
        )
        .unwrap();
        NurbsSurface::extrude(&line, Vector3::new(0.0, 0.0, z_hi - z_lo)).unwrap()
    }

    /// Horizontal parallelogram patch at `z = z0`.
    fn z_patch(z0: f64, x_lo: f64, x_hi: f64, y_lo: f64, y_hi: f64) -> NurbsSurface {
        let line = NurbsCurve3D::from_unweighted(
            vec![Point3::new(x_lo, y_lo, z0), Point3::new(x_hi, y_lo, z0)],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
        )
        .unwrap();
        NurbsSurface::extrude(&line, Vector3::new(0.0, y_hi - y_lo, 0.0)).unwrap()
    }

    /// Wall-style surface: a half circle (radius 2, XY plane) extruded up.
    fn arc_wall(height: f64) -> NurbsSurface {
        let arc = NurbsCurve3D::arc(
            Point3::origin(),
            2.0,
            Vector3::z(),
            Vector3::x(),
            0.0,
            std::f64::consts::PI,
        )
        .unwrap();
        NurbsSurface::extrude(&arc, Vector3::new(0.0, 0.0, height)).unwrap()
    }

    /// Full closed cylinder (annulus-wall style), radius 1.
    fn full_cylinder(height: f64) -> NurbsSurface {
        let circle =
            NurbsCurve3D::circle(Point3::origin(), 1.0, Vector3::z(), Vector3::x()).unwrap();
        NurbsSurface::extrude(&circle, Vector3::new(0.0, 0.0, height)).unwrap()
    }

    /// Every branch point must satisfy BOTH surface equations through its
    /// UV traces — the strongest per-point correctness property.
    fn assert_on_both(a: &NurbsSurface, b: &NurbsSurface, branches: &[SurfaceIntersectionCurve]) {
        for br in branches {
            assert_eq!(br.points.len(), br.uv_a.len());
            assert_eq!(br.points.len(), br.uv_b.len());
            for ((p, ua), ub) in br.points.iter().zip(&br.uv_a).zip(&br.uv_b) {
                let pa = a.point_at(ua.x, ua.y).unwrap();
                let pb = b.point_at(ub.x, ub.y).unwrap();
                assert!(
                    (pa - p).norm() < 1e-9,
                    "uv_a off the curve: {pa:?} vs {p:?}"
                );
                assert!(
                    (pb - p).norm() < 1e-9,
                    "uv_b off the curve: {pb:?} vs {p:?}"
                );
            }
        }
    }

    /// The analytic branch set must match the marcher's on branch count,
    /// closure flags, and endpoint geometry (within one marching step).
    fn assert_matches_marcher(a: &NurbsSurface, b: &NurbsSurface) {
        let options = IntersectionOptions::default();
        let analytic = try_intersect(a, b, &options)
            .unwrap()
            .expect("pair must classify for the analytic path");
        let marched = intersect_surfaces_marching(a, b, &options).unwrap();
        assert_on_both(a, b, &analytic);

        assert_eq!(
            analytic.len(),
            marched.len(),
            "branch count: analytic {} vs marcher {}",
            analytic.len(),
            marched.len()
        );
        let step = (seed_leaf_extent(a, b) * options.step_factor).max(1e-6);
        for branch in &analytic {
            let has_counterpart = marched.iter().any(|counterpart| {
                if counterpart.closed != branch.closed {
                    return false;
                }
                if branch.closed {
                    // Closed loops have no distinguished endpoints; match by
                    // every marcher point lying within a step of the
                    // analytic polyline sampling.
                    return counterpart
                        .points
                        .iter()
                        .all(|mp| branch.points.iter().any(|ap| (ap - mp).norm() < 2.0 * step));
                }
                let (a0, a1) = (branch.points[0], branch.points[branch.points.len() - 1]);
                let (m0, m1) = (
                    counterpart.points[0],
                    counterpart.points[counterpart.points.len() - 1],
                );
                let fwd = (a0 - m0).norm() < 2.0 * step && (a1 - m1).norm() < 2.0 * step;
                let rev = (a0 - m1).norm() < 2.0 * step && (a1 - m0).norm() < 2.0 * step;
                fwd || rev
            });
            assert!(
                has_counterpart,
                "analytic branch (closed={}, {:?}..{:?}) has no marcher counterpart",
                branch.closed,
                branch.points.first(),
                branch.points.last()
            );
        }
    }

    #[test]
    fn straight_wall_times_vertical_jamb_matches_marcher() {
        // Parallel case: both planes vertical, iso-line branch.
        let wall = x_patch(0.0, -3.0, 3.0, 0.0, 2.4); // wall face in x=0 plane
        let jamb = z_patch(1.2, -1.0, 1.0, -0.5, 0.5); // horizontal plane crossing it
        assert_matches_marcher(&wall, &jamb);
    }

    #[test]
    fn arc_wall_times_sill_plane_matches_marcher() {
        // Transversal graph case: horizontal sill plane × curved wall.
        let wall = arc_wall(2.4);
        let sill = z_patch(0.9, -1.5, 1.5, 0.5, 2.5);
        assert_matches_marcher(&wall, &sill);
    }

    #[test]
    fn arc_wall_times_vertical_jamb_matches_marcher() {
        // Parallel case on a curved wall: vertical jamb plane cuts the arc
        // at two azimuths -> two iso-line branches.
        let wall = arc_wall(2.4);
        let jamb = x_patch(0.7, -3.0, 3.0, 0.4, 1.9);
        assert_matches_marcher(&wall, &jamb);
    }

    #[test]
    fn closed_cylinder_full_period_sill_is_closed_loop() {
        // A sill plane crossing the ENTIRE closed cylinder: one closed loop
        // spanning the full u period.
        let wall = full_cylinder(2.0);
        let sill = z_patch(1.0, -2.0, 2.0, -2.0, 2.0);
        let options = IntersectionOptions::default();
        let analytic = try_intersect(&wall, &sill, &options)
            .unwrap()
            .expect("pair must classify");
        assert_on_both(&wall, &sill, &analytic);
        assert_eq!(analytic.len(), 1, "one full-period loop");
        assert!(analytic[0].closed, "full-period branch must be closed");
        // Seam representation: first sample at u_min, last at u_max, same
        // 3D point (the marcher's seam-sample convention).
        let ((u0, u1), _) = wall.parameter_domain();
        let br = &analytic[0];
        assert!((br.uv_a[0].x - u0).abs() < 1e-12);
        assert!((br.uv_a[br.uv_a.len() - 1].x - u1).abs() < 1e-12);
        assert!(
            (br.points[0] - br.points[br.points.len() - 1]).norm() < 1e-9,
            "seam endpoints must be one 3D point"
        );
    }

    #[test]
    fn closed_cylinder_seam_straddling_patch_merges_runs() {
        // A sill patch clipped to straddle the +X seam azimuth: the two
        // sweep runs merge into ONE branch across the seam.
        let wall = full_cylinder(2.0);
        let sill = z_patch(1.0, 0.5, 2.0, -0.8, 0.8); // x in [0.5,2] straddles azimuth 0
        let options = IntersectionOptions::default();
        let analytic = try_intersect(&wall, &sill, &options)
            .unwrap()
            .expect("pair must classify");
        assert_on_both(&wall, &sill, &analytic);
        assert_eq!(
            analytic.len(),
            1,
            "seam-straddling arc must be one merged branch"
        );
        assert!(!analytic[0].closed);
        // The seam junction carries the same 3D point in both parametric
        // representations (u = u1 sample adjacent to u = u0 sample).
        let br = &analytic[0];
        let ((u0, u1), _) = wall.parameter_domain();
        let seam_idx = br
            .uv_a
            .windows(2)
            .position(|w| (w[0].x - u1).abs() < 1e-9 && (w[1].x - u0).abs() < 1e-9)
            .expect("seam sample pair present");
        assert!(
            (br.points[seam_idx] - br.points[seam_idx + 1]).norm() < 1e-9,
            "seam pair must be one 3D point"
        );
    }

    #[test]
    fn disjoint_pair_yields_empty() {
        let wall = arc_wall(2.4);
        let far = z_patch(0.9, 10.0, 12.0, 10.0, 12.0);
        let options = IntersectionOptions::default();
        let analytic = try_intersect(&wall, &far, &options)
            .unwrap()
            .expect("pair must classify");
        assert!(analytic.is_empty());
    }

    #[test]
    fn coincident_planes_terminate_cleanly_with_no_branches() {
        let a = z_patch(1.0, -1.0, 1.0, -1.0, 1.0);
        let b = z_patch(1.0, -0.5, 0.5, -0.5, 0.5);
        let options = IntersectionOptions::default();
        let analytic = try_intersect(&a, &b, &options)
            .unwrap()
            .expect("pair must classify");
        assert!(
            analytic.is_empty(),
            "coincident overlap: clean, no branches"
        );
    }
}
