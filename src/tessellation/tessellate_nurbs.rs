//! Adaptive tessellation of NURBS curves and surfaces.
//!
//! Curves are sampled by recursive chord-deviation bisection, seeded at
//! distinct interior knots so degree-reducing discontinuities are honored.
//! Surfaces are meshed onto a uniformly refined grid whose division count is
//! doubled per direction until adjacent-sample normal deviation falls under
//! tolerance, also seeded with interior knot lines.

use crate::error::{Result, TessellationError};
use crate::geometry::nurbs::{NurbsCurve3D, NurbsSurface};
use crate::geometry::surface::Surface;
use crate::math::{Point2, Point3, Vector3};

use super::TriangleMesh;

/// Options for adaptive NURBS curve tessellation.
#[derive(Debug, Clone, Copy)]
pub struct CurveTessellationOptions {
    /// Maximum chord deviation from the true curve, in model units.
    pub chord_tolerance: f64,
    /// Maximum recursion depth per seed interval.
    pub max_depth: usize,
}

impl Default for CurveTessellationOptions {
    fn default() -> Self {
        Self {
            chord_tolerance: 1e-3,
            max_depth: 16,
        }
    }
}

/// Adaptively samples a NURBS curve into a polyline.
///
/// The parameter domain is first split at every distinct interior knot, then
/// each sub-interval is recursively bisected while the curve's midpoint
/// deviates from the chord midpoint by more than `chord_tolerance` (capped at
/// `max_depth`). Endpoints are always emitted exactly once; the returned points
/// are ordered and contain no consecutive duplicates.
///
/// # Errors
///
/// Returns [`TessellationError::InvalidParameters`] if `chord_tolerance` is not
/// strictly positive, or propagates any evaluation error from the curve.
pub fn tessellate_nurbs_curve(
    curve: &NurbsCurve3D,
    options: &CurveTessellationOptions,
) -> Result<Vec<Point3>> {
    let params = tessellate_nurbs_curve_params(curve, options)?;
    params.iter().map(|&t| curve.point_at(t)).collect()
}

/// Adaptively samples a NURBS curve, returning the ordered **parameter** values
/// (endpoints included) rather than the evaluated points.
///
/// This is the geometry-only sampling used to make adjacent faces conform: two
/// faces that share a boundary curve sample it at the identical parameter set,
/// so their boundary polylines coincide exactly (see [`conforming_boundary_uv`]).
///
/// # Errors
///
/// Returns [`TessellationError::InvalidParameters`] if `chord_tolerance` is not
/// strictly positive, or propagates any evaluation error from the curve.
pub(crate) fn tessellate_nurbs_curve_params(
    curve: &NurbsCurve3D,
    options: &CurveTessellationOptions,
) -> Result<Vec<f64>> {
    if options.chord_tolerance <= 0.0 {
        return Err(TessellationError::InvalidParameters(
            "chord_tolerance must be strictly positive".to_owned(),
        )
        .into());
    }

    let (t_min, t_max) = curve.parameter_domain();
    let seeds = seed_parameters(curve, t_min, t_max);

    // Start with the first seed parameter, then append each interval's interior +
    // end parameters so endpoints are shared exactly once.
    let mut params = Vec::new();
    let mut prev_end = Sample {
        t: seeds[0],
        point: curve.point_at(seeds[0])?,
    };
    params.push(prev_end.t);
    for window in seeds.windows(2) {
        let end = Sample {
            t: window[1],
            point: curve.point_at(window[1])?,
        };
        bisect_curve(curve, prev_end, end, options, 0, &mut params)?;
        params.push(end.t);
        prev_end = end;
    }

    Ok(params)
}

/// A curve sample: parameter plus its evaluated point.
#[derive(Clone, Copy)]
struct Sample {
    t: f64,
    point: Point3,
}

/// Distinct parameter seeds: domain endpoints plus every interior knot that
/// lies strictly inside the domain. Returned sorted and de-duplicated, always
/// containing at least the two endpoints.
fn seed_parameters(curve: &NurbsCurve3D, t_min: f64, t_max: f64) -> Vec<f64> {
    let mut seeds = vec![t_min, t_max];
    let span = (t_max - t_min).abs();
    let eps = span * 1e-9;
    for &k in curve.knots().as_slice() {
        if k > t_min + eps && k < t_max - eps {
            seeds.push(k);
        }
    }
    seeds.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    seeds.dedup_by(|a, b| (*a - *b).abs() <= eps);
    seeds
}

/// Recursively appends interior sample **parameters** for the open interval
/// `(start.t, end.t)`.
///
/// `start` / `end` are the already-evaluated endpoints. The midpoint parameter
/// is appended only when the curve deviates from the chord beyond tolerance and
/// depth remains; the `end` endpoint is appended by the caller.
fn bisect_curve(
    curve: &NurbsCurve3D,
    start: Sample,
    end: Sample,
    options: &CurveTessellationOptions,
    depth: usize,
    out: &mut Vec<f64>,
) -> Result<()> {
    let mid_t = 0.5 * (start.t + end.t);
    let mid = Sample {
        t: mid_t,
        point: curve.point_at(mid_t)?,
    };

    let chord_mid = Point3::from(0.5 * (start.point.coords + end.point.coords));
    let deviation = (mid.point - chord_mid).norm();

    if depth >= options.max_depth || deviation <= options.chord_tolerance {
        return Ok(());
    }

    bisect_curve(curve, start, mid, options, depth + 1, out)?;
    out.push(mid.t);
    bisect_curve(curve, mid, end, options, depth + 1, out)?;
    Ok(())
}

/// Chord tolerance (model units) for boundary-curve sampling shared by adjacent
/// faces. Because both faces sample the *same* boundary curve at the *same*
/// parameters, the inter-face deviation is exact regardless of this value; it
/// only governs how finely the silhouette itself is approximated.
pub(crate) const BOUNDARY_CHORD_TOLERANCE: f64 = 1e-3;

/// Builds the curve-intrinsic UV boundary polyline for a NURBS face whose outer
/// boundary is the full parameter rectangle.
///
/// Each of the four boundary isocurves is sampled at its own chord-adaptive
/// parameters ([`tessellate_nurbs_curve_params`]) — a function of the curve
/// geometry alone. Any other face that shares one of these boundary curves (for
/// example a ruled side wall extruded from it) samples the identical curve at
/// the identical parameters, so the two faces emit coincident 3D boundary
/// vertices and no silhouette sliver forms. The returned loop is CCW in UV:
/// `v_min` edge (u increasing), `u_max` edge (v increasing), `v_max` edge (u
/// decreasing), `u_min` edge (v decreasing).
///
/// # Errors
///
/// Propagates isocurve extraction or curve-sampling errors.
pub(crate) fn conforming_boundary_uv(
    surface: &NurbsSurface,
    chord_tolerance: f64,
) -> Result<Vec<Point2>> {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_domain();
    let options = CurveTessellationOptions {
        chord_tolerance,
        max_depth: 16,
    };

    // Boundary isocurves with their parameter axis:
    //   v_min / v_max edges vary u (knots_u domain [u_min, u_max]).
    //   u_min / u_max edges vary v (knots_v domain [v_min, v_max]).
    let u_at_vmin = tessellate_nurbs_curve_params(&surface.isocurve_v(v_min)?, &options)?;
    let v_at_umax = tessellate_nurbs_curve_params(&surface.isocurve_u(u_max)?, &options)?;
    let u_at_vmax = tessellate_nurbs_curve_params(&surface.isocurve_v(v_max)?, &options)?;
    let v_at_umin = tessellate_nurbs_curve_params(&surface.isocurve_u(u_min)?, &options)?;

    let mut pts: Vec<Point2> = Vec::new();
    for &u in &u_at_vmin {
        pts.push(Point2::new(u, v_min));
    }
    for &v in &v_at_umax {
        pts.push(Point2::new(u_max, v));
    }
    for &u in u_at_vmax.iter().rev() {
        pts.push(Point2::new(u, v_max));
    }
    for &v in v_at_umin.iter().rev() {
        pts.push(Point2::new(u_min, v));
    }

    // Drop consecutive near-duplicate points (corner joins) and the closing
    // wrap-around duplicate so the CDT never receives a zero-length constraint.
    pts.dedup_by(|a, b| (*a - *b).norm() < 1e-9);
    while pts.len() >= 2 && (pts[0] - pts[pts.len() - 1]).norm() < 1e-9 {
        pts.pop();
    }
    Ok(pts)
}

/// Reports whether an untrimmed NURBS surface's parameter rectangle maps to a
/// non-degenerate quadrilateral boundary (four distinct corners).
///
/// Closed/seam surfaces (a revolved wall, an extruded tube whose `u` wraps) or
/// pole-collapsed patches fail this test: their opposite boundary edges
/// coincide, so the rectangle-outline CDT path is inapplicable and the caller
/// falls back to the tensor-grid tessellator.
pub(crate) fn nurbs_surface_is_open(surface: &NurbsSurface) -> bool {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_domain();
    let corners = [
        surface.point_at(u_min, v_min),
        surface.point_at(u_max, v_min),
        surface.point_at(u_max, v_max),
        surface.point_at(u_min, v_max),
    ];
    let corners: Vec<Point3> = match corners.into_iter().collect::<Result<_>>() {
        Ok(c) => c,
        Err(_) => return false,
    };
    for i in 0..corners.len() {
        for j in (i + 1)..corners.len() {
            if (corners[i] - corners[j]).norm() < crate::math::TOLERANCE {
                return false;
            }
        }
    }
    true
}

/// Options for adaptive NURBS surface tessellation.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceTessellationOptions {
    /// Maximum normal deviation (radians) tolerated between adjacent grid
    /// samples before the offending direction is refined.
    pub normal_tolerance: f64,
    /// Minimum subdivision count per direction.
    pub min_divisions: usize,
    /// Maximum subdivision count per direction.
    pub max_divisions: usize,
}

impl Default for SurfaceTessellationOptions {
    fn default() -> Self {
        Self {
            normal_tolerance: 0.05,
            min_divisions: 4,
            max_divisions: 128,
        }
    }
}

/// Adaptively tessellates a NURBS surface into a watertight triangle mesh.
///
/// The parameter rectangle is sampled on a uniform grid that is independently
/// refined in `u` and `v`: each direction starts at `min_divisions` intervals
/// (merged with that direction's interior knot lines) and its interval count is
/// doubled while any adjacent-sample normal pair deviates by more than
/// `normal_tolerance`, capped at `max_divisions`. Because all triangles share
/// grid-vertex indices, every interior edge is referenced by exactly two
/// triangles.
///
/// Emitted vertices carry positions, unit normals, and `[u, v]` UVs normalized
/// to `[0, 1]^2`.
///
/// # Errors
///
/// Returns [`TessellationError::InvalidParameters`] if `normal_tolerance` is
/// not strictly positive or `min_divisions` is zero or exceeds
/// `max_divisions`, and propagates any surface evaluation error.
pub fn tessellate_nurbs_surface(
    surface: &NurbsSurface,
    options: &SurfaceTessellationOptions,
) -> Result<TriangleMesh> {
    validate_surface_options(options)?;

    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_domain();
    let (u_params, v_params) = adaptive_grid_parameters(surface, options);
    build_grid_mesh(surface, &u_params, &v_params, u_min, u_max, v_min, v_max)
}

/// Computes the adaptive `u`/`v` sampling-parameter lists for `surface` under
/// `options`: each axis starts at `min_divisions` uniform intervals (merged with
/// that axis's interior knot lines) and is doubled while the adjacent-sample
/// normal deviation exceeds `normal_tolerance`, capped at `max_divisions`.
///
/// Shared by the full-domain surface tessellator and the trimmed tessellator so
/// both use the identical refinement; the caller is responsible for validating
/// `options`.
pub(crate) fn adaptive_grid_parameters(
    surface: &NurbsSurface,
    options: &SurfaceTessellationOptions,
) -> (Vec<f64>, Vec<f64>) {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_domain();
    let interior_u = interior_knots(surface.knots_u().as_slice(), u_min, u_max);
    let interior_v = interior_knots(surface.knots_v().as_slice(), v_min, v_max);
    let u_params = refine_axis(surface, options, u_min, u_max, &interior_u, Axis::U);
    let v_params = refine_axis(surface, options, v_min, v_max, &interior_v, Axis::V);
    (u_params, v_params)
}

/// Direction along which an adjacent-sample normal check is taken.
#[derive(Clone, Copy)]
enum Axis {
    /// Vary `u` while holding `v` fixed at several stations.
    U,
    /// Vary `v` while holding `u` fixed at several stations.
    V,
}

fn validate_surface_options(options: &SurfaceTessellationOptions) -> Result<()> {
    if options.normal_tolerance <= 0.0 {
        return Err(TessellationError::InvalidParameters(
            "normal_tolerance must be strictly positive".to_owned(),
        )
        .into());
    }
    if options.min_divisions == 0 || options.min_divisions > options.max_divisions {
        return Err(TessellationError::InvalidParameters(
            "require 1 <= min_divisions <= max_divisions".to_owned(),
        )
        .into());
    }
    Ok(())
}

/// Distinct interior knot values strictly inside `(lo, hi)`.
fn interior_knots(knots: &[f64], lo: f64, hi: f64) -> Vec<f64> {
    let eps = (hi - lo).abs() * 1e-9;
    let mut out: Vec<f64> = knots
        .iter()
        .copied()
        .filter(|&k| k > lo + eps && k < hi - eps)
        .collect();
    out.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    out.dedup_by(|a, b| (*a - *b).abs() <= eps);
    out
}

/// Build a sorted, de-duplicated parameter list of `divisions` uniform steps
/// over `[lo, hi]` merged with the supplied interior knot lines.
fn axis_parameters(lo: f64, hi: f64, divisions: usize, interior: &[f64]) -> Vec<f64> {
    let eps = (hi - lo).abs() * 1e-9;
    let mut params = Vec::with_capacity(divisions + 1 + interior.len());
    for i in 0..=divisions {
        #[allow(clippy::cast_precision_loss)]
        let frac = i as f64 / divisions as f64;
        params.push(lo + frac * (hi - lo));
    }
    params.extend_from_slice(interior);
    params.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    params.dedup_by(|a, b| (*a - *b).abs() <= eps);
    params
}

/// Maximum normal deviation (radians) between adjacent samples taken along
/// `axis` over the given parameter stations, holding the other parameter fixed
/// at a few representative cross-stations.
fn max_normal_deviation_along(
    surface: &NurbsSurface,
    params: &[f64],
    cross_lo: f64,
    cross_hi: f64,
    axis: Axis,
) -> f64 {
    // A handful of fixed cross-stations is enough to catch curvature; the
    // interior point avoids degenerate normals at rational seams.
    let cross_stations = [cross_lo, 0.5 * (cross_lo + cross_hi), cross_hi];

    let mut max_dev = 0.0_f64;
    for &cross in &cross_stations {
        let mut prev: Option<Vector3> = None;
        for &p in params {
            let normal = match axis {
                Axis::U => surface.normal(p, cross),
                Axis::V => surface.normal(cross, p),
            };
            // Skip degenerate normals (e.g. a collapsed pole) rather than fail
            // the whole tessellation.
            let Ok(n) = normal else {
                prev = None;
                continue;
            };
            if let Some(pn) = prev {
                let dev = angle_between(&pn, &n);
                if dev > max_dev {
                    max_dev = dev;
                }
            }
            prev = Some(n);
        }
    }
    max_dev
}

/// Angle in radians between two unit vectors, clamped for numerical safety.
fn angle_between(a: &Vector3, b: &Vector3) -> f64 {
    a.dot(b).clamp(-1.0, 1.0).acos()
}

/// Determine the sampling parameters for one axis by doubling the uniform
/// division count until the adjacent-sample normal deviation falls under
/// tolerance or the cap is hit.
fn refine_axis(
    surface: &NurbsSurface,
    options: &SurfaceTessellationOptions,
    lo: f64,
    hi: f64,
    interior: &[f64],
    axis: Axis,
) -> Vec<f64> {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_domain();
    let (cross_lo, cross_hi) = match axis {
        Axis::U => (v_min, v_max),
        Axis::V => (u_min, u_max),
    };

    let mut divisions = options.min_divisions;
    loop {
        let params = axis_parameters(lo, hi, divisions, interior);
        let dev = max_normal_deviation_along(surface, &params, cross_lo, cross_hi, axis);
        if dev <= options.normal_tolerance || divisions >= options.max_divisions {
            return params;
        }
        divisions = (divisions * 2).min(options.max_divisions);
    }
}

/// Build a watertight triangle mesh from a tensor grid of `u`/`v` parameters.
#[allow(clippy::cast_possible_truncation, clippy::too_many_arguments)]
fn build_grid_mesh(
    surface: &NurbsSurface,
    u_params: &[f64],
    v_params: &[f64],
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
) -> Result<TriangleMesh> {
    let nu = u_params.len();
    let nv = v_params.len();

    let u_span = u_max - u_min;
    let v_span = v_max - v_min;

    let mut vertices = Vec::with_capacity(nu * nv);
    let mut normals = Vec::with_capacity(nu * nv);
    let mut uvs = Vec::with_capacity(nu * nv);

    for &u in u_params {
        for &v in v_params {
            let p = surface.point_at(u, v)?;
            // A collapsed pole yields a zero normal; fall back to +Z so the
            // mesh stays well-formed (callers can recompute if needed).
            let n = surface.normal(u, v).unwrap_or_else(|_| Vector3::z());
            vertices.push(p);
            normals.push(n);
            let su = if u_span.abs() > f64::EPSILON {
                (u - u_min) / u_span
            } else {
                0.0
            };
            let sv = if v_span.abs() > f64::EPSILON {
                (v - v_min) / v_span
            } else {
                0.0
            };
            uvs.push(Point2::new(su, sv));
        }
    }

    let mut indices = Vec::with_capacity((nu.saturating_sub(1)) * (nv.saturating_sub(1)) * 2);
    let idx = |i: usize, j: usize| -> u32 { (i * nv + j) as u32 };
    for i in 0..nu.saturating_sub(1) {
        for j in 0..nv.saturating_sub(1) {
            let v00 = idx(i, j);
            let v10 = idx(i + 1, j);
            let v01 = idx(i, j + 1);
            let v11 = idx(i + 1, j + 1);
            // Consistent winding: (v00, v10, v11) then (v00, v11, v01).
            indices.push([v00, v10, v11]);
            indices.push([v00, v11, v01]);
        }
    }

    Ok(TriangleMesh {
        vertices,
        normals,
        uvs,
        indices,
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod curve_tests {
    use super::*;
    use crate::geometry::nurbs::{KnotVector, NurbsCurve3D};
    use crate::math::Vector3;

    /// Straight degree-1 curve from (0,0,0) to (5,0,0).
    fn line_curve() -> NurbsCurve3D {
        NurbsCurve3D::from_unweighted(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(5.0, 0.0, 0.0)],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
        )
        .unwrap()
    }

    fn unit_circle() -> NurbsCurve3D {
        NurbsCurve3D::circle(
            Point3::new(0.0, 0.0, 0.0),
            1.0,
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 0.0),
        )
        .unwrap()
    }

    #[test]
    fn straight_line_tessellates_to_two_endpoints() {
        let curve = line_curve();
        let pts = tessellate_nurbs_curve(&curve, &CurveTessellationOptions::default()).unwrap();
        assert_eq!(pts.len(), 2, "degree-1 line must not oversample");
        assert!((pts[0] - Point3::new(0.0, 0.0, 0.0)).norm() < 1e-12);
        assert!((pts[1] - Point3::new(5.0, 0.0, 0.0)).norm() < 1e-12);
    }

    #[test]
    fn rational_circle_within_chord_tolerance() {
        let curve = unit_circle();
        let options = CurveTessellationOptions {
            chord_tolerance: 1e-3,
            max_depth: 16,
        };
        let pts = tessellate_nurbs_curve(&curve, &options).unwrap();

        // Every emitted point lies on the unit circle.
        for p in &pts {
            let r = (p.x * p.x + p.y * p.y).sqrt();
            assert!((r - 1.0).abs() < 1e-9, "point off circle: r = {r}");
            assert!(p.z.abs() < 1e-9, "point off plane: z = {}", p.z);
        }

        // Adjacent-chord deviation: sample the true curve midpoint of each
        // emitted segment and measure its distance to the chord midpoint.
        for w in pts.windows(2) {
            let chord_mid = Point3::from(0.5 * (w[0].coords + w[1].coords));
            // Invert each endpoint's angle to find the bracketing parameters,
            // then evaluate the true midpoint by angle.
            let a0 = w[0].y.atan2(w[0].x);
            let a1 = w[1].y.atan2(w[1].x);
            let mut da = a1 - a0;
            // Normalize to the short arc.
            while da > std::f64::consts::PI {
                da -= std::f64::consts::TAU;
            }
            while da < -std::f64::consts::PI {
                da += std::f64::consts::TAU;
            }
            let mid_angle = a0 + 0.5 * da;
            let true_mid = Point3::new(mid_angle.cos(), mid_angle.sin(), 0.0);
            let deviation = (true_mid - chord_mid).norm();
            assert!(
                deviation < options.chord_tolerance + 1e-9,
                "segment deviation {deviation} exceeds tolerance"
            );
        }

        assert!(pts.len() < 200, "unreasonable point count: {}", pts.len());
        assert!(pts.len() > 4, "circle should be refined: {}", pts.len());
    }

    #[test]
    fn tightening_tolerance_increases_point_count() {
        let curve = unit_circle();
        let coarse = tessellate_nurbs_curve(
            &curve,
            &CurveTessellationOptions {
                chord_tolerance: 1e-2,
                max_depth: 20,
            },
        )
        .unwrap();
        let fine = tessellate_nurbs_curve(
            &curve,
            &CurveTessellationOptions {
                chord_tolerance: 1e-3,
                max_depth: 20,
            },
        )
        .unwrap();
        assert!(
            fine.len() > coarse.len(),
            "tighter tolerance must add points: {} vs {}",
            fine.len(),
            coarse.len()
        );
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod surface_tests {
    use super::*;
    use crate::geometry::nurbs::KnotVector;
    use std::collections::HashMap;

    /// 2x2 bilinear (planar) patch spanning [0,2]x[0,2] in the XY plane.
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

    /// Quarter-cylinder shell: rational quadratic quarter circle in u (radius 1
    /// about the origin in XY) extruded linearly in v along +Z by 2.
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

    /// Count how many triangles reference each undirected edge.
    fn edge_use_counts(mesh: &TriangleMesh) -> HashMap<(u32, u32), usize> {
        let mut counts: HashMap<(u32, u32), usize> = HashMap::new();
        for tri in &mesh.indices {
            for &(a, b) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let key = if a < b { (a, b) } else { (b, a) };
                *counts.entry(key).or_insert(0) += 1;
            }
        }
        counts
    }

    #[test]
    fn planar_patch_stays_at_min_divisions_and_is_watertight() {
        let surface = bilinear_patch();
        let options = SurfaceTessellationOptions::default();
        let mesh = tessellate_nurbs_surface(&surface, &options).unwrap();

        // A flat patch never triggers refinement: a (min+1)^2 vertex grid.
        let expected_dim = options.min_divisions + 1;
        assert_eq!(
            mesh.vertices.len(),
            expected_dim * expected_dim,
            "flat patch must stay at min_divisions"
        );

        // All normals identical (planar surface) and unit length.
        let n0 = mesh.normals[0];
        for n in &mesh.normals {
            assert!((n - n0).norm() < 1e-9, "planar normals must be identical");
            assert!((n.norm() - 1.0).abs() < 1e-9, "normal must be unit length");
        }

        // Watertight: every interior edge shared by exactly 2 triangles, every
        // boundary edge by exactly 1.
        let counts = edge_use_counts(&mesh);
        let interior = counts.values().filter(|&&c| c == 2).count();
        let boundary = counts.values().filter(|&&c| c == 1).count();
        assert_eq!(
            interior + boundary,
            counts.len(),
            "every edge used by 1 or 2 triangles"
        );
        assert!(interior > 0 && boundary > 0);
    }

    #[test]
    fn quarter_cylinder_refines_u_only_and_lies_on_cylinder() {
        let surface = quarter_cylinder_patch();
        let options = SurfaceTessellationOptions::default();
        let mesh = tessellate_nurbs_surface(&surface, &options).unwrap();

        // Refinement happens in u (curved) but not v (straight extrusion).
        let interior_u = interior_knots(surface.knots_u().as_slice(), 0.0, 1.0);
        let interior_v = interior_knots(surface.knots_v().as_slice(), 0.0, 1.0);
        let u_params = refine_axis(&surface, &options, 0.0, 1.0, &interior_u, Axis::U);
        let v_params = refine_axis(&surface, &options, 0.0, 1.0, &interior_v, Axis::V);
        assert!(
            u_params.len() > options.min_divisions + 1,
            "curved u must refine: {}",
            u_params.len()
        );
        assert_eq!(
            v_params.len(),
            options.min_divisions + 1,
            "straight v must not refine: {}",
            v_params.len()
        );

        // Every vertex lies on the unit cylinder (x^2 + y^2 == 1).
        for p in &mesh.vertices {
            let r = (p.x * p.x + p.y * p.y).sqrt();
            assert!((r - 1.0).abs() < 1e-9, "vertex off cylinder: r = {r}");
        }

        // Vertex normals match the analytic outward cylinder normal (radial in
        // XY, zero Z).
        for (p, n) in mesh.vertices.iter().zip(mesh.normals.iter()) {
            let analytic = Vector3::new(p.x, p.y, 0.0).normalize();
            let dev = angle_between(&analytic, n);
            assert!(
                dev < options.normal_tolerance,
                "normal deviates by {dev} rad"
            );
        }
    }

    #[test]
    fn indices_in_bounds_and_no_degenerate_triangles() {
        for surface in [bilinear_patch(), quarter_cylinder_patch()] {
            let mesh =
                tessellate_nurbs_surface(&surface, &SurfaceTessellationOptions::default()).unwrap();
            let n = mesh.vertices.len();
            for tri in &mesh.indices {
                for &i in tri {
                    assert!((i as usize) < n, "index {i} out of bounds ({n} vertices)");
                }
                let a = mesh.vertices[tri[0] as usize];
                let b = mesh.vertices[tri[1] as usize];
                let c = mesh.vertices[tri[2] as usize];
                let area = (b - a).cross(&(c - a)).norm() * 0.5;
                assert!(area > 1e-12, "degenerate triangle, area = {area}");
            }
        }
    }
}
