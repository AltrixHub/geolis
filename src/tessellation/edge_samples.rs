//! Per-edge boundary sampling — the single source of truth for "the polyline
//! of this edge".
//!
//! Shared-edge `BRep` topology makes boundary conformance structural: adjacent
//! faces reference the same [`EdgeId`], the edge is sampled into a polyline
//! exactly once per [`EdgeSampleCache`], and every consumer (planar cap CDT,
//! NURBS boundary CDT via face pcurves) looks the samples up instead of
//! re-deriving them from geometry. Two faces meeting at an edge therefore emit
//! bit-identical 3D boundary vertices by construction.
//!
//! Sampling rules per curve kind:
//! - NURBS: chord-adaptive curve-intrinsic parameters
//!   ([`tessellate_nurbs_curve_params`] at [`BOUNDARY_CHORD_TOLERANCE`]) — the
//!   same algorithm design (a) used, so pre-existing geometric conformance is
//!   preserved bit-for-bit.
//! - Line: the two endpoints.
//! - Arc / Circle / Ellipse: the sagitta-bounded segment count previously
//!   local to the planar face path.

use std::collections::HashMap;

use crate::error::Result;
use crate::geometry::curve::Curve;
use crate::math::Point3;
use crate::topology::{EdgeCurve, EdgeId, TopologyStore};

use super::tessellate_nurbs::{
    tessellate_nurbs_curve_params, CurveTessellationOptions, BOUNDARY_CHORD_TOLERANCE,
};
use super::TessellationParams;

/// Chord-adaptive samples of one edge: parameters on the edge curve plus the
/// evaluated 3D points, in the edge's natural (`t_start` → `t_end`) order.
#[derive(Debug, Clone)]
pub(crate) struct EdgeSamples {
    /// Curve parameters, ascending from `t_start` to `t_end` (both included).
    pub params: Vec<f64>,
    /// `curve.point_at(params[i])`, synchronized with `params`.
    pub points: Vec<Point3>,
}

/// Per-solid cache: every edge is sampled exactly once, so every face that
/// references the edge consumes the identical polyline.
#[derive(Debug, Default)]
pub(crate) struct EdgeSampleCache {
    params: TessellationParams,
    samples: HashMap<EdgeId, EdgeSamples>,
}

impl EdgeSampleCache {
    /// Creates a cache sampling analytic edges under `params`.
    pub fn new(params: TessellationParams) -> Self {
        Self {
            params,
            samples: HashMap::new(),
        }
    }

    /// Returns the samples for `edge`, computing them on first access.
    ///
    /// # Errors
    ///
    /// Propagates curve evaluation errors.
    pub fn get(&mut self, store: &TopologyStore, edge: EdgeId) -> Result<&EdgeSamples> {
        if !self.samples.contains_key(&edge) {
            let computed = sample_edge(store, edge, &self.params)?;
            self.samples.insert(edge, computed);
        }
        // The entry was just inserted (or already present); index is safe.
        Ok(&self.samples[&edge])
    }
}

/// Samples one edge according to the per-kind rules in the module docs.
fn sample_edge(
    store: &TopologyStore,
    edge: EdgeId,
    params: &TessellationParams,
) -> Result<EdgeSamples> {
    let data = store.edge(edge)?;
    let (t_start, t_end) = (data.t_start, data.t_end);

    let ts: Vec<f64> = match &data.curve {
        EdgeCurve::Line(_) => vec![t_start, t_end],
        EdgeCurve::Arc(arc) => uniform_params(
            t_start,
            t_end,
            sagitta_segments(arc.radius(), t_start, t_end, params),
        ),
        EdgeCurve::Circle(circle) => uniform_params(
            t_start,
            t_end,
            sagitta_segments(circle.radius(), t_start, t_end, params),
        ),
        EdgeCurve::Ellipse(ellipse) => uniform_params(
            t_start,
            t_end,
            sagitta_segments(ellipse.semi_major(), t_start, t_end, params),
        ),
        // A degree-1 NURBS curve IS its control polygon: sampling exactly at
        // the knot breakpoints reproduces the curve losslessly (the same rule
        // the trim-loop sampler applies to degree-1 pcurves). This keeps SSI
        // trace edges — degree-1 polylines through the marcher's samples —
        // conformal with the faces whose trim loops carry the identical
        // sample points.
        EdgeCurve::Nurbs(nurbs) if nurbs.degree() == 1 => {
            clip_params(breakpoint_params(nurbs.knots()), t_start, t_end)
        }
        EdgeCurve::Nurbs(nurbs) => {
            let options = CurveTessellationOptions {
                chord_tolerance: BOUNDARY_CHORD_TOLERANCE,
                max_depth: 16,
            };
            clip_params(
                tessellate_nurbs_curve_params(nurbs, &options)?,
                t_start,
                t_end,
            )
        }
    };

    let mut points = Vec::with_capacity(ts.len());
    for &t in &ts {
        points.push(evaluate(&data.curve, t)?);
    }
    Ok(EdgeSamples { params: ts, points })
}

/// Evaluates the edge curve at `t` regardless of the curve kind.
fn evaluate(curve: &EdgeCurve, t: f64) -> Result<Point3> {
    match curve {
        EdgeCurve::Line(c) => c.evaluate(t),
        EdgeCurve::Arc(c) => c.evaluate(t),
        EdgeCurve::Circle(c) => c.evaluate(t),
        EdgeCurve::Ellipse(c) => c.evaluate(t),
        EdgeCurve::Nurbs(c) => c.point_at(t),
    }
}

/// The distinct knot values of a knot vector — the breakpoints of a degree-1
/// curve, where its control polygon vertices sit.
fn breakpoint_params(knots: &crate::geometry::nurbs::KnotVector) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();
    for &k in knots.as_slice() {
        if out.last().is_none_or(|&last| k > last) {
            out.push(k);
        }
    }
    out
}

/// `n + 1` uniformly spaced parameters from `t_start` to `t_end` inclusive.
fn uniform_params(t_start: f64, t_end: f64, segments: usize) -> Vec<f64> {
    let n = segments.max(1);
    let mut ts = Vec::with_capacity(n + 1);
    for i in 0..=n {
        #[allow(clippy::cast_precision_loss)]
        let frac = i as f64 / n as f64;
        ts.push(t_start + (t_end - t_start) * frac);
    }
    ts
}

/// Sagitta-bounded segment count for circular-ish edges (the rule previously
/// local to the planar face path): the chord deviation of each segment stays
/// below `params.tolerance`.
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn sagitta_segments(radius: f64, t_start: f64, t_end: f64, params: &TessellationParams) -> usize {
    let sweep = (t_end - t_start).abs();
    if radius > params.tolerance {
        let half_angle = (1.0 - params.tolerance / radius).acos();
        let computed = (sweep / (2.0 * half_angle)).ceil() as usize;
        computed.clamp(params.min_segments, params.max_segments)
    } else {
        params.min_segments
    }
}

/// Restricts full-domain chord-adaptive parameters to `[t_start, t_end]`,
/// keeping the exact endpoints. Full-domain edges (the common case from the
/// creation ops) pass through unchanged.
fn clip_params(full: Vec<f64>, t_start: f64, t_end: f64) -> Vec<f64> {
    let (lo, hi) = if t_start <= t_end {
        (t_start, t_end)
    } else {
        (t_end, t_start)
    };
    let eps = 1e-12 * (hi - lo).abs().max(1.0);
    let mut ts = vec![lo];
    ts.extend(full.into_iter().filter(|&t| t > lo + eps && t < hi - eps));
    ts.push(hi);
    ts
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::curve::Line;
    use crate::geometry::nurbs::NurbsCurve3D;
    use crate::math::Vector3;
    use crate::topology::{EdgeData, VertexData};

    fn store_with_edge(curve: EdgeCurve, t_start: f64, t_end: f64) -> (TopologyStore, EdgeId) {
        let mut store = TopologyStore::new();
        let p = evaluate(&curve, t_start).unwrap();
        let q = evaluate(&curve, t_end).unwrap();
        let start = store.add_vertex(VertexData { point: p });
        let end = store.add_vertex(VertexData { point: q });
        let edge = store.add_edge(EdgeData {
            start,
            end,
            curve,
            t_start,
            t_end,
        });
        (store, edge)
    }

    #[test]
    fn nurbs_edge_matches_curve_intrinsic_params() {
        // The cache must reproduce the design-(a) parameters exactly so
        // pre-existing geometric conformance is preserved bit-for-bit.
        let circle =
            NurbsCurve3D::circle(Point3::origin(), 0.8, Vector3::z(), Vector3::x()).unwrap();
        let (t0, t1) = circle.parameter_domain();
        let expected = tessellate_nurbs_curve_params(
            &circle,
            &CurveTessellationOptions {
                chord_tolerance: BOUNDARY_CHORD_TOLERANCE,
                max_depth: 16,
            },
        )
        .unwrap();

        let (store, edge) = store_with_edge(EdgeCurve::Nurbs(circle), t0, t1);
        let mut cache = EdgeSampleCache::new(TessellationParams::default());
        let samples = cache.get(&store, edge).unwrap();
        assert_eq!(
            samples.params, expected,
            "cache must reuse the design-(a) params"
        );
        assert_eq!(samples.params.len(), samples.points.len());
    }

    #[test]
    fn repeated_get_returns_identical_samples() {
        let circle =
            NurbsCurve3D::circle(Point3::origin(), 1.0, Vector3::z(), Vector3::x()).unwrap();
        let (t0, t1) = circle.parameter_domain();
        let (store, edge) = store_with_edge(EdgeCurve::Nurbs(circle), t0, t1);
        let mut cache = EdgeSampleCache::new(TessellationParams::default());
        let first: Vec<f64> = cache.get(&store, edge).unwrap().params.clone();
        let second: Vec<f64> = cache.get(&store, edge).unwrap().params.clone();
        assert_eq!(first, second);
    }

    /// A degree-1 NURBS polyline edge samples exactly at its breakpoints
    /// (the control polygon vertices) — losslessly and deterministically, so
    /// SSI trace edges conform with trim loops built from the same points.
    #[test]
    fn degree_one_nurbs_edge_samples_at_breakpoints() {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.5, 0.0),
            Point3::new(2.0, 0.2, 0.3),
            Point3::new(3.0, 1.0, 0.1),
        ];
        let poly = NurbsCurve3D::polyline(&pts).unwrap();
        let (t0, t1) = poly.parameter_domain();
        let (store, edge) = store_with_edge(EdgeCurve::Nurbs(poly), t0, t1);
        let mut cache = EdgeSampleCache::new(TessellationParams::default());
        let samples = cache.get(&store, edge).unwrap();
        assert_eq!(
            samples.points.len(),
            pts.len(),
            "one sample per polyline vertex"
        );
        for (sample, expected) in samples.points.iter().zip(&pts) {
            assert!(
                (*sample - *expected).norm() < 1e-12,
                "breakpoint sample must reproduce the polyline vertex"
            );
        }
    }

    #[test]
    fn line_edge_yields_exactly_its_endpoints() {
        // `Line` normalizes its direction, so `t` is arc length.
        let line = Line::new(Point3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)).unwrap();
        let (store, edge) = store_with_edge(EdgeCurve::Line(line), 0.0, 2.0);
        let mut cache = EdgeSampleCache::new(TessellationParams::default());
        let samples = cache.get(&store, edge).unwrap();
        assert_eq!(samples.points.len(), 2);
        assert!((samples.points[0] - Point3::new(0.0, 0.0, 0.0)).norm() < 1e-12);
        assert!((samples.points[1] - Point3::new(2.0, 0.0, 0.0)).norm() < 1e-12);
    }

    #[test]
    fn circle_edge_respects_sagitta_bound() {
        use crate::geometry::curve::Circle;
        let circle = Circle::new(Point3::origin(), 1.0, Vector3::z(), Vector3::x()).unwrap();
        let (store, edge) = store_with_edge(EdgeCurve::Circle(circle), 0.0, std::f64::consts::TAU);
        let params = TessellationParams::default();
        let expected = sagitta_segments(1.0, 0.0, std::f64::consts::TAU, &params) + 1;
        let mut cache = EdgeSampleCache::new(params);
        let samples = cache.get(&store, edge).unwrap();
        assert_eq!(samples.points.len(), expected);
    }
}
