use crate::error::Result;
use crate::geometry::curve::Curve;
use crate::math::Point3;
use crate::topology::{EdgeCurve, EdgeId, TopologyStore};

/// Result of a curve-curve intersection.
#[derive(Debug, Clone, Copy)]
pub struct IntersectionResult {
    /// The intersection point.
    pub point: Point3,
    /// Parameter on the first curve.
    pub t1: f64,
    /// Parameter on the second curve.
    pub t2: f64,
}

/// Computes intersections between two curves.
pub struct CurveCurveIntersect {
    edge_a: EdgeId,
    edge_b: EdgeId,
}

impl CurveCurveIntersect {
    /// Creates a new `CurveCurveIntersect` query.
    #[must_use]
    pub fn new(edge_a: EdgeId, edge_b: EdgeId) -> Self {
        Self { edge_a, edge_b }
    }

    /// Executes the query, returning all intersection points.
    ///
    /// # Errors
    ///
    /// Returns an error if any edge is not found.
    #[allow(clippy::similar_names)]
    pub fn execute(&self, store: &TopologyStore) -> Result<Vec<IntersectionResult>> {
        let edge_a = store.edge(self.edge_a)?;
        let edge_b = store.edge(self.edge_b)?;

        let curve_a = &edge_a.curve;
        let curve_b = &edge_b.curve;
        let (ta_start, ta_end) = (edge_a.t_start, edge_a.t_end);
        let (tb_start, tb_end) = (edge_b.t_start, edge_b.t_end);

        match (curve_a, curve_b) {
            (EdgeCurve::Line(la), EdgeCurve::Line(lb)) => {
                intersect_line_line(la, ta_start, ta_end, lb, tb_start, tb_end)
            }
            (EdgeCurve::Line(la), EdgeCurve::Arc(ab)) => {
                intersect_line_arc(la, ta_start, ta_end, ab, tb_start, tb_end)
            }
            (EdgeCurve::Arc(aa), EdgeCurve::Line(lb)) => {
                // Swap and reverse parameters in results
                let results = intersect_line_arc(lb, tb_start, tb_end, aa, ta_start, ta_end)?;
                Ok(results
                    .into_iter()
                    .map(|r| IntersectionResult {
                        point: r.point,
                        t1: r.t2,
                        t2: r.t1,
                    })
                    .collect())
            }
            (EdgeCurve::Arc(aa), EdgeCurve::Arc(ab)) => {
                Ok(intersect_arc_arc(aa, ta_start, ta_end, ab, tb_start, tb_end))
            }
        }
    }
}

/// Line-Line intersection (segments).
#[allow(clippy::similar_names)]
fn intersect_line_line(
    la: &crate::geometry::curve::Line,
    ta_start: f64,
    ta_end: f64,
    lb: &crate::geometry::curve::Line,
    tb_start: f64,
    tb_end: f64,
) -> Result<Vec<IntersectionResult>> {
    use crate::math::intersect_2d::segment_segment_intersect_2d;

    let a0 = la.evaluate(ta_start)?;
    let a1 = la.evaluate(ta_end)?;
    let b0 = lb.evaluate(tb_start)?;
    let b1 = lb.evaluate(tb_end)?;

    let mut results = Vec::new();
    if let Some((pt, t, u)) = segment_segment_intersect_2d(&a0, &a1, &b0, &b1) {
        // Map from [0, 1] back to edge parameter space
        let t1 = ta_start + t * (ta_end - ta_start);
        let t2 = tb_start + u * (tb_end - tb_start);
        results.push(IntersectionResult {
            point: pt,
            t1,
            t2,
        });
    }
    Ok(results)
}

/// Line-Arc intersection (bounded).
#[allow(clippy::similar_names)]
fn intersect_line_arc(
    line: &crate::geometry::curve::Line,
    tl_start: f64,
    tl_end: f64,
    arc: &crate::geometry::curve::Arc,
    ta_start: f64,
    ta_end: f64,
) -> Result<Vec<IntersectionResult>> {
    use crate::math::intersect_2d::line_arc_intersect_2d;

    let l0 = line.evaluate(tl_start)?;
    let l1 = line.evaluate(tl_end)?;

    let center = arc.center();
    let radius = arc.radius();
    let sweep = ta_end - ta_start;

    let hits = line_arc_intersect_2d(
        l0.x, l0.y, l1.x, l1.y, center.x, center.y, radius, ta_start, sweep,
    );

    let mut results = Vec::new();
    for ((px, py), t_seg, t_arc) in hits {
        let t1 = tl_start + t_seg * (tl_end - tl_start);
        let t2 = ta_start + t_arc * sweep;
        results.push(IntersectionResult {
            point: Point3::new(px, py, l0.z),
            t1,
            t2,
        });
    }
    Ok(results)
}

/// Arc-Arc intersection.
#[allow(clippy::similar_names)]
fn intersect_arc_arc(
    aa: &crate::geometry::curve::Arc,
    ta_start: f64,
    ta_end: f64,
    ab: &crate::geometry::curve::Arc,
    tb_start: f64,
    tb_end: f64,
) -> Vec<IntersectionResult> {
    use crate::math::intersect_2d::arc_arc_intersect_2d;

    let ca = aa.center();
    let cb = ab.center();
    let ra = aa.radius();
    let rb = ab.radius();
    let sweep_a = ta_end - ta_start;
    let sweep_b = tb_end - tb_start;

    let hits = arc_arc_intersect_2d(
        ca.x, ca.y, ra, ta_start, sweep_a, cb.x, cb.y, rb, tb_start, sweep_b,
    );

    hits.into_iter()
        .map(|((px, py), t1_frac, t2_frac)| IntersectionResult {
            point: Point3::new(px, py, ca.z),
            t1: ta_start + t1_frac * sweep_a,
            t2: tb_start + t2_frac * sweep_b,
        })
        .collect()
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use crate::operations::creation::MakeWire;
    use crate::topology::TopologyStore;

    #[test]
    fn two_crossing_lines() {
        let mut store = TopologyStore::new();

        // Line A: (0, 0) -> (2, 2) diagonal
        let wire_a = MakeWire::new(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 2.0, 0.0)],
            false,
        )
        .execute(&mut store)
        .unwrap();
        let edge_a = store.wire(wire_a).unwrap().edges[0].edge;

        // Line B: (0, 2) -> (2, 0) anti-diagonal
        let wire_b = MakeWire::new(
            vec![Point3::new(0.0, 2.0, 0.0), Point3::new(2.0, 0.0, 0.0)],
            false,
        )
        .execute(&mut store)
        .unwrap();
        let edge_b = store.wire(wire_b).unwrap().edges[0].edge;

        let results = CurveCurveIntersect::new(edge_a, edge_b)
            .execute(&store)
            .unwrap();

        assert_eq!(results.len(), 1);
        assert!((results[0].point.x - 1.0).abs() < 1e-6);
        assert!((results[0].point.y - 1.0).abs() < 1e-6);
    }

    #[test]
    fn parallel_lines_no_intersection() {
        let mut store = TopologyStore::new();

        let wire_a = MakeWire::new(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)],
            false,
        )
        .execute(&mut store)
        .unwrap();
        let edge_a = store.wire(wire_a).unwrap().edges[0].edge;

        let wire_b = MakeWire::new(
            vec![Point3::new(0.0, 1.0, 0.0), Point3::new(2.0, 1.0, 0.0)],
            false,
        )
        .execute(&mut store)
        .unwrap();
        let edge_b = store.wire(wire_b).unwrap().edges[0].edge;

        let results = CurveCurveIntersect::new(edge_a, edge_b)
            .execute(&store)
            .unwrap();

        assert!(results.is_empty());
    }
}
