use crate::error::Result;
use crate::math::Point3;
use crate::topology::{EdgeCurve, EdgeId, TopologyStore};

/// Result of a closest point query.
#[derive(Debug, Clone, Copy)]
pub struct ClosestPointResult {
    /// The closest point on the curve.
    pub point: Point3,
    /// The parameter value at the closest point.
    pub parameter: f64,
    /// The distance from the query point to the closest point.
    pub distance: f64,
}

/// Finds the closest point on a curve to a given point.
pub struct ClosestPointOnCurve {
    edge: EdgeId,
    point: Point3,
}

impl ClosestPointOnCurve {
    /// Creates a new `ClosestPointOnCurve` query.
    #[must_use]
    pub fn new(edge: EdgeId, point: Point3) -> Self {
        Self { edge, point }
    }

    /// Executes the query, returning the closest point result.
    ///
    /// For `Line` edges, computes the analytical projection.
    /// For `Arc` edges, finds the closest angle and clamps to the arc domain.
    ///
    /// # Errors
    ///
    /// Returns an error if the edge is not found.
    pub fn execute(&self, store: &TopologyStore) -> Result<ClosestPointResult> {
        let edge = store.edge(self.edge)?;
        match &edge.curve {
            EdgeCurve::Line(line) => {
                closest_point_on_line(line, edge.t_start, edge.t_end, &self.point)
            }
            EdgeCurve::Arc(arc) => {
                closest_point_on_arc(arc, edge.t_start, edge.t_end, &self.point)
            }
            EdgeCurve::Circle(_) | EdgeCurve::Ellipse(_) => {
                todo!("ClosestPointOnCurve for Circle/Ellipse")
            }
        }
    }
}

/// Finds the closest point on a bounded line segment.
fn closest_point_on_line(
    line: &crate::geometry::curve::Line,
    t_start: f64,
    t_end: f64,
    point: &Point3,
) -> Result<ClosestPointResult> {
    use crate::geometry::curve::Curve;

    let origin = line.origin();
    let dir = line.direction();

    // Project point onto line: t = dot(point - origin, dir)
    let to_point = point - origin;
    let t = to_point.dot(dir).clamp(t_start, t_end);

    let closest = line.evaluate(t)?;
    let distance = (point - closest).norm();

    Ok(ClosestPointResult {
        point: closest,
        parameter: t,
        distance,
    })
}

/// Finds the closest point on a bounded arc.
fn closest_point_on_arc(
    arc: &crate::geometry::curve::Arc,
    t_start: f64,
    t_end: f64,
    point: &Point3,
) -> Result<ClosestPointResult> {
    use crate::geometry::curve::Curve;

    let center = arc.center();
    let normal = arc.normal();

    // Project the query point onto the arc's plane
    let to_point = point - center;
    let in_plane = to_point - normal * to_point.dot(normal);
    let dist_from_center = in_plane.norm();

    if dist_from_center < crate::math::TOLERANCE {
        // Point is at the center of the arc; the entire arc is equidistant.
        // Return the start point.
        let closest = arc.evaluate(t_start)?;
        let distance = (point - closest).norm();
        return Ok(ClosestPointResult {
            point: closest,
            parameter: t_start,
            distance,
        });
    }

    // Find the angle of the projection
    // Use atan2 in the arc's local coordinate system
    // Sample the arc at multiple points to find the closest
    let n_samples = 64;
    let mut best_t = t_start;
    let mut best_dist = f64::INFINITY;

    for i in 0..=n_samples {
        #[allow(clippy::cast_precision_loss)]
        let frac = f64::from(i) / f64::from(n_samples);
        let t = t_start + frac * (t_end - t_start);
        let pt = arc.evaluate(t)?;
        let d = (point - pt).norm();
        if d < best_dist {
            best_dist = d;
            best_t = t;
        }
    }

    // Refine with bisection around the best sample
    #[allow(clippy::cast_precision_loss)]
    let dt = (t_end - t_start) / f64::from(n_samples);
    let mut lo = (best_t - dt).max(t_start);
    let mut hi = (best_t + dt).min(t_end);

    for _ in 0..50 {
        let mid1 = lo + (hi - lo) / 3.0;
        let mid2 = hi - (hi - lo) / 3.0;
        let d1 = (point - arc.evaluate(mid1)?).norm();
        let d2 = (point - arc.evaluate(mid2)?).norm();
        if d1 < d2 {
            hi = mid2;
        } else {
            lo = mid1;
        }
    }

    #[allow(clippy::manual_midpoint)]
    let best_t = (lo + hi) / 2.0;
    let closest = arc.evaluate(best_t)?;
    let distance = (point - closest).norm();

    Ok(ClosestPointResult {
        point: closest,
        parameter: best_t,
        distance,
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use crate::operations::creation::MakeWire;
    use crate::topology::TopologyStore;

    #[test]
    fn closest_point_on_line_perpendicular() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 0.0, 0.0)],
            false,
        )
        .execute(&mut store)
        .unwrap();
        let edge_id = store.wire(wire).unwrap().edges[0].edge;

        let result = ClosestPointOnCurve::new(edge_id, Point3::new(5.0, 3.0, 0.0))
            .execute(&store)
            .unwrap();

        assert!((result.point.x - 5.0).abs() < 1e-10);
        assert!(result.point.y.abs() < 1e-10);
        assert!((result.distance - 3.0).abs() < 1e-10);
    }

    #[test]
    fn closest_point_clamps_to_start() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 0.0, 0.0)],
            false,
        )
        .execute(&mut store)
        .unwrap();
        let edge_id = store.wire(wire).unwrap().edges[0].edge;

        let result = ClosestPointOnCurve::new(edge_id, Point3::new(-5.0, 0.0, 0.0))
            .execute(&store)
            .unwrap();

        assert!(result.point.x.abs() < 1e-10);
        assert!((result.distance - 5.0).abs() < 1e-10);
    }

    #[test]
    fn closest_point_clamps_to_end() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 0.0, 0.0)],
            false,
        )
        .execute(&mut store)
        .unwrap();
        let edge_id = store.wire(wire).unwrap().edges[0].edge;

        let result = ClosestPointOnCurve::new(edge_id, Point3::new(15.0, 0.0, 0.0))
            .execute(&store)
            .unwrap();

        assert!((result.point.x - 10.0).abs() < 1e-10);
        assert!((result.distance - 5.0).abs() < 1e-10);
    }
}
