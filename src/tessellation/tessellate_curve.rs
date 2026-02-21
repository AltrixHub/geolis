use crate::error::Result;
use crate::geometry::curve::Curve;
use crate::topology::{EdgeCurve, EdgeId, TopologyStore};

use super::{Polyline, TessellationParams};

/// Tessellates a curve (edge) into a polyline.
pub struct TessellateCurve {
    edge: EdgeId,
    params: TessellationParams,
}

impl TessellateCurve {
    /// Creates a new `TessellateCurve` operation.
    #[must_use]
    pub fn new(edge: EdgeId, params: TessellationParams) -> Self {
        Self { edge, params }
    }

    /// Executes the tessellation, returning a polyline.
    ///
    /// For `Line` edges, returns just the two endpoints.
    /// For `Arc` edges, subdivides based on tolerance and segment limits.
    ///
    /// # Errors
    ///
    /// Returns an error if the edge is not found or evaluation fails.
    pub fn execute(&self, store: &TopologyStore) -> Result<Polyline> {
        let edge = store.edge(self.edge)?;
        match &edge.curve {
            EdgeCurve::Line(line) => {
                let p0 = line.evaluate(edge.t_start)?;
                let p1 = line.evaluate(edge.t_end)?;
                Ok(Polyline {
                    points: vec![p0, p1],
                })
            }
            EdgeCurve::Arc(arc) => tessellate_arc(arc, edge.t_start, edge.t_end, &self.params),
            EdgeCurve::Circle(circle) => {
                tessellate_circular(circle.radius(), circle, edge.t_start, edge.t_end, &self.params)
            }
            EdgeCurve::Ellipse(ellipse) => {
                // Approximate with the semi-major axis for chord error calculation
                tessellate_circular(ellipse.semi_major(), ellipse, edge.t_start, edge.t_end, &self.params)
            }
        }
    }
}

/// Tessellates an arc into a polyline with adaptive subdivision.
fn tessellate_arc(
    arc: &crate::geometry::curve::Arc,
    t_start: f64,
    t_end: f64,
    params: &TessellationParams,
) -> Result<Polyline> {
    tessellate_circular(arc.radius(), arc, t_start, t_end, params)
}

/// Tessellates any circular/elliptical curve into a polyline with adaptive subdivision.
///
/// `approx_radius` is used for chord error calculation (use radius for circles,
/// semi-major axis for ellipses).
fn tessellate_circular(
    approx_radius: f64,
    curve: &dyn Curve,
    t_start: f64,
    t_end: f64,
    params: &TessellationParams,
) -> Result<Polyline> {
    let sweep = (t_end - t_start).abs();

    // Compute number of segments based on tolerance:
    // chord error = r * (1 - cos(theta/2)) where theta = sweep/n
    // So n = sweep / (2 * acos(1 - tolerance/r))
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let n = if approx_radius > params.tolerance {
        let half_angle = (1.0 - params.tolerance / approx_radius).acos();
        let computed = (sweep / (2.0 * half_angle)).ceil() as usize;
        computed.clamp(params.min_segments, params.max_segments)
    } else {
        params.min_segments
    };

    let mut points = Vec::with_capacity(n + 1);
    for i in 0..=n {
        #[allow(clippy::cast_precision_loss)]
        let frac = i as f64 / n as f64;
        let t = t_start + frac * (t_end - t_start);
        points.push(curve.evaluate(t)?);
    }

    Ok(Polyline { points })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::curve::Arc;
    use crate::math::{Point3, Vector3};
    use crate::operations::creation::MakeWire;
    use crate::topology::{EdgeCurve, EdgeData, TopologyStore, VertexData};

    #[test]
    fn line_tessellation_returns_two_points() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(5.0, 0.0, 0.0)],
            false,
        )
        .execute(&mut store)
        .unwrap();
        let edge_id = store.wire(wire).unwrap().edges[0].edge;

        let polyline = TessellateCurve::new(edge_id, TessellationParams::default())
            .execute(&store)
            .unwrap();

        assert_eq!(polyline.points.len(), 2);
        assert!((polyline.points[0].x).abs() < 1e-10);
        assert!((polyline.points[1].x - 5.0).abs() < 1e-10);
    }

    #[test]
    fn arc_tessellation_produces_multiple_points() {
        let mut store = TopologyStore::new();

        // Create a semicircular arc edge manually
        let start = Point3::new(1.0, 0.0, 0.0);
        let end = Point3::new(-1.0, 0.0, 0.0);
        let v_start = store.add_vertex(VertexData::new(start));
        let v_end = store.add_vertex(VertexData::new(end));

        let arc = Arc::new(
            Point3::new(0.0, 0.0, 0.0),
            1.0,
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 0.0),
            0.0,
            std::f64::consts::PI,
        )
        .unwrap();

        let edge_id = store.add_edge(EdgeData {
            start: v_start,
            end: v_end,
            curve: EdgeCurve::Arc(arc),
            t_start: 0.0,
            t_end: std::f64::consts::PI,
        });

        let polyline = TessellateCurve::new(edge_id, TessellationParams::default())
            .execute(&store)
            .unwrap();

        assert!(polyline.points.len() >= 5); // At least min_segments + 1
        // First point should be (1, 0, 0)
        assert!((polyline.points[0].x - 1.0).abs() < 1e-10);
        // Last point should be (-1, 0, 0)
        assert!((polyline.points.last().unwrap().x + 1.0).abs() < 1e-10);
    }
}
