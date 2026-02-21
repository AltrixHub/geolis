use crate::error::{OperationError, Result};
use crate::geometry::curve::{Curve, Line};
use crate::math::{Point3, Vector3};
use crate::topology::{EdgeCurve, EdgeData, EdgeId, TopologyStore, VertexData};

/// Offsets a 2D curve (edge) by a given distance.
///
/// For line edges, creates a parallel line offset perpendicular to the edge
/// direction. Positive distance = left side, negative = right side.
pub struct CurveOffset2D {
    edge: EdgeId,
    distance: f64,
}

impl CurveOffset2D {
    /// Creates a new `CurveOffset2D` operation.
    #[must_use]
    pub fn new(edge: EdgeId, distance: f64) -> Self {
        Self { edge, distance }
    }

    /// Executes the offset, creating the result edge in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if the edge cannot be offset (e.g. arc collapses).
    pub fn execute(&self, store: &mut TopologyStore) -> Result<EdgeId> {
        let edge = store.edge(self.edge)?;
        let curve = edge.curve.clone();
        let t_start = edge.t_start;
        let t_end = edge.t_end;

        match &curve {
            EdgeCurve::Line(line) => offset_line(store, line, t_start, t_end, self.distance),
            EdgeCurve::Arc(arc) => offset_arc(store, arc, t_start, t_end, self.distance),
            EdgeCurve::Circle(_) | EdgeCurve::Ellipse(_) => {
                todo!("CurveOffset2D for Circle/Ellipse")
            }
        }
    }
}

/// Offsets a line edge by computing a perpendicular displacement.
fn offset_line(
    store: &mut TopologyStore,
    line: &crate::geometry::curve::Line,
    t_start: f64,
    t_end: f64,
    distance: f64,
) -> Result<EdgeId> {
    let dir = line.direction();
    // Left-perpendicular in XY plane: rotate 90° CCW
    let offset_dir = Vector3::new(-dir.y, dir.x, 0.0);
    let offset = offset_dir * distance;

    let start_point = line.evaluate(t_start)? + offset;
    let end_point = line.evaluate(t_end)? + offset;

    let start_v = store.add_vertex(VertexData::new(start_point));
    let end_v = store.add_vertex(VertexData::new(end_point));

    let new_dir = end_point - start_point;
    let new_t_end = new_dir.norm();
    let new_line = Line::new(start_point, new_dir)?;

    Ok(store.add_edge(EdgeData {
        start: start_v,
        end: end_v,
        curve: EdgeCurve::Line(new_line),
        t_start: 0.0,
        t_end: new_t_end,
    }))
}

/// Offsets an arc edge by adjusting the radius.
#[allow(clippy::similar_names)]
fn offset_arc(
    store: &mut TopologyStore,
    arc: &crate::geometry::curve::Arc,
    t_start: f64,
    t_end: f64,
    distance: f64,
) -> Result<EdgeId> {
    use crate::math::arc_2d::offset_arc_segment;

    let center = arc.center();
    let start_pt = arc.evaluate(t_start)?;
    let end_pt = arc.evaluate(t_end)?;

    // Compute bulge from the arc parameters
    let sweep = t_end - t_start;
    let bulge = (sweep / 4.0).tan();

    let result = offset_arc_segment(start_pt.x, start_pt.y, end_pt.x, end_pt.y, bulge, distance);

    match result {
        Some((ox0, oy0, ox1, oy1, _new_bulge)) => {
            let new_start = Point3::new(ox0, oy0, start_pt.z);
            let new_end = Point3::new(ox1, oy1, end_pt.z);
            let start_v = store.add_vertex(VertexData::new(new_start));
            let end_v = store.add_vertex(VertexData::new(new_end));

            // Reconstruct arc from offset endpoints
            let to_start = new_start - center;
            let new_radius = to_start.norm();
            let ref_dir = to_start / new_radius;

            let new_arc = crate::geometry::curve::Arc::new(
                *center,
                new_radius,
                *arc.normal(),
                ref_dir,
                t_start,
                t_end,
            )?;

            Ok(store.add_edge(EdgeData {
                start: start_v,
                end: end_v,
                curve: EdgeCurve::Arc(new_arc),
                t_start,
                t_end,
            }))
        }
        None => Err(
            OperationError::Failed("arc offset collapsed (radius <= 0)".into()).into(),
        ),
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use crate::operations::creation::MakeWire;
    use crate::topology::TopologyStore;

    #[test]
    fn offset_horizontal_line_upward() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 0.0, 0.0)],
            false,
        )
        .execute(&mut store)
        .unwrap();
        let edge_id = store.wire(wire).unwrap().edges[0].edge;

        let result = CurveOffset2D::new(edge_id, 2.0)
            .execute(&mut store)
            .unwrap();

        let edge = store.edge(result).unwrap();
        let start = store.vertex(edge.start).unwrap().point;
        let end = store.vertex(edge.end).unwrap().point;
        // Line going right: left offset = upward (positive Y)
        assert!((start.x).abs() < 1e-10);
        assert!((start.y - 2.0).abs() < 1e-10);
        assert!((end.x - 10.0).abs() < 1e-10);
        assert!((end.y - 2.0).abs() < 1e-10);
    }

    #[test]
    fn offset_line_right() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 0.0, 0.0)],
            false,
        )
        .execute(&mut store)
        .unwrap();
        let edge_id = store.wire(wire).unwrap().edges[0].edge;

        let result = CurveOffset2D::new(edge_id, -2.0)
            .execute(&mut store)
            .unwrap();

        let edge = store.edge(result).unwrap();
        let start = store.vertex(edge.start).unwrap().point;
        // Right offset = downward (negative Y)
        assert!((start.y + 2.0).abs() < 1e-10);
    }
}
