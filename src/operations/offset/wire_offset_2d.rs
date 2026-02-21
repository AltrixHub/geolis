use crate::error::{OperationError, Result};
use crate::geometry::pline::Pline;
use crate::math::Point3;
use crate::operations::creation::MakeWire;
use crate::topology::{TopologyStore, WireId};

use super::PlineOffset2D;

/// Offsets a 2D wire by a given distance.
///
/// Converts the wire to a [`Pline`], delegates to [`PlineOffset2D`],
/// then converts the result back to a wire.
pub struct WireOffset2D {
    wire: WireId,
    distance: f64,
}

impl WireOffset2D {
    /// Creates a new `WireOffset2D` operation.
    #[must_use]
    pub fn new(wire: WireId, distance: f64) -> Self {
        Self { wire, distance }
    }

    /// Executes the offset, creating the result wire in the topology store.
    ///
    /// Returns the first (largest) offset result. For closed wires that split
    /// into multiple loops, only the first is returned.
    ///
    /// # Errors
    ///
    /// Returns an error if the wire cannot be offset (e.g. collapses entirely).
    pub fn execute(&self, store: &mut TopologyStore) -> Result<WireId> {
        let wire = store.wire(self.wire)?;
        let is_closed = wire.is_closed;
        let edges = wire.edges.clone();

        // Collect 2D points from the wire
        let mut points = Vec::with_capacity(edges.len() + 1);
        for oe in &edges {
            let edge = store.edge(oe.edge)?;
            let vid = if oe.forward { edge.start } else { edge.end };
            points.push(store.vertex(vid)?.point);
        }

        // For open wires, add the last vertex
        if !is_closed {
            if let Some(last_oe) = edges.last() {
                let edge = store.edge(last_oe.edge)?;
                let vid = if last_oe.forward { edge.end } else { edge.start };
                points.push(store.vertex(vid)?.point);
            }
        }

        // Convert to Pline
        let pline = Pline::from_points(&points, is_closed);

        // Execute offset
        let results = PlineOffset2D::new(pline, self.distance).execute()?;

        if results.is_empty() {
            return Err(
                OperationError::Failed("wire offset collapsed entirely".into()).into(),
            );
        }

        // Convert first result back to wire
        let result_pline = &results[0];
        let mut result_points: Vec<Point3> = result_pline.to_points(0.01);

        // For closed wires, the Pline output may include a closing point that
        // duplicates the first point — remove it to avoid MakeWire error.
        if is_closed && result_points.len() >= 3 {
            let first = result_points[0];
            if let Some(last) = result_points.last() {
                let dist = (last - first).norm();
                if dist < crate::math::TOLERANCE * 100.0 {
                    result_points.pop();
                }
            }
        }

        if result_points.len() < 2 {
            return Err(
                OperationError::Failed("offset result has too few points".into()).into(),
            );
        }

        MakeWire::new(result_points, is_closed).execute(store)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use crate::operations::creation::MakeWire;
    use crate::topology::TopologyStore;

    fn p(x: f64, y: f64) -> Point3 {
        Point3::new(x, y, 0.0)
    }

    #[test]
    fn offset_square_inward() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![p(0.0, 0.0), p(10.0, 0.0), p(10.0, 10.0), p(0.0, 10.0)],
            true,
        )
        .execute(&mut store)
        .unwrap();

        let result = WireOffset2D::new(wire, 1.0)
            .execute(&mut store)
            .unwrap();

        let result_wire = store.wire(result).unwrap();
        assert!(result_wire.is_closed);
        // An inward offset of 1.0 on a 10x10 square should give 4 edges (8x8)
        assert_eq!(result_wire.edges.len(), 4);
    }
}
