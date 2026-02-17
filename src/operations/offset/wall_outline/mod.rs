mod decompose;
mod junction;
mod offset_edges;
mod trace;

use crate::error::{OperationError, Result};
use crate::geometry::pline::Pline;

/// Generates wall outlines from a centerline network.
///
/// Given a `Pline` representing wall centerlines (potentially with
/// self-intersecting paths), produces closed outline polygons at the
/// specified half-width distance.
#[derive(Debug)]
pub struct WallOutline2D {
    pline: Pline,
    half_width: f64,
}

impl WallOutline2D {
    /// Creates a new wall outline operation.
    #[must_use]
    pub fn new(pline: Pline, half_width: f64) -> Self {
        Self { pline, half_width }
    }

    /// Executes the wall outline generation.
    ///
    /// # Errors
    ///
    /// Returns `OperationError::InvalidInput` if the polyline has fewer than
    /// 2 vertices, or `OperationError::Failed` if no outline can be generated.
    pub fn execute(&self) -> Result<Vec<Pline>> {
        if self.pline.vertices.len() < 2 {
            return Err(OperationError::InvalidInput(
                "at least 2 vertices required for wall outline".to_owned(),
            )
            .into());
        }

        if self.half_width.abs() < crate::math::TOLERANCE {
            return Ok(vec![self.pline.clone()]);
        }

        // Step 1: Decompose into unique segments.
        let segments = decompose::decompose(&self.pline);
        if segments.is_empty() {
            return Err(
                OperationError::Failed("no valid segments in pline".to_owned()).into(),
            );
        }

        // Step 2: Detect junctions and split segments.
        let network = junction::build_network(&segments);

        // Step 3: Generate offset edges with junction resolution.
        let edges = offset_edges::build(&network, self.half_width);

        // Step 4: Trace outer boundaries.
        let outlines = trace::trace_boundaries(&edges);

        if outlines.is_empty() {
            return Err(
                OperationError::Failed("wall outline trace produced no results".to_owned())
                    .into(),
            );
        }

        Ok(outlines)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::pline::PlineVertex;

    fn double_cross_pline() -> Pline {
        Pline {
            vertices: vec![
                PlineVertex::line(3.0, 0.0),
                PlineVertex::line(3.0, 10.0),
                PlineVertex::line(3.0, 7.0),
                PlineVertex::line(0.0, 7.0),
                PlineVertex::line(10.0, 7.0),
                PlineVertex::line(7.0, 7.0),
                PlineVertex::line(7.0, 10.0),
                PlineVertex::line(7.0, 0.0),
                PlineVertex::line(7.0, 3.0),
                PlineVertex::line(10.0, 3.0),
                PlineVertex::line(0.0, 3.0),
            ],
            closed: false,
        }
    }

    #[test]
    fn debug_double_cross_wall_outline() {
        let wall = WallOutline2D::new(double_cross_pline(), 0.3);
        let result = wall.execute().unwrap();
        eprintln!("=== WallOutline2D double-cross d=0.3 ===");
        eprintln!("Result: {} polygons", result.len());
        for (i, p) in result.iter().enumerate() {
            eprintln!("  Poly {i}: {} verts, closed={}", p.vertices.len(), p.closed);
            for (j, v) in p.vertices.iter().enumerate() {
                eprintln!("    [{j}] ({:.4}, {:.4})", v.x, v.y);
            }
        }
    }
}
