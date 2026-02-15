use crate::error::Result;
use crate::math::Point3;
use crate::topology::{TopologyStore, WireId};

/// Creates a wire from a sequence of 3D points.
pub struct MakeWire {
    points: Vec<Point3>,
    close: bool,
}

impl MakeWire {
    /// Creates a new `MakeWire` operation.
    #[must_use]
    pub fn new(points: Vec<Point3>, close: bool) -> Self {
        Self { points, close }
    }

    /// Executes the operation, creating the wire in the topology store.
    ///
    /// # Errors
    ///
    /// Returns an error if the operation fails.
    pub fn execute(&self, _store: &mut TopologyStore) -> Result<WireId> {
        let _ = (&self.points, self.close);
        todo!()
    }
}
