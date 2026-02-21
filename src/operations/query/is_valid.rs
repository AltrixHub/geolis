use crate::topology::{SolidId, TopologyStore};

/// Validates the topological and geometric consistency of a solid.
pub struct IsValid {
    solid: SolidId,
}

impl IsValid {
    /// Creates a new `IsValid` query.
    #[must_use]
    pub fn new(solid: SolidId) -> Self {
        Self { solid }
    }

    /// Executes the validation, returning `true` if the solid is valid.
    ///
    /// Checks:
    /// - All referenced entities exist in the store
    /// - All wires in faces are closed
    /// - Outer shell has at least one face
    /// - Each edge in a closed shell is used exactly twice across all faces
    #[must_use]
    pub fn execute(&self, store: &TopologyStore) -> bool {
        self.validate(store).is_ok()
    }

    fn validate(&self, store: &TopologyStore) -> Result<(), &'static str> {
        use std::collections::HashMap;

        let solid = store.solid(self.solid).map_err(|_| "solid not found")?;
        let shell = store
            .shell(solid.outer_shell)
            .map_err(|_| "outer shell not found")?;

        if shell.faces.is_empty() {
            return Err("shell has no faces");
        }

        let mut edge_usage: HashMap<crate::topology::EdgeId, usize> = HashMap::new();

        for &face_id in &shell.faces {
            let face = store.face(face_id).map_err(|_| "face not found")?;

            // Validate outer wire
            let outer_wire = store
                .wire(face.outer_wire)
                .map_err(|_| "outer wire not found")?;
            if !outer_wire.is_closed {
                return Err("outer wire is not closed");
            }

            for oe in &outer_wire.edges {
                let edge = store.edge(oe.edge).map_err(|_| "edge not found")?;
                store
                    .vertex(edge.start)
                    .map_err(|_| "start vertex not found")?;
                store
                    .vertex(edge.end)
                    .map_err(|_| "end vertex not found")?;
                *edge_usage.entry(oe.edge).or_insert(0) += 1;
            }

            // Validate inner wires
            for &inner_wire_id in &face.inner_wires {
                let inner_wire = store
                    .wire(inner_wire_id)
                    .map_err(|_| "inner wire not found")?;
                if !inner_wire.is_closed {
                    return Err("inner wire is not closed");
                }

                for oe in &inner_wire.edges {
                    store.edge(oe.edge).map_err(|_| "edge not found")?;
                    *edge_usage.entry(oe.edge).or_insert(0) += 1;
                }
            }
        }

        // In a closed shell, every edge should be shared by exactly 2 faces
        if shell.is_closed {
            for count in edge_usage.values() {
                if *count != 2 {
                    return Err("edge not shared by exactly 2 faces in closed shell");
                }
            }
        }

        // Validate inner shells
        for &inner_shell_id in &solid.inner_shells {
            let _shell = store
                .shell(inner_shell_id)
                .map_err(|_| "inner shell not found")?;
        }

        Ok(())
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::{Point3, Vector3};
    use crate::operations::creation::{MakeFace, MakeWire};
    use crate::operations::shaping::Extrude;
    use crate::topology::TopologyStore;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    #[test]
    fn extruded_cube_is_valid() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(0.0, 1.0, 0.0)],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(wire, vec![]).execute(&mut store).unwrap();
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 1.0))
            .execute(&mut store)
            .unwrap();

        assert!(IsValid::new(solid).execute(&store));
    }

    #[test]
    fn extruded_prism_is_valid() {
        let mut store = TopologyStore::new();
        let wire = MakeWire::new(
            vec![p(0.0, 0.0, 0.0), p(3.0, 0.0, 0.0), p(1.5, 2.0, 0.0)],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(wire, vec![]).execute(&mut store).unwrap();
        let solid = Extrude::new(face, Vector3::new(0.0, 0.0, 3.0))
            .execute(&mut store)
            .unwrap();

        assert!(IsValid::new(solid).execute(&store));
    }
}
