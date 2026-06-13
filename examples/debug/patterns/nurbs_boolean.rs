//! NURBS through-cut boolean showcase (P6/P7) — the program target.
//!
//! Builds a curved NURBS slab and a NURBS tube that passes fully through it,
//! then runs `boolean_execute(slab, tube, Subtract)`. The result is the curved
//! slab with a real hole punched through both curved faces and a band face
//! forming the hole wall. The whole result solid is tessellated through the
//! standard pipeline (trimmed constrained-Delaunay for the punched faces and
//! the band), so the through-cut path is exercised end-to-end.
//!
//! The revolved-solid variant (a horizontal-axis tube whose band surface is
//! geometrically closed but parametrically non-periodic) is deferred: it fails
//! the v1 seam-closure criterion in the SSI marcher (see
//! `geometry::nurbs::intersect::surface_surface`).

use geolis::math::Point3;
use geolis::operations::boolean::Subtract;
use geolis::operations::creation::{MakeCurvedSlab, MakeNurbsTube};
use geolis::tessellation::{TessellateSolid, TessellationParams};
use geolis::topology::TopologyStore;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_edges, register_face, register_label, SceneBounds};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const SLAB_COLOR: Color = Color::rgb(120, 190, 210);
const EDGE_COLOR: Color = Color::rgb(255, 255, 255);

pub fn register(storage: &MeshStorage, bounds: &mut SceneBounds) {
    register_label(storage, bounds, -1.5, 8.0, "1", LABEL_SIZE, LABEL_COLOR);

    let mut store = TopologyStore::new();

    // Curved slab spanning [0,6]^2 in XY, peaking 1.5 above z=0, 1.0 thick.
    let Ok(slab) = MakeCurvedSlab::new(6.0, 0.0, 1.5, 1.0).execute(&mut store) else {
        return;
    };
    // Tube through the slab center, rising from below to above.
    let Ok(tube) = MakeNurbsTube::new(Point3::new(3.0, 3.0, -1.5), 0.8, 5.0).execute(&mut store)
    else {
        return;
    };

    let Ok(result) = Subtract::new(slab, tube).execute(&mut store) else {
        return;
    };

    let Ok(mesh) = TessellateSolid::new(result, TessellationParams::default()).execute(&store)
    else {
        return;
    };
    register_face(storage, bounds, mesh, SLAB_COLOR);

    if let Ok(solid) = store.solid(result) {
        register_edges(storage, bounds, &store, solid.outer_shell, EDGE_COLOR);
    }
}
