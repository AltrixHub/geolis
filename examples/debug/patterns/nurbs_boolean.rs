//! NURBS through-cut boolean showcase (P6/P7 + F1) — the program target.
//!
//! Variant 1 builds a curved NURBS slab and a NURBS tube that passes fully
//! through it, then runs `Subtract`. The result is the curved slab with a real
//! hole punched through both curved faces and a band face forming the hole
//! wall. The whole result solid is tessellated through the standard pipeline
//! (trimmed constrained-Delaunay for the punched faces and the band), so the
//! through-cut path is exercised end-to-end.
//!
//! Variant 2 is the F1 periodic-wrap showcase: a revolved solid (vase) cut by
//! a HORIZONTAL tube. Both holes land on the same geometrically closed wall
//! face and the tool's periodic direction wraps during SSI — the case that was
//! deferred while the marcher still terminated at UV seams. The tube runs
//! along +Y so its holes sit at wall azimuths safely away from the wall's
//! parametric seam at +X (a seam-straddling hole is a typed error until
//! general boolean face splitting).

use geolis::geometry::nurbs::NurbsCurve3D;
use geolis::math::{Point3, Vector3};
use geolis::operations::boolean::Subtract;
use geolis::operations::creation::{
    MakeCurvedSlab, MakeNurbsPrism, MakeNurbsTube, MakeRevolvedSolid,
};
use geolis::tessellation::{TessellateSolid, TessellationParams};
use geolis::topology::TopologyStore;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_edges, register_face, register_label, SceneBounds};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const SLAB_COLOR: Color = Color::rgb(120, 190, 210);
const VASE_COLOR: Color = Color::rgb(210, 160, 120);
const EDGE_COLOR: Color = Color::rgb(255, 255, 255);

pub fn register(storage: &MeshStorage, bounds: &mut SceneBounds) {
    register_slab_variant(storage, bounds);
    register_revolved_variant(storage, bounds);
}

/// Variant 1: curved slab − vertical tube.
fn register_slab_variant(storage: &MeshStorage, bounds: &mut SceneBounds) {
    register_label(storage, -1.5, 8.0, "1", LABEL_SIZE, LABEL_COLOR);

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

/// Variant 2: revolved vase − horizontal tube (F1 periodic-wrap showcase),
/// placed beside variant 1.
fn register_revolved_variant(storage: &MeshStorage, bounds: &mut SceneBounds) {
    const OFFSET_X: f64 = 12.0;
    register_label(storage, OFFSET_X - 4.5, 8.0, "2", LABEL_SIZE, LABEL_COLOR);

    let mut store = TopologyStore::new();

    // Vase-like profile revolved about the Z axis through the origin (the
    // revolve axis is fixed); the tessellated mesh is shifted beside variant 1
    // afterwards, which is why this variant registers the mesh without the
    // store-based edge overlay.
    let Ok(vase) = MakeRevolvedSolid::new(vec![(2.0, 0.0), (2.4, 1.2), (2.1, 2.4), (2.6, 3.6)])
        .execute(&mut store)
    else {
        return;
    };
    // Horizontal tube along +Y through both walls at mid-height. Azimuths
    // ±π/2 avoid the revolved wall's parametric seam at +X.
    let Ok(circle) =
        NurbsCurve3D::circle(Point3::new(0.0, -4.0, 1.8), 0.5, Vector3::y(), Vector3::x())
    else {
        return;
    };
    let Ok(tube) = MakeNurbsPrism::new(circle, Vector3::new(0.0, 8.0, 0.0)).execute(&mut store)
    else {
        return;
    };

    let Ok(result) = Subtract::new(vase, tube).execute(&mut store) else {
        return;
    };

    let Ok(mut mesh) = TessellateSolid::new(result, TessellationParams::default()).execute(&store)
    else {
        return;
    };
    // Shift the whole variant beside the slab showcase.
    for v in &mut mesh.vertices {
        v.x += OFFSET_X;
    }
    register_face(storage, bounds, mesh, VASE_COLOR);
}
