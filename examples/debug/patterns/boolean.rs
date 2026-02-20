use geolis::math::{Point3, Vector3};
use geolis::operations::boolean::{Intersect, Subtract, Union};
use geolis::operations::creation::{MakeFace, MakeWire};
use geolis::operations::shaping::Extrude;
use geolis::tessellation::{TessellateSolid, TessellationParams};
use geolis::topology::TopologyStore;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_edges, register_face, register_label};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const GREEN: Color = Color::rgb(100, 200, 100);
const BLUE: Color = Color::rgb(100, 150, 255);
const RED: Color = Color::rgb(255, 100, 100);

fn make_box(
    store: &mut TopologyStore,
    x: f64,
    y: f64,
    z: f64,
    dx: f64,
    dy: f64,
    dz: f64,
) -> Option<geolis::topology::SolidId> {
    let pts = vec![
        Point3::new(x, y, z),
        Point3::new(x + dx, y, z),
        Point3::new(x + dx, y + dy, z),
        Point3::new(x, y + dy, z),
    ];
    let wire = MakeWire::new(pts, true).execute(store).ok()?;
    let face = MakeFace::new(wire, vec![]).execute(store).ok()?;
    Extrude::new(face, Vector3::new(0.0, 0.0, dz))
        .execute(store)
        .ok()
}

fn render_solid(
    storage: &MeshStorage,
    store: &TopologyStore,
    solid: geolis::topology::SolidId,
    mesh_color: Color,
    edge_color: Color,
) {
    if let Ok(mesh) =
        TessellateSolid::new(solid, TessellationParams::default()).execute(store)
    {
        register_face(storage, mesh, mesh_color);
    }
    if let Ok(solid_data) = store.solid(solid) {
        register_edges(storage, store, solid_data.outer_shell, edge_color);
    }
}

pub fn register(storage: &MeshStorage) {
    // Case 1: Subtract — window hole in wall
    // Window extends beyond wall in z to avoid coplanar faces
    let bx = -14.0;
    let by = 0.0;
    register_label(storage, bx - 2.0, by + 8.0, "1", LABEL_SIZE, LABEL_COLOR);
    {
        let mut store = TopologyStore::new();
        let wall = make_box(&mut store, bx, by, 0.0, 6.0, 6.0, 4.0);
        let window = make_box(&mut store, bx + 1.5, by + 1.5, -0.5, 3.0, 3.0, 5.0);
        if let (Some(a), Some(b)) = (wall, window) {
            if let Ok(result) = Subtract::new(a, b).execute(&mut store) {
                let edge_color = Color::rgb(60, 60, 60);
                render_solid(storage, &store, result, GREEN, edge_color);
            }
        }
    }

    // Case 2: Union — two overlapping boxes
    let bx = -4.0;
    let by = 0.0;
    register_label(storage, bx - 2.0, by + 8.0, "2", LABEL_SIZE, LABEL_COLOR);
    {
        let mut store = TopologyStore::new();
        let a = make_box(&mut store, bx, by, 0.0, 4.0, 4.0, 3.0);
        let b = make_box(&mut store, bx + 2.0, by + 2.0, 1.0, 4.0, 4.0, 3.0);
        if let (Some(a), Some(b)) = (a, b) {
            if let Ok(result) = Union::new(a, b).execute(&mut store) {
                let edge_color = Color::rgb(60, 60, 60);
                render_solid(storage, &store, result, BLUE, edge_color);
            }
        }
    }

    // Case 3: Intersect — overlap region
    let bx = 8.0;
    let by = 0.0;
    register_label(storage, bx - 2.0, by + 8.0, "3", LABEL_SIZE, LABEL_COLOR);
    {
        let mut store = TopologyStore::new();
        let a = make_box(&mut store, bx, by, 0.0, 4.0, 4.0, 3.0);
        let b = make_box(&mut store, bx + 2.0, by + 2.0, 1.0, 4.0, 4.0, 3.0);
        if let (Some(a), Some(b)) = (a, b) {
            if let Ok(result) = Intersect::new(a, b).execute(&mut store) {
                let edge_color = Color::rgb(60, 60, 60);
                render_solid(storage, &store, result, RED, edge_color);
            }
        }
    }
}
