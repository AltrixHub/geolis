use geolis::math::{Point3, Vector3};
use geolis::operations::creation::{MakeFace, MakeWire};
use geolis::operations::shaping::Extrude;
use geolis::tessellation::{StrokeStyle, TessellateSolid, TessellationParams};
use geolis::topology::TopologyStore;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_edges, register_face, register_label, register_stroke};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const GRAY: Color = Color::rgb(150, 150, 150);
const GREEN: Color = Color::rgb(100, 200, 100);
const BLUE: Color = Color::rgb(100, 150, 255);
const RED: Color = Color::rgb(255, 100, 100);
const EDGE_COLOR: Color = Color::rgb(60, 60, 60);

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

fn stroke(storage: &MeshStorage, points: &[Point3], closed: bool, color: Color) {
    if let Ok(style) = StrokeStyle::new(0.05) {
        register_stroke(storage, points, style, closed, color);
    }
}

pub fn register(storage: &MeshStorage) {
    // Case 1: Subtract — wall(6x6x4) - window(3x3x5)
    // Show input boxes as 3D solids
    let bx = -14.0;
    let by = 0.0;
    register_label(storage, bx - 2.0, by + 8.0, "1", LABEL_SIZE, LABEL_COLOR);
    {
        let mut store = TopologyStore::new();
        // Wall (gray solid)
        if let Some(wall) = make_box(&mut store, bx, by, 0.0, 6.0, 6.0, 4.0) {
            render_solid(storage, &store, wall, GRAY, EDGE_COLOR);
        }
        // Window opening (green wireframe outline at z=0 and z=4)
        let hole_bottom = [
            Point3::new(bx + 1.5, by + 1.5, 0.0),
            Point3::new(bx + 4.5, by + 1.5, 0.0),
            Point3::new(bx + 4.5, by + 4.5, 0.0),
            Point3::new(bx + 1.5, by + 4.5, 0.0),
        ];
        let hole_top = [
            Point3::new(bx + 1.5, by + 1.5, 4.0),
            Point3::new(bx + 4.5, by + 1.5, 4.0),
            Point3::new(bx + 4.5, by + 4.5, 4.0),
            Point3::new(bx + 1.5, by + 4.5, 4.0),
        ];
        stroke(storage, &hole_bottom, true, GREEN);
        stroke(storage, &hole_top, true, GREEN);
        // Vertical edges of the hole
        for i in 0..4 {
            stroke(storage, &[hole_bottom[i], hole_top[i]], false, GREEN);
        }
    }

    // Case 2: Union — two overlapping boxes
    // Box A: (bx..bx+4, by..by+4, 0..3)
    // Box B: (bx+2..bx+6, by+2..by+6, 1..4)
    let bx = -4.0;
    let by = 0.0;
    register_label(storage, bx - 2.0, by + 8.0, "2", LABEL_SIZE, LABEL_COLOR);
    {
        let mut store = TopologyStore::new();
        if let Some(a) = make_box(&mut store, bx, by, 0.0, 4.0, 4.0, 3.0) {
            render_solid(storage, &store, a, GRAY, EDGE_COLOR);
        }
        if let Some(b) = make_box(&mut store, bx + 2.0, by + 2.0, 1.0, 4.0, 4.0, 3.0) {
            render_solid(storage, &store, b, BLUE, EDGE_COLOR);
        }
    }

    // Case 3: Intersect — overlap region is a 2x2x2 box
    // Box A: (bx..bx+4, by..by+4, 0..3)
    // Box B: (bx+2..bx+6, by+2..by+6, 1..4)
    // Expected: (bx+2..bx+4, by+2..by+4, 1..3)
    let bx = 8.0;
    let by = 0.0;
    register_label(storage, bx - 2.0, by + 8.0, "3", LABEL_SIZE, LABEL_COLOR);
    {
        let mut store = TopologyStore::new();
        // Input boxes as gray wireframe
        let a_bottom = [
            Point3::new(bx, by, 0.0),
            Point3::new(bx + 4.0, by, 0.0),
            Point3::new(bx + 4.0, by + 4.0, 0.0),
            Point3::new(bx, by + 4.0, 0.0),
        ];
        let a_top = [
            Point3::new(bx, by, 3.0),
            Point3::new(bx + 4.0, by, 3.0),
            Point3::new(bx + 4.0, by + 4.0, 3.0),
            Point3::new(bx, by + 4.0, 3.0),
        ];
        stroke(storage, &a_bottom, true, GRAY);
        stroke(storage, &a_top, true, GRAY);
        for i in 0..4 {
            stroke(storage, &[a_bottom[i], a_top[i]], false, GRAY);
        }

        let b_bottom = [
            Point3::new(bx + 2.0, by + 2.0, 1.0),
            Point3::new(bx + 6.0, by + 2.0, 1.0),
            Point3::new(bx + 6.0, by + 6.0, 1.0),
            Point3::new(bx + 2.0, by + 6.0, 1.0),
        ];
        let b_top = [
            Point3::new(bx + 2.0, by + 2.0, 4.0),
            Point3::new(bx + 6.0, by + 2.0, 4.0),
            Point3::new(bx + 6.0, by + 6.0, 4.0),
            Point3::new(bx + 2.0, by + 6.0, 4.0),
        ];
        stroke(storage, &b_bottom, true, GRAY);
        stroke(storage, &b_top, true, GRAY);
        for i in 0..4 {
            stroke(storage, &[b_bottom[i], b_top[i]], false, GRAY);
        }

        // Expected intersection: 2x2x2 box at (bx+2, by+2, 1) to (bx+4, by+4, 3)
        if let Some(expected) = make_box(&mut store, bx + 2.0, by + 2.0, 1.0, 2.0, 2.0, 2.0) {
            render_solid(storage, &store, expected, RED, EDGE_COLOR);
        }
    }
}
