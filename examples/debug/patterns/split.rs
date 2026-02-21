use geolis::math::{Point3, Vector3};
use geolis::operations::creation::MakeBox;
use geolis::operations::modification::Split;
use geolis::tessellation::{TessellateSolid, TessellationParams};
use geolis::topology::TopologyStore;
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_edges, register_face, register_label};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const GREEN: Color = Color::rgb(100, 200, 100);
const BLUE: Color = Color::rgb(100, 150, 255);

fn render_solid(
    storage: &MeshStorage,
    topo: &TopologyStore,
    solid: geolis::topology::SolidId,
    mesh_color: Color,
    edge_color: Color,
) {
    if let Ok(mesh) = TessellateSolid::new(solid, TessellationParams::default()).execute(topo) {
        register_face(storage, mesh, mesh_color);
    }
    if let Ok(solid_data) = topo.solid(solid) {
        register_edges(storage, topo, solid_data.outer_shell, edge_color);
    }
}

pub fn register(storage: &MeshStorage) {
    let spacing = 14.0;
    let edge_color = Color::rgb(60, 60, 60);

    // Case 1: Box split horizontally (z=2)
    {
        let bx = 0.0;
        let by = 0.0;
        register_label(storage, bx - 2.0, by + 8.0, "1", LABEL_SIZE, LABEL_COLOR);

        let mut topo = TopologyStore::new();
        if let Ok(solid) =
            MakeBox::new(Point3::new(bx, by, 0.0), Point3::new(bx + 4.0, by + 4.0, 4.0))
                .execute(&mut topo)
        {
            if let Ok((top, bottom)) =
                Split::new(solid, Point3::new(0.0, 0.0, 2.0), Vector3::z()).execute(&mut topo)
            {
                render_solid(storage, &topo, top, GREEN, edge_color);
                render_solid(storage, &topo, bottom, BLUE, edge_color);
            }
        }
    }

    // Case 2: Box split vertically (x=2)
    {
        let bx = spacing;
        let by = 0.0;
        register_label(storage, bx - 2.0, by + 8.0, "2", LABEL_SIZE, LABEL_COLOR);

        let mut topo = TopologyStore::new();
        if let Ok(solid) =
            MakeBox::new(Point3::new(bx, by, 0.0), Point3::new(bx + 4.0, by + 4.0, 4.0))
                .execute(&mut topo)
        {
            if let Ok((left, right)) =
                Split::new(solid, Point3::new(bx + 2.0, 0.0, 0.0), Vector3::x()).execute(&mut topo)
            {
                render_solid(storage, &topo, left, GREEN, edge_color);
                render_solid(storage, &topo, right, BLUE, edge_color);
            }
        }
    }

    // Case 3: Box split diagonally (45 deg in XZ plane)
    {
        let bx = spacing * 2.0;
        let by = 0.0;
        register_label(storage, bx - 2.0, by + 8.0, "3", LABEL_SIZE, LABEL_COLOR);

        let mut topo = TopologyStore::new();
        if let Ok(solid) =
            MakeBox::new(Point3::new(bx, by, 0.0), Point3::new(bx + 4.0, by + 4.0, 4.0))
                .execute(&mut topo)
        {
            let normal = Vector3::new(1.0, 0.0, 1.0).normalize();
            if let Ok((a, b)) =
                Split::new(solid, Point3::new(bx + 2.0, by + 2.0, 2.0), normal).execute(&mut topo)
            {
                render_solid(storage, &topo, a, GREEN, edge_color);
                render_solid(storage, &topo, b, BLUE, edge_color);
            }
        }
    }
}
