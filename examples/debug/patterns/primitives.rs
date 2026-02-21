use geolis::math::{Point3, Vector3};
use geolis::operations::creation::{MakeCone, MakeCylinder, MakeSphere};
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

#[allow(clippy::too_many_lines)]
pub fn register(storage: &MeshStorage) {
    let spacing = 12.0;
    let edge_color = Color::rgb(60, 60, 60);

    // Case 1: MakeCylinder (r=3, h=6)
    {
        let bx = 0.0;
        let by = 0.0;
        register_label(storage, bx - 2.0, by + 8.0, "1", LABEL_SIZE, LABEL_COLOR);

        let mut topo = TopologyStore::new();
        if let Ok(solid) =
            MakeCylinder::new(Point3::new(bx, by, 0.0), 3.0, Vector3::z(), 6.0).execute(&mut topo)
        {
            render_solid(storage, &topo, solid, GREEN, edge_color);
        }
    }

    // Case 2: MakeSphere (r=3)
    {
        let bx = spacing;
        let by = 0.0;
        register_label(storage, bx - 2.0, by + 8.0, "2", LABEL_SIZE, LABEL_COLOR);

        let mut topo = TopologyStore::new();
        if let Ok(solid) = MakeSphere::new(Point3::new(bx, by, 0.0), 3.0).execute(&mut topo) {
            render_solid(storage, &topo, solid, BLUE, edge_color);
        }
    }

    // Case 3: MakeCone full (r=3, h=6, pointed)
    {
        let bx = spacing * 2.0;
        let by = 0.0;
        register_label(storage, bx - 2.0, by + 8.0, "3", LABEL_SIZE, LABEL_COLOR);

        let mut topo = TopologyStore::new();
        if let Ok(solid) =
            MakeCone::new(Point3::new(bx, by, 0.0), 3.0, 0.0, Vector3::z(), 6.0)
                .execute(&mut topo)
        {
            render_solid(storage, &topo, solid, GREEN, edge_color);
        }
    }

    // Case 4: MakeCone truncated (r1=3, r2=1.5, h=6)
    {
        let bx = spacing * 3.0;
        let by = 0.0;
        register_label(storage, bx - 2.0, by + 8.0, "4", LABEL_SIZE, LABEL_COLOR);

        let mut topo = TopologyStore::new();
        if let Ok(solid) =
            MakeCone::new(Point3::new(bx, by, 0.0), 3.0, 1.5, Vector3::z(), 6.0)
                .execute(&mut topo)
        {
            render_solid(storage, &topo, solid, BLUE, edge_color);
        }
    }
}
