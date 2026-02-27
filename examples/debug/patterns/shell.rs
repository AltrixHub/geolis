use geolis::math::Point3;
use geolis::operations::creation::MakeBox;
use geolis::operations::modification::Shell;
use geolis::tessellation::{TessellateSolid, TessellationParams};
use geolis::topology::{FaceSurface, TopologyStore};
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_edges, register_face, register_label};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const GREEN: Color = Color::rgb(100, 200, 100);
const BLUE: Color = Color::rgb(100, 150, 255);
const RED: Color = Color::rgb(230, 100, 100);

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

/// Gets the face with the highest component in the given axis direction.
fn get_face_by_normal(
    store: &TopologyStore,
    solid: geolis::topology::SolidId,
    axis: fn(&geolis::math::Vector3) -> f64,
) -> Option<geolis::topology::FaceId> {
    let solid_data = store.solid(solid).ok()?;
    let shell = store.shell(solid_data.outer_shell).ok()?;
    let mut best = *shell.faces.first()?;
    let mut best_val = f64::NEG_INFINITY;
    for &face_id in &shell.faces {
        let Ok(face) = store.face(face_id) else {
            continue;
        };
        if let FaceSurface::Plane(plane) = &face.surface {
            let n = if face.same_sense {
                *plane.plane_normal()
            } else {
                -*plane.plane_normal()
            };
            let val = axis(&n);
            if val > best_val {
                best_val = val;
                best = face_id;
            }
        }
    }
    Some(best)
}

pub fn register(storage: &MeshStorage) {
    let spacing = 14.0;
    let edge_color = Color::rgb(60, 60, 60);

    // Case 1: Shell with top removed (thickness=0.5)
    {
        let bx = 0.0;
        let by = 0.0;
        register_label(storage, bx - 2.0, by + 8.0, "1", LABEL_SIZE, LABEL_COLOR);

        let mut topo = TopologyStore::new();
        if let Ok(solid) =
            MakeBox::new(Point3::new(bx, by, 0.0), Point3::new(bx + 4.0, by + 4.0, 4.0))
                .execute(&mut topo)
        {
            if let Some(top) = get_face_by_normal(&topo, solid, |n| n.z) {
                if let Ok(result) = Shell::new(solid, 0.5, vec![top]).execute(&mut topo) {
                    render_solid(storage, &topo, result, GREEN, edge_color);
                }
            }
        }
    }

    // Case 2: Shell with side (+x) removed (thickness=0.5)
    {
        let bx = spacing;
        let by = 0.0;
        register_label(storage, bx - 2.0, by + 8.0, "2", LABEL_SIZE, LABEL_COLOR);

        let mut topo = TopologyStore::new();
        if let Ok(solid) =
            MakeBox::new(Point3::new(bx, by, 0.0), Point3::new(bx + 4.0, by + 4.0, 4.0))
                .execute(&mut topo)
        {
            if let Some(side) = get_face_by_normal(&topo, solid, |n| n.x) {
                if let Ok(result) = Shell::new(solid, 0.5, vec![side]).execute(&mut topo) {
                    render_solid(storage, &topo, result, BLUE, edge_color);
                }
            }
        }
    }

    // Case 3: Shell with thick walls (thickness=1.5, 6x6x6 box, top removed)
    {
        let bx = spacing * 2.0;
        let by = 0.0;
        register_label(storage, bx - 2.0, by + 8.0, "3", LABEL_SIZE, LABEL_COLOR);

        let mut topo = TopologyStore::new();
        if let Ok(solid) =
            MakeBox::new(Point3::new(bx, by, 0.0), Point3::new(bx + 6.0, by + 6.0, 6.0))
                .execute(&mut topo)
        {
            if let Some(top) = get_face_by_normal(&topo, solid, |n| n.z) {
                if let Ok(result) = Shell::new(solid, 1.5, vec![top]).execute(&mut topo) {
                    render_solid(storage, &topo, result, RED, edge_color);
                }
            }
        }
    }
}
