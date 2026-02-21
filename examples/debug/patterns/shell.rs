use geolis::math::Point3;
use geolis::operations::creation::MakeBox;
use geolis::operations::modification::Shell;
use geolis::tessellation::{TessellateFace, TessellateSolid, TessellationParams};
use geolis::topology::{FaceSurface, TopologyStore};
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_edges, register_face, register_label};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const GREEN: Color = Color::rgb(100, 200, 100);
const BLUE: Color = Color::rgb(100, 150, 255);
const RED: Color = Color::rgb(230, 100, 80);

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

/// Renders each face of a shell's solid individually with alternating colors
/// to make the structure visible (outer, inner, side faces).
fn render_shell_faces(
    storage: &MeshStorage,
    topo: &TopologyStore,
    solid: geolis::topology::SolidId,
    outer_color: Color,
    inner_color: Color,
    side_color: Color,
    n_kept: usize,
    edge_color: Color,
) {
    let params = TessellationParams::default();
    let Ok(solid_data) = topo.solid(solid) else {
        return;
    };
    let Ok(shell) = topo.shell(solid_data.outer_shell) else {
        return;
    };

    for (i, &face_id) in shell.faces.iter().enumerate() {
        let color = if i < n_kept * 2 {
            // First n_kept*2 faces are outer/inner pairs
            if i % 2 == 0 { outer_color } else { inner_color }
        } else {
            side_color
        };
        if let Ok(mesh) = TessellateFace::new(face_id, params).execute(topo) {
            register_face(storage, mesh, color);
        }
    }
    register_edges(storage, topo, solid_data.outer_shell, edge_color);
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

    // Case 1: Shell with top removed — rendered as full solid
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
                if let Ok(result) = Shell::new(solid, 1.5, vec![top]).execute(&mut topo) {
                    render_solid(storage, &topo, result, GREEN, edge_color);
                }
            }
        }
    }

    // Case 2: Shell with top removed — faces colored by type
    // Outer=green, Inner=red, Side=blue
    {
        let bx = spacing;
        let by = 0.0;
        register_label(storage, bx - 2.0, by + 8.0, "2", LABEL_SIZE, LABEL_COLOR);

        let mut topo = TopologyStore::new();
        if let Ok(solid) =
            MakeBox::new(Point3::new(bx, by, 0.0), Point3::new(bx + 4.0, by + 4.0, 4.0))
                .execute(&mut topo)
        {
            if let Some(top) = get_face_by_normal(&topo, solid, |n| n.z) {
                // 5 kept faces (6 - 1 removed)
                if let Ok(result) = Shell::new(solid, 1.5, vec![top]).execute(&mut topo) {
                    render_shell_faces(
                        storage, &topo, result,
                        GREEN, RED, BLUE,
                        5, edge_color,
                    );
                }
            }
        }
    }

    // Case 3: Shell with side removed — faces colored by type
    {
        let bx = spacing * 2.0;
        let by = 0.0;
        register_label(storage, bx - 2.0, by + 8.0, "3", LABEL_SIZE, LABEL_COLOR);

        let mut topo = TopologyStore::new();
        if let Ok(solid) =
            MakeBox::new(Point3::new(bx, by, 0.0), Point3::new(bx + 4.0, by + 4.0, 4.0))
                .execute(&mut topo)
        {
            if let Some(side) = get_face_by_normal(&topo, solid, |n| n.x) {
                if let Ok(result) = Shell::new(solid, 1.5, vec![side]).execute(&mut topo) {
                    render_shell_faces(
                        storage, &topo, result,
                        GREEN, RED, BLUE,
                        5, edge_color,
                    );
                }
            }
        }
    }
}
