//! NURBS surface tessellation showcase.
//!
//! Renders adaptively tessellated NURBS surfaces as shaded meshes:
//! 1. a free-form bicubic bump patch, and
//! 2. a rational quarter-cylinder shell.

use geolis::geometry::nurbs::{KnotVector, NurbsSurface};
use geolis::math::Point3;
use geolis::tessellation::{tessellate_nurbs_surface, SurfaceTessellationOptions, TriangleMesh};
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_face, register_label};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const GREEN: Color = Color::rgb(120, 200, 140);
const BLUE: Color = Color::rgb(120, 170, 230);

/// Tessellate a surface, translate it by `(bx, by)`, and register the shaded mesh.
fn register_surface(storage: &MeshStorage, surface: &NurbsSurface, bx: f64, by: f64, color: Color) {
    let options = SurfaceTessellationOptions::default();
    let Ok(mut mesh) = tessellate_nurbs_surface(surface, &options) else {
        return;
    };
    translate(&mut mesh, bx, by);
    register_face(storage, mesh, color);
}

/// Offset every vertex in the XY plane (normals/UVs are unaffected).
fn translate(mesh: &mut TriangleMesh, bx: f64, by: f64) {
    for v in &mut mesh.vertices {
        v.x += bx;
        v.y += by;
    }
}

/// Free-form bicubic patch with a central bump (4x4 control grid).
fn bump_patch() -> Option<NurbsSurface> {
    let heights = [
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 2.5, 2.5, 0.0],
        [0.0, 2.5, 2.5, 0.0],
        [0.0, 0.0, 0.0, 0.0],
    ];
    let mut control = Vec::with_capacity(16);
    for (i, row) in heights.iter().enumerate() {
        for (j, &h) in row.iter().enumerate() {
            #[allow(clippy::cast_precision_loss)]
            control.push(Point3::new(i as f64 * 2.0, j as f64 * 2.0, h));
        }
    }
    let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]).ok()?;
    NurbsSurface::from_unweighted(control, 4, 4, knots.clone(), knots, 3, 3).ok()
}

/// Rational quarter-cylinder shell (radius 3, height 6 along +Z).
fn quarter_cylinder() -> Option<NurbsSurface> {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let r = 3.0;
    let h = 6.0;
    NurbsSurface::new(
        vec![
            Point3::new(r, 0.0, 0.0),
            Point3::new(r, 0.0, h),
            Point3::new(r, r, 0.0),
            Point3::new(r, r, h),
            Point3::new(0.0, r, 0.0),
            Point3::new(0.0, r, h),
        ],
        vec![1.0, 1.0, w, w, 1.0, 1.0],
        3,
        2,
        KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]).ok()?,
        KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).ok()?,
        2,
        1,
    )
    .ok()
}

pub fn register(storage: &MeshStorage) {
    // Case 1: free-form bump patch.
    {
        let bx = 0.0;
        let by = 0.0;
        register_label(storage, bx - 1.5, by + 8.0, "1", LABEL_SIZE, LABEL_COLOR);
        if let Some(surface) = bump_patch() {
            register_surface(storage, &surface, bx, by, GREEN);
        }
    }

    // Case 2: rational quarter-cylinder shell.
    {
        let bx = 10.0;
        let by = 0.0;
        register_label(storage, bx - 1.5, by + 8.0, "2", LABEL_SIZE, LABEL_COLOR);
        if let Some(surface) = quarter_cylinder() {
            register_surface(storage, &surface, bx, by, BLUE);
        }
    }
}
