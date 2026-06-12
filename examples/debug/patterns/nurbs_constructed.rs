//! NURBS surface-constructor showcase (P2 constructors).
//!
//! Renders four adaptively tessellated surfaces, one per constructor:
//! 1. extrude  — a circle swept linearly into a cylinder shell,
//! 2. revolve  — a profile line revolved into a cone,
//! 3. loft     — three stacked section curves blended into a surface,
//! 4. sweep    — a circle profile swept along a curved rail.

use geolis::geometry::nurbs::{NurbsCurve3D, NurbsSurface};
use geolis::math::{Point3, Vector3};
use geolis::tessellation::{tessellate_nurbs_surface, SurfaceTessellationOptions, TriangleMesh};
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_face, register_label};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const GREEN: Color = Color::rgb(120, 200, 140);
const BLUE: Color = Color::rgb(120, 170, 230);
const ORANGE: Color = Color::rgb(230, 160, 90);
const PURPLE: Color = Color::rgb(180, 140, 220);

fn register_surface(storage: &MeshStorage, surface: &NurbsSurface, bx: f64, by: f64, color: Color) {
    let options = SurfaceTessellationOptions::default();
    let Ok(mut mesh) = tessellate_nurbs_surface(surface, &options) else {
        return;
    };
    translate(&mut mesh, bx, by);
    register_face(storage, mesh, color);
}

fn translate(mesh: &mut TriangleMesh, bx: f64, by: f64) {
    for v in &mut mesh.vertices {
        v.x += bx;
        v.y += by;
    }
}

/// Cylinder shell: a unit circle extruded along +Z.
fn extruded_cylinder() -> Option<NurbsSurface> {
    let circle = NurbsCurve3D::circle(
        Point3::new(0.0, 0.0, 0.0),
        2.0,
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(1.0, 0.0, 0.0),
    )
    .ok()?;
    NurbsSurface::extrude(&circle, Vector3::new(0.0, 0.0, 5.0)).ok()
}

/// Cone: a slanted profile line revolved a full turn about the Z axis.
fn revolved_cone() -> Option<NurbsSurface> {
    let profile =
        NurbsCurve3D::polyline(&[Point3::new(0.0, 0.0, 5.0), Point3::new(3.0, 0.0, 0.0)]).ok()?;
    NurbsSurface::revolve(
        &profile,
        Point3::new(0.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
        std::f64::consts::TAU,
    )
    .ok()
}

/// Lofted surface from three interpolated section curves at rising heights.
fn lofted_surface() -> Option<NurbsSurface> {
    let section = |z: f64, bow: f64| -> Option<NurbsCurve3D> {
        let (curve, _) = NurbsCurve3D::interpolate(
            &[
                Point3::new(0.0, 0.0, z),
                Point3::new(2.0, bow, z),
                Point3::new(4.0, 0.0, z),
            ],
            2,
        )
        .ok()?;
        Some(curve)
    };
    let sections = [section(0.0, 0.0)?, section(2.5, 2.0)?, section(5.0, 0.0)?];
    NurbsSurface::loft(&sections, None).ok()
}

/// Swept surface: a circular profile carried along a curved rail.
fn swept_surface() -> Option<NurbsSurface> {
    let profile = NurbsCurve3D::circle(
        Point3::new(0.0, 0.0, 0.0),
        0.8,
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(1.0, 0.0, 0.0),
    )
    .ok()?;
    let (rail, _) = NurbsCurve3D::interpolate(
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 3.0),
            Point3::new(5.0, 0.0, 4.0),
            Point3::new(7.0, 0.0, 1.0),
        ],
        3,
    )
    .ok()?;
    NurbsSurface::sweep(&profile, &rail).ok()
}

pub fn register(storage: &MeshStorage) {
    let cases: [(f64, f64, &str, Color, Option<NurbsSurface>); 4] = [
        (0.0, 0.0, "1", GREEN, extruded_cylinder()),
        (10.0, 0.0, "2", BLUE, revolved_cone()),
        (20.0, 0.0, "3", ORANGE, lofted_surface()),
        (30.0, 0.0, "4", PURPLE, swept_surface()),
    ];
    for (bx, by, label, color, surface) in cases {
        register_label(storage, bx - 1.5, by + 8.0, label, LABEL_SIZE, LABEL_COLOR);
        if let Some(surface) = surface {
            register_surface(storage, &surface, bx, by, color);
        }
    }
}
