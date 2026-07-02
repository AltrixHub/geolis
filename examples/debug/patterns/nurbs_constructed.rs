//! NURBS surface-constructor showcase (P2 constructors).
//!
//! Renders four adaptively tessellated surfaces, one per constructor:
//! 1. extrude  — a circle swept linearly into a cylinder shell,
//! 2. revolve  — a profile line revolved into a cone,
//! 3. loft     — six sinusoidally sized circle sections skinned into a
//!    bulging/necking vase silhouette,
//! 4. sweep    — a circular profile carried along a helical coil rail
//!    (radius 2, height 6, 2.5 turns).

use geolis::geometry::nurbs::{NurbsCurve3D, NurbsSurface};
use geolis::math::{Point3, Vector3};
use geolis::tessellation::{tessellate_nurbs_surface, SurfaceTessellationOptions, TriangleMesh};
use revion_ui::value_objects::Color;
use revion_ui::MeshStorage;

use super::{register_face, register_label, SceneBounds};

const LABEL_SIZE: f64 = 1.2;
const LABEL_COLOR: Color = Color::rgb(255, 220, 80);
const GREEN: Color = Color::rgb(120, 200, 140);
const BLUE: Color = Color::rgb(120, 170, 230);
const ORANGE: Color = Color::rgb(230, 160, 90);
const PURPLE: Color = Color::rgb(180, 140, 220);

fn register_surface(
    storage: &MeshStorage,
    bounds: &mut SceneBounds,
    surface: &NurbsSurface,
    bx: f64,
    by: f64,
    color: Color,
) {
    let options = SurfaceTessellationOptions::default();
    let Ok(mut mesh) = tessellate_nurbs_surface(surface, &options) else {
        return;
    };
    translate(&mut mesh, bx, by);
    register_face(storage, bounds, mesh, color);
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

/// Wavy vase: six exact NURBS circle sections stacked along +Z, with radii
/// modulated by a sine so the silhouette bulges and necks.
fn lofted_surface() -> Option<NurbsSurface> {
    let mut sections = Vec::with_capacity(6);
    for k in 0..6 {
        let z = 1.2 * f64::from(k);
        let radius = 1.6 + 0.6 * (f64::from(k) * 1.4).sin();
        let circle =
            NurbsCurve3D::circle(Point3::new(0.0, 0.0, z), radius, Vector3::z(), Vector3::x())
                .ok()?;
        sections.push(circle);
    }
    // Circles share degree/knots, so the loft interpolates them exactly.
    NurbsSurface::loft(&sections, None).ok()
}

/// Helical coil: a small circular profile swept along a helix rail.
///
/// The rail is interpolated through helix samples `(2cosθ, 2sinθ, h·θ/θ_max)`
/// for θ ∈ [0, 2.5·2π]. The sweep rigidly transports the profile from the rail
/// start with rotation-minimizing frames, so the profile is centred at the rail
/// start point in the plane perpendicular to the rail start tangent.
fn swept_surface() -> Option<NurbsSurface> {
    const HELIX_RADIUS: f64 = 2.0;
    const HELIX_HEIGHT: f64 = 6.0;
    const PROFILE_RADIUS: f64 = 0.5;
    const TURNS: f64 = 2.5;
    let theta_max = TURNS * std::f64::consts::TAU;

    // ~24 samples along the helix for the interpolated rail.
    let samples = 24_u32;
    let mut points = Vec::with_capacity(samples as usize);
    for i in 0..samples {
        let theta = theta_max * f64::from(i) / f64::from(samples - 1);
        points.push(Point3::new(
            HELIX_RADIUS * theta.cos(),
            HELIX_RADIUS * theta.sin(),
            HELIX_HEIGHT * theta / theta_max,
        ));
    }
    let (rail, _) = NurbsCurve3D::interpolate(&points, 3).ok()?;

    // Centre the profile at the rail start, plane perpendicular to the start
    // tangent (sweep transports control points relative to the rail start).
    let (t0, _) = rail.parameter_domain();
    let ders = rail.derivatives(t0, 1).ok()?;
    let start = Point3::from(ders[0]);
    let tangent = ders[1].normalize();
    // Any direction perpendicular to the tangent works as the circle's ref dir.
    let seed = if tangent.dot(&Vector3::z()).abs() < 0.9 {
        Vector3::z()
    } else {
        Vector3::x()
    };
    let perp = (seed - tangent * tangent.dot(&seed)).normalize();
    let profile = NurbsCurve3D::circle(start, PROFILE_RADIUS, tangent, perp).ok()?;

    NurbsSurface::sweep(&profile, &rail).ok()
}

pub fn register(storage: &MeshStorage, bounds: &mut SceneBounds) {
    let cases: [(f64, f64, &str, Color, Option<NurbsSurface>); 4] = [
        (0.0, 0.0, "1", GREEN, extruded_cylinder()),
        (10.0, 0.0, "2", BLUE, revolved_cone()),
        (20.0, 0.0, "3", ORANGE, lofted_surface()),
        (30.0, 0.0, "4", PURPLE, swept_surface()),
    ];
    for (bx, by, label, color, surface) in cases {
        register_label(storage, bx - 1.5, by + 8.0, label, LABEL_SIZE, LABEL_COLOR);
        if let Some(surface) = surface {
            register_surface(storage, bounds, &surface, bx, by, color);
        }
    }
}
