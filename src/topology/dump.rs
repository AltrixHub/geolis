//! Diagnostic helpers — render topology entities as full-precision strings
//! suitable for logging and for copy-pasting into regression tests.
//!
//! These are debug-only utilities: callers route the returned string through
//! their own logger (e.g. `tracing::info!`). The format intentionally uses
//! `{:.17e}` so a coordinate round-trips back to the exact `f64` it came from
//! — that is what lets us pin a failing runtime placement as a `Point3::new`
//! literal in a unit test.

use crate::math::Point3;
use crate::topology::{FaceId, SolidId, TopologyStore};

/// Render every face of `solid_id` as `outer_wire` (and inner wires when
/// present) using full-precision (`{:.17e}`) coordinates.
///
/// # Errors
///
/// Returns `Err(String)` if any topology lookup fails (the solid /
/// shell / face / wire / edge / vertex is missing from `store`).
pub fn dump_solid_full_precision(
    store: &TopologyStore,
    solid_id: SolidId,
) -> Result<String, String> {
    let solid = store
        .solid(solid_id)
        .map_err(|e| format!("dump_solid: solid lookup failed: {e}"))?;
    let shell = store
        .shell(solid.outer_shell)
        .map_err(|e| format!("dump_solid: shell lookup failed: {e}"))?;

    let mut out = String::new();
    out.push_str(&format!(
        "SOLID {solid_id:?} face_count={}\n",
        shell.faces.len()
    ));
    for (face_idx, &face_id) in shell.faces.iter().enumerate() {
        out.push_str(&format!("  face[{face_idx}] {face_id:?}\n"));
        let face_dump = dump_face_full_precision(store, face_id)?;
        for line in face_dump.lines() {
            out.push_str("    ");
            out.push_str(line);
            out.push('\n');
        }
    }
    Ok(out)
}

/// Render a single face as `outer_wire` (and inner wires when present)
/// using full-precision (`{:.17e}`) coordinates. Used both standalone and
/// from `dump_solid_full_precision`.
///
/// # Errors
///
/// Returns `Err(String)` if any topology lookup fails (the face / wire /
/// edge / vertex is missing from `store`).
pub fn dump_face_full_precision(store: &TopologyStore, face_id: FaceId) -> Result<String, String> {
    let face = store
        .face(face_id)
        .map_err(|e| format!("dump_face: face lookup failed: {e}"))?;
    let mut out = String::new();
    out.push_str(&format!("FACE {face_id:?}\n"));
    out.push_str("  outer_wire:\n");
    let outer_pts = collect_wire_points(store, face.outer_wire)?;
    for (i, p) in outer_pts.iter().enumerate() {
        out.push_str(&format!("    [{i}] {}\n", format_point(p)));
    }
    for (hi, &inner_wire) in face.inner_wires.iter().enumerate() {
        out.push_str(&format!("  inner_wire[{hi}]:\n"));
        let inner_pts = collect_wire_points(store, inner_wire)?;
        for (i, p) in inner_pts.iter().enumerate() {
            out.push_str(&format!("    [{i}] {}\n", format_point(p)));
        }
    }
    Ok(out)
}

/// Walk a wire's oriented edges and produce the ordered list of start
/// points (respecting each `OrientedEdge::forward`). The last edge's end
/// point is `outer_pts[0]` for a closed wire so we don't repeat it.
fn collect_wire_points(
    store: &TopologyStore,
    wire_id: crate::topology::WireId,
) -> Result<Vec<Point3>, String> {
    let wire = store
        .wire(wire_id)
        .map_err(|e| format!("dump: wire lookup failed: {e}"))?;
    let mut pts = Vec::with_capacity(wire.edges.len());
    for oe in &wire.edges {
        let edge = store
            .edge(oe.edge)
            .map_err(|e| format!("dump: edge lookup failed: {e}"))?;
        let start_vid = if oe.forward { edge.start } else { edge.end };
        let v = store
            .vertex(start_vid)
            .map_err(|e| format!("dump: vertex lookup failed: {e}"))?;
        pts.push(v.point);
    }
    Ok(pts)
}

fn format_point(p: &Point3) -> String {
    format!("({:.17e}, {:.17e}, {:.17e})", p.x, p.y, p.z)
}
