use std::collections::{HashMap, HashSet, VecDeque};
use std::f64::consts::TAU;

use spade::handles::FixedFaceHandle;
use spade::{
    ConstrainedDelaunayTriangulation, InsertionError, Point2 as SpadePoint2, Triangulation,
};

use crate::error::{Result, TessellationError};
use crate::geometry::surface::Surface;
use crate::math::{Point2, Vector3};
use crate::topology::{EdgeCurve, FaceId, FaceSurface, TopologyStore, WireId};

use super::{TessellationMode, TessellationParams, TriangleMesh};

/// Tessellates a face into a triangle mesh.
pub struct TessellateFace {
    face: FaceId,
    params: TessellationParams,
}

impl TessellateFace {
    /// Creates a new `TessellateFace` operation.
    #[must_use]
    pub fn new(face: FaceId, params: TessellationParams) -> Self {
        Self { face, params }
    }

    /// Executes the tessellation, returning a triangle mesh.
    ///
    /// # Errors
    ///
    /// Returns an error if the face cannot be tessellated.
    #[allow(clippy::cast_possible_truncation)]
    pub fn execute(&self, store: &TopologyStore) -> Result<TriangleMesh> {
        let face = store.face(self.face)?;
        let same_sense = face.same_sense;
        let outer_wire_id = face.outer_wire;

        let full_rev = wire_has_full_circle(store, outer_wire_id);

        match &face.surface {
            FaceSurface::Plane(plane) => {
                let plane = plane.clone();
                if full_rev {
                    // Annular disc (or full disc) from revolve — use polar grid
                    let (r_min, r_max, center) = extract_annular_radii(store, outer_wire_id)?;
                    tessellate_annular_disc(&plane, &center, r_min, r_max, same_sense, &self.params)
                } else {
                    let inner_wire_ids = face.inner_wires.clone();
                    tessellate_plane(store, &plane, same_sense, outer_wire_id, &inner_wire_ids)
                }
            }
            FaceSurface::Cylinder(cyl) => {
                let cyl = cyl.clone();
                let outer_3d = collect_wire_points_tessellated(store, outer_wire_id, &self.params)?;
                let (_, _, v_min, v_max) = compute_uv_bounds(&outer_3d, |p| cyl.inverse(p));
                let (u_min, u_max) = if full_rev {
                    (0.0, TAU)
                } else {
                    compute_unwrapped_u_bounds(&outer_3d, |p| cyl.inverse(p))
                };
                let n_u = adaptive_angular_segments(cyl.radius(), u_max - u_min, &self.params);
                let n_v = adaptive_linear_segments(v_max - v_min, &self.params);
                tessellate_surface(&cyl, u_min, u_max, v_min, v_max, n_u, n_v, same_sense, &self.params)
            }
            FaceSurface::Sphere(sph) => {
                let sph = sph.clone();
                let outer_3d = collect_wire_points_tessellated(store, outer_wire_id, &self.params)?;
                let (_, _, v_min, v_max) = compute_uv_bounds(&outer_3d, |p| sph.inverse(p));
                let (u_min, u_max) = if full_rev {
                    (0.0, TAU)
                } else {
                    compute_unwrapped_u_bounds(&outer_3d, |p| sph.inverse(p))
                };
                let n_u = adaptive_angular_segments(sph.radius(), u_max - u_min, &self.params);
                let n_v = adaptive_angular_segments(sph.radius(), v_max - v_min, &self.params);
                tessellate_surface(&sph, u_min, u_max, v_min, v_max, n_u, n_v, same_sense, &self.params)
            }
            FaceSurface::Cone(cone) => {
                let cone = cone.clone();
                let outer_3d = collect_wire_points_tessellated(store, outer_wire_id, &self.params)?;
                let (_, _, v_min, v_max) = compute_uv_bounds(&outer_3d, |p| cone.inverse(p));
                let (u_min, u_max) = if full_rev {
                    (0.0, TAU)
                } else {
                    compute_unwrapped_u_bounds(&outer_3d, |p| cone.inverse(p))
                };
                let max_radius = v_max * cone.half_angle().sin();
                let n_u = adaptive_angular_segments(max_radius, u_max - u_min, &self.params);
                let n_v = adaptive_linear_segments(v_max - v_min, &self.params);
                tessellate_surface(&cone, u_min, u_max, v_min, v_max, n_u, n_v, same_sense, &self.params)
            }
            FaceSurface::Torus(torus) => {
                let torus = torus.clone();
                let outer_3d = collect_wire_points_tessellated(store, outer_wire_id, &self.params)?;
                let (_, _, v_min, v_max) = compute_uv_bounds(&outer_3d, |p| torus.inverse(p));
                let (u_min, u_max) = if full_rev {
                    (0.0, TAU)
                } else {
                    compute_unwrapped_u_bounds(&outer_3d, |p| torus.inverse(p))
                };
                let n_u = adaptive_angular_segments(
                    torus.major_radius() + torus.minor_radius(),
                    u_max - u_min,
                    &self.params,
                );
                let n_v = adaptive_angular_segments(torus.minor_radius(), v_max - v_min, &self.params);
                tessellate_surface(&torus, u_min, u_max, v_min, v_max, n_u, n_v, same_sense, &self.params)
            }
        }
    }
}

/// Tessellates a planar face using CDT.
#[allow(clippy::cast_possible_truncation)]
fn tessellate_plane(
    store: &TopologyStore,
    plane: &crate::geometry::surface::Plane,
    same_sense: bool,
    outer_wire_id: crate::topology::WireId,
    inner_wire_ids: &[crate::topology::WireId],
) -> Result<TriangleMesh> {
    let params = TessellationParams::default();
    let outer_3d = collect_wire_points_tessellated(store, outer_wire_id, &params)?;
    let mut inner_3d_list = Vec::new();
    for &wire_id in inner_wire_ids {
        inner_3d_list.push(collect_wire_points_tessellated(store, wire_id, &params)?);
    }

    let origin = plane.origin();
    let u_dir = plane.u_dir();
    let v_dir = plane.v_dir();
    let normal = if same_sense {
        *plane.plane_normal()
    } else {
        -*plane.plane_normal()
    };

    let project = |p: &crate::math::Point3| -> SpadePoint2<f64> {
        let d = p - origin;
        SpadePoint2::new(d.dot(u_dir), d.dot(v_dir))
    };

    let outer_2d: Vec<_> = outer_3d.iter().map(&project).collect();
    let inner_2d_list: Vec<Vec<_>> = inner_3d_list
        .iter()
        .map(|pts| pts.iter().map(&project).collect())
        .collect();

    let mut cdt = ConstrainedDelaunayTriangulation::<SpadePoint2<f64>>::new();
    insert_constraint_loop(&mut cdt, &outer_2d)?;
    for inner_2d in &inner_2d_list {
        insert_constraint_loop(&mut cdt, inner_2d)?;
    }

    let interior_faces = classify_interior_faces(&cdt);

    let mut mesh = TriangleMesh::default();
    let mut vertex_map: HashMap<usize, u32> = HashMap::new();

    for face_handle in cdt.inner_faces() {
        let fix = face_handle.fix();
        if !interior_faces.contains(&fix.index()) {
            continue;
        }

        let verts = face_handle.vertices();
        let mut tri_indices = [0u32; 3];

        for (i, vh) in verts.iter().enumerate() {
            let idx = vh.fix().index();
            let mesh_idx = if let Some(&existing) = vertex_map.get(&idx) {
                existing
            } else {
                let pos = vh.position();
                let u = pos.x;
                let v = pos.y;
                let p3 = *origin + *u_dir * u + *v_dir * v;
                let new_idx = mesh.vertices.len() as u32;
                mesh.vertices.push(p3);
                mesh.normals.push(normal);
                mesh.uvs.push(Point2::new(u, v));
                vertex_map.insert(idx, new_idx);
                new_idx
            };
            tri_indices[i] = mesh_idx;
        }

        mesh.indices.push(tri_indices);
    }

    Ok(mesh)
}

/// Checks if a wire contains a full-circle edge (sweep ≈ TAU).
fn wire_has_full_circle(store: &TopologyStore, wire_id: WireId) -> bool {
    let Ok(wire) = store.wire(wire_id) else {
        return false;
    };
    for oe in &wire.edges {
        let Ok(edge) = store.edge(oe.edge) else {
            continue;
        };
        if matches!(&edge.curve, EdgeCurve::Circle(_)) {
            let sweep = (edge.t_end - edge.t_start).abs();
            if sweep > TAU - 0.01 {
                return true;
            }
        }
    }
    false
}

/// Computes u-bounds by unwrapping `atan2` values along the wire boundary.
///
/// The surface's `inverse()` returns `u` via `atan2`, which has a discontinuity
/// at ±π. This function tracks cumulative angular deltas to produce a continuous
/// `(u_min, u_max)` range that correctly handles sweeps beyond π (e.g., 270°).
///
/// This works regardless of whether the surface's angular direction matches the
/// Arc edge's direction (e.g., Cone with reversed axis vs Cylinder with aligned axis).
fn compute_unwrapped_u_bounds(
    points: &[crate::math::Point3],
    inverse: impl Fn(&crate::math::Point3) -> (f64, f64),
) -> (f64, f64) {
    if points.is_empty() {
        return (0.0, 0.0);
    }

    let (first_u, _) = inverse(&points[0]);
    let mut u_min = first_u;
    let mut u_max = first_u;
    let mut prev_raw = first_u;
    let mut running = first_u;

    for p in &points[1..] {
        let (raw_u, _) = inverse(p);
        let mut delta = raw_u - prev_raw;
        // Unwrap: keep delta in (-π, π]
        if delta > std::f64::consts::PI {
            delta -= TAU;
        } else if delta < -std::f64::consts::PI {
            delta += TAU;
        }
        running += delta;
        u_min = u_min.min(running);
        u_max = u_max.max(running);
        prev_raw = raw_u;
    }

    (u_min, u_max)
}

/// Extracts the min/max radii and center from circle edges in a wire.
///
/// Used for annular disc tessellation. If only one circle is found,
/// `r_min` is 0 (full disc).
fn extract_annular_radii(
    store: &TopologyStore,
    wire_id: WireId,
) -> Result<(f64, f64, crate::math::Point3)> {
    let wire = store.wire(wire_id)?;
    let mut radii = Vec::new();
    let mut center = None;

    for oe in &wire.edges {
        let edge = store.edge(oe.edge)?;
        if let EdgeCurve::Circle(circle) = &edge.curve {
            radii.push(circle.radius());
            if center.is_none() {
                center = Some(*circle.center());
            }
        }
    }

    let center = center
        .ok_or_else(|| TessellationError::Failed("no circle edges in annular disc wire".into()))?;

    let r_max = radii.iter().copied().fold(0.0_f64, f64::max);
    let r_min = if radii.len() >= 2 {
        radii.iter().copied().fold(f64::INFINITY, f64::min)
    } else {
        0.0
    };

    Ok((r_min, r_max, center))
}

/// Tessellates an annular disc (or full disc) on a plane using a polar grid.
///
/// Instead of CDT (which struggles with slit-annulus constraint polygons from
/// full-circle edges), this generates a regular grid in polar coordinates
/// `(θ, r)` and evaluates points directly on the plane.
#[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss, clippy::unnecessary_wraps)]
fn tessellate_annular_disc(
    plane: &crate::geometry::surface::Plane,
    center: &crate::math::Point3,
    r_min: f64,
    r_max: f64,
    same_sense: bool,
    params: &TessellationParams,
) -> Result<TriangleMesh> {
    let n_theta = adaptive_angular_segments(r_max, TAU, params);
    let n_r = adaptive_linear_segments(r_max - r_min, params).max(1);

    let normal = if same_sense {
        *plane.plane_normal()
    } else {
        -*plane.plane_normal()
    };
    let u_dir = plane.u_dir();
    let v_dir = plane.v_dir();

    let mut mesh = TriangleMesh::default();
    // n_theta columns (last column wraps to first — no +1)
    let cols = n_theta;
    let rows = n_r + 1;

    mesh.vertices.reserve(rows * cols);
    mesh.normals.reserve(rows * cols);
    mesh.uvs.reserve(rows * cols);
    mesh.indices.reserve(n_theta * n_r * 2);

    // Generate vertices in polar grid
    for ir in 0..rows {
        let r = r_min + (r_max - r_min) * ir as f64 / n_r as f64;
        for itheta in 0..cols {
            let theta = TAU * itheta as f64 / n_theta as f64;
            let pt = *center + *u_dir * (r * theta.cos()) + *v_dir * (r * theta.sin());
            mesh.vertices.push(pt);
            mesh.normals.push(normal);
            mesh.uvs.push(Point2::new(theta, r));
        }
    }

    // Generate triangles — wrap around in θ direction
    for ir in 0..n_r {
        for itheta in 0..n_theta {
            let next_theta = (itheta + 1) % n_theta;
            let i00 = (ir * cols + itheta) as u32;
            let i10 = (ir * cols + next_theta) as u32;
            let i01 = ((ir + 1) * cols + itheta) as u32;
            let i11 = ((ir + 1) * cols + next_theta) as u32;
            if same_sense {
                mesh.indices.push([i00, i10, i11]);
                mesh.indices.push([i00, i11, i01]);
            } else {
                mesh.indices.push([i00, i11, i10]);
                mesh.indices.push([i00, i01, i11]);
            }
        }
    }

    Ok(mesh)
}

/// Tessellates a parametric surface on a UV grid.
///
/// Generates `(n_u + 1) * (n_v + 1)` vertices via `surface.evaluate(u, v)`,
/// then splits each quad cell into two triangles.
#[allow(clippy::cast_possible_truncation)]
#[allow(clippy::too_many_arguments)]
fn tessellate_uv_grid(
    surface: &dyn Surface,
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
    n_u: usize,
    n_v: usize,
    same_sense: bool,
) -> Result<TriangleMesh> {
    let mut mesh = TriangleMesh::default();
    let rows = n_v + 1;
    let cols = n_u + 1;
    mesh.vertices.reserve(rows * cols);
    mesh.normals.reserve(rows * cols);
    mesh.uvs.reserve(rows * cols);
    mesh.indices.reserve(n_u * n_v * 2);

    // Generate vertices
    for iv in 0..rows {
        #[allow(clippy::cast_precision_loss)]
        let v = v_min + (v_max - v_min) * iv as f64 / n_v as f64;
        for iu in 0..cols {
            #[allow(clippy::cast_precision_loss)]
            let u = u_min + (u_max - u_min) * iu as f64 / n_u as f64;
            let pt = surface.evaluate(u, v)?;
            let n = surface.normal(u, v).unwrap_or(Vector3::z());
            let n = if same_sense { n } else { -n };
            mesh.vertices.push(pt);
            mesh.normals.push(n);
            mesh.uvs.push(Point2::new(u, v));
        }
    }

    // Generate triangles (two per quad cell)
    for iv in 0..n_v {
        for iu in 0..n_u {
            let i00 = (iv * cols + iu) as u32;
            let i10 = (iv * cols + iu + 1) as u32;
            let i01 = ((iv + 1) * cols + iu) as u32;
            let i11 = ((iv + 1) * cols + iu + 1) as u32;
            if same_sense {
                mesh.indices.push([i00, i10, i11]);
                mesh.indices.push([i00, i11, i01]);
            } else {
                mesh.indices.push([i00, i11, i10]);
                mesh.indices.push([i00, i01, i11]);
            }
        }
    }

    Ok(mesh)
}

/// Dispatches to either `tessellate_uv_grid` or `tessellate_uv_adaptive` based on mode.
///
/// In adaptive mode, the curvature-computed `n_u`/`n_v` are ignored; instead a
/// coarse base grid (`min_segments × min_segments`) is used, and cells are
/// recursively subdivided where the midpoint deviation exceeds the tolerance.
#[allow(clippy::too_many_arguments)]
fn tessellate_surface(
    surface: &dyn Surface,
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
    n_u: usize,
    n_v: usize,
    same_sense: bool,
    params: &TessellationParams,
) -> Result<TriangleMesh> {
    match params.mode {
        TessellationMode::Default => {
            tessellate_uv_grid(surface, u_min, u_max, v_min, v_max, n_u, n_v, same_sense)
        }
        TessellationMode::Adaptive => {
            let base = params.min_segments;
            tessellate_uv_adaptive(surface, u_min, u_max, v_min, v_max, base, base, same_sense, params.tolerance)
        }
    }
}

/// Maximum recursion depth for adaptive subdivision.
///
/// Each level doubles the effective resolution per base cell, so depth 6
/// yields up to 64× the base resolution in each direction.
const MAX_ADAPTIVE_DEPTH: usize = 6;

/// Tessellates a parametric surface using adaptive midpoint subdivision.
///
/// Starts with a `base_n_u × base_n_v` grid and recursively subdivides cells
/// whose midpoint deviation from bilinear interpolation exceeds `tolerance`.
#[allow(clippy::too_many_arguments, clippy::similar_names)]
fn tessellate_uv_adaptive(
    surface: &dyn Surface,
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
    base_n_u: usize,
    base_n_v: usize,
    same_sense: bool,
    tolerance: f64,
) -> Result<TriangleMesh> {
    let mut mesh = TriangleMesh::default();
    let mut vertex_cache: HashMap<(u64, u64), u32> = HashMap::new();

    #[allow(clippy::cast_precision_loss)]
    let du = (u_max - u_min) / base_n_u as f64;
    #[allow(clippy::cast_precision_loss)]
    let dv = (v_max - v_min) / base_n_v as f64;

    for iv in 0..base_n_v {
        #[allow(clippy::cast_precision_loss)]
        let cv0 = v_min + dv * iv as f64;
        #[allow(clippy::cast_precision_loss)]
        let cv1 = v_min + dv * (iv + 1) as f64;
        for iu in 0..base_n_u {
            #[allow(clippy::cast_precision_loss)]
            let cu0 = u_min + du * iu as f64;
            #[allow(clippy::cast_precision_loss)]
            let cu1 = u_min + du * (iu + 1) as f64;
            subdivide_cell(
                surface, cu0, cu1, cv0, cv1, same_sense, tolerance, 0,
                &mut mesh, &mut vertex_cache,
            )?;
        }
    }

    Ok(mesh)
}

/// Recursively subdivides a UV cell if its midpoint deviation exceeds the tolerance.
///
/// If the surface midpoint deviates from the bilinear interpolation of the 4 corners
/// by more than `tolerance`, the cell is split into 4 sub-cells. Otherwise, 2 triangles
/// are emitted for the cell.
#[allow(clippy::too_many_arguments)]
fn subdivide_cell(
    surface: &dyn Surface,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    same_sense: bool,
    tolerance: f64,
    depth: usize,
    mesh: &mut TriangleMesh,
    cache: &mut HashMap<(u64, u64), u32>,
) -> Result<()> {
    let mid_u = f64::midpoint(u0, u1);
    let mid_v = f64::midpoint(v0, v1);

    let p00 = surface.evaluate(u0, v0)?;
    let p10 = surface.evaluate(u1, v0)?;
    let p01 = surface.evaluate(u0, v1)?;
    let p11 = surface.evaluate(u1, v1)?;
    let actual_mid = surface.evaluate(mid_u, mid_v)?;

    let bilinear_mid = crate::math::Point3::new(
        (p00.x + p10.x + p01.x + p11.x) / 4.0,
        (p00.y + p10.y + p01.y + p11.y) / 4.0,
        (p00.z + p10.z + p01.z + p11.z) / 4.0,
    );

    let deviation = (actual_mid - bilinear_mid).norm();

    if deviation > tolerance && depth < MAX_ADAPTIVE_DEPTH {
        subdivide_cell(surface, u0, mid_u, v0, mid_v, same_sense, tolerance, depth + 1, mesh, cache)?;
        subdivide_cell(surface, mid_u, u1, v0, mid_v, same_sense, tolerance, depth + 1, mesh, cache)?;
        subdivide_cell(surface, u0, mid_u, mid_v, v1, same_sense, tolerance, depth + 1, mesh, cache)?;
        subdivide_cell(surface, mid_u, u1, mid_v, v1, same_sense, tolerance, depth + 1, mesh, cache)?;
    } else {
        let i00 = get_or_insert_vertex(mesh, cache, surface, u0, v0, same_sense)?;
        let i10 = get_or_insert_vertex(mesh, cache, surface, u1, v0, same_sense)?;
        let i01 = get_or_insert_vertex(mesh, cache, surface, u0, v1, same_sense)?;
        let i11 = get_or_insert_vertex(mesh, cache, surface, u1, v1, same_sense)?;

        if same_sense {
            mesh.indices.push([i00, i10, i11]);
            mesh.indices.push([i00, i11, i01]);
        } else {
            mesh.indices.push([i00, i11, i10]);
            mesh.indices.push([i00, i01, i11]);
        }
    }

    Ok(())
}

/// Gets an existing vertex from the cache or inserts a new one.
///
/// Uses `f64::to_bits()` for exact UV deduplication so adjacent cells
/// at the same subdivision level share vertices perfectly.
#[allow(clippy::cast_possible_truncation)]
fn get_or_insert_vertex(
    mesh: &mut TriangleMesh,
    cache: &mut HashMap<(u64, u64), u32>,
    surface: &dyn Surface,
    u: f64,
    v: f64,
    same_sense: bool,
) -> Result<u32> {
    let key = (u.to_bits(), v.to_bits());
    if let Some(&idx) = cache.get(&key) {
        return Ok(idx);
    }
    let pt = surface.evaluate(u, v)?;
    let n = surface.normal(u, v).unwrap_or(Vector3::z());
    let n = if same_sense { n } else { -n };
    let idx = mesh.vertices.len() as u32;
    mesh.vertices.push(pt);
    mesh.normals.push(n);
    mesh.uvs.push(Point2::new(u, v));
    cache.insert(key, idx);
    Ok(idx)
}

/// Computes UV bounds from wire vertex positions using an inverse parametrization function.
fn compute_uv_bounds(
    points: &[crate::math::Point3],
    inverse: impl Fn(&crate::math::Point3) -> (f64, f64),
) -> (f64, f64, f64, f64) {
    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;
    let mut v_min = f64::INFINITY;
    let mut v_max = f64::NEG_INFINITY;

    for pt in points {
        let (u, v) = inverse(pt);
        u_min = u_min.min(u);
        u_max = u_max.max(u);
        v_min = v_min.min(v);
        v_max = v_max.max(v);
    }

    (u_min, u_max, v_min, v_max)
}

/// Computes the number of segments for an angular (u) parameter range based on chord error.
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn adaptive_angular_segments(radius: f64, sweep: f64, params: &TessellationParams) -> usize {
    if radius > params.tolerance {
        let half_angle = (1.0 - params.tolerance / radius).acos();
        let computed = (sweep / (2.0 * half_angle)).ceil() as usize;
        computed.clamp(params.min_segments, params.max_segments)
    } else {
        params.min_segments
    }
}

/// Computes the number of segments for a linear (v) parameter range.
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn adaptive_linear_segments(extent: f64, params: &TessellationParams) -> usize {
    let computed = (extent / params.tolerance).ceil() as usize;
    computed.clamp(params.min_segments, params.max_segments)
}

/// Collects 3D points from a wire, tessellating curved edges into polylines.
///
/// For Line edges, only the start point is included (avoiding duplicates).
/// For Circle/Arc/Ellipse edges, intermediate points are sampled along the curve.
fn collect_wire_points_tessellated(
    store: &TopologyStore,
    wire_id: crate::topology::WireId,
    params: &TessellationParams,
) -> Result<Vec<crate::math::Point3>> {
    use crate::geometry::curve::Curve;
    use crate::topology::EdgeCurve;

    let edges = store.wire(wire_id)?.edges.clone();
    let mut points = Vec::new();

    for oe in &edges {
        let edge = store.edge(oe.edge)?;
        let (t_start, t_end) = if oe.forward {
            (edge.t_start, edge.t_end)
        } else {
            (edge.t_end, edge.t_start)
        };

        match &edge.curve {
            EdgeCurve::Line(line) => {
                points.push(line.evaluate(t_start)?);
            }
            EdgeCurve::Arc(arc) => {
                let n = tessellate_edge_segments(arc.radius(), t_start, t_end, params);
                add_curve_samples(&mut points, arc, t_start, t_end, n)?;
            }
            EdgeCurve::Circle(circle) => {
                let n = tessellate_edge_segments(circle.radius(), t_start, t_end, params);
                add_curve_samples(&mut points, circle, t_start, t_end, n)?;
            }
            EdgeCurve::Ellipse(ellipse) => {
                let n = tessellate_edge_segments(ellipse.semi_major(), t_start, t_end, params);
                add_curve_samples(&mut points, ellipse, t_start, t_end, n)?;
            }
        }
    }

    Ok(points)
}

/// Computes the number of segments for a curved edge.
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn tessellate_edge_segments(
    radius: f64,
    t_start: f64,
    t_end: f64,
    params: &TessellationParams,
) -> usize {
    let sweep = (t_end - t_start).abs();
    if radius > params.tolerance {
        let half_angle = (1.0 - params.tolerance / radius).acos();
        let computed = (sweep / (2.0 * half_angle)).ceil() as usize;
        computed.clamp(params.min_segments, params.max_segments)
    } else {
        params.min_segments
    }
}

/// Adds sample points from a curve (excluding the last point to avoid duplicates).
fn add_curve_samples(
    points: &mut Vec<crate::math::Point3>,
    curve: &dyn crate::geometry::curve::Curve,
    t_start: f64,
    t_end: f64,
    n: usize,
) -> Result<()> {
    for i in 0..n {
        #[allow(clippy::cast_precision_loss)]
        let frac = i as f64 / n as f64;
        let t = t_start + frac * (t_end - t_start);
        points.push(curve.evaluate(t)?);
    }
    Ok(())
}

/// Inserts a closed polygon as constraint edges into the CDT.
fn insert_constraint_loop(
    cdt: &mut ConstrainedDelaunayTriangulation<SpadePoint2<f64>>,
    points: &[SpadePoint2<f64>],
) -> Result<()> {
    if points.len() < 3 {
        return Err(
            TessellationError::Failed("constraint loop needs at least 3 points".into()).into(),
        );
    }

    let mut handles = Vec::with_capacity(points.len());
    for &pt in points {
        let h = cdt
            .insert(pt)
            .map_err(|e: InsertionError| TessellationError::Failed(format!("CDT insert: {e}")))?;
        handles.push(h);
    }

    for i in 0..handles.len() {
        let from = handles[i];
        let to = handles[(i + 1) % handles.len()];
        if from != to {
            cdt.add_constraint(from, to);
        }
    }

    Ok(())
}

/// Classifies which inner faces of the CDT are inside the polygon using flood-fill.
///
/// Starts from faces adjacent to the outer (infinite) face at depth 0. Each time
/// a constraint edge is crossed, depth increments. Odd depth = interior.
fn classify_interior_faces(
    cdt: &ConstrainedDelaunayTriangulation<SpadePoint2<f64>>,
) -> HashSet<usize> {
    let mut interior = HashSet::new();
    let mut depth_map: HashMap<usize, u32> = HashMap::new();
    let mut queue: VecDeque<(FixedFaceHandle<spade::handles::InnerTag>, u32)> = VecDeque::new();

    let outer_fix = cdt.outer_face().fix();

    // Seed: find inner faces adjacent to the outer face via directed edges
    for edge in cdt.directed_edges() {
        if edge.face().fix() == outer_fix {
            let rev_face = edge.rev().face();
            if let Some(inner) = rev_face.as_inner() {
                let idx = inner.fix().index();
                if depth_map.contains_key(&idx) {
                    continue;
                }
                let depth = u32::from(cdt.is_constraint_edge(edge.as_undirected().fix()));
                depth_map.insert(idx, depth);
                if depth % 2 == 1 {
                    interior.insert(idx);
                }
                queue.push_back((inner.fix(), depth));
            }
        }
    }

    // BFS flood-fill
    while let Some((face_fix, depth)) = queue.pop_front() {
        let face = cdt.face(face_fix);
        for edge in face.adjacent_edges() {
            let neighbor = edge.rev().face();
            if let Some(inner_neighbor) = neighbor.as_inner() {
                let n_idx = inner_neighbor.fix().index();
                if depth_map.contains_key(&n_idx) {
                    continue;
                }
                let new_depth = if cdt.is_constraint_edge(edge.as_undirected().fix()) {
                    depth + 1
                } else {
                    depth
                };
                depth_map.insert(n_idx, new_depth);
                if new_depth % 2 == 1 {
                    interior.insert(n_idx);
                }
                queue.push_back((inner_neighbor.fix(), new_depth));
            }
        }
    }

    interior
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::math::Point3;
    use crate::operations::creation::{MakeFace, MakeWire};

    fn p(x: f64, y: f64) -> Point3 {
        Point3::new(x, y, 0.0)
    }

    fn make_face_from_points(
        store: &mut crate::topology::TopologyStore,
        points: Vec<Point3>,
    ) -> FaceId {
        let wire = MakeWire::new(points, true).execute(store).unwrap();
        MakeFace::new(wire, vec![]).execute(store).unwrap()
    }

    #[test]
    fn triangle_produces_1_triangle() {
        let mut store = crate::topology::TopologyStore::new();
        let face = make_face_from_points(&mut store, vec![p(0.0, 0.0), p(4.0, 0.0), p(2.0, 3.0)]);
        let mesh = TessellateFace::new(face, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert_eq!(mesh.indices.len(), 1);
        assert_eq!(mesh.vertices.len(), 3);
        assert_eq!(mesh.normals.len(), 3);
        assert_eq!(mesh.uvs.len(), 3);
    }

    #[test]
    fn square_produces_2_triangles() {
        let mut store = crate::topology::TopologyStore::new();
        let face = make_face_from_points(
            &mut store,
            vec![p(0.0, 0.0), p(4.0, 0.0), p(4.0, 4.0), p(0.0, 4.0)],
        );
        let mesh = TessellateFace::new(face, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert_eq!(mesh.indices.len(), 2);
        assert_eq!(mesh.vertices.len(), 4);
    }

    #[test]
    fn l_shape_concave_tessellates() {
        let mut store = crate::topology::TopologyStore::new();
        let face = make_face_from_points(
            &mut store,
            vec![
                p(0.0, 0.0),
                p(4.0, 0.0),
                p(4.0, 2.0),
                p(2.0, 2.0),
                p(2.0, 4.0),
                p(0.0, 4.0),
            ],
        );
        let mesh = TessellateFace::new(face, TessellationParams::default())
            .execute(&store)
            .unwrap();
        // L-shape (6 vertices, concave) → should produce 4 triangles
        assert_eq!(mesh.indices.len(), 4);
        assert_eq!(mesh.vertices.len(), 6);
    }

    #[test]
    fn face_with_hole_excludes_interior() {
        let mut store = crate::topology::TopologyStore::new();
        let outer = MakeWire::new(
            vec![p(0.0, 0.0), p(10.0, 0.0), p(10.0, 10.0), p(0.0, 10.0)],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let inner = MakeWire::new(
            vec![p(3.0, 3.0), p(7.0, 3.0), p(7.0, 7.0), p(3.0, 7.0)],
            true,
        )
        .execute(&mut store)
        .unwrap();
        let face = MakeFace::new(outer, vec![inner])
            .execute(&mut store)
            .unwrap();
        let mesh = TessellateFace::new(face, TessellationParams::default())
            .execute(&store)
            .unwrap();

        // No triangle center should be inside the hole (3..7, 3..7)
        for tri in &mesh.indices {
            let cx = (mesh.vertices[tri[0] as usize].x
                + mesh.vertices[tri[1] as usize].x
                + mesh.vertices[tri[2] as usize].x)
                / 3.0;
            let cy = (mesh.vertices[tri[0] as usize].y
                + mesh.vertices[tri[1] as usize].y
                + mesh.vertices[tri[2] as usize].y)
                / 3.0;
            let in_hole = cx > 3.0 && cx < 7.0 && cy > 3.0 && cy < 7.0;
            assert!(
                !in_hole,
                "triangle centroid ({cx}, {cy}) is inside the hole"
            );
        }
    }

    #[test]
    fn normals_match_plane_normal() {
        let mut store = crate::topology::TopologyStore::new();
        let face = make_face_from_points(
            &mut store,
            vec![p(0.0, 0.0), p(4.0, 0.0), p(4.0, 4.0), p(0.0, 4.0)],
        );
        let mesh = TessellateFace::new(face, TessellationParams::default())
            .execute(&store)
            .unwrap();
        for n in &mesh.normals {
            assert!(n.z.abs() > 0.99, "normal z should be ~1.0, got {}", n.z);
        }
    }

    // ── Curved surface tessellation tests ──────────────────────

    use crate::geometry::curve::Circle;
    use crate::geometry::surface::{Cylinder, Sphere, Torus};
    use crate::math::Vector3;
    use crate::topology::{
        EdgeCurve, EdgeData, FaceData, OrientedEdge, VertexData, WireData,
    };
    use std::f64::consts::TAU;

    /// Helper: creates a cylindrical face with a wire of 4 vertices
    /// spanning the cylinder from u=0..TAU, v=0..height.
    fn make_cylinder_face(store: &mut crate::topology::TopologyStore, radius: f64, height: f64) -> FaceId {
        let cyl = Cylinder::new(
            Point3::origin(), radius, Vector3::z(), Vector3::x(),
        ).unwrap();

        // 4 corner vertices of the "unrolled" cylinder patch
        let v0 = store.add_vertex(VertexData::new(Point3::new(radius, 0.0, 0.0)));
        let v1 = store.add_vertex(VertexData::new(Point3::new(radius, 0.0, height)));

        // Circle edges at bottom and top (full circles, so start=end vertex)
        let bottom_circle = Circle::new(
            Point3::new(0.0, 0.0, 0.0), radius, Vector3::z(), Vector3::x(),
        ).unwrap();
        let top_circle = Circle::new(
            Point3::new(0.0, 0.0, height), radius, Vector3::z(), Vector3::x(),
        ).unwrap();

        let e_bottom = store.add_edge(EdgeData {
            start: v0, end: v0,
            curve: EdgeCurve::Circle(bottom_circle),
            t_start: 0.0, t_end: TAU,
        });
        let e_top = store.add_edge(EdgeData {
            start: v1, end: v1,
            curve: EdgeCurve::Circle(top_circle),
            t_start: 0.0, t_end: TAU,
        });
        let e_seam = store.add_edge(EdgeData {
            start: v0, end: v1,
            curve: EdgeCurve::Line(
                crate::geometry::curve::Line::new(
                    Point3::new(radius, 0.0, 0.0),
                    Vector3::new(0.0, 0.0, height),
                ).unwrap()
            ),
            t_start: 0.0, t_end: height,
        });

        let wire = store.add_wire(WireData {
            edges: vec![
                OrientedEdge::new(e_bottom, true),
                OrientedEdge::new(e_seam, true),
                OrientedEdge::new(e_top, false),
                OrientedEdge::new(e_seam, false),
            ],
            is_closed: true,
        });

        store.add_face(FaceData {
            surface: FaceSurface::Cylinder(cyl),
            outer_wire: wire,
            inner_wires: vec![],
            same_sense: true,
        })
    }

    #[test]
    fn cylinder_face_tessellates() {
        let mut store = crate::topology::TopologyStore::new();
        let face = make_cylinder_face(&mut store, 2.0, 5.0);
        let mesh = TessellateFace::new(face, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty());
        assert_eq!(mesh.vertices.len(), mesh.normals.len());
        assert_eq!(mesh.vertices.len(), mesh.uvs.len());
    }

    #[test]
    fn cylinder_normals_point_outward() {
        let mut store = crate::topology::TopologyStore::new();
        let face = make_cylinder_face(&mut store, 2.0, 5.0);
        let mesh = TessellateFace::new(face, TessellationParams::default())
            .execute(&store)
            .unwrap();
        for (i, n) in mesh.normals.iter().enumerate() {
            let v = &mesh.vertices[i];
            // For z-axis cylinder, normal should point radially outward in XY
            let radial = Vector3::new(v.x, v.y, 0.0);
            let radial_len = radial.norm();
            if radial_len > 1e-6 {
                let dot = n.dot(&(radial / radial_len));
                assert!(dot > 0.9, "normal at {v:?} not outward: dot={dot}");
            }
        }
    }

    /// Helper: creates a full sphere face (u=0..TAU, v=-PI/2..PI/2).
    fn make_sphere_face(store: &mut crate::topology::TopologyStore, radius: f64) -> FaceId {
        let sph = Sphere::new(
            Point3::origin(), radius, Vector3::z(), Vector3::x(),
        ).unwrap();

        // For a full sphere, we need a wire. Use two poles + two meridian seams.
        let south = Point3::new(0.0, 0.0, -radius);
        let north = Point3::new(0.0, 0.0, radius);
        let v_south = store.add_vertex(VertexData::new(south));
        let v_north = store.add_vertex(VertexData::new(north));

        // Seam edge along u=0 from south to north
        let seam_dir = north - south;
        let e_seam = store.add_edge(EdgeData {
            start: v_south, end: v_north,
            curve: EdgeCurve::Line(
                crate::geometry::curve::Line::new(south, seam_dir).unwrap()
            ),
            t_start: 0.0, t_end: seam_dir.norm(),
        });

        let wire = store.add_wire(WireData {
            edges: vec![
                OrientedEdge::new(e_seam, true),
                OrientedEdge::new(e_seam, false),
            ],
            is_closed: true,
        });

        store.add_face(FaceData {
            surface: FaceSurface::Sphere(sph),
            outer_wire: wire,
            inner_wires: vec![],
            same_sense: true,
        })
    }

    #[test]
    fn sphere_face_tessellates() {
        let mut store = crate::topology::TopologyStore::new();
        let face = make_sphere_face(&mut store, 3.0);
        let mesh = TessellateFace::new(face, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty());
        assert_eq!(mesh.vertices.len(), mesh.normals.len());
    }

    #[test]
    fn sphere_normals_point_outward() {
        let mut store = crate::topology::TopologyStore::new();
        let face = make_sphere_face(&mut store, 3.0);
        let mesh = TessellateFace::new(face, TessellationParams::default())
            .execute(&store)
            .unwrap();
        for (i, n) in mesh.normals.iter().enumerate() {
            let v = &mesh.vertices[i];
            let len = Vector3::new(v.x, v.y, v.z).norm();
            if len > 1e-6 {
                let expected = Vector3::new(v.x, v.y, v.z) / len;
                let dot = n.dot(&expected);
                assert!(dot > 0.9, "normal at {v:?} not outward: dot={dot}");
            }
        }
    }

    /// Helper: creates a full torus face.
    fn make_torus_face(
        store: &mut crate::topology::TopologyStore,
        major_r: f64,
        minor_r: f64,
    ) -> FaceId {
        let torus = Torus::new(
            Point3::origin(), major_r, minor_r, Vector3::z(), Vector3::x(),
        ).unwrap();

        // Outer point at u=0, v=0
        let pt = Point3::new(major_r + minor_r, 0.0, 0.0);
        let v0 = store.add_vertex(VertexData::new(pt));

        // Single seam edge (degenerate closed loop for full torus)
        let e_seam = store.add_edge(EdgeData {
            start: v0, end: v0,
            curve: EdgeCurve::Circle(
                Circle::new(
                    Point3::new(major_r, 0.0, 0.0),
                    minor_r,
                    Vector3::y(),
                    Vector3::x(),
                ).unwrap()
            ),
            t_start: 0.0, t_end: TAU,
        });

        let wire = store.add_wire(WireData {
            edges: vec![OrientedEdge::new(e_seam, true)],
            is_closed: true,
        });

        store.add_face(FaceData {
            surface: FaceSurface::Torus(torus),
            outer_wire: wire,
            inner_wires: vec![],
            same_sense: true,
        })
    }

    #[test]
    fn torus_face_tessellates() {
        let mut store = crate::topology::TopologyStore::new();
        let face = make_torus_face(&mut store, 3.0, 1.0);
        let mesh = TessellateFace::new(face, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty());
        assert_eq!(mesh.vertices.len(), mesh.normals.len());
        assert_eq!(mesh.vertices.len(), mesh.uvs.len());
    }

    // ── Adaptive tessellation tests ──────────────────────────────

    use super::TessellationMode;

    #[test]
    fn adaptive_sphere_tessellates() {
        let mut store = crate::topology::TopologyStore::new();
        let face = make_sphere_face(&mut store, 3.0);
        let params = TessellationParams {
            mode: TessellationMode::Adaptive,
            ..TessellationParams::default()
        };
        let mesh = TessellateFace::new(face, params)
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty());
        assert_eq!(mesh.vertices.len(), mesh.normals.len());
        assert_eq!(mesh.vertices.len(), mesh.uvs.len());
    }

    #[test]
    fn adaptive_cylinder_tessellates() {
        let mut store = crate::topology::TopologyStore::new();
        let face = make_cylinder_face(&mut store, 2.0, 5.0);
        let params = TessellationParams {
            mode: TessellationMode::Adaptive,
            ..TessellationParams::default()
        };
        let mesh = TessellateFace::new(face, params)
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty());
        assert_eq!(mesh.vertices.len(), mesh.normals.len());
    }

    #[test]
    fn adaptive_torus_tessellates() {
        let mut store = crate::topology::TopologyStore::new();
        let face = make_torus_face(&mut store, 3.0, 1.0);
        let params = TessellationParams {
            mode: TessellationMode::Adaptive,
            ..TessellationParams::default()
        };
        let mesh = TessellateFace::new(face, params)
            .execute(&store)
            .unwrap();
        assert!(!mesh.indices.is_empty());
        assert_eq!(mesh.vertices.len(), mesh.normals.len());
    }

    #[test]
    fn adaptive_cylinder_subdivides_coarse_cells() {
        // With a coarse tolerance, the default grid has large cells.
        // Adaptive should subdivide cells whose midpoint deviation exceeds
        // tolerance, producing more triangles than the base grid.
        let mut store = crate::topology::TopologyStore::new();
        let face = make_cylinder_face(&mut store, 2.0, 5.0);

        let coarse = TessellationParams {
            tolerance: 0.5,
            min_segments: 4,
            max_segments: 256,
            mode: TessellationMode::Default,
        };
        let default_mesh = TessellateFace::new(face, coarse)
            .execute(&store)
            .unwrap();

        let adaptive = TessellationParams {
            tolerance: 0.5,
            min_segments: 4,
            max_segments: 256,
            mode: TessellationMode::Adaptive,
        };
        let adaptive_mesh = TessellateFace::new(face, adaptive)
            .execute(&store)
            .unwrap();

        // With coarse tolerance on a cylinder, the initial 4-segment grid
        // has cells with significant midpoint deviation, so adaptive should
        // subdivide and produce more triangles.
        assert!(
            adaptive_mesh.indices.len() > default_mesh.indices.len(),
            "adaptive ({}) should produce more triangles than default ({}) at coarse tolerance",
            adaptive_mesh.indices.len(),
            default_mesh.indices.len(),
        );
    }

    #[test]
    fn adaptive_sphere_normals_outward() {
        let mut store = crate::topology::TopologyStore::new();
        let face = make_sphere_face(&mut store, 3.0);
        let params = TessellationParams {
            mode: TessellationMode::Adaptive,
            ..TessellationParams::default()
        };
        let mesh = TessellateFace::new(face, params)
            .execute(&store)
            .unwrap();
        for (i, n) in mesh.normals.iter().enumerate() {
            let v = &mesh.vertices[i];
            let len = Vector3::new(v.x, v.y, v.z).norm();
            if len > 1e-6 {
                let expected = Vector3::new(v.x, v.y, v.z) / len;
                let dot = n.dot(&expected);
                assert!(dot > 0.9, "normal at {v:?} not outward: dot={dot}");
            }
        }
    }

    #[test]
    fn default_mode_unchanged() {
        // Verify that TessellationMode::Default produces same result as before
        let mut store = crate::topology::TopologyStore::new();
        let face = make_cylinder_face(&mut store, 2.0, 5.0);
        let params_default = TessellationParams {
            mode: TessellationMode::Default,
            ..TessellationParams::default()
        };
        let mesh = TessellateFace::new(face, params_default)
            .execute(&store)
            .unwrap();

        // Should produce the same mesh as with no mode specified (same Default)
        let mesh2 = TessellateFace::new(face, TessellationParams::default())
            .execute(&store)
            .unwrap();
        assert_eq!(mesh.indices.len(), mesh2.indices.len());
        assert_eq!(mesh.vertices.len(), mesh2.vertices.len());
    }
}
