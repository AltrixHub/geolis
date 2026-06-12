use crate::error::{GeometryError, Result};
use crate::math::{Point3, Vector3, TOLERANCE};

use crate::geometry::nurbs::{NurbsCurve3D, NurbsSurface};

/// An orthonormal frame `(tangent, normal, binormal)` anchored at `origin`.
#[derive(Debug, Clone, Copy)]
struct Frame {
    origin: Point3,
    tangent: Vector3,
    normal: Vector3,
    binormal: Vector3,
}

impl NurbsSurface {
    /// Sweeps a profile curve along a rail curve (The NURBS Book, §10.4).
    ///
    /// The rail is sampled at uniform parameters and a rotation-minimizing
    /// frame (the double-reflection method of Wang et al.) is propagated along
    /// it. At each station the profile control points are rigidly transported
    /// from the start frame into the station frame, and the transported
    /// profiles are skinned with [`NurbsSurface::loft`].
    ///
    /// The sweep is geometrically approximate: the skin interpolates a finite
    /// set of stations rather than reproducing the exact swept surface.
    ///
    /// # Errors
    ///
    /// Returns an error if the rail tangent vanishes at a sampled station, the
    /// profile cannot be transported, or the loft fails.
    #[allow(clippy::cast_precision_loss)]
    pub fn sweep(profile: &NurbsCurve3D, rail: &NurbsCurve3D) -> Result<Self> {
        // Sample count: enough stations to resolve rail curvature (§10.4).
        let stations = (2 * rail.control_points().len()).max(8);
        let (t0, t1) = rail.parameter_domain();

        // Rail sample points and tangents at uniform parameters.
        let mut points = Vec::with_capacity(stations);
        let mut tangents = Vec::with_capacity(stations);
        for i in 0..stations {
            let t = t0 + (t1 - t0) * (i as f64) / ((stations - 1) as f64);
            let ders = rail.derivatives(t, 1)?;
            let tangent = ders[1];
            let len = tangent.norm();
            if len < TOLERANCE {
                return Err(GeometryError::ZeroVector.into());
            }
            points.push(Point3::from(ders[0]));
            tangents.push(tangent / len);
        }

        let frames = rotation_minimizing_frames(&points, &tangents)?;
        let start = frames[0];

        // Profile control points expressed in the start frame's local
        // coordinates; transported into each station frame.
        let local: Vec<(Vector3, f64)> = profile
            .control_points()
            .iter()
            .zip(profile.weights())
            .map(|(p, &w)| {
                let d = p - start.origin;
                (
                    Vector3::new(
                        d.dot(&start.tangent),
                        d.dot(&start.normal),
                        d.dot(&start.binormal),
                    ),
                    w,
                )
            })
            .collect();

        let mut sections = Vec::with_capacity(stations);
        for frame in &frames {
            let control_points: Vec<Point3> = local
                .iter()
                .map(|(c, _)| {
                    frame.origin + frame.tangent * c.x + frame.normal * c.y + frame.binormal * c.z
                })
                .collect();
            let weights: Vec<f64> = local.iter().map(|(_, w)| *w).collect();
            sections.push(NurbsCurve3D::new(
                control_points,
                weights,
                profile.knots().clone(),
                profile.degree(),
            )?);
        }

        // Skin the transported profiles; cubic in v where the station count
        // allows, matching the loft default.
        NurbsSurface::loft(&sections, None)
    }
}

/// Propagates a rotation-minimizing frame along the sampled rail using the
/// double-reflection method (Wang, Jüttler, Zheng & Liu, 2008).
///
/// # Errors
///
/// Returns an error if an initial frame normal cannot be chosen (vanishing
/// tangent).
fn rotation_minimizing_frames(points: &[Point3], tangents: &[Vector3]) -> Result<Vec<Frame>> {
    let n = points.len();
    let mut frames = Vec::with_capacity(n);

    // Initial normal: any unit vector perpendicular to the first tangent.
    let t0 = tangents[0];
    let normal0 = initial_normal(&t0)?;
    let binormal0 = t0.cross(&normal0);
    frames.push(Frame {
        origin: points[0],
        tangent: t0,
        normal: normal0,
        binormal: binormal0,
    });

    for i in 0..(n - 1) {
        let prev = frames[i];
        let xi1 = points[i + 1];
        let ti1 = tangents[i + 1];

        // First reflection: across the plane bisecting xi and xi+1.
        let v1 = xi1 - prev.origin;
        let c1 = v1.dot(&v1);
        let (r_l, t_l) = if c1 < TOLERANCE {
            (prev.normal, prev.tangent)
        } else {
            let r_l = prev.normal - v1 * (2.0 / c1 * v1.dot(&prev.normal));
            let t_l = prev.tangent - v1 * (2.0 / c1 * v1.dot(&prev.tangent));
            (r_l, t_l)
        };

        // Second reflection: across the plane bisecting t_l and ti+1.
        let v2 = ti1 - t_l;
        let c2 = v2.dot(&v2);
        let normal = if c2 < TOLERANCE {
            r_l
        } else {
            r_l - v2 * (2.0 / c2 * v2.dot(&r_l))
        };
        let binormal = ti1.cross(&normal);

        frames.push(Frame {
            origin: xi1,
            tangent: ti1,
            normal,
            binormal,
        });
    }
    Ok(frames)
}

/// Picks a unit normal perpendicular to `tangent` by reflecting away the
/// most-aligned coordinate axis.
///
/// # Errors
///
/// Returns an error if `tangent` is not unit-length (vanishing tangent).
fn initial_normal(tangent: &Vector3) -> Result<Vector3> {
    if (tangent.norm() - 1.0).abs() > 1e-6 {
        return Err(GeometryError::ZeroVector.into());
    }
    // Choose the axis least aligned with the tangent to avoid a near-parallel
    // cross product.
    let abs = Vector3::new(tangent.x.abs(), tangent.y.abs(), tangent.z.abs());
    let axis = if abs.x <= abs.y && abs.x <= abs.z {
        Vector3::x()
    } else if abs.y <= abs.z {
        Vector3::y()
    } else {
        Vector3::z()
    };
    let normal = (axis - tangent * tangent.dot(&axis)).normalize();
    Ok(normal)
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::cast_precision_loss)]
mod tests {
    use super::*;
    use crate::geometry::curve::Curve;
    use crate::geometry::nurbs::KnotVector;

    fn unit_circle_xy() -> NurbsCurve3D {
        NurbsCurve3D::circle(Point3::origin(), 1.0, Vector3::z(), Vector3::x()).unwrap()
    }

    #[test]
    fn sweep_circle_along_straight_line_matches_extrude() {
        let profile = unit_circle_xy();
        let rail = NurbsCurve3D::from_unweighted(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 4.0)],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
        )
        .unwrap();
        let swept = NurbsSurface::sweep(&profile, &rail).unwrap();

        // The swept tube must lie on the cylinder of radius 1 about the rail
        // (the Z axis), matching a straight extrusion of the circle.
        let ((u0, u1), (v0, v1)) = swept.parameter_domain();
        for i in 0..=12 {
            let u = u0 + (u1 - u0) * f64::from(i) / 12.0;
            for j in 0..=12 {
                let v = v0 + (v1 - v0) * f64::from(j) / 12.0;
                let p = swept.point_at(u, v).unwrap();
                let radial = (p.x * p.x + p.y * p.y).sqrt();
                assert!((radial - 1.0).abs() < 1e-9, "off cylinder at ({u},{v})");
            }
        }
    }

    #[test]
    fn sweep_follows_quarter_circle_rail() {
        // Profile: small circle centred at the rail start (2,0,0) in the XY
        // plane. Rail: quarter circle in the XZ plane from (2,0,0) toward +Z.
        // Centring the profile on the rail start makes its centroid coincide
        // with the rail station, so the swept ring centre tracks the rail.
        let profile =
            NurbsCurve3D::circle(Point3::new(2.0, 0.0, 0.0), 0.2, Vector3::z(), Vector3::x())
                .unwrap();
        let rail = NurbsCurve3D::arc(
            Point3::origin(),
            2.0,
            -Vector3::y(),
            Vector3::x(),
            0.0,
            std::f64::consts::FRAC_PI_2,
        )
        .unwrap();
        let swept = NurbsSurface::sweep(&profile, &rail).unwrap();

        // Reconstruct the station rail points and their chord-length v
        // parameters exactly as `sweep` -> `loft` does (the profile centroid
        // equals the rail station, so the loft section parameters are the
        // chord-length parameters of the rail stations). At those v parameters
        // the loft reproduces each section exactly, so the swept ring centre
        // must equal the rail station point within the sweep approximation
        // tolerance (1e-6, NOT 1e-12).
        let stations = (2 * rail.control_points().len()).max(8);
        let (rt0, rt1) = rail.parameter_domain();
        let denom = (stations - 1) as f64;
        let station_pts: Vec<Point3> = (0..stations)
            .map(|i| {
                let t = rt0 + (rt1 - rt0) * (i as f64) / denom;
                rail.evaluate(t).unwrap()
            })
            .collect();
        let chords: Vec<f64> = station_pts
            .windows(2)
            .map(|w| (w[1] - w[0]).norm())
            .collect();
        let total: f64 = chords.iter().sum();
        let mut v_params = vec![0.0];
        let mut acc = 0.0;
        for c in &chords {
            acc += c;
            v_params.push(acc / total);
        }
        *v_params.last_mut().unwrap() = 1.0;

        let ((u0, u1), (v0, v1)) = swept.parameter_domain();
        for (station_pt, &vp) in station_pts.iter().zip(&v_params) {
            let v = v0 + (v1 - v0) * vp;
            // Average the surface ring over u to recover its centre.
            let samples = 24;
            let mut centre = Vector3::zeros();
            for s in 0..samples {
                let u = u0 + (u1 - u0) * (f64::from(s) / f64::from(samples));
                centre += swept.point_at(u, v).unwrap().coords;
            }
            centre /= f64::from(samples);
            assert!(
                (Point3::from(centre) - station_pt).norm() < 1e-6,
                "ring centre off rail at v={v}"
            );
        }
    }

    #[test]
    fn frames_stay_orthonormal_along_3d_rail() {
        // Interpolated wavy 3D rail.
        let rail = NurbsCurve3D::interpolate(
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.5),
                Point3::new(2.0, 0.0, 1.0),
                Point3::new(3.0, 1.0, 1.5),
                Point3::new(4.0, 0.0, 2.0),
            ],
            3,
        )
        .unwrap()
        .0;

        let stations = (2 * rail.control_points().len()).max(8);
        let (t0, t1) = rail.parameter_domain();
        let mut points = Vec::new();
        let mut tangents = Vec::new();
        for i in 0..stations {
            let t = t0 + (t1 - t0) * (i as f64) / ((stations - 1) as f64);
            let ders = rail.derivatives(t, 1).unwrap();
            points.push(Point3::from(ders[0]));
            tangents.push(ders[1].normalize());
        }
        let frames = rotation_minimizing_frames(&points, &tangents).unwrap();
        for (i, f) in frames.iter().enumerate() {
            assert!((f.tangent.norm() - 1.0).abs() < 1e-9, "t not unit at {i}");
            assert!((f.normal.norm() - 1.0).abs() < 1e-9, "n not unit at {i}");
            assert!((f.binormal.norm() - 1.0).abs() < 1e-9, "b not unit at {i}");
            assert!(f.tangent.dot(&f.normal).abs() < 1e-9, "t.n != 0 at {i}");
            assert!(f.tangent.dot(&f.binormal).abs() < 1e-9, "t.b != 0 at {i}");
            assert!(f.normal.dot(&f.binormal).abs() < 1e-9, "n.b != 0 at {i}");
        }
    }
}
