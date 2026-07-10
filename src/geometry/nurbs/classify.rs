//! Structural classification of NURBS surfaces into analytic forms.
//!
//! [`NurbsSurface::extrusion_form`] recognizes surfaces built by
//! [`NurbsSurface::extrude`] (The NURBS Book §8.3): `degree_v == 1`, two
//! v-rows where row 1 is row 0 translated by one common vector, with
//! per-column equal weights — exactly `S(u, v) = P(u) + t(v) · d` where
//! `t` maps the v-knot range linearly onto `[0, 1]`.
//!
//! [`NurbsSurface::parallelogram_form`] further recognizes the planar
//! special case (the profile is a single degree-1 segment with equal
//! weights): an affine patch `S(u, v) = origin + fu(u)·e_u + fv(v)·e_v`
//! whose parameterization can be inverted exactly.
//!
//! Both checks are structural and conservative: a `None` never rejects a
//! valid surface from general algorithms — callers fall back to the
//! numerical path. They exist so the surface×surface intersection can
//! dispatch closed-form fast paths for the prism-vs-prism booleans that
//! dominate the BIM workload (wall solids × opening cutters).

use crate::math::{Point3, Vector3, TOLERANCE};

use super::curve::NurbsCurve3D;
use super::surface::NurbsSurface;

/// A surface recognized as a translational extrusion:
/// `S(u, v) = profile(u) + t(v) · direction`, `t` linear from the v-knot
/// range onto `[0, 1]`.
#[derive(Debug, Clone)]
pub(crate) struct ExtrusionForm {
    /// The `v = v_min` boundary isocurve (u parameterization identical to the
    /// surface's u direction).
    pub profile: NurbsCurve3D,
    /// Total translation across the full v range.
    pub direction: Vector3,
}

/// A surface recognized as an affine parallelogram patch:
/// `S(u, v) = origin + fu(u)·e_u + fv(v)·e_v` with `fu`, `fv` mapping the
/// knot ranges linearly onto `[0, 1]`. Inversion of an in-plane point onto
/// `(fu, fv)` is an exact 2×2 solve.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ParallelogramForm {
    /// Corner at `(u_min, v_min)`.
    pub origin: Point3,
    /// Edge vector spanning the full u range at `v = v_min`.
    pub e_u: Vector3,
    /// Edge vector spanning the full v range at `u = u_min`.
    pub e_v: Vector3,
    /// Unit normal (`e_u × e_v`, normalized).
    pub normal: Vector3,
}

impl NurbsSurface {
    /// Recognizes a translational extrusion (see module docs).
    ///
    /// Returns `None` for anything that is not structurally a two-row
    /// degree-1 v-direction with one uniform row translation and equal
    /// per-column weights.
    pub(crate) fn extrusion_form(&self) -> Option<ExtrusionForm> {
        let (nu, nv) = self.grid_size();
        let (_, degree_v) = self.degrees();
        if degree_v != 1 || nv != 2 {
            return None;
        }

        // Scale-aware tolerance: constructor-exact data differs from the
        // ideal by rounding proportional to the coordinate magnitude.
        let (lo, hi) = self.bounding_box();
        let scale = (hi - lo).norm().max(1.0);
        let tol = TOLERANCE * scale;

        let direction = self.control_point(0, 1) - self.control_point(0, 0);
        for i in 0..nu {
            let d_i = self.control_point(i, 1) - self.control_point(i, 0);
            if (d_i - direction).norm() > tol {
                return None;
            }
            if (self.weight(i, 1) - self.weight(i, 0)).abs() > TOLERANCE {
                return None;
            }
        }
        if direction.norm() < tol {
            return None;
        }

        let (degree_u, _) = self.degrees();
        let profile_points: Vec<Point3> = (0..nu).map(|i| *self.control_point(i, 0)).collect();
        let profile_weights: Vec<f64> = (0..nu).map(|i| self.weight(i, 0)).collect();
        let profile = NurbsCurve3D::new(
            profile_points,
            profile_weights,
            self.knots_u().clone(),
            degree_u,
        )
        .ok()?;

        Some(ExtrusionForm { profile, direction })
    }

    /// Recognizes an affine parallelogram patch (see module docs): an
    /// extrusion whose profile is a single degree-1 segment with equal
    /// weights.
    pub(crate) fn parallelogram_form(&self) -> Option<ParallelogramForm> {
        let (nu, _) = self.grid_size();
        let (degree_u, _) = self.degrees();
        if degree_u != 1 || nu != 2 {
            return None;
        }
        // Rational degree-1 with unequal weights is a non-affine
        // reparameterization of the segment; the exact inversion below
        // requires the affine case.
        if (self.weight(0, 0) - self.weight(1, 0)).abs() > TOLERANCE {
            return None;
        }
        let form = self.extrusion_form()?;

        let origin = *self.control_point(0, 0);
        let e_u = self.control_point(1, 0) - self.control_point(0, 0);
        let e_v = form.direction;
        let n = e_u.cross(&e_v);
        let scale = e_u.norm().max(e_v.norm()).max(1.0);
        if n.norm() < TOLERANCE * scale * scale {
            return None;
        }
        Some(ParallelogramForm {
            origin,
            e_u,
            e_v,
            normal: n / n.norm(),
        })
    }
}

impl ParallelogramForm {
    /// Inverts an in-plane 3D point onto the normalized `(fu, fv)`
    /// fractions via the exact 2×2 Gram solve. The caller guarantees the
    /// point lies on the plane (analytic construction); any off-plane
    /// component is projected out.
    #[allow(clippy::many_single_char_names)]
    pub(crate) fn invert(&self, p: Point3) -> (f64, f64) {
        let r = p - self.origin;
        let a = self.e_u.dot(&self.e_u);
        let b = self.e_u.dot(&self.e_v);
        let c = self.e_v.dot(&self.e_v);
        let ru = r.dot(&self.e_u);
        let rv = r.dot(&self.e_v);
        let det = a * c - b * b;
        // `det > 0` is guaranteed by the non-degenerate cross product in
        // `parallelogram_form`.
        ((c * ru - b * rv) / det, (a * rv - b * ru) / det)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::nurbs::KnotVector;

    fn line(p0: Point3, p1: Point3) -> NurbsCurve3D {
        NurbsCurve3D::from_unweighted(
            vec![p0, p1],
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
        )
        .unwrap()
    }

    fn quarter_arc() -> NurbsCurve3D {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        NurbsCurve3D::new(
            vec![
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            vec![1.0, w, 1.0],
            KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]).unwrap(),
            2,
        )
        .unwrap()
    }

    #[test]
    fn extruded_arc_classifies_with_profile_and_direction() {
        let dir = Vector3::new(0.0, 0.0, 2.4);
        let s = NurbsSurface::extrude(&quarter_arc(), dir).unwrap();
        let form = s.extrusion_form().unwrap();
        assert!((form.direction - dir).norm() < 1e-12);
        for i in 0..=8 {
            let u = f64::from(i) / 8.0;
            let p_prof = form.profile.point_at(u).unwrap();
            let p_surf = s.point_at(u, 0.0).unwrap();
            assert!((p_prof - p_surf).norm() < 1e-12);
        }
    }

    #[test]
    fn extruded_line_classifies_as_parallelogram() {
        let dir = Vector3::new(0.0, 0.0, 3.0);
        let s = NurbsSurface::extrude(
            &line(Point3::new(1.0, 2.0, 0.0), Point3::new(4.0, 2.0, 0.0)),
            dir,
        )
        .unwrap();
        let para = s.parallelogram_form().unwrap();
        assert!((para.origin - Point3::new(1.0, 2.0, 0.0)).norm() < 1e-12);
        assert!((para.e_u - Vector3::new(3.0, 0.0, 0.0)).norm() < 1e-12);
        assert!((para.e_v - dir).norm() < 1e-12);
        // Inversion round-trips surface evaluation.
        for &(fu, fv) in &[(0.0, 0.0), (1.0, 1.0), (0.25, 0.75), (0.5, 0.5)] {
            let p = s.point_at(fu, fv).unwrap();
            let (gu, gv) = para.invert(p);
            assert!((gu - fu).abs() < 1e-12 && (gv - fv).abs() < 1e-12);
        }
    }

    #[test]
    fn extruded_arc_is_not_a_parallelogram() {
        let s = NurbsSurface::extrude(&quarter_arc(), Vector3::new(0.0, 0.0, 1.0)).unwrap();
        assert!(s.extrusion_form().is_some());
        assert!(s.parallelogram_form().is_none());
    }

    #[test]
    fn revolved_surface_is_not_an_extrusion() {
        let profile = line(Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 1.0));
        let s = NurbsSurface::revolve(
            &profile,
            Point3::origin(),
            Vector3::new(0.0, 0.0, 1.0),
            std::f64::consts::FRAC_PI_2,
        )
        .unwrap();
        assert!(s.extrusion_form().is_none());
    }

    #[test]
    fn non_uniform_translation_is_rejected() {
        // A bilinear patch whose two rows are NOT translates (twisted).
        let s = NurbsSurface::from_unweighted(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 0.0, 1.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 2.0),
            ],
            2,
            2,
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]).unwrap(),
            1,
            1,
        )
        .unwrap();
        assert!(s.extrusion_form().is_none());
        assert!(s.parallelogram_form().is_none());
    }
}
