//! NURBS curve and surface constructors.
//!
//! Each submodule attaches associated constructor functions to [`NurbsCurve3D`]
//! or [`NurbsSurface`] (The NURBS Book, chapters 7-10):
//!
//! - Curves: [`NurbsCurve3D::circle`] / [`NurbsCurve3D::arc`] (A7.1),
//!   [`NurbsCurve3D::polyline`], [`NurbsCurve3D::interpolate`] (A9.1).
//! - Surfaces: [`NurbsSurface::extrude`] (§8.3), [`NurbsSurface::revolve`]
//!   (A8.1), [`NurbsSurface::loft`] (§10.3), [`NurbsSurface::sweep`] (§10.4).
//!
//! [`NurbsCurve3D`]: super::NurbsCurve3D
//! [`NurbsSurface`]: super::NurbsSurface

mod circle;
mod extrude;
mod interpolate;
mod loft;
mod polyline;
mod revolve;
mod sweep;
