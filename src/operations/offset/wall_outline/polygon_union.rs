//! Backwards-compatibility shim — the engine has moved to
//! [`crate::operations::boolean_2d`].
//!
//! This module re-exports the subset of symbols that the surrounding
//! `wall_outline` code (and its tests) still consume by the old path
//! `super::polygon_union::...`. Update call sites to import directly
//! from `crate::operations::boolean_2d` instead of touching this file.

pub(crate) use crate::operations::boolean_2d::{
    point_in_polygon_class, seg_seg_intersect, PointClass, Polygon, PolygonWithHoles, WALL_EPS,
    WALL_EPS_SQ,
};

// Test-only re-export: the P3 oracle in `wall_outline::tests` calls
// `polygon_union::signed_area` to detect CCW/CW boundary windings.
#[cfg(test)]
pub(crate) use crate::operations::boolean_2d::signed_area;
