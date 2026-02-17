pub mod arc_2d;
pub mod distance_2d;
pub mod intersect_2d;
pub mod polygon_2d;

/// 2D point type.
pub type Point2 = nalgebra::Point2<f64>;

/// 3D point type.
pub type Point3 = nalgebra::Point3<f64>;

/// 2D vector type.
pub type Vector2 = nalgebra::Vector2<f64>;

/// 3D vector type.
pub type Vector3 = nalgebra::Vector3<f64>;

/// 4x4 transformation matrix.
pub type Matrix4 = nalgebra::Matrix4<f64>;

/// Global geometric tolerance for floating-point comparisons.
pub const TOLERANCE: f64 = 1e-10;
