mod assemble;
mod classify;
mod engine;
mod face_intersection;
mod intersect_op;
mod select;
mod split;
mod subtract;
mod union;

pub use classify::{classify_point_in_solid, PointClassification};
pub use face_intersection::{intersect_face_face, FaceFaceIntersection};
pub use intersect_op::Intersect;
pub use select::BooleanOp;
pub use split::{FaceFragment, SolidSource};
pub use subtract::Subtract;
pub use union::Union;
