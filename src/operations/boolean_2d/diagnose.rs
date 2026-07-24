//! Health assessment for 2D boolean results.
//!
//! Pure detection over [`PolygonWithHoles`] outputs, surfaced through the
//! [`crate::diagnostics`] contract so the app can log a diagnostic without the
//! kernel owning any logging policy.
//!
//! The arrangement engine already guarantees simple, CDT-safe output *or* an
//! `Err` (see the module docs), so on `Ok` the residual signals are: an
//! unexpectedly empty result (the base was fully consumed by the cut) and
//! near-zero-area sliver faces. A defensive non-finite check guards against
//! future engine changes.

use super::types::{signed_area, Polygon, PolygonWithHoles, WALL_EPS};
use crate::diagnostics::{OpHealth, Reason};

/// Outer-ring area below this is treated as a collapsed sliver (m^2).
pub const MIN_FACE_AREA: f64 = WALL_EPS;

/// Assess a 2D boolean result. `expect_nonempty` = the op should have left at
/// least one face (e.g. subtracting a cut that must not consume the whole base).
#[must_use]
pub fn assess(faces: &[PolygonWithHoles], expect_nonempty: bool) -> OpHealth {
    let mut degenerate = Vec::new();
    let mut suspicious = Vec::new();

    if let Some(at) = first_nonfinite(faces) {
        degenerate.push(Reason::NonFinite { at });
    }
    if expect_nonempty && faces.is_empty() {
        degenerate.push(Reason::EmptyResult);
    }
    for face in faces {
        let area = signed_area(&face.outer).abs();
        if area < MIN_FACE_AREA {
            suspicious.push(Reason::ZeroAreaFace { scale: area });
        }
    }

    if !degenerate.is_empty() {
        OpHealth::Degenerate(degenerate)
    } else if !suspicious.is_empty() {
        OpHealth::Suspicious(suspicious)
    } else {
        OpHealth::Ok
    }
}

fn ring_finite(ring: &Polygon) -> bool {
    ring.iter().all(|&(x, y)| x.is_finite() && y.is_finite())
}

/// Where the first non-finite coordinate appears, if any.
fn first_nonfinite(faces: &[PolygonWithHoles]) -> Option<&'static str> {
    for face in faces {
        if !ring_finite(&face.outer) {
            return Some("outer ring");
        }
        if face.holes.iter().any(|h| !ring_finite(h)) {
            return Some("hole ring");
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A CCW square of side `side` (area = side^2).
    fn square(side: f64) -> PolygonWithHoles {
        PolygonWithHoles {
            outer: vec![(0.0, 0.0), (side, 0.0), (side, side), (0.0, side)],
            holes: Vec::new(),
        }
    }

    #[test]
    fn clean_result_is_ok() {
        assert_eq!(assess(&[square(1.0)], true), OpHealth::Ok);
    }

    #[test]
    fn expected_nonempty_but_empty_is_degenerate() {
        assert_eq!(
            assess(&[], true),
            OpHealth::Degenerate(vec![Reason::EmptyResult])
        );
        // Emptiness is not flagged when it is an acceptable outcome.
        assert_eq!(assess(&[], false), OpHealth::Ok);
    }

    #[test]
    fn sliver_face_is_suspicious() {
        // side 1e-4 -> area 1e-8 < MIN_FACE_AREA (1e-6).
        match assess(&[square(1e-4)], true) {
            OpHealth::Suspicious(reasons) => {
                assert!(matches!(reasons.as_slice(), [Reason::ZeroAreaFace { .. }]));
            }
            other => panic!("expected Suspicious, got {other:?}"),
        }
    }

    #[test]
    fn non_finite_is_degenerate() {
        let mut face = square(1.0);
        face.outer[1].0 = f64::NAN;
        assert_eq!(
            assess(&[face], true),
            OpHealth::Degenerate(vec![Reason::NonFinite { at: "outer ring" }])
        );
    }
}
