use super::classify::PointClassification;
use super::split::SolidSource;

/// The type of boolean operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanOp {
    Union,
    Subtract,
    Intersect,
}

/// Decision about whether to keep a fragment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KeepDecision {
    Keep,
    KeepFlipped,
    Discard,
}

/// Determines whether a fragment should be kept based on its classification
/// relative to the other solid and the boolean operation.
///
/// | Fragment | vs Other Solid | Union | Subtract(A-B) | Intersect |
/// |----------|---------------|-------|----------------|-----------|
/// | from A   | OUTSIDE B     | keep  | keep           | discard   |
/// | from A   | INSIDE B      | discard | discard      | keep      |
/// | from B   | OUTSIDE A     | keep  | discard        | discard   |
/// | from B   | INSIDE A      | discard | keep (flip)  | keep      |
#[allow(clippy::match_same_arms)]
#[must_use]
pub fn should_keep_fragment(
    source: SolidSource,
    classification: PointClassification,
    op: BooleanOp,
) -> KeepDecision {
    // Each arm is kept explicit for readability — the decision table matches the
    // algorithm specification directly, even though some arms share the same body.
    match (source, classification, op) {
        // Fragment from A, classified vs B
        (SolidSource::A, PointClassification::Outside, BooleanOp::Union) => KeepDecision::Keep,
        (SolidSource::A, PointClassification::Outside, BooleanOp::Subtract) => KeepDecision::Keep,
        (SolidSource::A, PointClassification::Outside, BooleanOp::Intersect) => {
            KeepDecision::Discard
        }

        (SolidSource::A, PointClassification::Inside, BooleanOp::Union) => KeepDecision::Discard,
        (SolidSource::A, PointClassification::Inside, BooleanOp::Subtract) => {
            KeepDecision::Discard
        }
        (SolidSource::A, PointClassification::Inside, BooleanOp::Intersect) => KeepDecision::Keep,

        // Fragment from B, classified vs A
        (SolidSource::B, PointClassification::Outside, BooleanOp::Union) => KeepDecision::Keep,
        (SolidSource::B, PointClassification::Outside, BooleanOp::Subtract) => {
            KeepDecision::Discard
        }
        (SolidSource::B, PointClassification::Outside, BooleanOp::Intersect) => {
            KeepDecision::Discard
        }

        (SolidSource::B, PointClassification::Inside, BooleanOp::Union) => KeepDecision::Discard,
        (SolidSource::B, PointClassification::Inside, BooleanOp::Subtract) => {
            KeepDecision::KeepFlipped
        }
        (SolidSource::B, PointClassification::Inside, BooleanOp::Intersect) => KeepDecision::Keep,

        // OnBoundary — keep from A to avoid duplicates, discard from B
        (SolidSource::A, PointClassification::OnBoundary, BooleanOp::Union) => KeepDecision::Keep,
        (SolidSource::A, PointClassification::OnBoundary, BooleanOp::Subtract) => {
            KeepDecision::Keep
        }
        (SolidSource::A, PointClassification::OnBoundary, BooleanOp::Intersect) => {
            KeepDecision::Keep
        }
        (SolidSource::B, PointClassification::OnBoundary, _) => KeepDecision::Discard,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn union_keeps_outside_fragments() {
        assert_eq!(
            should_keep_fragment(SolidSource::A, PointClassification::Outside, BooleanOp::Union),
            KeepDecision::Keep
        );
        assert_eq!(
            should_keep_fragment(SolidSource::B, PointClassification::Outside, BooleanOp::Union),
            KeepDecision::Keep
        );
    }

    #[test]
    fn union_discards_inside_fragments() {
        assert_eq!(
            should_keep_fragment(SolidSource::A, PointClassification::Inside, BooleanOp::Union),
            KeepDecision::Discard
        );
        assert_eq!(
            should_keep_fragment(SolidSource::B, PointClassification::Inside, BooleanOp::Union),
            KeepDecision::Discard
        );
    }

    #[test]
    fn subtract_keeps_a_outside_discards_b_outside() {
        assert_eq!(
            should_keep_fragment(
                SolidSource::A,
                PointClassification::Outside,
                BooleanOp::Subtract
            ),
            KeepDecision::Keep
        );
        assert_eq!(
            should_keep_fragment(
                SolidSource::B,
                PointClassification::Outside,
                BooleanOp::Subtract
            ),
            KeepDecision::Discard
        );
    }

    #[test]
    fn subtract_flips_b_inside() {
        assert_eq!(
            should_keep_fragment(
                SolidSource::B,
                PointClassification::Inside,
                BooleanOp::Subtract
            ),
            KeepDecision::KeepFlipped
        );
    }

    #[test]
    fn intersect_keeps_inside_fragments() {
        assert_eq!(
            should_keep_fragment(
                SolidSource::A,
                PointClassification::Inside,
                BooleanOp::Intersect
            ),
            KeepDecision::Keep
        );
        assert_eq!(
            should_keep_fragment(
                SolidSource::B,
                PointClassification::Inside,
                BooleanOp::Intersect
            ),
            KeepDecision::Keep
        );
    }
}
