use crate::error::{GeometryError, Result};
use crate::math::TOLERANCE;

/// Non-decreasing knot sequence for B-spline / NURBS bases.
///
/// Invariants enforced at construction: at least 4 knots (the minimum for a
/// degree-1 curve with 2 control points), all values finite, non-decreasing.
#[derive(Debug, Clone, PartialEq)]
pub struct KnotVector(Vec<f64>);

impl KnotVector {
    /// Creates a knot vector, validating the invariants.
    ///
    /// # Errors
    ///
    /// Returns an error if there are fewer than 4 knots, any knot is not
    /// finite, or the sequence is decreasing anywhere.
    pub fn new(knots: Vec<f64>) -> Result<Self> {
        if knots.len() < 4 {
            return Err(
                GeometryError::Degenerate("knot vector needs at least 4 knots".into()).into(),
            );
        }
        if knots.iter().any(|k| !k.is_finite()) {
            return Err(GeometryError::Degenerate("knot values must be finite".into()).into());
        }
        if knots.windows(2).any(|w| w[1] < w[0]) {
            return Err(
                GeometryError::Degenerate("knot vector must be non-decreasing".into()).into(),
            );
        }
        Ok(Self(knots))
    }

    /// Number of knots.
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Whether the knot vector is empty (never true for a valid instance).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Knot values as a slice.
    #[must_use]
    pub fn as_slice(&self) -> &[f64] {
        &self.0
    }

    /// Parameter domain `[knots[degree], knots[len - 1 - degree]]` of a
    /// curve of the given degree over this knot vector.
    #[must_use]
    pub fn domain(&self, degree: usize) -> (f64, f64) {
        (self.0[degree], self.0[self.0.len() - 1 - degree])
    }

    /// Finds the knot span index containing `u` (The NURBS Book, A2.1).
    ///
    /// `control_count` is the number of control points of the curve using
    /// this knot vector. `u` is assumed to lie within the domain.
    #[must_use]
    pub fn find_span(&self, degree: usize, control_count: usize, u: f64) -> usize {
        let n = control_count - 1;
        let knots = &self.0;
        if u >= knots[n + 1] {
            return n;
        }
        if u <= knots[degree] {
            return degree;
        }
        let mut low = degree;
        let mut high = n + 1;
        let mut mid = usize::midpoint(low, high);
        while u < knots[mid] || u >= knots[mid + 1] {
            if u < knots[mid] {
                high = mid;
            } else {
                low = mid;
            }
            mid = usize::midpoint(low, high);
        }
        mid
    }

    /// Multiplicity of `u` in the knot vector (tolerance-based equality).
    #[must_use]
    pub fn multiplicity(&self, u: f64) -> usize {
        self.0
            .iter()
            .filter(|&&k| (k - u).abs() < TOLERANCE)
            .count()
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    // Knots from The NURBS Book, Ex2.3: U = {0,0,0,1,2,3,4,4,5,5,5}, p = 2.
    fn book_knots() -> KnotVector {
        KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0]).unwrap()
    }

    #[test]
    fn rejects_decreasing_knots() {
        assert!(KnotVector::new(vec![0.0, 1.0, 0.5, 2.0]).is_err());
    }

    #[test]
    fn rejects_non_finite_knots() {
        assert!(KnotVector::new(vec![0.0, f64::NAN, 1.0, 2.0]).is_err());
    }

    #[test]
    fn rejects_too_short() {
        assert!(KnotVector::new(vec![0.0, 1.0, 2.0]).is_err());
    }

    #[test]
    fn find_span_book_example() {
        // 11 knots, p=2 -> control_count = 11 - 2 - 1 = 8; u = 2.5 lies in [knots[4], knots[5])
        assert_eq!(book_knots().find_span(2, 8, 2.5), 4);
    }

    #[test]
    fn find_span_at_domain_end_returns_last_span() {
        // u = 5.0 (domain end) must return n = 7, not run past the array
        assert_eq!(book_knots().find_span(2, 8, 5.0), 7);
    }

    #[test]
    fn domain_is_clamped_range() {
        let (t0, t1) = book_knots().domain(2);
        assert!((t0 - 0.0).abs() < TOLERANCE);
        assert!((t1 - 5.0).abs() < TOLERANCE);
    }

    #[test]
    fn multiplicity_counts_equal_knots() {
        let k = book_knots();
        assert_eq!(k.multiplicity(4.0), 2);
        assert_eq!(k.multiplicity(2.0), 1);
        assert_eq!(k.multiplicity(2.5), 0);
    }
}
