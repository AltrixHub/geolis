use super::knot::KnotVector;
#[cfg(test)]
use crate::math::TOLERANCE;

/// Evaluates the `degree + 1` non-vanishing B-spline basis functions at `u`
/// for the given knot span (The NURBS Book, A2.2).
///
/// Returns `[N_{span-degree,degree}(u), ..., N_{span,degree}(u)]`.
#[must_use]
pub fn basis_functions(knots: &KnotVector, span: usize, u: f64, degree: usize) -> Vec<f64> {
    let k = knots.as_slice();
    let mut n = vec![0.0; degree + 1];
    let mut left = vec![0.0; degree + 1];
    let mut right = vec![0.0; degree + 1];
    n[0] = 1.0;
    for j in 1..=degree {
        left[j] = u - k[span + 1 - j];
        right[j] = k[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            let temp = n[r] / (right[r + 1] + left[j - r]);
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    n
}

/// Binomial coefficient C(n, k) as f64.
#[must_use]
pub fn binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut result = 1.0;
    for i in 0..k {
        result = result * ((n - i) as f64) / ((i + 1) as f64);
    }
    result
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn book_knots() -> KnotVector {
        KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0]).unwrap()
    }

    #[test]
    fn book_example_2_3() {
        // The NURBS Book Ex2.3: p=2, u=5/2, span 4 -> N = [1/8, 6/8, 1/8]
        let k = book_knots();
        let n = basis_functions(&k, 4, 2.5, 2);
        assert!((n[0] - 0.125).abs() < TOLERANCE);
        assert!((n[1] - 0.75).abs() < TOLERANCE);
        assert!((n[2] - 0.125).abs() < TOLERANCE);
    }

    #[test]
    fn partition_of_unity() {
        let k = book_knots();
        for i in 0..=50 {
            let u = 5.0 * f64::from(i) / 50.0;
            let span = k.find_span(2, 8, u);
            let n = basis_functions(&k, span, u, 2);
            let sum: f64 = n.iter().sum();
            assert!((sum - 1.0).abs() < 1e-12, "sum {sum} at u={u}");
            assert!(n.iter().all(|&x| x >= -1e-12), "negative basis at u={u}");
        }
    }

    #[test]
    fn binomial_values() {
        assert!((binomial(4, 0) - 1.0).abs() < TOLERANCE);
        assert!((binomial(4, 2) - 6.0).abs() < TOLERANCE);
        assert!((binomial(5, 3) - 10.0).abs() < TOLERANCE);
    }
}
