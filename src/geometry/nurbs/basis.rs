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

/// Evaluates the non-vanishing basis functions and their derivatives up to
/// `num_derivatives` at `u` (The NURBS Book, A2.3).
///
/// Returns `ders` where `ders[k][j]` is the k-th derivative of
/// `N_{span-degree+j,degree}` at `u`. Derivatives of order greater than
/// `degree` are zero.
#[must_use]
pub fn basis_function_derivatives(
    knots: &KnotVector,
    span: usize,
    u: f64,
    degree: usize,
    num_derivatives: usize,
) -> Vec<Vec<f64>> {
    let k = knots.as_slice();
    let p = degree;
    let n = num_derivatives.min(p);

    // ndu[j][r] stores basis functions and knot differences (A2.3 layout).
    let mut ndu = vec![vec![0.0; p + 1]; p + 1];
    let mut left = vec![0.0; p + 1];
    let mut right = vec![0.0; p + 1];
    ndu[0][0] = 1.0;
    for j in 1..=p {
        left[j] = u - k[span + 1 - j];
        right[j] = k[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            ndu[j][r] = right[r + 1] + left[j - r];
            let temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    let mut ders = vec![vec![0.0; p + 1]; num_derivatives + 1];
    for j in 0..=p {
        ders[0][j] = ndu[j][p];
    }

    let mut a = [vec![0.0; p + 1], vec![0.0; p + 1]];
    for r in 0..=p {
        let mut s1 = 0usize;
        let mut s2 = 1usize;
        a[0].iter_mut().for_each(|x| *x = 0.0);
        a[1].iter_mut().for_each(|x| *x = 0.0);
        a[0][0] = 1.0;
        for kk in 1..=n {
            let mut d = 0.0;
            let rk = r as isize - kk as isize;
            let pk = p - kk;
            if r >= kk {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk as usize];
                d = a[s2][0] * ndu[rk as usize][pk];
            }
            let j1 = if rk >= -1 { 1 } else { (-rk) as usize };
            let j2 = if r <= pk + 1 { kk - 1 } else { p - r };
            for j in j1..=j2 {
                let col = (rk + j as isize) as usize;
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][col];
                d += a[s2][j] * ndu[col][pk];
            }
            if r <= pk {
                a[s2][kk] = -a[s1][kk - 1] / ndu[pk + 1][r];
                d += a[s2][kk] * ndu[r][pk];
            }
            ders[kk][r] = d;
            std::mem::swap(&mut s1, &mut s2);
        }
    }

    // Multiply by the factorial-style factors p!/(p-k)!
    let mut factor = p as f64;
    for kk in 1..=n {
        for j in 0..=p {
            ders[kk][j] *= factor;
        }
        factor *= (p - kk) as f64;
    }
    ders
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

    #[test]
    fn derivatives_row_zero_matches_basis() {
        let k = book_knots();
        let n = basis_functions(&k, 4, 2.5, 2);
        let d = basis_function_derivatives(&k, 4, 2.5, 2, 2);
        for j in 0..=2 {
            assert!((d[0][j] - n[j]).abs() < TOLERANCE);
        }
    }

    #[test]
    fn first_derivatives_match_central_difference() {
        let k = book_knots();
        let h = 1e-6;
        let u = 2.5;
        let span = k.find_span(2, 8, u);
        let d = basis_function_derivatives(&k, span, u, 2, 1);
        let np = basis_functions(&k, span, u + h, 2);
        let nm = basis_functions(&k, span, u - h, 2);
        for j in 0..=2 {
            let fd = (np[j] - nm[j]) / (2.0 * h);
            assert!((d[1][j] - fd).abs() < 1e-5, "j={j}: {} vs {fd}", d[1][j]);
        }
    }

    #[test]
    fn derivative_rows_sum_to_zero() {
        // d/du of partition of unity is 0
        let k = book_knots();
        let u = 1.7;
        let span = k.find_span(2, 8, u);
        let d = basis_function_derivatives(&k, span, u, 2, 2);
        let s1: f64 = d[1].iter().sum();
        let s2: f64 = d[2].iter().sum();
        assert!(s1.abs() < 1e-10);
        assert!(s2.abs() < 1e-10);
    }
}
