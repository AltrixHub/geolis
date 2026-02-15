use crate::error::{Result, TessellationError};

/// Style parameters for polyline stroke tessellation.
#[derive(Debug, Clone, Copy)]
pub struct StrokeStyle {
    width: f64,
}

impl StrokeStyle {
    /// Creates a new stroke style.
    ///
    /// # Errors
    ///
    /// Returns an error if `width` is not positive.
    pub fn new(width: f64) -> Result<Self> {
        if width <= 0.0 {
            return Err(TessellationError::InvalidParameters(
                "stroke width must be positive".to_owned(),
            )
            .into());
        }
        Ok(Self { width })
    }

    /// Returns the stroke width.
    #[must_use]
    pub fn width(&self) -> f64 {
        self.width
    }

    /// Returns half the stroke width.
    #[must_use]
    pub fn half_width(&self) -> f64 {
        self.width * 0.5
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn new_with_valid_width() {
        let style = StrokeStyle::new(2.0).unwrap();
        assert!((style.width() - 2.0).abs() < f64::EPSILON);
        assert!((style.half_width() - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn new_with_zero_width_fails() {
        let result = StrokeStyle::new(0.0);
        assert!(result.is_err());
    }

    #[test]
    fn new_with_negative_width_fails() {
        let result = StrokeStyle::new(-1.0);
        assert!(result.is_err());
    }
}
