use crate::error::{Result, TessellationError};

/// Line join style at polyline vertices.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum LineJoin {
    /// Extend edges to their intersection point (clamped to prevent spikes).
    Miter,
    /// Flat cap at every interior vertex — no miter extension.
    Bevel,
    /// Miter for gentle angles, bevel for sharp angles (default).
    #[default]
    Auto,
}

/// Style parameters for polyline stroke tessellation.
#[derive(Debug, Clone, Copy)]
pub struct StrokeStyle {
    width: f64,
    line_join: LineJoin,
}

impl StrokeStyle {
    /// Creates a new stroke style with the given width and default join (`Auto`).
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
        Ok(Self {
            width,
            line_join: LineJoin::default(),
        })
    }

    /// Sets the line join style.
    #[must_use]
    pub fn with_line_join(mut self, join: LineJoin) -> Self {
        self.line_join = join;
        self
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

    /// Returns the line join style.
    #[must_use]
    pub fn line_join(&self) -> LineJoin {
        self.line_join
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
        assert_eq!(style.line_join(), LineJoin::Auto);
    }

    #[test]
    fn with_line_join_sets_join() {
        let style = StrokeStyle::new(1.0).unwrap().with_line_join(LineJoin::Bevel);
        assert_eq!(style.line_join(), LineJoin::Bevel);
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
