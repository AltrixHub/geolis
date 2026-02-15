use thiserror::Error;

/// Top-level error type for the Geolis CAD kernel.
#[derive(Debug, Error)]
pub enum GeolisError {
    #[error(transparent)]
    Geometry(#[from] GeometryError),

    #[error(transparent)]
    Topology(#[from] TopologyError),

    #[error(transparent)]
    Operation(#[from] OperationError),

    #[error(transparent)]
    Tessellation(#[from] TessellationError),
}

/// Errors related to geometric computations.
#[derive(Debug, Error)]
pub enum GeometryError {
    #[error("parameter {parameter} = {value} is out of range [{min}, {max}]")]
    ParameterOutOfRange {
        parameter: &'static str,
        value: f64,
        min: f64,
        max: f64,
    },

    #[error("degenerate geometry: {0}")]
    Degenerate(String),

    #[error("zero-length vector")]
    ZeroVector,
}

/// Errors related to topological operations.
#[derive(Debug, Error)]
pub enum TopologyError {
    #[error("entity not found: {0}")]
    EntityNotFound(String),

    #[error("wire is not closed")]
    WireNotClosed,

    #[error("invalid topology: {0}")]
    InvalidTopology(String),
}

/// Errors related to CAD operations.
#[derive(Debug, Error)]
pub enum OperationError {
    #[error("invalid input: {0}")]
    InvalidInput(String),

    #[error("operation failed: {0}")]
    Failed(String),
}

/// Errors related to tessellation.
#[derive(Debug, Error)]
pub enum TessellationError {
    #[error("invalid tessellation parameters: {0}")]
    InvalidParameters(String),

    #[error("tessellation failed: {0}")]
    Failed(String),
}

/// Convenience type alias for results using [`GeolisError`].
pub type Result<T> = std::result::Result<T, GeolisError>;
