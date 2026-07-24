//! Operation diagnostics — the kernel *reports* health as data; the caller
//! (the app) owns logging policy and serialization.
//!
//! A modeling op returns its value together with an [`OpHealth`] verdict. The
//! value is always produced best-effort (even when degenerate); the optional
//! [`InputSnapshot`] carries a readable summary of the inputs and is populated
//! **only when health is not `Ok`**, so the normal path costs nothing.
//!
//! geolis stays serialization-free: [`OpHealth`] / [`Reason`] / [`InputSnapshot`]
//! are plain data the app layer (which owns `serde`) turns into a diagnostic log
//! line or a replay dump. See the modeling-diagnostics design doc.

/// A modeling-op result paired with a health verdict.
#[derive(Debug, Clone)]
pub struct OpDiagnostic<T> {
    /// The op output (best-effort even when degenerate).
    pub value: T,
    /// Health verdict for the op.
    pub health: OpHealth,
    /// Readable input summary, captured only when `health != Ok`.
    pub inputs: Option<InputSnapshot>,
}

impl<T> OpDiagnostic<T> {
    /// A clean result: `Ok` health, no captured inputs.
    #[must_use]
    pub fn ok(value: T) -> Self {
        Self {
            value,
            health: OpHealth::Ok,
            inputs: None,
        }
    }

    /// A flagged result: a non-`Ok` verdict with its input snapshot.
    #[must_use]
    pub fn flagged(value: T, health: OpHealth, inputs: InputSnapshot) -> Self {
        Self {
            value,
            health,
            inputs: Some(inputs),
        }
    }

    /// Whether the verdict is `Ok`.
    pub fn is_ok(&self) -> bool {
        self.health.is_ok()
    }

    /// Transform the value, preserving the verdict and snapshot.
    #[must_use]
    pub fn map<U>(self, f: impl FnOnce(T) -> U) -> OpDiagnostic<U> {
        OpDiagnostic {
            value: f(self.value),
            health: self.health,
            inputs: self.inputs,
        }
    }
}

/// Health verdict for an operation.
#[derive(Debug, Clone, PartialEq)]
pub enum OpHealth {
    /// Output is well-formed.
    Ok,
    /// Succeeded, but numerically suspect — worth a warning, not a failure.
    Suspicious(Vec<Reason>),
    /// Output is broken (empty / non-finite / self-intersecting / ...).
    Degenerate(Vec<Reason>),
    /// The op could not complete; carries the error message.
    Failed(String),
}

impl OpHealth {
    /// Whether the verdict is `Ok`.
    #[must_use]
    pub fn is_ok(&self) -> bool {
        matches!(self, OpHealth::Ok)
    }

    /// The reasons behind a non-`Ok` verdict (empty for `Ok` / `Failed`).
    #[must_use]
    pub fn reasons(&self) -> &[Reason] {
        match self {
            OpHealth::Suspicious(r) | OpHealth::Degenerate(r) => r,
            OpHealth::Ok | OpHealth::Failed(_) => &[],
        }
    }
}

/// A specific, numeric reason an op's output is suspect or degenerate.
///
/// Extended per op as detection is wired in; each variant carries the number
/// that made the call so a reader can judge severity without re-deriving it.
#[derive(Debug, Clone, PartialEq)]
pub enum Reason {
    /// No geometry was produced where output was expected.
    EmptyResult,
    /// A coordinate or scalar was NaN or infinite.
    NonFinite {
        /// Where the non-finite value appeared, e.g. `"outer vertex"`.
        at: &'static str,
    },
    /// The result self-intersects.
    SelfIntersection,
    /// A face collapsed to (near) zero area.
    ZeroAreaFace {
        /// The offending area / scale.
        scale: f64,
    },
    /// Two vertices sit within an unsafe distance of each other.
    NearCoincidentVertices {
        /// The measured separation.
        dist: f64,
    },
    /// A numerically ill-conditioned step (high condition number).
    HighConditionNumber {
        /// The measured condition number.
        value: f64,
    },
}

/// A readable snapshot of an op's inputs, for a diagnostic log line and (later)
/// a replay dump.
///
/// Serialization-free here: `summary` pairs are plain strings the app can print
/// or serialize. Reproduction-grade geometry payload is attached where each op
/// is wired (a follow-up once the app-side serialization boundary is settled).
#[derive(Debug, Clone, PartialEq)]
pub struct InputSnapshot {
    /// The op that produced this snapshot, e.g. `"boolean_2d::subtract"`.
    pub op: &'static str,
    /// Readable key facts (bbox, counts, params) — one `(key, value)` per fact.
    pub summary: Vec<(String, String)>,
}

impl InputSnapshot {
    /// An empty snapshot for `op`.
    #[must_use]
    pub fn new(op: &'static str) -> Self {
        Self {
            op,
            summary: Vec::new(),
        }
    }

    /// Add one readable fact (builder style).
    #[must_use]
    pub fn with(mut self, key: impl Into<String>, value: impl std::fmt::Display) -> Self {
        self.summary.push((key.into(), value.to_string()));
        self
    }
}

#[cfg(test)]
#[allow(
    clippy::expect_used,
    reason = "expect is the failure signal in these unit tests"
)]
mod tests {
    use super::*;

    #[test]
    fn ok_diagnostic_has_no_inputs() {
        let d = OpDiagnostic::ok(42);
        assert!(d.is_ok());
        assert!(d.inputs.is_none());
        assert_eq!(d.health.reasons(), &[]);
    }

    #[test]
    fn flagged_carries_reasons_and_snapshot() {
        let snap = InputSnapshot::new("boolean_2d::subtract")
            .with("base_area", 12.5)
            .with("holes", 3);
        let d = OpDiagnostic::flagged(
            Vec::<u8>::new(),
            OpHealth::Degenerate(vec![Reason::EmptyResult]),
            snap,
        );
        assert!(!d.is_ok());
        assert_eq!(d.health.reasons(), &[Reason::EmptyResult]);
        let inputs = d.inputs.expect("flagged snapshot present");
        assert_eq!(inputs.op, "boolean_2d::subtract");
        assert_eq!(inputs.summary.len(), 2);
        assert_eq!(
            inputs.summary[0],
            ("base_area".to_string(), "12.5".to_string())
        );
    }

    #[test]
    fn map_preserves_verdict_and_snapshot() {
        let d = OpDiagnostic::flagged(
            3u32,
            OpHealth::Suspicious(vec![Reason::ZeroAreaFace { scale: 1e-9 }]),
            InputSnapshot::new("x"),
        );
        let mapped = d.map(|v| u64::from(v) + 1);
        assert_eq!(mapped.value, 4u64);
        assert!(matches!(mapped.health, OpHealth::Suspicious(_)));
        assert!(mapped.inputs.is_some());
    }

    #[test]
    fn health_is_ok_only_for_ok() {
        assert!(OpHealth::Ok.is_ok());
        assert!(!OpHealth::Failed("boom".to_string()).is_ok());
        assert!(!OpHealth::Degenerate(vec![]).is_ok());
        assert!(!OpHealth::Suspicious(vec![Reason::SelfIntersection]).is_ok());
    }
}
