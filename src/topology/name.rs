//! Persistent topology names (topological naming).
//!
//! Slotmap ids are fresh on every rebuild, so cross-rebuild references (an
//! opening on a wall face, a per-face material) need rebuild-stable names.
//! Names are **derivational**: a pure function of the creating operation's
//! identity ([`OpId`], supplied by the caller — geolis never invents one),
//! the entity's role within that operation, and — for boolean products — the
//! parent names. Same inputs, same names, independent of allocation order.
//!
//! A boolean's result carries the target's names forward UNCHANGED (a punched
//! wall face is still "the wall's outer face"); new faces (band, pocket
//! floor) and new edges (hole rims) get names composed from their parents.
//! Resolution failure is `None` — no geometric best-match heuristics
//! (Kripac 1997 / Chen & Hoffmann 1995 / OCCT TNaming motivate the problem;
//! geolis's small deterministic op set lets derivational names replace their
//! heavyweight history machinery).

use std::collections::HashMap;
use std::fmt;
use std::sync::Arc;

use super::edge::EdgeId;
use super::face::FaceId;

/// Identity of the graph operation that created an entity, supplied by the
/// caller (revion passes its cognet node id). Rebuild-stable by construction.
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub struct OpId(Arc<str>);

impl OpId {
    /// Creates an operation id from the caller's stable identifier.
    pub fn new(id: impl Into<Arc<str>>) -> Self {
        Self(id.into())
    }

    /// The raw identifier.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for OpId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

/// A caller-supplied, opaque tag identifying one profile segment of a
/// segmented creation op (e.g. [`MakeSegmentedPrism`]). Like [`OpId`], geolis
/// never invents the identity: the caller derives tags from its own outline
/// provenance, so a face keeps its name when positional segment indices shift
/// (junction re-trims). Rebuild-stable by construction.
///
/// [`MakeSegmentedPrism`]: crate::operations::creation::MakeSegmentedPrism
#[derive(Clone, Debug, Hash, Eq, PartialEq, Ord, PartialOrd)]
pub struct SegmentTag(Arc<str>);

impl SegmentTag {
    /// Creates a segment tag from the caller's stable identifier.
    pub fn new(tag: impl Into<Arc<str>>) -> Self {
        Self(tag.into())
    }

    /// The raw identifier.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for SegmentTag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

/// The role of a face within its creation operation.
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub enum FaceRole {
    /// A side face, indexed by the op's deterministic side order (e.g. the
    /// curved wall: 0 = inner, 1 = outer, 2 = start end, 3 = end end).
    Side(u8),
    /// A side face identified by a caller-supplied segment tag
    /// (junction-stable — survives segment-count changes, unlike `Side(k)`).
    Tagged(SegmentTag),
    /// The cap at the extrusion start (`v0` end).
    CapStart,
    /// The cap at the extrusion end (`v1` end).
    CapEnd,
    /// The top face (slab / wall).
    Top,
    /// The bottom face (slab / wall).
    Bottom,
    /// The revolved wall surface.
    Wall,
}

/// The role of an edge within its creation operation.
#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq)]
pub enum EdgeRole {
    /// The shared ring edge at the extrusion start (`v0` / first profile end).
    RingStart,
    /// The shared ring edge at the extrusion end (`v1` / last profile end).
    RingEnd,
}

/// A persistent, rebuild-stable name for a face.
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub enum FaceName {
    /// A face born in a creation operation.
    Created {
        /// The creating operation.
        op: OpId,
        /// The face's role within that operation.
        role: FaceRole,
    },
    /// The band (hole / pocket wall) a boolean carved with a tool face.
    Band {
        /// The boolean operation.
        op: OpId,
        /// The tool side face the band lies on.
        tool_face: Box<FaceName>,
        /// Deterministic loop index (loops sorted by mean tool-`v`).
        loop_index: u32,
    },
    /// The pocket floor: the buried tool cap kept (sense-flipped) by a
    /// pocket subtract.
    Floor {
        /// The boolean operation.
        op: OpId,
        /// The buried cap's name.
        cap: Box<FaceName>,
    },
}

/// A persistent, rebuild-stable name for an edge.
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub enum EdgeName {
    /// An edge born in a creation operation.
    Created {
        /// The creating operation.
        op: OpId,
        /// The edge's role within that operation.
        role: EdgeRole,
    },
    /// A hole-rim ring a boolean cut into a target face.
    CutRim {
        /// The boolean operation.
        op: OpId,
        /// The punched target face's name.
        target: Box<FaceName>,
        /// Deterministic loop index (loops sorted by mean tool-`v`).
        loop_index: u32,
    },
}

/// Bidirectional registry `PersistentName ↔ current slotmap id`.
///
/// Both directions stay bijective: registering a name that is already bound
/// rebinds it (the previous holder drops out), which is exactly the boolean
/// move semantics — the newest result owns the name.
#[derive(Debug, Default)]
pub struct NameRegistry {
    face_names: HashMap<FaceId, FaceName>,
    faces_by_name: HashMap<FaceName, FaceId>,
    edge_names: HashMap<EdgeId, EdgeName>,
    edges_by_name: HashMap<EdgeName, EdgeId>,
}

impl NameRegistry {
    /// Binds `name` to `face`, unbinding any previous holder of the name and
    /// any previous name of the face.
    pub fn bind_face(&mut self, face: FaceId, name: FaceName) {
        if let Some(old_face) = self.faces_by_name.remove(&name) {
            self.face_names.remove(&old_face);
        }
        if let Some(old_name) = self.face_names.remove(&face) {
            self.faces_by_name.remove(&old_name);
        }
        self.face_names.insert(face, name.clone());
        self.faces_by_name.insert(name, face);
    }

    /// Binds `name` to `edge`, unbinding any previous holders.
    pub fn bind_edge(&mut self, edge: EdgeId, name: EdgeName) {
        if let Some(old_edge) = self.edges_by_name.remove(&name) {
            self.edge_names.remove(&old_edge);
        }
        if let Some(old_name) = self.edge_names.remove(&edge) {
            self.edges_by_name.remove(&old_name);
        }
        self.edge_names.insert(edge, name.clone());
        self.edges_by_name.insert(name, edge);
    }

    /// Moves the name of `from` (if any) onto `to` — the boolean carry-over.
    pub fn transfer_face(&mut self, from: FaceId, to: FaceId) {
        if let Some(name) = self.face_names.remove(&from) {
            self.faces_by_name.remove(&name);
            self.bind_face(to, name);
        }
    }

    /// Resolves a face name to the current face id.
    #[must_use]
    pub fn face(&self, name: &FaceName) -> Option<FaceId> {
        self.faces_by_name.get(name).copied()
    }

    /// The current name of a face, if registered.
    #[must_use]
    pub fn name_of_face(&self, face: FaceId) -> Option<&FaceName> {
        self.face_names.get(&face)
    }

    /// Resolves an edge name to the current edge id.
    #[must_use]
    pub fn edge(&self, name: &EdgeName) -> Option<EdgeId> {
        self.edges_by_name.get(name).copied()
    }

    /// The current name of an edge, if registered.
    #[must_use]
    pub fn name_of_edge(&self, edge: EdgeId) -> Option<&EdgeName> {
        self.edge_names.get(&edge)
    }
}

// ---------------------------------------------------------------------------
// Canonical string form (opaque to consumers; stable for graph storage).
//
// Grammar (all tokens ASCII; op ids and segment tags are percent-escaped for
// `%(:)`):
//   face := "created:" op ":" role
//         | "band:" op ":" index ":(" face ")"
//         | "floor:" op ":(" face ")"
//   edge := "ring:" op ":" ("start" | "end")
//         | "rim:"  op ":" index ":(" face ")"
//   role := "side" u8 | "side:" tag | "cap-start" | "cap-end" | "top"
//         | "bottom" | "wall"
// ---------------------------------------------------------------------------

/// Percent-escapes the characters the grammar reserves.
fn escape_component(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '%' => out.push_str("%25"),
            ':' => out.push_str("%3A"),
            '(' => out.push_str("%28"),
            ')' => out.push_str("%29"),
            _ => out.push(c),
        }
    }
    out
}

fn unescape_component(s: &str) -> Option<String> {
    let mut out = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if c == '%' {
            let hex: String = chars.by_ref().take(2).collect();
            let code = u8::from_str_radix(&hex, 16).ok()?;
            out.push(code as char);
        } else {
            out.push(c);
        }
    }
    Some(out)
}

fn escape_op(op: &OpId) -> String {
    escape_component(op.as_str())
}

fn unescape_op(s: &str) -> Option<OpId> {
    unescape_component(s).map(OpId::new)
}

impl fmt::Display for FaceRole {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Side(k) => write!(f, "side{k}"),
            Self::Tagged(tag) => write!(f, "side:{}", escape_component(tag.as_str())),
            Self::CapStart => f.write_str("cap-start"),
            Self::CapEnd => f.write_str("cap-end"),
            Self::Top => f.write_str("top"),
            Self::Bottom => f.write_str("bottom"),
            Self::Wall => f.write_str("wall"),
        }
    }
}

impl std::str::FromStr for FaceRole {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, ()> {
        match s {
            "cap-start" => Ok(Self::CapStart),
            "cap-end" => Ok(Self::CapEnd),
            "top" => Ok(Self::Top),
            "bottom" => Ok(Self::Bottom),
            "wall" => Ok(Self::Wall),
            _ => {
                // "side:" (tagged) must be checked before the bare "side"
                // prefix of the positional form.
                if let Some(tag) = s.strip_prefix("side:") {
                    return unescape_component(tag)
                        .map(|t| Self::Tagged(SegmentTag::new(t)))
                        .ok_or(());
                }
                s.strip_prefix("side")
                    .and_then(|k| k.parse::<u8>().ok())
                    .map(Self::Side)
                    .ok_or(())
            }
        }
    }
}

impl fmt::Display for FaceName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Created { op, role } => write!(f, "created:{}:{role}", escape_op(op)),
            Self::Band {
                op,
                tool_face,
                loop_index,
            } => write!(f, "band:{}:{loop_index}:({tool_face})", escape_op(op)),
            Self::Floor { op, cap } => write!(f, "floor:{}:({cap})", escape_op(op)),
        }
    }
}

impl fmt::Display for EdgeName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Created { op, role } => {
                let end = match role {
                    EdgeRole::RingStart => "start",
                    EdgeRole::RingEnd => "end",
                };
                write!(f, "ring:{}:{end}", escape_op(op))
            }
            Self::CutRim {
                op,
                target,
                loop_index,
            } => write!(f, "rim:{}:{loop_index}:({target})", escape_op(op)),
        }
    }
}

/// Splits `"op:rest"` at the first unescaped `:`.
fn split_op(s: &str) -> Option<(OpId, &str)> {
    let idx = s.find(':')?;
    Some((unescape_op(&s[..idx])?, &s[idx + 1..]))
}

/// Parses `"index:(inner)"` into the index and the inner text.
fn split_indexed_inner(s: &str) -> Option<(u32, &str)> {
    let idx = s.find(':')?;
    let n = s[..idx].parse::<u32>().ok()?;
    let inner = s[idx + 1..].strip_prefix('(')?.strip_suffix(')')?;
    Some((n, inner))
}

impl std::str::FromStr for FaceName {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, ()> {
        if let Some(rest) = s.strip_prefix("created:") {
            let (op, role) = split_op(rest).ok_or(())?;
            return Ok(Self::Created {
                op,
                role: role.parse()?,
            });
        }
        if let Some(rest) = s.strip_prefix("band:") {
            let (op, rest) = split_op(rest).ok_or(())?;
            let (loop_index, inner) = split_indexed_inner(rest).ok_or(())?;
            return Ok(Self::Band {
                op,
                tool_face: Box::new(inner.parse()?),
                loop_index,
            });
        }
        if let Some(rest) = s.strip_prefix("floor:") {
            let (op, rest) = split_op(rest).ok_or(())?;
            let inner = rest.strip_prefix('(').and_then(|r| r.strip_suffix(')'));
            return Ok(Self::Floor {
                op,
                cap: Box::new(inner.ok_or(())?.parse()?),
            });
        }
        Err(())
    }
}

impl std::str::FromStr for EdgeName {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, ()> {
        if let Some(rest) = s.strip_prefix("ring:") {
            let (op, end) = split_op(rest).ok_or(())?;
            let role = match end {
                "start" => EdgeRole::RingStart,
                "end" => EdgeRole::RingEnd,
                _ => return Err(()),
            };
            return Ok(Self::Created { op, role });
        }
        if let Some(rest) = s.strip_prefix("rim:") {
            let (op, rest) = split_op(rest).ok_or(())?;
            let (loop_index, inner) = split_indexed_inner(rest).ok_or(())?;
            return Ok(Self::CutRim {
                op,
                target: Box::new(inner.parse()?),
                loop_index,
            });
        }
        Err(())
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::geometry::surface::Plane;
    use crate::math::{Point3, Vector3};
    use crate::topology::{FaceData, FaceSurface, TopologyStore, WireId};

    /// Adds a minimal placeholder face (registry tests never dereference the
    /// wire, so a null wire id is fine).
    fn dummy_face(store: &mut TopologyStore) -> FaceId {
        let plane = Plane::new(Point3::origin(), Vector3::z(), Vector3::x()).unwrap();
        store.add_face(FaceData {
            surface: FaceSurface::Plane(plane),
            outer_wire: WireId::default(),
            inner_wires: vec![],
            same_sense: true,
            trim: None,
            pcurves: Vec::new(),
        })
    }

    fn outer(op: &str) -> FaceName {
        FaceName::Created {
            op: OpId::new(op),
            role: FaceRole::Side(1),
        }
    }

    #[test]
    fn bind_resolves_both_directions() {
        let mut store = TopologyStore::new();
        let face = dummy_face(&mut store);
        let mut reg = NameRegistry::default();
        reg.bind_face(face, outer("wall1"));
        assert_eq!(reg.face(&outer("wall1")), Some(face));
        assert_eq!(reg.name_of_face(face), Some(&outer("wall1")));
        assert_eq!(reg.face(&outer("wall2")), None);
    }

    #[test]
    fn rebinding_a_name_moves_it_off_the_old_face() {
        let mut store = TopologyStore::new();
        let old = dummy_face(&mut store);
        let new = dummy_face(&mut store);
        let mut reg = NameRegistry::default();
        reg.bind_face(old, outer("wall1"));
        reg.bind_face(new, outer("wall1"));
        assert_eq!(reg.face(&outer("wall1")), Some(new));
        assert_eq!(reg.name_of_face(old), None, "old holder must be unbound");
    }

    #[test]
    fn canonical_string_round_trips() {
        let cases: Vec<FaceName> = vec![
            FaceName::Created {
                op: OpId::new("wall1"),
                role: FaceRole::Side(1),
            },
            FaceName::Created {
                op: OpId::new("node:with(specials)%"),
                role: FaceRole::CapEnd,
            },
            FaceName::Band {
                op: OpId::new("cut1"),
                tool_face: Box::new(FaceName::Created {
                    op: OpId::new("win1"),
                    role: FaceRole::Side(0),
                }),
                loop_index: 0,
            },
            FaceName::Floor {
                op: OpId::new("cut1"),
                cap: Box::new(FaceName::Created {
                    op: OpId::new("win1"),
                    role: FaceRole::CapEnd,
                }),
            },
            // Tagged side face (F5): opaque caller tag, plain.
            FaceName::Created {
                op: OpId::new("wall1"),
                role: FaceRole::Tagged(SegmentTag::new("centerline-3/outer")),
            },
            // Tagged side face with every reserved character in the tag.
            FaceName::Created {
                op: OpId::new("op:with(specials)%"),
                role: FaceRole::Tagged(SegmentTag::new("tag:with(specials)%25")),
            },
            // Tagged side face nested inside a band name (paren escaping).
            FaceName::Band {
                op: OpId::new("cut1"),
                tool_face: Box::new(FaceName::Created {
                    op: OpId::new("wall1"),
                    role: FaceRole::Tagged(SegmentTag::new("seg(0)")),
                }),
                loop_index: 2,
            },
            // Nested: a band cut into a band (future splits compose too).
            FaceName::Band {
                op: OpId::new("cut2"),
                tool_face: Box::new(FaceName::Band {
                    op: OpId::new("cut1"),
                    tool_face: Box::new(FaceName::Created {
                        op: OpId::new("win1"),
                        role: FaceRole::Side(0),
                    }),
                    loop_index: 0,
                }),
                loop_index: 3,
            },
        ];
        for name in cases {
            let text = name.to_string();
            let parsed: FaceName = text.parse().unwrap_or_else(|()| panic!("parse {text}"));
            assert_eq!(parsed, name, "round trip failed for {text}");
        }

        let edges: Vec<EdgeName> = vec![
            EdgeName::Created {
                op: OpId::new("tube1"),
                role: EdgeRole::RingStart,
            },
            EdgeName::CutRim {
                op: OpId::new("cut1"),
                target: Box::new(FaceName::Created {
                    op: OpId::new("slab1"),
                    role: FaceRole::Top,
                }),
                loop_index: 1,
            },
        ];
        for name in edges {
            let text = name.to_string();
            let parsed: EdgeName = text.parse().unwrap_or_else(|()| panic!("parse {text}"));
            assert_eq!(parsed, name, "round trip failed for {text}");
        }

        assert!("garbage".parse::<FaceName>().is_err());
        assert!("band:cut1:zero:(created:a:top)"
            .parse::<FaceName>()
            .is_err());
    }

    /// The positional `side<k>` and tagged `side:<tag>` forms stay distinct:
    /// a tag that is itself a digit string must not collapse into `Side(k)`.
    #[test]
    fn positional_and_tagged_side_forms_are_distinct() {
        let positional: FaceRole = "side0".parse().unwrap();
        assert_eq!(positional, FaceRole::Side(0));

        let tagged: FaceRole = "side:0".parse().unwrap();
        assert_eq!(tagged, FaceRole::Tagged(SegmentTag::new("0")));
        assert_ne!(positional, tagged);

        // Round trip preserves the distinction.
        assert_eq!(FaceRole::Side(0).to_string(), "side0");
        assert_eq!(FaceRole::Tagged(SegmentTag::new("0")).to_string(), "side:0");
    }

    #[test]
    fn transfer_moves_the_name_to_the_copy() {
        let mut store = TopologyStore::new();
        let original = dummy_face(&mut store);
        let copy = dummy_face(&mut store);
        let mut reg = NameRegistry::default();
        reg.bind_face(original, outer("wall1"));
        reg.transfer_face(original, copy);
        assert_eq!(reg.face(&outer("wall1")), Some(copy));
        assert_eq!(reg.name_of_face(original), None);
        // Transferring from an unnamed face is a no-op.
        let unrelated = dummy_face(&mut store);
        reg.transfer_face(unrelated, original);
        assert_eq!(reg.face(&outer("wall1")), Some(copy));
    }
}
