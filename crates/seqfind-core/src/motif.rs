//! Named motif presets -a thin lookup table over the `search` machinery.
//!
//! Names are case-insensitive. Keep the list small and conservative; this is
//! convenience sugar, not a curated motif database.

pub struct Motif {
    pub name: &'static str,
    pub pattern: &'static str,
    pub description: &'static str,
}

pub const PRESETS: &[Motif] = &[
    Motif {
        name: "EcoRI",
        pattern: "GAATTC",
        description: "EcoRI restriction site",
    },
    Motif {
        name: "BamHI",
        pattern: "GGATCC",
        description: "BamHI restriction site",
    },
    Motif {
        name: "HindIII",
        pattern: "AAGCTT",
        description: "HindIII restriction site",
    },
    Motif {
        name: "NotI",
        pattern: "GCGGCCGC",
        description: "NotI restriction site",
    },
    Motif {
        name: "TATA",
        pattern: "TATAAA",
        description: "TATA box (core)",
    },
    Motif {
        name: "Kozak",
        pattern: "GCCACCATGG",
        description: "Kozak consensus (vertebrate)",
    },
    Motif {
        name: "PolyA",
        pattern: "AATAAA",
        description: "Canonical polyadenylation signal",
    },
    Motif {
        name: "ShineDalgarno",
        pattern: "AGGAGG",
        description: "Shine-Dalgarno ribosome binding (E. coli)",
    },
];

pub fn lookup(name: &str) -> Option<&'static Motif> {
    PRESETS.iter().find(|m| m.name.eq_ignore_ascii_case(name))
}
