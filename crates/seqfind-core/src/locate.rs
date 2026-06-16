//! Regex-based motif locator -SeqKit's `seqkit locate` with `--pattern` / `--pattern-file`.
//!
//! Unlike [`crate::matcher`] which expects fixed strings with `?` wildcards,
//! this module accepts full PCRE-lite regex syntax (via the `regex` crate)
//! and supports IUPAC ambiguity codes via [`iupac_to_regex`].

use crate::matcher;
use crate::seq::SeqRecord;
use anyhow::{Context, Result};
use regex::bytes::Regex;
use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub struct LocateHit {
    pub record: String,
    pub pattern_name: String,
    pub pattern: String,
    pub start: usize,
    pub end: usize,
    pub strand: String,
    pub matched: String,
}

#[derive(Debug, Clone)]
pub struct LocatePattern {
    pub name: String,
    pub raw: String,
    pub regex: Regex,
}

impl LocatePattern {
    pub fn new(name: impl Into<String>, raw: impl Into<String>) -> Result<Self> {
        let raw = raw.into();
        let name = name.into();
        let regex = Regex::new(&raw)
            .with_context(|| format!("invalid regex for pattern '{name}': {raw}"))?;
        Ok(Self { name, raw, regex })
    }

    pub fn from_iupac(name: impl Into<String>, iupac: impl Into<String>) -> Result<Self> {
        let raw = iupac.into();
        let pattern = iupac_to_regex(&raw);
        let regex = Regex::new(&pattern)
            .with_context(|| format!("failed to compile IUPAC pattern '{raw}'"))?;
        Ok(Self {
            name: name.into(),
            raw,
            regex,
        })
    }
}

/// Translate IUPAC ambiguity codes (RYSWKM…BDHVN) into a character-class regex.
/// Useful for non-regex callers who want SeqKit-style "ambiguous nucleotide" matching.
pub fn iupac_to_regex(s: &str) -> String {
    let mut out = String::with_capacity(s.len() * 4);
    for ch in s.chars() {
        match ch.to_ascii_uppercase() {
            'A' => out.push('A'),
            'C' => out.push('C'),
            'G' => out.push('G'),
            'T' => out.push('T'),
            'U' => out.push('T'),
            'R' => out.push_str("[AG]"),
            'Y' => out.push_str("[CT]"),
            'S' => out.push_str("[GC]"),
            'W' => out.push_str("[AT]"),
            'K' => out.push_str("[GT]"),
            'M' => out.push_str("[AC]"),
            'B' => out.push_str("[CGT]"),
            'D' => out.push_str("[AGT]"),
            'H' => out.push_str("[ACT]"),
            'V' => out.push_str("[ACG]"),
            'N' => out.push_str("[ACGT]"),
            other => out.push(other), // pass through regex metacharacters
        }
    }
    out
}

pub fn locate_in_record(
    record: &SeqRecord,
    patterns: &[LocatePattern],
    reverse_complement: bool,
) -> Vec<LocateHit> {
    let mut hits = Vec::new();
    let seq = record.sequence.as_bytes();
    for p in patterns {
        for m in p.regex.find_iter(seq) {
            hits.push(LocateHit {
                record: record.name.clone(),
                pattern_name: p.name.clone(),
                pattern: p.raw.clone(),
                start: m.start(),
                end: m.end(),
                strand: "+".to_string(),
                matched: String::from_utf8_lossy(m.as_bytes()).into_owned(),
            });
        }
    }
    if reverse_complement {
        let rc = matcher::reverse_complement(&record.sequence);
        let rc_bytes = rc.as_bytes();
        for p in patterns {
            for m in p.regex.find_iter(rc_bytes) {
                // Translate rc-frame coordinates back to forward strand
                // (same convention as matcher: forward index of the RC occurrence).
                let forward_start = record.sequence.len() - m.end();
                hits.push(LocateHit {
                    record: record.name.clone(),
                    pattern_name: p.name.clone(),
                    pattern: p.raw.clone(),
                    start: forward_start,
                    end: forward_start + (m.end() - m.start()),
                    strand: "-".to_string(),
                    matched: String::from_utf8_lossy(m.as_bytes()).into_owned(),
                });
            }
        }
    }
    hits
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(name: &str, seq: &str) -> SeqRecord {
        SeqRecord {
            name: name.into(),
            sequence: seq.into(),
            quality: None,
        }
    }

    #[test]
    fn locates_forward_regex() {
        let r = rec("c", "ATGCGTAAATGCAT");
        let p = LocatePattern::new("atg", "ATG").unwrap();
        let hits = locate_in_record(&r, &[p], false);
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0].start, 0);
        assert_eq!(hits[1].start, 8);
    }

    #[test]
    fn locates_iupac() {
        let r = rec("c", "ATGCNGTAATG");
        // R = [AG]
        let p = LocatePattern::from_iupac("ambig", "RTG").unwrap();
        let hits = locate_in_record(&r, &[p], false);
        // matches at "ATG" (positions 0, 8) and "GTA" no -RTG needs R then T then G
        // ATG matches R=A T G ✓ at 0 and 8
        assert_eq!(hits.len(), 2);
    }

    #[test]
    fn locates_reverse_complement() {
        let r = rec("c", "AAATGC");
        // search for "GCAT" -reverse complement of "ATGC" appears at forward 2
        let p = LocatePattern::new("gcat", "GCAT").unwrap();
        let hits = locate_in_record(&r, &[p], true);
        // Forward "GCAT" doesn't appear, but RC pattern "ATGC" found at fwd index 2.
        let neg: Vec<_> = hits.iter().filter(|h| h.strand == "-").collect();
        assert_eq!(neg.len(), 1);
        assert_eq!(neg[0].start, 2);
    }
}
