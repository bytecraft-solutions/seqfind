//! Identifier-level edits -SeqKit's `replace`, `rename`, `concat`, `restart`, `mutate`.

use crate::seq::SeqRecord;
use anyhow::{Context, Result, anyhow, bail};
use regex::Regex;
use std::collections::HashMap;

/// Replace occurrences of `pattern` in record IDs with `replacement` (regex semantics).
pub fn replace_id(
    records: &[SeqRecord],
    pattern: &str,
    replacement: &str,
) -> Result<Vec<SeqRecord>> {
    let re = Regex::new(pattern).with_context(|| format!("invalid regex pattern: {pattern}"))?;
    Ok(records
        .iter()
        .map(|r| SeqRecord {
            name: re.replace_all(&r.name, replacement).into_owned(),
            sequence: r.sequence.clone(),
            quality: r.quality.clone(),
        })
        .collect())
}

/// Replace occurrences of `pattern` in record sequences. Operates on the
/// uppercase sequence text (sequences are already upper-cased on read).
pub fn replace_seq(
    records: &[SeqRecord],
    pattern: &str,
    replacement: &str,
) -> Result<Vec<SeqRecord>> {
    let re = Regex::new(pattern).with_context(|| format!("invalid regex pattern: {pattern}"))?;
    Ok(records
        .iter()
        .map(|r| SeqRecord {
            name: r.name.clone(),
            sequence: re.replace_all(&r.sequence, replacement).into_owned(),
            quality: r.quality.clone(),
        })
        .collect())
}

/// Disambiguate duplicate IDs by appending `_2`, `_3`, … to repeats.
/// Returns the new record list plus the number of records renamed.
pub fn rename_duplicates(records: &[SeqRecord]) -> (Vec<SeqRecord>, usize) {
    let mut counts: HashMap<String, usize> = HashMap::new();
    let mut renamed = 0usize;
    let out = records
        .iter()
        .map(|r| {
            let entry = counts.entry(r.name.clone()).or_insert(0);
            *entry += 1;
            let name = if *entry == 1 {
                r.name.clone()
            } else {
                renamed += 1;
                format!("{}_{}", r.name, *entry)
            };
            SeqRecord {
                name,
                sequence: r.sequence.clone(),
                quality: r.quality.clone(),
            }
        })
        .collect();
    (out, renamed)
}

/// How to align records from multiple files for [`concat`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConcatBy {
    /// Pair records by record ID across files. Records with no counterpart in
    /// every file are dropped (warning emitted by caller).
    Id,
    /// Pair records by index order across files. Stops at the shortest file.
    Order,
}

/// Concatenate corresponding records from multiple `files` into one sequence per
/// record. Names are taken from the first file. Returns a single record per
/// pairing.
pub fn concat(files: &[Vec<SeqRecord>], by: ConcatBy) -> Result<Vec<SeqRecord>> {
    if files.is_empty() {
        return Ok(Vec::new());
    }
    match by {
        ConcatBy::Order => {
            let count = files.iter().map(|f| f.len()).min().unwrap_or(0);
            let mut out = Vec::with_capacity(count);
            for i in 0..count {
                let first = &files[0][i];
                let mut seq =
                    String::with_capacity(files.iter().map(|f| f[i].sequence.len()).sum());
                for f in files {
                    seq.push_str(&f[i].sequence);
                }
                out.push(SeqRecord {
                    name: first.name.clone(),
                    sequence: seq,
                    quality: None,
                });
            }
            Ok(out)
        }
        ConcatBy::Id => {
            let mut out = Vec::with_capacity(files[0].len());
            let by_id: Vec<HashMap<&str, &SeqRecord>> = files
                .iter()
                .map(|f| f.iter().map(|r| (r.name.as_str(), r)).collect())
                .collect();
            for first in &files[0] {
                let mut seq = String::with_capacity(first.sequence.len());
                let mut all_present = true;
                for map in &by_id {
                    match map.get(first.name.as_str()) {
                        Some(r) => seq.push_str(&r.sequence),
                        None => {
                            all_present = false;
                            break;
                        }
                    }
                }
                if all_present {
                    out.push(SeqRecord {
                        name: first.name.clone(),
                        sequence: seq,
                        quality: None,
                    });
                }
            }
            Ok(out)
        }
    }
}

/// Rotate a circular sequence so it starts at 1-based `position`. `position`
/// must be in `[1, len]`. Quality bytes are rotated in lockstep.
pub fn restart(records: &[SeqRecord], position: usize) -> Result<Vec<SeqRecord>> {
    if position == 0 {
        bail!("restart position must be 1-based (got 0)");
    }
    let mut out = Vec::with_capacity(records.len());
    for r in records {
        let n = r.sequence.len();
        if n == 0 {
            out.push(r.clone());
            continue;
        }
        if position > n {
            bail!(
                "restart position {position} exceeds length {n} for record {}",
                r.name
            );
        }
        let cut = position - 1;
        let mut seq = String::with_capacity(n);
        seq.push_str(&r.sequence[cut..]);
        seq.push_str(&r.sequence[..cut]);
        let quality = r.quality.as_ref().map(|q| {
            let mut v = Vec::with_capacity(q.len());
            v.extend_from_slice(&q[cut..]);
            v.extend_from_slice(&q[..cut]);
            v
        });
        out.push(SeqRecord {
            name: r.name.clone(),
            sequence: seq,
            quality,
        });
    }
    Ok(out)
}

/// A single mutation applied to a sequence. Coordinates are 1-based inclusive.
#[derive(Debug, Clone)]
pub enum Mutation {
    /// `pos:BASE` -replace the base at `pos` with `base`.
    Point { pos: usize, base: char },
    /// `start:end` -delete the inclusive range (1-based).
    Deletion { start: usize, end: usize },
    /// `pos:BASES` -insert `bases` before `pos`. Use `pos = len + 1` to append.
    Insertion { pos: usize, bases: String },
}

impl Mutation {
    /// Parse `5:G`, `10:20`, `5:ATGC` based on which `kind` is requested.
    pub fn parse_point(s: &str) -> Result<Self> {
        let (a, b) = split_colon(s)?;
        let pos: usize = a
            .parse()
            .map_err(|e| anyhow!("invalid point pos '{a}': {e}"))?;
        if b.chars().count() != 1 {
            bail!("point mutation expects a single base, got '{b}'");
        }
        Ok(Self::Point {
            pos,
            base: b.chars().next().unwrap().to_ascii_uppercase(),
        })
    }

    pub fn parse_deletion(s: &str) -> Result<Self> {
        let (a, b) = split_colon(s)?;
        let start: usize = a
            .parse()
            .map_err(|e| anyhow!("invalid deletion start '{a}': {e}"))?;
        let end: usize = b
            .parse()
            .map_err(|e| anyhow!("invalid deletion end '{b}': {e}"))?;
        if start == 0 || end < start {
            bail!("deletion start:end must be 1-based with start<=end, got '{s}'");
        }
        Ok(Self::Deletion { start, end })
    }

    pub fn parse_insertion(s: &str) -> Result<Self> {
        let (a, b) = split_colon(s)?;
        let pos: usize = a
            .parse()
            .map_err(|e| anyhow!("invalid insertion pos '{a}': {e}"))?;
        if pos == 0 {
            bail!("insertion position must be >=1");
        }
        Ok(Self::Insertion {
            pos,
            bases: b.to_ascii_uppercase(),
        })
    }

    fn sort_key(&self) -> usize {
        match self {
            Mutation::Point { pos, .. } => *pos,
            Mutation::Deletion { start, .. } => *start,
            Mutation::Insertion { pos, .. } => *pos,
        }
    }
}

fn split_colon(s: &str) -> Result<(&str, &str)> {
    s.split_once(':')
        .ok_or_else(|| anyhow!("expected START:VALUE, got '{s}'"))
}

/// Apply a list of mutations to every record. Mutations operate on the original
/// coordinates: we sort by descending position and apply in that order so the
/// indices for earlier mutations stay valid. Quality bytes are kept aligned
/// (insertions inherit `I`, Phred 40).
pub fn mutate(records: &[SeqRecord], mutations: &[Mutation]) -> Result<Vec<SeqRecord>> {
    let mut ordered: Vec<Mutation> = mutations.to_vec();
    ordered.sort_by_key(|m| std::cmp::Reverse(m.sort_key()));

    let mut out = Vec::with_capacity(records.len());
    for r in records {
        let mut seq = r.sequence.as_bytes().to_vec();
        let mut qual = r.quality.clone();
        let n = seq.len();

        for m in &ordered {
            match m {
                Mutation::Point { pos, base } => {
                    if *pos == 0 || *pos > n {
                        bail!("point pos {pos} out of range 1..={n} for record {}", r.name);
                    }
                    seq[*pos - 1] = *base as u8;
                }
                Mutation::Deletion { start, end } => {
                    if *start == 0 || *end > n {
                        bail!(
                            "deletion {start}:{end} out of range 1..={n} for record {}",
                            r.name
                        );
                    }
                    seq.drain(*start - 1..*end);
                    if let Some(q) = qual.as_mut() {
                        q.drain(*start - 1..*end);
                    }
                }
                Mutation::Insertion { pos, bases } => {
                    if *pos == 0 || *pos > n + 1 {
                        bail!(
                            "insertion pos {pos} out of range 1..={} for record {}",
                            n + 1,
                            r.name
                        );
                    }
                    let idx = *pos - 1;
                    seq.splice(idx..idx, bases.bytes());
                    if let Some(q) = qual.as_mut() {
                        let pad = vec![b'I'; bases.len()];
                        q.splice(idx..idx, pad);
                    }
                }
            }
        }
        let new_seq = String::from_utf8(seq)
            .with_context(|| format!("mutated sequence not UTF-8 for {}", r.name))?;
        out.push(SeqRecord {
            name: r.name.clone(),
            sequence: new_seq,
            quality: qual,
        });
    }
    Ok(out)
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
    fn replace_id_simple() {
        let out = replace_id(&[rec("chr1", "A"), rec("chr2", "T")], r"chr", "scaffold_").unwrap();
        assert_eq!(out[0].name, "scaffold_1");
        assert_eq!(out[1].name, "scaffold_2");
    }

    #[test]
    fn replace_id_captures() {
        let out = replace_id(&[rec("gene_42_v1", "A")], r"gene_(\d+)_v\d+", "g$1").unwrap();
        assert_eq!(out[0].name, "g42");
    }

    #[test]
    fn replace_seq_basic() {
        let out = replace_seq(&[rec("a", "ATGN")], r"N", "A").unwrap();
        assert_eq!(out[0].sequence, "ATGA");
    }

    #[test]
    fn concat_by_order_joins_sequences() {
        let a = vec![rec("x", "AAA"), rec("y", "CCC")];
        let b = vec![rec("x", "TTT"), rec("y", "GGG")];
        let out = concat(&[a, b], ConcatBy::Order).unwrap();
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].name, "x");
        assert_eq!(out[0].sequence, "AAATTT");
        assert_eq!(out[1].sequence, "CCCGGG");
    }

    #[test]
    fn concat_by_id_drops_unmatched() {
        let a = vec![rec("x", "AA"), rec("y", "CC")];
        let b = vec![rec("y", "GG"), rec("x", "TT")];
        let out = concat(&[a, b], ConcatBy::Id).unwrap();
        assert_eq!(out[0].name, "x");
        assert_eq!(out[0].sequence, "AATT");
        assert_eq!(out[1].sequence, "CCGG");
    }

    #[test]
    fn restart_rotates_circular() {
        let out = restart(&[rec("c", "ABCDEF")], 3).unwrap();
        assert_eq!(out[0].sequence, "CDEFAB");
    }

    #[test]
    fn restart_position_out_of_range_errors() {
        assert!(restart(&[rec("c", "AAA")], 10).is_err());
    }

    #[test]
    fn mutate_point_substitution() {
        let m = Mutation::parse_point("2:G").unwrap();
        let out = mutate(&[rec("r", "AAAA")], &[m]).unwrap();
        assert_eq!(out[0].sequence, "AGAA");
    }

    #[test]
    fn mutate_deletion_and_insertion_compose() {
        // original "AAAGCCC"; insertion at 4 "TT", deletion 1:2
        let muts = vec![
            Mutation::parse_insertion("4:TT").unwrap(),
            Mutation::parse_deletion("1:2").unwrap(),
        ];
        let out = mutate(&[rec("r", "AAAGCCC")], &muts).unwrap();
        // Apply insertion first (highest pos), then deletion: "AAATTGCCC" → del 1..2 → "ATTGCCC"
        assert_eq!(out[0].sequence, "ATTGCCC");
    }

    #[test]
    fn rename_duplicates_appends_suffix() {
        let (out, n) = rename_duplicates(&[rec("a", ""), rec("a", ""), rec("b", ""), rec("a", "")]);
        let names: Vec<_> = out.iter().map(|r| r.name.clone()).collect();
        assert_eq!(names, vec!["a", "a_2", "b", "a_3"]);
        assert_eq!(n, 2);
    }
}
