//! Identifier-level edits -SeqKit's `replace` and `rename`.

use crate::seq::SeqRecord;
use anyhow::{Context, Result};
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
    fn rename_duplicates_appends_suffix() {
        let (out, n) = rename_duplicates(&[rec("a", ""), rec("a", ""), rec("b", ""), rec("a", "")]);
        let names: Vec<_> = out.iter().map(|r| r.name.clone()).collect();
        assert_eq!(names, vec!["a", "a_2", "b", "a_3"]);
        assert_eq!(n, 2);
    }
}
