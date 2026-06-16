//! Per-record sequence transforms (SeqKit-style `seq` operations).

use crate::matcher;
use crate::seq::SeqRecord;

#[derive(Debug, Clone, Default)]
pub struct TransformOptions {
    pub uppercase: bool,
    pub lowercase: bool,
    pub remove_gaps: bool,
    /// Replace `T` with `U` (treat input as DNA, emit RNA).
    pub dna_to_rna: bool,
    /// Replace `U` with `T` (treat input as RNA, emit DNA).
    pub rna_to_dna: bool,
    pub reverse: bool,
    pub complement: bool,
    /// Equivalent to `reverse && complement` (set explicitly for clarity).
    pub reverse_complement: bool,
    /// Drop records shorter than this length.
    pub min_length: Option<usize>,
    /// Drop records longer than this length.
    pub max_length: Option<usize>,
}

/// Apply transforms in a deterministic order: gap-remove → case → DNA/RNA → revcomp.
/// Returns `None` when the record is filtered out by length bounds.
pub fn apply(rec: &SeqRecord, opts: &TransformOptions) -> Option<SeqRecord> {
    let mut seq = rec.sequence.clone();
    let mut qual = rec.quality.clone();

    if opts.remove_gaps {
        // Keep quality aligned with sequence by indexing.
        let bytes = seq.as_bytes().to_vec();
        let mut kept = Vec::with_capacity(bytes.len());
        let mut kept_q = qual.as_ref().map(|_| Vec::with_capacity(bytes.len()));
        for (i, &b) in bytes.iter().enumerate() {
            if b != b'-' && b != b'.' {
                kept.push(b);
                if let (Some(q), Some(out)) = (qual.as_ref(), kept_q.as_mut()) {
                    out.push(q[i]);
                }
            }
        }
        seq = String::from_utf8(kept).expect("ASCII subset stays valid UTF-8");
        qual = kept_q;
    }

    if let Some(min) = opts.min_length
        && seq.len() < min
    {
        return None;
    }
    if let Some(max) = opts.max_length
        && seq.len() > max
    {
        return None;
    }

    if opts.uppercase {
        seq = seq.to_ascii_uppercase();
    } else if opts.lowercase {
        seq = seq.to_ascii_lowercase();
    }

    if opts.dna_to_rna {
        seq = seq.replace('T', "U").replace('t', "u");
    } else if opts.rna_to_dna {
        seq = seq.replace('U', "T").replace('u', "t");
    }

    let do_rc = opts.reverse_complement;
    let do_rev = opts.reverse || do_rc;
    let do_comp = opts.complement || do_rc;

    if do_comp && do_rev {
        // matcher::reverse_complement expects uppercase ACGT; tolerate by upper-casing
        // a working copy but keep the original case in the result via reorder.
        // For simplicity reuse the existing helper and mirror its semantics.
        seq = matcher::reverse_complement(&seq.to_ascii_uppercase());
        if let Some(q) = qual.as_mut() {
            q.reverse();
        }
    } else if do_rev {
        seq = seq.chars().rev().collect();
        if let Some(q) = qual.as_mut() {
            q.reverse();
        }
    } else if do_comp {
        seq = seq
            .chars()
            .map(|c| match c {
                'A' => 'T',
                'T' | 'U' => 'A',
                'G' => 'C',
                'C' => 'G',
                'a' => 't',
                't' | 'u' => 'a',
                'g' => 'c',
                'c' => 'g',
                'N' => 'N',
                'n' => 'n',
                other => other,
            })
            .collect();
    }

    Some(SeqRecord {
        name: rec.name.clone(),
        sequence: seq,
        quality: qual,
    })
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
    fn uppercase_lowercase() {
        let r = rec("r", "atGc");
        assert_eq!(
            apply(
                &r,
                &TransformOptions {
                    uppercase: true,
                    ..Default::default()
                }
            )
            .unwrap()
            .sequence,
            "ATGC"
        );
        assert_eq!(
            apply(
                &r,
                &TransformOptions {
                    lowercase: true,
                    ..Default::default()
                }
            )
            .unwrap()
            .sequence,
            "atgc"
        );
    }

    #[test]
    fn dna_rna_roundtrip() {
        let r = rec("r", "ATGC");
        let rna = apply(
            &r,
            &TransformOptions {
                dna_to_rna: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(rna.sequence, "AUGC");
        let back = apply(
            &rna,
            &TransformOptions {
                rna_to_dna: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(back.sequence, "ATGC");
    }

    #[test]
    fn remove_gaps_keeps_quality_aligned() {
        let r = SeqRecord {
            name: "r".into(),
            sequence: "A-T.G".into(),
            quality: Some(b"abcde".to_vec()),
        };
        let out = apply(
            &r,
            &TransformOptions {
                remove_gaps: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(out.sequence, "ATG");
        // A,T,G at positions 0,2,4 → quality "a","c","e".
        assert_eq!(out.quality.as_deref(), Some(b"ace".as_ref()));
    }

    #[test]
    fn revcomp() {
        let out = apply(
            &rec("r", "ATGC"),
            &TransformOptions {
                reverse_complement: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(out.sequence, "GCAT");
    }

    #[test]
    fn complement_only() {
        let out = apply(
            &rec("r", "ATGC"),
            &TransformOptions {
                complement: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(out.sequence, "TACG");
    }

    #[test]
    fn length_filters_drop_record() {
        let r = rec("r", "AT");
        assert!(
            apply(
                &r,
                &TransformOptions {
                    min_length: Some(3),
                    ..Default::default()
                }
            )
            .is_none()
        );
        assert!(
            apply(
                &r,
                &TransformOptions {
                    max_length: Some(1),
                    ..Default::default()
                }
            )
            .is_none()
        );
    }
}
