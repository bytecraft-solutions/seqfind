//! In-silico PCR -extract sequences between a forward and reverse primer pair.
//!
//! Both primers are 5'→3'. The reverse primer's reverse-complement is searched
//! on the forward strand; the amplicon spans from the forward primer start to
//! the matched reverse-primer-RC end. IUPAC ambiguity codes are supported via
//! the same regex translation used by [`crate::locate`].

use crate::locate::iupac_to_regex;
use crate::matcher;
use crate::seq::SeqRecord;
use anyhow::{Context, Result, bail};
use regex::bytes::Regex;
use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub struct Amplicon {
    pub record: String,
    pub start: usize,
    pub end: usize,
    pub strand: String,
    pub length: usize,
    pub sequence: String,
    pub forward_primer: String,
    pub reverse_primer: String,
}

#[derive(Debug, Clone)]
pub struct AmpliconOptions {
    pub forward: String,
    pub reverse: String,
    /// Also scan the reverse-complement strand (swap roles of the two primers).
    pub both_strands: bool,
    /// Optional maximum amplicon length. Drop hits longer than this.
    pub max_length: Option<usize>,
}

pub fn find_amplicons(records: &[SeqRecord], opts: &AmpliconOptions) -> Result<Vec<Amplicon>> {
    if opts.forward.is_empty() || opts.reverse.is_empty() {
        bail!("both --forward and --reverse primers must be non-empty");
    }
    let forward_re = compile_iupac(&opts.forward).context("forward primer")?;
    let reverse_rc = matcher::reverse_complement(&opts.reverse.to_ascii_uppercase());
    let reverse_re = compile_iupac(&reverse_rc).context("reverse primer (RC)")?;

    let mut out = Vec::new();
    for r in records {
        scan_strand(
            &r.name,
            &r.sequence,
            &forward_re,
            &reverse_re,
            "+",
            opts,
            &mut out,
        );
        if opts.both_strands {
            let rc = matcher::reverse_complement(&r.sequence);
            let len = r.sequence.len();
            let mut tmp = Vec::new();
            scan_strand(&r.name, &rc, &forward_re, &reverse_re, "-", opts, &mut tmp);
            for mut a in tmp {
                // Translate rc-frame coordinates back to forward strand.
                let new_start = len - a.end;
                let new_end = len - a.start;
                a.start = new_start;
                a.end = new_end;
                out.push(a);
            }
        }
    }
    Ok(out)
}

fn scan_strand(
    name: &str,
    seq: &str,
    forward_re: &Regex,
    reverse_re: &Regex,
    strand: &str,
    opts: &AmpliconOptions,
    out: &mut Vec<Amplicon>,
) {
    let bytes = seq.as_bytes();
    for fm in forward_re.find_iter(bytes) {
        let fwd_start = fm.start();
        // Reverse-primer match must start at or after the forward primer match.
        for rm in reverse_re.find_iter(&bytes[fwd_start..]) {
            let rev_end = fwd_start + rm.end();
            let length = rev_end - fwd_start;
            if let Some(max) = opts.max_length
                && length > max
            {
                continue;
            }
            out.push(Amplicon {
                record: name.to_string(),
                start: fwd_start,
                end: rev_end,
                strand: strand.to_string(),
                length,
                sequence: seq[fwd_start..rev_end].to_string(),
                forward_primer: opts.forward.clone(),
                reverse_primer: opts.reverse.clone(),
            });
            break; // shortest amplicon per forward hit
        }
    }
}

fn compile_iupac(s: &str) -> Result<Regex> {
    let pattern = iupac_to_regex(s);
    Regex::new(&pattern).with_context(|| format!("failed to compile primer '{s}' as regex"))
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
    fn finds_simple_amplicon() {
        // forward primer "ATGC", reverse primer "GGGG" (RC = "CCCC")
        // sequence:        ATGCNNNNCCCC
        let recs = [rec("c", "ATGCNNNNCCCC")];
        let opts = AmpliconOptions {
            forward: "ATGC".into(),
            reverse: "GGGG".into(),
            both_strands: false,
            max_length: None,
        };
        let out = find_amplicons(&recs, &opts).unwrap();
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].sequence, "ATGCNNNNCCCC");
        assert_eq!(out[0].start, 0);
        assert_eq!(out[0].end, 12);
    }

    #[test]
    fn no_amplicon_when_primers_missing() {
        let recs = [rec("c", "AAAAAA")];
        let opts = AmpliconOptions {
            forward: "ATGC".into(),
            reverse: "GGGG".into(),
            both_strands: false,
            max_length: None,
        };
        assert!(find_amplicons(&recs, &opts).unwrap().is_empty());
    }

    #[test]
    fn max_length_filter() {
        let recs = [rec("c", "ATGCNNNNNNNNNNNNNNNNCCCC")];
        let opts = AmpliconOptions {
            forward: "ATGC".into(),
            reverse: "GGGG".into(),
            both_strands: false,
            max_length: Some(10),
        };
        assert!(find_amplicons(&recs, &opts).unwrap().is_empty());
    }
}
