//! Region extraction -SeqKit's `subseq --region`. 1-based inclusive coordinates
//! with negative-index support (counted from the end), matching SeqKit's semantics.

use crate::matcher;
use crate::seq::SeqRecord;
use anyhow::{Result, bail};

/// SeqKit-style 1-based inclusive coordinates. Negative values count from the end:
/// `-1` is the last base, `-30:-1` is the last 30 bases.
#[derive(Debug, Clone, Copy)]
pub struct Region {
    pub start: i64,
    pub end: i64,
}

impl Region {
    /// Parse `"5:10"`, `"-30:-1"`, `"100:-50"`.
    pub fn parse(s: &str) -> Result<Self> {
        let (a, b) = s
            .split_once(':')
            .ok_or_else(|| anyhow::anyhow!("region must be START:END, got '{s}'"))?;
        let start: i64 = a
            .trim()
            .parse()
            .map_err(|e| anyhow::anyhow!("invalid region start '{a}': {e}"))?;
        let end: i64 = b
            .trim()
            .parse()
            .map_err(|e| anyhow::anyhow!("invalid region end '{b}': {e}"))?;
        Ok(Self { start, end })
    }

    /// Resolve 1-based inclusive coordinates to a `[start, end)` byte slice for
    /// a sequence of length `len`. Returns `None` if the resolved range is empty
    /// after clamping (e.g. start > sequence length).
    pub fn resolve(&self, len: usize) -> Option<(usize, usize)> {
        let len = len as i64;
        if len == 0 {
            return None;
        }
        let resolve_one = |v: i64| -> i64 { if v < 0 { len + v + 1 } else { v } };
        let s = resolve_one(self.start).max(1);
        let e = resolve_one(self.end).min(len);
        if e < s {
            return None;
        }
        Some(((s - 1) as usize, e as usize))
    }
}

#[derive(Debug, Clone)]
pub struct SubseqOptions {
    pub region: Region,
    /// Take the region from the reverse-complement strand instead of forward.
    pub reverse_complement: bool,
    /// Bases of flanking context to include upstream of the region.
    pub up_stream: usize,
    /// Bases of flanking context to include downstream of the region.
    pub down_stream: usize,
}

pub fn subseq_record(record: &SeqRecord, opts: &SubseqOptions) -> Result<Option<SeqRecord>> {
    let Some((mut s, mut e)) = opts.region.resolve(record.sequence.len()) else {
        return Ok(None);
    };
    s = s.saturating_sub(opts.up_stream);
    e = (e + opts.down_stream).min(record.sequence.len());
    if e <= s {
        bail!("computed subseq window is empty for record {}", record.name);
    }
    let mut slice = record.sequence[s..e].to_string();
    let mut qual = record.quality.as_ref().map(|q| q[s..e].to_vec());
    if opts.reverse_complement {
        slice = matcher::reverse_complement(&slice);
        if let Some(q) = qual.as_mut() {
            q.reverse();
        }
    }
    let name = format!(
        "{}:{}-{}{}",
        record.name,
        s + 1,
        e,
        if opts.reverse_complement { ":-" } else { "" }
    );
    Ok(Some(SeqRecord {
        name,
        sequence: slice,
        quality: qual,
    }))
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
    fn region_positive_inclusive() {
        let r = Region::parse("2:5").unwrap();
        assert_eq!(r.resolve(10), Some((1, 5)));
    }

    #[test]
    fn region_negative_from_end() {
        let r = Region::parse("-3:-1").unwrap();
        // len 10 → -1 → 10, -3 → 8 → start=8, end=10 → 0-based [7, 10)
        assert_eq!(r.resolve(10), Some((7, 10)));
    }

    #[test]
    fn subseq_basic() {
        let out = subseq_record(
            &rec("c", "ACGTACGT"),
            &SubseqOptions {
                region: Region::parse("2:5").unwrap(),
                reverse_complement: false,
                up_stream: 0,
                down_stream: 0,
            },
        )
        .unwrap()
        .unwrap();
        assert_eq!(out.sequence, "CGTA");
        assert!(out.name.contains(":2-5"));
    }

    #[test]
    fn subseq_with_flanks_and_rc() {
        let out = subseq_record(
            &rec("c", "AAACGTAAA"),
            &SubseqOptions {
                region: Region::parse("4:6").unwrap(), // CGT
                reverse_complement: true,
                up_stream: 1,
                down_stream: 1,
            },
        )
        .unwrap()
        .unwrap();
        // window = ACGTA → RC = TACGT
        assert_eq!(out.sequence, "TACGT");
        assert!(out.name.ends_with(":-"));
    }
}
