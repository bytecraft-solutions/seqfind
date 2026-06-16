//! Set/ordering operations on collections of records: head/range/sample/rmdup/
//! sort/shuffle/split -SeqKit's `head`/`range`/`sample`/`rmdup`/`sort`/`shuffle`/`split`.

use crate::seq::SeqRecord;
use anyhow::{Result, bail};
use rand::SeedableRng;
use rand::prelude::*;
use rand_chacha::ChaCha20Rng;
use std::collections::HashSet;

pub fn head(records: &[SeqRecord], n: usize) -> Vec<SeqRecord> {
    records.iter().take(n).cloned().collect()
}

/// 1-based inclusive range, SeqKit style. `--range 5:10` → indices 5..=10.
pub fn range(records: &[SeqRecord], start: usize, end: usize) -> Result<Vec<SeqRecord>> {
    if start == 0 || end == 0 {
        bail!("range indices are 1-based; got start={start} end={end}");
    }
    if end < start {
        bail!("range end ({end}) is before start ({start})");
    }
    Ok(records
        .iter()
        .skip(start - 1)
        .take(end - start + 1)
        .cloned()
        .collect())
}

#[derive(Debug, Clone, Copy)]
pub enum SampleKind {
    /// Sample exactly N records uniformly without replacement.
    Number(usize),
    /// Keep each record with probability `p`.
    Proportion(f64),
}

pub fn sample(records: &[SeqRecord], kind: SampleKind, seed: u64) -> Vec<SeqRecord> {
    let mut rng = ChaCha20Rng::seed_from_u64(seed);
    match kind {
        SampleKind::Proportion(p) => records
            .iter()
            .filter(|_| rng.r#gen::<f64>() < p)
            .cloned()
            .collect(),
        SampleKind::Number(n) => {
            if n >= records.len() {
                return records.to_vec();
            }
            // Indices reservoir sample to preserve original order.
            let mut idx: Vec<usize> = (0..records.len()).collect();
            idx.shuffle(&mut rng);
            idx.truncate(n);
            idx.sort_unstable();
            idx.into_iter().map(|i| records[i].clone()).collect()
        }
    }
}

pub fn shuffle(records: &[SeqRecord], seed: u64) -> Vec<SeqRecord> {
    let mut out = records.to_vec();
    let mut rng = ChaCha20Rng::seed_from_u64(seed);
    out.shuffle(&mut rng);
    out
}

#[derive(Debug, Clone, Copy)]
pub enum RmdupKey {
    Id,
    Sequence,
}

pub fn rmdup(records: &[SeqRecord], key: RmdupKey) -> (Vec<SeqRecord>, usize) {
    let mut seen = HashSet::new();
    let mut out = Vec::with_capacity(records.len());
    let mut dropped = 0usize;
    for r in records {
        let k = match key {
            RmdupKey::Id => r.name.clone(),
            RmdupKey::Sequence => r.sequence.clone(),
        };
        if seen.insert(k) {
            out.push(r.clone());
        } else {
            dropped += 1;
        }
    }
    (out, dropped)
}

#[derive(Debug, Clone, Copy)]
pub enum SortKey {
    Id,
    Length,
    Sequence,
}

pub fn sort(records: &[SeqRecord], key: SortKey, reverse: bool) -> Vec<SeqRecord> {
    let mut out = records.to_vec();
    match key {
        SortKey::Id => out.sort_by(|a, b| a.name.cmp(&b.name)),
        SortKey::Sequence => out.sort_by(|a, b| a.sequence.cmp(&b.sequence)),
        SortKey::Length => out.sort_by_key(|r| r.sequence.len()),
    }
    if reverse {
        out.reverse();
    }
    out
}

/// Split into `n` near-equal chunks (round-robin to preserve uniform sizing).
pub fn split_into_parts(records: &[SeqRecord], n: usize) -> Vec<Vec<SeqRecord>> {
    if n == 0 {
        return vec![records.to_vec()];
    }
    let mut buckets: Vec<Vec<SeqRecord>> = (0..n).map(|_| Vec::new()).collect();
    let chunk_size = records.len().div_ceil(n).max(1);
    for (i, r) in records.iter().enumerate() {
        let bucket = (i / chunk_size).min(n - 1);
        buckets[bucket].push(r.clone());
    }
    buckets
}

/// Split into chunks of size `size`.
pub fn split_by_size(records: &[SeqRecord], size: usize) -> Vec<Vec<SeqRecord>> {
    if size == 0 {
        return vec![records.to_vec()];
    }
    records.chunks(size).map(|c| c.to_vec()).collect()
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

    fn make(n: usize) -> Vec<SeqRecord> {
        (0..n).map(|i| rec(&format!("r{i}"), "A")).collect()
    }

    #[test]
    fn head_truncates() {
        assert_eq!(head(&make(5), 2).len(), 2);
        assert_eq!(head(&make(5), 99).len(), 5);
    }

    #[test]
    fn range_inclusive_one_based() {
        let r = range(&make(5), 2, 4).unwrap();
        assert_eq!(
            r.iter().map(|x| x.name.clone()).collect::<Vec<_>>(),
            ["r1", "r2", "r3"]
        );
    }

    #[test]
    fn range_rejects_zero() {
        assert!(range(&make(5), 0, 2).is_err());
    }

    #[test]
    fn sample_seed_is_deterministic() {
        let a = sample(&make(20), SampleKind::Number(5), 42);
        let b = sample(&make(20), SampleKind::Number(5), 42);
        assert_eq!(
            a.iter().map(|x| x.name.clone()).collect::<Vec<_>>(),
            b.iter().map(|x| x.name.clone()).collect::<Vec<_>>()
        );
        assert_eq!(a.len(), 5);
    }

    #[test]
    fn sample_proportion_seed_is_deterministic() {
        let a = sample(&make(50), SampleKind::Proportion(0.5), 7);
        let b = sample(&make(50), SampleKind::Proportion(0.5), 7);
        assert_eq!(a.len(), b.len());
        assert_eq!(
            a.iter().map(|x| x.name.clone()).collect::<Vec<_>>(),
            b.iter().map(|x| x.name.clone()).collect::<Vec<_>>()
        );
    }

    #[test]
    fn shuffle_seed_is_deterministic() {
        let a = shuffle(&make(10), 11);
        let b = shuffle(&make(10), 11);
        assert_eq!(
            a.iter().map(|x| x.name.clone()).collect::<Vec<_>>(),
            b.iter().map(|x| x.name.clone()).collect::<Vec<_>>()
        );
    }

    #[test]
    fn rmdup_by_id_and_seq() {
        let recs = vec![
            rec("a", "ATG"),
            rec("a", "TTT"),
            rec("b", "ATG"),
            rec("c", "GGG"),
        ];
        let (by_id, dropped_id) = rmdup(&recs, RmdupKey::Id);
        assert_eq!(by_id.len(), 3);
        assert_eq!(dropped_id, 1);
        let (by_seq, dropped_seq) = rmdup(&recs, RmdupKey::Sequence);
        assert_eq!(by_seq.len(), 3);
        assert_eq!(dropped_seq, 1);
    }

    #[test]
    fn sort_by_length() {
        let recs = vec![rec("a", "AA"), rec("b", "A"), rec("c", "AAA")];
        let s = sort(&recs, SortKey::Length, false);
        assert_eq!(
            s.iter().map(|r| r.name.clone()).collect::<Vec<_>>(),
            ["b", "a", "c"]
        );
        let s = sort(&recs, SortKey::Length, true);
        assert_eq!(
            s.iter().map(|r| r.name.clone()).collect::<Vec<_>>(),
            ["c", "a", "b"]
        );
    }

    #[test]
    fn split_into_parts_balances() {
        let recs = make(10);
        let parts = split_into_parts(&recs, 3);
        assert_eq!(parts.len(), 3);
        assert_eq!(parts.iter().map(|p| p.len()).sum::<usize>(), 10);
    }

    #[test]
    fn split_by_size_chunks() {
        let recs = make(7);
        let parts = split_by_size(&recs, 3);
        assert_eq!(parts.len(), 3);
        assert_eq!(parts[0].len(), 3);
        assert_eq!(parts[2].len(), 1);
    }
}
