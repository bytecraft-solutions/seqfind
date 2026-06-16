//! Aggregate sequence stats -SeqKit-style `stats` summary (N50, sum/min/avg/max,
//! GC%, optional FASTQ Q20/Q30). Computed per-record stats stay in [`crate::output`].

use crate::seq::SeqRecord;
use serde::Serialize;

#[derive(Debug, Clone, Default, Serialize)]
pub struct StatsSummary {
    pub num_seqs: usize,
    pub sum_len: usize,
    pub min_len: usize,
    pub avg_len: f64,
    pub max_len: usize,
    pub n50: usize,
    pub gc_percent: f64,
    /// Total Ns across all records.
    pub n_count: usize,
    /// Percent bases at or above Phred 20 (FASTQ only). `None` for FASTA inputs.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub q20_percent: Option<f64>,
    /// Percent bases at or above Phred 30 (FASTQ only). `None` for FASTA inputs.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub q30_percent: Option<f64>,
}

/// Assume Phred+33 quality encoding (Sanger / Illumina 1.8+). That's the
/// modern default and the only encoding SeqKit treats as canonical too.
const PHRED_OFFSET: u8 = 33;

pub fn compute(records: &[SeqRecord]) -> StatsSummary {
    if records.is_empty() {
        return StatsSummary::default();
    }
    let mut lengths: Vec<usize> = records.iter().map(|r| r.sequence.len()).collect();
    let sum_len: usize = lengths.iter().sum();
    let min_len = *lengths.iter().min().unwrap_or(&0);
    let max_len = *lengths.iter().max().unwrap_or(&0);
    let avg_len = sum_len as f64 / records.len() as f64;

    lengths.sort_unstable_by(|a, b| b.cmp(a)); // descending
    let half = sum_len / 2;
    let mut acc = 0usize;
    let mut n50 = 0usize;
    for l in &lengths {
        acc += l;
        if acc >= half {
            n50 = *l;
            break;
        }
    }

    let mut gc_count = 0usize;
    let mut n_count = 0usize;
    for r in records {
        for ch in r.sequence.chars() {
            match ch {
                'G' | 'C' | 'g' | 'c' => gc_count += 1,
                'N' | 'n' => n_count += 1,
                _ => {}
            }
        }
    }
    let gc_percent = if sum_len == 0 {
        0.0
    } else {
        gc_count as f64 / sum_len as f64 * 100.0
    };

    let has_quality = records.iter().any(|r| r.quality.is_some());
    let (q20_percent, q30_percent) = if has_quality {
        let mut total_q = 0usize;
        let mut q20 = 0usize;
        let mut q30 = 0usize;
        for r in records {
            if let Some(q) = &r.quality {
                for &b in q {
                    let phred = b.saturating_sub(PHRED_OFFSET);
                    total_q += 1;
                    if phred >= 20 {
                        q20 += 1;
                    }
                    if phred >= 30 {
                        q30 += 1;
                    }
                }
            }
        }
        if total_q == 0 {
            (Some(0.0), Some(0.0))
        } else {
            (
                Some(q20 as f64 / total_q as f64 * 100.0),
                Some(q30 as f64 / total_q as f64 * 100.0),
            )
        }
    } else {
        (None, None)
    };

    StatsSummary {
        num_seqs: records.len(),
        sum_len,
        min_len,
        avg_len,
        max_len,
        n50,
        gc_percent,
        n_count,
        q20_percent,
        q30_percent,
    }
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
    fn empty_input() {
        let s = compute(&[]);
        assert_eq!(s.num_seqs, 0);
        assert_eq!(s.sum_len, 0);
        assert_eq!(s.n50, 0);
    }

    #[test]
    fn n50_basic() {
        // lengths 9, 5, 3 → sum 17, half 8 → cumulative descending: 9 (≥8), n50=9
        let recs = [rec("a", "AAAAAAAAA"), rec("b", "AAAAA"), rec("c", "AAA")];
        let s = compute(&recs);
        assert_eq!(s.num_seqs, 3);
        assert_eq!(s.sum_len, 17);
        assert_eq!(s.min_len, 3);
        assert_eq!(s.max_len, 9);
        assert_eq!(s.n50, 9);
    }

    #[test]
    fn n50_three_equal() {
        let recs = [rec("a", "AAAA"), rec("b", "AAAA"), rec("c", "AAAA")];
        let s = compute(&recs);
        // sum 12, half 6 → after 4 acc=4 (<6), after 8 acc=8 (≥6) → n50=4
        assert_eq!(s.n50, 4);
    }

    #[test]
    fn gc_and_n() {
        let recs = [rec("a", "ATGCNNN"), rec("b", "GGGG")];
        let s = compute(&recs);
        assert_eq!(s.sum_len, 11);
        assert_eq!(s.n_count, 3);
        // 2 (G+C from first) + 4 (GGGG) = 6 GC, 6/11 ≈ 54.5
        assert!((s.gc_percent - (6.0 / 11.0 * 100.0)).abs() < 1e-9);
    }

    #[test]
    fn q20_q30_phred33() {
        // "!" = Phred 0, "5" = Phred 20, "?" = Phred 30
        let recs = [SeqRecord {
            name: "r".into(),
            sequence: "AAAA".into(),
            quality: Some(b"!5?I".to_vec()),
        }];
        let s = compute(&recs);
        // 3 of 4 are ≥20, 2 of 4 are ≥30
        assert_eq!(s.q20_percent, Some(75.0));
        assert_eq!(s.q30_percent, Some(50.0));
    }

    #[test]
    fn fasta_has_no_quality_percent() {
        let s = compute(&[rec("a", "ATGC")]);
        assert!(s.q20_percent.is_none());
        assert!(s.q30_percent.is_none());
    }
}
