use serde::Serialize;

use crate::matcher::reverse_complement;

/// One open reading frame.
///
/// Coordinates are always **forward-strand**, 0-based, half-open:
/// `start` is the index of the first base of the start codon, `end` is one
/// past the last base of the stop codon, so `end - start` is the ORF length
/// in bases (including start and stop codons).
#[derive(Debug, Clone, Serialize)]
pub struct Orf {
    pub record: String,
    pub start: usize,
    pub end: usize,
    pub strand: String,
    pub start_codon: String,
    pub stop_codon: String,
    pub length: usize,
}

const START_CODON: &[u8] = b"ATG";

fn is_stop(codon: &[u8]) -> bool {
    matches!(codon, b"TAA" | b"TAG" | b"TGA")
}

/// Scan one strand of `seq` in all three reading frames for ORFs of at least
/// `min_length` bases (including start and stop codons). `strand` is the
/// label to emit. `to_forward` converts a local index on `seq` into the
/// forward-strand index for the record (identity for `+`, mirror for `-`).
fn scan_strand<F>(
    record: &str,
    seq: &[u8],
    strand: &str,
    min_length: usize,
    to_forward: F,
    out: &mut Vec<Orf>,
) where
    F: Fn(usize) -> usize,
{
    for frame in 0..3 {
        let mut i = frame;
        let mut current_start: Option<usize> = None;

        while i + 3 <= seq.len() {
            let codon = &seq[i..i + 3];

            if current_start.is_none() {
                if codon == START_CODON {
                    current_start = Some(i);
                }
            } else if is_stop(codon) {
                let s = current_start.take().unwrap();
                let e = i + 3;
                if e - s >= min_length {
                    let (fwd_start, fwd_end) = if strand == "+" {
                        (to_forward(s), to_forward(e))
                    } else {
                        // Mirror: end on the forward strand corresponds to start
                        // on the RC strand, and vice versa.
                        (to_forward(e), to_forward(s))
                    };
                    out.push(Orf {
                        record: record.to_string(),
                        start: fwd_start,
                        end: fwd_end,
                        strand: strand.to_string(),
                        start_codon: String::from_utf8_lossy(&seq[s..s + 3]).into_owned(),
                        stop_codon: String::from_utf8_lossy(&seq[i..i + 3]).into_owned(),
                        length: e - s,
                    });
                }
            }

            i += 3;
        }
    }
}

pub fn find_orfs(record: &str, sequence: &str, min_length: usize, include_rc: bool) -> Vec<Orf> {
    let fwd = sequence.as_bytes();
    let mut orfs = Vec::new();

    scan_strand(record, fwd, "+", min_length, |i| i, &mut orfs);

    if include_rc {
        let rc = reverse_complement(sequence);
        let rc_bytes = rc.as_bytes().to_vec();
        let seq_len = fwd.len();
        scan_strand(
            record,
            &rc_bytes,
            "-",
            min_length,
            |i| seq_len - i,
            &mut orfs,
        );
    }

    orfs.sort_by(|a, b| {
        a.record
            .cmp(&b.record)
            .then(a.start.cmp(&b.start))
            .then(a.strand.cmp(&b.strand))
    });
    orfs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_simple_forward_orf() {
        // ATG CCC TAA → ORF length 9
        let orfs = find_orfs("seq1", "ATGCCCTAA", 9, false);
        assert_eq!(orfs.len(), 1);
        let o = &orfs[0];
        assert_eq!(o.start, 0);
        assert_eq!(o.end, 9);
        assert_eq!(o.strand, "+");
        assert_eq!(o.start_codon, "ATG");
        assert_eq!(o.stop_codon, "TAA");
    }

    #[test]
    fn respects_min_length() {
        let orfs = find_orfs("seq1", "ATGCCCTAA", 12, false);
        assert!(orfs.is_empty());
    }

    #[test]
    fn finds_orf_in_offset_frame() {
        // ATG starts at position 1, in frame 1
        let orfs = find_orfs("seq1", "GATGCCCTAA", 9, false);
        assert_eq!(orfs.len(), 1);
        assert_eq!(orfs[0].start, 1);
        assert_eq!(orfs[0].end, 10);
    }

    #[test]
    fn rc_orf_coordinates_are_forward_strand() {
        // Forward seq: TTAGGGCAT (RC = ATGCCCTAA → ORF on minus strand)
        let orfs = find_orfs("seq1", "TTAGGGCAT", 9, true);
        let minus: Vec<&Orf> = orfs.iter().filter(|o| o.strand == "-").collect();
        assert_eq!(minus.len(), 1);
        // RC ORF spans the full 9 bp → forward start=0, end=9
        assert_eq!(minus[0].start, 0);
        assert_eq!(minus[0].end, 9);
        assert_eq!(minus[0].start_codon, "ATG");
        assert_eq!(minus[0].stop_codon, "TAA");
    }

    #[test]
    fn skips_start_without_stop() {
        let orfs = find_orfs("seq1", "ATGCCCCCC", 6, false);
        assert!(orfs.is_empty());
    }

    #[test]
    fn finds_consecutive_orfs_in_same_frame() {
        // ATG-AAA-TAA-ATG-CCC-TAA → two ORFs
        let orfs = find_orfs("seq1", "ATGAAATAAATGCCCTAA", 9, false);
        assert_eq!(orfs.len(), 2);
        assert_eq!(orfs[0].start, 0);
        assert_eq!(orfs[1].start, 9);
    }
}
