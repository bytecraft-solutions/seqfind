use bio::pattern_matching::horspool::Horspool;
use serde::Serialize;

/// One sequence-search hit.
///
/// Coordinates are forward-strand, 0-based, half-open: `start` is the index of
/// the first matching base on the forward strand, `end = start + matched.len()`.
#[derive(Debug, Clone, Serialize)]
pub struct Hit {
    pub record: String,
    pub start: usize,
    pub end: usize,
    pub strand: String,
    pub matched: String,
    pub pattern: String,
}

fn complement(base: char) -> char {
    match base {
        'A' => 'T',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
        'N' => 'N',
        other => other,
    }
}

pub fn reverse_complement(seq: &str) -> String {
    seq.chars().rev().map(complement).collect()
}

fn wildcard_match(pattern: &str, text: &str) -> bool {
    if pattern.len() != text.len() {
        return false;
    }
    pattern
        .chars()
        .zip(text.chars())
        .all(|(p, t)| p == '?' || p == t)
}

fn find_in_strand(record: &str, sequence: &str, pattern: &str, strand: &str) -> Vec<Hit> {
    let pat_len = pattern.len();
    let seq_len = sequence.len();
    let mut hits = Vec::new();

    if pat_len > seq_len {
        return hits;
    }

    if pattern.contains('?') {
        for i in 0..=(seq_len - pat_len) {
            let window = &sequence[i..i + pat_len];
            if wildcard_match(pattern, window) {
                hits.push(Hit {
                    record: record.to_string(),
                    start: i,
                    end: i + pat_len,
                    strand: strand.to_string(),
                    matched: window.to_string(),
                    pattern: pattern.to_string(),
                });
            }
        }
    } else {
        let matcher = Horspool::new(pattern.as_bytes());
        for i in matcher.find_all(sequence.as_bytes()) {
            hits.push(Hit {
                record: record.to_string(),
                start: i,
                end: i + pat_len,
                strand: strand.to_string(),
                matched: sequence[i..i + pat_len].to_string(),
                pattern: pattern.to_string(),
            });
        }
    }

    hits
}

pub fn find_matches(record: &str, sequence: &str, pattern: &str, include_rc: bool) -> Vec<Hit> {
    let pattern = pattern.to_uppercase();
    let mut hits = find_in_strand(record, sequence, &pattern, "+");

    if include_rc {
        let rc_seq = reverse_complement(sequence);
        let rc_hits = find_in_strand(record, &rc_seq, &pattern, "-");
        // Mirror RC indices back onto forward-strand coordinates.
        let seq_len = sequence.len();
        let pat_len = pattern.len();
        for mut hit in rc_hits {
            let fwd_start = seq_len - hit.start - pat_len;
            hit.start = fwd_start;
            hit.end = fwd_start + pat_len;
            hits.push(hit);
        }
    }

    hits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exact_match() {
        let hits = find_matches("seq1", "ATGCGTAAC", "ATGCGT", false);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].start, 0);
        assert_eq!(hits[0].end, 6);
        assert_eq!(hits[0].strand, "+");
        assert_eq!(hits[0].pattern, "ATGCGT");
    }

    #[test]
    fn test_wildcard_match() {
        let hits = find_matches("seq1", "ATGCGTAAC", "ATG???", false);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].matched, "ATGCGT");
    }

    #[test]
    fn test_reverse_complement() {
        // ATGCGT rc = ACGCAT
        let rc = reverse_complement("ATGCGT");
        assert_eq!(rc, "ACGCAT");
    }

    #[test]
    fn test_rc_search() {
        // ACGCAT is the rc of ATGCGT; search for ATGCGT on both strands in a seq that has ACGCAT
        let hits = find_matches("seq1", "TTTACGCATGGG", "ATGCGT", true);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].strand, "-");
    }

    #[test]
    fn test_no_match() {
        let hits = find_matches("seq1", "AAAAAAAA", "ATGCGT", false);
        assert!(hits.is_empty());
    }
}
