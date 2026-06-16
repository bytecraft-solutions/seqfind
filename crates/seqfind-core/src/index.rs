//! FM-index over reference sequences, used by `compare --fm-index` to amortize
//! index construction across many query files / many query records against a
//! chromosome-scale reference. For small refs (E. coli, SARS-CoV-2) Horspool
//! or BNDM are faster end-to-end; the FM-index pays off when the same
//! reference is hit by many queries.

use crate::alphabets::{fm_index_alphabet, normalize_for_fm_index};
use crate::matcher::reverse_complement;
use bio::data_structures::bwt::{BWT, Less, Occ, bwt, less};
use bio::data_structures::fmindex::{BackwardSearchResult, FMIndex, FMIndexable};
use bio::data_structures::suffix_array::{RawSuffixArray, suffix_array};

pub struct ReferenceIndex {
    pub name: String,
    pub length: usize,
    sa: RawSuffixArray,
    bwt: BWT,
    less: Less,
    occ: Occ,
}

impl ReferenceIndex {
    pub fn build(name: String, sequence: &str) -> Self {
        let mut text = normalize_for_fm_index(sequence);
        text.push(b'$');
        let alphabet = fm_index_alphabet();
        let sa = suffix_array(&text);
        let bwt_data = bwt(&text, &sa);
        let less_data = less(&bwt_data, &alphabet);
        let occ_data = Occ::new(&bwt_data, 3, &alphabet);
        let length = sequence.len();
        Self {
            name,
            length,
            sa,
            bwt: bwt_data,
            less: less_data,
            occ: occ_data,
        }
    }

    fn fm(&self) -> FMIndex<&BWT, &Less, &Occ> {
        FMIndex::new(&self.bwt, &self.less, &self.occ)
    }

    fn search_normalized(&self, pattern: &[u8]) -> Vec<usize> {
        if pattern.is_empty() || pattern.len() > self.length {
            return Vec::new();
        }
        let result = self.fm().backward_search(pattern.iter());
        let interval = match result {
            BackwardSearchResult::Complete(i) => i,
            _ => return Vec::new(),
        };
        let mut positions = interval.occ(&self.sa);
        // Drop matches that bleed into the `$` sentinel (defensive -should
        // not happen since the sentinel is not in the query alphabet).
        positions.retain(|p| p + pattern.len() <= self.length);
        positions.sort_unstable();
        positions
    }

    /// Forward-strand hit positions of `pattern` in this reference record.
    pub fn find_forward(&self, pattern: &str) -> Vec<usize> {
        let normalized = normalize_for_fm_index(pattern);
        self.search_normalized(&normalized)
    }

    /// Forward-strand positions where the reverse complement of `pattern`
    /// occurs in this reference -i.e. BLAST-style minus-strand hits.
    pub fn find_reverse(&self, pattern: &str) -> Vec<usize> {
        let rc = reverse_complement(pattern);
        let normalized = normalize_for_fm_index(&rc);
        self.search_normalized(&normalized)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_forward_and_rc_hits() {
        let index = ReferenceIndex::build("ref".into(), "TTACGCATGGGCCCAAATTTATGCGTAAA");
        let fwd = index.find_forward("ATGCGT");
        assert_eq!(fwd, vec![20]);
        // ACGCAT is at index 2 -minus-strand hit for ATGCGT.
        let rc = index.find_reverse("ATGCGT");
        assert_eq!(rc, vec![2]);
    }

    #[test]
    fn no_hits_for_absent_pattern() {
        let index = ReferenceIndex::build("ref".into(), "AAAAAAAA");
        assert!(index.find_forward("CCCC").is_empty());
    }
}
