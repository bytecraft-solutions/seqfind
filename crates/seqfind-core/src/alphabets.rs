//! DNA alphabet validation / normalization helpers.
//!
//! Uses Rust-Bio's IUPAC DNA alphabet to detect bases that fall outside the
//! standard nucleotide set. Today this is a soft warning emitted at read time
//! -it does not fail the run, since real-world FASTA files routinely include
//! ambiguity codes and soft-masked lowercase bases.

use bio::alphabets::Alphabet;
use bio::alphabets::dna;
use std::sync::OnceLock;

fn iupac() -> &'static Alphabet {
    static ALPHA: OnceLock<Alphabet> = OnceLock::new();
    ALPHA.get_or_init(dna::iupac_alphabet)
}

fn fmindex_alphabet_bytes() -> &'static [u8] {
    // Sequences are uppercased on read; add `$` for the FM-index sentinel and
    // `N` for the catch-all bucket every non-ACGT base gets normalized to.
    b"$ACGTN"
}

/// Normalize a sequence into the strict {A,C,G,T,N,$} alphabet that the
/// FM-index expects. Any IUPAC ambiguity code or stray character collapses to
/// `N`; U (RNA) is translated to T. The result owns its bytes.
pub fn normalize_for_fm_index(sequence: &str) -> Vec<u8> {
    sequence
        .bytes()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' | b'C' | b'G' | b'T' | b'N' => b.to_ascii_uppercase(),
            b'U' => b'T',
            _ => b'N',
        })
        .collect()
}

/// Alphabet used to build the FM-index over a normalized sequence + sentinel.
pub fn fm_index_alphabet() -> Alphabet {
    Alphabet::new(fmindex_alphabet_bytes())
}

/// Returns the count of bases that lie outside the IUPAC DNA alphabet.
pub fn count_invalid_bases(sequence: &str) -> usize {
    let alpha = iupac();
    sequence
        .bytes()
        .filter(|b| !alpha.is_word(std::slice::from_ref(b)))
        .count()
}

/// Emit a single stderr warning if `sequence` contains characters outside the
/// IUPAC DNA alphabet. Quiet when the record is clean.
pub fn warn_if_invalid(record_name: &str, sequence: &str) {
    let invalid = count_invalid_bases(sequence);
    if invalid > 0 {
        eprintln!(
            "warning: record '{record_name}' has {invalid} base(s) outside the IUPAC DNA alphabet"
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_invalid_bases() {
        assert_eq!(count_invalid_bases("ACGTN"), 0);
        // IUPAC codes like R/Y/W/S are valid DNA ambiguity codes.
        assert_eq!(count_invalid_bases("ACGTRYWS"), 0);
        assert_eq!(count_invalid_bases("ACGT*?"), 2);
    }

    #[test]
    fn normalizes_to_fm_alphabet() {
        let bytes = normalize_for_fm_index("AcGtNuRYWX");
        assert_eq!(bytes, b"ACGTNTNNNN");
    }
}
