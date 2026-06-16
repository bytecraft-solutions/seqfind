//! Per-record digests -SeqKit's `seqkit sum`. SHA-256 hex over the uppercased
//! sequence (canonicalised so the digest survives case/wrap differences).

use crate::seq::SeqRecord;
use serde::Serialize;
use sha2::{Digest, Sha256};

#[derive(Debug, Clone, Serialize)]
pub struct RecordDigest {
    pub record: String,
    pub length: usize,
    pub sha256: String,
}

pub fn sha256_hex(seq: &str) -> String {
    let mut h = Sha256::new();
    for &b in seq.as_bytes() {
        h.update([b.to_ascii_uppercase()]);
    }
    let bytes = h.finalize();
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

pub fn digest_records(records: &[SeqRecord]) -> Vec<RecordDigest> {
    records
        .iter()
        .map(|r| RecordDigest {
            record: r.name.clone(),
            length: r.sequence.len(),
            sha256: sha256_hex(&r.sequence),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sha256_is_case_insensitive() {
        assert_eq!(sha256_hex("ATGC"), sha256_hex("atgc"));
    }

    #[test]
    fn sha256_distinguishes_different_sequences() {
        assert_ne!(sha256_hex("ATGC"), sha256_hex("ATGCG"));
        assert_ne!(sha256_hex("ATGC"), sha256_hex("ATGT"));
    }

    #[test]
    fn sha256_hex_length_is_64() {
        assert_eq!(sha256_hex("ATGC").len(), 64);
        assert!(sha256_hex("").chars().all(|c| c.is_ascii_hexdigit()));
    }
}
