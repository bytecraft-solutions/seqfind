use crate::alphabets;
use anyhow::{Context, Result};
use noodles_fasta as fasta;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct FastaRecord {
    pub name: String,
    pub sequence: String,
}

pub fn read_fasta(path: &Path) -> Result<Vec<FastaRecord>> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    read_fasta_reader(BufReader::new(file), &path.display().to_string())
}

/// Parse FASTA from an in-memory string. Used by the HTTP API which receives
/// FASTA payloads in request bodies rather than from the filesystem.
pub fn read_fasta_str(input: &str) -> Result<Vec<FastaRecord>> {
    read_fasta_reader(BufReader::new(input.as_bytes()), "<input>")
}

fn read_fasta_reader<R: std::io::BufRead>(reader: R, source: &str) -> Result<Vec<FastaRecord>> {
    let mut reader = fasta::io::Reader::new(reader);
    let mut records = Vec::new();
    for result in reader.records() {
        let record = result.with_context(|| format!("while reading {source}"))?;

        let name = String::from_utf8(record.name().to_vec())
            .with_context(|| format!("invalid FASTA record name in {source}"))?;

        let sequence: String = record
            .sequence()
            .as_ref()
            .iter()
            .map(|b| b.to_ascii_uppercase() as char)
            .collect();

        alphabets::warn_if_invalid(&name, &sequence);
        records.push(FastaRecord { name, sequence });
    }

    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn parses_multiline_fasta() {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, ">seq1 description here\nATGC\nGTAA\n>seq2\nCCC").unwrap();
        let recs = read_fasta(f.path()).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].name, "seq1");
        assert_eq!(recs[0].sequence, "ATGCGTAA");
        assert_eq!(recs[1].name, "seq2");
        assert_eq!(recs[1].sequence, "CCC");
    }
}
