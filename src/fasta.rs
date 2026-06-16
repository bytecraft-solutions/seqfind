use anyhow::{Context, Result};
use seq_io::fasta::{Reader, Record};
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum FastaError {
    #[error("invalid FASTA record: {0}")]
    Parse(String),
}

#[derive(Debug, Clone)]
pub struct FastaRecord {
    pub name: String,
    pub sequence: String,
}

pub fn read_fasta(path: &Path) -> Result<Vec<FastaRecord>> {
    let mut reader =
        Reader::from_path(path).with_context(|| format!("failed to open {}", path.display()))?;

    let mut records = Vec::new();
    while let Some(result) = reader.next() {
        let rec = result
            .map_err(|e| FastaError::Parse(e.to_string()))
            .with_context(|| format!("while reading {}", path.display()))?;

        let id = rec
            .id()
            .map_err(|e| FastaError::Parse(e.to_string()))?
            .to_string();

        let mut sequence = String::with_capacity(rec.seq().len());
        for line in rec.seq_lines() {
            for &b in line {
                sequence.push(b.to_ascii_uppercase() as char);
            }
        }

        records.push(FastaRecord { name: id, sequence });
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
