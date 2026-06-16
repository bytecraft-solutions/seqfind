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
    let mut reader = fasta::io::Reader::new(BufReader::new(file));

    let mut records = Vec::new();
    for result in reader.records() {
        let record = result.with_context(|| format!("while reading {}", path.display()))?;

        let name = String::from_utf8(record.name().to_vec())
            .with_context(|| format!("invalid FASTA record name in {}", path.display()))?;

        let sequence = record
            .sequence()
            .as_ref()
            .iter()
            .map(|b| b.to_ascii_uppercase() as char)
            .collect();

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
