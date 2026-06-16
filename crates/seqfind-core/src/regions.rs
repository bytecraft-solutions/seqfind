//! Lightweight BED / GFF / GTF region readers used by `seqfind subseq`.
//!
//! We deliberately parse only the columns the deterministic core needs:
//! `seqid`, `start`, `end`, `strand`, and `name` (for BED) / `attributes`
//! (for GFF/GTF). Full feature graphs aren't required for extraction.

use anyhow::{Context, Result, anyhow, bail};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomicRegion {
    pub seqid: String,
    /// 1-based inclusive start.
    pub start: u64,
    /// 1-based inclusive end.
    pub end: u64,
    /// '+' / '-' / '.'.
    pub strand: char,
    /// Display label (BED name column; GFF/GTF feature type plus first attribute).
    pub label: String,
}

pub fn read_bed(path: &Path) -> Result<Vec<GenomicRegion>> {
    let f = File::open(path).with_context(|| format!("opening BED {}", path.display()))?;
    let mut out = Vec::new();
    for (i, line) in BufReader::new(f).lines().enumerate() {
        let line = line?;
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with("track") {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() < 3 {
            bail!("BED line {} has fewer than 3 columns: {trimmed}", i + 1);
        }
        let start_0: u64 = cols[1]
            .parse()
            .map_err(|e| anyhow!("BED line {}: bad start '{}': {e}", i + 1, cols[1]))?;
        let end_0: u64 = cols[2]
            .parse()
            .map_err(|e| anyhow!("BED line {}: bad end '{}': {e}", i + 1, cols[2]))?;
        if end_0 <= start_0 {
            // Skip zero-width intervals silently.
            continue;
        }
        let label = cols.get(3).copied().unwrap_or(".").to_string();
        let strand = cols
            .get(5)
            .and_then(|s| s.chars().next())
            .filter(|c| matches!(c, '+' | '-' | '.'))
            .unwrap_or('+');
        out.push(GenomicRegion {
            seqid: cols[0].to_string(),
            // BED is 0-based half-open; convert to 1-based inclusive.
            start: start_0 + 1,
            end: end_0,
            strand,
            label,
        });
    }
    Ok(out)
}

/// Read a GFF3-or-GTF-style file. Both share the same 9-column tab layout and
/// 1-based inclusive coordinates; we ignore the dialect-specific attribute
/// formatting and just keep the raw attribute field as the label.
pub fn read_gff_or_gtf(path: &Path, kind: &str) -> Result<Vec<GenomicRegion>> {
    let f = File::open(path).with_context(|| format!("opening {kind} {}", path.display()))?;
    let mut out = Vec::new();
    for (i, line) in BufReader::new(f).lines().enumerate() {
        let line = line?;
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() < 8 {
            bail!("{kind} line {} has fewer than 8 columns", i + 1);
        }
        let start: u64 = cols[3]
            .parse()
            .map_err(|e| anyhow!("{kind} line {}: bad start '{}': {e}", i + 1, cols[3]))?;
        let end: u64 = cols[4]
            .parse()
            .map_err(|e| anyhow!("{kind} line {}: bad end '{}': {e}", i + 1, cols[4]))?;
        if end < start {
            continue;
        }
        let strand = cols
            .get(6)
            .and_then(|s| s.chars().next())
            .filter(|c| matches!(c, '+' | '-' | '.'))
            .unwrap_or('.');
        let feature_type = cols.get(2).copied().unwrap_or(".");
        let attrs = cols.get(8).copied().unwrap_or("");
        let label = if attrs.is_empty() {
            feature_type.to_string()
        } else {
            format!("{feature_type}|{attrs}")
        };
        out.push(GenomicRegion {
            seqid: cols[0].to_string(),
            start,
            end,
            strand,
            label,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn parses_bed_with_strand() {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            f,
            "# comment\nchr1\t0\t10\tfeatA\t0\t+\nchr1\t20\t30\tfeatB\t.\t-"
        )
        .unwrap();
        let regs = read_bed(f.path()).unwrap();
        assert_eq!(regs.len(), 2);
        assert_eq!(regs[0].start, 1); // 0-based 0 → 1-based 1
        assert_eq!(regs[0].end, 10);
        assert_eq!(regs[1].strand, '-');
    }

    #[test]
    fn parses_gff_columns() {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            f,
            "##gff-version 3\nchr1\tsrc\tCDS\t100\t200\t.\t-\t.\tID=cds1"
        )
        .unwrap();
        let regs = read_gff_or_gtf(f.path(), "GFF").unwrap();
        assert_eq!(regs.len(), 1);
        assert_eq!(regs[0].seqid, "chr1");
        assert_eq!(regs[0].start, 100);
        assert_eq!(regs[0].end, 200);
        assert_eq!(regs[0].strand, '-');
        assert!(regs[0].label.starts_with("CDS|"));
    }
}
