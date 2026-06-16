use anyhow::Result;
use serde::Serialize;

use crate::matcher::Hit;
use crate::orf::Orf;
use crate::seq::SeqRecord;

#[derive(Debug, Clone, Serialize)]
pub struct RecordStats {
    pub record: String,
    pub length: usize,
    pub gc_percent: f64,
    pub a: usize,
    pub t: usize,
    pub g: usize,
    pub c: usize,
    pub n: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct CompareHit {
    pub reference_file: String,
    pub reference_record: String,
    pub query_file: String,
    pub query_record: String,
    pub start: usize,
    pub end: usize,
    pub strand: String,
    pub matched: String,
    pub query_length: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub edits: Option<usize>,
}

#[derive(Debug, Clone, Serialize)]
pub struct Variant {
    pub reference_file: String,
    pub reference_record: String,
    pub query_file: String,
    pub query_record: String,
    /// 1-based reference coordinate.
    pub position: usize,
    #[serde(rename = "ref")]
    pub reference: String,
    pub alt: String,
    #[serde(rename = "type")]
    pub variant_type: String,
}

pub fn compute_stats(record: &SeqRecord) -> RecordStats {
    let mut a = 0usize;
    let mut t = 0usize;
    let mut g = 0usize;
    let mut c = 0usize;
    let mut n = 0usize;

    for ch in record.sequence.chars() {
        match ch {
            'A' => a += 1,
            'T' => t += 1,
            'G' => g += 1,
            'C' => c += 1,
            _ => n += 1,
        }
    }

    let len = record.sequence.len();
    let gc_percent = if len == 0 {
        0.0
    } else {
        (g + c) as f64 / len as f64 * 100.0
    };

    RecordStats {
        record: record.name.clone(),
        length: len,
        gc_percent,
        a,
        t,
        g,
        c,
        n,
    }
}

// --- hits output --------------------------------------------------------

pub fn print_hits_text(hits: &[Hit]) {
    if hits.is_empty() {
        println!("No matches found.");
        return;
    }
    for hit in hits {
        let mut line = format!(
            "{}\tstart={}\tend={}\tstrand={}\tmatched={}\tpattern={}",
            hit.record, hit.start, hit.end, hit.strand, hit.matched, hit.pattern
        );
        if let Some(edits) = hit.edits {
            line.push_str(&format!("\tedits={edits}"));
        }
        println!("{line}");
    }
}

pub fn print_hits_json(hits: &[Hit]) -> Result<()> {
    println!("{}", serde_json::to_string_pretty(hits)?);
    Ok(())
}

fn write_hits_delimited(hits: &[Hit], delimiter: u8) -> Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(delimiter)
        .from_writer(std::io::stdout());
    wtr.write_record([
        "record", "start", "end", "strand", "matched", "pattern", "edits",
    ])?;
    for h in hits {
        let edits = h.edits.map(|e| e.to_string()).unwrap_or_default();
        wtr.write_record([
            h.record.as_str(),
            &h.start.to_string(),
            &h.end.to_string(),
            h.strand.as_str(),
            h.matched.as_str(),
            h.pattern.as_str(),
            &edits,
        ])?;
    }
    wtr.flush()?;
    Ok(())
}

pub fn print_hits_csv(hits: &[Hit]) -> Result<()> {
    write_hits_delimited(hits, b',')
}

pub fn print_hits_tsv(hits: &[Hit]) -> Result<()> {
    write_hits_delimited(hits, b'\t')
}

/// BED3+3: chrom, start (0-based), end (exclusive), name, score, strand.
pub fn print_hits_bed(hits: &[Hit]) -> Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(std::io::stdout());
    for h in hits {
        let score = h.edits.map(|e| e.to_string()).unwrap_or_else(|| "0".into());
        wtr.write_record([
            h.record.as_str(),
            &h.start.to_string(),
            &h.end.to_string(),
            h.matched.as_str(),
            &score,
            h.strand.as_str(),
        ])?;
    }
    wtr.flush()?;
    Ok(())
}

// --- compare output -----------------------------------------------------

pub fn print_compare_text(hits: &[CompareHit]) {
    if hits.is_empty() {
        println!("No compare matches found.");
        return;
    }
    for hit in hits {
        let mut line = format!(
            "{}:{}\tstart={}\tend={}\tstrand={}\tquery_file={}\tquery_record={}\tquery_length={}\tmatched={}",
            hit.reference_file,
            hit.reference_record,
            hit.start,
            hit.end,
            hit.strand,
            hit.query_file,
            hit.query_record,
            hit.query_length,
            hit.matched
        );
        if let Some(edits) = hit.edits {
            line.push_str(&format!("\tedits={edits}"));
        }
        println!("{line}");
    }
}

pub fn print_compare_json(hits: &[CompareHit]) -> Result<()> {
    println!("{}", serde_json::to_string_pretty(hits)?);
    Ok(())
}

fn write_compare_delimited(hits: &[CompareHit], delimiter: u8) -> Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(delimiter)
        .from_writer(std::io::stdout());
    wtr.write_record([
        "reference_file",
        "reference_record",
        "query_file",
        "query_record",
        "start",
        "end",
        "strand",
        "matched",
        "query_length",
        "edits",
    ])?;
    for h in hits {
        let edits = h.edits.map(|e| e.to_string()).unwrap_or_default();
        wtr.write_record([
            h.reference_file.as_str(),
            h.reference_record.as_str(),
            h.query_file.as_str(),
            h.query_record.as_str(),
            &h.start.to_string(),
            &h.end.to_string(),
            h.strand.as_str(),
            h.matched.as_str(),
            &h.query_length.to_string(),
            &edits,
        ])?;
    }
    wtr.flush()?;
    Ok(())
}

pub fn print_compare_csv(hits: &[CompareHit]) -> Result<()> {
    write_compare_delimited(hits, b',')
}

pub fn print_compare_tsv(hits: &[CompareHit]) -> Result<()> {
    write_compare_delimited(hits, b'\t')
}

/// BED3+3 over reference coordinates. Name is `query_file:query_record`;
/// score is the query length.
pub fn print_compare_bed(hits: &[CompareHit]) -> Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(std::io::stdout());
    for h in hits {
        let name = format!("{}:{}", h.query_file, h.query_record);
        let score = h
            .edits
            .map(|e| e.to_string())
            .unwrap_or_else(|| h.query_length.to_string());
        wtr.write_record([
            h.reference_record.as_str(),
            &h.start.to_string(),
            &h.end.to_string(),
            &name,
            &score,
            h.strand.as_str(),
        ])?;
    }
    wtr.flush()?;
    Ok(())
}

// --- variant output -----------------------------------------------------

pub fn print_variants_text(variants: &[Variant]) {
    if variants.is_empty() {
        println!("No variants found.");
        return;
    }
    for v in variants {
        println!(
            "{}:{}\tposition={}\tref={}\talt={}\ttype={}\tquery_file={}\tquery_record={}",
            v.reference_file,
            v.reference_record,
            v.position,
            v.reference,
            v.alt,
            v.variant_type,
            v.query_file,
            v.query_record
        );
    }
}

pub fn print_variants_json(variants: &[Variant]) -> Result<()> {
    println!("{}", serde_json::to_string_pretty(variants)?);
    Ok(())
}

fn write_variants_delimited(variants: &[Variant], delimiter: u8) -> Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(delimiter)
        .from_writer(std::io::stdout());
    wtr.write_record([
        "reference_file",
        "reference_record",
        "query_file",
        "query_record",
        "position",
        "ref",
        "alt",
        "type",
    ])?;
    for v in variants {
        wtr.write_record([
            v.reference_file.as_str(),
            v.reference_record.as_str(),
            v.query_file.as_str(),
            v.query_record.as_str(),
            &v.position.to_string(),
            v.reference.as_str(),
            v.alt.as_str(),
            v.variant_type.as_str(),
        ])?;
    }
    wtr.flush()?;
    Ok(())
}

pub fn print_variants_csv(variants: &[Variant]) -> Result<()> {
    write_variants_delimited(variants, b',')
}

pub fn print_variants_tsv(variants: &[Variant]) -> Result<()> {
    write_variants_delimited(variants, b'\t')
}

pub fn print_variants_bed(variants: &[Variant]) -> Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(std::io::stdout());
    for v in variants {
        let start = v.position.saturating_sub(1);
        let end = if v.reference.is_empty() {
            start
        } else {
            start + v.reference.len()
        };
        let name = format!(
            "{}:{}:{}>{}",
            v.query_file, v.query_record, v.reference, v.alt
        );
        wtr.write_record([
            v.reference_record.as_str(),
            &start.to_string(),
            &end.to_string(),
            &name,
            "0",
            "+",
        ])?;
    }
    wtr.flush()?;
    Ok(())
}

// --- stats output -------------------------------------------------------

pub fn print_stats_text(stats: &[RecordStats]) {
    for s in stats {
        println!(
            "{}\tlength={}\tGC={:.1}%\tA={} T={} G={} C={} N={}",
            s.record, s.length, s.gc_percent, s.a, s.t, s.g, s.c, s.n
        );
    }
}

pub fn print_stats_json(stats: &[RecordStats]) -> Result<()> {
    println!("{}", serde_json::to_string_pretty(stats)?);
    Ok(())
}

fn write_stats_delimited(stats: &[RecordStats], delimiter: u8) -> Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(delimiter)
        .from_writer(std::io::stdout());
    wtr.write_record(["record", "length", "gc_percent", "a", "t", "g", "c", "n"])?;
    for s in stats {
        wtr.write_record([
            s.record.as_str(),
            &s.length.to_string(),
            &format!("{:.4}", s.gc_percent),
            &s.a.to_string(),
            &s.t.to_string(),
            &s.g.to_string(),
            &s.c.to_string(),
            &s.n.to_string(),
        ])?;
    }
    wtr.flush()?;
    Ok(())
}

pub fn print_stats_csv(stats: &[RecordStats]) -> Result<()> {
    write_stats_delimited(stats, b',')
}

pub fn print_stats_tsv(stats: &[RecordStats]) -> Result<()> {
    write_stats_delimited(stats, b'\t')
}

// --- ORF output --------------------------------------------------------

pub fn print_orfs_text(orfs: &[Orf]) {
    if orfs.is_empty() {
        println!("No ORFs found.");
        return;
    }
    for o in orfs {
        println!(
            "{}\tstart={}\tend={}\tlength={}\tstrand={}\tstart_codon={}\tstop_codon={}",
            o.record, o.start, o.end, o.length, o.strand, o.start_codon, o.stop_codon,
        );
    }
}

pub fn print_orfs_json(orfs: &[Orf]) -> Result<()> {
    println!("{}", serde_json::to_string_pretty(orfs)?);
    Ok(())
}

fn write_orfs_delimited(orfs: &[Orf], delimiter: u8) -> Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(delimiter)
        .from_writer(std::io::stdout());
    wtr.write_record([
        "record",
        "start",
        "end",
        "length",
        "strand",
        "start_codon",
        "stop_codon",
    ])?;
    for o in orfs {
        wtr.write_record([
            o.record.as_str(),
            &o.start.to_string(),
            &o.end.to_string(),
            &o.length.to_string(),
            o.strand.as_str(),
            o.start_codon.as_str(),
            o.stop_codon.as_str(),
        ])?;
    }
    wtr.flush()?;
    Ok(())
}

pub fn print_orfs_csv(orfs: &[Orf]) -> Result<()> {
    write_orfs_delimited(orfs, b',')
}

pub fn print_orfs_tsv(orfs: &[Orf]) -> Result<()> {
    write_orfs_delimited(orfs, b'\t')
}

pub fn print_orfs_bed(orfs: &[Orf]) -> Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(std::io::stdout());
    for o in orfs {
        let name = format!("ORF_{}_{}", o.start_codon, o.stop_codon);
        wtr.write_record([
            o.record.as_str(),
            &o.start.to_string(),
            &o.end.to_string(),
            &name,
            &o.length.to_string(),
            o.strand.as_str(),
        ])?;
    }
    wtr.flush()?;
    Ok(())
}

// --- FASTA output (rc --file) ------------------------------------------

/// Write a FASTA record to stdout with sequence lines wrapped at `width` bp.
pub fn print_fasta_record(name: &str, sequence: &str, width: usize) {
    println!(">{name}");
    if width == 0 {
        println!("{sequence}");
        return;
    }
    let bytes = sequence.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let end = (i + width).min(bytes.len());
        // SAFETY: input is ASCII (we uppercase on read).
        println!("{}", std::str::from_utf8(&bytes[i..end]).unwrap());
        i = end;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::SeqRecord;

    fn rec(name: &str, sequence: &str) -> SeqRecord {
        SeqRecord {
            name: name.to_string(),
            sequence: sequence.to_string(),
            quality: None,
        }
    }

    #[test]
    fn gc_percent_basic() {
        // 4 G/C out of 8 = 50%
        let s = compute_stats(&rec("seq1", "ATGCATGC"));
        assert_eq!(s.length, 8);
        assert_eq!(s.g, 2);
        assert_eq!(s.c, 2);
        assert!((s.gc_percent - 50.0).abs() < 1e-9);
    }

    #[test]
    fn gc_percent_all_at() {
        let s = compute_stats(&rec("seq2", "AAAATTTT"));
        assert_eq!(s.gc_percent, 0.0);
    }

    #[test]
    fn gc_percent_all_gc() {
        let s = compute_stats(&rec("seq3", "GGGCCC"));
        assert!((s.gc_percent - 100.0).abs() < 1e-9);
    }

    #[test]
    fn gc_percent_empty_sequence() {
        let s = compute_stats(&rec("seq4", ""));
        assert_eq!(s.length, 0);
        assert_eq!(s.gc_percent, 0.0);
    }

    #[test]
    fn gc_percent_counts_n_separately() {
        // 1 G out of 4 known bases, but length is 4, so GC% = 25%.
        let s = compute_stats(&rec("seq5", "ANGN"));
        assert_eq!(s.n, 2);
        assert_eq!(s.length, 4);
        assert!((s.gc_percent - 25.0).abs() < 1e-9);
    }
}
