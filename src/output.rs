use anyhow::Result;
use serde::Serialize;

use crate::fasta::FastaRecord;
use crate::matcher::Hit;

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

pub fn compute_stats(record: &FastaRecord) -> RecordStats {
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
        println!(
            "{}\tpos={}\tstrand={}\tmatched={}",
            hit.record, hit.position, hit.strand, hit.matched
        );
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
    wtr.write_record(["record", "position", "strand", "matched"])?;
    for h in hits {
        wtr.write_record([
            h.record.as_str(),
            &h.position.to_string(),
            h.strand.as_str(),
            h.matched.as_str(),
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
        let end = h.position + h.matched.len();
        wtr.write_record([
            h.record.as_str(),
            &h.position.to_string(),
            &end.to_string(),
            h.matched.as_str(),
            "0",
            h.strand.as_str(),
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
