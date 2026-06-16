//! Unified FASTA + FASTQ sequence reader with transparent gzip support.
//!
//! Format and compression are auto-detected from the input itself:
//! - gzip magic bytes (`0x1f 0x8b`) → decompress through [`flate2::read::MultiGzDecoder`]
//! - first non-whitespace byte `>` → FASTA, `@` → FASTQ
//!
//! All records are returned as [`SeqRecord`]; the `quality` field is `None`
//! for FASTA inputs and populated with raw quality bytes (typically Phred+33)
//! for FASTQ inputs.

use crate::alphabets;
use anyhow::{Context, Result, bail};
use bio::io::fastq;
use flate2::read::MultiGzDecoder;
use noodles_fasta as nfasta;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Cursor, Read, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub struct SeqRecord {
    pub name: String,
    pub sequence: String,
    /// Raw quality bytes as stored in the source (no Phred offset adjustment).
    /// `None` for FASTA inputs, `Some` for FASTQ.
    pub quality: Option<Vec<u8>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Fasta,
    Fastq,
}

/// Read all sequences from `path`. Format and gzip compression are auto-detected.
/// `path == "-"` reads from stdin.
pub fn read_seqs(path: &Path) -> Result<Vec<SeqRecord>> {
    if path.as_os_str() == "-" {
        // Leak Stdin (a tiny process-singleton handle) so the lock has a
        // 'static lifetime and slots into the gzip/format dispatch path.
        let stdin: &'static std::io::Stdin = Box::leak(Box::new(std::io::stdin()));
        let locked = StdinHolder(stdin.lock());
        let reader = decode_gzip(BufReader::new(locked))?;
        return read_seqs_reader(reader, "<stdin>");
    }
    let source = path.display().to_string();
    let file = File::open(path).with_context(|| format!("failed to open {source}"))?;
    let reader = decode_gzip(BufReader::new(file))?;
    read_seqs_reader(reader, &source)
}

/// Open `path` for writing. `path == None` or `path == Some("-")` writes to stdout.
pub fn open_output(path: Option<&Path>) -> Result<Box<dyn Write>> {
    match path {
        None => Ok(Box::new(BufWriter::new(std::io::stdout()))),
        Some(p) if p.as_os_str() == "-" => Ok(Box::new(BufWriter::new(std::io::stdout()))),
        Some(p) => {
            let f = File::create(p).with_context(|| format!("failed to create {}", p.display()))?;
            Ok(Box::new(BufWriter::new(f)))
        }
    }
}

/// Path helper for `clap`: parse `"-"` as stdin/stdout sentinel.
pub fn parse_path_or_stdio(s: &str) -> PathBuf {
    PathBuf::from(s)
}

struct StdinHolder(std::io::StdinLock<'static>);

impl Read for StdinHolder {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.0.read(buf)
    }
}

/// Write a single FASTA record. `wrap` of 0 disables line wrapping.
pub fn write_fasta_record<W: Write>(w: &mut W, name: &str, seq: &str, wrap: usize) -> Result<()> {
    writeln!(w, ">{name}")?;
    if wrap == 0 {
        writeln!(w, "{seq}")?;
    } else {
        let bytes = seq.as_bytes();
        for chunk in bytes.chunks(wrap) {
            w.write_all(chunk)?;
            writeln!(w)?;
        }
    }
    Ok(())
}

/// Write a single FASTQ record. Quality defaults to all `I` (Phred 40) when missing.
pub fn write_fastq_record<W: Write>(w: &mut W, rec: &SeqRecord) -> Result<()> {
    writeln!(w, "@{}", rec.name)?;
    writeln!(w, "{}", rec.sequence)?;
    writeln!(w, "+")?;
    if let Some(q) = &rec.quality {
        w.write_all(q)?;
        writeln!(w)?;
    } else {
        // Default quality so any FASTQ writer call with a FASTA-origin record
        // still produces well-formed output (matches SeqKit's default).
        for _ in 0..rec.sequence.len() {
            w.write_all(b"I")?;
        }
        writeln!(w)?;
    }
    Ok(())
}

#[derive(Debug, Clone, Copy)]
pub enum WriteFormat {
    Fasta,
    Fastq,
}

/// Stream a slice of records to `w` in either format. Default FASTA wrap = 70.
pub fn write_records<W: Write>(
    w: &mut W,
    records: &[SeqRecord],
    fmt: WriteFormat,
    fasta_wrap: usize,
) -> Result<()> {
    match fmt {
        WriteFormat::Fasta => {
            for r in records {
                write_fasta_record(w, &r.name, &r.sequence, fasta_wrap)?;
            }
        }
        WriteFormat::Fastq => {
            for r in records {
                write_fastq_record(w, r)?;
            }
        }
    }
    Ok(())
}

/// Read all sequences from an in-memory string. Used by the HTTP API which
/// receives uncompressed FASTA/FASTQ payloads in request bodies.
pub fn read_seqs_str(input: &str) -> Result<Vec<SeqRecord>> {
    let reader: Box<dyn BufRead> = Box::new(BufReader::new(Cursor::new(input.as_bytes().to_vec())));
    read_seqs_reader(reader, "<input>")
}

fn decode_gzip<R: BufRead + 'static>(mut reader: R) -> Result<Box<dyn BufRead>> {
    let head = reader.fill_buf().context("failed to peek input header")?;
    if head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(
            BufReadAsRead(reader),
        ))))
    } else {
        Ok(Box::new(reader))
    }
}

/// Adapter: `MultiGzDecoder<R>` requires `R: Read`. `BufRead: Read` already,
/// but we wrap to keep the trait-object construction unambiguous.
struct BufReadAsRead<R: BufRead>(R);

impl<R: BufRead> Read for BufReadAsRead<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.0.read(buf)
    }
}

fn read_seqs_reader(mut reader: Box<dyn BufRead>, source: &str) -> Result<Vec<SeqRecord>> {
    let head = reader
        .fill_buf()
        .with_context(|| format!("reading {source}"))?;
    // Skip leading whitespace/BOM-style junk so an accidental leading newline
    // doesn't trip detection.
    let first = head.iter().find(|&&b| !b.is_ascii_whitespace()).copied();
    match first {
        None => Ok(Vec::new()),
        Some(b'>') => read_fasta_records(reader, source),
        Some(b'@') => read_fastq_records(reader, source),
        Some(b) => bail!(
            "unrecognized leading byte 0x{b:02x} in {source}: expected '>' (FASTA) or '@' (FASTQ)"
        ),
    }
}

fn read_fasta_records(reader: Box<dyn BufRead>, source: &str) -> Result<Vec<SeqRecord>> {
    let mut reader = nfasta::io::Reader::new(reader);
    let mut records = Vec::new();
    for result in reader.records() {
        let record = result.with_context(|| format!("while reading FASTA {source}"))?;

        let name = String::from_utf8(record.name().to_vec())
            .with_context(|| format!("invalid FASTA record name in {source}"))?;

        let sequence: String = record
            .sequence()
            .as_ref()
            .iter()
            .map(|b| b.to_ascii_uppercase() as char)
            .collect();

        alphabets::warn_if_invalid(&name, &sequence);
        records.push(SeqRecord {
            name,
            sequence,
            quality: None,
        });
    }

    Ok(records)
}

fn read_fastq_records(reader: Box<dyn BufRead>, source: &str) -> Result<Vec<SeqRecord>> {
    let reader = fastq::Reader::new(reader);
    let mut records = Vec::new();
    for result in reader.records() {
        let record = result.with_context(|| format!("while reading FASTQ {source}"))?;

        let name = record.id().to_string();
        let sequence: String = record
            .seq()
            .iter()
            .map(|b| b.to_ascii_uppercase() as char)
            .collect();
        let quality = record.qual().to_vec();

        alphabets::warn_if_invalid(&name, &sequence);
        records.push(SeqRecord {
            name,
            sequence,
            quality: Some(quality),
        });
    }

    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    #[test]
    fn parses_multiline_fasta() {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, ">seq1 description here\nATGC\nGTAA\n>seq2\nCCC").unwrap();
        let recs = read_seqs(f.path()).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].name, "seq1");
        assert_eq!(recs[0].sequence, "ATGCGTAA");
        assert!(recs[0].quality.is_none());
        assert_eq!(recs[1].name, "seq2");
        assert_eq!(recs[1].sequence, "CCC");
    }

    #[test]
    fn parses_fastq() {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "@r1 desc\nATGC\n+\nIIII\n@r2\nGG\n+\n!#").unwrap();
        let recs = read_seqs(f.path()).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].name, "r1");
        assert_eq!(recs[0].sequence, "ATGC");
        assert_eq!(recs[0].quality.as_deref(), Some(b"IIII".as_ref()));
        assert_eq!(recs[1].name, "r2");
        assert_eq!(recs[1].sequence, "GG");
        assert_eq!(recs[1].quality.as_deref(), Some(b"!#".as_ref()));
    }

    #[test]
    fn parses_gzipped_fasta() {
        let mut buf = Vec::new();
        {
            let mut enc = GzEncoder::new(&mut buf, Compression::default());
            enc.write_all(b">a\nACGT\n>b\nTTTT\n").unwrap();
            enc.finish().unwrap();
        }
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(&buf).unwrap();
        let recs = read_seqs(f.path()).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].sequence, "ACGT");
        assert_eq!(recs[1].sequence, "TTTT");
    }

    #[test]
    fn parses_gzipped_fastq() {
        let mut buf = Vec::new();
        {
            let mut enc = GzEncoder::new(&mut buf, Compression::default());
            enc.write_all(b"@x\nAAAA\n+\n####\n").unwrap();
            enc.finish().unwrap();
        }
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(&buf).unwrap();
        let recs = read_seqs(f.path()).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].name, "x");
        assert_eq!(recs[0].sequence, "AAAA");
        assert_eq!(recs[0].quality.as_deref(), Some(b"####".as_ref()));
    }

    #[test]
    fn empty_input_is_empty() {
        let f = tempfile::NamedTempFile::new().unwrap();
        let recs = read_seqs(f.path()).unwrap();
        assert!(recs.is_empty());
    }

    #[test]
    fn unrecognized_format_errors() {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "not a real format").unwrap();
        let err = read_seqs(f.path()).unwrap_err();
        assert!(err.to_string().contains("unrecognized leading byte"));
    }

    #[test]
    fn parses_fasta_from_str() {
        let recs = read_seqs_str(">a\nACGT\n>b\nTT\n").unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].name, "a");
        assert_eq!(recs[1].sequence, "TT");
    }
}
