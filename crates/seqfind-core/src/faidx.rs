//! Plain-text FASTA random access via a `.fai` index -SeqKit's `faidx`.
//!
//! The on-disk format matches samtools / `noodles-fasta`:
//! `name\tlength\toffset\tlinebases\tlinewidth`
//!
//! We build the index by scanning the file once. Subsequent region lookups
//! `chr:start-end` (1-based inclusive) seek directly into the file without
//! reparsing the whole record. Gzipped inputs are not supported (FAI assumes
//! random-access byte offsets in the uncompressed file).

use anyhow::{Context, Result, anyhow, bail};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub struct FaiEntry {
    pub name: String,
    pub length: u64,
    /// Byte offset in the FASTA file of the first base of this record.
    pub offset: u64,
    /// Number of bases on each (non-final) sequence line.
    pub line_bases: u64,
    /// Number of bytes on each (non-final) sequence line (includes line-terminator).
    pub line_width: u64,
}

#[derive(Debug, Default)]
pub struct FaiIndex {
    pub entries: Vec<FaiEntry>,
}

impl FaiIndex {
    pub fn find(&self, name: &str) -> Option<&FaiEntry> {
        self.entries.iter().find(|e| e.name == name)
    }

    pub fn write_to<W: Write>(&self, w: &mut W) -> Result<()> {
        for e in &self.entries {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}",
                e.name, e.length, e.offset, e.line_bases, e.line_width
            )?;
        }
        Ok(())
    }

    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let f =
            File::open(&path).with_context(|| format!("opening {}", path.as_ref().display()))?;
        let mut entries = Vec::new();
        for line in BufReader::new(f).lines() {
            let line = line?;
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 5 {
                bail!("malformed .fai line: {line}");
            }
            entries.push(FaiEntry {
                name: parts[0].to_string(),
                length: parts[1].parse()?,
                offset: parts[2].parse()?,
                line_bases: parts[3].parse()?,
                line_width: parts[4].parse()?,
            });
        }
        Ok(Self { entries })
    }
}

/// Build a `.fai` index for `fasta_path` by scanning the file once.
pub fn build_index(fasta_path: &Path) -> Result<FaiIndex> {
    let file =
        File::open(fasta_path).with_context(|| format!("opening {}", fasta_path.display()))?;
    let mut reader = BufReader::new(file);
    let mut entries: Vec<FaiEntry> = Vec::new();
    let mut offset: u64 = 0;
    let mut line = Vec::new();
    let mut current: Option<FaiEntry> = None;
    let mut line_bases: u64 = 0;
    let mut line_width: u64 = 0;
    let mut seen_short_line = false;

    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        // Detect line terminator length: \n (1), \r\n (2), or no terminator at EOF.
        let has_newline = line.last() == Some(&b'\n');
        let mut content_len = line.len();
        if has_newline {
            content_len -= 1;
            if content_len > 0 && line[content_len - 1] == b'\r' {
                content_len -= 1;
            }
        }
        let content = &line[..content_len];
        let term_len = line.len() - content_len;

        if !content.is_empty() && content[0] == b'>' {
            if let Some(prev) = current.take() {
                entries.push(prev);
            }
            // Header: parse name = first whitespace-delimited token after '>'.
            let after = &content[1..];
            let name_end = after
                .iter()
                .position(|&b| b == b' ' || b == b'\t')
                .unwrap_or(after.len());
            let name = std::str::from_utf8(&after[..name_end])
                .map_err(|e| anyhow!("invalid UTF-8 in FASTA header: {e}"))?
                .to_string();
            offset += line.len() as u64;
            current = Some(FaiEntry {
                name,
                length: 0,
                offset,
                line_bases: 0,
                line_width: 0,
            });
            line_bases = 0;
            line_width = 0;
            seen_short_line = false;
            continue;
        }

        if let Some(entry) = current.as_mut() {
            let bases = content.len() as u64;
            let width = bases + term_len as u64;
            if entry.line_bases == 0 {
                entry.line_bases = bases;
                entry.line_width = width;
                line_bases = bases;
                line_width = width;
            } else if !seen_short_line {
                // Subsequent lines must match previous unless this is the final line.
                if bases != line_bases || width != line_width {
                    seen_short_line = true;
                }
            } else {
                bail!(
                    "FASTA record '{}' has inconsistent line widths; cannot build .fai",
                    entry.name
                );
            }
            entry.length += bases;
        }
        offset += line.len() as u64;
    }
    if let Some(prev) = current.take() {
        entries.push(prev);
    }
    Ok(FaiIndex { entries })
}

/// Default `.fai` location: append `.fai` to the FASTA path.
pub fn default_fai_path(fasta: &Path) -> PathBuf {
    let mut p = fasta.as_os_str().to_os_string();
    p.push(".fai");
    PathBuf::from(p)
}

/// Resolved 1-based-inclusive coordinates for a region lookup.
#[derive(Debug, Clone)]
pub struct Region {
    pub name: String,
    pub start: u64,
    pub end: u64,
}

impl Region {
    /// Parse `chr:start-end` (1-based inclusive). Bare `chr` selects the
    /// whole record (1..=length, resolved later).
    pub fn parse(s: &str) -> Result<Self> {
        if let Some((name, range)) = s.split_once(':') {
            let (a, b) = range
                .split_once('-')
                .ok_or_else(|| anyhow!("region must be NAME:START-END, got '{s}'"))?;
            Ok(Self {
                name: name.to_string(),
                start: a
                    .parse()
                    .map_err(|e| anyhow!("invalid region start '{a}': {e}"))?,
                end: b
                    .parse()
                    .map_err(|e| anyhow!("invalid region end '{b}': {e}"))?,
            })
        } else {
            Ok(Self {
                name: s.to_string(),
                start: 1,
                end: 0, // sentinel -caller fills in with entry.length
            })
        }
    }
}

/// Fetch the sequence bytes for `region` from the indexed FASTA.
pub fn fetch_region(fasta: &Path, index: &FaiIndex, region: &Region) -> Result<String> {
    let entry = index
        .find(&region.name)
        .ok_or_else(|| anyhow!("record '{}' not in index", region.name))?;
    let end = if region.end == 0 {
        entry.length
    } else {
        region.end
    };
    if region.start == 0 || end > entry.length || end < region.start {
        bail!(
            "region {}:{}-{} out of range for record length {}",
            region.name,
            region.start,
            end,
            entry.length
        );
    }
    let count = end - region.start + 1;
    if entry.line_bases == 0 {
        return Ok(String::new());
    }
    let start_idx = region.start - 1;
    let line_idx = start_idx / entry.line_bases;
    let col = start_idx % entry.line_bases;
    let byte_offset = entry.offset + line_idx * entry.line_width + col;

    let mut f = File::open(fasta).with_context(|| format!("opening {}", fasta.display()))?;
    f.seek(SeekFrom::Start(byte_offset))?;

    let mut out = String::with_capacity(count as usize);
    let mut remaining = count;
    // Walk through (line_bases - col) chunk first, then full lines.
    let mut chunk_cap = entry.line_bases - col;
    while remaining > 0 {
        let take = chunk_cap.min(remaining);
        let mut buf = vec![0u8; take as usize];
        f.read_exact(&mut buf)
            .with_context(|| format!("reading record {} bytes", region.name))?;
        out.push_str(
            std::str::from_utf8(&buf)
                .map_err(|e| anyhow!("non-UTF-8 sequence byte in {}: {e}", region.name))?,
        );
        remaining -= take;
        if remaining == 0 {
            break;
        }
        // Skip line terminator (line_width - line_bases bytes).
        let skip = entry.line_width - entry.line_bases;
        let mut skip_buf = vec![0u8; skip as usize];
        f.read_exact(&mut skip_buf)?;
        chunk_cap = entry.line_bases;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn temp_fasta(body: &[u8]) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(body).unwrap();
        f.flush().unwrap();
        f
    }

    #[test]
    fn build_and_query_simple() {
        // 2 records, line_bases=4
        let f = temp_fasta(b">chr1\nACGT\nACGT\nAC\n>chr2\nTTTT\nTT\n");
        let idx = build_index(f.path()).unwrap();
        assert_eq!(idx.entries.len(), 2);
        assert_eq!(idx.entries[0].name, "chr1");
        assert_eq!(idx.entries[0].length, 10);
        assert_eq!(idx.entries[0].line_bases, 4);
        assert_eq!(idx.entries[1].length, 6);

        let r = Region::parse("chr1:1-10").unwrap();
        assert_eq!(fetch_region(f.path(), &idx, &r).unwrap(), "ACGTACGTAC");

        let r = Region::parse("chr1:5-8").unwrap();
        assert_eq!(fetch_region(f.path(), &idx, &r).unwrap(), "ACGT");

        let r = Region::parse("chr2:3-5").unwrap();
        assert_eq!(fetch_region(f.path(), &idx, &r).unwrap(), "TTT");
    }

    #[test]
    fn full_record_when_no_range() {
        let f = temp_fasta(b">x\nACGT\nACGT\n");
        let idx = build_index(f.path()).unwrap();
        let r = Region::parse("x").unwrap();
        assert_eq!(fetch_region(f.path(), &idx, &r).unwrap(), "ACGTACGT");
    }
}
