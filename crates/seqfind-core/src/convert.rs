//! Format converters -`fq2fa` (drop quality), `fx2tab` (records → TSV), `tab2fx`
//! (TSV → records).

use crate::seq::SeqRecord;
use anyhow::{Context, Result, bail};
use std::io::{BufRead, Write};

/// Drop quality strings; result is suitable for the FASTA writer.
pub fn fq2fa(records: &[SeqRecord]) -> Vec<SeqRecord> {
    records
        .iter()
        .map(|r| SeqRecord {
            name: r.name.clone(),
            sequence: r.sequence.clone(),
            quality: None,
        })
        .collect()
}

/// Serialize as `id<TAB>sequence[<TAB>quality]`. Quality column is omitted
/// when no record has it; preserves a stable column order across calls.
pub fn fx2tab<W: Write>(w: &mut W, records: &[SeqRecord]) -> Result<()> {
    let any_quality = records.iter().any(|r| r.quality.is_some());
    if any_quality {
        writeln!(w, "id\tsequence\tquality")?;
    } else {
        writeln!(w, "id\tsequence")?;
    }
    for r in records {
        if any_quality {
            let q = r.quality.as_deref().unwrap_or(&[]);
            writeln!(
                w,
                "{}\t{}\t{}",
                r.name,
                r.sequence,
                std::str::from_utf8(q).unwrap_or("")
            )?;
        } else {
            writeln!(w, "{}\t{}", r.name, r.sequence)?;
        }
    }
    Ok(())
}

/// Parse a `fx2tab`-style table back into records.
pub fn tab2fx<R: BufRead>(reader: R) -> Result<Vec<SeqRecord>> {
    let mut out = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("reading tab2fx input at line {}", i + 1))?;
        if line.is_empty() {
            continue;
        }
        // Skip a single header row matching our fx2tab format.
        if i == 0 && (line == "id\tsequence" || line == "id\tsequence\tquality") {
            continue;
        }
        let mut cols = line.split('\t');
        let id = cols.next().unwrap_or("").to_string();
        let seq = cols.next().unwrap_or("").to_string();
        let q = cols.next();
        if id.is_empty() {
            bail!("tab2fx: missing id at line {}", i + 1);
        }
        out.push(SeqRecord {
            name: id,
            sequence: seq,
            quality: q.map(|s| s.as_bytes().to_vec()),
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn rec_q(name: &str, seq: &str, q: Option<&[u8]>) -> SeqRecord {
        SeqRecord {
            name: name.into(),
            sequence: seq.into(),
            quality: q.map(|b| b.to_vec()),
        }
    }

    #[test]
    fn fq2fa_drops_quality() {
        let out = fq2fa(&[rec_q("a", "ATGC", Some(b"IIII"))]);
        assert!(out[0].quality.is_none());
    }

    #[test]
    fn fx2tab_with_quality() {
        let mut buf = Vec::new();
        fx2tab(
            &mut buf,
            &[
                rec_q("a", "ATGC", Some(b"IIII")),
                rec_q("b", "GG", Some(b"!!")),
            ],
        )
        .unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("id\tsequence\tquality\n"));
        assert!(s.contains("a\tATGC\tIIII"));
    }

    #[test]
    fn fx2tab_no_quality() {
        let mut buf = Vec::new();
        fx2tab(&mut buf, &[rec_q("a", "ATGC", None)]).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("id\tsequence\n"));
    }

    #[test]
    fn tab2fx_roundtrip() {
        let mut buf = Vec::new();
        fx2tab(
            &mut buf,
            &[
                rec_q("a", "ATGC", Some(b"IIII")),
                rec_q("b", "GG", Some(b"!!")),
            ],
        )
        .unwrap();
        let recs = tab2fx(Cursor::new(buf)).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].name, "a");
        assert_eq!(recs[0].sequence, "ATGC");
        assert_eq!(recs[0].quality.as_deref(), Some(b"IIII".as_ref()));
    }
}
