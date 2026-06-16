//! Thin NCBI Efetch client.
//!
//! Streams FASTA responses directly to disk so large records, such as human
//! chromosomes, do not need to fit in memory.

use anyhow::{Context, Result, bail};
use std::fs::{self, File};
use std::io::{self, BufWriter, Read, Write};
use std::path::Path;
use std::time::Duration;

const EFETCH_URL: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";

/// Download a FASTA record for `accession` from NCBI Nucleotide and write it
/// to `out_path`. Returns the number of bytes written.
pub fn fetch_fasta(accession: &str, out_path: &Path, timeout: Duration) -> Result<u64> {
    if accession.trim().is_empty() {
        bail!("accession is empty");
    }
    if timeout.is_zero() {
        bail!("fetch timeout must be greater than zero seconds");
    }

    let agent = ureq::AgentBuilder::new()
        .timeout_connect(Duration::from_secs(30))
        .timeout_read(timeout)
        .timeout_write(Duration::from_secs(30))
        .user_agent("seqfind/0.1 (+https://github.com/binyuga/seqfind)")
        .build();

    let resp = agent
        .get(EFETCH_URL)
        .query("db", "nuccore")
        .query("id", accession)
        .query("rettype", "fasta")
        .query("retmode", "text")
        .call()
        .with_context(|| format!("NCBI Efetch request failed for {accession}"))?;

    let status = resp.status();
    if status != 200 {
        bail!("NCBI returned HTTP {status} for {accession}");
    }

    let mut reader = resp.into_reader();

    let mut first = [0_u8; 1];
    let n = reader
        .read(&mut first)
        .with_context(|| format!("failed to read first byte from NCBI response for {accession}"))?;

    if n == 0 {
        bail!("NCBI returned an empty response for {accession}");
    }

    if first[0] != b'>' {
        // Efetch returns a plain-text error message with no FASTA header on bad IDs.
        let mut preview_bytes = vec![first[0]];
        let mut limited = reader.take(159);
        limited
            .read_to_end(&mut preview_bytes)
            .with_context(|| format!("failed to read NCBI error preview for {accession}"))?;
        let preview = String::from_utf8_lossy(&preview_bytes);
        bail!("NCBI did not return a FASTA record for {accession}: {preview}");
    }

    if let Some(parent) = out_path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    let file = File::create(out_path)
        .with_context(|| format!("failed to create {}", out_path.display()))?;
    let mut writer = BufWriter::new(file);

    writer
        .write_all(&first)
        .with_context(|| format!("failed to write {}", out_path.display()))?;
    let copied = io::copy(&mut reader, &mut writer)
        .with_context(|| format!("failed to stream NCBI response body for {accession}"))?;
    writer
        .flush()
        .with_context(|| format!("failed to flush {}", out_path.display()))?;

    Ok(copied + 1)
}
