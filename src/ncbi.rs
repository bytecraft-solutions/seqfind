//! Thin NCBI Efetch client. One blocking GET, no streaming or caching yet -
//! reference genomes are small enough that simple wins.

use anyhow::{Context, Result, bail};
use std::fs;
use std::path::Path;
use std::time::Duration;

const EFETCH_URL: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";

/// Download a FASTA record for `accession` from NCBI Nucleotide and write it
/// to `out_path`. Returns the number of bytes written.
pub fn fetch_fasta(accession: &str, out_path: &Path) -> Result<u64> {
    if accession.trim().is_empty() {
        bail!("accession is empty");
    }

    let agent = ureq::AgentBuilder::new()
        .timeout(Duration::from_secs(60))
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

    let body = resp
        .into_string()
        .with_context(|| format!("failed to read NCBI response body for {accession}"))?;

    if !body.starts_with('>') {
        // Efetch returns a plain-text error message with no FASTA header on bad IDs.
        let preview: String = body.chars().take(160).collect();
        bail!("NCBI did not return a FASTA record for {accession}: {preview}");
    }

    if let Some(parent) = out_path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    fs::write(out_path, &body)
        .with_context(|| format!("failed to write {}", out_path.display()))?;

    Ok(body.len() as u64)
}
