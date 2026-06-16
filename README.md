# seqfind

`seqfind` is a Rust CLI for deterministic sequence search and analysis on
FASTA files. It is built as the core of a future bioinformatics agent: the
binary performs the actual work, and later AI/MCP layers can route requests to
these tested operations.

## Current Features

- Fetch FASTA records from NCBI Nucleotide with streamed writes for large files.
- Read FASTA using `noodles-fasta`.
- Exact sequence search using Rust-Bio's Horspool matcher.
- Approximate sequence search using Rust-Bio's Myers matcher via
  `--max-edits`.
- Wildcard sequence search with `?` as a single-base wildcard.
- Reverse-complement search and whole-record reverse-complement output.
- Sequence stats: length, GC%, and base counts.
- ORF search across frames and strands.
- Named motif presets such as EcoRI, BamHI, HindIII, TATA, Kozak, and PolyA.
- Compare one reference FASTA against one or more query FASTA files.
- Variant-style FASTA comparison using Rust-Bio pairwise alignment via
  `compare --mode variants`.
- Output formats: text, JSON, CSV, TSV, and BED where applicable.
- PostgreSQL run history with `runs list`, `runs show`, and `runs export`.
- Timing metadata on CLI operations.

## Build

```bash
cargo build --release
./target/release/seqfind --help
```

For local development:

```bash
cargo run -- --help
cargo test
```

Install the current workspace binary:

```bash
cargo install --path . --locked
```

## Quick Start

Fetch a small reference:

```bash
seqfind fetch --accession NC_045512.2 --out data/sars-cov-2.fna
```

For larger chromosomes, increase the read timeout if needed:

```bash
seqfind fetch --accession NC_000001.11 --out data/chr1.fna --timeout-seconds 3600
```

Compute stats:

```bash
seqfind stats --file data/sars-cov-2.fna --format json --time
```

Search an exact pattern:

```bash
seqfind search \
  --file data/sars-cov-2.fna \
  --pattern ATGCGT \
  --reverse-complement \
  --format bed \
  --time
```

Search with wildcards:

```bash
seqfind search \
  --file data/sars-cov-2.fna \
  --pattern 'ATG???TAA' \
  --format json
```

Search with up to one edit:

```bash
seqfind search \
  --file data/sars-cov-2.fna \
  --pattern ATGCGT \
  --max-edits 1 \
  --format json
```

Find ORFs:

```bash
seqfind orf \
  --file data/sars-cov-2.fna \
  --min-length 300 \
  --format bed
```

Use motif presets:

```bash
seqfind motif --list
seqfind motif --file data/sars-cov-2.fna --name EcoRI --format bed
```

Compare a reference against multiple query FASTA files:

```bash
seqfind compare \
  --reference data/chr1.fna \
  --query data/query_marker_1.fna \
  --query data/query_marker_2.fna \
  --reverse-complement \
  --format json \
  --time
```

By default, `compare` performs containment search: every query FASTA record is
searched inside every reference FASTA record. Use `--max-edits <N>` for mutated
markers/contigs:

```bash
seqfind compare \
  --reference examples/compare_batch/ref.fna \
  --query examples/compare_batch/ecoli_marker.fna \
  --max-edits 1 \
  --format json \
  --time
```

For SNP/indel-style record-to-record comparison, use variant mode. Records are
paired by FASTA ID by default:

```bash
seqfind compare \
  --reference examples/demo_ref.fna \
  --query examples/demo_sample_mut.fna \
  --mode variants \
  --format json \
  --time
```

If the files have different record names but the same order, add
`--pair-by-order`. Variant mode uses global pairwise alignment, so use it for
gene/contig-sized paired records rather than whole chromosomes.

## PostgreSQL Run History

Set a database URL:

```bash
export SEQFIND_DB_URL=postgresql://postgres:passw0rd@127.0.0.1/seqfind
```

Persist supported results:

```bash
seqfind search --file data/sars-cov-2.fna --pattern GAATTC --save-db
seqfind runs list
seqfind runs show 1
seqfind runs export 1 --format json
```

## Implementation Notes

Current core crates:

- `clap` for CLI parsing.
- `serde` and `serde_json` for structured output.
- `csv` for CSV and TSV output.
- `noodles-fasta` for FASTA parsing.
- `bio` for sequence matching algorithms.
- `ureq` for NCBI EFetch.
- `postgres` for run persistence.

## Roadmap

Near-term algorithm work:

- Rust-Bio BNDM or Shift-And for short motif/primer matching.
- Rust-Bio FM-index or suffix-array based repeated search against large references.
- `noodles-fasta` FAI support for indexed region queries.
- `noodles-gff`, `noodles-gtf`, and `noodles-bed` for annotation and interval operations.
- `noodles-vcf` for variant output once variant compare lands.

Agent/MCP work:

- Natural-language `ask` command that maps prompts to validated task JSON.
- Reproducible run folders with inputs, plan, command, outputs, and logs.
- Local MCP server exposing deterministic tools such as `search_exact`,
  `search_approx`, `find_orfs`, `read_gff_features`, and `write_bed`.
