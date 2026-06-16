# Example data

`sample.fasta` is a tiny hand-written file used by the unit tests.

For larger, real-world data, download reference genomes from NCBI into the
gitignored `data/` directory at the repo root. Two ways:

## Option A -direct `curl` (no extra tools)

```bash
mkdir -p data

# SARS-CoV-2 reference (NC_045512.2, ~30 KB)
curl -sS -o data/sars-cov-2.fna \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text"

# E. coli K-12 MG1655 complete genome (NC_000913.3, ~4.7 MB)
curl -sS -o data/ecoli.fna \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text"
```

## Option B -NCBI `datasets` CLI

Install once:

```bash
# macOS
brew install ncbi-datasets-cli
```

Then:

```bash
mkdir -p data && cd data

# SARS-CoV-2
datasets download genome accession GCF_009858895.2 --include genome
unzip -o ncbi_dataset.zip
mv ncbi_dataset/data/GCF_009858895.2/*.fna sars-cov-2.fna

# E. coli K-12 MG1655
datasets download genome accession GCF_000005845.2 --include genome
unzip -o ncbi_dataset.zip
mv ncbi_dataset/data/GCF_000005845.2/*.fna ecoli.fna
```

## Test commands

A few representative runs against the real data:

```bash
# Build the release binary first for realistic timing
cargo build --release
SF=./target/release/seqfind

# SARS-CoV-2 -sequence stats
$SF stats --file data/sars-cov-2.fna
# → NC_045512.2  length=29903  GC=38.0%  A=8954 T=9594 G=5863 C=5492 N=0

# SARS-CoV-2 -find ATG…TAA-like 9mers (start to stop) on both strands as BED
$SF search --file data/sars-cov-2.fna --pattern 'ATG???TAA' \
  --reverse-complement --format bed | head

# E. coli -full genome stats (4.64 Mb in ~40 ms on M-series Mac)
time $SF stats --file data/ecoli.fna

# E. coli -count ATGCGT hits on both strands (~2500 hits in ~50 ms)
time $SF search --file data/ecoli.fna --pattern ATGCGT \
  --reverse-complement --format bed | wc -l

# E. coli -wildcard search ATG???TAA
time $SF search --file data/ecoli.fna --pattern 'ATG???TAA' \
  --format bed | wc -l

# Pipe BED into other tools (sort, bedtools, etc.)
$SF search --file data/ecoli.fna --pattern GAATTC --format bed \
  | sort -k1,1 -k2,2n > data/ecoli_ecori_sites.bed
```

## Reference benchmarks

Recorded on an M-series Mac, `cargo build --release`, single thread:

| Run                                                       | Wall time |
| --------------------------------------------------------- | --------- |
| `stats` on E. coli (4.64 Mb)                              | ~40 ms    |
| `search ATGCGT --reverse-complement --format bed` E. coli | ~60 ms    |
| `search 'ATG???TAA' --format bed` E. coli                 | ~20 ms    |
| `stats` on SARS-CoV-2 (30 Kb)                             | <10 ms    |

Numbers will grow roughly linearly with file size and pattern count; the
matcher is a naive sliding window, which is plenty for MVP scale.
