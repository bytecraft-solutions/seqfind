use clap::{Args, Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "seqfind", about = "Sequence search and analysis CLI")]
pub struct Cli {
    /// Print timing information. Timing is currently always printed for CLI operations.
    #[arg(long, global = true)]
    pub time: bool,

    #[command(subcommand)]
    pub command: Command,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum SearchFormat {
    Text,
    Json,
    Csv,
    Tsv,
    Bed,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum CompareMode {
    /// Search query FASTA records as subsequences inside reference records.
    Exact,

    /// Compare paired records and emit SNP/indel-style variants.
    Variants,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum StatsFormat {
    Text,
    Json,
    Csv,
    Tsv,
}

/// Output formats for record-emitting subcommands (`seq`, `convert`, `sample`,
/// `head`, `rmdup`, `sort`, `shuffle`, `split`, `subseq`).
#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum RecordFormat {
    /// Detect from input format and `--out` extension; falls back to FASTA.
    Auto,
    Fasta,
    Fastq,
    /// Tab-separated (id, sequence, optional quality).
    Tab,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum SortKeyArg {
    Id,
    Length,
    Sequence,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum RmdupKeyArg {
    Id,
    Sequence,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum LocateFormat {
    Text,
    Json,
    Csv,
    Tsv,
    Bed,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum DigestFormat {
    Text,
    Json,
    Csv,
    Tsv,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum OrfFormat {
    Text,
    Json,
    Csv,
    Tsv,
    Bed,
}

#[derive(Debug, Args)]
#[group(required = true, multiple = false)]
pub struct RcInput {
    /// DNA sequence to reverse-complement
    pub sequence: Option<String>,

    /// Reverse-complement every record in a FASTA file, output as FASTA
    #[arg(short, long)]
    pub file: Option<PathBuf>,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    /// Search for a pattern in FASTA sequences
    Search {
        /// Path to FASTA file
        #[arg(short, long)]
        file: PathBuf,

        /// Pattern to search (use ? as single-base wildcard)
        #[arg(short, long)]
        pattern: String,

        /// Also search the reverse complement strand
        #[arg(long)]
        reverse_complement: bool,

        /// Allow up to N substitutions/insertions/deletions using Rust-Bio Myers
        #[arg(long)]
        max_edits: Option<usize>,

        /// Output format
        #[arg(long, value_enum, default_value_t = SearchFormat::Text)]
        format: SearchFormat,

        /// Persist results into PostgreSQL. URL comes from SEQFIND_DB_URL env,
        /// or pass --save-db <URL> to override.
        #[arg(long, num_args = 0..=1, default_missing_value = "")]
        save_db: Option<String>,
    },

    /// Print sequence statistics for a FASTA/FASTQ file
    Stats {
        /// Path to FASTA/FASTQ file (use `-` for stdin)
        #[arg(short, long)]
        file: PathBuf,

        /// Print a SeqKit-style aggregate summary (N50, sum/min/avg/max, GC%,
        /// plus Q20/Q30 for FASTQ) instead of per-record rows.
        #[arg(long)]
        summary: bool,

        /// Output format
        #[arg(long, value_enum, default_value_t = StatsFormat::Text)]
        format: StatsFormat,

        /// Persist results into PostgreSQL. URL comes from SEQFIND_DB_URL env,
        /// or pass --save-db <URL> to override.
        #[arg(long, num_args = 0..=1, default_missing_value = "")]
        save_db: Option<String>,
    },

    /// Reverse complement: pass a sequence inline, or `--file` for whole-record FASTA
    Rc {
        #[command(flatten)]
        input: RcInput,
    },

    /// Find open reading frames (ATG → stop) on both strands
    Orf {
        /// Path to FASTA file
        #[arg(short, long)]
        file: PathBuf,

        /// Minimum ORF length in bases (includes start and stop codons)
        #[arg(long, default_value_t = 100)]
        min_length: usize,

        /// Only scan the forward strand
        #[arg(long)]
        forward_only: bool,

        /// Output format
        #[arg(long, value_enum, default_value_t = OrfFormat::Text)]
        format: OrfFormat,

        /// Persist results into PostgreSQL. URL comes from SEQFIND_DB_URL env,
        /// or pass --save-db <URL> to override.
        #[arg(long, num_args = 0..=1, default_missing_value = "")]
        save_db: Option<String>,
    },

    /// Search a named motif preset (run `seqfind motif --list` to see them)
    Motif {
        /// Path to FASTA file (omit when using --list)
        #[arg(short, long)]
        file: Option<PathBuf>,

        /// Preset name, e.g. EcoRI, TATA, Kozak (case-insensitive)
        #[arg(short, long)]
        name: Option<String>,

        /// Also search the reverse complement strand
        #[arg(long)]
        reverse_complement: bool,

        /// Output format
        #[arg(long, value_enum, default_value_t = SearchFormat::Text)]
        format: SearchFormat,

        /// List the available motif presets and exit
        #[arg(long)]
        list: bool,

        /// Persist results into PostgreSQL. URL comes from SEQFIND_DB_URL env,
        /// or pass --save-db <URL> to override.
        #[arg(long, num_args = 0..=1, default_missing_value = "")]
        save_db: Option<String>,
    },

    /// Search one or more query FASTA files inside a reference FASTA file
    Compare {
        /// Reference FASTA file to search within
        #[arg(short, long)]
        reference: PathBuf,

        /// Query FASTA file(s) to search for. Repeat --query for multiple files.
        #[arg(short, long, required = true, num_args = 1..)]
        query: Vec<PathBuf>,

        /// Also search the reverse complement strand
        #[arg(long)]
        reverse_complement: bool,

        /// Allow up to N substitutions/insertions/deletions for containment search
        #[arg(long)]
        max_edits: Option<usize>,

        /// Compare mode: exact containment search or paired variant calling
        #[arg(long, value_enum, default_value_t = CompareMode::Exact)]
        mode: CompareMode,

        /// Build an FM-index per reference record (exact mode only). Pays off
        /// when many queries are searched against a chromosome-scale reference.
        #[arg(long)]
        fm_index: bool,

        /// In variant mode, pair records by file order instead of FASTA record ID
        #[arg(long)]
        pair_by_order: bool,

        /// Output format
        #[arg(long, value_enum, default_value_t = SearchFormat::Text)]
        format: SearchFormat,
    },

    /// Download a FASTA record from NCBI by accession (e.g. NC_045512.2)
    Fetch {
        /// NCBI Nucleotide accession (e.g. NC_000913.3 for E. coli)
        #[arg(short, long)]
        accession: String,

        /// Output FASTA path (default: data/<accession>.fna)
        #[arg(short, long)]
        out: Option<PathBuf>,

        /// Read timeout in seconds for large NCBI downloads
        #[arg(long, default_value_t = 1_800)]
        timeout_seconds: u64,
    },

    /// Inspect past runs stored in the database
    Runs {
        #[command(subcommand)]
        action: RunsAction,
    },

    /// Apply SeqKit-style per-record transforms (case, RNA↔DNA, length filters, RC, gap removal)
    Seq {
        /// Path to FASTA/FASTQ file (use `-` for stdin)
        #[arg(short, long)]
        file: PathBuf,
        /// Output path (omit or use `-` for stdout)
        #[arg(short, long)]
        out: Option<PathBuf>,
        #[arg(long)]
        upper: bool,
        #[arg(long)]
        lower: bool,
        /// Remove `-` and `.` from sequences (alignment gaps).
        #[arg(long)]
        remove_gaps: bool,
        /// Replace T with U (DNA → RNA).
        #[arg(long)]
        rna: bool,
        /// Replace U with T (RNA → DNA).
        #[arg(long)]
        dna: bool,
        #[arg(long)]
        reverse: bool,
        #[arg(long)]
        complement: bool,
        /// Reverse + complement.
        #[arg(long)]
        reverse_complement: bool,
        /// Drop records shorter than N bases.
        #[arg(long)]
        min_length: Option<usize>,
        /// Drop records longer than N bases.
        #[arg(long)]
        max_length: Option<usize>,
        /// Print only the record IDs.
        #[arg(long, conflicts_with_all = ["only_seq", "only_name"])]
        only_id: bool,
        /// Print only the FASTA names (id + description).
        #[arg(long, conflicts_with_all = ["only_seq", "only_id"])]
        only_name: bool,
        /// Print only the sequences.
        #[arg(long, conflicts_with_all = ["only_id", "only_name"])]
        only_seq: bool,
        /// FASTA line wrap width (0 = no wrap). Ignored for FASTQ output.
        #[arg(long, default_value_t = 70)]
        wrap: usize,
        /// Output format
        #[arg(long, value_enum, default_value_t = RecordFormat::Auto)]
        format: RecordFormat,
    },

    /// Convert between formats: FASTQ→FASTA, fx2tab, tab2fx.
    Convert {
        #[command(subcommand)]
        action: ConvertAction,
    },

    /// Take the first N records (SeqKit's `head`).
    Head {
        #[arg(short, long)]
        file: PathBuf,
        #[arg(short, long, default_value_t = 10)]
        number: usize,
        #[arg(short, long)]
        out: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = RecordFormat::Auto)]
        format: RecordFormat,
    },

    /// Take a 1-based inclusive record range (SeqKit's `range`).
    Range {
        #[arg(short, long)]
        file: PathBuf,
        /// `START:END` (1-based inclusive). Negative values count from the end.
        #[arg(short, long)]
        range: String,
        #[arg(short, long)]
        out: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = RecordFormat::Auto)]
        format: RecordFormat,
    },

    /// Sample N records uniformly without replacement, with a deterministic seed.
    Sample {
        #[arg(short, long)]
        file: PathBuf,
        /// Number of records to sample (mutually exclusive with --proportion).
        #[arg(short, long, conflicts_with = "proportion")]
        number: Option<usize>,
        /// Sampling proportion 0.0..=1.0.
        #[arg(short, long, conflicts_with = "number")]
        proportion: Option<f64>,
        /// PRNG seed (default 11).
        #[arg(short, long, default_value_t = 11)]
        seed: u64,
        #[arg(short, long)]
        out: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = RecordFormat::Auto)]
        format: RecordFormat,
    },

    /// Shuffle records with a deterministic seed (SeqKit's `shuffle`).
    Shuffle {
        #[arg(short, long)]
        file: PathBuf,
        #[arg(short, long, default_value_t = 23)]
        seed: u64,
        #[arg(short, long)]
        out: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = RecordFormat::Auto)]
        format: RecordFormat,
    },

    /// Drop duplicate records by ID or sequence (SeqKit's `rmdup`).
    Rmdup {
        #[arg(short, long)]
        file: PathBuf,
        #[arg(long, value_enum, default_value_t = RmdupKeyArg::Id)]
        by: RmdupKeyArg,
        #[arg(short, long)]
        out: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = RecordFormat::Auto)]
        format: RecordFormat,
    },

    /// Sort records by id, length, or sequence (SeqKit's `sort`).
    Sort {
        #[arg(short, long)]
        file: PathBuf,
        #[arg(long, value_enum, default_value_t = SortKeyArg::Id)]
        by: SortKeyArg,
        #[arg(long)]
        reverse: bool,
        #[arg(short, long)]
        out: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = RecordFormat::Auto)]
        format: RecordFormat,
    },

    /// Split records into N parts or chunks of size K (SeqKit's `split`).
    Split {
        #[arg(short, long)]
        file: PathBuf,
        /// Split into this many parts (mutually exclusive with --size).
        #[arg(long, conflicts_with = "size")]
        parts: Option<usize>,
        /// Split into chunks of this many records.
        #[arg(long, conflicts_with = "parts")]
        size: Option<usize>,
        /// Directory to write chunks into.
        #[arg(short = 'O', long, default_value = ".")]
        out_dir: PathBuf,
        /// Filename prefix for the chunks (default: input stem).
        #[arg(long)]
        prefix: Option<String>,
        #[arg(long, value_enum, default_value_t = RecordFormat::Auto)]
        format: RecordFormat,
    },

    /// Locate motifs/regex/IUPAC patterns inside records (SeqKit's `locate`).
    Locate {
        #[arg(short, long)]
        file: PathBuf,
        /// One or more patterns. Use `--regex` to disable IUPAC translation.
        #[arg(short, long, num_args = 1..)]
        pattern: Vec<String>,
        /// Path to a file with one pattern per line (use with or instead of --pattern).
        #[arg(long)]
        pattern_file: Option<PathBuf>,
        /// Treat patterns as raw regex instead of IUPAC ambiguity codes.
        #[arg(long)]
        regex: bool,
        /// Also locate on the reverse-complement strand.
        #[arg(long)]
        reverse_complement: bool,
        #[arg(long, value_enum, default_value_t = LocateFormat::Text)]
        format: LocateFormat,
    },

    /// Extract a 1-based inclusive region from each record (SeqKit's `subseq --region`).
    Subseq {
        #[arg(short, long)]
        file: PathBuf,
        /// `START:END` (1-based inclusive). Negatives count from the end.
        /// Required unless one of --bed/--gff/--gtf is set.
        #[arg(short, long, conflicts_with_all = ["bed", "gff", "gtf"])]
        region: Option<String>,
        /// BED region file (chrom\tstart\tend; 0-based half-open).
        #[arg(long, conflicts_with_all = ["gff", "gtf"])]
        bed: Option<PathBuf>,
        /// GFF3 region file (1-based inclusive).
        #[arg(long, conflicts_with_all = ["bed", "gtf"])]
        gff: Option<PathBuf>,
        /// GTF region file (1-based inclusive).
        #[arg(long, conflicts_with_all = ["bed", "gff"])]
        gtf: Option<PathBuf>,
        /// Bases of upstream context to include.
        #[arg(long, default_value_t = 0)]
        up_stream: usize,
        /// Bases of downstream context to include.
        #[arg(long, default_value_t = 0)]
        down_stream: usize,
        /// Take the region from the reverse-complement strand.
        #[arg(long)]
        reverse_complement: bool,
        #[arg(short, long)]
        out: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = RecordFormat::Auto)]
        format: RecordFormat,
    },

    /// Print per-record SHA-256 digests (case-insensitive) -SeqKit's `sum`.
    Sum {
        #[arg(short, long)]
        file: PathBuf,
        #[arg(long, value_enum, default_value_t = DigestFormat::Text)]
        format: DigestFormat,
    },

    /// Regex-replace record IDs or sequences (SeqKit's `replace`).
    Replace {
        #[arg(short, long)]
        file: PathBuf,
        /// Regex pattern to match.
        #[arg(short, long)]
        pattern: String,
        /// Replacement string. Supports `$1`-style backreferences.
        #[arg(short, long)]
        replacement: String,
        /// Operate on sequences instead of IDs.
        #[arg(long)]
        by_seq: bool,
        #[arg(short, long)]
        out: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = RecordFormat::Auto)]
        format: RecordFormat,
    },

    /// Disambiguate duplicate record IDs by appending `_2`, `_3`, … (SeqKit's `rename`).
    Rename {
        #[arg(short, long)]
        file: PathBuf,
        #[arg(short, long)]
        out: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = RecordFormat::Auto)]
        format: RecordFormat,
    },

    /// Concatenate corresponding records from multiple files (SeqKit's `concat`).
    Concat {
        /// Two or more input files. Records are paired by ID by default.
        /// Repeat the flag (`--file a --file b`) or pass several values
        /// (`--file a b`).
        #[arg(short, long, required = true, num_args = 1.., action = clap::ArgAction::Append)]
        file: Vec<PathBuf>,
        /// Pair records by index order instead of ID.
        #[arg(long)]
        by_order: bool,
        #[arg(short, long)]
        out: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = RecordFormat::Auto)]
        format: RecordFormat,
    },

    /// Rotate circular genomes so each record starts at 1-based POSITION (SeqKit's `restart`).
    Restart {
        #[arg(short, long)]
        file: PathBuf,
        /// 1-based position to rotate to the start of the record.
        #[arg(short, long)]
        position: usize,
        #[arg(short, long)]
        out: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = RecordFormat::Auto)]
        format: RecordFormat,
    },

    /// Apply explicit point / deletion / insertion mutations (SeqKit's `mutate`).
    Mutate {
        #[arg(short, long)]
        file: PathBuf,
        /// Point mutation `POS:BASE` (1-based). Repeat for multiple.
        #[arg(long = "point")]
        point: Vec<String>,
        /// Deletion `START:END` (1-based inclusive). Repeat for multiple.
        #[arg(long = "deletion")]
        deletion: Vec<String>,
        /// Insertion `POS:BASES` (1-based). Repeat for multiple.
        #[arg(long = "insertion")]
        insertion: Vec<String>,
        #[arg(short, long)]
        out: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = RecordFormat::Auto)]
        format: RecordFormat,
    },

    /// Extract amplicons defined by a forward+reverse primer pair (SeqKit's `amplicon`).
    Amplicon {
        #[arg(short, long)]
        file: PathBuf,
        /// Forward primer (5'→3'). IUPAC codes allowed.
        #[arg(short = 'F', long)]
        forward: String,
        /// Reverse primer (5'→3'). Its RC is searched on the forward strand.
        #[arg(short = 'R', long)]
        reverse: String,
        /// Scan the reverse strand as well.
        #[arg(long)]
        both_strands: bool,
        /// Drop amplicons longer than this many bases.
        #[arg(long)]
        max_length: Option<usize>,
        /// Output extracted amplicons as records.
        #[arg(short, long)]
        out: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = RecordFormat::Auto)]
        format: RecordFormat,
    },

    /// Build a .fai index and fetch random-access regions (SeqKit's `faidx`).
    Faidx {
        #[arg(short, long)]
        file: PathBuf,
        /// One or more `chr` or `chr:start-end` (1-based inclusive) regions.
        /// If omitted, just builds/refreshes the .fai index.
        #[arg()]
        region: Vec<String>,
        #[arg(short, long)]
        out: Option<PathBuf>,
        /// Wrap width for output FASTA. 0 disables wrapping.
        #[arg(long, default_value_t = 70)]
        wrap: usize,
        /// Always rebuild the .fai index even if one exists.
        #[arg(long)]
        force: bool,
    },
}

#[derive(Debug, Subcommand)]
pub enum ConvertAction {
    /// FASTQ → FASTA (drop quality).
    Fq2fa {
        #[arg(short, long)]
        file: PathBuf,
        #[arg(short, long)]
        out: Option<PathBuf>,
        #[arg(long, default_value_t = 70)]
        wrap: usize,
    },
    /// Records → TSV (`id\tsequence[\tquality]`).
    Fx2tab {
        #[arg(short, long)]
        file: PathBuf,
        #[arg(short, long)]
        out: Option<PathBuf>,
    },
    /// TSV → records (inverse of `fx2tab`).
    Tab2fx {
        #[arg(short, long)]
        file: PathBuf,
        #[arg(short, long)]
        out: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = RecordFormat::Auto)]
        format: RecordFormat,
    },
}

#[derive(Debug, Subcommand)]
pub enum RunsAction {
    /// List recent runs
    List {
        /// PostgreSQL connection URL (default: $SEQFIND_DB_URL)
        #[arg(long, env = seqfind_core::db::URL_ENV)]
        db: String,

        /// Max number of rows to show
        #[arg(long, default_value_t = 20)]
        limit: usize,
    },

    /// Show metadata and result counts for one run
    Show {
        /// PostgreSQL connection URL (default: $SEQFIND_DB_URL)
        #[arg(long, env = seqfind_core::db::URL_ENV)]
        db: String,

        /// Run id
        id: i64,
    },

    /// Re-emit a run's results in any supported format
    Export {
        /// PostgreSQL connection URL (default: $SEQFIND_DB_URL)
        #[arg(long, env = seqfind_core::db::URL_ENV)]
        db: String,

        /// Run id
        id: i64,

        /// Output format (uses the search/orf/stats format set for that kind)
        #[arg(long, value_enum, default_value_t = SearchFormat::Json)]
        format: SearchFormat,
    },
}
