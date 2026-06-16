use clap::{Args, Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "seqfind", about = "Sequence search and analysis CLI")]
pub struct Cli {
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
pub enum StatsFormat {
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

        /// Output format
        #[arg(long, value_enum, default_value_t = SearchFormat::Text)]
        format: SearchFormat,

        /// Persist results into PostgreSQL. URL comes from SEQFIND_DB_URL env,
        /// or pass --save-db <URL> to override.
        #[arg(long, num_args = 0..=1, default_missing_value = "")]
        save_db: Option<String>,
    },

    /// Print sequence statistics for a FASTA file
    Stats {
        /// Path to FASTA file
        #[arg(short, long)]
        file: PathBuf,

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

    /// Download a FASTA record from NCBI by accession (e.g. NC_045512.2)
    Fetch {
        /// NCBI Nucleotide accession (e.g. NC_000913.3 for E. coli)
        #[arg(short, long)]
        accession: String,

        /// Output FASTA path (default: data/<accession>.fna)
        #[arg(short, long)]
        out: Option<PathBuf>,
    },

    /// Inspect past runs stored in the database
    Runs {
        #[command(subcommand)]
        action: RunsAction,
    },
}

#[derive(Debug, Subcommand)]
pub enum RunsAction {
    /// List recent runs
    List {
        /// PostgreSQL connection URL (default: $SEQFIND_DB_URL)
        #[arg(long, env = crate::db::URL_ENV)]
        db: String,

        /// Max number of rows to show
        #[arg(long, default_value_t = 20)]
        limit: usize,
    },

    /// Show metadata and result counts for one run
    Show {
        /// PostgreSQL connection URL (default: $SEQFIND_DB_URL)
        #[arg(long, env = crate::db::URL_ENV)]
        db: String,

        /// Run id
        id: i64,
    },

    /// Re-emit a run's results in any supported format
    Export {
        /// PostgreSQL connection URL (default: $SEQFIND_DB_URL)
        #[arg(long, env = crate::db::URL_ENV)]
        db: String,

        /// Run id
        id: i64,

        /// Output format (uses the search/orf/stats format set for that kind)
        #[arg(long, value_enum, default_value_t = SearchFormat::Json)]
        format: SearchFormat,
    },
}
