use clap::{Parser, Subcommand, ValueEnum};
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
    },

    /// Print sequence statistics for a FASTA file
    Stats {
        /// Path to FASTA file
        #[arg(short, long)]
        file: PathBuf,

        /// Output format
        #[arg(long, value_enum, default_value_t = StatsFormat::Text)]
        format: StatsFormat,
    },

    /// Print reverse complement of a sequence
    Rc {
        /// DNA sequence to reverse-complement
        sequence: String,
    },
}
