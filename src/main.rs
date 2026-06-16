mod cli;
mod fasta;
mod matcher;
mod output;

use anyhow::Result;
use clap::Parser;

use cli::{Cli, Command, SearchFormat, StatsFormat};

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Command::Search {
            file,
            pattern,
            reverse_complement,
            format,
        } => {
            let records = fasta::read_fasta(&file)?;
            let mut all_hits = Vec::new();

            for record in &records {
                let hits = matcher::find_matches(
                    &record.name,
                    &record.sequence,
                    &pattern,
                    reverse_complement,
                );
                all_hits.extend(hits);
            }

            match format {
                SearchFormat::Text => output::print_hits_text(&all_hits),
                SearchFormat::Json => output::print_hits_json(&all_hits)?,
                SearchFormat::Csv => output::print_hits_csv(&all_hits)?,
                SearchFormat::Tsv => output::print_hits_tsv(&all_hits)?,
                SearchFormat::Bed => output::print_hits_bed(&all_hits)?,
            }
        }

        Command::Stats { file, format } => {
            let records = fasta::read_fasta(&file)?;
            let stats: Vec<_> = records.iter().map(output::compute_stats).collect();

            match format {
                StatsFormat::Text => output::print_stats_text(&stats),
                StatsFormat::Json => output::print_stats_json(&stats)?,
                StatsFormat::Csv => output::print_stats_csv(&stats)?,
                StatsFormat::Tsv => output::print_stats_tsv(&stats)?,
            }
        }

        Command::Rc { sequence } => {
            println!("{}", matcher::reverse_complement(&sequence.to_uppercase()));
        }
    }

    Ok(())
}
