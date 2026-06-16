mod cli;
mod db;
mod fasta;
mod matcher;
mod motif;
mod ncbi;
mod orf;
mod output;

use anyhow::{Result, bail};
use clap::Parser;
use serde_json::json;
use std::path::{Path, PathBuf};

use cli::{Cli, Command, OrfFormat, RunsAction, SearchFormat, StatsFormat};

/// Restore the default SIGPIPE handler so that piping output to `head`, `less`,
/// etc. terminates the process cleanly instead of panicking inside `println!`.
#[cfg(unix)]
fn reset_sigpipe() {
    // SAFETY: setting SIGPIPE to SIG_DFL is the standard, well-defined Unix
    // behavior for CLI tools that produce streaming output.
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }
}

#[cfg(not(unix))]
fn reset_sigpipe() {}

fn maybe_save_hits(
    save_db: &Option<String>,
    subcommand: &str,
    input_file: &Path,
    params: serde_json::Value,
    hits: &[matcher::Hit],
) -> Result<()> {
    let Some(override_url) = save_db else {
        return Ok(());
    };
    let url = db::resolve_url(Some(override_url.as_str()))?;
    let mut db = db::Database::open(&url)?;
    let run_id = db.record_run(
        subcommand,
        Some(&input_file.to_string_lossy()),
        &params.to_string(),
    )?;
    db.insert_hits(run_id, hits)?;
    eprintln!("✓ saved run #{run_id} ({} hits) to postgres", hits.len());
    Ok(())
}

fn run_search(
    file: &Path,
    pattern: &str,
    reverse_complement: bool,
    format: SearchFormat,
    save_db: &Option<String>,
    subcommand_label: &str,
) -> Result<()> {
    let records = fasta::read_fasta(file)?;
    let mut all_hits = Vec::new();
    for record in &records {
        let hits =
            matcher::find_matches(&record.name, &record.sequence, pattern, reverse_complement);
        all_hits.extend(hits);
    }
    match format {
        SearchFormat::Text => output::print_hits_text(&all_hits),
        SearchFormat::Json => output::print_hits_json(&all_hits)?,
        SearchFormat::Csv => output::print_hits_csv(&all_hits)?,
        SearchFormat::Tsv => output::print_hits_tsv(&all_hits)?,
        SearchFormat::Bed => output::print_hits_bed(&all_hits)?,
    }
    let params = json!({
        "pattern": pattern,
        "reverse_complement": reverse_complement,
        "format": format!("{format:?}").to_lowercase(),
    });
    maybe_save_hits(save_db, subcommand_label, file, params, &all_hits)
}

fn main() -> Result<()> {
    reset_sigpipe();
    let cli = Cli::parse();

    match cli.command {
        Command::Search {
            file,
            pattern,
            reverse_complement,
            format,
            save_db,
        } => run_search(
            &file,
            &pattern,
            reverse_complement,
            format,
            &save_db,
            "search",
        )?,

        Command::Stats {
            file,
            format,
            save_db,
        } => {
            let records = fasta::read_fasta(&file)?;
            let stats: Vec<_> = records.iter().map(output::compute_stats).collect();

            match format {
                StatsFormat::Text => output::print_stats_text(&stats),
                StatsFormat::Json => output::print_stats_json(&stats)?,
                StatsFormat::Csv => output::print_stats_csv(&stats)?,
                StatsFormat::Tsv => output::print_stats_tsv(&stats)?,
            }

            if let Some(override_url) = &save_db {
                let url = db::resolve_url(Some(override_url.as_str()))?;
                let mut db = db::Database::open(&url)?;
                let params = json!({ "format": format!("{format:?}").to_lowercase() });
                let run_id =
                    db.record_run("stats", Some(&file.to_string_lossy()), &params.to_string())?;
                db.insert_stats(run_id, &stats)?;
                eprintln!(
                    "✓ saved run #{run_id} ({} records) to postgres",
                    stats.len()
                );
            }
        }

        Command::Rc { input } => match (input.sequence, input.file) {
            (Some(seq), None) => {
                println!("{}", matcher::reverse_complement(&seq.to_uppercase()));
            }
            (None, Some(path)) => {
                let records = fasta::read_fasta(&path)?;
                for r in &records {
                    let rc = matcher::reverse_complement(&r.sequence);
                    let name = format!("{}_rc", r.name);
                    output::print_fasta_record(&name, &rc, 70);
                }
            }
            _ => unreachable!("clap arg group guarantees exactly one"),
        },

        Command::Orf {
            file,
            min_length,
            forward_only,
            format,
            save_db,
        } => {
            let records = fasta::read_fasta(&file)?;
            let mut all_orfs = Vec::new();
            for r in &records {
                let orfs = orf::find_orfs(&r.name, &r.sequence, min_length, !forward_only);
                all_orfs.extend(orfs);
            }
            match format {
                OrfFormat::Text => output::print_orfs_text(&all_orfs),
                OrfFormat::Json => output::print_orfs_json(&all_orfs)?,
                OrfFormat::Csv => output::print_orfs_csv(&all_orfs)?,
                OrfFormat::Tsv => output::print_orfs_tsv(&all_orfs)?,
                OrfFormat::Bed => output::print_orfs_bed(&all_orfs)?,
            }

            if let Some(override_url) = &save_db {
                let url = db::resolve_url(Some(override_url.as_str()))?;
                let mut db = db::Database::open(&url)?;
                let params = json!({
                    "min_length": min_length,
                    "forward_only": forward_only,
                    "format": format!("{format:?}").to_lowercase(),
                });
                let run_id =
                    db.record_run("orf", Some(&file.to_string_lossy()), &params.to_string())?;
                db.insert_orfs(run_id, &all_orfs)?;
                eprintln!(
                    "✓ saved run #{run_id} ({} ORFs) to postgres",
                    all_orfs.len()
                );
            }
        }

        Command::Motif {
            file,
            name,
            reverse_complement,
            format,
            list,
            save_db,
        } => {
            if list {
                println!("Available motif presets:");
                for m in motif::PRESETS {
                    println!("  {:<14} {:<12} {}", m.name, m.pattern, m.description);
                }
                return Ok(());
            }
            let Some(name) = name else {
                bail!("--name is required unless --list is set");
            };
            let Some(file) = file else {
                bail!("--file is required when searching a motif");
            };
            let Some(preset) = motif::lookup(&name) else {
                bail!(
                    "unknown motif '{}'. Run `seqfind motif --list` to see presets.",
                    name
                );
            };
            run_search(
                &file,
                preset.pattern,
                reverse_complement,
                format,
                &save_db,
                "motif",
            )?;
        }

        Command::Fetch { accession, out } => {
            let out_path = out.unwrap_or_else(|| {
                PathBuf::from(format!("data/{}.fna", accession.replace('/', "_")))
            });
            let bytes = ncbi::fetch_fasta(&accession, &out_path)?;
            eprintln!(
                "✓ fetched {} ({} bytes) → {}",
                accession,
                bytes,
                out_path.display()
            );
        }

        Command::Runs { action } => match action {
            RunsAction::List { db, limit } => {
                let mut database = db::Database::open(&db)?;
                let runs = database.list_runs(limit)?;
                if runs.is_empty() {
                    println!(
                        "No runs yet. Use --save-db <path> on search/orf/stats to record one."
                    );
                    return Ok(());
                }
                println!(
                    "{:>4}  {:<19}  {:<10}  {:>8}  input",
                    "id", "ts (UTC)", "kind", "rows"
                );
                for r in runs {
                    let input = r.input_file.as_deref().unwrap_or("-");
                    println!(
                        "{:>4}  {:<19}  {:<10}  {:>8}  {}",
                        r.id, r.ts, r.subcommand, r.result_count, input
                    );
                }
            }
            RunsAction::Show { db, id } => {
                let mut database = db::Database::open(&db)?;
                let Some(meta) = database.get_run(id)? else {
                    bail!("no run with id {id}");
                };
                println!("id:        {}", meta.id);
                println!("ts (UTC):  {}", meta.ts);
                println!("kind:      {}", meta.subcommand);
                println!("input:     {}", meta.input_file.as_deref().unwrap_or("-"));
                println!("params:    {}", meta.params);
                let hits = database.get_hits(id)?.len();
                let orfs = database.get_orfs(id)?.len();
                let stats = database.get_stats(id)?.len();
                println!("results:   hits={hits} orfs={orfs} stats={stats}");
            }
            RunsAction::Export { db, id, format } => {
                let mut database = db::Database::open(&db)?;
                let Some(meta) = database.get_run(id)? else {
                    bail!("no run with id {id}");
                };
                match meta.subcommand.as_str() {
                    "search" | "motif" => {
                        let hits = database.get_hits(id)?;
                        match format {
                            SearchFormat::Text => output::print_hits_text(&hits),
                            SearchFormat::Json => output::print_hits_json(&hits)?,
                            SearchFormat::Csv => output::print_hits_csv(&hits)?,
                            SearchFormat::Tsv => output::print_hits_tsv(&hits)?,
                            SearchFormat::Bed => output::print_hits_bed(&hits)?,
                        }
                    }
                    "orf" => {
                        let orfs = database.get_orfs(id)?;
                        match format {
                            SearchFormat::Text => output::print_orfs_text(&orfs),
                            SearchFormat::Json => output::print_orfs_json(&orfs)?,
                            SearchFormat::Csv => output::print_orfs_csv(&orfs)?,
                            SearchFormat::Tsv => output::print_orfs_tsv(&orfs)?,
                            SearchFormat::Bed => output::print_orfs_bed(&orfs)?,
                        }
                    }
                    "stats" => {
                        let stats = database.get_stats(id)?;
                        match format {
                            SearchFormat::Text => output::print_stats_text(&stats),
                            SearchFormat::Json => output::print_stats_json(&stats)?,
                            SearchFormat::Csv => output::print_stats_csv(&stats)?,
                            SearchFormat::Tsv => output::print_stats_tsv(&stats)?,
                            SearchFormat::Bed => bail!("BED is not supported for stats runs"),
                        }
                    }
                    other => bail!("unknown subcommand recorded for run {id}: {other}"),
                }
            }
        },
    }

    Ok(())
}
