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
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

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

struct CliTimer {
    label: &'static str,
    started: Instant,
}

impl CliTimer {
    fn start(label: &'static str) -> Self {
        let started = Instant::now();
        let started_at = SystemTime::now();
        eprintln!("{label} started: {}", format_system_time(started_at));
        Self { label, started }
    }

    fn finish(&self, success: bool) {
        let ended_at = SystemTime::now();
        let status = if success { "ok" } else { "failed" };
        eprintln!(
            "{label} ended:   {}",
            format_system_time(ended_at),
            label = self.label
        );
        eprintln!(
            "{label} took:    {} ({status})",
            format_duration(self.started.elapsed()),
            label = self.label
        );
    }
}

fn format_system_time(time: SystemTime) -> String {
    match time.duration_since(UNIX_EPOCH) {
        Ok(duration) => format_unix_timestamp(duration.as_secs()),
        Err(err) => format!(
            "{} before 1970-01-01 00:00:00 UTC",
            format_duration(err.duration())
        ),
    }
}

fn format_unix_timestamp(seconds: u64) -> String {
    let days = (seconds / 86_400) as i64;
    let seconds_of_day = seconds % 86_400;
    let (year, month, day) = civil_from_days(days);
    let hour = seconds_of_day / 3_600;
    let minute = (seconds_of_day % 3_600) / 60;
    let second = seconds_of_day % 60;

    format!("{year:04}-{month:02}-{day:02} {hour:02}:{minute:02}:{second:02} UTC")
}

fn civil_from_days(days_since_unix_epoch: i64) -> (i32, u32, u32) {
    let z = days_since_unix_epoch + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = mp + if mp < 10 { 3 } else { -9 };
    let year = y + if month <= 2 { 1 } else { 0 };

    (year as i32, month as u32, day as u32)
}

fn format_duration(duration: Duration) -> String {
    let total_seconds = duration.as_secs();
    let subsecond_millis = duration.subsec_millis();

    if total_seconds == 0 {
        let nanos = duration.subsec_nanos();
        if nanos < 1_000 {
            return format!("{nanos} ns");
        }
        if nanos < 1_000_000 {
            return format!("{:.3} us", nanos as f64 / 1_000.0);
        }
        return format!("{:.3} ms", nanos as f64 / 1_000_000.0);
    }

    if total_seconds < 60 {
        return format!(
            "{:.3} s",
            total_seconds as f64 + subsecond_millis as f64 / 1_000.0
        );
    }

    let hours = total_seconds / 3_600;
    let minutes = (total_seconds % 3_600) / 60;
    let seconds = total_seconds % 60;
    let seconds_with_fraction = seconds as f64 + subsecond_millis as f64 / 1_000.0;

    if hours == 0 {
        format!("{minutes}m {seconds_with_fraction:06.3}s")
    } else {
        format!("{hours}h {minutes:02}m {seconds_with_fraction:06.3}s")
    }
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

fn run_compare(
    reference: &Path,
    query_files: &[PathBuf],
    reverse_complement: bool,
    format: SearchFormat,
) -> Result<()> {
    if query_files.is_empty() {
        bail!("at least one --query FASTA file is required");
    }

    let reference_records = fasta::read_fasta(reference)?;
    let reference_file = reference.to_string_lossy().into_owned();
    let mut all_hits = Vec::new();

    for query_file in query_files {
        let query_records = fasta::read_fasta(query_file)?;
        let query_file_label = query_file.to_string_lossy().into_owned();

        for query_record in &query_records {
            if query_record.sequence.is_empty() {
                bail!(
                    "query record '{}' in {} is empty",
                    query_record.name,
                    query_file.display()
                );
            }

            for reference_record in &reference_records {
                let hits = matcher::find_matches(
                    &reference_record.name,
                    &reference_record.sequence,
                    &query_record.sequence,
                    reverse_complement,
                );

                for hit in hits {
                    all_hits.push(output::CompareHit {
                        reference_file: reference_file.clone(),
                        reference_record: hit.record,
                        query_file: query_file_label.clone(),
                        query_record: query_record.name.clone(),
                        start: hit.start,
                        end: hit.end,
                        strand: hit.strand,
                        matched: hit.matched,
                        query_length: query_record.sequence.len(),
                    });
                }
            }
        }
    }

    match format {
        SearchFormat::Text => output::print_compare_text(&all_hits),
        SearchFormat::Json => output::print_compare_json(&all_hits)?,
        SearchFormat::Csv => output::print_compare_csv(&all_hits)?,
        SearchFormat::Tsv => output::print_compare_tsv(&all_hits)?,
        SearchFormat::Bed => output::print_compare_bed(&all_hits)?,
    }

    Ok(())
}

fn command_label(command: &Command) -> &'static str {
    match command {
        Command::Search { .. } => "search",
        Command::Stats { .. } => "stats",
        Command::Rc { .. } => "rc",
        Command::Orf { .. } => "orf",
        Command::Motif { .. } => "motif",
        Command::Compare { .. } => "compare",
        Command::Fetch { .. } => "fetch",
        Command::Runs {
            action: RunsAction::List { .. },
        } => "runs list",
        Command::Runs {
            action: RunsAction::Show { .. },
        } => "runs show",
        Command::Runs {
            action: RunsAction::Export { .. },
        } => "runs export",
    }
}

fn main() -> Result<()> {
    reset_sigpipe();
    let cli = Cli::parse();
    let Cli { command, time: _ } = cli;
    let timer = CliTimer::start(command_label(&command));
    let result = run_command(command);
    timer.finish(result.is_ok());
    result
}

fn run_command(command: Command) -> Result<()> {
    match command {
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

        Command::Compare {
            reference,
            query,
            reverse_complement,
            format,
        } => run_compare(&reference, &query, reverse_complement, format)?,

        Command::Fetch {
            accession,
            out,
            timeout_seconds,
        } => {
            let out_path = out.unwrap_or_else(|| {
                PathBuf::from(format!("data/{}.fna", accession.replace('/', "_")))
            });
            let bytes =
                ncbi::fetch_fasta(&accession, &out_path, Duration::from_secs(timeout_seconds))?;
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
