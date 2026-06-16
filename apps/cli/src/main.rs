mod cli;

#[allow(unused_imports)]
use seqfind_core::{
    alphabets, convert as core_convert, db, digest, edit as core_edit, index, locate, matcher,
    motif, ncbi, orf, output, select, seq, subseq as core_subseq, summary, transform,
};

use anyhow::{Result, bail};
use bio::alignment::AlignmentOperation;
use bio::alignment::pairwise::Aligner;
use clap::Parser;
use serde_json::json;
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use cli::{
    Cli, Command, CompareMode, ConvertAction, DigestFormat, LocateFormat, OrfFormat, RecordFormat,
    RmdupKeyArg, RunsAction, SearchFormat, SortKeyArg, StatsFormat,
};

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
    max_edits: Option<usize>,
    format: SearchFormat,
    save_db: &Option<String>,
    subcommand_label: &str,
) -> Result<()> {
    if max_edits.is_some_and(|edits| edits > 0) && pattern.contains('?') {
        bail!("--max-edits cannot currently be combined with ? wildcard patterns");
    }

    let records = seq::read_seqs(file)?;
    let mut all_hits = Vec::new();
    for record in &records {
        let hits = match max_edits {
            Some(edits) => matcher::find_matches_with_max_edits(
                &record.name,
                &record.sequence,
                pattern,
                reverse_complement,
                edits,
            ),
            None => {
                matcher::find_matches(&record.name, &record.sequence, pattern, reverse_complement)
            }
        };
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
        "max_edits": max_edits,
        "format": format!("{format:?}").to_lowercase(),
    });
    maybe_save_hits(save_db, subcommand_label, file, params, &all_hits)
}

fn fm_index_hits(
    reference: &index::ReferenceIndex,
    query: &str,
    include_rc: bool,
) -> Vec<matcher::Hit> {
    let pat_len = query.len();
    let mut hits = Vec::new();
    for start in reference.find_forward(query) {
        let end = start + pat_len;
        hits.push(matcher::Hit {
            record: reference.name.clone(),
            start,
            end,
            strand: "+".to_string(),
            matched: query.to_string(),
            pattern: query.to_string(),
            edits: None,
        });
    }
    if include_rc {
        // `matched` mirrors the non-FM path: report the (uppercased) pattern,
        // not the forward-strand substring. The start index already encodes
        // where the reverse complement sits on the forward strand.
        for start in reference.find_reverse(query) {
            let end = start + pat_len;
            hits.push(matcher::Hit {
                record: reference.name.clone(),
                start,
                end,
                strand: "-".to_string(),
                matched: query.to_string(),
                pattern: query.to_string(),
                edits: None,
            });
        }
    }
    hits
}

#[allow(clippy::too_many_arguments)]
fn run_compare(
    reference: &Path,
    query_files: &[PathBuf],
    reverse_complement: bool,
    max_edits: Option<usize>,
    mode: CompareMode,
    pair_by_order: bool,
    fm_index: bool,
    format: SearchFormat,
) -> Result<()> {
    if query_files.is_empty() {
        bail!("at least one --query FASTA file is required");
    }
    match mode {
        CompareMode::Exact => run_compare_exact(
            reference,
            query_files,
            reverse_complement,
            max_edits,
            pair_by_order,
            fm_index,
            format,
        ),
        CompareMode::Variants => {
            if fm_index {
                bail!("--fm-index is only supported with exact compare mode");
            }
            run_compare_variants(
                reference,
                query_files,
                reverse_complement,
                max_edits,
                pair_by_order,
                format,
            )
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn run_compare_exact(
    reference: &Path,
    query_files: &[PathBuf],
    reverse_complement: bool,
    max_edits: Option<usize>,
    pair_by_order: bool,
    fm_index: bool,
    format: SearchFormat,
) -> Result<()> {
    if pair_by_order {
        bail!("--pair-by-order is only supported with --mode variants");
    }
    if fm_index && max_edits.is_some_and(|edits| edits > 0) {
        bail!("--fm-index is exact-only; remove --max-edits or drop --fm-index");
    }

    let reference_records = seq::read_seqs(reference)?;
    let reference_file = reference.to_string_lossy().into_owned();

    let indexes: Vec<index::ReferenceIndex> = if fm_index {
        reference_records
            .iter()
            .map(|r| index::ReferenceIndex::build(r.name.clone(), &r.sequence))
            .collect()
    } else {
        Vec::new()
    };

    let mut all_hits = Vec::new();

    for query_file in query_files {
        let query_records = seq::read_seqs(query_file)?;
        let query_file_label = query_file.to_string_lossy().into_owned();

        for query_record in &query_records {
            if query_record.sequence.is_empty() {
                bail!(
                    "query record '{}' in {} is empty",
                    query_record.name,
                    query_file.display()
                );
            }

            for (ref_idx, reference_record) in reference_records.iter().enumerate() {
                let hits = if fm_index {
                    fm_index_hits(
                        &indexes[ref_idx],
                        &query_record.sequence,
                        reverse_complement,
                    )
                } else {
                    match max_edits {
                        Some(edits) => matcher::find_matches_with_max_edits(
                            &reference_record.name,
                            &reference_record.sequence,
                            &query_record.sequence,
                            reverse_complement,
                            edits,
                        ),
                        None => matcher::find_matches(
                            &reference_record.name,
                            &reference_record.sequence,
                            &query_record.sequence,
                            reverse_complement,
                        ),
                    }
                };

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
                        edits: hit.edits,
                    });
                }
            }
        }
    }

    print_compare_hits(&all_hits, format)?;
    Ok(())
}

fn run_compare_variants(
    reference: &Path,
    query_files: &[PathBuf],
    reverse_complement: bool,
    max_edits: Option<usize>,
    pair_by_order: bool,
    format: SearchFormat,
) -> Result<()> {
    if reverse_complement {
        bail!("--reverse-complement is only supported with exact compare mode");
    }
    if max_edits.is_some() {
        bail!("--max-edits is only supported with exact compare mode");
    }

    let reference_records = seq::read_seqs(reference)?;
    let reference_file = reference.to_string_lossy().into_owned();
    let reference_by_name: HashMap<&str, &seq::SeqRecord> = reference_records
        .iter()
        .map(|record| (record.name.as_str(), record))
        .collect();
    let mut variants = Vec::new();

    for query_file in query_files {
        let query_records = seq::read_seqs(query_file)?;
        let query_file_label = query_file.to_string_lossy().into_owned();

        if pair_by_order {
            let count = reference_records.len().min(query_records.len());
            for i in 0..count {
                variants.extend(compare_record_variants(
                    &reference_file,
                    &reference_records[i],
                    &query_file_label,
                    &query_records[i],
                ));
            }
            if reference_records.len() > query_records.len() {
                for record in &reference_records[query_records.len()..] {
                    eprintln!(
                        "warning: reference record {} has no paired query record in {}",
                        record.name,
                        query_file.display()
                    );
                }
            } else if query_records.len() > reference_records.len() {
                for record in &query_records[reference_records.len()..] {
                    eprintln!(
                        "warning: query record {} in {} has no paired reference record",
                        record.name,
                        query_file.display()
                    );
                }
            }
            continue;
        }

        let mut seen_references = HashSet::new();
        for query_record in &query_records {
            let Some(reference_record) = reference_by_name.get(query_record.name.as_str()) else {
                eprintln!(
                    "warning: query record {} in {} has no matching reference record",
                    query_record.name,
                    query_file.display()
                );
                continue;
            };
            seen_references.insert(reference_record.name.as_str());
            variants.extend(compare_record_variants(
                &reference_file,
                reference_record,
                &query_file_label,
                query_record,
            ));
        }

        for reference_record in &reference_records {
            if !seen_references.contains(reference_record.name.as_str()) {
                eprintln!(
                    "warning: reference record {} has no matching query record in {}",
                    reference_record.name,
                    query_file.display()
                );
            }
        }
    }

    print_variants(&variants, format)?;
    Ok(())
}

fn compare_record_variants(
    reference_file: &str,
    reference_record: &seq::SeqRecord,
    query_file: &str,
    query_record: &seq::SeqRecord,
) -> Vec<output::Variant> {
    let context = VariantContext {
        reference_file,
        reference_record,
        query_file,
        query_record,
    };

    if reference_record.sequence.len() == query_record.sequence.len() {
        return compare_equal_length_records(context);
    }

    let score = |a: u8, b: u8| if a == b { 1 } else { -1 };
    let mut aligner = Aligner::with_capacity(
        query_record.sequence.len(),
        reference_record.sequence.len(),
        -5,
        -1,
        &score,
    );
    let alignment = aligner.global(
        query_record.sequence.as_bytes(),
        reference_record.sequence.as_bytes(),
    );
    let query = query_record.sequence.as_bytes();
    let reference = reference_record.sequence.as_bytes();
    let mut query_i = alignment.xstart;
    let mut reference_i = alignment.ystart;
    let mut variants = Vec::new();

    for op in alignment.operations {
        match op {
            AlignmentOperation::Match => {
                query_i += 1;
                reference_i += 1;
            }
            AlignmentOperation::Subst => {
                variants.push(make_variant(
                    context,
                    reference_i + 1,
                    byte_to_string(reference[reference_i]),
                    byte_to_string(query[query_i]),
                    "snp",
                ));
                query_i += 1;
                reference_i += 1;
            }
            AlignmentOperation::Ins => {
                variants.push(make_variant(
                    context,
                    reference_i + 1,
                    String::new(),
                    byte_to_string(query[query_i]),
                    "insertion",
                ));
                query_i += 1;
            }
            AlignmentOperation::Del => {
                variants.push(make_variant(
                    context,
                    reference_i + 1,
                    byte_to_string(reference[reference_i]),
                    String::new(),
                    "deletion",
                ));
                reference_i += 1;
            }
            AlignmentOperation::Xclip(len) => {
                query_i += len;
            }
            AlignmentOperation::Yclip(len) => {
                reference_i += len;
            }
        }
    }

    variants
}

#[derive(Clone, Copy)]
struct VariantContext<'a> {
    reference_file: &'a str,
    reference_record: &'a seq::SeqRecord,
    query_file: &'a str,
    query_record: &'a seq::SeqRecord,
}

fn compare_equal_length_records(context: VariantContext<'_>) -> Vec<output::Variant> {
    context
        .reference_record
        .sequence
        .bytes()
        .zip(context.query_record.sequence.bytes())
        .enumerate()
        .filter_map(|(i, (reference_base, query_base))| {
            if reference_base == query_base {
                return None;
            }
            Some(make_variant(
                context,
                i + 1,
                byte_to_string(reference_base),
                byte_to_string(query_base),
                "snp",
            ))
        })
        .collect()
}

fn byte_to_string(base: u8) -> String {
    char::from(base).to_string()
}

fn make_variant(
    context: VariantContext<'_>,
    position: usize,
    reference: String,
    alt: String,
    variant_type: &str,
) -> output::Variant {
    output::Variant {
        reference_file: context.reference_file.to_string(),
        reference_record: context.reference_record.name.clone(),
        query_file: context.query_file.to_string(),
        query_record: context.query_record.name.clone(),
        position,
        reference,
        alt,
        variant_type: variant_type.to_string(),
    }
}

fn print_compare_hits(hits: &[output::CompareHit], format: SearchFormat) -> Result<()> {
    match format {
        SearchFormat::Text => output::print_compare_text(hits),
        SearchFormat::Json => output::print_compare_json(hits)?,
        SearchFormat::Csv => output::print_compare_csv(hits)?,
        SearchFormat::Tsv => output::print_compare_tsv(hits)?,
        SearchFormat::Bed => output::print_compare_bed(hits)?,
    }

    Ok(())
}

fn print_variants(variants: &[output::Variant], format: SearchFormat) -> Result<()> {
    match format {
        SearchFormat::Text => output::print_variants_text(variants),
        SearchFormat::Json => output::print_variants_json(variants)?,
        SearchFormat::Csv => output::print_variants_csv(variants)?,
        SearchFormat::Tsv => output::print_variants_tsv(variants)?,
        SearchFormat::Bed => output::print_variants_bed(variants)?,
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
        Command::Seq { .. } => "seq",
        Command::Convert { .. } => "convert",
        Command::Head { .. } => "head",
        Command::Range { .. } => "range",
        Command::Sample { .. } => "sample",
        Command::Shuffle { .. } => "shuffle",
        Command::Rmdup { .. } => "rmdup",
        Command::Sort { .. } => "sort",
        Command::Split { .. } => "split",
        Command::Locate { .. } => "locate",
        Command::Subseq { .. } => "subseq",
        Command::Sum { .. } => "sum",
        Command::Replace { .. } => "replace",
        Command::Rename { .. } => "rename",
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
            max_edits,
            format,
            save_db,
        } => run_search(
            &file,
            &pattern,
            reverse_complement,
            max_edits,
            format,
            &save_db,
            "search",
        )?,

        Command::Stats {
            file,
            summary: want_summary,
            format,
            save_db,
        } => {
            let records = seq::read_seqs(&file)?;
            let stats: Vec<_> = records.iter().map(output::compute_stats).collect();

            if want_summary {
                let agg = summary::compute(&records);
                print_summary(&agg, format)?;
            } else {
                match format {
                    StatsFormat::Text => output::print_stats_text(&stats),
                    StatsFormat::Json => output::print_stats_json(&stats)?,
                    StatsFormat::Csv => output::print_stats_csv(&stats)?,
                    StatsFormat::Tsv => output::print_stats_tsv(&stats)?,
                }
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
                let records = seq::read_seqs(&path)?;
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
            let records = seq::read_seqs(&file)?;
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
                None,
                format,
                &save_db,
                "motif",
            )?;
        }

        Command::Compare {
            reference,
            query,
            reverse_complement,
            max_edits,
            mode,
            pair_by_order,
            fm_index,
            format,
        } => run_compare(
            &reference,
            &query,
            reverse_complement,
            max_edits,
            mode,
            pair_by_order,
            fm_index,
            format,
        )?,

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

        Command::Seq {
            file,
            out,
            upper,
            lower,
            remove_gaps,
            rna,
            dna,
            reverse,
            complement,
            reverse_complement,
            min_length,
            max_length,
            only_id,
            only_name,
            only_seq,
            wrap,
            format,
        } => {
            let records = seq::read_seqs(&file)?;
            let opts = transform::TransformOptions {
                uppercase: upper,
                lowercase: lower,
                remove_gaps,
                dna_to_rna: rna,
                rna_to_dna: dna,
                reverse,
                complement,
                reverse_complement,
                min_length,
                max_length,
            };
            let mut writer = seq::open_output(out.as_deref())?;
            let any_quality = records.iter().any(|r| r.quality.is_some());
            let wfmt = resolve_record_format(format, any_quality, out.as_deref());
            for r in &records {
                let Some(t) = transform::apply(r, &opts) else {
                    continue;
                };
                if only_id || only_name {
                    writeln!(writer, "{}", t.name)?;
                } else if only_seq {
                    writeln!(writer, "{}", t.sequence)?;
                } else {
                    match wfmt {
                        seq::WriteFormat::Fasta => {
                            seq::write_fasta_record(&mut writer, &t.name, &t.sequence, wrap)?
                        }
                        seq::WriteFormat::Fastq => seq::write_fastq_record(&mut writer, &t)?,
                    }
                }
            }
            writer.flush()?;
        }

        Command::Convert { action } => match action {
            ConvertAction::Fq2fa { file, out, wrap } => {
                let records = seq::read_seqs(&file)?;
                let mut writer = seq::open_output(out.as_deref())?;
                for r in core_convert::fq2fa(&records) {
                    seq::write_fasta_record(&mut writer, &r.name, &r.sequence, wrap)?;
                }
                writer.flush()?;
            }
            ConvertAction::Fx2tab { file, out } => {
                let records = seq::read_seqs(&file)?;
                let mut writer = seq::open_output(out.as_deref())?;
                core_convert::fx2tab(&mut writer, &records)?;
                writer.flush()?;
            }
            ConvertAction::Tab2fx { file, out, format } => {
                let reader = std::io::BufReader::new(std::fs::File::open(&file)?);
                let records = core_convert::tab2fx(reader)?;
                let any_quality = records.iter().any(|r| r.quality.is_some());
                let wfmt = resolve_record_format(format, any_quality, out.as_deref());
                let mut writer = seq::open_output(out.as_deref())?;
                seq::write_records(&mut writer, &records, wfmt, 70)?;
                writer.flush()?;
            }
        },

        Command::Head {
            file,
            number,
            out,
            format,
        } => {
            let records = seq::read_seqs(&file)?;
            let picked = select::head(&records, number);
            emit_records(&picked, out.as_deref(), format)?;
        }

        Command::Range {
            file,
            range,
            out,
            format,
        } => {
            let (a, b) = range
                .split_once(':')
                .ok_or_else(|| anyhow::anyhow!("--range must be START:END, got '{range}'"))?;
            let start: usize = a.parse()?;
            let end: usize = b.parse()?;
            let records = seq::read_seqs(&file)?;
            let picked = select::range(&records, start, end)?;
            emit_records(&picked, out.as_deref(), format)?;
        }

        Command::Sample {
            file,
            number,
            proportion,
            seed,
            out,
            format,
        } => {
            let records = seq::read_seqs(&file)?;
            let kind = match (number, proportion) {
                (Some(n), None) => select::SampleKind::Number(n),
                (None, Some(p)) => {
                    if !(0.0..=1.0).contains(&p) {
                        bail!("--proportion must be between 0.0 and 1.0");
                    }
                    select::SampleKind::Proportion(p)
                }
                (None, None) => bail!("pass --number or --proportion"),
                (Some(_), Some(_)) => unreachable!("clap conflicts_with"),
            };
            let picked = select::sample(&records, kind, seed);
            emit_records(&picked, out.as_deref(), format)?;
        }

        Command::Shuffle {
            file,
            seed,
            out,
            format,
        } => {
            let records = seq::read_seqs(&file)?;
            let picked = select::shuffle(&records, seed);
            emit_records(&picked, out.as_deref(), format)?;
        }

        Command::Rmdup {
            file,
            by,
            out,
            format,
        } => {
            let records = seq::read_seqs(&file)?;
            let key = match by {
                RmdupKeyArg::Id => select::RmdupKey::Id,
                RmdupKeyArg::Sequence => select::RmdupKey::Sequence,
            };
            let (kept, dropped) = select::rmdup(&records, key);
            eprintln!("rmdup: kept {}, dropped {}", kept.len(), dropped);
            emit_records(&kept, out.as_deref(), format)?;
        }

        Command::Sort {
            file,
            by,
            reverse,
            out,
            format,
        } => {
            let records = seq::read_seqs(&file)?;
            let key = match by {
                SortKeyArg::Id => select::SortKey::Id,
                SortKeyArg::Length => select::SortKey::Length,
                SortKeyArg::Sequence => select::SortKey::Sequence,
            };
            let sorted = select::sort(&records, key, reverse);
            emit_records(&sorted, out.as_deref(), format)?;
        }

        Command::Split {
            file,
            parts,
            size,
            out_dir,
            prefix,
            format,
        } => {
            let records = seq::read_seqs(&file)?;
            let chunks = match (parts, size) {
                (Some(n), None) => select::split_into_parts(&records, n),
                (None, Some(s)) => select::split_by_size(&records, s),
                (None, None) => bail!("pass --parts or --size"),
                (Some(_), Some(_)) => unreachable!("clap conflicts_with"),
            };
            std::fs::create_dir_all(&out_dir)?;
            let stem = prefix.unwrap_or_else(|| {
                file.file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("part")
                    .to_string()
            });
            let any_quality = records.iter().any(|r| r.quality.is_some());
            let wfmt = resolve_record_format(format, any_quality, None);
            let ext = match wfmt {
                seq::WriteFormat::Fasta => "fasta",
                seq::WriteFormat::Fastq => "fastq",
            };
            for (i, chunk) in chunks.iter().enumerate() {
                let path = out_dir.join(format!("{stem}.part_{:03}.{ext}", i + 1));
                let mut w = seq::open_output(Some(&path))?;
                seq::write_records(&mut w, chunk, wfmt, 70)?;
                w.flush()?;
                eprintln!("wrote {} records → {}", chunk.len(), path.display());
            }
        }

        Command::Locate {
            file,
            pattern,
            pattern_file,
            regex,
            reverse_complement,
            format,
        } => {
            let mut raw_patterns: Vec<(String, String)> = pattern
                .iter()
                .enumerate()
                .map(|(i, p)| (format!("p{}", i + 1), p.clone()))
                .collect();
            if let Some(pf) = pattern_file {
                let text = std::fs::read_to_string(&pf)?;
                for (i, line) in text.lines().enumerate() {
                    let line = line.trim();
                    if line.is_empty() || line.starts_with('#') {
                        continue;
                    }
                    raw_patterns.push((format!("file_p{}", i + 1), line.to_string()));
                }
            }
            if raw_patterns.is_empty() {
                bail!("provide --pattern <P> and/or --pattern-file <FILE>");
            }
            let patterns: Vec<locate::LocatePattern> = raw_patterns
                .into_iter()
                .map(|(name, raw)| {
                    if regex {
                        locate::LocatePattern::new(name, raw)
                    } else {
                        locate::LocatePattern::from_iupac(name, raw)
                    }
                })
                .collect::<Result<_>>()?;
            let records = seq::read_seqs(&file)?;
            let mut hits = Vec::new();
            for r in &records {
                hits.extend(locate::locate_in_record(r, &patterns, reverse_complement));
            }
            print_locate_hits(&hits, format)?;
        }

        Command::Subseq {
            file,
            region,
            up_stream,
            down_stream,
            reverse_complement,
            out,
            format,
        } => {
            let region = core_subseq::Region::parse(&region)?;
            let opts = core_subseq::SubseqOptions {
                region,
                reverse_complement,
                up_stream,
                down_stream,
            };
            let records = seq::read_seqs(&file)?;
            let mut extracted = Vec::with_capacity(records.len());
            for r in &records {
                if let Some(out_rec) = core_subseq::subseq_record(r, &opts)? {
                    extracted.push(out_rec);
                }
            }
            emit_records(&extracted, out.as_deref(), format)?;
        }

        Command::Sum { file, format } => {
            let records = seq::read_seqs(&file)?;
            let digests = digest::digest_records(&records);
            print_digests(&digests, format)?;
        }

        Command::Replace {
            file,
            pattern,
            replacement,
            by_seq,
            out,
            format,
        } => {
            let records = seq::read_seqs(&file)?;
            let edited = if by_seq {
                core_edit::replace_seq(&records, &pattern, &replacement)?
            } else {
                core_edit::replace_id(&records, &pattern, &replacement)?
            };
            emit_records(&edited, out.as_deref(), format)?;
        }

        Command::Rename { file, out, format } => {
            let records = seq::read_seqs(&file)?;
            let (renamed, n) = core_edit::rename_duplicates(&records);
            eprintln!("rename: {n} record(s) renamed");
            emit_records(&renamed, out.as_deref(), format)?;
        }
    }

    Ok(())
}

/// Pick FASTA vs FASTQ output based on:
/// 1. an explicit `--format` flag (Fasta/Fastq), 2. the destination extension
///    if `--out` was given, 3. whether the input carries quality data.
fn resolve_record_format(
    flag: RecordFormat,
    input_has_quality: bool,
    out: Option<&Path>,
) -> seq::WriteFormat {
    match flag {
        RecordFormat::Fasta | RecordFormat::Tab => seq::WriteFormat::Fasta,
        RecordFormat::Fastq => seq::WriteFormat::Fastq,
        RecordFormat::Auto => {
            if let Some(p) = out
                && let Some(ext) = p.extension().and_then(|s| s.to_str())
            {
                let ext = ext.to_ascii_lowercase();
                if matches!(ext.as_str(), "fq" | "fastq") {
                    return seq::WriteFormat::Fastq;
                }
                if matches!(ext.as_str(), "fa" | "fna" | "fasta") {
                    return seq::WriteFormat::Fasta;
                }
            }
            if input_has_quality {
                seq::WriteFormat::Fastq
            } else {
                seq::WriteFormat::Fasta
            }
        }
    }
}

fn emit_records(
    records: &[seq::SeqRecord],
    out: Option<&Path>,
    format: RecordFormat,
) -> Result<()> {
    let any_quality = records.iter().any(|r| r.quality.is_some());
    let wfmt = resolve_record_format(format, any_quality, out);
    let mut w = seq::open_output(out)?;
    if matches!(format, RecordFormat::Tab) {
        core_convert::fx2tab(&mut w, records)?;
    } else {
        seq::write_records(&mut w, records, wfmt, 70)?;
    }
    w.flush()?;
    Ok(())
}

fn print_summary(s: &summary::StatsSummary, format: StatsFormat) -> Result<()> {
    match format {
        StatsFormat::Text => {
            println!("num_seqs:    {}", s.num_seqs);
            println!("sum_len:     {}", s.sum_len);
            println!("min_len:     {}", s.min_len);
            println!("avg_len:     {:.2}", s.avg_len);
            println!("max_len:     {}", s.max_len);
            println!("N50:         {}", s.n50);
            println!("GC%:         {:.2}", s.gc_percent);
            println!("N_count:     {}", s.n_count);
            if let Some(q20) = s.q20_percent {
                println!("Q20%:        {:.2}", q20);
            }
            if let Some(q30) = s.q30_percent {
                println!("Q30%:        {:.2}", q30);
            }
        }
        StatsFormat::Json => {
            println!("{}", serde_json::to_string_pretty(s)?);
        }
        StatsFormat::Csv => {
            let mut w = csv::Writer::from_writer(std::io::stdout());
            w.serialize(s)?;
            w.flush()?;
        }
        StatsFormat::Tsv => {
            let mut w = csv::WriterBuilder::new()
                .delimiter(b'\t')
                .from_writer(std::io::stdout());
            w.serialize(s)?;
            w.flush()?;
        }
    }
    Ok(())
}

fn print_locate_hits(hits: &[locate::LocateHit], format: LocateFormat) -> Result<()> {
    match format {
        LocateFormat::Text => {
            for h in hits {
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    h.record, h.start, h.end, h.strand, h.pattern_name, h.pattern, h.matched
                );
            }
        }
        LocateFormat::Json => println!("{}", serde_json::to_string_pretty(hits)?),
        LocateFormat::Csv => {
            let mut w = csv::Writer::from_writer(std::io::stdout());
            for h in hits {
                w.serialize(h)?;
            }
            w.flush()?;
        }
        LocateFormat::Tsv => {
            let mut w = csv::WriterBuilder::new()
                .delimiter(b'\t')
                .from_writer(std::io::stdout());
            for h in hits {
                w.serialize(h)?;
            }
            w.flush()?;
        }
        LocateFormat::Bed => {
            for h in hits {
                println!(
                    "{}\t{}\t{}\t{}\t.\t{}",
                    h.record, h.start, h.end, h.pattern_name, h.strand
                );
            }
        }
    }
    Ok(())
}

fn print_digests(rows: &[digest::RecordDigest], format: DigestFormat) -> Result<()> {
    match format {
        DigestFormat::Text => {
            for r in rows {
                println!("{}\t{}\t{}", r.sha256, r.length, r.record);
            }
        }
        DigestFormat::Json => println!("{}", serde_json::to_string_pretty(rows)?),
        DigestFormat::Csv => {
            let mut w = csv::Writer::from_writer(std::io::stdout());
            for r in rows {
                w.serialize(r)?;
            }
            w.flush()?;
        }
        DigestFormat::Tsv => {
            let mut w = csv::WriterBuilder::new()
                .delimiter(b'\t')
                .from_writer(std::io::stdout());
            for r in rows {
                w.serialize(r)?;
            }
            w.flush()?;
        }
    }
    Ok(())
}
