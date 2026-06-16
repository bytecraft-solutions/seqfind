use assert_cmd::Command;
use predicates::prelude::*;
use std::io::Write;
use tempfile::NamedTempFile;

// ---------------------------------------------------------------------------
// Postgres helpers -DB tests use a real, shared Postgres instance and silently
// skip when one isn't reachable.

const DEFAULT_PG_URL: &str = "postgresql://postgres:passw0rd@127.0.0.1/seqfind";

fn pg_url() -> String {
    std::env::var("SEQFIND_TEST_PG_URL").unwrap_or_else(|_| DEFAULT_PG_URL.to_string())
}

/// Returns Some(url) if a server is reachable on the test URL, None otherwise.
fn pg_available() -> Option<String> {
    let url = pg_url();
    let ok = Command::cargo_bin("seqfind")
        .unwrap()
        .env("SEQFIND_DB_URL", &url)
        .args(["runs", "list", "--limit", "1"])
        .ok()
        .is_ok();
    if ok { Some(url) } else { None }
}

/// Pull `saved run #N` out of the binary's stderr.
fn parse_run_id(stderr: &[u8]) -> i64 {
    let s = std::str::from_utf8(stderr).unwrap();
    let after = s.split("saved run #").nth(1).expect("no run id in stderr");
    let digits: String = after.chars().take_while(|c| c.is_ascii_digit()).collect();
    digits.parse().expect("non-numeric run id")
}

fn fixture() -> NamedTempFile {
    // chr1 (33 bp): ATGCGT at positions 0 and 11.
    // chr2 (29 bp): ATGCGT at position 20; ACGCAT at position 2 (RC of ATGCGT).
    let mut f = NamedTempFile::new().unwrap();
    writeln!(
        f,
        ">chr1\nATGCGTAACCC\nATGCGTAAAAAAGGGGGCCCCC\n>chr2\nTTACGCATGGGCCCAAATTTATGCGTAAA"
    )
    .unwrap();
    f
}

fn bin() -> Command {
    Command::cargo_bin("seqfind").unwrap()
}

#[test]
fn search_text_format() {
    let f = fixture();
    bin()
        .args(["search", "--file"])
        .arg(f.path())
        .args(["--pattern", "ATGCGT"])
        .assert()
        .success()
        .stdout(predicate::str::contains(
            "chr1\tstart=0\tend=6\tstrand=+\tmatched=ATGCGT\tpattern=ATGCGT",
        ))
        .stdout(predicate::str::contains(
            "chr1\tstart=11\tend=17\tstrand=+\tmatched=ATGCGT\tpattern=ATGCGT",
        ))
        .stdout(predicate::str::contains(
            "chr2\tstart=20\tend=26\tstrand=+\tmatched=ATGCGT\tpattern=ATGCGT",
        ));
}

#[test]
fn search_reports_timing_to_stderr() {
    let f = fixture();
    bin()
        .args(["search", "--file"])
        .arg(f.path())
        .args(["--pattern", "ATGCGT"])
        .assert()
        .success()
        .stderr(
            predicate::str::contains("search started:")
                .and(predicate::str::contains("search ended:"))
                .and(predicate::str::contains("search took:"))
                .and(predicate::str::contains("(ok)")),
        );
}

#[test]
fn search_json_format() {
    let f = fixture();
    let out = bin()
        .args(["search", "--file"])
        .arg(f.path())
        .args(["--pattern", "ATGCGT", "--format", "json"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let v: serde_json::Value = serde_json::from_slice(&out).expect("valid JSON");
    let arr = v.as_array().unwrap();
    assert_eq!(arr.len(), 3);
    assert_eq!(arr[0]["record"], "chr1");
    assert_eq!(arr[0]["start"], 0);
    assert_eq!(arr[0]["end"], 6);
    assert_eq!(arr[0]["strand"], "+");
    assert_eq!(arr[0]["matched"], "ATGCGT");
    assert_eq!(arr[0]["pattern"], "ATGCGT");
}

#[test]
fn search_csv_format_has_header_and_rows() {
    let f = fixture();
    bin()
        .args(["search", "--file"])
        .arg(f.path())
        .args(["--pattern", "ATGCGT", "--format", "csv"])
        .assert()
        .success()
        .stdout(predicate::str::starts_with(
            "record,start,end,strand,matched,pattern\n",
        ))
        .stdout(predicate::str::contains("chr1,0,6,+,ATGCGT,ATGCGT"))
        .stdout(predicate::str::contains("chr2,20,26,+,ATGCGT,ATGCGT"));
}

#[test]
fn search_tsv_format_has_header_and_rows() {
    let f = fixture();
    bin()
        .args(["search", "--file"])
        .arg(f.path())
        .args(["--pattern", "ATGCGT", "--format", "tsv"])
        .assert()
        .success()
        .stdout(predicate::str::starts_with(
            "record\tstart\tend\tstrand\tmatched\tpattern\n",
        ))
        .stdout(predicate::str::contains("chr1\t0\t6\t+\tATGCGT\tATGCGT"));
}

#[test]
fn search_bed_format_no_header_with_intervals() {
    let f = fixture();
    bin()
        .args(["search", "--file"])
        .arg(f.path())
        .args(["--pattern", "ATGCGT", "--format", "bed"])
        .assert()
        .success()
        // BED3+6 has no header. chr1 pos 0 → end 6, score 0, strand +.
        .stdout(predicate::str::starts_with("chr1\t0\t6\tATGCGT\t0\t+\n"))
        .stdout(predicate::str::contains("chr1\t11\t17\tATGCGT\t0\t+"))
        .stdout(predicate::str::contains("chr2\t20\t26\tATGCGT\t0\t+"));
}

#[test]
fn search_rc_strand_uses_forward_coordinates_in_bed() {
    // chr2 has ACGCAT at index 2; ACGCAT is the RC of ATGCGT.
    // RC hit reports the forward-strand index of the RC match -i.e., 2.
    let f = fixture();
    bin()
        .args(["search", "--file"])
        .arg(f.path())
        .args([
            "--pattern",
            "ATGCGT",
            "--reverse-complement",
            "--format",
            "bed",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("chr2\t2\t8\tATGCGT\t0\t-"));
}

#[test]
fn stats_text_format() {
    let f = fixture();
    bin()
        .args(["stats", "--file"])
        .arg(f.path())
        .assert()
        .success()
        .stdout(predicate::str::contains("chr1\tlength=33"))
        .stdout(predicate::str::contains("chr2\tlength=29"));
}

#[test]
fn stats_accepts_time_flag_after_subcommand() {
    let f = fixture();
    bin()
        .args(["stats", "--file"])
        .arg(f.path())
        .arg("--time")
        .assert()
        .success()
        .stdout(predicate::str::contains("chr1\tlength=33"))
        .stderr(predicate::str::contains("stats took:"));
}

#[test]
fn stats_json_format() {
    let f = fixture();
    let out = bin()
        .args(["stats", "--file"])
        .arg(f.path())
        .args(["--format", "json"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let v: serde_json::Value = serde_json::from_slice(&out).expect("valid JSON");
    let arr = v.as_array().unwrap();
    assert_eq!(arr.len(), 2);
    assert_eq!(arr[0]["record"], "chr1");
    assert_eq!(arr[0]["length"], 33);
}

#[test]
fn stats_csv_format() {
    let f = fixture();
    bin()
        .args(["stats", "--file"])
        .arg(f.path())
        .args(["--format", "csv"])
        .assert()
        .success()
        .stdout(predicate::str::starts_with(
            "record,length,gc_percent,a,t,g,c,n\n",
        ))
        .stdout(predicate::str::contains("chr1,33,"));
}

#[test]
fn stats_tsv_format() {
    let f = fixture();
    bin()
        .args(["stats", "--file"])
        .arg(f.path())
        .args(["--format", "tsv"])
        .assert()
        .success()
        .stdout(predicate::str::starts_with(
            "record\tlength\tgc_percent\ta\tt\tg\tc\tn\n",
        ))
        .stdout(predicate::str::contains("chr1\t33\t"));
}

#[test]
fn rc_subcommand() {
    bin()
        .args(["rc", "ATGCGT"])
        .assert()
        .success()
        .stdout("ACGCAT\n");
}

#[test]
fn rc_file_emits_fasta_per_record() {
    let mut f = NamedTempFile::new().unwrap();
    writeln!(f, ">a\nATGC\n>b\nGGGG").unwrap();
    bin()
        .args(["rc", "--file"])
        .arg(f.path())
        .assert()
        .success()
        // RC of ATGC = GCAT, RC of GGGG = CCCC
        .stdout(predicate::str::contains(">a_rc\nGCAT"))
        .stdout(predicate::str::contains(">b_rc\nCCCC"));
}

#[test]
fn rc_rejects_both_inline_and_file() {
    let mut f = NamedTempFile::new().unwrap();
    writeln!(f, ">a\nATGC").unwrap();
    bin()
        .args(["rc", "ATGC", "--file"])
        .arg(f.path())
        .assert()
        .failure();
}

fn orf_fixture() -> NamedTempFile {
    // chr1 = ATG CCC TAA (frame 0, +) = 9 bp ORF starting at 0
    // Tail keeps it unambiguous.
    let mut f = NamedTempFile::new().unwrap();
    writeln!(f, ">chr1\nATGCCCTAAGGG\n>chr2\nTTAGGGCAT").unwrap();
    f
}

#[test]
fn orf_text_finds_forward_orf() {
    let f = orf_fixture();
    bin()
        .args(["orf", "--file"])
        .arg(f.path())
        .args(["--min-length", "9", "--forward-only"])
        .assert()
        .success()
        .stdout(predicate::str::contains(
            "chr1\tstart=0\tend=9\tlength=9\tstrand=+\tstart_codon=ATG\tstop_codon=TAA",
        ));
}

#[test]
fn orf_json_schema_matches_mvp_doc() {
    let f = orf_fixture();
    let out = bin()
        .args(["orf", "--file"])
        .arg(f.path())
        .args(["--min-length", "9", "--forward-only", "--format", "json"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let v: serde_json::Value = serde_json::from_slice(&out).unwrap();
    let arr = v.as_array().unwrap();
    let o = &arr[0];
    // Schema fields per docs/MVP_PHASES.md
    for key in [
        "record",
        "start",
        "end",
        "strand",
        "start_codon",
        "stop_codon",
        "length",
    ] {
        assert!(o.get(key).is_some(), "missing key {key}");
    }
    assert_eq!(o["record"], "chr1");
    assert_eq!(o["start"], 0);
    assert_eq!(o["end"], 9);
    assert_eq!(o["strand"], "+");
    assert_eq!(o["start_codon"], "ATG");
    assert_eq!(o["stop_codon"], "TAA");
    assert_eq!(o["length"], 9);
}

#[test]
fn orf_bed_format() {
    let f = orf_fixture();
    bin()
        .args(["orf", "--file"])
        .arg(f.path())
        .args(["--min-length", "9", "--forward-only", "--format", "bed"])
        .assert()
        .success()
        .stdout(predicate::str::starts_with(
            "chr1\t0\t9\tORF_ATG_TAA\t9\t+\n",
        ));
}

#[test]
fn orf_finds_reverse_strand_orf() {
    let f = orf_fixture();
    // chr2 = TTAGGGCAT; RC = ATGCCCTAA → ORF on minus strand covering 0..9
    let out = bin()
        .args(["orf", "--file"])
        .arg(f.path())
        .args(["--min-length", "9", "--format", "json"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let v: serde_json::Value = serde_json::from_slice(&out).unwrap();
    let arr = v.as_array().unwrap();
    let minus: Vec<&serde_json::Value> = arr
        .iter()
        .filter(|o| o["strand"] == "-" && o["record"] == "chr2")
        .collect();
    assert_eq!(minus.len(), 1);
    assert_eq!(minus[0]["start"], 0);
    assert_eq!(minus[0]["end"], 9);
}

#[test]
fn motif_list_prints_presets() {
    bin()
        .args(["motif", "--list"])
        .assert()
        .success()
        .stdout(predicate::str::contains("EcoRI"))
        .stdout(predicate::str::contains("GAATTC"))
        .stdout(predicate::str::contains("TATA"));
}

#[test]
fn motif_search_resolves_preset() {
    let mut f = NamedTempFile::new().unwrap();
    writeln!(f, ">chr1\nAAGAATTCGG").unwrap();
    bin()
        .args(["motif", "--file"])
        .arg(f.path())
        .args(["--name", "ecori", "--format", "tsv"])
        .assert()
        .success()
        .stdout(predicate::str::contains("chr1\t2\t8\t+\tGAATTC\tGAATTC"));
}

#[test]
fn motif_unknown_name_errors() {
    let mut f = NamedTempFile::new().unwrap();
    writeln!(f, ">a\nA").unwrap();
    bin()
        .args(["motif", "--file"])
        .arg(f.path())
        .args(["--name", "NotAReal"])
        .assert()
        .failure();
}

#[test]
fn search_save_db_via_env() {
    let Some(url) = pg_available() else {
        eprintln!("skipping: no Postgres at {}", pg_url());
        return;
    };
    let f = fixture();

    // No URL on the CLI -`--save-db` alone reads from $SEQFIND_DB_URL.
    let out = bin()
        .env("SEQFIND_DB_URL", &url)
        .args(["search", "--file"])
        .arg(f.path())
        .args(["--pattern", "ATGCGT", "--save-db"])
        .assert()
        .success()
        .stderr(predicate::str::contains("saved run #"))
        .get_output()
        .clone();
    let run_id = parse_run_id(&out.stderr);

    // Export reads $SEQFIND_DB_URL via clap's `env` attribute on --db.
    let out = bin()
        .env("SEQFIND_DB_URL", &url)
        .args(["runs", "export", &run_id.to_string(), "--format", "json"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let v: serde_json::Value = serde_json::from_slice(&out).unwrap();
    let arr = v.as_array().unwrap();
    assert_eq!(arr.len(), 3);
    assert_eq!(arr[0]["pattern"], "ATGCGT");
}

#[test]
fn save_db_url_flag_overrides_env() {
    let Some(url) = pg_available() else {
        eprintln!("skipping: no Postgres at {}", pg_url());
        return;
    };
    let f = orf_fixture();

    // Env points to a bogus URL; --save-db <URL> overrides and succeeds.
    let out = bin()
        .env(
            "SEQFIND_DB_URL",
            "postgresql://postgres:wrong@127.0.0.1:1/seqfind",
        )
        .args(["orf", "--file"])
        .arg(f.path())
        .args(["--min-length", "9", "--forward-only", "--save-db", &url])
        .assert()
        .success()
        .stderr(predicate::str::contains("ORFs"))
        .get_output()
        .clone();
    let run_id = parse_run_id(&out.stderr);

    bin()
        .env("SEQFIND_DB_URL", &url)
        .args(["runs", "show", &run_id.to_string()])
        .assert()
        .success()
        .stdout(predicate::str::contains("kind:      orf"))
        .stdout(predicate::str::contains("orfs="));
}

#[test]
fn stats_save_db_persists_records() {
    let Some(url) = pg_available() else {
        eprintln!("skipping: no Postgres at {}", pg_url());
        return;
    };
    let f = fixture();

    let out = bin()
        .env("SEQFIND_DB_URL", &url)
        .args(["stats", "--file"])
        .arg(f.path())
        .args(["--save-db"])
        .assert()
        .success()
        .get_output()
        .clone();
    let run_id = parse_run_id(&out.stderr);

    let out = bin()
        .env("SEQFIND_DB_URL", &url)
        .args(["runs", "export", &run_id.to_string(), "--format", "json"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let v: serde_json::Value = serde_json::from_slice(&out).unwrap();
    let arr = v.as_array().unwrap();
    assert_eq!(arr.len(), 2);
    assert_eq!(arr[0]["record"], "chr1");
    assert_eq!(arr[0]["length"], 33);
}

#[test]
fn runs_list_requires_env_or_flag() {
    // No SEQFIND_DB_URL and no --db → clap errors before we touch Postgres.
    bin()
        .env_remove("SEQFIND_DB_URL")
        .args(["runs", "list"])
        .assert()
        .failure();
}

#[test]
fn save_db_without_env_or_url_errors() {
    let f = fixture();
    bin()
        .env_remove("SEQFIND_DB_URL")
        .args(["search", "--file"])
        .arg(f.path())
        .args(["--pattern", "ATGCGT", "--save-db"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("SEQFIND_DB_URL"));
}

#[test]
fn fetch_rejects_garbage_accession() {
    // Definitely-invalid accession; NCBI returns a plain-text error.
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("garbage.fna");
    bin()
        .args([
            "fetch",
            "--accession",
            "DEFINITELY_NOT_AN_ACCESSION_XYZ",
            "--out",
        ])
        .arg(&out)
        .assert()
        .failure();
}

#[test]
fn missing_file_errors() {
    bin()
        .args([
            "search",
            "--file",
            "/nonexistent/path.fasta",
            "--pattern",
            "A",
        ])
        .assert()
        .failure();
}
