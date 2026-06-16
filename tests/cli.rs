use assert_cmd::Command;
use predicates::prelude::*;
use std::io::Write;
use tempfile::NamedTempFile;

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
            "chr1\tpos=0\tstrand=+\tmatched=ATGCGT",
        ))
        .stdout(predicate::str::contains(
            "chr1\tpos=11\tstrand=+\tmatched=ATGCGT",
        ))
        .stdout(predicate::str::contains(
            "chr2\tpos=20\tstrand=+\tmatched=ATGCGT",
        ));
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
    assert_eq!(arr[0]["position"], 0);
    assert_eq!(arr[0]["strand"], "+");
    assert_eq!(arr[0]["matched"], "ATGCGT");
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
            "record,position,strand,matched\n",
        ))
        .stdout(predicate::str::contains("chr1,0,+,ATGCGT"))
        .stdout(predicate::str::contains("chr2,20,+,ATGCGT"));
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
            "record\tposition\tstrand\tmatched\n",
        ))
        .stdout(predicate::str::contains("chr1\t0\t+\tATGCGT"));
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
