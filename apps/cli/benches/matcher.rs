//! Benchmark suite for the matching pipelines.
//!
//! Covers exact (BNDM/Horspool dispatch), wildcard, approximate (Myers), and
//! repeated-search (FM-index) workloads on synthetic SARS-CoV-2-sized,
//! E. coli-sized, and human-chromosome-sized inputs. Real FASTA files under
//! `data/` are not required -sequences are seeded deterministically so the
//! benchmark is reproducible from a clean clone.
//!
//! Run with: `cargo bench --bench matcher`

use criterion::{Criterion, criterion_group, criterion_main};
use seqfind_core::index::ReferenceIndex;
use seqfind_core::matcher::{
    exact_match_positions, find_matches, find_matches_with_max_edits, reverse_complement,
};

/// Deterministic pseudo-random DNA sequence of `n` bases. xorshift keeps the
/// generator dependency-free.
fn synthetic_dna(n: usize, seed: u64) -> String {
    let mut state = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let bases = [b'A', b'C', b'G', b'T'];
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        out.push(bases[(state & 0b11) as usize]);
    }
    String::from_utf8(out).unwrap()
}

fn bench_exact(c: &mut Criterion) {
    let mut group = c.benchmark_group("exact_search");

    // SARS-CoV-2 sized (~30 kb)
    let sars = synthetic_dna(30_000, 1);
    group.bench_function("sars_cov2_30kb_pattern20", |b| {
        let pat = &sars[12_345..12_365];
        b.iter(|| exact_match_positions(sars.as_bytes(), pat.as_bytes()).len());
    });

    // E. coli sized (~4.6 Mb)
    let ecoli = synthetic_dna(4_641_652, 2);
    group.bench_function("ecoli_4mb_pattern20", |b| {
        let pat = &ecoli[1_234_567..1_234_587];
        b.iter(|| exact_match_positions(ecoli.as_bytes(), pat.as_bytes()).len());
    });

    // Long pattern (>64) falls back to Horspool.
    group.bench_function("ecoli_4mb_pattern128_horspool", |b| {
        let pat = &ecoli[1_234_567..1_234_695];
        b.iter(|| exact_match_positions(ecoli.as_bytes(), pat.as_bytes()).len());
    });

    // Human-chromosome-sized slice (~10 Mb keeps the benchmark wall time bounded
    // while still exercising the chromosome-scale code path).
    let chr_slice = synthetic_dna(10_000_000, 3);
    group.bench_function("chr_slice_10mb_pattern30", |b| {
        let pat = &chr_slice[5_000_000..5_000_030];
        b.iter(|| exact_match_positions(chr_slice.as_bytes(), pat.as_bytes()).len());
    });

    group.finish();
}

fn bench_wildcard(c: &mut Criterion) {
    let mut group = c.benchmark_group("wildcard_search");
    let sars = synthetic_dna(30_000, 4);
    group.bench_function("sars_cov2_wildcard", |b| {
        b.iter(|| find_matches("ref", &sars, "ATG???TAA", false).len());
    });
    group.finish();
}

fn bench_approximate(c: &mut Criterion) {
    let mut group = c.benchmark_group("approximate_search");
    let sars = synthetic_dna(30_000, 5);
    let pat = "ATGCGTAACCCATGCGT";
    group.bench_function("sars_cov2_myers_edits2", |b| {
        b.iter(|| find_matches_with_max_edits("ref", &sars, pat, false, 2).len());
    });
    group.finish();
}

fn bench_compare_repeated(c: &mut Criterion) {
    let mut group = c.benchmark_group("compare_repeated");

    // Simulate a chromosome-scale reference with 50 short queries -the case
    // FM-index is meant to win.
    let reference = synthetic_dna(10_000_000, 6);
    let queries: Vec<String> = (0..50)
        .map(|i| reference[(100_000 + i * 50_000)..(100_000 + i * 50_000 + 30)].to_string())
        .collect();

    group.bench_function("horspool_bndm_dispatch", |b| {
        b.iter(|| {
            let mut total = 0usize;
            for q in &queries {
                total += find_matches("ref", &reference, q, false).len();
            }
            total
        });
    });

    group.bench_function("fm_index_build_and_query", |b| {
        b.iter(|| {
            let idx = ReferenceIndex::build("ref".into(), &reference);
            let mut total = 0usize;
            for q in &queries {
                total += idx.find_forward(q).len();
            }
            total
        });
    });

    group.bench_function("fm_index_query_only_prebuilt", |b| {
        let idx = ReferenceIndex::build("ref".into(), &reference);
        b.iter(|| {
            let mut total = 0usize;
            for q in &queries {
                total += idx.find_forward(q).len();
            }
            total
        });
    });

    group.finish();
}

fn bench_reverse_complement(c: &mut Criterion) {
    let mut group = c.benchmark_group("reverse_complement");
    let sars = synthetic_dna(30_000, 7);
    group.bench_function("rc_30kb", |b| {
        b.iter(|| reverse_complement(&sars).len());
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_exact,
    bench_wildcard,
    bench_approximate,
    bench_compare_repeated,
    bench_reverse_complement,
);
criterion_main!(benches);
