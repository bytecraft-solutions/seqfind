//! Library facade -exposes the deterministic core modules so the CLI binary,
//! the HTTP API, integration tests, and benchmarks can all call them directly.

pub mod alphabets;
pub mod db;
pub mod fasta;
pub mod index;
pub mod matcher;
pub mod motif;
pub mod ncbi;
pub mod orf;
pub mod output;
