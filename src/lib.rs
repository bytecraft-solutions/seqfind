//! Library facade -exposes the deterministic core modules so integration
//! tests and benchmarks can call them directly without going through the CLI.

pub mod alphabets;
pub mod fasta;
pub mod index;
pub mod matcher;
pub mod motif;
pub mod orf;
