//! Library facade -exposes the deterministic core modules so the CLI binary,
//! the HTTP API, integration tests, and benchmarks can all call them directly.

pub mod alphabets;
pub mod convert;
pub mod db;
pub mod digest;
pub mod edit;
pub mod index;
pub mod locate;
pub mod matcher;
pub mod motif;
pub mod ncbi;
pub mod orf;
pub mod output;
pub mod select;
pub mod seq;
pub mod subseq;
pub mod summary;
pub mod transform;
