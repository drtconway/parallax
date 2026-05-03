//! Parallax - Sequence indexing and alignment utilities
//!
//! This library provides alignment algorithms and utilities for DNA sequence analysis.

use crate::{index::Index, reference::InMemoryReference, writer::AlignmentWriter};

pub mod align;
pub mod annotate;
pub mod cluster;
pub mod config;
pub mod error;
pub mod index;
pub mod kmers;
pub mod metrics;
pub mod reads;
pub mod reference;
pub mod scores;
pub mod seeding;
pub mod utils;
pub mod validation;
pub mod writer;

pub trait AlignerBuilder<'a, const K: usize, const S: usize> {
    type AlignerType: Aligner<'a, K, S>;

    fn new(reference: &'a InMemoryReference, index: &'a Index<K, S>, writer: &'a AlignmentWriter) -> Self;

    fn build(self) -> Self::AlignerType;
}

pub trait Aligner<'a, const K: usize, const S: usize> {
    fn align(&mut self, name: &str, query: &[u8], quality: &[u8]) -> std::io::Result<()>;

    fn finish(self) -> std::io::Result<()>;
}

pub mod explanatory;
