use parallax::{index::Index, reference::InMemoryReference};
use crate::writer::AlignmentWriter;

pub trait AlignerBuilder<'a, const K: usize, const S: usize> {
    type AlignerType: Aligner<'a, K, S>;

    fn new(reference: &'a InMemoryReference, index: &'a Index<K, S>, writer: &'a AlignmentWriter) -> Self;

    fn build(self) -> Self::AlignerType;
}

pub trait Aligner<'a, const K: usize, const S: usize> {
    fn align(&mut self, name: &str, query: &[u8], quality: &[u8]) -> std::io::Result<()>;

    fn finish(self) -> std::io::Result<()>;
}
