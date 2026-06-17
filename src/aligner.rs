use parallax::{index::Index, reference::InMemoryReference};
use crate::writer::RecordWriter;

pub trait AlignerBuilder<'a> {
    type AlignerType: Aligner<'a>;

    fn new(reference: &'a InMemoryReference, index: &'a dyn Index, writer: &'a dyn RecordWriter) -> Self;

    fn build(self) -> Self::AlignerType;
}

pub trait Aligner<'a> {
    fn align(&mut self, name: &str, query: &[u8], quality: &[u8]) -> std::io::Result<()>;

    fn finish(self) -> std::io::Result<()>;
}
