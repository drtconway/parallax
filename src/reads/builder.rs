use noodles::sam::alignment::{
    record::Flags,
    record::MappingQuality,
    record_buf::{Cigar, Data, QualityScores, RecordBuf, Sequence},
};

/// Build a RecordBuf for a mapped alignment segment.
///
/// Quality scores are expected as Phred+33 ASCII (FASTQ convention);
/// they are converted to raw Phred (0-based) for noodles.
#[allow(clippy::too_many_arguments)]
pub fn build_record(
    name: &str,
    flags: Flags,
    reference_sequence_id: usize,
    pos: usize, // 1-based
    mapq: u8,
    cigar: Cigar,
    mate_ref_id: Option<usize>,
    mate_pos: Option<usize>, // 1-based
    seq: &[u8],              // ASCII nucleotides
    qual: &[u8],             // Phred+33 ASCII
    data: Data,
) -> RecordBuf {
    let mut builder = RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_reference_sequence_id(reference_sequence_id)
        .set_alignment_start(
            noodles::core::Position::try_from(pos).expect("alignment position must be >= 1"),
        )
        .set_mapping_quality(
            MappingQuality::try_from(mapq.min(254)).expect("mapping quality must be <= 254"),
        )
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(
            qual.iter().map(|&q| q.saturating_sub(33)).collect::<Vec<u8>>(),
        ))
        .set_data(data);

    if let Some(mate_id) = mate_ref_id {
        builder = builder.set_mate_reference_sequence_id(mate_id);
    }
    if let Some(mp) = mate_pos {
        builder = builder.set_mate_alignment_start(
            noodles::core::Position::try_from(mp).expect("mate position must be >= 1"),
        );
    }

    builder.build()
}

/// Build a RecordBuf for an unmapped read.
///
/// Quality scores are expected as Phred+33 ASCII.
pub fn build_unmapped_record(name: &str, seq: &[u8], qual: &[u8]) -> RecordBuf {
    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::UNMAPPED)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(
            qual.iter().map(|&q| q.saturating_sub(33)).collect::<Vec<u8>>(),
        ))
        .build()
}

