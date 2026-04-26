use noodles::sam::alignment::{
    record::{
        Flags,
        cigar::{Op, op::Kind},
        data::field::Tag,
    },
    record_buf::{Data, data::field::Value},
};

use crate::{
    Aligner, AlignerBuilder,
    align::Alignment,
    index::Index,
    reads::{
        builder::{build_record, build_unmapped_record},
        extended::ExtendedSeed,
    },
    reference::InMemoryReference,
    seeding::SeedCollector,
    utils::sequence::{complement, reverse_complement_into},
    writer::AlignmentWriter,
};

pub struct ExplanatoryAlignerBuilder<'a, const K: usize, const S: usize> {
    reference: &'a InMemoryReference,
    index: &'a Index<K, S>,
    writer: &'a AlignmentWriter,
}

impl<'a, const K: usize, const S: usize> AlignerBuilder<'a, K, S>
    for ExplanatoryAlignerBuilder<'a, K, S>
{
    type AlignerType = ExplanatoryAligner<'a, K, S>;

    fn new(
        reference: &'a InMemoryReference,
        index: &'a Index<K, S>,
        writer: &'a AlignmentWriter,
    ) -> Self {
        Self {
            reference,
            index,
            writer,
        }
    }

    fn build(self) -> Self::AlignerType {
        ExplanatoryAligner {
            reference: self.reference,
            index: self.index,
            writer: self.writer,
            seeder: SeedCollector::new(),
            aligner: crate::align::DpAligner::new(),
            all_seeds: Vec::new(),
        }
    }
}

pub struct ExplanatoryAligner<'a, const K: usize, const S: usize> {
    reference: &'a InMemoryReference,
    index: &'a Index<K, S>,
    writer: &'a AlignmentWriter,
    seeder: SeedCollector,
    aligner: crate::align::DpAligner,
    all_seeds: Vec<ExtendedSeed>,
}

impl<'a, const K: usize, const S: usize> Aligner<'a, K, S> for ExplanatoryAligner<'a, K, S> {
    fn align(&mut self, name: &str, query: &[u8], quality: &[u8]) -> std::io::Result<()> {
        let query_len = query.len();
        let mut query_rc = Vec::with_capacity(query_len);
        reverse_complement_into(query, &mut query_rc);

        let rc_quality: Vec<u8> = quality.iter().rev().copied().collect();

        self.all_seeds.clear();

        self.seeder
            .gather_seeds_batched::<K, S>(query, false, self.index, self.reference, name);
        self.all_seeds.extend(
            self.seeder
                .hits
                .iter()
                .map(|seed| ExtendedSeed::from_seed_hit(seed, false, query_len)),
        );
        self.seeder
            .gather_seeds_batched::<K, S>(&query_rc, true, self.index, self.reference, name);
        self.all_seeds.extend(
            self.seeder
                .hits
                .iter()
                .map(|seed| ExtendedSeed::from_seed_hit(seed, true, query_len)),
        );

        // Simplify seeds by merging overlapping ones on the same diagonal
        ExtendedSeed::simplify_seeds(&mut self.all_seeds);

        let mut groups = ExtendedSeed::form_explanatory_groups(&self.all_seeds);

        if groups.is_empty() {
            let record = build_unmapped_record(name, query, quality);
            self.writer.write_record(&record).expect("write failed");
            return Ok(());
        }

        // Assemble segments: each segment is a maximal run of colinear seeds.
        // A None gap (or end of group) terminates the current segment.
        struct Segment {
            first_seed: ExtendedSeed,
            last_seed: ExtendedSeed,
            alignment: Alignment,
        }

        let mut explanations: Vec<Vec<Segment>> = Vec::new();

        for i in 0..groups.len() {
            let group = &mut groups[i];

            ExtendedSeed::extend_and_trim(group, query, self.reference);

            let mut gaps =
                ExtendedSeed::align_gaps(group, query, self.reference, &mut self.aligner);

            let n = group.len();
            let mut segments: Vec<Segment> = Vec::new();
            let mut current_parts: Vec<Alignment> = Vec::new();
            let mut segment_start = 0;
            for j in 0..n {
                if current_parts.is_empty() {
                    segment_start = j;
                }
                current_parts.push(group[j].to_alignment());
                match gaps.get_mut(j) {
                    Some(Some(aln)) => {
                        let mut tmp = Alignment::default();
                        std::mem::swap(aln, &mut tmp);
                        current_parts.push(tmp);
                    }
                    None | Some(None) => {
                        segments.push(Segment {
                            first_seed: group[segment_start].clone(),
                            last_seed: group[j].clone(),
                            alignment: Alignment::concat(&std::mem::take(&mut current_parts)),
                        });
                    }
                }
            }
            explanations.push(segments);
        }

        if false {
            for (i, segmentss) in explanations.iter().enumerate() {
                let mut query_coverage = 0usize;
                let mut total_score = 0.0f64;
                for segment in segmentss.iter() {
                    query_coverage += segment.alignment.query_length();
                    total_score += segment.alignment.divergence.0;
                }
                // Treat the uncovered portion of the query as a deletion
                let missing_coverage = query_len - query_coverage;
                total_score += missing_coverage as f64;
                let coverage_pct = 100.0 * (query_coverage as f64) / (query_len as f64);
                println!(
                    "Group {}: {} segments, total score {}, query coverage {:.1}%",
                    i,
                    segmentss.len(),
                    total_score,
                    coverage_pct
                );
            }
        }

        for (i, segments) in explanations.iter().enumerate() {
            // Build SA tag summaries for each segment so we can cross-reference.
            // Format per SAM spec: rname,pos,strand,CIGAR,mapQ,NM
            let sa_entries: Vec<String> = segments
                .iter()
                .map(|segment| {
                    let first = &segment.first_seed;
                    let last = &segment.last_seed;
                    let is_reverse = first.is_reverse();
                    let chrom_id = first.ref_chrom_id();
                    let chrom_name = self.reference.chrom_name(chrom_id);

                    let ref_pos = if is_reverse {
                        last.ref_start() + 1
                    } else {
                        first.ref_start() + 1
                    };
                    let strand = if is_reverse { "-" } else { "+" };
                    let summary_cigar = segment.alignment.summary_cigar(
                        first.read_start(),
                        last.read_end(),
                        query_len,
                        is_reverse,
                    );
                    let nm = segment.alignment.mismatch_count();

                    format!(
                        "{},{},{},{},255,{}",
                        chrom_name, ref_pos, strand, summary_cigar, nm
                    )
                })
                .collect();

            // Pick the best segment (longest query span) as the representative.
            let best_seg_idx = segments
                .iter()
                .enumerate()
                .max_by_key(|(_, seg)| seg.alignment.query_length())
                .map(|(idx, _)| idx)
                .unwrap_or(0);

            for (seg_idx, segment) in segments.iter().enumerate() {
                let first = &segment.first_seed;
                let last = &segment.last_seed;
                let is_reverse = first.is_reverse();
                let chrom_id = first.ref_chrom_id();

                // SAM POS: leftmost reference position (1-based).
                // Forward: first seed has the leftmost ref position.
                // Reverse: last seed has the leftmost ref position (ref
                // decreases as read advances in a colinear reverse segment).
                let ref_pos = if is_reverse {
                    last.ref_start() + 1
                } else {
                    first.ref_start() + 1
                };

                // Read range covered by this segment.
                let seg_read_start = first.read_start();
                let seg_read_end = last.read_end();

                // Validate the segment alignment against the reference and query.
                if false {
                    let (ref_begin, ref_end) = if is_reverse {
                        (last.ref_start(), first.ref_start() + first.length())
                    } else {
                        (first.ref_start(), last.ref_start() + last.length())
                    };

                    let ref_slice: Vec<u8> = if is_reverse {
                        self.reference
                            .get_seq(chrom_id, ref_begin, ref_end)
                            .iter()
                            .rev()
                            .map(|&b| complement(b))
                            .collect()
                    } else {
                        self.reference
                            .get_seq(chrom_id, ref_begin, ref_end)
                            .to_vec()
                    };

                    // The alignment was built against seq (forward read),
                    // regardless of strand.
                    let query_seq = &query[seg_read_start..seg_read_end];

                    if let Err(e) = segment.alignment.validate(&ref_slice, query_seq, 0) {
                        let chrom_name = self.reference.chrom_name(chrom_id);
                        let strand = if is_reverse { "-" } else { "+" };
                        log::error!(
                            "VALIDATION FAILED: group {} seg {} ({} {}:{}-{} {}): {}",
                            i,
                            seg_idx,
                            name,
                            chrom_name,
                            ref_begin,
                            ref_end,
                            strand,
                            e
                        );
                    }
                }

                // Flags
                let mut flags = Flags::empty();
                if is_reverse {
                    flags |= Flags::REVERSE_COMPLEMENTED;
                }
                if i > 0 {
                    flags |= Flags::SECONDARY;
                }
                if seg_idx != best_seg_idx {
                    flags |= Flags::SUPPLEMENTARY;
                }

                let is_primary = i == 0 && seg_idx == best_seg_idx;

                // Build CIGAR: primary gets soft clips, secondary/supplementary
                // get hard clips (and a truncated SEQ/QUAL).
                let clip_kind = if is_primary {
                    Kind::SoftClip
                } else {
                    Kind::HardClip
                };
                let mut cigar = Vec::new();
                if is_reverse {
                    // The alignment was built as seq vs rc(ref), but SAM
                    // convention is rc_seq vs forward_ref.  We reverse the
                    // CIGAR and use rc_seq coordinates for clipping:
                    //   rc_start = seq_len - seg_read_end
                    //   rc_end   = seq_len - seg_read_start
                    let rc_start = query_len - seg_read_end;
                    let rc_end = query_len - seg_read_start;
                    if rc_start > 0 {
                        cigar.push(Op::new(clip_kind, rc_start));
                    }
                    for &op in segment.alignment.cigar.iter().rev() {
                        cigar.push(op);
                    }
                    if rc_end < query_len {
                        cigar.push(Op::new(clip_kind, query_len - rc_end));
                    }
                } else {
                    if seg_read_start > 0 {
                        cigar.push(Op::new(clip_kind, seg_read_start));
                    }
                    cigar.extend_from_slice(&segment.alignment.cigar);
                    if seg_read_end < query_len {
                        cigar.push(Op::new(clip_kind, query_len - seg_read_end));
                    }
                }
                let noodles_cigar: noodles::sam::alignment::record_buf::Cigar =
                    cigar.iter().copied().collect();

                // SEQ/QUAL: primary emits the full read; secondary/supplementary
                // emit only the aligned portion.
                let (strand_seq, strand_qual) = if is_reverse {
                    (&query_rc[..], &rc_quality[..])
                } else {
                    (query, quality)
                };
                let (out_seq, out_qual) = if is_primary {
                    (strand_seq, strand_qual)
                } else if is_reverse {
                    // Non-primary reverse: use rc_seq coordinates.
                    let rc_start = query_len - seg_read_end;
                    let rc_end = query_len - seg_read_start;
                    (&query_rc[rc_start..rc_end], &rc_quality[rc_start..rc_end])
                } else {
                    (
                        &strand_seq[seg_read_start..seg_read_end],
                        &strand_qual[seg_read_start..seg_read_end],
                    )
                };

                // Build SA tag: list all OTHER segments in this group.
                let sa_value: String = sa_entries
                    .iter()
                    .enumerate()
                    .filter(|&(k, _)| k != seg_idx)
                    .map(|(_, entry)| entry.as_str())
                    .collect::<Vec<_>>()
                    .join(";");

                let data: Data = if segments.len() > 1 {
                    vec![(
                        Tag::try_from(*b"SA").unwrap(),
                        Value::from(sa_value.as_str()),
                    )]
                    .into_iter()
                    .collect()
                } else {
                    Data::default()
                };

                let record = build_record(
                    name,
                    flags,
                    chrom_id,
                    ref_pos,
                    255, // mapq placeholder
                    noodles_cigar,
                    None, // mate_ref_id
                    None, // mate_pos
                    out_seq,
                    out_qual,
                    data,
                );
                self.writer.write_record(&record).expect("write failed");
            }

            if i > 2 {
                break;
            }
        }

        Ok(())
    }

    fn finish(self) -> std::io::Result<()> {
        self.writer.finish()?;
        Ok(())
    }
}
