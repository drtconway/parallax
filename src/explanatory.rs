use std::sync::{Mutex, OnceLock};

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
    config,
    index::Index,
    reads::{
        builder::{build_record, build_unmapped_record},
        extended::{ExtendedSeed, ExtendedSeedDumpItem, TagValue},
    },
    reference::InMemoryReference,
    seeding::SeedCollector,
    utils::{
        dump::DumpItem,
        sequence::{complement, reverse_complement_into},
    },
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

        let quality_rc: Vec<u8> = quality.iter().rev().copied().collect();

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

        if !config::get().seeding.debug_seeds_sam.is_empty() {
            static SEED_DUMPER: OnceLock<Mutex<(std::fs::File, usize)>> = OnceLock::new();

            let path = config::get().seeding.debug_seeds_sam.clone();
            let mut dump = SEED_DUMPER
                .get_or_init(|| {
                    Mutex::new((
                        std::fs::File::create(path).expect("failed to create dumpt file"),
                        0,
                    ))
                })
                .lock()
                .unwrap();

            if dump.1 == 0 {
                ExtendedSeedDumpItem::write_header(self.reference, &mut dump.0);
            }

            let query_str: String = query.iter().map(|c| *c as char).collect();

            let n = query.len();

            for (i, seed) in self.all_seeds.iter().enumerate() {
                // SEQ is always taken from the forward-strand query at the seed's
                // forward-strand coordinates. For FLAG=16 records IGV RC's the SEQ
                // to compare against the reference, and RC(query_rc[b..e]) =
                // query[read_start..read_end], so this is the correct orientation.
                let b = seed.read_start();
                let e = seed.read_end();
                let q_str: String = quality[b..e].iter().map(|q| *q as char).collect();
                let item = ExtendedSeedDumpItem::from((
                    self.reference,
                    name,
                    n,
                    i,
                    seed,
                    &query_str[b..e],
                    &q_str as &str,
                ));
                item.write(&mut dump.0);
                dump.1 += 1;
            }
        }

        // Simplify seeds by merging overlapping ones on the same diagonal
        ExtendedSeed::simplify_seeds(&mut self.all_seeds);

        let mut groups = ExtendedSeed::form_explanatory_groups(&self.all_seeds);

        if groups.is_empty() {
            let record = build_unmapped_record(name, query, quality);
            self.writer.write_record(&record).expect("write failed");
            return Ok(());
        }

        if !config::get().seeding.debug_chains_sam.is_empty() {
            static CHAIN_DUMPER: OnceLock<Mutex<(std::fs::File, usize)>> = OnceLock::new();

            let path = config::get().seeding.debug_chains_sam.clone();
            let mut dump = CHAIN_DUMPER
                .get_or_init(|| {
                    Mutex::new((
                        std::fs::File::create(path).expect("failed to create dumpt file"),
                        0,
                    ))
                })
                .lock()
                .unwrap();

            if dump.1 == 0 {
                ExtendedSeedDumpItem::write_header(self.reference, &mut dump.0);
            }

            let query_str: String = query.iter().map(|c| *c as char).collect();

            let n = query.len();
            let mut k = 0;
            for (j, group) in groups.iter().enumerate() {
                let alts: Vec<String> = group
                    .iter()
                    .map(|seed| {
                        // SA tag clips must be in strand space (RC space for reverse seeds).
                        let (left_clip, right_clip) = if seed.is_reverse() {
                            (n - seed.read_end(), seed.read_start())
                        } else {
                            (seed.read_start(), n - seed.read_end())
                        };
                        let left = if left_clip > 0 {
                            format!("{}S", left_clip)
                        } else {
                            String::new()
                        };
                        let right = if right_clip > 0 {
                            format!("{}S", right_clip)
                        } else {
                            String::new()
                        };
                        let chrom = self.reference.chrom_name(seed.ref_chrom_id());
                        let strand = if seed.is_reverse() { "-" } else { "+" };
                        let mapq = (seed.weight().floor() as i32).min(200);
                        format!(
                            "{},{},{},{}{}={},{},0;",
                            chrom,
                            seed.ref_start() + 1,
                            strand,
                            left,
                            seed.length(),
                            right,
                            mapq
                        )
                    })
                    .collect();
                let g = alts.len();
                for (i, seed) in group.iter().enumerate() {
                    k += 1;
                    let b = seed.read_start();
                    let e = seed.read_end();
                    let q_str: String = quality[b..e].iter().map(|q| *q as char).collect();
                    let sa_parts: Vec<String> = (0..g)
                        .filter(|v| *v != i)
                        .map(|v| alts[v].clone())
                        .collect();
                    let tags = vec![
                        (String::from("XG"), TagValue::Int(j as i64)),
                        (String::from("XS"), TagValue::Int(i as i64)),
                        (String::from("SA"), TagValue::Str(sa_parts.join(""))),
                    ];
                    let item = ExtendedSeedDumpItem::from((
                        self.reference,
                        name,
                        n,
                        k,
                        seed,
                        &query_str[b..e],
                        &q_str as &str,
                        tags,
                    ));
                    item.write(&mut dump.0);
                    dump.1 += 1;
                }
            }
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

                // ── X-drop end extensions ─────────────────────────────────────
                // Only the outer ends of the primary group are extended; secondary
                // groups' clipped ends are already represented by the primary
                // alignment, and the cost is not justified.
                let is_group_start = seg_idx == 0;
                let is_group_end = seg_idx == segments.len() - 1;

                let fwd_left_budget = if is_group_start && i == 0 {
                    seg_read_start
                } else {
                    0
                };
                let fwd_right_budget = if is_group_end && i == 0 {
                    query_len - seg_read_end
                } else {
                    0
                };

                // Translate to strand space: for reverse, forward-left is
                // strand-right and vice versa (mirrors seed_cluster.rs:135).
                let (strand_left_budget, strand_right_budget) = if is_reverse {
                    (fwd_right_budget, fwd_left_budget)
                } else {
                    (fwd_left_budget, fwd_right_budget)
                };

                // Strand-space read/ref boundaries for this segment.
                let (strand_read_start, strand_read_end, strand_ref_start, strand_ref_end) =
                    if is_reverse {
                        (
                            query_len - seg_read_end,
                            query_len - seg_read_start,
                            last.ref_start(),
                            first.ref_start() + first.length(),
                        )
                    } else {
                        (
                            seg_read_start,
                            seg_read_end,
                            first.ref_start(),
                            last.ref_start() + last.length(),
                        )
                    };

                let strand_seq_ext: &[u8] = if is_reverse { &query_rc } else { query };
                let chrom_len = self.reference.chrom_length(chrom_id) as usize;

                const REF_EXTENSION_PAD: usize = 10_000;

                let left_ext: Option<Alignment> = if strand_left_budget > 0 && strand_read_start > 0
                {
                    let available = strand_read_start.min(strand_left_budget);
                    let read_slice =
                        &strand_seq_ext[strand_read_start - available..strand_read_start];
                    let ref_start = strand_ref_start.saturating_sub(available + REF_EXTENSION_PAD);
                    let ref_slice = self
                        .reference
                        .get_seq(chrom_id, ref_start, strand_ref_start);
                    self.aligner.extend_left(read_slice, ref_slice).ok()
                } else {
                    None
                };

                let right_ext: Option<Alignment> = if strand_right_budget > 0
                    && strand_read_end < query_len
                {
                    let available = (query_len - strand_read_end).min(strand_right_budget);
                    let read_slice = &strand_seq_ext[strand_read_end..strand_read_end + available];
                    let ref_end = (strand_ref_end + available + REF_EXTENSION_PAD).min(chrom_len);
                    let ref_slice = self.reference.get_seq(chrom_id, strand_ref_end, ref_end);
                    self.aligner.extend_right(read_slice, ref_slice).ok()
                } else {
                    None
                };

                let ext_left_qlen = left_ext.as_ref().map_or(0, |e| e.query_length());
                let ext_right_qlen = right_ext.as_ref().map_or(0, |e| e.query_length());
                let ref_start_adj = left_ext.as_ref().map_or(0, |e| e.reference_consumed());

                // Apply the left-extension reference adjustment to SAM POS.
                let ref_pos = ref_pos - ref_start_adj;

                // Validate the segment alignment against the reference and query.
                if true {
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
                // Extended strand-space boundaries after x-drop.
                let new_strand_start = strand_read_start - ext_left_qlen;
                let new_strand_end = strand_read_end + ext_right_qlen;

                // Both extend_left and extend_right return CIGARs in ref
                // left-to-right order, so they slot directly around the main
                // alignment without additional reversal — even for reverse strand
                // where the main CIGAR is reversed.
                let mut cigar = Vec::new();
                if new_strand_start > 0 {
                    cigar.push(Op::new(clip_kind, new_strand_start));
                }
                if let Some(e) = &left_ext {
                    cigar.extend_from_slice(&e.cigar);
                }
                if is_reverse {
                    for &op in segment.alignment.cigar.iter().rev() {
                        cigar.push(op);
                    }
                } else {
                    cigar.extend_from_slice(&segment.alignment.cigar);
                }
                if let Some(e) = &right_ext {
                    cigar.extend_from_slice(&e.cigar);
                }
                if new_strand_end < query_len {
                    cigar.push(Op::new(clip_kind, query_len - new_strand_end));
                }
                let noodles_cigar: noodles::sam::alignment::record_buf::Cigar =
                    cigar.iter().copied().collect();

                let mapq = if i > 0 {
                    0
                } else {
                    100 // XXX properly compute mapq
                };

                // SEQ/QUAL: primary emits the full strand sequence (soft clips
                // are included in SEQ); non-primary emits only the aligned
                // portion in strand-space coordinates, expanded to include any
                // x-drop extension.
                let (strand_seq, strand_qual) = if is_reverse {
                    (&query_rc[..], &quality_rc[..])
                } else {
                    (query, quality)
                };
                let (out_seq, out_qual) = if is_primary {
                    (strand_seq, strand_qual)
                } else {
                    (
                        &strand_seq[new_strand_start..new_strand_end],
                        &strand_qual[new_strand_start..new_strand_end],
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
                    mapq,
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
