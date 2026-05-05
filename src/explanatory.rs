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

        for (group, sv_breaks) in groups.iter_mut() {
            ExtendedSeed::prune_repetitive_seeds(group, sv_breaks, 10, Some(self.reference));
            ExtendedSeed::extend_and_trim(group, sv_breaks, query, self.reference);
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
            let mut s = 0;
            for (j, (group, sv_breaks)) in groups.iter().enumerate() {
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
                let mut segment_score = 0.0;
                for (i, seed) in group.iter().enumerate() {
                    segment_score += seed.weight();
                    if i < sv_breaks.len() {
                        if !sv_breaks[i] {
                            if let Some((weight, _)) = seed.edge_penalty(&group[i + 1]) {
                                segment_score += weight;
                            }
                        }
                    }
                    if i < sv_breaks.len() && sv_breaks[i] {
                        log::info!("group {}, segment {}: score {:.1}", j, s, segment_score);
                        s += 1;
                        segment_score = 0.0;
                    }
                    k += 1;

                    log::info!(
                        "group {}, segment {}, seed {}: length: {}, weight {:.1}, diagonal {}",
                        j,
                        s,
                        k,
                        seed.length(),
                        seed.weight(),
                        seed.diagonal()
                    );

                    let b = seed.read_start();
                    let e = seed.read_end();
                    let q_str: String = quality[b..e].iter().map(|q| *q as char).collect();
                    let sa_parts: Vec<String> = (0..g)
                        .filter(|v| *v != i)
                        .map(|v| alts[v].clone())
                        .collect();
                    let tags = vec![
                        (String::from("XG"), TagValue::Int(j as i64)),
                        (String::from("XR"), TagValue::Int(s as i64)),
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
                log::info!("group {}, segment {}: score {:.1}", j, s, segment_score);
                s += 1;
            }
        }

        let mut explanations: Vec<Vec<Segment>> = Vec::new();

        for (i, (group, sv_breaks)) in groups.iter_mut().enumerate() {

            let mut gaps = ExtendedSeed::align_gaps(
                group,
                sv_breaks,
                query,
                self.reference,
                &mut self.aligner,
            );

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
                        let first = &group[segment_start];
                        let last = &group[j];
                        let is_reverse = first.is_reverse();
                        // Forward-strand ref range: [leftmost, rightmost).
                        // For reverse seeds, last (highest read pos) has the lowest ref coord.
                        let (ref_start, ref_end) = if is_reverse {
                            (last.ref_start(), first.ref_start() + first.length())
                        } else {
                            (first.ref_start(), last.ref_start() + last.length())
                        };
                        segments.push(Segment {
                            alignment: Alignment::concat(&std::mem::take(&mut current_parts)),
                            chrom_id: first.ref_chrom_id(),
                            is_reverse,
                            fwd_read_start: first.read_start(),
                            fwd_read_end: last.read_end(),
                            ref_start,
                            ref_end,
                        });
                    }
                }
            }

            // Second pass: compute x-drop extensions for the outer ends of the
            // primary group (group 0) only, stitch them into the alignment, and
            // update the segment bounds so they remain self-consistent.
            let n_segs = segments.len();
            for (seg_idx, segment) in segments.iter_mut().enumerate() {
                if i != 0 {
                    break;
                }
                let is_reverse = segment.is_reverse;
                let chrom_id = segment.chrom_id;

                let fwd_left_budget = if seg_idx == 0 {
                    segment.fwd_read_start
                } else {
                    0
                };
                let fwd_right_budget = if seg_idx == n_segs - 1 {
                    query_len - segment.fwd_read_end
                } else {
                    0
                };
                let (strand_left_budget, strand_right_budget) = if is_reverse {
                    (fwd_right_budget, fwd_left_budget)
                } else {
                    (fwd_left_budget, fwd_right_budget)
                };

                let chrom_len = self.reference.chrom_length(chrom_id) as usize;
                const REF_EXTENSION_PAD: usize = 10_000;

                // All extensions are computed in the same coordinate system as the
                // main CIGAR: forward query vs RC ref (for reverse) or forward ref
                // (for forward). This keeps the stitched CIGAR internally consistent
                // so validate() can walk query and ref with a single forward pass.
                //
                // For forward strand:
                //   left ext:  query[fwd_read_start-avail..fwd_read_start] vs fwd ref to the left
                //   right ext: query[fwd_read_end..fwd_read_end+avail]     vs fwd ref to the right
                //
                // For reverse strand (main CIGAR uses query[fwd_read_start..fwd_read_end] vs
                // RC(ref[ref_start..ref_end]), so RC-ref position 0 == genome[ref_end-1]):
                //   right ext (prepended): query[fwd_read_start-avail..fwd_read_start]
                //                         vs RC(ref[ref_end..ref_end+avail+pad])  — extends RC-ref leftward
                //   left ext  (appended):  query[fwd_read_end..fwd_read_end+avail]
                //                         vs RC(ref[ref_start-avail-pad..ref_start]) — extends RC-ref rightward

                let (fwd_left_ext, fwd_right_ext) = if is_reverse {
                    // fwd_left_budget  = unclipped bases to the LEFT of fwd_read_start (strand-right)
                    // fwd_right_budget = unclipped bases to the RIGHT of fwd_read_end  (strand-left)
                    let right_ext = if fwd_left_budget > 0 && segment.fwd_read_start > 0 {
                        let available = segment.fwd_read_start.min(fwd_left_budget);
                        let read_slice =
                            &query[segment.fwd_read_start - available..segment.fwd_read_start];
                        let ref_end_ext =
                            (segment.ref_end + available + REF_EXTENSION_PAD).min(chrom_len);
                        let rc_ref: Vec<u8> = self
                            .reference
                            .get_seq(chrom_id, segment.ref_end, ref_end_ext)
                            .iter()
                            .rev()
                            .map(|&b| complement(b))
                            .collect();
                        self.aligner.extend_left(read_slice, &rc_ref).ok()
                    } else {
                        None
                    };
                    let left_ext = if fwd_right_budget > 0 && segment.fwd_read_end < query_len {
                        let available = (query_len - segment.fwd_read_end).min(fwd_right_budget);
                        let read_slice =
                            &query[segment.fwd_read_end..segment.fwd_read_end + available];
                        let ref_start_ext = segment
                            .ref_start
                            .saturating_sub(available + REF_EXTENSION_PAD);
                        let rc_ref: Vec<u8> = self
                            .reference
                            .get_seq(chrom_id, ref_start_ext, segment.ref_start)
                            .iter()
                            .rev()
                            .map(|&b| complement(b))
                            .collect();
                        self.aligner.extend_right(read_slice, &rc_ref).ok()
                    } else {
                        None
                    };
                    (left_ext, right_ext)
                } else {
                    let left_ext = if strand_left_budget > 0 && segment.fwd_read_start > 0 {
                        let available = segment.fwd_read_start.min(strand_left_budget);
                        let read_slice =
                            &query[segment.fwd_read_start - available..segment.fwd_read_start];
                        let ref_start_ext = segment
                            .ref_start
                            .saturating_sub(available + REF_EXTENSION_PAD);
                        let ref_slice =
                            self.reference
                                .get_seq(chrom_id, ref_start_ext, segment.ref_start);
                        self.aligner.extend_left(read_slice, ref_slice).ok()
                    } else {
                        None
                    };
                    let right_ext = if strand_right_budget > 0 && segment.fwd_read_end < query_len {
                        let available = (query_len - segment.fwd_read_end).min(strand_right_budget);
                        let read_slice =
                            &query[segment.fwd_read_end..segment.fwd_read_end + available];
                        let ref_end_ext =
                            (segment.ref_end + available + REF_EXTENSION_PAD).min(chrom_len);
                        let ref_slice =
                            self.reference
                                .get_seq(chrom_id, segment.ref_end, ref_end_ext);
                        self.aligner.extend_right(read_slice, ref_slice).ok()
                    } else {
                        None
                    };
                    (left_ext, right_ext)
                };

                // Update segment bounds. For both strands:
                //   left ext:  prepends to fwd_read_start (shrinks it) and shrinks ref_start
                //   right ext: appends  to fwd_read_end   (grows it)   and grows   ref_end
                // For reverse the ref directions are:
                //   fwd_right_ext (stored as right_ext): shrinks fwd_read_start, grows ref_end
                //   fwd_left_ext  (stored as left_ext):  grows fwd_read_end,     shrinks ref_start
                let left_qlen = fwd_left_ext.as_ref().map_or(0, |e| e.query_length());
                let left_ref = fwd_left_ext.as_ref().map_or(0, |e| e.reference_consumed());
                let right_qlen = fwd_right_ext.as_ref().map_or(0, |e| e.query_length());
                let right_ref = fwd_right_ext.as_ref().map_or(0, |e| e.reference_consumed());

                if is_reverse {
                    segment.fwd_read_end += left_qlen;
                    segment.ref_start -= left_ref;
                    segment.fwd_read_start -= right_qlen;
                    segment.ref_end += right_ref;
                } else {
                    segment.fwd_read_start -= left_qlen;
                    segment.ref_start -= left_ref;
                    segment.fwd_read_end += right_qlen;
                    segment.ref_end += right_ref;
                }

                // Stitch extensions into the alignment in forward-query order:
                // [right_ext (prepended), main, left_ext (appended)] for reverse,
                // [left_ext,              main, right_ext]            for forward.
                let mut parts: Vec<Alignment> = Vec::new();
                if is_reverse {
                    if let Some(e) = fwd_right_ext {
                        parts.push(e);
                    }
                    parts.push(std::mem::take(&mut segment.alignment));
                    if let Some(e) = fwd_left_ext {
                        parts.push(e);
                    }
                } else {
                    if let Some(e) = fwd_left_ext {
                        parts.push(e);
                    }
                    parts.push(std::mem::take(&mut segment.alignment));
                    if let Some(e) = fwd_right_ext {
                        parts.push(e);
                    }
                }
                segment.alignment = Alignment::concat(&parts);
                segment.alignment.normalize();
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
                    "{}\t{}\t{}\t{}\t{}\t{}\t{:.1}",
                    name,
                    i,
                    segmentss.len(),
                    total_score,
                    query_coverage,
                    query_len,
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
                    let chrom_name = self.reference.chrom_name(segment.chrom_id);
                    let strand = if segment.is_reverse { "-" } else { "+" };
                    let summary_cigar = segment.alignment.summary_cigar(
                        segment.fwd_read_start,
                        segment.fwd_read_end,
                        query_len,
                        segment.is_reverse,
                    );
                    let nm = segment.alignment.mismatch_count();
                    format!(
                        "{},{},{},{},255,{}",
                        chrom_name,
                        segment.sam_pos(),
                        strand,
                        summary_cigar,
                        nm
                    )
                })
                .collect();

            // Pick the best segment (longest query span) as the representative.
            let best_seg_idx = segments
                .iter()
                .enumerate()
                .max_by_key(|(_, seg)| seg.fwd_read_end - seg.fwd_read_start)
                .map(|(idx, _)| idx)
                .unwrap_or(0);

            for (seg_idx, segment) in segments.iter().enumerate() {
                let is_reverse = segment.is_reverse;
                let chrom_id = segment.chrom_id;

                let strand_read_start = segment.strand_read_start(query_len);
                let strand_read_end = segment.strand_read_end(query_len);

                // Validate: all CIGARs (including extensions) use forward query vs
                // RC ref (reverse) or forward ref (forward), so validation is uniform.
                {
                    let ref_slice: Vec<u8> = if is_reverse {
                        self.reference
                            .get_seq(chrom_id, segment.ref_start, segment.ref_end)
                            .iter()
                            .rev()
                            .map(|&b| complement(b))
                            .collect()
                    } else {
                        self.reference
                            .get_seq(chrom_id, segment.ref_start, segment.ref_end)
                            .to_vec()
                    };
                    let query_seq = &query[segment.fwd_read_start..segment.fwd_read_end];
                    if let Err(e) = segment.alignment.validate(&ref_slice, query_seq, 0) {
                        let chrom_name = self.reference.chrom_name(chrom_id);
                        let strand = if is_reverse { "-" } else { "+" };
                        log::error!(
                            "VALIDATION FAILED: group {} seg {} ({} {}:{}-{} {}): {}",
                            i,
                            seg_idx,
                            name,
                            chrom_name,
                            segment.ref_start,
                            segment.ref_end,
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

                // Build CIGAR: primary gets soft clips, supplementary/secondary
                // get hard clips (and a truncated SEQ/QUAL).
                let clip_kind = if is_primary {
                    Kind::SoftClip
                } else {
                    Kind::HardClip
                };

                let mut cigar = Vec::new();
                if strand_read_start > 0 {
                    cigar.push(Op::new(clip_kind, strand_read_start));
                }
                if is_reverse {
                    for &op in segment.alignment.cigar.iter().rev() {
                        cigar.push(op);
                    }
                } else {
                    cigar.extend_from_slice(&segment.alignment.cigar);
                }
                if strand_read_end < query_len {
                    cigar.push(Op::new(clip_kind, query_len - strand_read_end));
                }
                let noodles_cigar: noodles::sam::alignment::record_buf::Cigar =
                    cigar.iter().copied().collect();

                let mapq = if i == 0 {
                    let alt = explanations.get(1).map(|s| s.as_slice()).unwrap_or(&[]);
                    compute_mapq(segment, alt)
                } else {
                    0u8
                };

                // SEQ/QUAL: primary emits the full strand sequence; non-primary
                // emits only the aligned portion in strand-space coordinates.
                let (strand_seq, strand_qual) = if is_reverse {
                    (&query_rc[..], &quality_rc[..])
                } else {
                    (query, quality)
                };
                let (out_seq, out_qual) = if is_primary {
                    (strand_seq, strand_qual)
                } else {
                    (
                        &strand_seq[strand_read_start..strand_read_end],
                        &strand_qual[strand_read_start..strand_read_end],
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
                    segment.sam_pos(),
                    mapq,
                    noodles_cigar,
                    None,
                    None,
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

// Assemble segments: each segment is a maximal run of colinear seeds.
// A None gap (or end of group) terminates the current segment.
// X-drop extensions are computed in a second pass and the bounds updated
// in place, so the struct is always self-consistent for validation and
// SAM emission without referencing the original seeds.
struct Segment {
    alignment: Alignment,
    chrom_id: usize,
    is_reverse: bool,
    // Forward-strand query range covered by the alignment (half-open).
    fwd_read_start: usize,
    fwd_read_end: usize,
    // Forward-strand ref range covered by the alignment (half-open).
    ref_start: usize,
    ref_end: usize,
}

impl Segment {
    // Strand-space read start (index into query_rc for reverse, query for forward).
    fn strand_read_start(&self, query_len: usize) -> usize {
        if self.is_reverse {
            query_len - self.fwd_read_end
        } else {
            self.fwd_read_start
        }
    }
    // Strand-space read end.
    fn strand_read_end(&self, query_len: usize) -> usize {
        if self.is_reverse {
            query_len - self.fwd_read_start
        } else {
            self.fwd_read_end
        }
    }
    // SAM POS: leftmost 1-based ref position.
    fn sam_pos(&self) -> usize {
        self.ref_start + 1
    }
    // Aligned query length.
    fn aligned_len(&self) -> usize {
        self.fwd_read_end - self.fwd_read_start
    }
}

// Takes two slices of segments (sorted by fwd_read_start) and returns all pairwise
// intersections as (lhs_index, rhs_index, overlap_start, overlap_end).
fn intersections(lhs: &[Segment], rhs: &[Segment]) -> Vec<(usize, usize, usize, usize)> {
    let mut result = Vec::new();
    let mut rhs_start = 0; // lowest rhs index that can still overlap the current lhs

    for (i, l) in lhs.iter().enumerate() {
        // Advance rhs_start past segments that end before l begins.
        while rhs_start < rhs.len() && rhs[rhs_start].fwd_read_end <= l.fwd_read_start {
            rhs_start += 1;
        }
        // Walk rhs from rhs_start while segments can still overlap l.
        let mut j = rhs_start;
        while j < rhs.len() && rhs[j].fwd_read_start < l.fwd_read_end {
            let overlap_start = l.fwd_read_start.max(rhs[j].fwd_read_start);
            let overlap_end = l.fwd_read_end.min(rhs[j].fwd_read_end);
            if overlap_end > overlap_start {
                result.push((i, j, overlap_start, overlap_end));
            }
            j += 1;
        }
    }

    result
}

// Compute MAPQ for a group-0 segment given the group-1 alternative segments.
//
// MAPQ = 10 * (alt_divergence - seg0_divergence) / ln(10)  ≈ 4.343 * score_diff
//
// This is the log-likelihood ratio between the best and alternative alignments
// under a model where each error costs 1 unit.  Using raw divergence differences
// (not rates) means longer alignments accumulate larger score differences at the
// same per-base error rate, reflecting that they are stronger evidence for the
// correct mapping.  Uncovered bases in the alternative are treated as fully
// divergent (1 unit each).
fn compute_mapq(seg0: &Segment, alt_segments: &[Segment]) -> u8 {
    let len = seg0.aligned_len();
    if len == 0 {
        return 0;
    }

    let mut alt_divergence = 0.0f64;
    let mut alt_covered = 0usize;

    for (_, j, overlap_start, overlap_end) in
        intersections(std::slice::from_ref(seg0), alt_segments)
    {
        let overlap_len = overlap_end - overlap_start;
        let alt_len = alt_segments[j].aligned_len();
        if alt_len == 0 {
            continue;
        }
        let scaled = alt_segments[j].alignment.divergence.0 * (overlap_len as f64 / alt_len as f64);
        alt_divergence += scaled;
        alt_covered += overlap_len;
    }

    alt_divergence += len.saturating_sub(alt_covered) as f64;

    let score_diff = alt_divergence - seg0.alignment.divergence.0;
    ((score_diff * 10.0 / std::f64::consts::LN_10)
        .clamp(0.0, 60.0)
        .round()) as u8
}
