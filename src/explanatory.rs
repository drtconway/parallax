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
    align::Alignment,
    aligner::{Aligner, AlignerBuilder},
    reads::{
        builder::{build_record, build_unmapped_record},
        extended::{
            ExtendedSeed, ExtendedSeedDumpItem, SeedFilter, ShortSingleSeedSegmentFilter, TagValue,
        },
    },
    seeding::SeedCollector,
    writer::AlignmentWriter,
};
use parallax::{
    config::{self, FilteringConfig, SeedingConfig},
    index::Index,
    reference::InMemoryReference,
    utils::{
        dump::DumpItem,
        sequence::{complement, reverse_complement_into},
        telemetry::{
            Recorder, RecorderExt, histogram::HistogramRecorder, summary::SimpleSummaryRecorder,
        },
    },
};

pub struct ExplanatoryAlignerBuilder<'a, const K: usize, const S: usize> {
    reference: &'a InMemoryReference,
    index: &'a Index<K, S>,
    writer: &'a AlignmentWriter,
    no_secondary: bool,
}

impl<'a, const K: usize, const S: usize> ExplanatoryAlignerBuilder<'a, K, S> {
    pub fn no_secondary(mut self, no_secondary: bool) -> Self {
        self.no_secondary = no_secondary;
        self
    }
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
            no_secondary: false,
        }
    }

    fn build(self) -> Self::AlignerType {
        let cfg = config::get();
        ExplanatoryAligner {
            reference: self.reference,
            index: self.index,
            writer: self.writer,
            seeder: SeedCollector::new(),
            aligner: crate::align::DpAligner::from_config(&cfg.alignment, &cfg.block_aligner),
            all_seeds: Vec::new(),
            no_secondary: self.no_secondary,
            seeding_cfg: cfg.seeding.clone(),
            filtering_cfg: cfg.filtering.clone(),
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
    no_secondary: bool,
    seeding_cfg: SeedingConfig,
    filtering_cfg: FilteringConfig,
}

impl<'a, const K: usize, const S: usize> Aligner<'a, K, S> for ExplanatoryAligner<'a, K, S> {
    fn align(&mut self, name: &str, query: &[u8], quality: &[u8]) -> std::io::Result<()> {
        let start = std::time::Instant::now();

        let query_len = query.len();
        let mut query_rc = Vec::with_capacity(query_len);
        reverse_complement_into(query, &mut query_rc);

        let quality_rc: Vec<u8> = quality.iter().rev().copied().collect();

        self.all_seeds.clear();

        self.seeder.gather_seeds_batched::<K, S>(
            query,
            false,
            self.index,
            self.reference,
            name,
            &self.seeding_cfg,
        );
        self.all_seeds.extend(
            self.seeder
                .hits
                .iter()
                .map(|seed| ExtendedSeed::from_seed_hit(seed, false, query_len)),
        );
        self.seeder.gather_seeds_batched::<K, S>(
            &query_rc,
            true,
            self.index,
            self.reference,
            name,
            &self.seeding_cfg,
        );
        self.all_seeds.extend(
            self.seeder
                .hits
                .iter()
                .map(|seed| ExtendedSeed::from_seed_hit(seed, true, query_len)),
        );

        if true {
            for seed in self.all_seeds.iter() {
                seed_length_recorder().record(seed.length());
            }
        }

        if false {
            let mut permutation = (0..self.all_seeds.len()).collect::<Vec<_>>();
            permutation.sort_by_key(|&i| {
                let seed = &self.all_seeds[i];
                (
                    seed.read_start(),
                    seed.read_end(),
                    seed.kmer_uniqueness(),
                    seed.ref_chrom_id(),
                    seed.ref_start(),
                )
            });

            let columns = [
                "read_start",
                "read_end",
                "length",
                "ref_chrom",
                "ref_start",
                "ref_end",
                "strand",
                "weight",
                "uniqueness",
                "frequency",
            ];
            println!("{}", columns.join("\t"));
            for i in permutation.into_iter() {
                let seed = &self.all_seeds[i];
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{}",
                    seed.read_start(),
                    seed.read_end(),
                    seed.length(),
                    self.reference.chrom_name(seed.ref_chrom_id()),
                    seed.ref_start(),
                    seed.ref_start() + seed.length(),
                    (if seed.is_reverse() { "-" } else { "+" }),
                    seed.weight(),
                    seed.kmer_uniqueness(),
                    seed.read_frequency(),
                );
            }
        }

        if !self.seeding_cfg.debug_seeds_sam.is_empty() {
            static SEED_DUMPER: OnceLock<Mutex<(std::fs::File, usize)>> = OnceLock::new();

            let path = self.seeding_cfg.debug_seeds_sam.clone();
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

        let mut groups = ExtendedSeed::form_explanatory_groups(&self.all_seeds, &self.seeding_cfg);

        if groups.is_empty() {
            let record = build_unmapped_record(name, query, quality);
            self.writer.write_record(&record).expect("write failed");
            return Ok(());
        }

        // Resolve ref and read overlaps introduced by the chaining DP before
        // any filtering or extension runs.
        for (group, sv_breaks) in groups.iter_mut() {
            ExtendedSeed::resolve_ref_overlaps(group, sv_breaks);
            ExtendedSeed::resolve_read_overlaps(group, sv_breaks);
        }

        let short_segment_filter = ShortSingleSeedSegmentFilter {
            min_length: self.seeding_cfg.min_single_seed_length,
        };
        for (group, sv_breaks) in groups.iter_mut() {
            ExtendedSeed::prune_repetitive_seeds(group, sv_breaks, 10, &self.seeding_cfg);
            if let Err(e) = ExtendedSeed::validate_chain(group, sv_breaks) {
                log::error!("{name}: chain invalid after prune_repetitive_seeds: {e}");
            }
            ExtendedSeed::extend_and_trim(
                name,
                group,
                sv_breaks,
                query,
                self.reference,
                &self.seeding_cfg,
            );
            if let Err(e) = ExtendedSeed::validate_chain(group, sv_breaks) {
                log::error!("{name}: chain invalid after extend_and_trim: {e}");
            }
            short_segment_filter.apply_filter(group, sv_breaks, &self.seeding_cfg);
            if group.is_empty() {
                log::debug!("{name}: group emptied by short_segment_filter");
            }
            if let Err(e) = ExtendedSeed::validate_chain(group, sv_breaks) {
                log::error!("{name}: chain invalid after short_segment_filter: {e}");
            }
        }

        if !self.seeding_cfg.debug_chains_sam.is_empty() {
            static CHAIN_DUMPER: OnceLock<Mutex<(std::fs::File, usize)>> = OnceLock::new();

            let path = self.seeding_cfg.debug_chains_sam.clone();
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
            let seeding_cfg = &self.seeding_cfg;
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
                            if let Some((weight, _)) = seed.edge_penalty(&group[i + 1], seeding_cfg)
                            {
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

                    log::debug!(
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
                        (
                            String::from("XJ"),
                            TagValue::Int(seed.read_frequency() as i64),
                        ),
                        (
                            String::from("XK"),
                            TagValue::Int(seed.kmer_uniqueness() as i64),
                        ),
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
                log::debug!("group {}, segment {}: score {:.1}", j, s, segment_score);
                s += 1;
            }
        }

        let mut explanations: Vec<Vec<Segment>> = Vec::new();

        for (i, (group, sv_breaks)) in groups.iter_mut().enumerate() {
            if group.is_empty() {
                continue;
            }

            let mut gaps = ExtendedSeed::align_gaps(
                name,
                group,
                sv_breaks,
                query,
                self.reference,
                &mut self.aligner,
            );

            // Identify spans where the direct bridging SW alignment scores better
            // than the segmented representation, and collapse them.
            {
                let align_params = crate::align::AlignParams::default();

                let is_colinear_pair = |a: &ExtendedSeed, b: &ExtendedSeed| -> bool {
                    if a.ref_chrom_id() != b.ref_chrom_id() || a.is_reverse() != b.is_reverse() {
                        return false;
                    }
                    if a.is_reverse() {
                        b.ref_start() + b.length() <= a.ref_start()
                    } else {
                        b.ref_start() >= a.ref_start() + a.length()
                    }
                };

                // Collect candidate spans: (l, r, bridging_alignment, score_improvement).
                // l = index of seed before the first sv_break in the span.
                // r = index of seed after the last sv_break in the span.
                // score_improvement = segmented_score - bridging_score (lower is better).
                struct Span {
                    l: usize,
                    r: usize,
                    bridging: Option<crate::align::Alignment>,
                    improvement: f64,
                }

                let sv_spans_tsv = !self.seeding_cfg.debug_sv_spans_tsv.is_empty() && i == 0;
                let mut tsv_file: Option<std::sync::MutexGuard<std::fs::File>> = if sv_spans_tsv {
                    static SV_SPAN_DUMPER: OnceLock<Mutex<std::fs::File>> = OnceLock::new();
                    let path = self.seeding_cfg.debug_sv_spans_tsv.clone();
                    let guard = SV_SPAN_DUMPER
                        .get_or_init(|| {
                            let mut f = std::fs::File::create(path).expect("failed to create sv spans file");
                            use std::io::Write;
                            writeln!(f, "read_name\tanchor_before_read_start\tanchor_before_read_end\tchrom\tanchor_before_ref_start\tanchor_before_ref_end\tnum_sv_breaks\tanchor_after_read_start\tanchor_after_read_end\tanchor_after_ref_start\tanchor_after_ref_end\tstrand\tread_gap\tref_gap\tsegmented_score\tbridging_score\tcollapsed").unwrap();
                            Mutex::new(f)
                        })
                        .lock()
                        .unwrap();
                    Some(guard)
                } else {
                    None
                };

                let mut candidates: Vec<Span> = Vec::new();
                let mut j = 0;
                while j < sv_breaks.len() {
                    if !sv_breaks[j] {
                        j += 1;
                        continue;
                    }
                    let l = j;
                    let mut r = j + 1;
                    while r < sv_breaks.len() && sv_breaks[r] {
                        r += 1;
                    }
                    let num_sv = sv_breaks[l..r].iter().filter(|&&b| b).count();
                    // group[l] is before, group[r] is after.
                    let before = &group[l];
                    let after = &group[r];
                    if num_sv > 1 && is_colinear_pair(before, after) {
                        let read_gap = after
                            .read_start()
                            .saturating_sub(before.read_start() + before.length());
                        let ref_gap = if before.is_reverse() {
                            before
                                .ref_start()
                                .saturating_sub(after.ref_start() + after.length())
                        } else {
                            after
                                .ref_start()
                                .saturating_sub(before.ref_start() + before.length())
                        };
                        const MAX_BRIDGE_GAP: usize = 10_000;
                        if read_gap > MAX_BRIDGE_GAP || ref_gap > MAX_BRIDGE_GAP {
                            j = r;
                            continue;
                        }
                        let segmented_score: f64 = (l + 1..r)
                            .map(|k| group[k].to_alignment().quality(&align_params).0)
                            .sum::<f64>()
                            + gaps[l..r]
                                .iter()
                                .filter_map(|g| g.as_ref())
                                .map(|aln| aln.quality(&align_params).0)
                                .sum::<f64>();
                        let bridging = ExtendedSeed::align_gap(
                            before,
                            after,
                            query,
                            self.reference,
                            &mut self.aligner,
                        );
                        let bridging_score = bridging
                            .as_ref()
                            .map(|aln| aln.quality(&align_params).0)
                            .unwrap_or(f64::NEG_INFINITY);
                        let improvement = segmented_score - bridging_score;
                        let will_collapse =
                            improvement < (num_sv as f64) * self.seeding_cfg.sv_penalty;
                        if let Some(ref mut file) = tsv_file.as_deref_mut() {
                            use std::io::Write;
                            let strand = if before.is_reverse() { "-" } else { "+" };
                            writeln!(
                                file,
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{}",
                                name,
                                before.read_start(), before.read_end(),
                                self.reference.chrom_name(before.ref_chrom_id()),
                                before.ref_start(), before.ref_start() + before.length(),
                                num_sv,
                                after.read_start(), after.read_end(),
                                after.ref_start(), after.ref_start() + after.length(),
                                strand,
                                read_gap as isize, ref_gap as isize,
                                segmented_score, bridging_score,
                                will_collapse,
                            ).unwrap();
                        }
                        if will_collapse {
                            candidates.push(Span {
                                l,
                                r,
                                bridging,
                                improvement,
                            });
                        }
                    }
                    j = r;
                }

                // Resolve overlapping candidates: linear scan, keep the better of each
                // overlapping pair (lower improvement = bigger gain).
                let mut accepted: Vec<Span> = Vec::new();
                for span in candidates {
                    if let Some(prev) = accepted.last_mut() {
                        // Spans overlap if span.l < prev.r (indices share seeds).
                        if span.l < prev.r {
                            if span.improvement < prev.improvement {
                                *prev = span;
                            }
                            continue;
                        }
                    }
                    accepted.push(span);
                }

                // Apply accepted spans right-to-left using pop-and-rebuild.
                for span in accepted.into_iter().rev() {
                    let l = span.l;
                    let r = span.r;

                    // Pop group[l+1..=r] and gaps[l..=r-1] off into temporaries,
                    // then push back only group[r] and the bridging gap.
                    // We handle group and gaps in lock-step.

                    // Drain seeds l+1..=r and gaps l..=r from the ends.
                    // Since l..r are all interior to the current vectors and we go
                    // right-to-left across spans, the tail indices are stable.
                    let tail_seeds: Vec<_> = group.drain(l + 1..).collect();
                    let tail_gaps: Vec<_> = gaps.drain(l..).collect();
                    let tail_breaks: Vec<_> = sv_breaks.drain(l..).collect();

                    // From the tails, keep only what comes after the span.
                    // tail_seeds[0..r-l-1] = interior seeds (discard)
                    // tail_seeds[r-l-1..] = seeds[r..] (keep)
                    // tail_gaps[0..r-l] = gaps[l..r] (discard)
                    // tail_gaps[r-l..] = gaps[r..] (keep)
                    // tail_breaks[0..r-l] = sv_breaks[l..r] (discard)
                    // tail_breaks[r-l..] = sv_breaks[r..] (keep)
                    let interior = r - l;
                    group.push(tail_seeds[interior - 1].clone()); // group[r] = after
                    group.extend_from_slice(&tail_seeds[interior..]);
                    gaps.push(span.bridging);
                    gaps.extend_from_slice(&tail_gaps[interior..]);
                    sv_breaks.push(false); // the collapsed span is now a colinear gap
                    sv_breaks.extend_from_slice(&tail_breaks[interior..]);
                }
            }

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

            merge_overlapping_segments(
                &mut segments,
                self.seeding_cfg.overlap_merge_max_identity_ratio,
                self.seeding_cfg.overlap_merge_min_forced,
                query_len,
            );

            segments.retain(|seg| {
                let read_bases = seg.fwd_read_end - seg.fwd_read_start;
                if read_bases < self.filtering_cfg.min_aligned_length {
                    log::debug!(
                        "{}: dropping segment with read span {} below threshold {}",
                        name,
                        read_bases,
                        self.filtering_cfg.min_aligned_length
                    );
                    false
                } else {
                    true
                }
            });
            if segments.is_empty() {
                log::debug!(
                    "{name}: all segments in group {i} dropped by min_aligned_length filter"
                );
                continue;
            }

            explanations.push(segments);
        }

        let mut segment_written = false;

        for (i, segments) in explanations.iter().enumerate() {
            assert!(!segments.is_empty());

            if i > 0 && self.no_secondary {
                break;
            }

            let mut query_covered = 0;
            let mut identities = 0;
            for segment in segments {
                let aln = &segment.alignment;
                query_covered += aln.query_consumed();
                identities += aln
                    .cigar
                    .iter()
                    .filter_map(|op| match op.kind() {
                        Kind::Match | Kind::SequenceMatch => Some(op.len()),
                        _ => None,
                    })
                    .sum::<usize>();
            }

            let query_coverage = (query_covered as f64) / (query_len as f64);
            if query_coverage < self.filtering_cfg.min_read_coverage {
                log::debug!(
                    "{}: skipping group {} with query coverage {:.1}% below threshold {:.1}%",
                    name,
                    i,
                    query_coverage * 100.0,
                    self.filtering_cfg.min_read_coverage * 100.0
                );
                continue;
            }

            let identity_fraction = if query_covered > 0 {
                (identities as f64) / (query_covered as f64)
            } else {
                0.0
            };
            if identity_fraction < self.filtering_cfg.min_identity {
                log::debug!(
                    "{}: skipping group {} with identity {:.1}% below threshold {:.1}%",
                    name,
                    i,
                    identity_fraction * 100.0,
                    self.filtering_cfg.min_identity * 100.0
                );
                continue;
            }

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

            // Scan the segments in this group for overlaps in the reference.
            // For each overlap, record the read span (in forward-read coords) for
            // each of the two segments involved, tagged with the ref region.
            // xo_entries[seg_idx] accumulates "read_start,read_end,chrom,ref_start,ref_end" strings.
            let mut xo_entries: Vec<Vec<String>> = vec![Vec::new(); segments.len()];
            for (a_idx, a) in segments.iter().enumerate() {
                for (b_idx, b) in segments.iter().enumerate().skip(a_idx + 1) {
                    if a.chrom_id != b.chrom_id {
                        continue;
                    }
                    let overlap_start = a.ref_start.max(b.ref_start);
                    let overlap_end = a.ref_end.min(b.ref_end);
                    if overlap_start >= overlap_end {
                        continue;
                    }
                    let chrom_name = self.reference.chrom_name(a.chrom_id);
                    for (seg_idx, seg) in [(a_idx, a), (b_idx, b)] {
                        if let Some((rs, re)) =
                            seg.read_range_for_ref_overlap(overlap_start, overlap_end)
                        {
                            xo_entries[seg_idx].push(format!(
                                "{},{},{},{},{}",
                                rs, re, chrom_name, overlap_start, overlap_end
                            ));
                        }
                    }
                }
            }

            // Pick the best segment (longest query span) as the representative.
            let best_seg_idx = segments
                .iter()
                .enumerate()
                .max_by_key(|(_, seg)| seg.fwd_read_end - seg.fwd_read_start)
                .map(|(idx, _)| idx)
                .unwrap_or(0);

            let num_segs = segments.len();
            if i == 0 {
                segment_count_recorder()
                    .record_value(parallax::utils::telemetry::Value::from(num_segs));
            }

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
                            "VALIDATION FAILED: group {} seg {} ({} {}:{}-{} {}): {}\nCIGAR: {}\nREF:   {}\nQUERY: {}",
                            i,
                            seg_idx,
                            name,
                            chrom_name,
                            segment.ref_start,
                            segment.ref_end,
                            strand,
                            e,
                            segment.alignment.cigar_string(),
                            String::from_utf8_lossy(&ref_slice),
                            String::from_utf8_lossy(query_seq)
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

                let mut tags = Vec::new();

                // Build SA tag: list all OTHER segments in this group.
                let sa_value: String = sa_entries
                    .iter()
                    .enumerate()
                    .filter(|&(k, _)| k != seg_idx)
                    .map(|(_, entry)| entry.as_str())
                    .collect::<Vec<_>>()
                    .join(";");

                tags.push((
                    Tag::try_from(*b"SA").unwrap(),
                    Value::from(sa_value.as_str()),
                ));

                tags.push((Tag::try_from(*b"XG").unwrap(), Value::from(num_segs as i32)));

                if seg_idx > 0 {
                    tags.push((
                        Tag::try_from(*b"XP").unwrap(),
                        Value::from(sa_entries[seg_idx - 1].as_str()),
                    ));
                }

                if seg_idx < num_segs - 1 {
                    tags.push((
                        Tag::try_from(*b"XN").unwrap(),
                        Value::from(sa_entries[seg_idx + 1].as_str()),
                    ));
                }

                let xo_value = xo_entries[seg_idx].join(";");
                if !xo_value.is_empty() {
                    tags.push((
                        Tag::try_from(*b"XO").unwrap(),
                        Value::from(xo_value.as_str()),
                    ));
                }

                let data: Data = if segments.len() > 1 {
                    tags.into_iter().collect()
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
                segment_written = true;
            }

            if i > 2 {
                break;
            }
        }

        if !segment_written {
            let record = build_unmapped_record(name, query, quality);
            self.writer.write_record(&record).expect("write failed");
            return Ok(());
        }

        let elapsed = start.elapsed().as_secs_f64();
        align_time_recorder().record(elapsed);

        Ok(())
    }

    fn finish(self) -> std::io::Result<()> {
        self.writer.finish()?;
        Ok(())
    }
}

fn align_time_recorder() -> &'static SimpleSummaryRecorder {
    static RECORDER: OnceLock<&'static SimpleSummaryRecorder> = OnceLock::new();
    RECORDER.get_or_init(|| SimpleSummaryRecorder::new_registered("read_time"))
}

// Assemble segments: each segment is a maximal run of colinear seeds.
// A None gap (or end of group) terminates the current segment.
// X-drop extensions are computed in a second pass and the bounds updated
// in place, so the struct is always self-consistent for validation and
// SAM emission without referencing the original seeds.
#[derive(Debug)]
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

    // Given a genome ref interval [target_ref_start, target_ref_end) that overlaps
    // this segment, return the forward-read interval [read_start, read_end) that
    // the CIGAR maps to that ref region.
    //
    // The internal CIGAR is always forward-query vs RC-ref for reverse segments:
    // CIGAR offset 0 corresponds to fwd_read_start and genome ref_end respectively.
    // For forward segments, CIGAR offset 0 = fwd_read_start and genome ref_start.
    //
    // Returns None if the overlap falls entirely in a deletion/skip (no read bases).
    fn read_range_for_ref_overlap(
        &self,
        target_ref_start: usize,
        target_ref_end: usize,
    ) -> Option<(usize, usize)> {
        // Express the target as an offset range within the CIGAR's internal ref axis.
        // For forward:  cigar_ref_offset = genome_pos - ref_start
        // For reverse:  cigar_ref_offset = ref_end - genome_pos  (RC-ref is mirrored)
        let (cigar_target_start, cigar_target_end) = if self.is_reverse {
            // Clamp to the segment's genome ref range first.
            let gs = target_ref_start.max(self.ref_start);
            let ge = target_ref_end.min(self.ref_end);
            if gs >= ge {
                return None;
            }
            // Mirrored: larger genome pos → smaller cigar ref offset.
            (self.ref_end - ge, self.ref_end - gs)
        } else {
            let gs = target_ref_start.max(self.ref_start);
            let ge = target_ref_end.min(self.ref_end);
            if gs >= ge {
                return None;
            }
            (gs - self.ref_start, ge - self.ref_start)
        };

        let mut cigar_ref_pos: usize = 0; // ref bases consumed so far
        let mut read_pos: usize = self.fwd_read_start; // read position in forward coords

        let mut result_start: Option<usize> = None;
        let mut result_end: usize = self.fwd_read_start;

        for &op in &self.alignment.cigar {
            if cigar_ref_pos >= cigar_target_end {
                break;
            }
            let n = op.len();
            let consumes_ref = op.kind().consumes_reference();
            let consumes_read = op.kind().consumes_read();

            if consumes_ref {
                let op_ref_start = cigar_ref_pos;
                let op_ref_end = cigar_ref_pos + n;

                // Clamp the op to the target window.
                let clip_start = op_ref_start.max(cigar_target_start);
                let clip_end = op_ref_end.min(cigar_target_end);

                if clip_start < clip_end && consumes_read {
                    // How many read bases before the clip start within this op?
                    let before = clip_start - op_ref_start;
                    let within = clip_end - clip_start;
                    let seg_start = read_pos + before;
                    let seg_end = seg_start + within;
                    if result_start.is_none() {
                        result_start = Some(seg_start);
                    }
                    result_end = seg_end;
                } else if clip_start < clip_end && !consumes_read {
                    // Deletion/skip covering part of the target — read_pos unchanged.
                    // result boundaries not extended, but don't return None yet.
                    if result_start.is_none() {
                        // The overlap starts in a deletion; advance to note we entered it.
                        result_start = Some(read_pos);
                        result_end = read_pos;
                    }
                }
            }

            if consumes_ref {
                cigar_ref_pos += n;
            }
            if consumes_read {
                read_pos += n;
            }
        }

        result_start.map(|s| (s, result_end))
    }

    // Walk the CIGAR over the genome ref interval [target_ref_start, target_ref_end)
    // and return the fraction of identity within that window:
}

// For each pair of adjacent segments (in fwd_read order) that overlap on the
// reference, check whether one segment aligns significantly worse in the
// overlapping ref region than the other.  When the divergence ratio
// (worse / better) exceeds `divergence_ratio_threshold`, the poorly-aligning
// end is trimmed to the overlap boundary: the trimmed query bases plus any
// inter-segment read gap are replaced by a single INS op, and the two
// segments are merged into one.
//
// Only adjacent pairs in the vector are considered — non-adjacent ref overlaps
// are left to the existing XO-tag reporting machinery.
fn merge_overlapping_segments(
    segments: &mut Vec<Segment>,
    divergence_ratio_threshold: f64,
    min_forced_overlap: usize,
    read_len: usize,
) {
    // Build the result into a fresh vector.  For each incoming segment we
    // attempt to merge it with the last segment already in `out`.  If they
    // merge, the last element of `out` is updated in place and we continue
    // (the merged segment may again be mergeable with the next incoming one).
    // This is O(n) in the number of segments.

    let mut out: Vec<Segment> = Vec::with_capacity(segments.len());

    for incoming in segments.drain(..) {
        if let Some(tail) = out.pop() {
            match try_merge(
                tail,
                incoming,
                divergence_ratio_threshold,
                min_forced_overlap,
                read_len,
            ) {
                Ok(merged) => out.push(merged),
                Err((a, b)) => {
                    out.push(a);
                    out.push(b);
                }
            }
        } else {
            out.push(incoming);
        }
    }

    *segments = out;
}

// Attempt to merge two segments that overlap on the reference.
//
// Returns Ok(merged) if the merge happened, or Err((prev, next)) if the pair
// does not qualify (different chrom/strand, no simple overlap, similar quality,
// or inconsistent CIGAR spans).
//
// For reverse-strand segments the problem is translated to an equivalent
// forward-strand problem (CIGAR reversed, roles of prev/next swapped), merged
// using the forward-only path, then translated back.
fn try_merge(
    prev: Segment,
    next: Segment,
    divergence_ratio_threshold: f64,
    min_forced_overlap: usize,
    _read_len: usize,
) -> Result<Segment, (Segment, Segment)> {
    if prev.chrom_id != next.chrom_id || prev.is_reverse != next.is_reverse {
        return Err((prev, next));
    }

    if prev.is_reverse {
        try_merge_rev(prev, next, divergence_ratio_threshold, min_forced_overlap)
    } else {
        try_merge_fwd(prev, next, divergence_ratio_threshold, min_forced_overlap)
    }
}

// Reverse-strand merge. Both segments must have is_reverse == true.
//
// Segments are ordered by fwd_read_start. On reverse strand, lower fwd_read_start
// means higher ref position, so prev covers a HIGHER ref range than next:
//
//   ref:  [next.ref_start .. next.ref_end)
//                        [prev.ref_start .. prev.ref_end)
//   overlap:             [prev.ref_start .. next.ref_end)
//
// This is exactly symmetric with try_merge_fwd, just with the ref roles of prev
// and next exchanged. The CIGAR for the merged segment (high→low ref) is:
//
//   [prev_right | INS? | winner | INS? | next_tail]
//
// where prev_right is prev's non-overlapping high portion (CIGAR start),
// next_tail is next's non-overlapping low portion (CIGAR end), and
// merged ref = [next.ref_start .. prev.ref_end).
fn try_merge_rev(
    prev: Segment,
    next: Segment,
    divergence_ratio_threshold: f64,
    min_forced_overlap: usize,
) -> Result<Segment, (Segment, Segment)> {
    // Simple overlap geometry (no containment):
    // next.ref_start < prev.ref_start < next.ref_end <= prev.ref_end
    // Equality prev.ref_end == next.ref_end is allowed (prev is a high-end suffix of next).
    if prev.ref_start <= next.ref_start
        || prev.ref_end < next.ref_end
        || prev.ref_start >= next.ref_end
    {
        return Err((prev, next));
    }
    let ref_overlap_len = next.ref_end - prev.ref_start;

    let prev_cigar_ref: usize = prev
        .alignment
        .cigar
        .iter()
        .filter(|op| op.kind().consumes_reference())
        .map(|op| op.len())
        .sum();
    let next_cigar_ref: usize = next
        .alignment
        .cigar
        .iter()
        .filter(|op| op.kind().consumes_reference())
        .map(|op| op.len())
        .sum();
    if prev_cigar_ref != prev.ref_end - prev.ref_start
        || next_cigar_ref != next.ref_end - next.ref_start
    {
        log::debug!(
            "try_merge_rev: skipping — CIGAR ref span mismatch: \
             prev cigar_ref={prev_cigar_ref} ref_span={}, \
             next cigar_ref={next_cigar_ref} ref_span={}",
            prev.ref_end - prev.ref_start,
            next.ref_end - next.ref_start,
        );
        return Err((prev, next));
    }

    // Split each alignment at the overlap boundary.
    // CIGAR pos 0 = ref_end (high genome end) for reverse strand.
    //   prev: split at (prev_cigar_ref - ref_overlap_len) → (prev_right | prev_overlap)
    //   next: split at ref_overlap_len                    → (next_overlap | next_tail)
    let (prev_right, prev_overlap) = prev
        .alignment
        .split_at_ref_pos(prev_cigar_ref - ref_overlap_len);
    let (next_overlap, next_tail) = next.alignment.split_at_ref_pos(ref_overlap_len);

    let prev_identity = prev_overlap.identity();
    let next_identity = next_overlap.identity();

    let forced = min_forced_overlap > 0 && ref_overlap_len < min_forced_overlap;
    let prev_is_worse = forced && prev_identity <= next_identity
        || prev_identity < next_identity * divergence_ratio_threshold;
    let next_is_worse = forced && next_identity < prev_identity
        || next_identity < prev_identity * divergence_ratio_threshold;

    if !prev_is_worse && !next_is_worse {
        return Err((prev, next));
    }
    overlap_size_recorder().record(ref_overlap_len);

    // Query bases consumed by the losing overlap piece become an INS.
    let prev_overlap_query: usize = prev_overlap
        .cigar
        .iter()
        .filter(|op| op.kind().consumes_read())
        .map(|op| op.len())
        .sum();
    let next_overlap_query: usize = next_overlap
        .cigar
        .iter()
        .filter(|op| op.kind().consumes_read())
        .map(|op| op.len())
        .sum();

    assert!(
        next.fwd_read_start >= prev.fwd_read_end,
        "try_merge_rev: unexpected read overlap: prev.fwd_read_end={} > next.fwd_read_start={}",
        prev.fwd_read_end,
        next.fwd_read_start,
    );
    let inter_segment_read_gap = next.fwd_read_start - prev.fwd_read_end;

    // Symmetric with try_merge_fwd: [prev_right | INS? | winner | INS? | next_tail]
    let (winner, ins_near_prev_right, ins_near_next_tail) = if prev_is_worse {
        (next_overlap, prev_overlap_query + inter_segment_read_gap, 0)
    } else {
        (prev_overlap, 0, next_overlap_query + inter_segment_read_gap)
    };

    let ins = |n: usize| -> Option<Alignment> {
        if n > 0 {
            Some(Alignment::from(vec![Op::new(Kind::Insertion, n)]))
        } else {
            None
        }
    };

    let parts: Vec<Alignment> = [
        Some(prev_right),
        ins(ins_near_prev_right),
        Some(winner),
        ins(ins_near_next_tail),
        Some(next_tail),
    ]
    .into_iter()
    .flatten()
    .collect();

    let mut merged = prev;
    merged.alignment = Alignment::concat(&parts);
    merged.ref_start = next.ref_start;
    merged.fwd_read_end = next.fwd_read_end;

    log::debug!(
        "  merged: ref=[{},{}), read=[{},{})",
        merged.ref_start,
        merged.ref_end,
        merged.fwd_read_start,
        merged.fwd_read_end,
    );

    Ok(merged)
}

// Forward-strand-only merge logic. Both segments must have is_reverse == false
// and satisfy the simple overlap geometry (prev.ref_start < next.ref_start <
// prev.ref_end < next.ref_end).
fn try_merge_fwd(
    prev: Segment,
    next: Segment,
    divergence_ratio_threshold: f64,
    min_forced_overlap: usize,
) -> Result<Segment, (Segment, Segment)> {
    // Reject non-overlapping and non-simple-overlap (containment) geometries.
    // Requires: prev.ref_start < next.ref_start < prev.ref_end <= next.ref_end
    // The equality next.ref_end == prev.ref_end is allowed: next_right is empty
    // and next_overlap covers all of next (next is a suffix of prev on the ref).
    if next.ref_start <= prev.ref_start
        || next.ref_end < prev.ref_end
        || next.ref_start >= prev.ref_end
    {
        return Err((prev, next));
    }
    let ref_overlap_len = prev.ref_end - next.ref_start;

    // Guard: CIGAR must exactly span the declared ref range for both segments.
    // A mismatch means a previous merge left a segment in an inconsistent state.
    let prev_cigar_ref: usize = prev
        .alignment
        .cigar
        .iter()
        .filter(|op| op.kind().consumes_reference())
        .map(|op| op.len())
        .sum();
    let next_cigar_ref: usize = next
        .alignment
        .cigar
        .iter()
        .filter(|op| op.kind().consumes_reference())
        .map(|op| op.len())
        .sum();
    if prev_cigar_ref != prev.ref_end - prev.ref_start
        || next_cigar_ref != next.ref_end - next.ref_start
    {
        log::debug!(
            "try_merge_fwd: skipping — CIGAR ref span mismatch: \
             prev cigar_ref={prev_cigar_ref} ref_span={}, \
             next cigar_ref={next_cigar_ref} ref_span={}",
            prev.ref_end - prev.ref_start,
            next.ref_end - next.ref_start,
        );
        return Err((prev, next));
    }

    // Split each alignment at the overlap boundary.
    // CIGAR pos 0 = ref_start (low genome end) for forward strand.
    //   prev: split at (prev_cigar_ref - ref_overlap_len) → (prev_left | prev_overlap)
    //   next: split at ref_overlap_len                    → (next_overlap | next_right)
    let (prev_left, prev_overlap) = prev
        .alignment
        .split_at_ref_pos(prev_cigar_ref - ref_overlap_len);
    let (next_overlap, next_right) = next.alignment.split_at_ref_pos(ref_overlap_len);

    let prev_identity = prev_overlap.identity();
    let next_identity = next_overlap.identity();

    let forced = min_forced_overlap > 0 && ref_overlap_len < min_forced_overlap;
    let prev_is_worse = forced && prev_identity <= next_identity
        || prev_identity < next_identity * divergence_ratio_threshold;
    let next_is_worse = forced && next_identity < prev_identity
        || next_identity < prev_identity * divergence_ratio_threshold;

    if !prev_is_worse && !next_is_worse {
        return Err((prev, next));
    }
    overlap_size_recorder().record(ref_overlap_len);

    // Query bases consumed by the losing overlap piece become an INS.
    let prev_overlap_query: usize = prev_overlap
        .cigar
        .iter()
        .filter(|op| op.kind().consumes_read())
        .map(|op| op.len())
        .sum();
    let next_overlap_query: usize = next_overlap
        .cigar
        .iter()
        .filter(|op| op.kind().consumes_read())
        .map(|op| op.len())
        .sum();

    // Any read bases between the two segments are folded into the INS.
    assert!(
        next.fwd_read_start >= prev.fwd_read_end,
        "try_merge: unexpected read overlap: prev.fwd_read_end={} > next.fwd_read_start={}",
        prev.fwd_read_end,
        next.fwd_read_start,
    );
    let inter_segment_read_gap = next.fwd_read_start - prev.fwd_read_end;

    let (winner, ins_near_prev_left, ins_near_next_right) = if prev_is_worse {
        (next_overlap, prev_overlap_query + inter_segment_read_gap, 0)
    } else {
        (prev_overlap, 0, next_overlap_query + inter_segment_read_gap)
    };

    let ins = |n: usize| -> Option<Alignment> {
        if n > 0 {
            Some(Alignment::from(vec![Op::new(Kind::Insertion, n)]))
        } else {
            None
        }
    };

    let parts: Vec<Alignment> = [
        Some(prev_left),
        ins(ins_near_prev_left),
        Some(winner),
        ins(ins_near_next_right),
        Some(next_right),
    ]
    .into_iter()
    .flatten()
    .collect();

    let mut merged = prev;
    merged.alignment = Alignment::concat(&parts);
    merged.ref_end = next.ref_end;
    merged.fwd_read_end = next.fwd_read_end;

    log::debug!(
        "  merged: ref=[{},{}), read=[{},{})",
        merged.ref_start,
        merged.ref_end,
        merged.fwd_read_start,
        merged.fwd_read_end,
    );

    Ok(merged)
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

fn seed_length_recorder() -> &'static HistogramRecorder {
    static RECORDER: OnceLock<&'static HistogramRecorder> = OnceLock::new();
    RECORDER.get_or_init(|| HistogramRecorder::new_registered("seed_length"))
}

fn segment_count_recorder() -> &'static HistogramRecorder {
    static RECORDER: OnceLock<&'static HistogramRecorder> = OnceLock::new();
    RECORDER.get_or_init(|| HistogramRecorder::new_registered("segment_count"))
}

fn overlap_size_recorder() -> &'static HistogramRecorder {
    static RECORDER: OnceLock<&'static HistogramRecorder> = OnceLock::new();
    RECORDER.get_or_init(|| HistogramRecorder::new_registered("segment_merge_size"))
}

#[cfg(test)]
mod segment_tests {
    use super::*;

    fn make_aln(cigar_str: &str) -> Alignment {
        use crate::align::{Kind, Op};
        let mut cigar = Vec::new();
        let mut count = 0usize;
        for c in cigar_str.chars() {
            if let Some(d) = c.to_digit(10) {
                count = count * 10 + d as usize;
            } else {
                let kind = match c {
                    '=' => Kind::SequenceMatch,
                    'X' => Kind::SequenceMismatch,
                    'I' => Kind::Insertion,
                    'D' => Kind::Deletion,
                    _ => panic!("unknown cigar op {c}"),
                };
                cigar.push(Op::new(kind, count));
                count = 0;
            }
        }
        Alignment {
            divergence: parallax::scores::DivergenceScore::ZERO,
            cigar,
        }
    }

    fn fwd_seg(
        ref_start: usize,
        ref_end: usize,
        fwd_read_start: usize,
        fwd_read_end: usize,
        cigar: &str,
    ) -> Segment {
        Segment {
            alignment: make_aln(cigar),
            chrom_id: 0,
            is_reverse: false,
            fwd_read_start,
            fwd_read_end,
            ref_start,
            ref_end,
        }
    }

    fn rev_seg(
        ref_start: usize,
        ref_end: usize,
        fwd_read_start: usize,
        fwd_read_end: usize,
        cigar: &str,
    ) -> Segment {
        Segment {
            alignment: make_aln(cigar),
            chrom_id: 0,
            is_reverse: true,
            fwd_read_start,
            fwd_read_end,
            ref_start,
            ref_end,
        }
    }

    // ── read_range_for_ref_overlap ─────────────────────────────────────────

    #[test]
    fn read_range_fwd_exact() {
        // Forward segment: ref [100,110), cigar 10=, read [0,10)
        let seg = fwd_seg(100, 110, 0, 10, "10=");
        // Query overlap for ref [102,107)
        let r = seg.read_range_for_ref_overlap(102, 107);
        assert_eq!(r, Some((2, 7)));
    }

    #[test]
    fn read_range_fwd_with_deletion() {
        // ref [100,115), cigar 5=5D5=, read [0,10)
        // Forward: ref offset 0-4 = match, 5-9 = deletion (no read), 10-14 = match
        let seg = fwd_seg(100, 115, 0, 10, "5=5D5=");
        // Overlap in the deletion [105,110) — read range is a zero-width point
        let r = seg.read_range_for_ref_overlap(105, 110);
        assert_eq!(r, Some((5, 5))); // both sides of the gap map to same read pos
        // Overlap spanning deletion into second match [107,113)
        let r2 = seg.read_range_for_ref_overlap(107, 113);
        assert_eq!(r2, Some((5, 8))); // 3 bases from second match block
    }

    #[test]
    fn read_range_fwd_no_overlap() {
        let seg = fwd_seg(100, 110, 0, 10, "10=");
        assert_eq!(seg.read_range_for_ref_overlap(200, 210), None);
    }

    #[test]
    fn read_range_rev_exact() {
        // Reverse segment: ref [100,110), cigar 10=, fwd_read [0,10)
        // CIGAR offset 0 = ref_end-1=109, offset 9 = 100.
        // Ref [103,108) → cigar offsets [2,7) → read [2,7)
        let seg = rev_seg(100, 110, 0, 10, "10=");
        let r = seg.read_range_for_ref_overlap(103, 108);
        assert_eq!(r, Some((2, 7)));
    }

    #[test]
    fn read_range_rev_with_deletion() {
        // Reverse segment: ref [100,115), cigar 5=5D5= (same CIGAR, reverse orientation)
        // CIGAR walks 109→105 (5=), 104→100 (5D), 99→95 (5=).
        // Wait — ref_end=115, so offset 0=114, 1=113,... 4=110 (5=), 5-9=del(110-106?).
        // Actually offset k = ref_end - 1 - k in genome.
        // 5= covers cigar offsets 0-4 → genome 114..110 (exclusive)
        // 5D covers cigar offsets 5-9 → genome 109..105 (no read)
        // 5= covers cigar offsets 10-14 → genome 104..100
        //
        // Overlap ref [100,105) → cigar offsets [ref_end-105, ref_end-100) = [10,15)
        // → falls in the second 5= block → read [5,10)
        let seg = rev_seg(100, 115, 0, 10, "5=5D5=");
        let r = seg.read_range_for_ref_overlap(100, 105);
        assert_eq!(r, Some((5, 10)));
        // Overlap in the deletion [105,110) → cigar offsets [5,10)
        let r2 = seg.read_range_for_ref_overlap(105, 110);
        assert_eq!(r2, Some((5, 5)));
    }

    // ── try_merge ──────────────────────────────────────────────────────────

    #[test]
    fn try_merge_no_overlap_returns_err() {
        let prev = fwd_seg(100, 200, 0, 100, "100=");
        let next = fwd_seg(200, 300, 100, 200, "100=");
        assert!(try_merge(prev, next, 0.5, 0, 100000).is_err());
    }

    #[test]
    fn try_merge_similar_quality_returns_err() {
        // Both sides have equal identity in the overlap — below threshold
        let prev = fwd_seg(100, 210, 0, 110, "110=");
        let next = fwd_seg(200, 300, 100, 200, "100=");
        // ref overlap [200,210) = 10bp; both have identity 1.0
        // 1.0 < 1.0 * 0.5 is false — neither qualifies
        assert!(try_merge(prev, next, 0.5, 0, 100000).is_err());
    }

    #[test]
    fn try_merge_fwd_prev_worse_merges() {
        // prev: ref [100,210), cigar 100=10X (100 matches then 10 mismatches)
        // next: ref [200,300), cigar 10=90=  (first 10 = perfect over the overlap)
        // overlap [200,210): prev identity = 0/10 = 0.0, next identity = 1.0
        // 0.0 < 1.0 * 0.5 → prev_is_worse → trim prev suffix by 10
        let prev = fwd_seg(100, 210, 0, 110, "100=10X");
        let next = fwd_seg(200, 300, 110, 210, "10=90=");
        let merged = try_merge(prev, next, 0.5, 0, 100000).expect("should merge");
        assert_eq!(merged.ref_start, 100);
        assert_eq!(merged.ref_end, 300);
        // merged: 100= + 10=90= → 200= (concat merges adjacent same-kind ops)
        let ref_consumed: usize = merged
            .alignment
            .cigar
            .iter()
            .filter(|op| op.kind().consumes_reference())
            .map(|op| op.len())
            .sum();
        assert_eq!(ref_consumed, 200);
    }

    #[test]
    fn try_merge_fwd_next_worse_merges() {
        // prev: ref [100,210), cigar 100=10= (perfect)
        // next: ref [200,300), cigar 10X90= (10 mismatches over overlap then 90 matches)
        // overlap [200,210): prev identity = 1.0, next identity = 0/10 = 0.0
        // 0.0 < 1.0 * 0.5 → next_is_worse → trim next prefix by 10
        let prev = fwd_seg(100, 210, 0, 110, "100=10=");
        let next = fwd_seg(200, 300, 110, 210, "10X90=");
        let merged = try_merge(prev, next, 0.5, 0, 100000).expect("should merge");
        assert_eq!(merged.ref_start, 100);
        assert_eq!(merged.ref_end, 300);
        let ref_consumed: usize = merged
            .alignment
            .cigar
            .iter()
            .filter(|op| op.kind().consumes_reference())
            .map(|op| op.len())
            .sum();
        assert_eq!(ref_consumed, 200);
    }

    #[test]
    fn try_merge_guard_overlap_exceeds_cigar() {
        // Use a case where next.ref_end > prev.ref_end but cigar is short.
        let prev = fwd_seg(100, 210, 0, 110, "110=");
        let next = fwd_seg(130, 250, 80, 200, "50="); // cigar_ref=50 < overlap=80
        // overlap = [130,210) = 80, but next cigar_ref = 50 → guard bails
        assert!(try_merge(prev, next, 0.5, 0, 100000).is_err());
    }

    #[test]
    fn try_merge_rev_prev_worse_merges() {
        // Reverse strand: segments ordered by fwd_read_start, so prev covers HIGHER ref.
        // prev: ref [200,300), fwd_read [0,100), cigar 90=10X — mismatches at low ref end (overlap)
        //   CIGAR pos 0 = genome 299; mismatches cover genome [200,210)
        //   cigar ref span = 100 ✓, cigar read span = 100 ✓
        // next: ref [100,210), fwd_read [100,210), cigar 10=100= — all matches
        //   CIGAR pos 0 = genome 209; cigar ref span = 110 ✓
        // overlap [200,210) = 10bp; prev_identity=0.0, next_identity=1.0 → prev_is_worse
        //
        // Assembly (high→low ref): [prev_right(90=) | INS(10) | winner(10=) | next_tail(100=)]
        // concat merges adjacent 10= and 100= → [90= | 10I | 110=]
        // merged ref = [next.ref_start .. prev.ref_end) = [100, 300), ref_consumed = 200
        let prev = rev_seg(200, 300, 0, 100, "90=10X");
        let next = rev_seg(100, 210, 100, 210, "10=100=");
        let merged = try_merge(prev, next, 0.5, 0, 100000).expect("should merge");
        assert_eq!(merged.ref_start, 100);
        assert_eq!(merged.ref_end, 300);
        assert_eq!(merged.fwd_read_start, 0);
        assert_eq!(merged.fwd_read_end, 210);
        let ref_consumed: usize = merged
            .alignment
            .cigar
            .iter()
            .filter(|op| op.kind().consumes_reference())
            .map(|op| op.len())
            .sum();
        assert_eq!(ref_consumed, 200);
        // INS covers prev_overlap query bytes (10) + inter-segment gap (0)
        assert_eq!(merged.alignment.cigar[1], Op::new(Kind::Insertion, 10));
    }

    #[test]
    fn try_merge_different_chrom_returns_err() {
        let prev = fwd_seg(100, 210, 0, 110, "110=");
        let mut next = fwd_seg(200, 300, 110, 210, "100=");
        next.chrom_id = 1;
        assert!(try_merge(prev, next, 0.5, 0, 100000).is_err());
    }

    #[test]
    fn try_merge_mixed_strand_returns_err() {
        let prev = fwd_seg(100, 210, 0, 110, "110=");
        let next = rev_seg(200, 300, 110, 210, "100=");
        assert!(try_merge(prev, next, 0.5, 0, 100000).is_err());
    }

    // ── regression tests from real failing merges ──────────────────────────
    //
    // These are built directly from the debug log output of try_merge at the
    // point of the validation failures seen in production.  The pre-trim
    // segment state is used (the first ref=[] in the log line).

    // chr12 case 1 (seg 11): rev strand, prev_is_worse, overlap=190 > prev.cigar_ref=48
    // Guard must bail rather than attempting a trim that would exhaust the CIGAR.
    #[test]
    fn regression_chr12_prev_worse_overlap_exceeds_cigar() {
        let prev = rev_seg(2255803, 2256041, 16829, 16947, "2D46=");
        // cigar_ref = 48, ref_overlap = 190 → guard should fire
        let next = rev_seg(
            2255851,
            2256072,
            16948,
            17139,
            "25=1X13=1X37=1X21=1X1=1X27=1X1=1X5=30D22=1X2=1X28=",
        );
        assert!(
            try_merge(prev, next, 0.5, 0, 100000).is_err(),
            "must not merge when overlap ({}) exceeds prev cigar_ref ({})",
            190,
            48
        );
    }

    // chr12 case 2 (seg 13): rev strand, next_is_worse, overlap=179 > next.cigar_ref=43
    #[test]
    fn regression_chr12_next_worse_overlap_exceeds_cigar() {
        let prev = rev_seg(2255819, 2256029, 17321, 17531, "107=1X16=1X85=");
        // cigar_ref = 210, next.cigar_ref = 43, overlap = 179 → guard fires on next
        let next = rev_seg(2255850, 2256072, 17549, 17681, "25=1X17=");
        assert!(
            try_merge(prev, next, 0.5, 0, 100000).is_err(),
            "must not merge when overlap ({}) exceeds next cigar_ref ({})",
            179,
            43
        );
    }

    // chr13 case: rev strand, next_is_worse, overlap=96 > next.cigar_ref=91
    #[test]
    fn regression_chr13_next_worse_overlap_exceeds_cigar() {
        let prev = rev_seg(25159116, 25159221, 2823, 2933, "57=5I48=");
        // cigar_ref = 105, next.cigar_ref = 91, overlap = 96 → guard fires on next
        let next = rev_seg(25159125, 25159312, 2965, 3037, "23=1X4=63D");
        assert!(
            try_merge(prev, next, 0.5, 0, 100000).is_err(),
            "must not merge when overlap ({}) exceeds next cigar_ref ({})",
            96,
            91
        );
    }

    // chr6 case: rev strand, next_is_worse, overlap=146, next.cigar_ref=290 (guard
    // does NOT fire on size alone), but next.ref_end - next.cigar_ref = 167169535
    // which differs from next.ref_start = 167169389 — the CIGAR doesn't anchor at
    // ref_start, so the merge would produce an invalid joined CIGAR.
    // The correct behaviour is to bail; this test documents the expected fix.
    #[test]
    fn regression_chr6_next_cigar_does_not_span_ref_start() {
        let prev = rev_seg(
            167169085,
            167169535,
            12250,
            12389,
            "42=1D15=24D3=1X1=1X22=259D14=1X1=27D38=",
        );
        // next.ref_end - next.cigar_ref = 167169825 - 290 = 167169535 ≠ next.ref_start 167169389
        let next = rev_seg(
            167169389,
            167169825,
            12401,
            12522,
            "32=1X2=2I8=1X7=1X3=1X34=200D",
        );
        assert!(
            try_merge(prev, next, 0.5, 0, 100000).is_err(),
            "must not merge when next CIGAR does not span next.ref_start \
             (cigar anchors at {} not {})",
            167169825u32.saturating_sub(290),
            167169389
        );
    }

    // chr6 new cases (reads 3559, 3560, 3561): rev strand, next_is_worse, overlap=411.
    // next: ref=[167169124,167169775) span=651, cigar=31=2X1=1D5=175D25= cigar_ref=240.
    // 651 ≠ 240 → new span-mismatch guard must fire.
    #[test]
    fn regression_chr6_next_cigar_span_mismatch_overlap_411() {
        let prev = rev_seg(167169085, 167169535, 9772, 9912, "58=283D5=1X36=1X1=27D38=");
        // next ref span 651 but cigar_ref=240 — span mismatch guard must fire
        let next = rev_seg(167169124, 167169775, 9918, 10019, "31=2X1=1D5=175D25=");
        assert!(
            try_merge(prev, next, 0.5, 0, 100000).is_err(),
            "must not merge when next cigar_ref (240) != next ref_span (651)"
        );
    }

    // Full merge_overlapping_segments regression for read SRR29147690.2771 on chr12 (-).
    // Segments taken verbatim from the "before overlap merge" debug log.
    // The invariant checked: every output segment has cigar_ref == ref_end - ref_start.
    #[test]
    fn regression_chr12_2771_full_merge() {
        let mut segs = vec![
            rev_seg(
                2255820,
                2270751,
                0,
                15040,
                "2753=1D10=1D252=1X2004=1D2245=1D2161=1D27=1D357=1D279=1I178=2D1518=1D1314=1D326=1X577=1X280=1D365=180I43=1X43=1X29=1X29=1X7=60D22=1X2=1X29=",
            ),
            rev_seg(
                2255802,
                2256090,
                15040,
                15298,
                "43=1X13=1X29=1X32=1X17=1X8=1X1=1X5=1X11=1X9=1X29=1X2=1X5=1X12=30D29=",
            ),
            rev_seg(
                2255806,
                2256072,
                15328,
                15474,
                "25=1X3=2X38=1X29=1X21=120D25=",
            ),
            rev_seg(2255882, 2256046, 15474, 15608, "81=1X1=1X29=1X30D20="),
            rev_seg(2255802, 2256032, 15608, 15778, "37=1X29=1X60D102="),
            rev_seg(
                2255851,
                2256046,
                15804,
                16028,
                "23=1X22=1X6=1X40=1X16=1X6=1D22=1X24=30I29=",
            ),
            rev_seg(2255853, 2256030, 16029, 16177, "86=29D62="),
            rev_seg(2255802, 2256032, 16178, 16348, "110=1X30=60D29="),
            rev_seg(
                2255803,
                2256032,
                16358,
                16527,
                "32=1X26=1X1=1X35=1X11=1X9=1X21=60D28=",
            ),
            rev_seg(2255831, 2256072, 16528, 16679, "25=1X34=1X76=90D12=1X1="),
            rev_seg(2255831, 2256041, 16679, 16829, "76=1X29=30D22=30D20=1X1="),
            rev_seg(2255803, 2256041, 16829, 16947, "46=1X25=120D46="),
            rev_seg(
                2255851,
                2256072,
                16948,
                17139,
                "25=1X13=1X37=1X21=1X1=1X27=1X1=1X5=30D22=1X2=1X28=",
            ),
            rev_seg(
                2255803,
                2256030,
                17140,
                17307,
                "27=1X1=1X27=1X1=1X27=1X1=1X5=30D44=30D28=",
            ),
            rev_seg(2255819, 2256029, 17321, 17531, "107=1X16=1X85="),
            rev_seg(2255850, 2256072, 17549, 17681, "25=1X51=1X90D54="),
            rev_seg(2255803, 2256030, 17681, 17849, "35=1X2=29D29=30D101="),
            rev_seg(2255850, 2256072, 17850, 17982, "25=1X51=1X90D54="),
            rev_seg(2255803, 2256030, 17982, 18179, "57=1X30D139="),
            rev_seg(2255790, 2256030, 18192, 18342, "95=90D12=1X42="),
            rev_seg(
                2255806,
                2256084,
                18378,
                18566,
                "21=1X7=1X7=30D22=1X72=1X28=1X60D26=",
            ),
            rev_seg(2255874, 2256046, 18566, 18708, "34=1X46=1X13=30D47="),
            rev_seg(2255853, 2256046, 18716, 18879, "22=1X80=1X7=30D1=1X50="),
            rev_seg(2255882, 2256032, 18880, 19030, "29=1X1=1X78=1X8=1X8=1X21="),
            rev_seg(2255866, 2256032, 19030, 19166, "59=1X21=30D55="),
            rev_seg(2255853, 2256046, 19166, 19329, "22=1X88=30D1=1X50="),
            rev_seg(2255866, 2256032, 19330, 19526, "29=1X1=1X78=1X30I55="),
            rev_seg(
                2255850,
                2256046,
                19526,
                19722,
                "22=1X50=1X1=1X35=1X12=1X71=",
            ),
            rev_seg(2255820, 2256024, 19728, 19902, "29=1X29=1X42=1X58=1X30D12="),
            rev_seg(2255803, 2256090, 19902, 20098, "43=1X3=2X38=1X57=91D51="),
            rev_seg(2253889, 2256072, 20099, 22279, "786=1D182=1D1164=1D48="),
        ];

        merge_overlapping_segments(&mut segs, 0.5, 0, 22279);

        for (i, seg) in segs.iter().enumerate() {
            let cigar_ref: usize = seg
                .alignment
                .cigar
                .iter()
                .filter(|op| op.kind().consumes_reference())
                .map(|op| op.len())
                .sum();
            let cigar_read: usize = seg
                .alignment
                .cigar
                .iter()
                .filter(|op| op.kind().consumes_read())
                .map(|op| op.len())
                .sum();
            let ref_span = seg.ref_end - seg.ref_start;
            let read_span = seg.fwd_read_end - seg.fwd_read_start;
            assert_eq!(
                cigar_ref,
                ref_span,
                "seg {i}: cigar_ref={cigar_ref} != ref_span={ref_span} \
                 (ref=[{},{}), read=[{},{}), cigar={:?})",
                seg.ref_start,
                seg.ref_end,
                seg.fwd_read_start,
                seg.fwd_read_end,
                seg.alignment.cigar,
            );
            assert_eq!(
                cigar_read,
                read_span,
                "seg {i}: cigar_read={cigar_read} != read_span={read_span} \
                 (ref=[{},{}), read=[{},{}), cigar={:?})",
                seg.ref_start,
                seg.ref_end,
                seg.fwd_read_start,
                seg.fwd_read_end,
                seg.alignment.cigar,
            );
        }
    }

    // Isolated single-pair merge matching seg 11/12 from the chr12 2771 sequence,
    // with prev/next in the correct reverse-strand order (prev = higher ref range).
    // overlap = [2255851..2256041) = 190bp; prev_is_worse → winner is next_overlap.
    // merged ref = [next.ref_start .. prev.ref_end) = [2255803..2256072) = 269bp
    // read span = 17139 - 16829 = 310
    #[test]
    fn regression_chr12_2771_seg11_seg12_cigar_order() {
        let prev = rev_seg(
            2255851,
            2256072,
            16829,
            16947,
            "25=1X13=1X37=1X21=1X1=1X27=1X1=1X5=30D22=1X2=1X28=",
        );
        let next = rev_seg(2255803, 2256041, 16948, 17139, "46=1X25=120D46=");
        let merged = try_merge(prev, next, 0.5, 0, 100000).expect("should merge");
        assert_eq!(merged.ref_start, 2255803);
        assert_eq!(merged.ref_end, 2256072);
        assert_eq!(merged.fwd_read_start, 16829);
        assert_eq!(merged.fwd_read_end, 17139);
        let cigar_ref: usize = merged
            .alignment
            .cigar
            .iter()
            .filter(|op| op.kind().consumes_reference())
            .map(|op| op.len())
            .sum();
        let cigar_read: usize = merged
            .alignment
            .cigar
            .iter()
            .filter(|op| op.kind().consumes_read())
            .map(|op| op.len())
            .sum();
        assert_eq!(
            cigar_ref, 269,
            "CIGAR ref span mismatch: {:?}",
            merged.alignment.cigar
        );
        assert_eq!(
            cigar_read, 310,
            "CIGAR read span mismatch: {:?}",
            merged.alignment.cigar
        );
    }
}
