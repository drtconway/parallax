use std::{cmp::Reverse, sync::{Arc, OnceLock}};

use crate::{
    align::Alignment,
    aligner::{Aligner, AlignerBuilder},
    reads::{
        builder::{build_record, build_unmapped_record},
        compound::{
            EdgeType, Seed, SeedCollection, TagValue, Weighted, prune_isolated_seeds,
            seed_to_record,
        },
        extended::ExtendedSeed,
        segments::{BandedSegmentScheme, SegmentConfig, extract_chains, partition_by_diagonal},
    },
    seeding::SeedCollector,
    writer::{AlignmentWriter, OutputFormat, RecordWriter},
};
use noodles::sam::alignment::{
    RecordBuf,
    record::{
        Flags,
        cigar::{Op, op::Kind},
        data::field::Tag,
    },
    record_buf::{Cigar, Data, data::field::Value},
};
use ordered_float::OrderedFloat;
use parallax::{
    config::{self, FilteringConfig, SeedingConfig},
    index::Index,
    reference::InMemoryReference,
    utils::{
        sequence::{complement, reverse_complement_into},
        telemetry::{RecorderExt, summary::SimpleSummaryRecorder},
    },
};

static SEED_WRITER: OnceLock<Arc<AlignmentWriter>> = OnceLock::new();
static CHAIN_WRITER: OnceLock<Arc<AlignmentWriter>> = OnceLock::new();

pub struct ExplanatoryAlignerBuilder<'a> {
    reference: &'a InMemoryReference,
    index: &'a dyn Index,
    writer: &'a dyn RecordWriter,
    no_secondary: bool,
}

impl<'a> ExplanatoryAlignerBuilder<'a> {
    pub fn no_secondary(mut self, no_secondary: bool) -> Self {
        self.no_secondary = no_secondary;
        self
    }
}

impl<'a> AlignerBuilder<'a> for ExplanatoryAlignerBuilder<'a> {
    type AlignerType = ExplanatoryAligner<'a>;

    fn new(
        reference: &'a InMemoryReference,
        index: &'a dyn Index,
        writer: &'a dyn RecordWriter,
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

        let make_writer = |path: &str,
                           cell: &'static OnceLock<Arc<AlignmentWriter>>|
         -> Option<Arc<AlignmentWriter>> {
            if path.is_empty() {
                return None;
            }
            let path = path.to_string();
            let reference = self.reference;
            Some(Arc::clone(cell.get_or_init(|| {
                let format = OutputFormat::from_path(std::path::Path::new(&path))
                    .unwrap_or(OutputFormat::Sam);
                let file = std::fs::File::create(&path)
                    .unwrap_or_else(|e| panic!("failed to create {path}: {e}"));
                let repo = noodles::fasta::Repository::default();
                Arc::new(
                    AlignmentWriter::builder(Box::new(file), format, repo)
                        .add_contigs(
                            reference
                                .all_chrom_info()
                                .iter()
                                .map(|c| (c.name.as_str(), c.length as u64)),
                        )
                        .build()
                        .unwrap_or_else(|e| panic!("failed to write header to {path}: {e}")),
                )
            })))
        };

        ExplanatoryAligner {
            reference: self.reference,
            index: self.index,
            writer: self.writer,
            seed_writer: make_writer(&cfg.seeding.debug_seeds_sam, &SEED_WRITER),
            chain_writer: make_writer(&cfg.seeding.debug_chains_sam, &CHAIN_WRITER),
            seeder: SeedCollector::new(),
            aligner: crate::align::DpAligner::from_config(&cfg.alignment, &cfg.block_aligner),
            _all_seeds: Vec::new(),
            no_secondary: self.no_secondary,
            seeding_cfg: cfg.seeding.clone(),
            _filtering_cfg: cfg.filtering.clone(),
        }
    }
}

pub struct ExplanatoryAligner<'a> {
    reference: &'a InMemoryReference,
    index: &'a dyn Index,
    writer: &'a dyn RecordWriter,
    seed_writer: Option<Arc<AlignmentWriter>>,
    chain_writer: Option<Arc<AlignmentWriter>>,
    seeder: SeedCollector,
    aligner: crate::align::DpAligner,
    _all_seeds: Vec<ExtendedSeed>,
    no_secondary: bool,
    seeding_cfg: SeedingConfig,
    _filtering_cfg: FilteringConfig,
}

impl<'a> Aligner<'a> for ExplanatoryAligner<'a> {
    fn align(&mut self, name: &str, query: &[u8], quality: &[u8]) -> std::io::Result<()> {
        log::info!("aligning {} ({}bp)", name, query.len());
        let start = std::time::Instant::now();

        let k = self.index.k();

        let query_len = query.len();
        let mut query_rc = Vec::with_capacity(query_len);
        reverse_complement_into(query, &mut query_rc);

        let mut fwd_seeds = vec![];
        let mut rev_seeds = vec![];

        self.seeder.gather_seeds_batched2(
            query,
            false,
            self.index,
            name,
            &self.seeding_cfg,
            &mut fwd_seeds,
        );
        log::info!("read {}: {} fwd seeds", name, fwd_seeds.len());
        self.seeder.gather_seeds_batched2(
            &query_rc,
            true,
            self.index,
            name,
            &self.seeding_cfg,
            &mut rev_seeds,
        );

        let mut all_seeds: Vec<crate::reads::compound::AtomicSeed> = fwd_seeds;
        all_seeds.extend(rev_seeds);
        all_seeds.sort_unstable_by_key(|s| (s.read_start(), s.chrom_id(), s.is_reverse()));

        log::info!("read {}: {} seeds", name, all_seeds.len());

        let collection = SeedCollection::new(k, all_seeds);
        let mut compounds = collection.compound_seeds();

        let full_count = compounds.len();
        prune_isolated_seeds(&mut compounds, k, 400, 800, 100);
        compounds.sort_unstable_by_key(|s| s.read_pos());
        let pruned_count = compounds.len();
        log::info!("pruned seeds from {} to {}", full_count, pruned_count);

        if let Some(ref seed_writer) = self.seed_writer {
            log::info!("dumping {} seeds to SAM/BAM", compounds.len());
            for (i, seed) in compounds.iter().enumerate() {
                // SEQ is always taken from the forward-strand query at the seed's
                // forward-strand coordinates. For FLAG=16 records IGV RC's the SEQ
                // to compare against the reference, and RC(query_rc[b..e]) =
                // query[read_start..read_end], so this is the correct orientation.
                let b = seed.read_start() as usize;
                let e = seed.read_end(k) as usize;
                let tags = vec![
                    (String::from("XN"), TagValue::Int(i as i64)),
                    (
                        String::from("XF"),
                        TagValue::Int(seed.multiplicity() as i64),
                    ),
                    (String::from("XW"), TagValue::Flt(seed.weight())),
                ];
                let record =
                    seed_to_record(name, k, query_len, seed, &query[b..e], &quality[b..e], tags);
                seed_writer
                    .write_record(&record)
                    .expect("seed write failed");
            }
            seed_writer.finish().expect("finish failed");
        }

        // Diagonal bands.
        let mut bands = partition_by_diagonal(&compounds, k, 800, query_len);
        let total_bands = bands.len();
        bands.retain(|b| b.coverage >= 0.33);
        bands.sort_by_key(|b| Reverse(OrderedFloat(b.coverage)));

        log::info!("Bands recovered: {}, keeping {}", total_bands, bands.len());
        log::info!("band\tchrom\tdiag\tstrand\tcount\tspan\tcoverage\tdiag_var");
        for (i, band) in bands.iter().enumerate() {
            log::info!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.1}",
                i,
                &self.reference.chrom_name(band.chrom_id as usize),
                band.central_diagonal.floor() as i64 + 1,
                if band.is_reverse { "-" } else { "+" },
                band.members.len(),
                band.ref_max.max(band.ref_min) - band.ref_min.min(band.ref_max),
                band.coverage,
                band.diagonal_variance,
            );
        }

        let mut segment_written = false;

        for (band_rank, band) in bands.iter().enumerate() {
            let band_scheme = BandedSegmentScheme::new(
                SegmentConfig {
                    read_gap_cost_per_base: 0.003,
                    max_read_gap: 600,
                    ..SegmentConfig::default_for_k(k)
                },
                band.chrom_id,
                band.is_reverse,
                band.central_diagonal,
                band.diagonal_variance,
                1.0,    // diag_lambda
                15.0,   // sv_break_penalty
                600,
            );
            let mut band_chains = extract_chains(&compounds, k, &band_scheme);
            band_chains.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());

            for (chain_num, chain) in band_chains.iter().enumerate() {
                if let Some(ref chain_writer) = self.chain_writer {
                    let read_len = query_len as u32;
                    let alts: Vec<String> = chain
                        .chain
                        .iter()
                        .map(|&seed_num| {
                            let seed = &compounds[seed_num];
                            let (left_clip, right_clip) = if seed.is_reverse() {
                                (read_len - seed.read_end(k), seed.read_start())
                            } else {
                                (seed.read_start(), read_len - seed.read_end(k))
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
                            let chrom = self.reference.chrom_name(seed.chrom_id() as usize);
                            let strand = if seed.is_reverse() { "-" } else { "+" };
                            let mapq = (seed.weight().floor() as i32).min(200);
                            format!(
                                "{},{},{},{}{}={},{},0;",
                                chrom,
                                seed.ref_start() + 1,
                                strand,
                                left,
                                seed.length(k),
                                right,
                                mapq
                            )
                        })
                        .collect();
                    let g = alts.len();
                    for (i, &seed_num) in chain.chain.iter().enumerate() {
                        let seed = &compounds[seed_num];

                        let b = seed.read_start() as usize;
                        let e = seed.read_end(k) as usize;
                        let sa_parts: Vec<String> = (0..g)
                            .filter(|v| *v != i)
                            .map(|v| alts[v].clone())
                            .collect();
                        let tags = vec![
                            (String::from("XG"), TagValue::Int(band_rank as i64)),
                            (String::from("XR"), TagValue::Int(chain_num as i64)),
                            (String::from("XS"), TagValue::Int(i as i64)),
                            (
                                String::from("XK"),
                                TagValue::Int(seed.multiplicity() as i64),
                            ),
                            (String::from("SA"), TagValue::Str(sa_parts.join(""))),
                        ];
                        let record = seed_to_record(
                            name,
                            k,
                            query_len,
                            seed,
                            &query[b..e],
                            &quality[b..e],
                            tags,
                        );
                        chain_writer
                            .write_record(&record)
                            .expect("chain write failed");
                    }
                }

                let n = chain.chain.len();
                let mut gap_alignments: Vec<Option<Alignment>> = vec![None; n];
                let mut gap_trims: Vec<usize> = vec![0; n];
                for rhs_rank in 1..n {
                    let lhs_rank = rhs_rank - 1;
                    if chain.edge_type[lhs_rank] != EdgeType::Continuation {
                        continue;
                    }
                    let lhs_seed = &compounds[chain.chain[lhs_rank]];
                    let rhs_seed = &compounds[chain.chain[rhs_rank]];

                    let ga = self.align_gap(k, lhs_seed, rhs_seed, query);
                    if let Some((aln, trim)) = ga {
                        //log::info!("gap alignment yielded: {}", aln.cigar_string());
                        gap_alignments[lhs_rank] = Some(aln);
                        gap_trims[rhs_rank] = trim;
                    }
                }

                let mut segments: Vec<AlignedSegment> = vec![];
                let mut segment_start_rank = 0;
                let mut current_segment: Vec<Alignment> = vec![];
                for i in 0..n {
                    let j = chain.chain[i];
                    let seed = &compounds[j];
                    let seed_alignment =
                        Alignment::from_perfect_match(seed.length(k) - gap_trims[i]);
                    current_segment.push(seed_alignment);
                    // Close the segment when there is no continuation gap after
                    // this seed (either a non-continuation edge, or the last seed).
                    let has_gap_after = gap_alignments[i].is_some();
                    if has_gap_after {
                        current_segment.push(gap_alignments[i].take().unwrap());
                    } else {
                        let first = &compounds[chain.chain[segment_start_rank]];
                        let last = &compounds[j];
                        let is_reverse = first.is_reverse();
                        let combined = Alignment::concat(&current_segment);
                        current_segment.clear();
                        // For reverse strand the pieces were accumulated in descending
                        // ref order; a single CIGAR reversal gives ref-ascending order.
                        let combined = if is_reverse { combined.reversed() } else { combined };
                        let chrom_id = first.chrom_id();
                        let read_start = first.read_start();
                        let read_end = last.read_end(k);
                        let (ref_start, ref_end) = if is_reverse {
                            (last.ref_start(), first.ref_end(k))
                        } else {
                            (first.ref_start(), last.ref_end(k))
                        };
                        segments.push(AlignedSegment {
                            chrom_id,
                            read_start,
                            read_end,
                            ref_start,
                            ref_end,
                            is_reverse,
                            alignment: combined,
                        });
                        segment_start_rank = i + 1;
                    }
                }
                assert!(current_segment.is_empty());

                let n = segments.len();
                let alts: Vec<String> = segments
                    .iter()
                    .map(|seg| seg.sa_tag_value(query_len as u32, self.reference))
                    .collect();

                for (seg_num, seg) in segments.iter().enumerate() {
                    let seg_alts: Vec<&str> = (0..n)
                        .filter(|&i| i != seg_num)
                        .map(|i| (&alts[i]) as &str)
                        .collect();
                    let seg_alt = seg_alts.join(";");

                    let mut tags: Vec<(String, TagValue)> = vec![];
                    if seg_alt.len() > 0 {
                        tags.push(("SA".to_string(), TagValue::Str(seg_alt)));
                    }

                    if let Err(e) = seg.validate() {
                        log::error!("segment validation failed: {}", e);
                        log::error!("  cigar: {}", seg.cigar_string(query_len as u32, false));
                        continue;
                    }

                    if let Err(e) = seg.validate_sequences(query, self.reference) {
                        let chrom_name = self.reference.chrom_name(seg.chrom_id as usize);
                        let strand = if seg.is_reverse { "-" } else { "+" };
                        log::error!(
                            "sequence validation failed: {} {}:{}-{} {}: {}",
                            name, chrom_name, seg.ref_start, seg.ref_end, strand, e
                        );
                        log::error!("  cigar: {}", seg.cigar_string(query_len as u32, false));
                        log::error!("  read: {}..{}", seg.read_start, seg.read_end);
                        continue;
                    }

                    let rec = seg.to_record(name, query, quality, tags, seg_num > 0, chain_num > 0, self.reference);

                    self.writer.write_record(&rec)?;
                    segment_written = true;
                }
                if self.no_secondary {
                    break;
                }
            }
        }

        if let Some(ref chain_writer) = self.chain_writer {
            chain_writer.finish().expect("finish failed");
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
        // seed_writer and chain_writer are shared across threads via Arc+OnceLock;
        // finishing them here would race with other threads still writing.
        // They are flushed when the last Arc is dropped at process exit.
        Ok(())
    }
}

/// Compute the query and reference slice boundaries for the gap between two consecutive seeds,
/// returning `(qs, qe, rs, re, trim)` where:
/// - `query[qs..qe]` is the query gap sequence
/// - `ref[rs..re]`   is the reference gap sequence (forward-strand coords; caller RC's for reverse)
/// - `trim`          is the number of bases trimmed from the leading edge of rhs_seed to resolve
///                   any overlap, used to shorten `rhs_seed.from_perfect_match` by the caller
///
/// Returns `None` if the gap is empty in both query and ref (nothing to align).
pub fn gap_regions<S: Seed>(lhs: &S, rhs: &S, k: usize) -> Option<(usize, usize, usize, usize, usize)> {
    let qs = lhs.read_end(k) as usize;
    let qe = rhs.read_start() as usize;
    let qg = (qe as isize) - (qs as isize);

    // For forward strand: gap ref region is [lhs.ref_end .. rhs.ref_start].
    // For reverse strand: seeds run ref-descending as read_pos increases, so
    //   lhs (lower read_pos) has higher ref and rhs has lower ref.
    //   Gap ref region is [rhs.ref_end .. lhs.ref_start].
    let (rs, re) = if lhs.is_reverse() {
        (rhs.ref_end(k) as usize, lhs.ref_start() as usize)
    } else {
        (lhs.ref_end(k) as usize, rhs.ref_start() as usize)
    };
    let rg = (re as isize) - (rs as isize);

    // If both gaps are zero there is nothing to align (seeds abut perfectly).
    if qg == 0 && rg == 0 {
        return None;
    }

    // Trim resolves overlap: the number of bases to advance the rhs seed's leading
    // edge in *read* space (qe += trim) and to correspondingly shrink the ref gap.
    let trim = -(qg.min(rg).min(0)) as usize;

    // Advance query rhs boundary.
    let qe = qe + trim;

    // Advance the ref boundary closest to the rhs seed by `trim` to eliminate the overlap.
    // Forward: gap is [lhs.ref_end .. rhs.ref_start]; rhs boundary is `re`; re += trim.
    // Reverse: gap is [rhs.ref_end .. lhs.ref_start]; rhs seed sits at the rs end, but the
    //   boundary that needs to advance to close the overlap is also `re` (lhs.ref_start),
    //   because moving re upward tightens the upper edge of the gap to match the query trim.
    let re = re + trim;

    Some((qs, qe, rs, re, trim))
}

impl<'a> ExplanatoryAligner<'a> {
    fn align_gap<S: Seed>(
        &mut self,
        k: usize,
        lhs_seed: &S,
        rhs_seed: &S,
        query: &[u8],
    ) -> Option<(Alignment, usize)> {
        let (qs, qe, rs, re, trim) = gap_regions(lhs_seed, rhs_seed, k)?;

        let q = &query[qs..qe];

        let mut ref_rc: Option<Vec<u8>> = None;
        let r = if lhs_seed.is_reverse() {
            let fwd = self.reference.get_seq(lhs_seed.chrom_id() as usize, rs, re);
            ref_rc = Some(fwd.iter().rev().map(|&base| complement(base)).collect());
            ref_rc.as_ref().unwrap()
        } else {
            let _ = ref_rc;
            self.reference.get_seq(lhs_seed.chrom_id() as usize, rs, re)
        };

        // Fast path for short identical sequences.
        let qg = qe as isize - qs as isize;
        let rg = re as isize - rs as isize;
        if qg == rg && qg < 64 {
            if q == r {
                return Some((Alignment::from_perfect_match(q.len()), trim));
            }
        }

        self.aligner.align(q, r).map(|aln| (aln, trim))
    }
}

/// Finalise any debug writers (seed/chain SAM/BAM files) that were opened during
/// alignment.  Must be called after all aligner threads have joined, before
/// process exit, to ensure the BGZF EOF block is written for BAM output.
pub fn finish_debug_writers() -> std::io::Result<()> {
    if let Some(w) = SEED_WRITER.get() {
        w.finish()?;
    }
    if let Some(w) = CHAIN_WRITER.get() {
        w.finish()?;
    }
    Ok(())
}

fn align_time_recorder() -> &'static SimpleSummaryRecorder {
    static RECORDER: OnceLock<&'static SimpleSummaryRecorder> = OnceLock::new();
    RECORDER.get_or_init(|| SimpleSummaryRecorder::new_registered("read_time"))
}

pub struct AlignedSegment {
    read_start: u32,
    read_end: u32,
    chrom_id: u32,
    ref_start: u32,
    ref_end: u32,
    is_reverse: bool,
    alignment: Alignment,
}

impl AlignedSegment {
    pub fn to_record(
        &self,
        name: &str,
        query: &[u8],
        qual: &[u8],
        tags: Vec<(String, TagValue)>,
        is_supplementary: bool,
        is_secondary: bool,
        _reference: &InMemoryReference,
    ) -> RecordBuf {
        let read_length = query.len() as u32;
        let left_clip = if self.is_reverse {
            read_length - self.read_end
        } else {
            self.read_start
        } as usize;
        let right_clip = if self.is_reverse {
            self.read_start
        } else {
            read_length - self.read_end
        } as usize;

        let mut cigar_ops: Vec<Op> = Vec::with_capacity(self.alignment.cigar.len() + 2);
        if left_clip > 0 {
            cigar_ops.push(Op::new(Kind::HardClip, left_clip));
        }
        cigar_ops.extend_from_slice(&self.alignment.cigar);
        if right_clip > 0 {
            cigar_ops.push(Op::new(Kind::HardClip, right_clip));
        }
        let cigar: Cigar = cigar_ops.into_iter().collect();

        // SEQ is always taken from the forward-strand query slice.
        // For reverse-strand records IGV RC's the SEQ to compare against the ref.
        let seq_slice = &query[self.read_start as usize..self.read_end as usize];
        let qual_slice = &qual[self.read_start as usize..self.read_end as usize];
        let mut seq_buf: Option<Vec<u8>> = None;
        let mut qual_buf: Option<Vec<u8>> = None;
        let (out_seq, out_qual): (&[u8], &[u8]) = if self.is_reverse {
            seq_buf = Some(seq_slice.iter().rev().map(|&b| complement(b)).collect());
            qual_buf = Some(qual_slice.iter().rev().copied().collect());
            (seq_buf.as_ref().unwrap(), qual_buf.as_ref().unwrap())
        } else {
            let _ = seq_buf;
            let _ = qual_buf;
            (seq_slice, qual_slice)
        };

        let mut flags = Flags::empty();
        if self.is_reverse {
            flags |= Flags::REVERSE_COMPLEMENTED;
        }
        if is_supplementary {
            flags |= Flags::SUPPLEMENTARY;
        }
        if is_secondary {
            flags |= Flags::SECONDARY;
        }

        let mapq = 0u8;

        let mut data_tags: Vec<(Tag, Value)> = Vec::with_capacity(tags.len());
        for (key, value) in tags {
            let bytes = key.as_bytes();
            if bytes.len() == 2 {
                let tag = Tag::from([bytes[0], bytes[1]]);
                let v = match value {
                    TagValue::Str(s) => Value::from(s.as_str()),
                    TagValue::Int(i) => Value::from(i as i32),
                    TagValue::Flt(f) => Value::from(f as f32),
                };
                data_tags.push((tag, v));
            }
        }
        let data: Data = data_tags.into_iter().collect();

        build_record(
            name,
            flags,
            self.chrom_id as usize,
            (self.ref_start + 1) as usize,
            mapq,
            cigar,
            None,
            None,
            &out_seq,
            &out_qual,
            data,
        )
    }

    /// Check that the alignment CIGAR is consistent with the stored read and
    /// reference extents.  Returns `Ok(())` or an error string describing the
    /// mismatch, suitable for logging before a panic or skip.
    pub fn validate(&self) -> Result<(), String> {
        let expected_query = (self.read_end - self.read_start) as usize;
        let actual_query = self.alignment.query_length();
        if actual_query != expected_query {
            return Err(format!(
                "CIGAR query length {} != read extent {} (read {}..{})",
                actual_query, expected_query, self.read_start, self.read_end
            ));
        }
        let expected_ref = (self.ref_end - self.ref_start) as u64;
        let actual_ref = self.alignment.reference_span();
        if actual_ref != expected_ref {
            return Err(format!(
                "CIGAR ref span {} != ref extent {} (ref {}..{})",
                actual_ref, expected_ref, self.ref_start, self.ref_end
            ));
        }
        Ok(())
    }

    /// Validate the alignment against the actual query and reference sequences.
    /// `query` is always the forward-strand read; the ref slice is RC'd internally
    /// for reverse-strand segments, matching the convention used by the aligner.
    pub fn validate_sequences(
        &self,
        query: &[u8],
        reference: &InMemoryReference,
    ) -> Result<(), String> {
        // The CIGAR (after reversal for reverse-strand segments) goes ref-ascending,
        // matching the SAM convention: SEQ = RC(query[read_start..read_end]) aligned
        // against the forward reference.  Validate using the same orientation.
        let ref_slice = reference
            .get_seq(self.chrom_id as usize, self.ref_start as usize, self.ref_end as usize)
            .to_vec();
        let query_slice: Vec<u8> = if self.is_reverse {
            query[self.read_start as usize..self.read_end as usize]
                .iter()
                .rev()
                .map(|&b| complement(b))
                .collect()
        } else {
            query[self.read_start as usize..self.read_end as usize].to_vec()
        };
        self.alignment.validate(&ref_slice, &query_slice, 0)
    }

    pub fn sa_tag_value(&self, read_length: u32, reference: &InMemoryReference) -> String {
        let chrom = reference.chrom_name(self.chrom_id as usize);
        let cigar = self.summary_cigar_string(read_length, false);
        let mapq = 0;
        let nm = self.alignment.edit_distance();
        format!(
            "{},{},{},{},{},{}",
            chrom,
            self.ref_start,
            if self.is_reverse { "-" } else { "+" },
            cigar,
            mapq,
            nm
        )
    }

    pub fn cigar_string(&self, read_length: u32, soft_clip: bool) -> String {
        let clip_char = if soft_clip { 'S' } else { 'H' };
        let mut parts = vec![];
        if self.read_start > 0 {
            parts.push(format!("{}{}", self.read_start, clip_char));
        }
        parts.push(self.alignment.cigar_string());
        if self.read_end < read_length {
            parts.push(format!("{}{}", read_length - self.read_end, clip_char));
        }
        parts.join("")
    }

    pub fn basic_cigar_string(&self, read_length: u32, soft_clip: bool) -> String {
        let clip_char = if soft_clip { 'S' } else { 'H' };
        let mut parts = vec![];
        if self.read_start > 0 {
            parts.push(format!("{}{}", self.read_start, clip_char));
        }
        parts.push(self.alignment.basic_cigar_string());
        if self.read_end < read_length {
            parts.push(format!("{}{}", read_length - self.read_end, clip_char));
        }
        parts.join("")
    }

    pub fn summary_cigar_string(&self, read_length: u32, soft_clip: bool) -> String {
        let clip_char = if soft_clip { 'S' } else { 'H' };
        let mut parts = vec![];
        if self.read_start > 0 {
            parts.push(format!("{}{}", self.read_start, clip_char));
        }
        parts.push(self.alignment.summary_cigar_string());
        if self.read_end < read_length {
            parts.push(format!("{}{}", read_length - self.read_end, clip_char));
        }
        parts.join("")
    }
}

#[cfg(test)]
mod tests {
    use super::gap_regions;
    use crate::reads::compound::Seed;

    struct TestSeed {
        read_pos: u32,
        ref_pos: u32,
        is_reverse: bool,
        len: usize,
    }

    impl Seed for TestSeed {
        fn read_pos(&self) -> u32 { self.read_pos }
        fn ref_pos(&self) -> u32 { self.ref_pos }
        fn chrom_id(&self) -> u32 { 0 }
        fn is_reverse(&self) -> bool { self.is_reverse }
        fn length(&self, _k: usize) -> usize { self.len }
        fn multiplicity(&self) -> u32 { 1 }
        fn to_string(&self, _k: usize) -> String { String::new() }
    }

    fn fwd(read_pos: u32, ref_pos: u32, len: usize) -> TestSeed {
        TestSeed { read_pos, ref_pos, is_reverse: false, len }
    }

    fn rev(read_pos: u32, ref_pos: u32, len: usize) -> TestSeed {
        // For reverse seeds, ref_pos is the lower ref coordinate of the k-mer.
        // ref_start() = ref_pos, ref_end(k) = ref_pos + len (via default trait impls).
        TestSeed { read_pos, ref_pos, is_reverse: true, len }
    }

    const K: usize = 0; // length comes from TestSeed.len, k unused in these seeds

    // Forward strand: clean gap, no overlap
    #[test]
    fn fwd_clean_gap() {
        // lhs: read[100..110] ref[200..210], rhs: read[115..125] ref[215..225]
        // query gap = 5, ref gap = 5, trim = 0
        let lhs = fwd(100, 200, 10);
        let rhs = fwd(115, 215, 10);
        let (qs, qe, rs, re, trim) = gap_regions(&lhs, &rhs, K).unwrap();
        assert_eq!(qs, 110);
        assert_eq!(qe, 115);
        assert_eq!(rs, 210);
        assert_eq!(re, 215);
        assert_eq!(trim, 0);
        assert_eq!(re - rs, 5); // gap ref span matches query gap span
    }

    // Forward strand: 2-base ref overlap resolved by trim
    #[test]
    fn fwd_ref_overlap() {
        // lhs: read[100..110] ref[200..210], rhs: read[112..122] ref[208..218]
        // query gap = 2, ref gap = -2, trim = 2
        let lhs = fwd(100, 200, 10);
        let rhs = fwd(112, 208, 10);
        let (qs, qe, rs, re, trim) = gap_regions(&lhs, &rhs, K).unwrap();
        assert_eq!(trim, 2);
        assert_eq!(qs, 110);
        assert_eq!(qe, 114); // 112 + trim=2
        assert_eq!(rs, 210);
        assert_eq!(re, 210); // 208 + trim=2
        assert_eq!(re as isize - rs as isize, 0); // zero ref gap after trim
    }

    // Forward strand: 2-base query overlap resolved by trim
    #[test]
    fn fwd_query_overlap() {
        // lhs: read[100..110] ref[200..210], rhs: read[108..118] ref[212..222]
        // query gap = -2, ref gap = 2, trim = 2
        let lhs = fwd(100, 200, 10);
        let rhs = fwd(108, 212, 10);
        let (qs, qe, rs, re, trim) = gap_regions(&lhs, &rhs, K).unwrap();
        assert_eq!(trim, 2);
        assert_eq!(qs, 110); // lhs.read_end = 100 + 10, unchanged by trim
        assert_eq!(qe, 110); // 108 + 2
        assert_eq!(rs, 210); // lhs.ref_end = 200 + 10, unchanged by trim
        assert_eq!(re, 214); // 212 + 2
    }

    // Reverse strand: clean gap (lhs has higher ref, rhs has lower ref)
    #[test]
    fn rev_clean_gap() {
        // lhs: read[100..110] ref[300..310] (ref_pos=300, ref_end=310)
        // rhs: read[115..125] ref[285..295] (ref_pos=285, ref_end=295)
        // gap ref region: [rhs.ref_end..lhs.ref_start] = [295..300], span=5
        // gap query: [lhs.read_end..rhs.read_start] = [110..115], span=5
        let lhs = rev(100, 300, 10);
        let rhs = rev(115, 285, 10);
        let (qs, qe, rs, re, trim) = gap_regions(&lhs, &rhs, K).unwrap();
        assert_eq!(qs, 110);
        assert_eq!(qe, 115);
        assert_eq!(rs, 295); // rhs.ref_end = 285+10
        assert_eq!(re, 300); // lhs.ref_start = 300
        assert_eq!(trim, 0);
        assert_eq!(re - rs, 5);
    }

    // Reverse strand: 1-base ref overlap, trim=1
    #[test]
    fn rev_ref_overlap() {
        // lhs: read[100..110] ref[300..310]
        // rhs: read[112..122] ref[299..309] → rhs.ref_end = 309 > lhs.ref_start=300? No.
        // rhs.ref_end = 309, lhs.ref_start = 300: rg = 300-309 = -9 — too much.
        // Try: rhs.ref_end = 301, lhs.ref_start = 300 → rg = -1, trim=1
        // rhs: read[112..122] ref[291..301]
        let lhs = rev(100, 300, 10);
        let rhs = rev(112, 291, 10); // ref_end = 301
        // query gap = 112-110=2, ref gap = 300-301=-1, trim=1
        let (qs, qe, rs, re, trim) = gap_regions(&lhs, &rhs, K).unwrap();
        assert_eq!(trim, 1);
        assert_eq!(qs, 110); // lhs.read_end = 100 + 10
        assert_eq!(qe, 113); // 112 + 1
        assert_eq!(rs, 301); // rhs.ref_end = 291 + 10, unchanged by trim
        assert_eq!(re, 301); // lhs.ref_start=300 + trim=1 → 301
        assert_eq!(re as isize - rs as isize, 0); // zero gap after trim
    }
}
