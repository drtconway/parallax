use std::sync::{Arc, OnceLock};

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
            all_seeds: Vec::new(),
            no_secondary: self.no_secondary,
            seeding_cfg: cfg.seeding.clone(),
            filtering_cfg: cfg.filtering.clone(),
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
    all_seeds: Vec<ExtendedSeed>,
    no_secondary: bool,
    seeding_cfg: SeedingConfig,
    filtering_cfg: FilteringConfig,
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

impl<'a> ExplanatoryAligner<'a> {
    fn align_gap<S: Seed>(
        &mut self,
        k: usize,
        lhs_seed: &S,
        rhs_seed: &S,
        query: &[u8],
    ) -> Option<(Alignment, usize)> {
        let qs = lhs_seed.read_end(k) as usize;
        let qe = rhs_seed.read_start() as usize;
        let qg = (qe as isize) - (qs as isize);

        let (rs, re) = if lhs_seed.is_reverse() {
            (rhs_seed.ref_end(k) as usize, lhs_seed.ref_start() as usize)
        } else {
            (lhs_seed.ref_end(k) as usize, rhs_seed.ref_start() as usize)
        };
        let rg = (re as isize) - (rs as isize);

        let trim = -(qg.min(rg).min(0)) as usize;

        // Advance the rhs boundary by `trim` to eliminate any overlap.
        // On the query, rhs starts at qe — move it forward.
        // On the reference, the rhs boundary is re (forward) or rs (reverse).
        let qe = qe + trim;
        let (rs, re) = if lhs_seed.is_reverse() {
            (rs + trim, re)
        } else {
            (rs, re + trim)
        };

        let q = &query[qs..qe];

        let mut ref_rc: Option<Vec<u8>> = None; // storage for the rc if required.
        let r = if lhs_seed.is_reverse() {
            let fwd = self.reference.get_seq(lhs_seed.chrom_id() as usize, rs, re);
            ref_rc = Some(fwd.iter().rev().map(|&base| complement(base)).collect());
            ref_rc.as_ref().unwrap()
        } else {
            let _ = ref_rc;
            self.reference.get_seq(lhs_seed.chrom_id() as usize, rs, re)
        };

        // Fast path for the common case where the two sequences are short and equal.
        if qg == rg && qg < 64 {
            if q == r {
                let aln = Alignment::from_perfect_match(q.len());
                return Some((aln, trim));
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
        reference: &InMemoryReference,
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
