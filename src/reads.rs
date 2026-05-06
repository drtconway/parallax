#[cfg(feature = "conventional")]
use core::panic;
#[cfg(feature = "conventional")]
use std::io::Write;
use std::{sync::Arc, usize};

#[cfg(feature = "conventional")]
use ordered_float::OrderedFloat;

use crate::{
    Aligner, AlignerBuilder,
    error::Result,
    explanatory,
    index::Index,
    reference::InMemoryReference,
    utils::{debug, sequence::reverse_complement_into},
    writer::{AlignmentWriter, OutputFormat},
};

#[cfg(feature = "conventional")]
use crate::reads::seeds::analyze_gap_fills;
#[cfg(feature = "conventional")]
use crate::reads::seeds::seed_cluster::SeedCluster;
#[cfg(feature = "conventional")]
use crate::scores::compute_mapq_from_diff;
#[cfg(feature = "conventional")]
use crate::utils::GroupByTrait;
#[cfg(feature = "conventional")]
use crate::utils::heap::{Heap, HeapOrdering, Heapable};
#[cfg(feature = "conventional")]
use crate::utils::range_set::RangeSet;

pub mod builder;
#[cfg(feature = "conventional")]
pub mod chains;
pub mod extended;
pub mod seeds;

#[cfg(feature = "conventional")]
#[derive(Debug)]
enum AlignmentError {
    #[allow(dead_code)]
    NoClusters,
    #[allow(dead_code)]
    LowQuality,
}

/// Detected input file format.
enum InputFormat {
    Fastq,
    Bam,
}

/// Detect input format from file extension.
fn detect_input_format(path: &str) -> InputFormat {
    let path_lower = path.to_lowercase();
    if path_lower.ends_with(".bam") || path_lower.ends_with(".ubam") {
        InputFormat::Bam
    } else {
        // Default to FASTQ for .fq, .fastq, .fq.gz, .fastq.gz, etc.
        InputFormat::Fastq
    }
}

/// A read to be processed by a worker thread
struct ReadWork {
    name: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

/// Process reads from a FASTQ or unaligned BAM file using multiple threads.
///
/// The input format is auto-detected from the file extension:
/// - `.bam` → unaligned BAM
/// - anything else → FASTQ (with optional gzip/bzip2/xz compression)
///
/// Reads are distributed to worker threads via a channel. The InMemoryReference
/// is shared across all threads via Arc (no per-thread cloning needed).
pub fn process_reads_parallel<const K: usize, const S: usize>(
    index: &Index<K, S>,
    reference: &InMemoryReference,
    reads: &str,
    sam: Option<&str>,
    num_threads: usize,
    command_line: &str,
    read_group_header: Option<&str>,
    output_format: OutputFormat,
) -> Result<()> {
    use crossbeam::channel::bounded;

    let format = detect_input_format(reads);
    log::info!(
        "Processing reads from {} ({}) using {} threads, output format: {}",
        reads,
        match format {
            InputFormat::Fastq => "FASTQ",
            InputFormat::Bam => "BAM",
        },
        num_threads,
        output_format,
    );

    let now = std::time::Instant::now();
    let mut num_records = 0;

    // Store reference chromosome info for debug SAM headers
    debug::set_reference_info(reference.chromosomes());

    // Create writer with headers - either to file or stdout
    let output: Box<dyn std::io::Write + Send> = match sam {
        Some(path) => {
            log::info!("Writing output to {}", path);
            Box::new(std::fs::File::create(path)?)
        }
        None => Box::new(std::io::stdout()),
    };
    let writer = Arc::new(
        AlignmentWriter::builder(output, output_format, reference.to_fasta_repository())
            .add_contigs(reference.chromosomes())
            .read_group(read_group_header.map(String::from))
            .command_line(command_line)
            .build()?,
    );

    // Create a bounded channel for backpressure
    let (sender, receiver) = bounded::<ReadWork>(num_threads * 100);

    // Use crossbeam's scoped threads to safely borrow index, reference, and writer
    crossbeam::scope(|scope| {
        // Spawn worker threads
        for _ in 0..num_threads {
            let receiver = receiver.clone();
            let writer = writer.clone();
            scope.spawn(move |_| {
                let mut aligner =
                    explanatory::ExplanatoryAlignerBuilder::new(&reference, index, &writer).build();
                while let Ok(work) = receiver.recv() {
                    aligner
                        .align(&work.name, &work.seq, &work.qual)
                        .expect("alignment failed");
                }
            });
        }

        match format {
            InputFormat::Fastq => {
                // Read FASTQ and send to workers
                // Use niffler for transparent decompression (gzip, bzip2, xz)
                let (decompressed_reader, compression) =
                    niffler::from_path(std::path::Path::new(reads))
                        .expect("Failed to open FASTQ file");
                if compression != niffler::Format::No {
                    log::info!("Detected {:?} compression", compression);
                }
                let reader = std::io::BufReader::new(decompressed_reader);
                let mut reader = noodles::fastq::io::Reader::new(reader);

                for record in reader.records() {
                    let record = record.expect("Failed to read FASTQ record");
                    let seq: &[u8] = record.sequence().as_ref();
                    let qual: &[u8] = record.quality_scores().as_ref();
                    let work = ReadWork {
                        name: String::from_utf8_lossy(record.name()).into_owned(),
                        seq: seq.to_vec(),
                        qual: qual.to_vec(),
                    };
                    sender.send(work).expect("Failed to send work to thread");
                    num_records += 1;
                }
            }
            InputFormat::Bam => {
                // Read unaligned BAM and send to workers
                let file = std::fs::File::open(reads).expect("Failed to open BAM file");
                let mut reader = noodles::bam::io::Reader::new(file);
                let header = reader.read_header().expect("Failed to read BAM header");

                let mut rc_buf = Vec::new();

                for result in reader.record_bufs(&header) {
                    let record = result.expect("Failed to read BAM record");

                    let raw_seq: Vec<u8> = record.sequence().as_ref().iter().cloned().collect();
                    if raw_seq.is_empty() {
                        continue;
                    }

                    let name = record
                        .name()
                        .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                        .unwrap_or_else(|| format!("unnamed_{}", num_records));

                    let is_reverse = record.flags().is_reverse_complemented();
                    let raw_qual: Vec<u8> = record.quality_scores().as_ref().to_vec();

                    let (seq, qual) = if is_reverse {
                        // Undo reverse complement applied by previous aligner:
                        // reverse-complement the sequence and reverse the quality.
                        reverse_complement_into(&raw_seq, &mut rc_buf);
                        let seq = rc_buf.clone();
                        let qual: Vec<u8> = raw_qual
                            .iter()
                            .rev()
                            .map(|&q| q.saturating_add(33))
                            .collect();
                        (seq, qual)
                    } else {
                        // Convert quality from raw Phred (BAM) to Phred+33 (SAM/FASTQ).
                        let qual: Vec<u8> =
                            raw_qual.iter().map(|&q| q.saturating_add(33)).collect();
                        (raw_seq, qual)
                    };

                    // Handle missing quality scores (all 0xFF in BAM → empty after decode)
                    let qual = if qual.is_empty() {
                        vec![b'!'; seq.len()] // Phred 0 + 33 = '!' as placeholder
                    } else {
                        qual
                    };

                    let work = ReadWork { name, seq, qual };
                    sender.send(work).expect("Failed to send work to thread");
                    num_records += 1;
                }
            }
        }

        // Signal completion by dropping sender
        drop(sender);

        // Scoped threads automatically join when scope ends
    })
    .expect("Scoped thread panicked");

    writer.finish()?;

    let elapsed = now.elapsed();
    log::info!(
        "Completed processing {} reads from {} in {:.2?}",
        num_records,
        reads,
        elapsed
    );

    Ok(())
}

/// Compute the "core" read-end of a cluster by trimming minor tail seeds.
///
/// Cluster boundaries (`read_end`) can be inflated by small noisy seeds beyond
/// the main alignment (e.g., false-positive matches in a deleted region). This
/// function finds the rightmost seed that's part of the cluster's main mass by
/// trimming up to 10% of total seed length from the tail.
#[cfg(feature = "conventional")]
fn core_read_end(cluster: &SeedCluster) -> usize {
    let total: usize = cluster.chain.iter().map(|s| s.match_len).sum();
    if total == 0 {
        return cluster.read_end;
    }
    let trim_budget = total / 10;
    let mut trimmed = 0;
    for seed in cluster.chain.iter().rev() {
        if trimmed + seed.match_len > trim_budget {
            return seed.read_end();
        }
        trimmed += seed.match_len;
    }
    cluster.read_end
}

/// Compute the "core" read-start of a cluster by trimming minor head seeds.
#[cfg(feature = "conventional")]
fn core_read_start(cluster: &SeedCluster) -> usize {
    let total: usize = cluster.chain.iter().map(|s| s.match_len).sum();
    if total == 0 {
        return cluster.read_start;
    }
    let trim_budget = total / 10;
    let mut trimmed = 0;
    for seed in &cluster.chain {
        if trimmed + seed.match_len > trim_budget {
            return seed.read_pos;
        }
        trimmed += seed.match_len;
    }
    cluster.read_start
}

/// Merge clusters on the same chromosome and strand that are colinear and
/// nearly abutting on the read, indicating a deletion that the chainer
/// could not bridge (due to `MAX_DIAGONAL_DIST`).
///
/// For a genuine deletion the read is continuous across the breakpoint, so
/// the two clusters should be almost exactly adjacent on the read. Small
/// overlaps (≤ `del_merge_max_read_overlap`) are permitted for microhomology
/// at breakpoints; small gaps (≤ `del_merge_max_read_gap`) for seed placement
/// granularity.
///
/// Cluster boundaries (`read_start`/`read_end`) can be inflated by noisy
/// tail seeds that mapped into the deleted region, so the adjacency check
/// uses **core boundaries** — the extent of seeds near the cluster's
/// length-weighted median diagonal — rather than outermost seed positions.
///
/// After merging seed lists, `SeedCluster::new()` re-resolves overlaps and
/// filters noisy seeds, so the merged cluster is clean. The resulting gap
/// between the two seed groups will be aligned by `align_gaps()` and typically
/// produces a pure deletion (`D`) CIGAR operation.
#[cfg(feature = "conventional")]
fn merge_deletion_clusters(
    clusters: &mut Vec<SeedCluster>,
    read_len: usize,
    min_seed_length: usize,
) {
    let cfg = config::get();
    let max_read_overlap = cfg.seeding.del_merge_max_read_overlap as i64;
    let max_read_gap = cfg.seeding.del_merge_max_read_gap as i64;
    let max_ref_gap = cfg.seeding.del_merge_max_ref_gap as i64;

    if clusters.len() < 2 {
        return;
    }

    let mut merged_any = true;
    while merged_any {
        merged_any = false;

        // Sort by (chrom_id, is_reverse, read_start) in strand coordinates
        clusters.sort_by_key(|c| (c.chrom_id, c.is_reverse, c.read_start));

        let mut i = 0;
        while i + 1 < clusters.len() {
            let can_merge = {
                let a = &clusters[i];
                let b = &clusters[i + 1];

                if a.chrom_id != b.chrom_id
                    || a.is_reverse != b.is_reverse
                    || a.ref_end() > b.ref_start()
                    || (b.ref_start() as i64 - a.ref_end() as i64) > max_ref_gap
                {
                    false
                } else {
                    // Use core boundaries: the read extent of seeds close to
                    // each cluster's median diagonal. This ignores noisy tail
                    // seeds that mapped into the deleted region.
                    let a_core_end = core_read_end(a);
                    let b_core_start = core_read_start(b);
                    let read_delta = b_core_start as i64 - a_core_end as i64;
                    read_delta >= -max_read_overlap && read_delta <= max_read_gap
                }
            };

            if can_merge {
                let is_reverse = clusters[i].is_reverse;
                let mut combined = clusters[i].chain.clone();
                combined.extend(clusters[i + 1].chain.iter().cloned());

                if let Some(new_cluster) = SeedCluster::new(combined, is_reverse, min_seed_length) {
                    log::info!(
                        "Merged colinear clusters: read [{}-{}]+[{}-{}] -> [{}-{}], \
                         ref [{}-{}]+[{}-{}], gap {}bp",
                        clusters[i].read_start,
                        clusters[i].read_end,
                        clusters[i + 1].read_start,
                        clusters[i + 1].read_end,
                        new_cluster.read_start,
                        new_cluster.read_end,
                        clusters[i].ref_start(),
                        clusters[i].ref_end(),
                        clusters[i + 1].ref_start(),
                        clusters[i + 1].ref_end(),
                        clusters[i + 1].ref_start() as i64 - clusters[i].ref_end() as i64,
                    );
                    clusters[i] = new_cluster;
                    clusters.remove(i + 1);
                    merged_any = true;
                    continue; // retry at same position
                }
            }
            i += 1;
        }
    }

    // Re-sort by fwd_read_range for the rest of the pipeline
    clusters.sort_by_key(|c| c.fwd_read_range(read_len));
}

#[cfg(test)]
mod tests {
    use crate::reads::seeds::SeedHit;

    use super::*;

    // Helper to create a SeedHit with a dummy kmer value
    fn make_hit(chrom_id: usize, ref_pos: usize, read_pos: usize, match_len: usize) -> SeedHit {
        SeedHit::new(chrom_id, ref_pos, read_pos, 0, 1, match_len)
    }

    #[test]
    fn test_seed_hit_new() {
        let hit = SeedHit::new(1, 100, 50, 12345, 1, 20);
        assert_eq!(hit.chrom_id, 1);
        assert_eq!(hit.ref_pos, 100);
        assert_eq!(hit.read_pos, 50);
        assert_eq!(hit.kmer, 12345);
        assert_eq!(hit.kmer_uniqueness, 1);
        assert_eq!(hit.match_len, 20);
        // diagonal = ref_pos - read_pos = 100 - 50 = 50
        assert_eq!(hit.diagonal, 50);
    }

    #[test]
    fn test_seed_hit_read_end() {
        let hit = make_hit(0, 100, 50, 20);
        assert_eq!(hit.read_end(), 70); // 50 + 20
    }

    #[test]
    fn test_seed_hit_ref_end() {
        let hit = make_hit(0, 100, 50, 20);
        assert_eq!(hit.ref_end(), 120); // 100 + 20
    }

    #[test]
    fn test_seed_hit_diagonal_calculation() {
        // Forward diagonal (ref ahead of read)
        let hit1 = make_hit(0, 1000, 100, 20);
        assert_eq!(hit1.diagonal, 900);

        // Negative diagonal (read ahead of ref)
        let hit2 = make_hit(0, 100, 1000, 20);
        assert_eq!(hit2.diagonal, -900);

        // Zero diagonal (same position)
        let hit3 = make_hit(0, 500, 500, 20);
        assert_eq!(hit3.diagonal, 0);
    }

    #[test]
    fn test_extend_overlapping_same_diagonal() {
        // First seed at read pos 0, ref pos 100, length 20
        let mut hit = make_hit(0, 100, 0, 20);

        // Second seed at read pos 10, ref pos 110, length 20
        // This is on the same diagonal (110-10 = 100-0 = 100)
        // And overlaps: ref 110 < ref_end 120, gap is 10 < k=20
        let k = 20;
        let result = hit.extend(0, 110, 10, 0, 1, 1, k);

        assert!(
            result.is_none(),
            "Should extend in place, not return new hit"
        );
        // New end should be: (110 - 100) + 20 = 30
        assert_eq!(hit.match_len, 30);
        assert_eq!(hit.ref_end(), 130);
        assert_eq!(hit.read_end(), 30);
    }

    #[test]
    fn test_extend_adjacent_same_diagonal() {
        // First seed at read pos 0, ref pos 100, length 20
        let mut hit = make_hit(0, 100, 0, 20);

        // Second seed starts exactly where first ends
        // read pos 20, ref pos 120, still same diagonal
        let k = 20;
        let result = hit.extend(0, 120, 20, 0, 1, 1, k);

        // Gap is exactly 20, which equals k, so should NOT extend
        // because condition is: chrom_pos - self.ref_pos < self.match_len + k
        // 120 - 100 = 20 < 20 + 20 = 40, so it SHOULD extend
        assert!(result.is_none(), "Should extend in place");
        assert_eq!(hit.match_len, 40);
    }

    #[test]
    fn test_extend_gap_too_large() {
        // First seed at read pos 0, ref pos 100, length 20
        let mut hit = make_hit(0, 100, 0, 20);
        let original_len = hit.match_len;

        // Second seed with large gap (beyond match_len + k)
        // ref pos 200, read pos 100 (same diagonal = 100)
        // Gap check: 200 - 100 = 100 >= 20 + 20 = 40
        let k = 20;
        let result = hit.extend(0, 200, 100, 999, 1, 1, k);

        assert!(result.is_some(), "Should return new hit due to large gap");
        assert_eq!(hit.match_len, original_len, "Original should be unchanged");

        let new_hit = result.unwrap();
        assert_eq!(new_hit.ref_pos, 200);
        assert_eq!(new_hit.read_pos, 100);
        assert_eq!(new_hit.kmer, 999);
    }

    #[test]
    fn test_extend_different_chromosome() {
        let mut hit = make_hit(0, 100, 0, 20);
        let original_len = hit.match_len;

        // Same positions but different chromosome
        let k = 20;
        let result = hit.extend(1, 110, 10, 0, 1, 1, k);

        assert!(
            result.is_some(),
            "Different chromosome should create new hit"
        );
        assert_eq!(hit.match_len, original_len);
        assert_eq!(result.unwrap().chrom_id, 1);
    }

    #[test]
    fn test_extend_different_diagonal() {
        let mut hit = make_hit(0, 100, 0, 20);
        let original_len = hit.match_len;

        // Different diagonal: ref_pos - read_pos = 111 - 10 = 101 != 100
        let k = 20;
        let result = hit.extend(0, 111, 10, 0, 1, 1, k);

        assert!(result.is_some(), "Different diagonal should create new hit");
        assert_eq!(hit.match_len, original_len);

        let new_hit = result.unwrap();
        assert_eq!(new_hit.diagonal, 101);
    }

    #[test]
    fn test_extend_backwards_ref_position() {
        let mut hit = make_hit(0, 100, 50, 20);
        let original_len = hit.match_len;

        // New ref_pos before current ref_pos
        let k = 20;
        let result = hit.extend(0, 90, 40, 0, 1, 1, k);

        assert!(
            result.is_some(),
            "Backwards ref position should create new hit"
        );
        assert_eq!(hit.match_len, original_len);
    }

    #[test]
    fn test_extend_backwards_read_position() {
        let mut hit = make_hit(0, 100, 50, 20);
        let original_len = hit.match_len;

        // New read_pos before current read_pos (even if same diagonal)
        let k = 20;
        let result = hit.extend(0, 90, 40, 0, 1, 1, k);

        assert!(
            result.is_some(),
            "Backwards read position should create new hit"
        );
        assert_eq!(hit.match_len, original_len);
    }

    #[test]
    fn test_extend_fully_contained() {
        // Seed covering positions 0-20 in read, 100-120 in ref
        let mut hit = make_hit(0, 100, 0, 20);

        // New seed at pos 5-25 overlaps significantly
        // ref 105, read 5, same diagonal (100)
        let k = 20;
        let result = hit.extend(0, 105, 5, 0, 1, 1, k);

        assert!(result.is_none(), "Overlapping hit should extend");
        // New end: (105 - 100) + 20 = 25
        assert_eq!(hit.match_len, 25);
    }

    #[test]
    fn test_extend_no_length_change_if_contained() {
        // Seed covering positions 0-30 in read
        let mut hit = make_hit(0, 100, 0, 30);

        // New seed fully contained within existing match
        // ref 110, read 10, k=20 means it ends at read 30, ref 130
        // That's exactly where the original ends, so no extension needed
        let k = 20;
        let result = hit.extend(0, 110, 10, 0, 1, 1, k);

        assert!(result.is_none(), "Contained hit should not create new hit");
        // (110 - 100) + 20 = 30, which equals original, so no change
        assert_eq!(hit.match_len, 30);
    }

    #[test]
    fn test_extend_sequence_of_hits() {
        let k = 20;
        let mut hit = make_hit(0, 100, 0, k);

        // Simulate a sequence of overlapping syncmers ~6 bases apart
        // All on the same diagonal
        for i in 1..10 {
            let read_pos = i * 6;
            let ref_pos = 100 + i * 6;
            let result = hit.extend(0, ref_pos, read_pos, 0, 1, 1, k);
            assert!(result.is_none(), "Hit {} should extend in place", i);
        }

        // Final length should cover from 0 to (9*6 + 20) = 74
        assert_eq!(hit.match_len, 9 * 6 + k);
        assert_eq!(hit.read_end(), 74);
        assert_eq!(hit.ref_end(), 174);
    }

    // =========================================================================
    // SeedCluster tests (conventional feature only)
    // =========================================================================

    #[cfg(feature = "conventional")]
    mod seed_cluster_tests {
        use super::*;

        #[test]
        fn test_seed_cluster_new_sorts_by_read_pos() {
            // Create seeds in reverse read_pos order
            let seeds = vec![
                make_hit(0, 300, 200, 20), // read_pos = 200
                make_hit(0, 100, 0, 20),   // read_pos = 0
                make_hit(0, 200, 100, 20), // read_pos = 100
            ];

            let cluster = SeedCluster::new(seeds, false, 1).unwrap();

            // Should be sorted by read_pos
            assert_eq!(cluster.chain[0].read_pos, 0);
            assert_eq!(cluster.chain[1].read_pos, 100);
            assert_eq!(cluster.chain[2].read_pos, 200);

            assert_eq!(cluster.read_start, 0);
            assert_eq!(cluster.read_end, 220); // 200 + 20
        }

        #[test]
        fn test_seed_cluster_empty_returns_none() {
            let seeds: Vec<SeedHit> = vec![];
            assert!(SeedCluster::new(seeds, false, 1).is_none());
        }

        #[test]
        fn test_seed_cluster_fwd_read_range_forward_strand() {
            let seeds = vec![make_hit(0, 100, 50, 20)];
            let cluster = SeedCluster::new(seeds, false, 1).unwrap();

            let (start, end) = cluster.fwd_read_range(1000);
            assert_eq!(start, 50);
            assert_eq!(end, 70);
        }

        #[test]
        fn test_seed_cluster_fwd_read_range_reverse_strand() {
            // For reverse strand, coordinates need to be flipped
            let seeds = vec![make_hit(0, 100, 50, 20)];
            let cluster = SeedCluster::new(seeds, true, 1).unwrap();

            // read_start=50, read_end=70, read_len=1000
            // fwd_start = 1000 - 70 = 930, fwd_end = 1000 - 50 = 950
            let (start, end) = cluster.fwd_read_range(1000);
            assert_eq!(start, 930);
            assert_eq!(end, 950);
        }

        #[test]
        fn test_seed_cluster_split_at_gap() {
            use crate::align::{Alignment, Kind, Op};
            use crate::scores::DivergenceScore;

            // Create a chain with a gap between seeds 1 and 2
            let seeds = vec![
                make_hit(0, 100, 0, 20),   // 0-20
                make_hit(0, 200, 50, 20),  // 50-70
                make_hit(0, 400, 200, 20), // 200-220 (gap here)
                make_hit(0, 500, 250, 20), // 250-270
            ];

            let mut cluster = SeedCluster::new(seeds, false, 1).unwrap();

            // Add dummy gap alignments (3 gaps for 4 seeds)
            cluster.gap_alignments = vec![
                Alignment {
                    divergence: DivergenceScore::new(10.0),
                    cigar: vec![Op::new(Kind::SequenceMatch, 10)],
                },
                Alignment {
                    divergence: DivergenceScore::new(20.0),
                    cigar: vec![Op::new(Kind::SequenceMatch, 20)],
                },
                Alignment {
                    divergence: DivergenceScore::new(15.0),
                    cigar: vec![Op::new(Kind::SequenceMatch, 15)],
                },
            ];

            assert_eq!(cluster.chain.len(), 4);

            // Split at gap between index 1 and 2
            let (tail, _dropped_alignment) = cluster.split_at_gap(1).unwrap();

            // Original cluster should have seeds 0, 1
            assert_eq!(cluster.chain.len(), 2);
            assert_eq!(cluster.read_start, 0);
            assert_eq!(cluster.read_end, 70);

            // Tail cluster should have seeds 2, 3
            assert_eq!(tail.chain.len(), 2);
            assert_eq!(tail.read_start, 200);
            assert_eq!(tail.read_end, 270);
            assert_eq!(tail.is_reverse, cluster.is_reverse);
        }

        #[test]
        fn test_seed_cluster_split_preserves_strand() {
            use crate::align::{Alignment, Kind, Op};
            use crate::scores::DivergenceScore;

            let seeds = vec![make_hit(0, 100, 0, 20), make_hit(0, 300, 100, 20)];

            let mut cluster = SeedCluster::new(seeds, true, 1).unwrap();

            // Add dummy gap alignment (1 gap for 2 seeds)
            cluster.gap_alignments = vec![Alignment {
                divergence: DivergenceScore::new(10.0),
                cigar: vec![Op::new(Kind::SequenceMatch, 10)],
            }];

            let (tail, _dropped_alignment) = cluster.split_at_gap(0).unwrap();

            assert!(cluster.is_reverse);
            assert!(tail.is_reverse);
        }

        #[test]
        fn test_seed_cluster_gaps_forward_strand() {
            let seeds = vec![
                make_hit(0, 100, 0, 20),   // 0-20
                make_hit(0, 200, 50, 20),  // 50-70, gap of 30
                make_hit(0, 400, 200, 20), // 200-220, gap of 130
            ];

            let cluster = SeedCluster::new(seeds, false, 1).unwrap();
            let gaps: Vec<_> = cluster.gaps(1000, 50).collect();

            // Only the gap of 130 should be returned (min_gap=50)
            assert_eq!(gaps.len(), 1);
            assert_eq!(gaps[0].0, (70, 200)); // gap from 70 to 200
            assert_eq!(gaps[0].1, 1); // seed index 1 before this gap
        }

        #[test]
        fn test_seed_cluster_gaps_reverse_strand() {
            // For reverse strand, the chain is in RC coordinates
            // but gaps() should return forward-strand coordinates
            let seeds = vec![
                make_hit(0, 100, 0, 20),   // RC pos 0-20
                make_hit(0, 200, 50, 20),  // RC pos 50-70, gap of 30
                make_hit(0, 400, 200, 20), // RC pos 200-220, gap of 130
            ];

            let cluster = SeedCluster::new(seeds, true, 1).unwrap();
            let read_len = 1000;
            let gaps: Vec<_> = cluster.gaps(read_len, 50).collect();

            // Gap in RC coords: 70-200
            // In forward coords: (1000-200, 1000-70) = (800, 930)
            assert_eq!(gaps.len(), 1);
            assert_eq!(gaps[0].0, (800, 930));
        }

        // =========================================================================
        // Chain colinearity tests
        // =========================================================================

        #[test]
        fn test_colinear_chain_both_dimensions_increasing() {
            // A proper colinear chain should have both ref_pos and read_pos increasing
            let seeds = vec![
                make_hit(0, 100, 10, 20),
                make_hit(0, 200, 50, 20),
                make_hit(0, 300, 100, 20),
                make_hit(0, 400, 150, 20),
            ];

            let cluster = SeedCluster::new(seeds, false, 1).unwrap();

            // Verify both dimensions are strictly increasing
            for i in 1..cluster.chain.len() {
                assert!(
                    cluster.chain[i].read_pos > cluster.chain[i - 1].read_end() - 1,
                    "read_pos should be increasing: {} vs {}",
                    cluster.chain[i].read_pos,
                    cluster.chain[i - 1].read_end()
                );
                assert!(
                    cluster.chain[i].ref_pos >= cluster.chain[i - 1].ref_end(),
                    "ref_pos should be increasing: {} vs {}",
                    cluster.chain[i].ref_pos,
                    cluster.chain[i - 1].ref_end()
                );
            }
        }

        #[test]
        fn test_chain_ref_pos_monotonic_after_read_sort() {
            // This tests the invariant that should hold after SeedCluster::new
            // Even if seeds come in arbitrary order, after sorting by read_pos,
            // ref_pos should also be increasing for a proper colinear chain
            let seeds = vec![
                make_hit(0, 400, 150, 20), // Will be last after sort
                make_hit(0, 100, 10, 20),  // Will be first after sort
                make_hit(0, 300, 100, 20), // Will be third after sort
                make_hit(0, 200, 50, 20),  // Will be second after sort
            ];

            let cluster = SeedCluster::new(seeds, false, 1).unwrap();

            // Verify ref_pos is monotonically increasing
            for i in 1..cluster.chain.len() {
                assert!(
                    cluster.chain[i].ref_pos >= cluster.chain[i - 1].ref_pos,
                    "ref_pos not monotonic at {}: {} < {}",
                    i,
                    cluster.chain[i].ref_pos,
                    cluster.chain[i - 1].ref_pos
                );
            }
        }

        #[test]
        fn test_filter_misplaced_seeds_removes_bad_region() {
            // Construct a chain where seeds 2 & 3 are misplaced:
            // They introduce a big insertion gap followed by a big deletion gap
            // (or vice versa), indicating contradictory anchor placement.
            //
            // Good seeds are on diagonal ~0 (ref_pos ≈ read_pos).
            // Misplaced seeds shift the diagonal creating simultaneous ins + del.
            //
            // Layout (read_pos, ref_pos) for 5 seeds of length 20:
            //   seed 0: read 0,   ref 0     → gap to seed 1: read_delta 30, ref_delta 30, gap = 0   (ok)
            //   seed 1: read 50,  ref 50    → gap to seed 2: read_delta 30, ref_delta 80, gap = -50 (del)
            //   seed 2: read 100, ref 150   → gap to seed 3: read_delta 80, ref_delta 30, gap = +50 (ins)
            //   seed 3: read 200, ref 200   → gap to seed 4: read_delta 30, ref_delta 30, gap = 0   (ok)
            //   seed 4: read 250, ref 250
            //
            // At the two long gaps: n_del = 50, n_ins = 50 → diff = 2*min(50,50) = 100 > 40.
            // Seeds 2 and 3 (indices 2..4 from long_gap at index 2 to long_gap at index 3)
            // should be removed, but the minimap2 algorithm marks seeds from K[st]..K[en],
            // which means seed at index 2 is removed but seed at index 3 is kept.
            let seeds = vec![
                make_hit(0, 0, 0, 20),
                make_hit(0, 50, 50, 20),
                make_hit(0, 150, 100, 20), // misplaced: shifted +100bp on ref
                make_hit(0, 200, 200, 20),
                make_hit(0, 250, 250, 20),
            ];

            let cluster = SeedCluster::new(seeds, false, 1).unwrap();

            // The misplaced seed (ref 150 @ read 100) should have been removed,
            // leaving 4 seeds. The remaining seeds should all be on the good diagonal.
            assert_eq!(cluster.chain.len(), 4);
            assert_eq!(cluster.chain[0].ref_pos, 0);
            assert_eq!(cluster.chain[1].ref_pos, 50);
            assert_eq!(cluster.chain[2].ref_pos, 200);
            assert_eq!(cluster.chain[3].ref_pos, 250);
        }

        #[test]
        fn test_filter_misplaced_seeds_good_chain_unchanged() {
            // A well-behaved chain (all seeds on the same diagonal) should not
            // have any seeds removed.
            let seeds = vec![
                make_hit(0, 100, 0, 20),
                make_hit(0, 150, 50, 20),
                make_hit(0, 200, 100, 20),
                make_hit(0, 250, 150, 20),
                make_hit(0, 300, 200, 20),
            ];

            let cluster = SeedCluster::new(seeds, false, 1).unwrap();
            assert_eq!(cluster.chain.len(), 5);
        }

        // ── Jittery seed filter tests ────────────────────────────────────────────

        #[test]
        fn test_filter_jittery_removes_bouncing_seeds() {
            // Simulate a chain with a stable region followed by jittery seeds
            // from different repeat copies. Diagonal = ref_pos - read_pos.
            //
            // Stable region (diagonal ~1000, spread ~2):
            //   seed 0: read 0,    ref 1000, len 100 → diag 1000
            //   seed 1: read 120,  ref 1121, len 80  → diag 1001
            //   seed 2: read 220,  ref 1220, len 60  → diag 1000
            //
            // Jittery region (diagonal bouncing by 40-80bp, short seeds):
            //   seed 3: read 310,  ref 1350, len 25  → diag 1040 (shift +40)
            //   seed 4: read 360,  ref 1340, len 25  → diag  980 (shift -60)
            //   seed 5: read 410,  ref 1460, len 25  → diag 1050 (shift +70)
            //   seed 6: read 460,  ref 1430, len 25  → diag  970 (shift -80)
            //
            // Stable again:
            //   seed 7: read 520,  ref 1520, len 100 → diag 1000
            //   seed 8: read 640,  ref 1641, len 80  → diag 1001
            //
            // In the jittery zone (seeds 3-6), the shifts are 40+60+70+80 = 250
            // over a ref span of ~1430+25 - 1350 = 105 → density = 250/105 ≈ 2.4
            // which far exceeds the default threshold of 0.15.
            let seeds = vec![
                make_hit(0, 1000, 0, 100),
                make_hit(0, 1121, 120, 80),
                make_hit(0, 1220, 220, 60),
                make_hit(0, 1350, 310, 25),
                make_hit(0, 1340, 360, 25),
                make_hit(0, 1460, 410, 25),
                make_hit(0, 1430, 460, 25),
                make_hit(0, 1520, 520, 100),
                make_hit(0, 1641, 640, 80),
            ];

            let mut chain = seeds.clone();
            chain.sort_by_key(|h| h.read_pos);
            SeedCluster::filter_jittery_seeds(&mut chain, 0.15, 4);

            // The jittery interior seeds (3-6) should be removed.
            // Boundary seeds and stable seeds should remain.
            let remaining_reads: Vec<usize> = chain.iter().map(|s| s.read_pos).collect();
            // Seeds 0,1,2 (stable), 7,8 (stable) should survive.
            // Some boundary seeds from the jittery zone might also survive
            // since we keep the first and last of each window.
            assert!(
                remaining_reads.len() >= 5,
                "should keep at least the 5 stable seeds, got {:?}",
                remaining_reads
            );
            // The core jittery seeds (interior of the window) should be gone.
            // Seeds at read_pos 360 and 410 are always interior to any 4-gap window
            // that covers the jittery region.
            assert!(
                !remaining_reads.contains(&360) || !remaining_reads.contains(&410),
                "at least some jittery seeds should be removed, got {:?}",
                remaining_reads
            );
        }

        #[test]
        fn test_filter_jittery_preserves_stable_chain() {
            // A chain with consistent diagonal (small wobble from real indels)
            // should not lose any seeds.
            //
            // Diagonal ~500, wobble ≤ 5bp.
            let seeds = vec![
                make_hit(0, 500, 0, 50),   // diag 500
                make_hit(0, 553, 50, 40),  // diag 503
                make_hit(0, 600, 98, 30),  // diag 502
                make_hit(0, 635, 130, 35), // diag 505
                make_hit(0, 670, 165, 45), // diag 505
                make_hit(0, 718, 215, 50), // diag 503
                make_hit(0, 770, 268, 40), // diag 502
            ];

            let mut chain = seeds.clone();
            chain.sort_by_key(|h| h.read_pos);
            SeedCluster::filter_jittery_seeds(&mut chain, 0.15, 4);

            assert_eq!(chain.len(), 7, "stable chain should keep all seeds");
        }

        #[test]
        fn test_filter_jittery_single_shift_preserved() {
            // A single large diagonal shift (real indel) shouldn't trigger
            // removal because one shift doesn't make a high density over the
            // full window.
            //
            // Seeds: stable diagonal 500, then shift to 550 (real 50bp deletion),
            // then stable at 550.
            let seeds = vec![
                make_hit(0, 500, 0, 80),    // diag 500
                make_hit(0, 590, 88, 60),   // diag 502
                make_hit(0, 660, 150, 50),  // diag 510 (slight shift)
                make_hit(0, 770, 210, 40),  // diag 560 (big shift — real indel)
                make_hit(0, 870, 310, 80),  // diag 560
                make_hit(0, 1010, 448, 60), // diag 562
                make_hit(0, 1130, 568, 50), // diag 562
            ];

            let mut chain = seeds.clone();
            chain.sort_by_key(|h| h.read_pos);
            SeedCluster::filter_jittery_seeds(&mut chain, 0.15, 4);

            assert_eq!(
                chain.len(),
                7,
                "single diagonal shift (real indel) should not trigger removal"
            );
        }

        #[test]
        fn test_filter_jittery_disabled_at_zero() {
            // When threshold is 0.0, filter should be disabled.
            let seeds = vec![
                make_hit(0, 1000, 0, 25),
                make_hit(0, 1100, 30, 25),  // big shift
                make_hit(0, 1000, 60, 25),  // shift back
                make_hit(0, 1100, 90, 25),  // shift again
                make_hit(0, 1000, 120, 25), // shift back
            ];

            let mut chain = seeds.clone();
            chain.sort_by_key(|h| h.read_pos);
            SeedCluster::filter_jittery_seeds(&mut chain, 0.0, 4);

            assert_eq!(chain.len(), 5, "filter disabled at threshold 0.0");
        }

        #[test]
        fn test_filter_jittery_too_few_seeds() {
            // Chain shorter than window_size+1 should be untouched.
            let seeds = vec![
                make_hit(0, 1000, 0, 25),
                make_hit(0, 1100, 30, 25),
                make_hit(0, 1000, 60, 25),
            ];

            let mut chain = seeds.clone();
            chain.sort_by_key(|h| h.read_pos);
            SeedCluster::filter_jittery_seeds(&mut chain, 0.15, 4);

            assert_eq!(chain.len(), 3, "short chain should be untouched");
        }
    } // mod seed_cluster_tests
}
