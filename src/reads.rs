
use std::{sync::Arc, usize};

use parallax::{
    error::Result,
    index::Index,
    reference::InMemoryReference,
    utils::{debug, sequence::reverse_complement_into},
};
use crate::writer::{AlignmentWriter, OutputFormat};
use crate::aligner::{Aligner, AlignerBuilder};
use crate::explanatory;

pub mod builder;
pub mod extended;
pub mod seeds;

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
    no_secondary: bool,
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
                    explanatory::ExplanatoryAlignerBuilder::new(&reference, index, &writer)
                        .no_secondary(no_secondary)
                        .build();
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
}
