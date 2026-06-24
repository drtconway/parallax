"""Alignment fingerprinting utilities.

A fingerprint is a stable string representation of all primary+supplementary
alignments for a single read, suitable for detecting whether parallax has
changed its alignment between runs.

Format per segment:  chrom:ref_start-ref_end:strand:query_start-query_end
Multiple segments are sorted and joined with ';' so the fingerprint is
independent of the order records appear in the BAM.

Query coordinates are derived from the CIGAR string (not pysam properties)
using the same logic as compare_alignments_2.py, since pysam's query_alignment_start
/ query_alignment_end do not account for hard-clipped supplementary records
correctly for our purposes.
"""

from pysam import AlignedSegment, AlignmentFile


def query_length(rec: AlignedSegment) -> int:
    total = 0
    for (op, length) in rec.cigartuples:
        if op != 2:  # everything except deletion contributes to query length
            total += length
    return total


def query_range(rec: AlignedSegment) -> tuple[int, int]:
    cigar = rec.cigartuples
    lead  = cigar[0][1]  if cigar[0][0]  in (4, 5) else 0
    trail = cigar[-1][1] if cigar[-1][0] in (4, 5) else 0
    length = query_length(rec)
    q_start, q_end = lead, length - trail
    if rec.is_reverse:
        q_start, q_end = length - q_end, length - q_start
    return (q_start, q_end)


def segment_fingerprint(rec: AlignedSegment) -> str:
    chrom = rec.reference_name
    ref_start = rec.reference_start
    ref_end = rec.reference_end
    strand = '-' if rec.is_reverse else '+'
    q_start, q_end = query_range(rec)
    return f"{chrom}:{ref_start}-{ref_end}:{strand}:{q_start}-{q_end}"


def read_fingerprint(records: list[AlignedSegment]) -> str:
    """Return a stable fingerprint for all aligned segments of one read.

    Skips unmapped, secondary records. Includes primary and supplementary.
    Segments are sorted so fingerprint is independent of BAM record order.
    An unmapped read returns the empty string.
    """
    segments = []
    for rec in records:
        if rec.is_unmapped or rec.is_secondary:
            continue
        segments.append(segment_fingerprint(rec))
    if not segments:
        return ''
    segments.sort()
    return ';'.join(segments)


def fingerprints_from_bam(bam_path: str) -> dict[str, str]:
    """Return {read_name: fingerprint} for every read in a name-sorted BAM."""
    result = {}
    with AlignmentFile(bam_path, 'rb') as bam:
        current_name = None
        current_records: list[AlignedSegment] = []
        for rec in bam:
            if rec.is_secondary:
                continue
            if rec.query_name != current_name:
                if current_name is not None:
                    result[current_name] = read_fingerprint(current_records)
                current_name = rec.query_name
                current_records = [rec]
            else:
                current_records.append(rec)
        if current_name is not None:
            result[current_name] = read_fingerprint(current_records)
    return result
