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


# ─── Agreement between two fingerprints ─────────────────────────────────────

# Read bases whose projected reference positions lie within this distance of
# each other count as placed identically.  Fingerprints carry no CIGAR, so a
# segment's read-to-reference mapping is taken as a straight line between its
# endpoints, and indels make the true mapping drift from that line.
AGREE_TOLERANCE = 50  # bp


def parse_fingerprint(fp: str) -> list[tuple[str, int, int, str, int, int]]:
    """Parse a fingerprint into (chrom, ref_start, ref_end, strand, q_start, q_end)."""
    segs = []
    if not fp:
        return segs
    for part in fp.split(';'):
        rest, q = part.rsplit(':', 1)
        rest, strand = rest.rsplit(':', 1)
        chrom, ref = rest.rsplit(':', 1)
        ref_start, ref_end = (int(x) for x in ref.split('-'))
        q_start, q_end = (int(x) for x in q.split('-'))
        segs.append((chrom, ref_start, ref_end, strand, q_start, q_end))
    return segs


def _projection(seg) -> tuple[float, float]:
    """Return (slope, intercept) of the segment's read-to-reference line."""
    _, ref_start, ref_end, strand, q_start, q_end = seg
    scale = (ref_end - ref_start) / max(1, q_end - q_start)
    if strand == '+':
        return scale, ref_start - scale * q_start
    # Reverse strand: read coordinates run backwards along the reference.
    return -scale, ref_end + scale * q_start


def _union_length(intervals: list[tuple[float, float]]) -> float:
    total = 0.0
    end = float('-inf')
    for a, b in sorted(intervals):
        if b <= end:
            continue
        total += b - max(a, end)
        end = b
    return total


def alignment_agreement(fp_a: str, fp_b: str, read_len: int) -> dict:
    """Compare two fingerprints of one read, base by base along the read.

    Returns bases placed identically by both alignments ('agree'), and bases
    aligned by at least one but not placed identically ('differ'), each in bp
    and as a percentage of read_len.
    """
    segs_a = parse_fingerprint(fp_a)
    segs_b = parse_fingerprint(fp_b)
    aligned = [(s[4], s[5]) for s in segs_a + segs_b]
    agreeing = []
    for sa in segs_a:
        for sb in segs_b:
            if sa[0] != sb[0] or sa[3] != sb[3]:
                continue
            lo, hi = max(sa[4], sb[4]), min(sa[5], sb[5])
            if lo >= hi:
                continue
            # The projections' difference is linear in the read position, so
            # the positions where it is within tolerance form an interval.
            ma, ca = _projection(sa)
            mb, cb = _projection(sb)
            dm, dc = ma - mb, ca - cb
            if abs(dm) < 1e-12:
                if abs(dc) <= AGREE_TOLERANCE:
                    agreeing.append((lo, hi))
                continue
            x1 = (-AGREE_TOLERANCE - dc) / dm
            x2 = (AGREE_TOLERANCE - dc) / dm
            a, b = max(lo, min(x1, x2)), min(hi, max(x1, x2))
            if a < b:
                agreeing.append((a, b))
    agree = round(_union_length(agreeing))
    differ = round(_union_length(aligned)) - agree
    denom = max(1, read_len)
    return {
        'read_len': read_len,
        'agree_bp': agree,
        'agree_pct': 100 * agree / denom,
        'differ_bp': differ,
        'differ_pct': 100 * differ / denom,
    }
