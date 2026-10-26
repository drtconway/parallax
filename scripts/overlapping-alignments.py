#!/usr/bin/env python3
"""Print the alignments of reads that have two or more overlapping records.

The input must have all the records for a read adjacent (e.g. aligner output,
or `samtools sort -n`).  For each read whose primary/supplementary records
include a pair overlapping by at least --min-overlap bases, every record of
that read is written as SAM (with header) to stdout.

Overlap is measured on the reference (same chromosome, either strand) by
default, or in forward-strand read coordinates with --space query.

    overlapping-alignments.py in.bam [--space ref|query] [--min-overlap N]
"""

import argparse
import itertools
import sys

import pysam

CLIPS = (pysam.CSOFT_CLIP, pysam.CHARD_CLIP)


def query_range(rec):
    """Forward-strand read interval [qs, qe) covered by the alignment."""
    cigar = rec.cigartuples
    lead = sum(n for op, n in itertools.takewhile(lambda c: c[0] in CLIPS, cigar))
    trail = sum(n for op, n in itertools.takewhile(lambda c: c[0] in CLIPS, reversed(cigar)))
    length = rec.infer_read_length()
    qs, qe = lead, length - trail
    if rec.is_reverse:
        qs, qe = length - qe, length - qs
    return qs, qe


def interval(rec, space):
    if space == "query":
        return (None,) + query_range(rec)
    return rec.reference_name, rec.reference_start, rec.reference_end


def has_overlap(records, space, min_overlap):
    ivs = sorted(interval(r, space) for r in records
                 if not (r.is_unmapped or r.is_secondary) and r.cigartuples)
    # Sorted by (chrom, start), each interval's largest overlap with an earlier
    # one is with the one reaching furthest right.
    furthest = {}
    for chrom, start, end in ivs:
        if chrom in furthest and min(furthest[chrom], end) - start >= min_overlap:
            return True
        furthest[chrom] = max(furthest.get(chrom, end), end)
    return False


def main():
    p = argparse.ArgumentParser(
        description="Print the alignments of reads that have two or more overlapping records.")
    p.add_argument("bam", help="SAM/BAM/CRAM with each read's records adjacent")
    p.add_argument("--space", choices=("ref", "query"), default="ref",
                   help="overlap on the reference or in read coordinates (default: ref)")
    p.add_argument("--min-overlap", type=int, default=1,
                   help="minimum overlap in bases (default: 1)")
    args = p.parse_args()

    reads = flagged = 0
    with pysam.AlignmentFile(args.bam) as f, \
            pysam.AlignmentFile("-", "w", template=f) as out:
        for _, group in itertools.groupby(f.fetch(until_eof=True), key=lambda r: r.query_name):
            records = list(group)
            reads += 1
            if has_overlap(records, args.space, args.min_overlap):
                flagged += 1
                for r in records:
                    out.write(r)
    print(f"{flagged} of {reads} reads have overlapping alignments", file=sys.stderr)


if __name__ == "__main__":
    main()
