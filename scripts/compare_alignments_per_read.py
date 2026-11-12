#!/usr/bin/env python3
"""
Per-read comparison of two name-sorted BAM files (e.g. parallax vs minimap2).

For each read produces one TSV row with 17 columns:
  read_id, read_length,
  plx_n_primary, plx_n_secondary,
  plx_primary_frac, plx_secondary_frac, plx_mapq, plx_chrom,
  mm2_n_primary, mm2_n_secondary,
  mm2_primary_frac, mm2_secondary_frac, mm2_mapq, mm2_chrom,
  plx_primary_vs_mm2_primary_frac,
  mm2_primary_vs_plx_primary_frac,
  plx_primary_vs_mm2_secondary_frac,
  mm2_primary_vs_plx_secondary_frac

Overlap fractions are strand-aware: intervals only overlap if on the same strand.
Denominators for overlap fractions are the total non-secondary (or secondary) mapped
length of the query read in the numerator's aligner.

Both BAMs must be sorted by read name (samtools sort -n).
"""
from __future__ import annotations
import argparse
import sys
from itertools import groupby
from typing import NamedTuple

try:
    import pysam
except ImportError:
    sys.exit("pysam is required: pip install pysam")


# ─── Data types ──────────────────────────────────────────────────────────────

class Seg(NamedTuple):
    """A single aligned segment on the query read."""
    q_start: int   # 0-based, query coordinates
    q_end:   int   # exclusive
    strand:  str   # '+' or '-'
    chrom:   str


class ReadAlns:
    """All alignment records for one read from one BAM."""

    def __init__(self, read_length: int):
        self.read_length = read_length
        self.primary_segs:   list[Seg] = []   # non-secondary
        self.secondary_segs: list[Seg] = []   # secondary
        self.mapq = 0                          # max MAPQ across non-secondary
        self.chrom: str | None = None          # chrom of FLAG-primary record

    def add(self, rec: pysam.AlignedSegment) -> None:
        if rec.is_unmapped:
            return
        strand = '-' if rec.is_reverse else '+'
        qs, qe = query_span(rec)
        if rec.is_secondary:
            self.secondary_segs.append(Seg(qs, qe, strand, rec.reference_name))
        else:
            self.primary_segs.append(Seg(qs, qe, strand, rec.reference_name))
            if rec.mapping_quality is not None:
                self.mapq = max(self.mapq, rec.mapping_quality)
            if not rec.is_supplementary:
                self.chrom = rec.reference_name


# ─── Interval helpers ────────────────────────────────────────────────────────

def query_span(rec: pysam.AlignedSegment) -> tuple[int, int]:
    """Query (read) coordinates consumed by this alignment, 0-based exclusive."""
    qs = rec.query_alignment_start
    qe = rec.query_alignment_end
    return qs, qe


def merge_intervals(segs: list[Seg]) -> list[tuple[int, int, str, str]]:
    """Merge overlapping segments sharing the same chrom+strand. Returns (start, end, strand, chrom)."""
    by_key: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for s in segs:
        by_key.setdefault((s.chrom, s.strand), []).append((s.q_start, s.q_end))
    merged = []
    for (chrom, strand), ivs in by_key.items():
        ivs.sort()
        cur_start, cur_end = ivs[0]
        for start, end in ivs[1:]:
            if start <= cur_end:
                cur_end = max(cur_end, end)
            else:
                merged.append((cur_start, cur_end, strand, chrom))
                cur_start, cur_end = start, end
        merged.append((cur_start, cur_end, strand, chrom))
    return merged


def total_length(merged: list[tuple[int, int, str, str]]) -> int:
    return sum(e - s for s, e, _, _ in merged)


def intersect_length(a: list[tuple[int, int, str, str]], b: list[tuple[int, int, str, str]]) -> int:
    """Total length of chrom+strand-aware intersection of two merged interval lists."""
    from collections import defaultdict
    # index b by (chrom, strand)
    b_idx: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    for s, e, strand, chrom in b:
        b_idx[(chrom, strand)].append((s, e))
    for ivs in b_idx.values():
        ivs.sort()

    total = 0
    for s, e, strand, chrom in a:
        for bs, be in b_idx.get((chrom, strand), []):
            lo = max(s, bs)
            hi = min(e, be)
            if lo < hi:
                total += hi - lo
    return total


# ─── Per-read stats ──────────────────────────────────────────────────────────

def compute_stats(a: ReadAlns, b: ReadAlns) -> dict:
    rl = a.read_length or b.read_length or 1

    a_pri = merge_intervals(a.primary_segs)
    a_sec = merge_intervals(a.secondary_segs)
    b_pri = merge_intervals(b.primary_segs)
    b_sec = merge_intervals(b.secondary_segs)

    a_pri_len = total_length(a_pri)
    a_sec_len = total_length(a_sec)
    b_pri_len = total_length(b_pri)
    b_sec_len = total_length(b_sec)

    def safe_frac(num: int, denom: int) -> float:
        return num / denom if denom > 0 else 0.0

    return {
        "plx_n_primary":   len(a.primary_segs),
        "plx_n_secondary": len(a.secondary_segs),
        "plx_primary_frac":   safe_frac(a_pri_len, rl),
        "plx_secondary_frac": safe_frac(a_sec_len, rl),
        "plx_mapq":  a.mapq,
        "plx_chrom": a.chrom or ".",
        "mm2_n_primary":   len(b.primary_segs),
        "mm2_n_secondary": len(b.secondary_segs),
        "mm2_primary_frac":   safe_frac(b_pri_len, rl),
        "mm2_secondary_frac": safe_frac(b_sec_len, rl),
        "mm2_mapq":  b.mapq,
        "mm2_chrom": b.chrom or ".",
        "plx_primary_vs_mm2_primary_frac":    safe_frac(intersect_length(a_pri, b_pri), a_pri_len),
        "mm2_primary_vs_plx_primary_frac":    safe_frac(intersect_length(b_pri, a_pri), b_pri_len),
        "plx_primary_vs_mm2_secondary_frac":  safe_frac(intersect_length(a_pri, b_sec), a_pri_len),
        "mm2_primary_vs_plx_secondary_frac":  safe_frac(intersect_length(b_pri, a_sec), b_pri_len),
    }


COLUMNS = [
    "read_id", "read_length",
    "plx_n_primary", "plx_n_secondary",
    "plx_primary_frac", "plx_secondary_frac", "plx_mapq", "plx_chrom",
    "mm2_n_primary", "mm2_n_secondary",
    "mm2_primary_frac", "mm2_secondary_frac", "mm2_mapq", "mm2_chrom",
    "plx_primary_vs_mm2_primary_frac",
    "mm2_primary_vs_plx_primary_frac",
    "plx_primary_vs_mm2_secondary_frac",
    "mm2_primary_vs_plx_secondary_frac",
]


# ─── BAM reading ─────────────────────────────────────────────────────────────

def iter_reads(bam_path: str):
    """Yield (read_name, ReadAlns) for each read in a name-sorted BAM."""
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for name, records in groupby(bam.fetch(until_eof=True), key=lambda r: r.query_name):
            records = list(records)
            rl = next((r.query_length for r in records if r.query_length), 0)
            ra = ReadAlns(rl)
            for rec in records:
                ra.add(rec)
            yield name, ra


def empty_read(read_length: int = 0) -> ReadAlns:
    return ReadAlns(read_length)


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plx_bam", help="parallax BAM (name-sorted)")
    ap.add_argument("mm2_bam", help="minimap2 BAM (name-sorted)")
    ap.add_argument("-o", "--output", default="-", help="output TSV (default: stdout)")
    args = ap.parse_args()

    out = open(args.output, "w") if args.output != "-" else sys.stdout

    out.write("\t".join(COLUMNS) + "\n")

    plx_iter = iter_reads(args.plx_bam)
    mm2_iter = iter_reads(args.mm2_bam)

    plx_cur = next(plx_iter, None)
    mm2_cur = next(mm2_iter, None)

    n_written = 0
    while plx_cur is not None or mm2_cur is not None:
        plx_name = plx_cur[0] if plx_cur else None
        mm2_name = mm2_cur[0] if mm2_cur else None

        if plx_name is not None and (mm2_name is None or plx_name < mm2_name):
            name, plx_ra = plx_cur
            mm2_ra = empty_read(plx_cur[1].read_length)
            plx_cur = next(plx_iter, None)
        elif mm2_name is not None and (plx_name is None or mm2_name < plx_name):
            name, mm2_ra = mm2_cur
            plx_ra = empty_read(mm2_cur[1].read_length)
            mm2_cur = next(mm2_iter, None)
        else:
            name, plx_ra = plx_cur
            _,    mm2_ra = mm2_cur
            plx_cur = next(plx_iter, None)
            mm2_cur = next(mm2_iter, None)

        rl = plx_ra.read_length or mm2_ra.read_length
        stats = compute_stats(plx_ra, mm2_ra)

        row = [name, str(rl)] + [
            f"{stats[c]:.4f}" if isinstance(stats[c], float) else str(stats[c])
            for c in COLUMNS[2:]
        ]
        out.write("\t".join(row) + "\n")
        n_written += 1

    if out is not sys.stdout:
        out.close()

    print(f"Wrote {n_written} rows", file=sys.stderr)


if __name__ == "__main__":
    main()
