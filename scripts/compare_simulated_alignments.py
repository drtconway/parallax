#!/usr/bin/env python3
"""
Compare two aligners' results on simulated reads.

Takes two SAM files (e.g. parallax output and minimap2 output) produced from
the same set of simulated reads, and reports reads where the outcome differs.

Read names must encode the true origin (same format as check_simulated_alignments.py):
    sim_XXXXXXXX:chrom_start_end_strand[:errors]

For each read, the outcome is classified as:
    OK          - primary alignment matches expected position
    MISMATCH    - primary alignment at wrong position
    SECONDARY   - primary wrong, but a secondary/supplementary matches
    UNMAPPED    - no primary alignment

The script prints a confusion matrix (summary) and optionally the per-read
detail for reads where the two aligners disagree.

Example:
    python compare_simulated_alignments.py parallax.sam minimap2.sam
    python compare_simulated_alignments.py -v parallax.sam minimap2.sam
"""

import argparse
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from typing import Optional, List


# ── Data types (shared with check_simulated_alignments.py) ──────────────


@dataclass
class ExpectedAlignment:
    chrom: str
    start: int
    end: int
    strand: str
    errors: Optional[str] = None


@dataclass
class SamAlignment:
    qname: str
    flag: int
    rname: str
    pos: int
    mapq: int
    cigar: str
    seq_len: int

    @property
    def is_unmapped(self) -> bool:
        return (self.flag & 0x4) != 0

    @property
    def is_secondary(self) -> bool:
        return (self.flag & 0x100) != 0

    @property
    def is_supplementary(self) -> bool:
        return (self.flag & 0x800) != 0

    @property
    def is_reverse(self) -> bool:
        return (self.flag & 0x10) != 0

    @property
    def is_primary(self) -> bool:
        return not self.is_secondary and not self.is_supplementary and not self.is_unmapped


class Outcome(Enum):
    OK = "OK"
    MISMATCH = "MISMATCH"
    SECONDARY = "SECONDARY"
    UNMAPPED = "UNMAPPED"


# ── Parsing helpers ─────────────────────────────────────────────────────


def parse_read_name(qname: str) -> Optional[ExpectedAlignment]:
    parts = qname.split(":")
    if len(parts) < 2:
        return None
    loc_part = parts[1]
    errors = ":".join(parts[2:]) if len(parts) > 2 else None
    match = re.match(r"^(.+)_(\d+)_(\d+)_([+-])$", loc_part)
    if not match:
        return None
    return ExpectedAlignment(
        match.group(1), int(match.group(2)), int(match.group(3)), match.group(4), errors
    )


def parse_cigar_ref_length(cigar: str) -> int:
    length = 0
    for m in re.finditer(r"(\d+)([MIDNSHP=X])", cigar):
        if m.group(2) in "MDN=X":
            length += int(m.group(1))
    return length


def check_alignment(
    expected: ExpectedAlignment, actual: SamAlignment, tolerance: int
) -> bool:
    if expected.chrom != actual.rname:
        return False
    if (expected.strand == "-") != actual.is_reverse:
        return False
    ref_len = parse_cigar_ref_length(actual.cigar)
    actual_end = actual.pos + ref_len - 1
    if abs(actual.pos - expected.start) > tolerance:
        return False
    if abs(actual_end - expected.end) > tolerance:
        return False
    return True


# ── SAM loading ─────────────────────────────────────────────────────────


def load_sam(path: str) -> dict[str, List[SamAlignment]]:
    by_read: dict[str, List[SamAlignment]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            by_read[fields[0]].append(
                SamAlignment(
                    qname=fields[0],
                    flag=int(fields[1]),
                    rname=fields[2],
                    pos=int(fields[3]),
                    mapq=int(fields[4]),
                    cigar=fields[5],
                    seq_len=len(fields[9]),
                )
            )
    return by_read


# ── Classify a single read ──────────────────────────────────────────────


@dataclass
class ReadResult:
    outcome: Outcome
    primary: Optional[SamAlignment] = None
    reason: str = ""


def classify_read(
    expected: ExpectedAlignment, alns: List[SamAlignment], tolerance: int
) -> ReadResult:
    primary = None
    secondaries: List[SamAlignment] = []
    for a in alns:
        if a.is_primary:
            primary = a
        elif not a.is_unmapped:
            secondaries.append(a)

    if primary is None or primary.is_unmapped:
        return ReadResult(Outcome.UNMAPPED)

    if check_alignment(expected, primary, tolerance):
        return ReadResult(Outcome.OK, primary)

    for sec in secondaries:
        if check_alignment(expected, sec, tolerance):
            return ReadResult(Outcome.SECONDARY, primary, reason=f"primary={primary.rname}:{primary.pos}")
            
    ref_len = parse_cigar_ref_length(primary.cigar)
    actual_end = primary.pos + ref_len - 1
    reason = (
        f"expected {expected.chrom}:{expected.start}-{expected.end}{expected.strand}, "
        f"got {primary.rname}:{primary.pos}-{actual_end}{'-' if primary.is_reverse else '+'} MAPQ={primary.mapq}"
    )
    return ReadResult(Outcome.MISMATCH, primary, reason=reason)


# ── Formatting helpers ──────────────────────────────────────────────────


def overlap_pct(expected: ExpectedAlignment, actual: SamAlignment) -> float:
    """Return the percentage of the actual mapped span that overlaps the expected locus.

    Returns 0.0 when the alignment is on a different chromosome or unmapped.
    """
    if actual.is_unmapped or actual.rname != expected.chrom:
        return 0.0
    ref_len = parse_cigar_ref_length(actual.cigar)
    if ref_len == 0:
        return 0.0
    actual_start = actual.pos
    actual_end = actual.pos + ref_len  # half-open
    overlap = max(0, min(expected.end, actual_end) - max(expected.start, actual_start))
    expected_len = expected.end - expected.start
    if expected_len == 0:
        return 0.0
    return 100.0 * overlap / expected_len


def format_result_columns(res: ReadResult, expected: Optional[ExpectedAlignment] = None) -> tuple[str, str, str, str, str, str, str]:
    """Return (outcome, chrom, start, end, strand, mapq, overlap_pct) as strings."""
    if res.primary is None:
        return (res.outcome.value, "", "", "", "", "", "")
    ref_len = parse_cigar_ref_length(res.primary.cigar)
    end = res.primary.pos + ref_len - 1
    strand = "-" if res.primary.is_reverse else "+"
    pct = f"{overlap_pct(expected, res.primary):.1f}" if expected is not None else ""
    return (
        res.outcome.value,
        res.primary.rname,
        str(res.primary.pos),
        str(end),
        strand,
        str(res.primary.mapq),
        pct,
    )


# ── Main ────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="Compare two aligners on simulated reads"
    )
    parser.add_argument("sam_a", help="First SAM file (e.g. parallax)")
    parser.add_argument("sam_b", help="Second SAM file (e.g. minimap2)")
    parser.add_argument("--name-a", default="A", help="Label for first aligner (default: A)")
    parser.add_argument("--name-b", default="B", help="Label for second aligner (default: B)")
    parser.add_argument(
        "-t", "--tolerance", type=int, default=50,
        help="Position tolerance in bp (default: 50)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Print per-read detail for every disagreement",
    )
    parser.add_argument(
        "--only", choices=["a-better", "b-better", "both-wrong"],
        help="Only show a specific category of disagreement",
    )
    args = parser.parse_args()

    print(f"Loading {args.name_a}: {args.sam_a}", file=sys.stderr)
    reads_a = load_sam(args.sam_a)
    print(f"Loading {args.name_b}: {args.sam_b}", file=sys.stderr)
    reads_b = load_sam(args.sam_b)

    all_reads = sorted(set(reads_a.keys()) | set(reads_b.keys()))

    # Confusion matrix: (outcome_a, outcome_b) -> count
    matrix: dict[tuple[Outcome, Outcome], int] = defaultdict(int)
    skipped = 0
    diffs: list[tuple[str, ReadResult, ReadResult, ExpectedAlignment]] = []

    for qname in all_reads:
        expected = parse_read_name(qname)
        if expected is None:
            skipped += 1
            continue

        alns_a = reads_a.get(qname, [])
        alns_b = reads_b.get(qname, [])

        res_a = classify_read(expected, alns_a, args.tolerance)
        res_b = classify_read(expected, alns_b, args.tolerance)

        matrix[(res_a.outcome, res_b.outcome)] += 1

        if res_a.outcome != res_b.outcome:
            diffs.append((qname, res_a, res_b, expected))

    # ── Confusion matrix ────────────────────────────────────────────
    outcomes = [Outcome.OK, Outcome.MISMATCH, Outcome.SECONDARY, Outcome.UNMAPPED]
    total = sum(matrix.values())

    # Column widths
    label_w = max(len(args.name_a), len(args.name_b), 10) + 2
    col_w = 10

    print(f"\n{'':>{label_w}} | ", end="", file=sys.stderr)
    for o in outcomes:
        print(f"{args.name_b + ':' + o.value:>{col_w}}", end="  ", file=sys.stderr)
    print(f"{'total':>{col_w}}", file=sys.stderr)
    print("-" * (label_w + 3 + (col_w + 2) * (len(outcomes) + 1)), file=sys.stderr)

    for oa in outcomes:
        row_total = sum(matrix[(oa, ob)] for ob in outcomes)
        print(f"{args.name_a + ':' + oa.value:>{label_w}} | ", end="", file=sys.stderr)
        for ob in outcomes:
            c = matrix[(oa, ob)]
            cell = str(c) if c > 0 else "."
            print(f"{cell:>{col_w}}", end="  ", file=sys.stderr)
        print(f"{row_total:>{col_w}}", file=sys.stderr)

    col_totals = [sum(matrix[(oa, ob)] for oa in outcomes) for ob in outcomes]
    print(f"{'total':>{label_w}} | ", end="", file=sys.stderr)
    for ct in col_totals:
        print(f"{ct:>{col_w}}", end="  ", file=sys.stderr)
    print(f"{total:>{col_w}}", file=sys.stderr)

    # ── Summary counts ──────────────────────────────────────────────
    a_ok = sum(matrix[(Outcome.OK, ob)] for ob in outcomes)
    b_ok = sum(matrix[(oa, Outcome.OK)] for oa in outcomes)
    both_ok = matrix[(Outcome.OK, Outcome.OK)]
    a_better = sum(
        matrix[(Outcome.OK, ob)]
        for ob in outcomes
        if ob != Outcome.OK
    )
    b_better = sum(
        matrix[(oa, Outcome.OK)]
        for oa in outcomes
        if oa != Outcome.OK
    )
    both_wrong = total - a_ok - b_ok + both_ok  # inclusion-exclusion

    print(f"\n--- Summary ---", file=sys.stderr)
    print(f"Total reads:  {total}", file=sys.stderr)
    print(f"{args.name_a} correct: {a_ok} ({100*a_ok/total:.2f}%)" if total else f"{args.name_a} correct: 0", file=sys.stderr)
    print(f"{args.name_b} correct: {b_ok} ({100*b_ok/total:.2f}%)" if total else f"{args.name_b} correct: 0", file=sys.stderr)
    print(f"Both correct: {both_ok} ({100*both_ok/total:.2f}%)" if total else "Both correct: 0", file=sys.stderr)
    print(f"Only {args.name_a} correct: {a_better}", file=sys.stderr)
    print(f"Only {args.name_b} correct: {b_better}", file=sys.stderr)
    print(f"Both wrong:   {both_wrong}", file=sys.stderr)
    print(f"Disagree:     {len(diffs)}", file=sys.stderr)
    if skipped:
        print(f"Skipped (unparseable): {skipped}", file=sys.stderr)

    # ── Per-read detail ─────────────────────────────────────────────
    if args.verbose and diffs:
        # TSV header with separate columns for each locus
        header = "\t".join([
            "read",
            "expected_chrom", "expected_start", "expected_end", "expected_strand",
            f"{args.name_a}_outcome", f"{args.name_a}_chrom", f"{args.name_a}_start",
            f"{args.name_a}_end", f"{args.name_a}_strand", f"{args.name_a}_mapq",
            f"{args.name_a}_overlap_pct",
            f"{args.name_b}_outcome", f"{args.name_b}_chrom", f"{args.name_b}_start",
            f"{args.name_b}_end", f"{args.name_b}_strand", f"{args.name_b}_mapq",
            f"{args.name_b}_overlap_pct",
        ])
        print(header)
        for qname, res_a, res_b, exp in diffs:
            if args.only == "a-better" and res_a.outcome != Outcome.OK:
                continue
            if args.only == "b-better" and res_b.outcome != Outcome.OK:
                continue
            if args.only == "both-wrong" and (res_a.outcome == Outcome.OK or res_b.outcome == Outcome.OK):
                continue
            cols_a = format_result_columns(res_a, exp)
            cols_b = format_result_columns(res_b, exp)
            row = "\t".join([
                qname,
                exp.chrom, str(exp.start), str(exp.end), exp.strand,
                *cols_a,
                *cols_b,
            ])
            print(row)


if __name__ == "__main__":
    main()
