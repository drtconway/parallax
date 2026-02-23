#!/usr/bin/env python3
"""
Convert a UCSC RepeatMasker TSV export to a BED file, filtering out
highly truncated repeat instances.

The coverage of each instance relative to the full repeat consensus is:

    coverage = min(1.0, (repEnd - repStart + 1) / (repEnd + |repLeft|))

where the denominator estimates the total consensus length from the right-hand
coordinates (repEnd: last matched position; repLeft: bases remaining after
match, stored negative by convention but handled via abs()).  Instances where
repStart <= 0 extend beyond the left edge of the consensus and are treated as
having full left-side coverage.

Filtering is applied per (repClass, repFamily, repName) group using a
configurable quantile cutoff (default: 75th percentile within each group),
so that the threshold adapts to the size distribution of each element family
rather than applying a single global number across elements ranging from
~300 bp (Alu) to ~6 kb (full-length L1).

Output BED columns:
    chrom, chromStart, chromEnd, name, score, strand

where name is "repName|repClass|repFamily" and score is the swScore.

Usage:
    python ucsc_repeats_to_bed.py -i repeats.tsv.gz -o repeats.bed
    python ucsc_repeats_to_bed.py -i repeats.tsv.gz -o repeats.bed --quantile 0.5
    python ucsc_repeats_to_bed.py -i repeats.tsv.gz -o repeats.bed --min-coverage 0.5
    python ucsc_repeats_to_bed.py -i repeats.tsv.gz -o repeats.bed --stats
"""
from __future__ import annotations

import argparse
import gzip
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class RepeatRecord:
    # Genomic coordinates (genoStart is 0-based as exported by UCSC)
    chrom: str
    start: int
    end: int
    strand: str

    # Alignment quality
    sw_score: int

    # Repeat annotation
    rep_name: str
    rep_class: str
    rep_family: str

    # Consensus coordinates
    rep_start: int   # 1-based; can be <= 0 when instance overhangs consensus left edge
    rep_end: int     # 1-based last matched position in consensus
    rep_left: int    # raw value from TSV; represents -#bases remaining (often negative)

    @property
    def group_key(self) -> tuple[str, str, str]:
        return (self.rep_class, self.rep_family, self.rep_name)

    @property
    def coverage(self) -> float:
        """Fraction of the repeat consensus covered by this instance.

        Consensus length is estimated as repEnd + |repLeft|.  When repStart <= 0
        the instance overhangs the left edge, so covered bases are clamped to
        [1, repEnd] and coverage is capped at 1.0.
        """
        consensus_length = self.rep_end + abs(self.rep_left)
        if consensus_length <= 0:
            return 0.0
        covered = self.rep_end - max(self.rep_start, 1) + 1
        return min(1.0, covered / consensus_length)

    @property
    def bed_name(self) -> str:
        return f"{self.rep_name}|{self.rep_class}|{self.rep_family}"

    def to_bed_line(self) -> str:
        return (
            f"{self.chrom}\t{self.start}\t{self.end}\t"
            f"{self.bed_name}\t{self.sw_score}\t{self.strand}\n"
        )


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def open_input(path: Path):
    """Open a plain or gzip-compressed file for text reading."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_tsv(path: Path) -> Iterator[RepeatRecord]:
    """Parse a UCSC RepeatMasker TSV export, yielding one RepeatRecord per row.

    Expected column order (0-based indices):
        0  #bin
        1  swScore
        2  milliDiv
        3  milliDel
        4  milliIns
        5  genoName
        6  genoStart
        7  genoEnd
        8  genoLeft
        9  strand
        10 repName
        11 repClass
        12 repFamily
        13 repStart
        14 repEnd
        15 repLeft
        16 id
    """
    with open_input(path) as fh:
        saw_header = False
        for lineno, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                saw_header = True
                continue
            if not saw_header:
                # Treat the first non-comment line as a header if we never saw '#'
                saw_header = True
                continue

            fields = line.split("\t")
            if len(fields) < 17:
                print(
                    f"Warning: line {lineno} has only {len(fields)} fields, skipping",
                    file=sys.stderr,
                )
                continue
            try:
                yield RepeatRecord(
                    chrom=fields[5],
                    start=int(fields[6]),
                    end=int(fields[7]),
                    strand=fields[9],
                    sw_score=int(fields[1]),
                    rep_name=fields[10],
                    rep_class=fields[11],
                    rep_family=fields[12],
                    rep_start=int(fields[13]),
                    rep_end=int(fields[14]),
                    rep_left=int(fields[15]),
                )
            except ValueError as exc:
                print(
                    f"Warning: line {lineno}: could not parse fields ({exc}), skipping",
                    file=sys.stderr,
                )


# ---------------------------------------------------------------------------
# Statistics / filtering
# ---------------------------------------------------------------------------

def _quantile(values: list[float], q: float) -> float:
    """Linear-interpolation quantile of a list (no numpy required)."""
    if not values:
        return 0.0
    sv = sorted(values)
    pos = q * (len(sv) - 1)
    lo = int(pos)
    hi = lo + 1
    if hi >= len(sv):
        return sv[-1]
    return sv[lo] + (pos - lo) * (sv[hi] - sv[lo])


def compute_group_cutoffs(
    records: list[RepeatRecord],
    q: float,
) -> dict[tuple[str, str, str], float]:
    """Compute per-group coverage quantile cutoffs."""
    group_coverages: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for rec in records:
        group_coverages[rec.group_key].append(rec.coverage)
    return {
        group: _quantile(covs, q)
        for group, covs in group_coverages.items()
    }


def print_stats(
    records: list[RepeatRecord],
    kept: list[RepeatRecord],
    cutoffs: dict[tuple[str, str, str], float],
) -> None:
    group_total: dict[tuple, int] = defaultdict(int)
    group_kept: dict[tuple, int] = defaultdict(int)
    for rec in records:
        group_total[rec.group_key] += 1
    for rec in kept:
        group_kept[rec.group_key] += 1

    header = f"{'repClass':<20} {'repFamily':<20} {'repName':<30} {'cutoff':>8} {'total':>8} {'kept':>8} {'%kept':>7}"
    print(header, file=sys.stderr)
    print("-" * len(header), file=sys.stderr)
    for group in sorted(group_total):
        rc, rf, rn = group
        total = group_total[group]
        k = group_kept.get(group, 0)
        pct = 100.0 * k / total if total else 0.0
        print(
            f"{rc:<20} {rf:<20} {rn:<30} "
            f"{cutoffs[group]:>8.3f} {total:>8,} {k:>8,} {pct:>6.1f}%",
            file=sys.stderr,
        )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        type=Path,
        metavar="TSV",
        help="Input UCSC RepeatMasker TSV (plain or .gz)",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=None,
        metavar="BED",
        help="Output BED file (default: stdout)",
    )
    parser.add_argument(
        "-q", "--quantile",
        type=float,
        default=0.75,
        metavar="Q",
        help=(
            "Per-group coverage quantile below which instances are dropped "
            "(default: 0.75, i.e. bottom 75%% of coverage within each group)"
        ),
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=None,
        metavar="COV",
        help=(
            "Override: apply a fixed minimum coverage fraction instead of a "
            "per-group quantile (0.0–1.0)"
        ),
    )
    parser.add_argument(
        "--stats",
        action="store_true",
        help="Print a per-group filtering summary table to stderr",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    print(f"Reading {args.input} ...", file=sys.stderr)
    records = list(parse_tsv(args.input))
    print(f"  {len(records):,} records loaded", file=sys.stderr)

    # Compute cutoffs
    if args.min_coverage is not None:
        cutoffs: dict[tuple[str, str, str], float] = {
            rec.group_key: args.min_coverage for rec in records
        }
        print(
            f"  Using fixed minimum coverage {args.min_coverage}",
            file=sys.stderr,
        )
    else:
        cutoffs = compute_group_cutoffs(records, args.quantile)
        print(
            f"  Per-group coverage cutoffs computed at q={args.quantile} "
            f"across {len(cutoffs):,} groups",
            file=sys.stderr,
        )

    # Filter
    kept = [rec for rec in records if rec.coverage >= cutoffs[rec.group_key]]
    dropped = len(records) - len(kept)
    print(
        f"  Kept {len(kept):,} / {len(records):,} records "
        f"({100.0 * dropped / len(records):.1f}% dropped)",
        file=sys.stderr,
    )

    if args.stats:
        print("", file=sys.stderr)
        print_stats(records, kept, cutoffs)

    # Write BED
    if args.output is None:
        out = sys.stdout
        close_out = False
    else:
        out = open(args.output, "w")
        close_out = True

    try:
        for rec in kept:
            out.write(rec.to_bed_line())
    finally:
        if close_out:
            out.close()

    if args.output is not None:
        print(f"Done. Wrote {len(kept):,} BED records to {args.output}.", file=sys.stderr)


if __name__ == "__main__":
    main()
