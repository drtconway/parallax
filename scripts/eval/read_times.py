"""Per-read alignment times from the XT:f tags of a parallax BAM.

Usage:
    read_times.py <aln.bam> --prefix P [--top-fraction 0.01]

Writes
  P.read_times.tsv.gz   read, read_len, mapped, time_s, time_per_base_s,
                        slowest first by time per base
  P.slow_reads.fq.gz    the top fraction of reads by time per base, as FASTQ
                        in their original orientation

parallax writes XT:f (seconds to align the read, the same on every record of
the read) when run with `[seeding] debug_read_time_tag = true`.  Each read is
taken from its primary or unmapped record.  If the BAM has no XT tags, both
outputs are written empty (TSV header only) and a warning is printed.
"""

import argparse
import gzip
import math
import sys

import numpy as np

from common import open_bam

SKIP = 0x100 | 0x800  # secondary, supplementary


def representative_records(path):
    """The primary or unmapped record of each read."""
    with open_bam(path) as bam:
        for rec in bam.fetch(until_eof=True):
            if not rec.flag & SKIP:
                yield rec


def collect_times(path):
    names, lengths, mapped, times = [], [], [], []
    no_tag = 0
    for rec in representative_records(path):
        if not rec.has_tag("XT"):
            no_tag += 1
            continue
        names.append(rec.query_name)
        lengths.append(rec.infer_read_length() or rec.query_length)
        mapped.append(0 if rec.is_unmapped else 1)
        times.append(rec.get_tag("XT"))
    return (names, np.array(lengths, dtype=np.int64), np.array(mapped, dtype=np.int8),
            np.array(times, dtype=np.float64), no_tag)


def write_fastq(path, wanted, out_path):
    """Write the reads named in `wanted` from their primary/unmapped records."""
    written, clipped = 0, 0
    with gzip.open(out_path, "wt", compresslevel=6) as out:
        for rec in representative_records(path):
            if rec.query_name not in wanted:
                continue
            seq = rec.get_forward_sequence()
            if not seq or len(seq) != rec.infer_read_length():
                clipped += 1  # hard-clipped or SEQ omitted: cannot recover the read
                continue
            quals = rec.get_forward_qualities()
            qual = "".join(chr(q + 33) for q in quals) if quals is not None else "I" * len(seq)
            out.write(f"@{rec.query_name}\n{seq}\n+\n{qual}\n")
            written += 1
    return written, clipped


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam")
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--top-fraction", type=float, default=0.01,
                    help="fraction of reads, slowest by time per base, to write as FASTQ")
    args = ap.parse_args()

    names, lengths, mapped, times, no_tag = collect_times(args.bam)
    per_base = np.divide(times, lengths, out=np.full_like(times, np.nan), where=lengths > 0)
    # slowest first; reads of length 0 (NaN) last
    order = np.argsort(np.where(np.isnan(per_base), -np.inf, -per_base), kind="stable")

    with gzip.open(f"{args.prefix}.read_times.tsv.gz", "wt", compresslevel=6) as out:
        out.write("read\tread_len\tmapped\ttime_s\ttime_per_base_s\n")
        for i in order:
            out.write(f"{names[i]}\t{lengths[i]}\t{mapped[i]}\t{times[i]:.6g}\t{per_base[i]:.6g}\n")

    n = len(names)
    n_top = math.ceil(args.top_fraction * n) if n else 0
    wanted = {names[i] for i in order[:n_top]}
    written, clipped = write_fastq(args.bam, wanted, f"{args.prefix}.slow_reads.fq.gz")

    if n == 0:
        print(f"warning: no XT tags in {args.bam}; was parallax run with "
              "[seeding] debug_read_time_tag = true?", file=sys.stderr)
    elif no_tag:
        print(f"warning: {no_tag} reads without an XT tag were ignored", file=sys.stderr)
    if clipped:
        print(f"warning: {clipped} slow reads had no full sequence in the BAM and were not written",
              file=sys.stderr)
    if n:
        cutoff = per_base[order[n_top - 1]] if n_top else float("nan")
        print(f"{n} reads timed, {times.sum():.1f} s total; top {n_top} by time per base "
              f"(>= {cutoff:.3g} s/base) written: {written}", file=sys.stderr)


if __name__ == "__main__":
    main()
