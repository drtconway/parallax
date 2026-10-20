"""Per-alignment and per-read statistics for truth-free aligner evaluation.

Usage:
    alignment_stats.py <name-grouped.bam> --prefix P [--reference FA] [--strata TSV]

Writes:
  P.alignments.tsv.gz  one row per mapped record: class, coordinates, MAPQ,
                       match/mismatch/indel counts, identities, stratum, and
                       the kind of each end (read_end, junction or clipped)
  P.reads.tsv.gz       one row per read: length, alignment counts, bases
                       aligned (union), overlap, interior gaps, end clips,
                       identity windows and outliers, split-implied indels
  P.indels.tsv         indel events and bases by size band (powers of two),
                       from CIGARs and implied by colinear split alignments
  P.errpos.tsv         error events and aligned positions by distance from
                       the nearest alignment end, per end kind
  P.outliers.tsv.gz    read windows whose identity is far below the read's median

Records of each read must be adjacent (name-sorted or straight from the
aligner).  CIGARs with =/X are used directly; M is resolved against the
reference if one is given, else from the NM tag.
"""

import argparse
import statistics
from collections import Counter, defaultdict

import pysam

from common import (D, EQ, H, I, M, N, S, X, Strata, clips, count_dist_bins, dist_bin,
                    group_by_name, open_bam, open_text, query_range, read_length,
                    record_class, size_band, union_length)

ALN_COLUMNS = [
    "read", "class", "read_len", "qs", "qe", "chrom", "rs", "re", "strand", "mapq",
    "matches", "mismatches", "ins_events", "ins_bases", "del_events", "del_bases",
    "identity", "gc_identity", "stratum", "end5", "end3",
]
READ_COLUMNS = [
    "read", "read_len", "mapped", "n_primary", "n_supplementary", "n_secondary",
    "aligned_union", "aligned_sum", "interior_gap", "clip5", "clip3",
    "n_windows", "n_outlier_windows", "split_events",
]


class Windows:
    """Per-window identity counts over one read's forward coordinates."""

    def __init__(self, length, size):
        self.size = size
        n = max(1, (length + size - 1) // size)
        # aligned query bases, matches, errors (mismatch bases + ins bases + del bases)
        self.aligned = [0] * n
        self.matches = [0] * n
        self.errors = [0] * n

    def add_run(self, q0, q1, matches, errors):
        """Spread a run of aligned query bases [q0, q1) over the windows."""
        w = q0 // self.size
        while q0 < q1:
            end = min(q1, (w + 1) * self.size)
            n = end - q0
            self.aligned[w] += n
            if matches:
                self.matches[w] += n
            if errors:
                self.errors[w] += n
            q0 = end
            w += 1

    def add_point_errors(self, q, n):
        w = min(q // self.size, len(self.errors) - 1)
        self.errors[w] += n

    def identities(self):
        out = []
        for i, (a, m, e) in enumerate(zip(self.aligned, self.matches, self.errors)):
            if 2 * a >= self.size and m + e > 0:
                out.append((i, m / (m + e)))
        return out


class Stats:
    def __init__(self, args):
        self.args = args
        self.strata = Strata(args.strata)
        self.fasta = pysam.FastaFile(args.reference) if args.reference else None
        self.indels = Counter()   # (class, stratum, kind, source, band) -> events
        self.indel_bases = Counter()
        self.err = Counter()      # (end_kind, bin, what) -> count
        self.aln_out = open_text(f"{args.prefix}.alignments.tsv.gz", "wt")
        self.read_out = open_text(f"{args.prefix}.reads.tsv.gz", "wt")
        self.outlier_out = open_text(f"{args.prefix}.outliers.tsv.gz", "wt")
        self.aln_out.write("\t".join(ALN_COLUMNS) + "\n")
        self.read_out.write("\t".join(READ_COLUMNS) + "\n")
        self.outlier_out.write("read\tread_len\twin_start\twin_end\tidentity\tread_median\n")

    def end_kind(self, clip, boundary, others):
        """Kind of an alignment end at forward coordinate `boundary`."""
        if clip < self.args.near_end:
            return "read_end"
        for qs, qe in others:
            if abs(qs - boundary) < self.args.near_end or abs(qe - boundary) < self.args.near_end:
                return "junction"
        return "clipped"

    def mismatches_in_m(self, rec, q_seq, r, n, ref_seq, ref_start):
        query = rec.query_sequence
        if query is None or ref_seq is None:
            return None
        a = query[q_seq:q_seq + n].upper()
        b = ref_seq[r - ref_start:r - ref_start + n]
        if a == b:
            return []
        return [i for i in range(n) if a[i] != b[i]]

    def alignment(self, rec, length, ranges, windows):
        """Walk one record; returns its output row (as a list)."""
        args = self.args
        cls = record_class(rec)
        qs, qe = query_range(rec, length)
        lead, trail = clips(rec)
        hard_lead = rec.cigartuples[0][1] if rec.cigartuples[0][0] == H else 0
        others = [rg for rg in ranges if rg != (qs, qe)]
        chrom = rec.reference_name
        stratum = self.strata.classify(chrom, rec.reference_start, rec.reference_end)

        ref_seq = None
        if self.fasta and any(op == M for op, _ in rec.cigartuples) and chrom in self.fasta:
            ref_seq = self.fasta.fetch(chrom, rec.reference_start, rec.reference_end).upper()

        matches = mismatches = ins_e = ins_b = del_e = del_b = 0
        unresolved_m = 0
        events = []  # (alignment-orientation q, kind)
        q = 0
        r = rec.reference_start
        rev = rec.is_reverse

        def fwd(q0, q1):
            return (length - q1, length - q0) if rev else (q0, q1)

        for op, n in rec.cigartuples:
            if op == EQ or op == X:
                f0, f1 = fwd(q, q + n)
                if windows:
                    windows.add_run(f0, f1, op == EQ, op == X)
                if op == EQ:
                    matches += n
                else:
                    mismatches += n
                    events.extend((q + i, "mismatch") for i in range(n))
                q += n
                r += n
            elif op == M:
                mm = self.mismatches_in_m(rec, q - hard_lead, r, n, ref_seq, rec.reference_start)
                f0, f1 = fwd(q, q + n)
                if mm is None:
                    unresolved_m += n
                    if windows:
                        windows.add_run(f0, f1, False, False)
                else:
                    matches += n - len(mm)
                    mismatches += len(mm)
                    events.extend((q + i, "mismatch") for i in mm)
                    if windows:
                        windows.add_run(f0, f1, True, False)
                        for i in mm:
                            fq = (length - 1 - (q + i)) if rev else q + i
                            windows.matches[fq // windows.size] -= 1
                            windows.errors[fq // windows.size] += 1
                q += n
                r += n
            elif op == I:
                ins_e += 1
                ins_b += n
                events.append((q, "ins"))
                self.indels[(cls, stratum, "ins", "cigar", size_band(n))] += 1
                self.indel_bases[(cls, stratum, "ins", "cigar", size_band(n))] += n
                if windows:
                    f0, f1 = fwd(q, q + n)
                    windows.add_run(f0, f1, False, True)
                q += n
            elif op == D:
                del_e += 1
                del_b += n
                events.append((q, "del"))
                self.indels[(cls, stratum, "del", "cigar", size_band(n))] += 1
                self.indel_bases[(cls, stratum, "del", "cigar", size_band(n))] += n
                if windows:
                    fq = (length - q) if rev else q
                    windows.add_point_errors(min(fq, length - 1), n)
                r += n
            elif op == N:
                r += n
            elif op in (S, H):
                q += n

        if unresolved_m:
            nm = rec.get_tag("NM") if rec.has_tag("NM") else None
            if nm is not None:
                mm = max(0, nm - ins_b - del_b)
                mismatches += mm
                matches += unresolved_m - mm
            else:
                matches += unresolved_m

        # Ends in alignment orientation: lead side then trail side.
        a0, a1 = lead, length - trail
        lead_boundary = (length - a0) if rev else a0
        trail_boundary = (length - a1) if rev else a1
        kind_lead = self.end_kind(lead, lead_boundary, others)
        kind_trail = self.end_kind(trail, trail_boundary, others)
        if cls != "secondary":
            span = a1 - a0
            m = (span - 1) // 2
            count_dist_bins(m, self.err, lambda b: (kind_lead, b, "positions"))
            count_dist_bins(span - 2 - m, self.err, lambda b: (kind_trail, b, "positions"))
            for eq, what in events:
                dl, dr = eq - a0, a1 - 1 - eq
                if dl <= dr:
                    self.err[(kind_lead, dist_bin(max(dl, 0)), what)] += 1
                else:
                    self.err[(kind_trail, dist_bin(max(dr, 0)), what)] += 1

        denom = matches + mismatches + ins_b + del_b
        gc_denom = matches + mismatches + ins_e + del_e
        end5, end3 = (kind_trail, kind_lead) if rev else (kind_lead, kind_trail)
        return [
            rec.query_name, cls, length, qs, qe, chrom, rec.reference_start, rec.reference_end,
            "-" if rev else "+", rec.mapping_quality, matches, mismatches, ins_e, ins_b, del_e,
            del_b, f"{matches / denom:.5f}" if denom else "NA",
            f"{matches / gc_denom:.5f}" if gc_denom else "NA", stratum, end5, end3,
        ]

    def split_events(self, recs, length):
        """Indels implied by consecutive colinear primary/supplementary records."""
        segs = []
        for rec in recs:
            qs, qe = query_range(rec, length)
            segs.append((qs, qe, rec))
        segs.sort(key=lambda s: s[0])
        count = 0
        for (aqs, aqe, a), (bqs, bqe, b) in zip(segs, segs[1:]):
            if a.reference_id != b.reference_id or a.is_reverse != b.is_reverse:
                continue
            gq = bqs - aqe
            gr = (a.reference_start - b.reference_end) if a.is_reverse else (b.reference_start - a.reference_end)
            slack = self.args.near_end
            if gq < -slack or gr < -slack or gq > self.args.max_split_gap or gr > self.args.max_split_gap:
                continue
            d = gr - gq
            if d == 0:
                continue
            kind = "del" if d > 0 else "ins"
            size = abs(d)
            lo = b.reference_end if a.is_reverse else a.reference_end
            stratum = self.strata.classify(a.reference_name, lo, lo + max(gr, 1))
            cls = record_class(a) if record_class(a) == record_class(b) else "supplementary"
            self.indels[(cls, stratum, kind, "split", size_band(size))] += 1
            self.indel_bases[(cls, stratum, kind, "split", size_band(size))] += size
            count += 1
        return count

    def read(self, name, recs):
        mapped = [r for r in recs if not r.is_unmapped and r.cigartuples]
        if mapped:
            length = read_length(next((r for r in mapped if not r.is_secondary), mapped[0]))
        else:
            length = max((r.query_length or 0) for r in recs)
        counts = Counter(record_class(r) for r in mapped)
        main = [r for r in mapped if not r.is_secondary]
        ranges = [query_range(r, length) for r in main]
        windows = Windows(length, self.args.window) if main and length else None

        for rec in mapped:
            row = self.alignment(rec, length, ranges, windows if not rec.is_secondary else None)
            self.aln_out.write("\t".join(map(str, row)) + "\n")

        union = union_length(ranges)
        aligned_sum = sum(qe - qs for qs, qe in ranges)
        if ranges:
            first = min(qs for qs, _ in ranges)
            last = max(qe for _, qe in ranges)
            clip5, clip3 = first, length - last
            interior = (last - first) - union
        else:
            clip5, clip3, interior = length, 0, 0

        n_windows = n_out = 0
        if windows:
            ids = windows.identities()
            n_windows = len(ids)
            if len(ids) >= 3:
                values = [v for _, v in ids]
                med = statistics.median(values)
                mad = statistics.median(abs(v - med) for v in values)
                threshold = med - max(self.args.outlier_k * 1.4826 * mad, self.args.min_drop)
                for i, v in ids:
                    if v < threshold:
                        n_out += 1
                        s = i * windows.size
                        self.outlier_out.write(
                            f"{name}\t{length}\t{s}\t{min(s + windows.size, length)}\t{v:.4f}\t{med:.4f}\n")

        splits = self.split_events(main, length) if len(main) > 1 else 0
        self.read_out.write("\t".join(map(str, [
            name, length, int(bool(mapped)), counts["primary"], counts["supplementary"],
            counts["secondary"], union, aligned_sum, interior, clip5, clip3, n_windows, n_out,
            splits,
        ])) + "\n")

    def finish(self):
        for f in (self.aln_out, self.read_out, self.outlier_out):
            f.close()
        with open(f"{self.args.prefix}.indels.tsv", "w") as f:
            f.write("class\tstratum\tkind\tsource\tband\tevents\tbases\n")
            for key in sorted(self.indels):
                f.write("\t".join(map(str, key)) + f"\t{self.indels[key]}\t{self.indel_bases[key]}\n")
        table = defaultdict(Counter)
        for (kind, b, what), n in self.err.items():
            table[(kind, b)][what] += n
        with open(f"{self.args.prefix}.errpos.tsv", "w") as f:
            f.write("end_kind\tdist_bin\tpositions\tmismatches\tins_events\tdel_events\n")
            for (kind, b) in sorted(table):
                c = table[(kind, b)]
                f.write(f"{kind}\t{b}\t{c['positions']}\t{c['mismatch']}\t{c['ins']}\t{c['del']}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam")
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--reference")
    ap.add_argument("--strata")
    ap.add_argument("--window", type=int, default=500, help="identity window size")
    ap.add_argument("--outlier-k", type=float, default=5.0, help="outlier threshold in robust SDs")
    ap.add_argument("--min-drop", type=float, default=0.05,
                    help="minimum identity drop below the read median for an outlier window")
    ap.add_argument("--near-end", type=int, default=50,
                    help="clips shorter than this count as read ends; junction tolerance")
    ap.add_argument("--max-split-gap", type=int, default=1_000_000,
                    help="largest query/reference gap between split alignments read as an indel")
    args = ap.parse_args()

    stats = Stats(args)
    with open_bam(args.bam) as bam:
        for name, recs in group_by_name(bam):
            stats.read(name, recs)
    stats.finish()


if __name__ == "__main__":
    main()
