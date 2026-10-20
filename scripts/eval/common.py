"""Shared helpers for the alignment evaluation scripts used by evaluate.nf.

Coordinates on the read are always *forward-strand* query coordinates
(0-based, half-open), whatever the strand of the alignment, and include
hard-clipped bases, so the records of one read can be compared directly.
"""

import bisect
import gzip
from itertools import groupby

import pysam

M, I, D, N, S, H, P, EQ, X = range(9)
QUERY_OPS = (M, I, S, EQ, X, H)
MATCH_OPS = (M, EQ, X)


def open_text(path, mode="rt"):
    """Open a plain or gzipped text file."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode.replace("t", ""))


def record_class(rec):
    if rec.is_secondary:
        return "secondary"
    if rec.is_supplementary:
        return "supplementary"
    return "primary"


def read_length(rec):
    """Full read length, counting hard clips."""
    if rec.cigartuples:
        return sum(n for op, n in rec.cigartuples if op in QUERY_OPS)
    return rec.query_length or 0


def clips(rec):
    """(leading, trailing) clip lengths in alignment orientation, S and H."""
    cig = rec.cigartuples
    lead = trail = 0
    for op, n in cig:
        if op not in (S, H):
            break
        lead += n
    for op, n in reversed(cig):
        if op not in (S, H):
            break
        trail += n
    return lead, trail


def query_range(rec, length=None):
    """Forward-strand query interval [qs, qe) covered by the alignment."""
    length = read_length(rec) if length is None else length
    lead, trail = clips(rec)
    qs, qe = lead, length - trail
    if rec.is_reverse:
        qs, qe = length - qe, length - qs
    return qs, qe


def size_band(n):
    """Lower bound of the power-of-two size band containing n >= 1."""
    return 1 << (n.bit_length() - 1)


def dist_bin(d):
    """Lower bound of the log2 distance bin containing d >= 0 (0, 1, 2, 4, ...)."""
    return 0 if d == 0 else 1 << (d.bit_length() - 1)


def count_dist_bins(m, into, key_fn):
    """Add one count to into[key_fn(bin)] for each distance 0..m inclusive."""
    if m < 0:
        return
    into[key_fn(0)] += 1
    lo = 1
    while lo <= m:
        hi = min(2 * lo - 1, m)
        into[key_fn(lo)] += hi - lo + 1
        lo *= 2


def group_by_name(bam):
    """Yield (name, [records]) from a name-grouped BAM."""
    for name, recs in groupby(bam.fetch(until_eof=True), key=lambda r: r.query_name):
        yield name, list(recs)


def union_length(intervals):
    total = 0
    end = None
    for s, e in sorted(intervals):
        if end is None or s > end:
            total += e - s
            end = e
        elif e > end:
            total += e - end
            end = e
    return total


def merge_intervals(intervals):
    merged = []
    for s, e in sorted(intervals):
        if merged and s <= merged[-1][1]:
            if e > merged[-1][1]:
                merged[-1][1] = e
        else:
            merged.append([s, e])
    return merged


class Strata:
    """Stratum lookup from a TSV of name<TAB>bed lines (no header needed).

    An interval is assigned to the first stratum listed that covers at least
    `min_frac` of it, otherwise to "other".
    """

    def __init__(self, tsv=None, min_frac=0.5):
        self.names = []
        self.tables = []
        self.min_frac = min_frac
        if tsv:
            with open(tsv) as f:
                for line in f:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 2 or line.startswith("#") or fields[:2] == ["name", "bed"]:
                        continue
                    self.names.append(fields[0])
                    self.tables.append(self._load_bed(fields[1]))

    @staticmethod
    def _load_bed(path):
        by_chrom = {}
        with open_text(path) as f:
            for line in f:
                if line.startswith(("#", "track", "browser")):
                    continue
                fields = line.split("\t", 3)
                if len(fields) < 3:
                    continue
                by_chrom.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
        table = {}
        for chrom, ivs in by_chrom.items():
            merged = merge_intervals(ivs)
            table[chrom] = ([s for s, _ in merged], [e for _, e in merged])
        return table

    def _overlap(self, table, chrom, start, end):
        if chrom not in table:
            return 0
        starts, ends = table[chrom]
        j = bisect.bisect_right(ends, start)
        total = 0
        while j < len(starts) and starts[j] < end:
            total += min(end, ends[j]) - max(start, starts[j])
            j += 1
        return total

    def classify(self, chrom, start, end):
        if end <= start:
            end = start + 1
        need = self.min_frac * (end - start)
        for name, table in zip(self.names, self.tables):
            if self._overlap(table, chrom, start, end) >= need:
                return name
        return "other"


# --- sampled query -> reference maps, for comparing alignments of one read ---

UNMAPPED = -1


def encode(tid, reverse, rpos):
    return (tid << 33) | (int(reverse) << 32) | rpos


def flip_strand(code):
    return code ^ (1 << 32) if code >= 0 else code


def decode(code):
    return code >> 33, bool((code >> 32) & 1), code & 0xFFFFFFFF


def map_positions(records, length, positions, flip=False, offset=0):
    """Map forward-strand query positions to encoded reference positions.

    `positions` are positions in the *target* coordinate system.  The records
    belong to a read whose forward coordinate q relates to a target position t
    by t = q + offset, or, if `flip` (the read is the reverse complement of the
    target), by t = length - 1 - q + offset (strand is flipped to match).
    Uses primary and supplementary records; returns a list aligned with
    `positions`, holding UNMAPPED where no record places the base.
    """
    out = [UNMAPPED] * len(positions)
    index = {t: i for i, t in enumerate(positions)}
    for rec in records:
        if rec.is_unmapped or rec.is_secondary or not rec.cigartuples:
            continue
        # Target positions -> this record's alignment-orientation coordinates.
        to_aln = []
        for t in positions:
            q = (length - 1 - (t - offset)) if flip else (t - offset)
            if 0 <= q < length:
                qa = (length - 1 - q) if rec.is_reverse else q
                to_aln.append((qa, t))
        to_aln.sort()
        keys = [qa for qa, _ in to_aln]
        reverse = rec.is_reverse != flip
        tid = rec.reference_id
        q = 0
        r = rec.reference_start
        for op, n in rec.cigartuples:
            if op in MATCH_OPS:
                j = bisect.bisect_left(keys, q)
                while j < len(keys) and keys[j] < q + n:
                    qa, t = to_aln[j]
                    i = index[t]
                    if out[i] == UNMAPPED:
                        out[i] = encode(tid, reverse, r + qa - q)
                    j += 1
                q += n
                r += n
            elif op in (I, S, H):
                q += n
            elif op in (D, N):
                r += n
    return out


def compare_maps(a, b):
    """(agree, either) counts over two equal-length encoded maps."""
    agree = either = 0
    for x, y in zip(a, b):
        if x == UNMAPPED and y == UNMAPPED:
            continue
        either += 1
        if x == y:
            agree += 1
    return agree, either


def open_bam(path):
    return pysam.AlignmentFile(path, "rb", check_sq=False)
