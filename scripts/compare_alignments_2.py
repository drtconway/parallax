"""Compare Alignments.

Usage:
    compare_alignments_2.py <lhs> <rhs>


"""

from collections import defaultdict
from collections.abc import Generator
import sys

import docopt
from pysam import AlignmentFile, AlignedSegment
from shapely.geometry import box, MultiPolygon

class ChromCoordinates(object):
    def __init__(self, bam_name):
        with AlignmentFile(bam_name, "rb") as bam:
            names = bam.references
            lengths = bam.lengths
            self.chroms = dict(zip(names, lengths))
            cumu = 0
            self.starts = {}
            for (name, length) in zip(names, lengths):
                self.starts[name] = cumu
                cumu += length

    def pos_to_glob(self, chrom: str, pos: int, is_reverse: bool) -> int:
        glob = self.starts[chrom] + pos
        if is_reverse:
            glob = -glob
        return glob

    def interval_to_glob(self, chrom: str, start: int, end: int, is_reverse: bool) -> tuple[int, int]:
        start = self.starts[chrom] + start
        end = self.starts[chrom] + end
        if is_reverse:
            return (-end, -start)
        else:
            return (start, end)

def query_length(rec: AlignedSegment) -> int:
    l = 0
    for (op, length) in rec.cigartuples:
        if op != 2:
            l += length
    return l

def query_range(rec: AlignedSegment) -> tuple[int, int]:
    cigar = rec.cigartuples
    q1 = 0
    if cigar[0][0] == 4 or cigar[0][0] == 5:
        q1 = cigar[0][1]
    q2 = 0
    if cigar[-1][0] == 4 or cigar[-1][0] == 5:
        q2 = cigar[-1][1]
    l = query_length(rec)
    return (q1, l - q2)

def alignment_groups(bam_name: str) -> Generator[tuple[str, list[AlignedSegment]]]:
    with AlignmentFile(bam_name, "rb") as bam:
        group = []
        group_name = ""
        for record in bam:
            if record.is_unmapped:
                continue
            if record.is_secondary:
                continue
            if record.query_name != group_name:
                if len(group) > 0:
                    yield (group_name, group)
                group = [record]
                group_name = record.query_name
            else:
                group.append(record)
        if len(group) > 0:
            yield (group_name, group)

def alignment_group_boxes(coord_map: ChromCoordinates, group: list[AlignedSegment]) -> MultiPolygon:
    boxes = []
    for rec in group:
        if rec.is_unmapped:
            boxes.append(box(0, 0, 0, 0))
            continue
        chrom = rec.reference_name
        ref_start = rec.reference_start
        ref_end = rec.reference_end
        ref_is_reverse = rec.is_reverse
        (r1, r2) = coord_map.interval_to_glob(chrom, ref_start, ref_end, ref_is_reverse)
        (q1, q2) = query_range(rec)
        boxes.append(box(q1, r1, q2, r2))
    return MultiPolygon(boxes).buffer(0)

def group_jaccard(coord_map: ChromCoordinates, lhs: list[AlignedSegment], rhs: list[AlignedSegment]) -> float:
    lhs_boxes = alignment_group_boxes(coord_map, lhs)
    rhs_boxes = alignment_group_boxes(coord_map, rhs)

    den = lhs_boxes.union(rhs_boxes).area
    if den == 0:
        return 0.0
    
    num = lhs_boxes.intersection(rhs_boxes).area
    return num/den

CLIP_THRESHOLD = 10  # bp


def cigar_ops_from_end(cigar: list[tuple[int, int]], from_right: bool) -> list[tuple[int, int]]:
    """Return cigar ops ordered inward from the specified end, skipping leading clips."""
    ops = list(reversed(cigar)) if from_right else list(cigar)
    # skip terminal hard/soft clips
    while ops and ops[0][0] in (4, 5):
        ops.pop(0)
    return ops


def query_consumed(ops: list[tuple[int, int]]) -> int:
    """Count query bases consumed by a list of cigar ops."""
    return sum(l for op, l in ops if op in (0, 1, 4, 5, 7, 8))  # M, I, S, H, =, X


def is_clip_difference(
    mm2: AlignedSegment,
    plx: AlignedSegment,
) -> bool:
    """Return True if the difference between two single-segment alignments is
    purely a terminal clipping difference satisfying the clip criterion."""

    if mm2.is_reverse != plx.is_reverse:
        return False
    if mm2.reference_name != plx.reference_name:
        return False

    mm2_q = query_range(mm2)
    plx_q = query_range(plx)

    # Check each end that differs
    for side in ('left', 'right'):
        if side == 'left':
            mm2_qpos = mm2_q[0]
            plx_qpos = plx_q[0]
            mm2_rpos = mm2.reference_start
            plx_rpos = plx.reference_start
        else:
            mm2_qpos = mm2_q[1]
            plx_qpos = plx_q[1]
            mm2_rpos = mm2.reference_end
            plx_rpos = plx.reference_end

        if mm2_qpos == plx_qpos:
            continue  # this end agrees — skip

        # Identify long (more query bases) and short
        if side == 'left':
            # long has smaller qStart
            if mm2_qpos < plx_qpos:
                long_rec, long_qpos, long_rpos = mm2, mm2_qpos, mm2_rpos
                short_qpos, short_rpos = plx_qpos, plx_rpos
            else:
                long_rec, long_qpos, long_rpos = plx, plx_qpos, plx_rpos
                short_qpos, short_rpos = mm2_qpos, mm2_rpos
            q_diff = short_qpos - long_qpos
            r_diff = abs(short_rpos - long_rpos)
            from_right = False
        else:
            # long has larger qEnd
            if mm2_qpos > plx_qpos:
                long_rec, long_qpos, long_rpos = mm2, mm2_qpos, mm2_rpos
                short_qpos, short_rpos = plx_qpos, plx_rpos
            else:
                long_rec, long_qpos, long_rpos = plx, plx_qpos, plx_rpos
                short_qpos, short_rpos = mm2_qpos, mm2_rpos
            q_diff = long_qpos - short_qpos
            r_diff = abs(long_rpos - short_rpos)
            from_right = True

        # Query and reference extent of the difference must both be within threshold
        if q_diff > CLIP_THRESHOLD or r_diff > CLIP_THRESHOLD:
            return False

        # Diagonal at short's endpoint must match long's diagonal at that point.
        # Diagonal = refPos - qPos (+ strand) or refPos + qPos (- strand).
        if not long_rec.is_reverse:
            long_diag  = long_rpos  - long_qpos
            short_diag = short_rpos - short_qpos
        else:
            long_diag  = long_rpos  + long_qpos
            short_diag = short_rpos + short_qpos
        if long_diag != short_diag:
            return False

        # Scan long's CIGAR inward from this end; the first non-matching event
        # (X, I, or D — treating M as matching) must account for short's endpoint.
        ops = cigar_ops_from_end(long_rec.cigartuples, from_right)
        consumed = 0
        found_mismatch = False
        for op, length in ops:
            if op in (0, 7):  # M or = : matching, consume and continue
                consumed += length
                continue
            if op in (8, 1, 2):  # X, I, D : non-matching event
                # short's endpoint should fall at or within this op
                if side == 'left':
                    pos_at_op = long_qpos + consumed
                    if pos_at_op <= short_qpos <= long_qpos + consumed + length:
                        found_mismatch = True
                else:
                    pos_at_op = long_qpos - consumed
                    if long_qpos - consumed - length <= short_qpos <= pos_at_op:
                        found_mismatch = True
                break
            # soft/hard clips already skipped above; anything else → fail
            break

        if not found_mismatch:
            return False

    return True


def get_next(itr):
    try:
        return itr.__next__()
    except StopIteration:
        return None

def main(args):
    lhs_bam = args['<lhs>']
    rhs_bam = args['<rhs>']

    coord_map = ChromCoordinates(lhs_bam)

    lhs_groups = alignment_groups(lhs_bam)
    lhs_curr = get_next(lhs_groups)

    rhs_groups = alignment_groups(rhs_bam)
    rhs_curr = get_next(rhs_groups)

    i = 0
    e = 0
    print('read_id\tstate\tjaccard\tlhs_segs\trhs_segs')
    while lhs_curr is not None and rhs_curr is not None:
        i += 1
        (lhs_name, lhs_group) = lhs_curr
        (rhs_name, rhs_group) = rhs_curr

        if lhs_name < rhs_name:
            print(f"warning: rhs has no alignments for {lhs_name}", file=sys.stderr)
            print(f'{lhs_name}\tlhs\t0.0\t{len(lhs_group)}\t0')
            lhs_curr = get_next(lhs_groups)
            continue

        if rhs_name < lhs_name:
            print(f"warning: lhs has no alignments for {rhs_name}", file=sys.stderr)
            print(f'{rhs_name}\trhs\t0.0\t0\t{len(rhs_group)}')
            rhs_curr = get_next(rhs_groups)
            continue

        j = group_jaccard(coord_map, lhs_group, rhs_group)

        if j < 1.0:
            if (len(lhs_group) == 1 and len(rhs_group) == 1
                    and is_clip_difference(lhs_group[0], rhs_group[0])):
                print(f'{lhs_name}\tclip\t{j:2.4f}\t{len(lhs_group)}\t{len(rhs_group)}')
            else:
                print(f'{lhs_name}\tboth\t{j:2.4f}\t{len(lhs_group)}\t{len(rhs_group)}')
            #lhs_items = set()
            #for lhs in lhs_group:
            #    chrom = lhs.reference_name
            #    ref_start = lhs.reference_start
            #    ref_end = lhs.reference_end
            #    (qry_start, qry_end) = query_range(lhs)
            #    is_reverse = lhs.is_reverse
            #    lhs_items.add((qry_start, qry_end, chrom, ref_start, ref_end, is_reverse))
            #rhs_items = set()
            #for rhs in rhs_group:
            #    chrom = rhs.reference_name
            #    ref_start = rhs.reference_start
            #    ref_end = rhs.reference_end
            #    (qry_start, qry_end) = query_range(rhs)
            #    is_reverse = rhs.is_reverse
            #    rhs_items.add((qry_start, qry_end, chrom, ref_start, ref_end, is_reverse))
            #
            #lhs_only = lhs_items - rhs_items
            #rhs_only = rhs_items - lhs_items
            #
            #for (qry_start, qry_end, chrom, ref_start, ref_end, is_reverse) in sorted(lhs_only):
            #    print(f"{lhs_name}\tlhs\t{qry_start}\t{qry_end}\t{chrom}\t{ref_start}\t{ref_end}\t{is_reverse}\t{j:2.4f}")
            #for (qry_start, qry_end, chrom, ref_start, ref_end, is_reverse) in sorted(rhs_only):
            #    print(f"{lhs_name}\trhs\t{qry_start}\t{qry_end}\t{chrom}\t{ref_start}\t{ref_end}\t{is_reverse}\t{j:2.4f}")
            #sys.stdout.flush()
        else:
            e += 1

        lhs_curr = get_next(lhs_groups)
        rhs_curr = get_next(rhs_groups)
    print(f'total records: {i}', file=sys.stderr)
    print(f'total matching: {e}', file=sys.stderr)
    print(f'total mismatching: {i - e}', file=sys.stderr)

if __name__ == '__main__':
    args = docopt.docopt(__doc__, version='Compare Alignments 2.0')
    main(args)