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
    while lhs_curr is not None and rhs_curr is not None:
        i += 1
        (lhs_name, lhs_group) = lhs_curr
        (rhs_name, rhs_group) = rhs_curr

        if lhs_name < rhs_name:
            print(f"warning: rhs has no alignments for {lhs_name}", file=sys.stderr)
            lhs_curr = get_next(lhs_groups)
            continue

        if rhs_name < lhs_name:
            print(f"warning: lhs has no alignments for {rhs_name}", file=sys.stderr)
            rhs_curr = get_next(rhs_groups)
            continue

        j = group_jaccard(coord_map, lhs_group, rhs_group)

        if j < 1.0:
            print(f'{lhs_name}')
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