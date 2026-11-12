"""Compare alignments of a read set with alignments of a transformed copy.

Usage:
    invariance.py <orig.bam> <variant.bam> --kind rc|frag --out OUT.tsv.gz [--strata TSV]

A good aligner should place the same bases in the same place however the
read is presented.  Both BAMs must come from the same aligner and reference.

--kind rc    the variant reads are reverse complements of the originals
             (same names).  Each read is classified as identical, near,
             partial, disjoint, orig_only, variant_only or both_unmapped.
--kind frag  the variant reads are tiles named <read>/<start>-<end>.  Each
             tile is compared with the part of its parent's alignment that
             covers it, and classified as agree, partial, disagree_confident
             (MAPQ >= --min-mapq), disagree_lowmapq, unaligned (the parent
             aligned it but the tile did not), parent_unaligned or
             both_unaligned.

Comparison samples every --step query bases and asks whether both alignments
put that base at the same reference position and strand.  The agreement
score is agree / positions placed by either.

Primary and supplementary records only.  Variant records of each read must be
adjacent (name-sorted or straight from the aligner); the original BAM is held
in memory in compact form.
"""

import argparse
from array import array

from common import (UNMAPPED, Strata, compare_maps, decode, group_by_name, map_positions,
                    open_bam, open_text, read_length)


def main_records(recs):
    return [r for r in recs if not r.is_unmapped and not r.is_secondary and r.cigartuples]


def load_orig(path, step, strata):
    """read -> (length, positions map as array, n records, stratum, max mapq)."""
    reads = {}
    with open_bam(path) as bam:
        names = bam.references
        for name, recs in group_by_name(bam):
            main = main_records(recs)
            if main:
                length = read_length(main[0])
            else:
                length = max((r.query_length or 0) for r in recs)
            positions = list(range(0, length, step))
            codes = map_positions(main, length, positions) if main else [UNMAPPED] * len(positions)
            primary = next((r for r in main if not r.is_supplementary), None)
            stratum = (strata.classify(primary.reference_name, primary.reference_start,
                                       primary.reference_end) if primary else "NA")
            reads[name] = (length, array("q", codes), len(main), stratum)
    return reads, names


def classify_rc(agree, either, n_orig, n_var):
    if n_orig == 0 and n_var == 0:
        return "both_unmapped"
    if n_var == 0:
        return "orig_only"
    if n_orig == 0:
        return "variant_only"
    if either == 0:
        return "disjoint"
    score = agree / either
    if score == 1.0 and n_orig == n_var:
        return "identical"
    if score >= 0.95:
        return "near"
    if score > 0:
        return "partial"
    return "disjoint"


def run_rc(args, orig, out):
    out.write("read\tread_len\tclass\tscore\tn_orig\tn_variant\tstratum\n")
    seen = set()
    with open_bam(args.variant) as bam:
        for name, recs in group_by_name(bam):
            seen.add(name)
            if name not in orig:
                continue
            length, codes, n_orig, stratum = orig[name]
            main = main_records(recs)
            positions = list(range(0, length, args.step))
            var = map_positions(main, length, positions, flip=True) if main else [UNMAPPED] * len(positions)
            agree, either = compare_maps(codes, var)
            cls = classify_rc(agree, either, n_orig, len(main))
            score = f"{agree / either:.4f}" if either else "NA"
            out.write(f"{name}\t{length}\t{cls}\t{score}\t{n_orig}\t{len(main)}\t{stratum}\n")
    for name in orig.keys() - seen:
        length, _, n_orig, stratum = orig[name]
        cls = "orig_only" if n_orig else "both_unmapped"
        out.write(f"{name}\t{length}\t{cls}\tNA\t{n_orig}\t0\t{stratum}\n")


def parse_tile(name):
    parent, _, span = name.rpartition("/")
    start, _, end = span.partition("-")
    return parent, int(start), int(end)


def run_frag(args, orig, ref_names, strata, out):
    out.write("tile\tparent\tstart\tend\tclass\tscore\tmapq\tstratum\n")
    step = args.step
    with open_bam(args.variant) as bam:
        for name, recs in group_by_name(bam):
            parent, start, end = parse_tile(name)
            if parent not in orig:
                continue
            _, codes, _, _ = orig[parent]
            # Parent samples inside the tile, as tile-relative positions.
            first = (start + step - 1) // step
            last = (end - 1) // step
            expected = codes[first:last + 1]
            positions = [i * step for i in range(first, last + 1)]
            main = main_records(recs)
            mapq = max((r.mapping_quality for r in main), default=0)
            got = (map_positions(main, end - start, positions, offset=start)
                   if main else [UNMAPPED] * len(positions))

            placed = [c for c in expected if c != UNMAPPED]
            if placed:
                tid, _, r0 = decode(placed[0])
                r1 = decode(placed[-1])[2]
                stratum = strata.classify(ref_names[tid], min(r0, r1), max(r0, r1) + 1)
            elif main:
                stratum = strata.classify(main[0].reference_name, main[0].reference_start,
                                          main[0].reference_end)
            else:
                stratum = "NA"

            agree, either = compare_maps(expected, got)
            has_expected = bool(placed)
            has_got = any(c != UNMAPPED for c in got)
            score = agree / either if either else 0.0
            if not has_expected and not has_got:
                cls = "both_unaligned"
            elif not has_expected:
                cls = "parent_unaligned"
            elif not has_got:
                cls = "unaligned"
            elif score >= 0.9:
                cls = "agree"
            elif score > 0:
                cls = "partial"
            elif mapq >= args.min_mapq:
                cls = "disagree_confident"
            else:
                cls = "disagree_lowmapq"
            score_s = f"{score:.4f}" if either else "NA"
            out.write(f"{name}\t{parent}\t{start}\t{end}\t{cls}\t{score_s}\t{mapq}\t{stratum}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("orig")
    ap.add_argument("variant")
    ap.add_argument("--kind", choices=["rc", "frag"], required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--strata")
    ap.add_argument("--step", type=int, default=25, help="sample every STEP query bases")
    ap.add_argument("--min-mapq", type=int, default=20)
    args = ap.parse_args()

    strata = Strata(args.strata)
    orig, ref_names = load_orig(args.orig, args.step, strata)
    with open_text(args.out, "wt") as out:
        if args.kind == "rc":
            run_rc(args, orig, out)
        else:
            run_frag(args, orig, ref_names, strata, out)


if __name__ == "__main__":
    main()
