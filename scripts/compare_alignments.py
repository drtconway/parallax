#!/usr/bin/env python3
"""
Compare primary alignments between two BAM files (e.g., parallax vs minimap2).

Outputs:
- Summary statistics (concordant, discordant, unique to each)
- Optional TSV with per-read details
"""
from __future__ import annotations
import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional

try:
    import pysam
except ImportError:
    print("Error: pysam is required. Install with: pip install pysam", file=sys.stderr)
    sys.exit(1)


@dataclass
class Alignment:
    """Alignment info for a read segment."""
    chrom: str
    pos: int  # 0-based
    end: int  # 0-based, exclusive
    strand: str  # '+' or '-'
    mapq: int
    cigar: str
    is_supplementary: bool
    is_unmapped: bool


def get_alignments(bam_path: str) -> dict[str, list[Alignment]]:
    """Extract primary and supplementary alignments from a BAM file.
    
    Returns a dict mapping read_name -> list of alignments.
    Secondary alignments are excluded.
    """
    alignments = {}
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            # Skip secondary alignments only
            if read.is_secondary:
                continue
            
            read_name = read.query_name
            
            if read_name not in alignments:
                alignments[read_name] = []
            
            if read.is_unmapped:
                # Only add unmapped if we have no other alignments for this read
                if not alignments[read_name]:
                    alignments[read_name].append(Alignment(
                        chrom="*",
                        pos=0,
                        end=0,
                        strand=".",
                        mapq=0,
                        cigar="*",
                        is_supplementary=False,
                        is_unmapped=True,
                    ))
            else:
                alignments[read_name].append(Alignment(
                    chrom=read.reference_name,
                    pos=read.reference_start,
                    end=read.reference_end,
                    strand="-" if read.is_reverse else "+",
                    mapq=read.mapping_quality,
                    cigar=read.cigarstring,
                    is_supplementary=read.is_supplementary,
                    is_unmapped=False,
                ))
    
    # Sort each read's alignments: primary first, then supplementary by position
    for read_name in alignments:
        alignments[read_name].sort(key=lambda a: (a.is_supplementary, a.chrom, a.pos))
    
    return alignments


def alignments_match(aln1: Alignment, aln2: Alignment, pos_tolerance: int) -> bool:
    """Check if two alignments match (same chrom, strand, overlapping position)."""
    if aln1.is_unmapped or aln2.is_unmapped:
        return aln1.is_unmapped and aln2.is_unmapped
    
    if aln1.chrom != aln2.chrom or aln1.strand != aln2.strand:
        return False
    
    # Check position overlap or within tolerance
    pos_diff = abs(aln1.pos - aln2.pos)
    return pos_diff <= pos_tolerance


def compare_alignments(
    alns1: list[Alignment],
    alns2: list[Alignment],
    pos_tolerance: int = 100,
) -> tuple[str, Optional[int], dict]:
    """
    Compare alignment sets for a read.
    
    Returns:
        (category, position_diff, details)
        
    Categories:
        - "both_unmapped": Neither tool mapped the read
        - "concordant": All alignments match between tools
        - "partial_match": Some alignments match, some don't
        - "discordant_chrom": Primary alignments on different chromosomes
        - "discordant_strand": Primary on same chrom, different strand
        - "discordant_pos": Primary on same chrom/strand, position differs
        - "only_1_mapped": Only first BAM mapped
        - "only_2_mapped": Only second BAM mapped
        - "different_count": Different number of alignments (chimeric vs non-chimeric)
    """
    # Handle unmapped cases
    unmapped1 = all(a.is_unmapped for a in alns1) if alns1 else True
    unmapped2 = all(a.is_unmapped for a in alns2) if alns2 else True
    
    details = {
        "n_alns1": len([a for a in alns1 if not a.is_unmapped]) if alns1 else 0,
        "n_alns2": len([a for a in alns2 if not a.is_unmapped]) if alns2 else 0,
    }
    
    if unmapped1 and unmapped2:
        return "both_unmapped", None, details
    
    if unmapped1:
        return "only_2_mapped", None, details
    
    if unmapped2:
        return "only_1_mapped", None, details
    
    # Both mapped - get non-unmapped alignments
    mapped1 = [a for a in alns1 if not a.is_unmapped]
    mapped2 = [a for a in alns2 if not a.is_unmapped]
    
    # Compare primary alignments first (first in the sorted list)
    primary1 = mapped1[0]
    primary2 = mapped2[0]
    
    if primary1.chrom != primary2.chrom:
        return "discordant_chrom", None, details
    
    if primary1.strand != primary2.strand:
        return "discordant_strand", None, details
    
    pos_diff = abs(primary1.pos - primary2.pos)
    
    if pos_diff > pos_tolerance:
        return "discordant_pos", pos_diff, details
    
    # Primary alignments match - check supplementary alignments
    if len(mapped1) != len(mapped2):
        return "different_count", pos_diff, details
    
    # Try to match all alignments
    matched2 = [False] * len(mapped2)
    all_matched = True
    
    for a1 in mapped1:
        found_match = False
        for j, a2 in enumerate(mapped2):
            if not matched2[j] and alignments_match(a1, a2, pos_tolerance):
                matched2[j] = True
                found_match = True
                break
        if not found_match:
            all_matched = False
            break
    
    if all_matched and all(matched2):
        return "concordant", pos_diff, details
    else:
        return "partial_match", pos_diff, details


def main():
    parser = argparse.ArgumentParser(
        description="Compare primary and supplementary alignments between two BAM files"
    )
    parser.add_argument("bam1", help="First BAM file (e.g., parallax)")
    parser.add_argument("bam2", help="Second BAM file (e.g., minimap2)")
    parser.add_argument(
        "--tolerance", "-t",
        type=int,
        default=100,
        help="Position tolerance for concordance (default: 100bp)"
    )
    parser.add_argument(
        "--details", "-d",
        help="Output TSV file with per-read details"
    )
    parser.add_argument(
        "--discordant-only",
        action="store_true",
        help="In details output, only include discordant reads"
    )
    parser.add_argument(
        "--name1",
        default="bam1",
        help="Label for first BAM in output (default: bam1)"
    )
    parser.add_argument(
        "--name2", 
        default="bam2",
        help="Label for second BAM in output (default: bam2)"
    )
    
    args = parser.parse_args()
    
    print(f"Loading {args.bam1}...", file=sys.stderr)
    alns1 = get_alignments(args.bam1)
    n_alns1 = sum(len(v) for v in alns1.values())
    print(f"  Found {len(alns1)} reads with {n_alns1} alignments", file=sys.stderr)
    
    print(f"Loading {args.bam2}...", file=sys.stderr)
    alns2 = get_alignments(args.bam2)
    n_alns2 = sum(len(v) for v in alns2.values())
    print(f"  Found {len(alns2)} reads with {n_alns2} alignments", file=sys.stderr)
    
    # Get all read names
    all_reads = set(alns1.keys()) | set(alns2.keys())
    print(f"Total unique reads: {len(all_reads)}", file=sys.stderr)
    
    # Compare
    categories = defaultdict(int)
    pos_diffs = []  # For concordant alignments
    details = []
    
    for read_name in sorted(all_reads):
        read_alns1 = alns1.get(read_name, [])
        read_alns2 = alns2.get(read_name, [])
        
        if not read_alns1:
            category = f"only_in_{args.name2}"
            pos_diff = None
            cmp_details = {"n_alns1": 0, "n_alns2": len(read_alns2)}
        elif not read_alns2:
            category = f"only_in_{args.name1}"
            pos_diff = None
            cmp_details = {"n_alns1": len(read_alns1), "n_alns2": 0}
        else:
            category, pos_diff, cmp_details = compare_alignments(read_alns1, read_alns2, args.tolerance)
        
        categories[category] += 1
        
        if category == "concordant" and pos_diff is not None:
            pos_diffs.append(pos_diff)
        
        # Collect details if requested
        if args.details:
            if not args.discordant_only or category not in ("concordant", "both_unmapped"):
                # Get primary alignment info (first non-unmapped, or first)
                primary1 = next((a for a in read_alns1 if not a.is_unmapped), None) if read_alns1 else None
                primary2 = next((a for a in read_alns2 if not a.is_unmapped), None) if read_alns2 else None
                
                details.append({
                    "read_name": read_name,
                    "category": category,
                    "pos_diff": pos_diff if pos_diff is not None else "",
                    "n_alns1": cmp_details["n_alns1"],
                    "n_alns2": cmp_details["n_alns2"],
                    "chrom1": primary1.chrom if primary1 else "",
                    "pos1": primary1.pos if primary1 else "",
                    "strand1": primary1.strand if primary1 else "",
                    "mapq1": primary1.mapq if primary1 else "",
                    "chrom2": primary2.chrom if primary2 else "",
                    "pos2": primary2.pos if primary2 else "",
                    "strand2": primary2.strand if primary2 else "",
                    "mapq2": primary2.mapq if primary2 else "",
                })
    
    # Print summary
    print("\n" + "=" * 60)
    print("ALIGNMENT COMPARISON SUMMARY")
    print("=" * 60)
    print(f"{args.name1}: {args.bam1}")
    print(f"{args.name2}: {args.bam2}")
    print(f"Position tolerance: {args.tolerance}bp")
    print("-" * 60)
    
    total = len(all_reads)
    
    # Order categories nicely
    category_order = [
        "concordant",
        "partial_match",
        "different_count",
        "discordant_pos",
        "discordant_strand",
        "discordant_chrom",
        "only_1_mapped",
        "only_2_mapped",
        f"only_in_{args.name1}",
        f"only_in_{args.name2}",
        "both_unmapped",
    ]
    
    for cat in category_order:
        if cat in categories:
            count = categories[cat]
            pct = 100.0 * count / total
            # Make labels more readable
            label = cat.replace("only_1_mapped", f"only_{args.name1}_mapped")
            label = label.replace("only_2_mapped", f"only_{args.name2}_mapped")
            print(f"{label:30s} {count:8d} ({pct:5.1f}%)")
    
    # Print any remaining categories not in the order list
    for cat in sorted(categories.keys()):
        if cat not in category_order:
            count = categories[cat]
            pct = 100.0 * count / total
            print(f"{cat:30s} {count:8d} ({pct:5.1f}%)")
    
    print("-" * 60)
    print(f"{'TOTAL':30s} {total:8d}")
    
    # Position diff stats for concordant
    if pos_diffs:
        import statistics
        print("\nConcordant position differences:")
        print(f"  Mean:   {statistics.mean(pos_diffs):.1f}bp")
        print(f"  Median: {statistics.median(pos_diffs):.1f}bp")
        print(f"  Max:    {max(pos_diffs)}bp")
        
        # Histogram buckets
        buckets = [0, 1, 5, 10, 20, 50, 100]
        print("\n  Distribution:")
        for i, threshold in enumerate(buckets):
            count = sum(1 for d in pos_diffs if d <= threshold)
            pct = 100.0 * count / len(pos_diffs)
            print(f"    <= {threshold:3d}bp: {count:8d} ({pct:5.1f}%)")
    
    # Write details if requested
    if args.details:
        print(f"\nWriting details to {args.details}...", file=sys.stderr)
        with open(args.details, "w") as f:
            header = ["read_name", "category", "pos_diff",
                      f"n_alns_{args.name1}", f"n_alns_{args.name2}",
                      f"chrom_{args.name1}", f"pos_{args.name1}", 
                      f"strand_{args.name1}", f"mapq_{args.name1}",
                      f"chrom_{args.name2}", f"pos_{args.name2}",
                      f"strand_{args.name2}", f"mapq_{args.name2}"]
            f.write("\t".join(header) + "\n")
            for d in details:
                row = [
                    d["read_name"], d["category"], str(d["pos_diff"]),
                    str(d["n_alns1"]), str(d["n_alns2"]),
                    d["chrom1"], str(d["pos1"]), d["strand1"], str(d["mapq1"]),
                    d["chrom2"], str(d["pos2"]), d["strand2"], str(d["mapq2"]),
                ]
                f.write("\t".join(row) + "\n")
        print(f"  Wrote {len(details)} records", file=sys.stderr)


if __name__ == "__main__":
    main()
