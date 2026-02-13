#!/usr/bin/env python3
"""
Check if alignments of simulated reads match their expected positions.

Reads a SAM file where read names encode the true origin:
    sim_XXXXXXXX:chrom_start_end_strand[:errors]

For example:
    sim_00000001:chr1_1000_2000_+
    sim_00000002:chr1_5000_6000_-:A50T,G100C

Reports primary alignments that don't match the expected position.
Also checks if any secondary/supplementary alignments match when primary doesn't.
"""

import argparse
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional, List


@dataclass
class ExpectedAlignment:
    """Expected alignment parsed from read name."""
    chrom: str
    start: int  # 1-based
    end: int    # 1-based, inclusive
    strand: str  # '+' or '-'
    errors: Optional[str] = None


@dataclass 
class SamAlignment:
    """Alignment from SAM file."""
    qname: str
    flag: int
    rname: str
    pos: int  # 1-based
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


def parse_read_name(qname: str) -> Optional[ExpectedAlignment]:
    """Parse expected alignment from simulated read name."""
    # Format: sim_XXXXXXXX:chrom_start_end_strand[:errors]
    # Example: sim_00000001:chr1_1000_2000_+
    
    # Handle potential errors suffix
    parts = qname.split(':')
    if len(parts) < 2:
        return None
    
    # The location part is: chrom_start_end_strand
    # But chrom can contain underscores (e.g., chr1_random)
    # So we parse from the end
    loc_part = parts[1]
    if len(parts) > 2:
        errors = ':'.join(parts[2:])
    else:
        errors = None
    
    # Parse strand (last character after last _)
    match = re.match(r'^(.+)_(\d+)_(\d+)_([+-])$', loc_part)
    if not match:
        return None
    
    chrom = match.group(1)
    start = int(match.group(2))
    end = int(match.group(3))
    strand = match.group(4)
    
    return ExpectedAlignment(chrom, start, end, strand, errors)


def parse_cigar_ref_length(cigar: str) -> int:
    """Calculate reference length consumed by CIGAR."""
    length = 0
    for match in re.finditer(r'(\d+)([MIDNSHP=X])', cigar):
        count = int(match.group(1))
        op = match.group(2)
        if op in 'MDN=X':  # Operations that consume reference
            length += count
    return length


def check_alignment(expected: ExpectedAlignment, actual: SamAlignment, tolerance: int) -> tuple[bool, str]:
    """
    Check if alignment matches expected position.
    
    Returns (matches, reason) tuple.
    """
    # Check chromosome
    if expected.chrom != actual.rname:
        return False, f"chrom mismatch: expected {expected.chrom}, got {actual.rname}"
    
    # Check strand
    expected_reverse = (expected.strand == '-')
    if expected_reverse != actual.is_reverse:
        return False, f"strand mismatch: expected {expected.strand}, got {'-' if actual.is_reverse else '+'}"
    
    # Check position
    # The expected start/end are 1-based inclusive coordinates of the original genomic region
    # For forward strand: SAM POS should equal expected start
    # For reverse strand: SAM POS should also equal expected start (the leftmost position)
    
    ref_len = parse_cigar_ref_length(actual.cigar)
    actual_end = actual.pos + ref_len - 1  # 1-based inclusive
    
    # Allow some tolerance for soft-clipping
    start_diff = abs(actual.pos - expected.start)
    end_diff = abs(actual_end - expected.end)
    
    if start_diff > tolerance or end_diff > tolerance:
        return False, f"position mismatch: expected {expected.start}-{expected.end}, got {actual.pos}-{actual_end} (diff: start={start_diff}, end={end_diff})"
    
    return True, "ok"


def main():
    parser = argparse.ArgumentParser(
        description="Check simulated read alignments against expected positions"
    )
    parser.add_argument("sam_file", help="Input SAM file (use - for stdin)")
    parser.add_argument(
        "-t", "--tolerance", type=int, default=50,
        help="Position tolerance in bp (default: 50)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Print all alignments, not just mismatches"
    )
    parser.add_argument(
        "-s", "--summary", action="store_true",
        help="Print summary statistics"
    )
    parser.add_argument(
        "--include-secondary", action="store_true",
        help="Also check secondary/supplementary alignments"
    )
    args = parser.parse_args()
    
    # Open input
    if args.sam_file == "-":
        infile = sys.stdin
    else:
        infile = open(args.sam_file, 'r')
    
    # Collect all alignments per read
    alignments_by_read: dict[str, List[SamAlignment]] = defaultdict(list)
    
    try:
        for line in infile:
            line = line.strip()
            
            # Skip header lines
            if line.startswith('@'):
                continue
            
            fields = line.split('\t')
            if len(fields) < 11:
                continue
            
            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3])
            mapq = int(fields[4])
            cigar = fields[5]
            seq = fields[9]
            
            alignment = SamAlignment(
                qname=qname,
                flag=flag,
                rname=rname,
                pos=pos,
                mapq=mapq,
                cigar=cigar,
                seq_len=len(seq)
            )
            
            alignments_by_read[qname].append(alignment)
    
    finally:
        if infile != sys.stdin:
            infile.close()
    
    # Statistics
    total = 0
    matched = 0
    mismatched = 0
    mismatched_but_secondary_ok = 0
    unmapped = 0
    skipped = 0
    
    # Process each read
    for qname, alns in alignments_by_read.items():
        # Parse expected position from read name
        expected = parse_read_name(qname)
        if expected is None:
            skipped += 1
            continue
        
        total += 1
        
        # Find primary alignment
        primary = None
        secondaries: List[SamAlignment] = []
        for aln in alns:
            if aln.is_primary:
                primary = aln
            elif not aln.is_unmapped:
                secondaries.append(aln)
        
        # Check if unmapped
        if primary is None or primary.is_unmapped:
            unmapped += 1
            print(f"UNMAPPED\t{qname}\texpected {expected.chrom}:{expected.start}-{expected.end}{expected.strand}")
            continue
        
        primary_matches, primary_reason = check_alignment(expected, primary, args.tolerance)
        
        if primary_matches:
            matched += 1
            if args.verbose:
                print(f"OK\t{qname}\t{primary.rname}:{primary.pos}\tMAPQ={primary.mapq}")
        else:
            # Check if any secondary/supplementary alignment matches
            matching_secondary = None
            for sec in secondaries:
                sec_matches, _ = check_alignment(expected, sec, args.tolerance)
                if sec_matches:
                    matching_secondary = sec
                    break
            
            mismatched += 1
            if matching_secondary:
                mismatched_but_secondary_ok += 1
                sec_type = "supplementary" if matching_secondary.is_supplementary else "secondary"
                print(f"MISMATCH_BUT_SECONDARY_OK\t{qname}\tprimary={primary.rname}:{primary.pos} MAPQ={primary.mapq}\t"
                      f"{sec_type}={matching_secondary.rname}:{matching_secondary.pos} MAPQ={matching_secondary.mapq}\t{primary_reason}")
            else:
                print(f"MISMATCH\t{qname}\t{primary.rname}:{primary.pos}\tMAPQ={primary.mapq}\t{primary_reason}")
                if args.include_secondary and secondaries:
                    for sec in secondaries:
                        _, sec_reason = check_alignment(expected, sec, args.tolerance)
                        sec_type = "supplementary" if sec.is_supplementary else "secondary"
                        print(f"  {sec_type}: {sec.rname}:{sec.pos} MAPQ={sec.mapq} - {sec_reason}")
    
    # Print summary
    if args.summary or total > 0:
        print(f"\n--- Summary ---", file=sys.stderr)
        print(f"Total checked: {total}", file=sys.stderr)
        print(f"Matched: {matched} ({100*matched/total:.1f}%)" if total > 0 else "Matched: 0", file=sys.stderr)
        print(f"Mismatched: {mismatched} ({100*mismatched/total:.1f}%)" if total > 0 else "Mismatched: 0", file=sys.stderr)
        if mismatched_but_secondary_ok > 0:
            print(f"  - Primary wrong but secondary/supplementary correct: {mismatched_but_secondary_ok} ({100*mismatched_but_secondary_ok/mismatched:.1f}% of mismatches)", file=sys.stderr)
        print(f"Unmapped: {unmapped}", file=sys.stderr)
        if skipped > 0:
            print(f"Skipped (unparseable names): {skipped}", file=sys.stderr)


if __name__ == "__main__":
    main()
