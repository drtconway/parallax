#!/usr/bin/env python3
"""
Convert dbVar structural variant CSV to well-formed VCF 4.3.

Reads the CSV format exported from dbVar (NCBI Database of Genomic Structural
Variation) and produces a sorted VCF with symbolic alleles for structural
variants.

Usage:
    python sv_csv_to_vcf.py input.csv -o output.vcf
    python sv_csv_to_vcf.py input.csv -o output.vcf -r reference.fa
    python sv_csv_to_vcf.py input.csv  # writes to stdout
"""

import argparse
import csv
from collections import namedtuple
from fileinput import filename
import sys
from datetime import date

import gzip

def smart_open(filename: str, mode: str = "rt"):
    """Open a file, using gzip if it ends with .gz or .bgz."""
    if filename.endswith(".gz") or filename.endswith(".bgz"):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


SV_TYPE_MAP = {
    "duplication": "DUP",
    "deletion": "DEL",
    "copy number variation": "CNV",
    "copy number gain": "DUP",
    "copy number loss": "DEL",
    "insertion": "INS",
    "inversion": "INV",
    "tandem duplication": "DUP",
    "mobile element insertion": "INS:ME",
    "mobile element deletion": "DEL:ME",
    "alu insertion": "INS:ME:ALU",
    "alu deletion": "DEL:ME:ALU",
    "line1 insertion": "INS:ME:LINE1",
    "line1 deletion": "DEL:ME:LINE1",
    "sva insertion": "INS:ME:SVA",
    "sva deletion": "DEL:ME:SVA",
    "herv deletion": "DEL:ME:HERV",
}

VcfRecord = namedtuple("VcfRecord", ["chrom", "pos", "id", "ref", "alt", "info", "chrom_raw"])

# Chromosome sort order: 1-22, X, Y, MT, then anything else
CHROM_ORDER = {}
for i in range(1, 23):
    CHROM_ORDER[str(i)] = i
CHROM_ORDER["X"] = 23
CHROM_ORDER["Y"] = 24
CHROM_ORDER["M"] = 25
CHROM_ORDER["MT"] = 25


def parse_int(value):
    """Parse an integer from a possibly empty/quoted string."""
    if value is None:
        return None
    value = value.strip().strip('"')
    if not value:
        return None
    return int(value)


def get_sv_type(call_type):
    """Map variant call type string to VCF SVTYPE."""
    return SV_TYPE_MAP.get(call_type.strip().strip('"').lower(), "CNV")


def load_reference(ref_path):
    """Load a reference FASTA using pysam. Returns None if unavailable."""
    try:
        import pysam
        return pysam.FastaFile(ref_path)
    except ImportError:
        print(
            "Warning: pysam not installed; using 'N' for REF bases. "
            "Install with: pip install pysam",
            file=sys.stderr,
        )
        return None
    except Exception as e:
        print(f"Warning: could not open reference: {e}", file=sys.stderr)
        return None


def get_ref_base(ref_fasta, chrom, pos):
    """Fetch the reference base at a 1-based position."""
    if ref_fasta is None:
        return "N"
    try:
        for name in [f"chr{chrom}", chrom, f"Chr{chrom}"]:
            if name in ref_fasta.references:
                base = ref_fasta.fetch(name, pos - 1, pos).upper()
                return base if base else "N"
    except Exception:
        pass
    return "N"


def determine_positions(row):
    """
    Determine POS, END, CIPOS, and CIEND from the six coordinate columns.

    dbVar coordinate model:
        Outer Start <= Start <= Inner Start  ... variant ...  Inner End <= End <= Outer End

    VCF model:
        POS with CIPOS = (outer_start - pos, inner_start - pos)
        END with CIEND = (inner_end - end, outer_end - end)
    """
    outer_start = parse_int(row.get("Outer Start", ""))
    start = parse_int(row.get("Start", ""))
    inner_start = parse_int(row.get("Inner Start", ""))
    inner_end = parse_int(row.get("Inner End", ""))
    end = parse_int(row.get("End", ""))
    outer_end = parse_int(row.get("Outer End", ""))

    # Best estimate of start: prefer Start, then Inner Start, then Outer Start
    pos = start if start is not None else (
        inner_start if inner_start is not None else outer_start
    )
    # Best estimate of end: prefer End, then Inner End, then Outer End
    end_pos = end if end is not None else (
        inner_end if inner_end is not None else outer_end
    )

    # Confidence interval around POS
    cipos = None
    if pos is not None:
        lo = (outer_start - pos) if outer_start is not None else None
        hi = (inner_start - pos) if inner_start is not None else None
        if lo is not None or hi is not None:
            cipos = (lo if lo is not None else 0, hi if hi is not None else 0)

    # Confidence interval around END
    ciend = None
    if end_pos is not None:
        lo = (inner_end - end_pos) if inner_end is not None else None
        hi = (outer_end - end_pos) if outer_end is not None else None
        if lo is not None or hi is not None:
            ciend = (lo if lo is not None else 0, hi if hi is not None else 0)

    return pos, end_pos, cipos, ciend


def write_vcf_header(out, assembly, source_file, chroms):
    """Write a compliant VCF 4.3 header."""
    out.write("##fileformat=VCFv4.3\n")
    out.write(f"##fileDate={date.today().strftime('%Y%m%d')}\n")
    out.write(f"##source=sv_csv_to_vcf.py\n")
    out.write(f"##inputFile={source_file}\n")
    if assembly:
        out.write(f"##reference={assembly}\n")

    # Contig lines for chromosomes seen in the data
    for c in sorted(chroms, key=lambda x: CHROM_ORDER.get(x.replace("chr", ""), 99)):
        out.write(f"##contig=<ID={c}>\n")

    # ALT descriptions
    out.write('##ALT=<ID=DEL,Description="Deletion">\n')
    out.write('##ALT=<ID=DUP,Description="Duplication">\n')
    out.write('##ALT=<ID=CNV,Description="Copy number variation">\n')
    out.write('##ALT=<ID=INS,Description="Insertion">\n')
    out.write('##ALT=<ID=INV,Description="Inversion">\n')
    out.write('##ALT=<ID=INS:ME,Description="Mobile element insertion">\n')
    out.write('##ALT=<ID=DEL:ME,Description="Mobile element deletion">\n')
    out.write('##ALT=<ID=INS:ME:ALU,Description="Alu element insertion">\n')
    out.write('##ALT=<ID=DEL:ME:ALU,Description="Alu element deletion">\n')
    out.write('##ALT=<ID=INS:ME:LINE1,Description="LINE1 element insertion">\n')
    out.write('##ALT=<ID=DEL:ME:LINE1,Description="LINE1 element deletion">\n')
    out.write('##ALT=<ID=INS:ME:SVA,Description="SVA element insertion">\n')
    out.write('##ALT=<ID=DEL:ME:SVA,Description="SVA element deletion">\n')
    out.write('##ALT=<ID=DEL:ME:HERV,Description="HERV element deletion">\n')

    # INFO fields
    out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n')
    out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the structural variant">\n')
    out.write('##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">\n')
    out.write('##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">\n')
    out.write('##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">\n')
    out.write('##INFO=<ID=DBVARID,Number=1,Type=String,Description="dbVar study accession">\n')
    out.write('##INFO=<ID=CALLTYPE,Number=1,Type=String,Description="Original variant call type">\n')

    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


def chrom_label(chrom_raw):
    """Normalise chromosome name to 'chrN' form."""
    c = chrom_raw.strip().strip('"')
    if not c.startswith("chr"):
        c = f"chr{c}"
    return c


def convert_row(row, ref_fasta):
    """Convert one CSV row into a VCF record dict (or None to skip)."""
    chrom_raw = row["Chromosome"].strip().strip('"')
    variant_id = row["Variant ID"].strip().strip('"')
    call_type = row["Variant Call type"].strip().strip('"')
    study_id = row["Study ID"].strip().strip('"')

    svtype = get_sv_type(call_type)
    pos, end_pos, cipos, ciend = determine_positions(row)

    if pos is None:
        print(f"Warning: skipping {variant_id} — no usable start coordinate", file=sys.stderr)
        return None

    vcf_chrom = chrom_label(chrom_raw)
    ref_base = get_ref_base(ref_fasta, chrom_raw, pos)
    alt = f"<{svtype}>"

    # Determine whether the call is imprecise
    imprecise = (cipos is not None and cipos != (0, 0)) or (ciend is not None and ciend != (0, 0))

    # Build INFO
    info_parts = [f"SVTYPE={svtype}"]
    if imprecise:
        info_parts.append("IMPRECISE")
    if end_pos is not None:
        info_parts.append(f"END={end_pos}")
        svlen = end_pos - pos
        if svtype.startswith("DEL"):
            svlen = -svlen
        info_parts.append(f"SVLEN={svlen}")
    if cipos is not None and cipos != (0, 0):
        info_parts.append(f"CIPOS={cipos[0]},{cipos[1]}")
    if ciend is not None and ciend != (0, 0):
        info_parts.append(f"CIEND={ciend[0]},{ciend[1]}")
    if study_id:
        info_parts.append(f"DBVARID={study_id}")
    if call_type:
        info_parts.append(f"CALLTYPE={call_type.replace(' ', '_')}")

    info = ";".join(info_parts)

    return VcfRecord(
        chrom=vcf_chrom,
        pos=pos,
        id=variant_id,
        ref=ref_base,
        alt=alt,
        info=info,
        chrom_raw=chrom_raw,
    )


def sort_key(rec):
    """Sort records by chromosome then position."""
    return (CHROM_ORDER.get(rec.chrom_raw, 99), rec.pos)


def main():
    parser = argparse.ArgumentParser(
        description="Convert dbVar structural variant CSV to VCF 4.3 format",
    )
    parser.add_argument("input", help="Input CSV file from dbVar")
    parser.add_argument(
        "-o", "--output", default="-",
        help="Output VCF file (default: stdout)",
    )
    parser.add_argument(
        "-r", "--reference",
        help="Reference FASTA (indexed) for REF base lookup; if omitted, 'N' is used",
    )
    args = parser.parse_args()

    ref_fasta = load_reference(args.reference) if args.reference else None

    # Read all rows
    with smart_open(args.input) as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)

    if not rows:
        print("Error: input CSV is empty", file=sys.stderr)
        sys.exit(1)

    # Convert rows, deduplicating on all fields except variant ID
    records = {}
    converted = 0
    for row in rows:
        rec = convert_row(row, ref_fasta)
        if rec is not None:
            converted += 1
            key = (rec.chrom, rec.pos, rec.ref, rec.alt, rec.info)
            if key not in records:
                records[key] = rec

    if not records:
        print("Error: no valid records produced", file=sys.stderr)
        sys.exit(1)

    # Sort by chromosome and position
    sorted_records = sorted(records.values(), key=sort_key)

    # Collect unique chromosomes in sort order
    seen = set()
    chroms = []
    for r in sorted_records:
        if r.chrom not in seen:
            seen.add(r.chrom)
            chroms.append(r.chrom)

    # Determine assembly from first row
    assembly = rows[0].get("Assembly Name", "").strip().strip('"')

    # Write output
    out = sys.stdout if args.output == "-" else smart_open(args.output, "w")
    try:
        write_vcf_header(out, assembly, args.input, chroms)
        for rec in sorted_records:
            out.write(
                f"{rec.chrom}\t{rec.pos}\t{rec.id}\t{rec.ref}\t"
                f"{rec.alt}\t.\t.\t{rec.info}\n"
            )
    finally:
        if out is not sys.stdout:
            out.close()

    if ref_fasta is not None:
        ref_fasta.close()

    n = len(records)
    deduped = converted - n
    skipped = len(rows) - converted
    print(f"Wrote {n} records to VCF", file=sys.stderr)
    if deduped:
        print(f"Deduplicated {deduped} duplicate rows", file=sys.stderr)
    if skipped:
        print(f"Skipped {skipped} rows (no usable coordinates)", file=sys.stderr)


if __name__ == "__main__":
    main()
