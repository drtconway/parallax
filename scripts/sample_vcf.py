#!/usr/bin/env python3
"""Sample variants from a VCF based on the AF info field probability."""

import argparse
import gzip
import random
import sys

alu = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGTGGTAGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
line1 = "AAAAATAGAAACTATACTAACACAAATCAAAAGAAAGCTGGGGTAGCTATATTAATTTCAGACAAAGCTGACTTCAGAAGGAAAATTGTCAAAGGCATTACTTAATGAGAAGAGCTCTATCCTCCAGGAAGACATAACAATCCTTAATGTGTATGTGCCTAACAAGAGAGTGTCAAAATACAGAGACAAAAACTAGTAGAAATGCAAGGAGAAATAAACAATGCCATTATTATAGTTGGAGACTTCAGCACACCTTTATCAATAATTGACAGATCTTGCAGGCAGAAAATCAGTAAAAATAGTTTAACTAAACAGAACCATCAGTAAACTGATTTAATTGATGTTTACAGAATACTTCATCAAACAACAGCAGAATATGCATTATTCTTAAGCTCATATGGAACAGTCACCATGACAGACAACATGCTGGACTATTACACATACCTCAACAAATTTAAAGGACTAAAATTGACACAAAATATGCCCGAGGACACAGTGAAATTAAACTAGAAATCAATACCAAGAAGACAGCTGGAAAACCCCAATGTATTTGAGATTGAACAACATAATTCTAAATAACACATGGCTCAAAGAGGAAAACTTAAAGATATTGTGAAGTGT"

def smart_open(filename: str, mode: str = "rt"):
    """Open a file, using gzip if it ends with .gz or .bgz."""
    if filename.endswith(".gz") or filename.endswith(".bgz"):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def parse_af(info: str) -> float | None:
    """Extract AF value from the INFO field."""
    for field in info.split(";"):
        if field.startswith("AF="):
            # AF can be comma-separated for multi-allelic; take the first
            return float(field[3:].split(",")[0])
    return None

class ReservoirSampler:
    """Reservoir sampler for streaming data."""
    def __init__(self, k: int, rng: random.Random):
        self.k = k
        self.rng = rng
        self.sample = []
        self.n = 0

    def add(self, item):
        self.n += 1
        if len(self.sample) < self.k:
            self.sample.append(item)
        else:
            s = self.rng.randint(0, self.n - 1)
            # Weighted reservoir sampling: replace with
            # probability proportional to the relative AF
            if s < self.k:
                u = self.rng.random()
                v = self.rng.random()
                if item[1] * u > self.sample[s][1] * v:
                    self.sample[s] = item

def main():
    parser = argparse.ArgumentParser(
        description="Sample VCF events using the AF info field probability."
    )
    parser.add_argument("input_vcf", help="Input VCF file (- for stdin)")
    parser.add_argument(
        "-o", "--output", default="-", help="Output VCF file (default: stdout)"
    )
    parser.add_argument(
        "-s", "--seed", type=int, default=42, help="Random seed"
    )
    parser.add_argument(
        "-n", "--num-variants", type=int, default=1000,
        help="Number of variants to sample (default: 1000)"
    )
    args = parser.parse_args()

    rng = random.Random(args.seed)
    n = args.num_variants
    samplers = {}

    inp = sys.stdin if args.input_vcf == "-" else smart_open(args.input_vcf)
    out = sys.stdout if args.output == "-" else smart_open(args.output, "w")

    try:
        prev_chrom = None

        for vn, line in enumerate(inp):
            if line.startswith("#"):
                out.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0]
            if chrom != prev_chrom:
                prev_chrom = chrom
                print(f"Processing chromosome {chrom}...", file=sys.stderr)

            alt = fields[4]
            info = fields[7]
            seq = None

            # Skip BND events
            if alt == "<BND>" or "SVTYPE=BND" in info:
                continue

            # Skip CNV events
            if alt == "<CNV>" or "SVTYPE=CNV" in info:
                continue

            # Skip CPX events
            if alt == "<CPX>" or "SVTYPE=CPX" in info:
                continue

            # Skip INS events (common, but we don't have the inserted sequence)
            if alt == "<INS>" or "SVTYPE=INS" in info:
                continue

            if alt == "<INS:ME:ALU>":
                seq = alu
            elif alt == "<INS:ME:LINE1>":
                seq = line1

            af = parse_af(info)
            if af is None:
                # No AF field — include the variant as-is
                out.write(line)
                continue

        if seq is None:
            fields[7] = "."
        else:
            # Add some random purturbation to the sequence to make it more realistic
            seq = "".join(rng.choice("ACGT") if rng.random() < 0.05 else c for c in seq)
            fields[7] = f"SVLEN={len(seq)};SEQ={seq}"
        line = "\t".join(fields) + "\n"
        if alt not in samplers:
            samplers[alt] = ReservoirSampler(n, rng)
        samplers[alt].add((vn, af, line))

        # Collect all sampled variants of each kind
        lines = []
        for _alt, sampler in samplers.items():
             lines.extend(sampler.sample)

        # Write the sampled variants
        for _, _, line in sorted(lines, key=lambda x: x[0]):
            out.write(line)

    finally:
        if inp is not sys.stdin:
            inp.close()
        if out is not sys.stdout:
            out.close()



if __name__ == "__main__":
    main()
