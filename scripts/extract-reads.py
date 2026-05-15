#!/usr/bin/env python3
"""Extract named reads from a FASTQ stream.

Usage:
    python extract-reads.py names.txt < input.fastq > output.fastq
"""

import sys

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} names.txt < input.fastq", file=sys.stderr)
        sys.exit(1)

    with open(sys.argv[1]) as f:
        names = {line.strip() for line in f if line.strip()}

    found = 0
    while True:
        header = sys.stdin.readline()
        if not header:
            break
        seq = sys.stdin.readline()
        plus = sys.stdin.readline()
        qual = sys.stdin.readline()

        # Read name is everything after '@' up to the first whitespace
        read_name = header[1:].split()[0]
        if read_name in names:
            sys.stdout.write(header + seq + plus + qual)
            found += 1
            if found == len(names):
                break

    print(f"Extracted {found}/{len(names)} reads", file=sys.stderr)

if __name__ == "__main__":
    main()
