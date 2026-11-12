"""Make the read variants used for invariance testing.

Usage:
    make_variants.py <reads.fq[.gz]> [--rc OUT] [--frag LEN:OUT ...]

--rc OUT        write the reverse complement of every read (same names)
--frag LEN:OUT  cut each read into non-overlapping tiles of LEN bases, named
                <read>/<start>-<end> (forward-strand coordinates); a trailing
                tile is kept if it is at least LEN/2 long

Outputs are gzipped FASTQ.  Read names are truncated at the first whitespace.
"""

import argparse
import gzip
import sys

from common import open_text

COMPLEMENT = bytes.maketrans(b"ACGTNacgtnRYKMBVDHrykmbvdh", b"TGCANtgcanYRMKVBHDyrmkvbhd")


def fastq_records(path):
    with open_text(path, "rb") as f:
        while True:
            header = f.readline()
            if not header:
                return
            seq = f.readline().rstrip(b"\r\n")
            f.readline()
            qual = f.readline().rstrip(b"\r\n")
            if not header.startswith(b"@") or len(seq) != len(qual):
                sys.exit(f"malformed FASTQ record near {header!r}")
            yield header[1:].split(None, 1)[0], seq, qual


def write(out, name, seq, qual):
    out.write(b"@" + name + b"\n" + seq + b"\n+\n" + qual + b"\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reads")
    ap.add_argument("--rc")
    ap.add_argument("--frag", action="append", default=[])
    args = ap.parse_args()

    rc_out = gzip.open(args.rc, "wb", compresslevel=1) if args.rc else None
    frags = []
    for spec in args.frag:
        length, path = spec.split(":", 1)
        frags.append((int(length), gzip.open(path, "wb", compresslevel=1)))

    for name, seq, qual in fastq_records(args.reads):
        if rc_out:
            write(rc_out, name, seq.translate(COMPLEMENT)[::-1], qual[::-1])
        for length, out in frags:
            for start in range(0, len(seq), length):
                end = min(start + length, len(seq))
                if end - start < length and 2 * (end - start) < length:
                    break
                tile = name + b"/" + f"{start}-{end}".encode()
                write(out, tile, seq[start:end], qual[start:end])

    if rc_out:
        rc_out.close()
    for _, out in frags:
        out.close()


if __name__ == "__main__":
    main()
