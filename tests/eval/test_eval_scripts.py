"""Tests for scripts/eval (run with: .venv/bin/python3 -m pytest tests/eval)."""

import gzip
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts" / "eval"
sys.path.insert(0, str(SCRIPTS))

import common  # noqa: E402

HEADER = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": "chr1", "LN": 10000}]}


def make_bam(path, records):
    """records: (name, flag, pos, cigar, mapq, seq) tuples."""
    with pysam.AlignmentFile(path, "wb", header=HEADER) as out:
        for name, flag, pos, cigar, mapq, seq in records:
            a = pysam.AlignedSegment(out.header)
            a.query_name = name
            a.flag = flag
            if not flag & 4:
                a.reference_id = 0
                a.reference_start = pos
                a.cigarstring = cigar
                a.mapping_quality = mapq
            a.query_sequence = seq
            out.write(a)
    return path


def read_tsv(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f]


def run(script, *args):
    subprocess.run([sys.executable, str(SCRIPTS / script), *map(str, args)], check=True)


def test_query_range_reverse_and_hard_clips(tmp_path):
    bam = make_bam(tmp_path / "a.bam", [
        ("r", 0, 100, "5S20=10H", 60, "A" * 25),
        ("r", 16 | 2048, 500, "10H20=5S", 60, "A" * 25),
    ])
    with pysam.AlignmentFile(bam) as f:
        fwd, rev = list(f)
    assert common.read_length(fwd) == 35
    assert common.query_range(fwd) == (5, 25)
    # Reverse: alignment-orientation [10, 30) -> forward [35-30, 35-10).
    assert common.query_range(rev) == (5, 25)


def test_rc_invariance_identical_and_disjoint(tmp_path):
    seq = "ACGT" * 25
    orig = make_bam(tmp_path / "orig.bam", [
        ("same", 0, 100, "100=", 60, seq),
        ("moved", 0, 100, "100=", 60, seq),
        ("lost", 0, 100, "100=", 60, seq),
    ])
    rc = make_bam(tmp_path / "rc.bam", [
        ("same", 16, 100, "100=", 60, seq),
        ("moved", 16, 5000, "100=", 60, seq),
        ("lost", 4, 0, None, 0, seq),
    ])
    out = tmp_path / "rc.tsv.gz"
    run("invariance.py", orig, rc, "--kind", "rc", "--out", out, "--step", 10)
    got = {r["read"]: r["class"] for r in read_tsv(out)}
    assert got == {"same": "identical", "moved": "disjoint", "lost": "orig_only"}


def test_rc_invariance_detects_shifted_indel(tmp_path):
    seq = "A" * 100
    # Same span, but the 1 bp deletion is placed after base 40 in one and
    # after base 60 in the other, so bases 40..59 land one position apart.
    orig = make_bam(tmp_path / "orig.bam", [("r", 0, 100, "40=1D60=", 60, seq)])
    rc = make_bam(tmp_path / "rc.bam", [("r", 16, 100, "60=1D40=", 60, seq)])
    out = tmp_path / "rc.tsv.gz"
    run("invariance.py", orig, rc, "--kind", "rc", "--out", out, "--step", 1)
    (row,) = read_tsv(out)
    assert row["class"] == "partial"
    assert float(row["score"]) == pytest.approx(0.8)


def test_fragment_invariance(tmp_path):
    seq = "ACGT" * 50
    orig = make_bam(tmp_path / "orig.bam", [("p", 0, 1000, "200=", 60, seq)])
    tiles = make_bam(tmp_path / "tiles.bam", [
        ("p/0-100", 0, 1000, "100=", 60, seq[:100]),
        ("p/100-200", 0, 7000, "100=", 60, seq[100:]),
    ])
    out = tmp_path / "frag.tsv.gz"
    run("invariance.py", orig, tiles, "--kind", "frag", "--out", out, "--step", 10)
    got = {r["tile"]: r["class"] for r in read_tsv(out)}
    assert got == {"p/0-100": "agree", "p/100-200": "disagree_confident"}


def test_alignment_stats_counts_and_split_indel(tmp_path):
    seq = "A" * 100
    bam = make_bam(tmp_path / "a.bam", [
        # Primary covers query [0, 48): 10= 1X 5= 2I 10= 3D 20=
        ("r", 0, 1000, "10=1X5=2I10=3D20=52S", 60, seq),
        # Supplementary covers [48, 100), 500 bp further on: a split deletion.
        ("r", 2048, 1049 + 500, "48H52=", 60, seq[48:]),
    ])
    prefix = tmp_path / "x"
    run("alignment_stats.py", bam, "--prefix", prefix)
    alns = read_tsv(f"{prefix}.alignments.tsv.gz")
    prim = next(a for a in alns if a["class"] == "primary")
    assert (prim["matches"], prim["mismatches"]) == ("45", "1")
    assert (prim["ins_events"], prim["ins_bases"], prim["del_events"], prim["del_bases"]) == ("1", "2", "1", "3")
    assert prim["end3"] == "junction"
    (read,) = read_tsv(f"{prefix}.reads.tsv.gz")
    assert read["aligned_union"] == "100"
    assert read["split_events"] == "1"
    indels = {(r["kind"], r["source"], r["band"]): int(r["events"]) for r in read_tsv(f"{prefix}.indels.tsv")}
    assert indels[("ins", "cigar", "2")] == 1
    assert indels[("del", "cigar", "2")] == 1
    assert indels[("del", "split", "256")] == 1
    err = read_tsv(f"{prefix}.errpos.tsv")
    assert sum(int(r["positions"]) for r in err) == 100
    assert sum(int(r["mismatches"]) for r in err) == 1


def test_make_variants(tmp_path):
    fq = tmp_path / "r.fq"
    fq.write_text("@r1 comment\nAACGTTTG\n+\nABCDEFGH\n")
    run("make_variants.py", fq, "--rc", tmp_path / "rc.fq.gz", "--frag", f"3:{tmp_path / 'f.fq.gz'}")
    with gzip.open(tmp_path / "rc.fq.gz", "rt") as f:
        assert f.read() == "@r1\nCAAACGTT\n+\nHGFEDCBA\n"
    with gzip.open(tmp_path / "f.fq.gz", "rt") as f:
        names = [line.strip() for i, line in enumerate(f) if i % 4 == 0]
    # 8 bases in tiles of 3: 0-3, 3-6, then a 2 bp tail (>= 3/2) is kept.
    assert names == ["@r1/0-3", "@r1/3-6", "@r1/6-8"]


@pytest.mark.parametrize("d,expected", [(0, 0), (1, 1), (2, 2), (3, 2), (4, 4), (1023, 512)])
def test_dist_bin(d, expected):
    assert common.dist_bin(d) == expected
