import subprocess
from pathlib import Path
from typing import TypedDict

import httpx
import pysam

from .config import PARALLAX_URL


class AlignmentLocus(TypedDict):
    name: str
    chrom: str
    start: int
    end: int


async def align(fastq_path: Path, result_dir: Path) -> tuple[Path, list[AlignmentLocus]]:
    with open(fastq_path, "rb") as f:
        data = f.read()

    async with httpx.AsyncClient() as client:
        response = await client.post(
            f"{PARALLAX_URL}/align",
            content=data,
            headers={"Content-Type": "application/octet-stream"},
            timeout=120.0,
        )
        response.raise_for_status()

    bam_path = result_dir / "alignment.bam"
    bam_path.write_bytes(response.content)
    _index_bam(bam_path)

    loci = _primary_loci(bam_path)
    return bam_path, loci


def accept(result_dir: Path) -> None:
    """Copy alignment.bam to expected.bam and index it."""
    bam_path = result_dir / "alignment.bam"
    expected_path = result_dir / "expected.bam"
    expected_path.write_bytes(bam_path.read_bytes())
    _index_bam(expected_path)


def _index_bam(bam_path: Path) -> None:
    result = subprocess.run(
        ["samtools", "index", str(bam_path)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"samtools index failed: {result.stderr.strip()}")


def primary_alignment_keys(bam_path: Path) -> list[tuple]:
    """Return (chrom, start, end, flag, cigar) tuples for primary alignments."""
    keys = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch():
            if read.flag & 0x904:
                continue
            keys.append((
                read.reference_name or "",
                read.reference_start,
                read.reference_end,
                read.flag,
                read.cigarstring or "",
            ))
    return keys


def _primary_loci(bam_path: Path) -> list[AlignmentLocus]:
    loci: list[AlignmentLocus] = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch():
            if read.flag & 0x904:
                continue
            start = (read.reference_start or 0) + 1  # 1-based
            end = read.reference_end or start
            loci.append({
                "name": read.query_name or "",
                "chrom": read.reference_name or "",
                "start": start,
                "end": end,
            })
    return loci
