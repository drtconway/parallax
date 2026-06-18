import gzip
import hashlib
import io
from pathlib import Path

from .config import RESULTS_DIR


def read_digest(fastq_record: bytes) -> str:
    return hashlib.md5(fastq_record).hexdigest()


def result_dir_for(digest: str) -> Path:
    return RESULTS_DIR / digest


def _decompress_if_needed(data: bytes) -> bytes:
    if data[:2] == b"\x1f\x8b":
        with gzip.open(io.BytesIO(data)) as f:
            return f.read()
    return data


def split_fastq(data: bytes) -> list[bytes]:
    """Split a FASTQ file into individual records (4 lines each)."""
    data = _decompress_if_needed(data)
    lines = [l for l in data.splitlines() if l.strip()]
    records = []
    for i in range(0, len(lines) - 3, 4):
        group = lines[i:i + 4]
        if len(group) == 4 and group[0].startswith(b"@"):
            record = b"\n".join(group) + b"\n"
            records.append(record)
    return records


def ensure_read(digest: str, record: bytes) -> Path:
    """Write read.fastq.gz if it doesn't already exist. Returns the result dir."""
    result_dir = result_dir_for(digest)
    result_dir.mkdir(parents=True, exist_ok=True)
    read_path = result_dir / "read.fastq.gz"
    if not read_path.exists():
        with gzip.open(read_path, "wb") as f:
            f.write(record)
    return result_dir


def alignment_status(digest: str) -> str:
    """
    Return the status of a result directory:
      'pending'  — no expected.bam
      'passing'  — expected.bam exists and matches alignment.bam
      'failing'  — expected.bam exists but differs from alignment.bam
      'missing'  — no alignment.bam yet
    """
    result_dir = result_dir_for(digest)
    bam_path = result_dir / "alignment.bam"
    expected_path = result_dir / "expected.bam"

    if not bam_path.exists():
        return "missing"
    if not expected_path.exists():
        return "pending"

    from .aligner import primary_alignment_keys
    current = set(primary_alignment_keys(bam_path))
    expected = set(primary_alignment_keys(expected_path))
    return "passing" if current == expected else "failing"
