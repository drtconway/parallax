import logging
from pathlib import Path

from fastapi import FastAPI, HTTPException

log = logging.getLogger(__name__)
from fastapi.responses import FileResponse
from pydantic import BaseModel

from .aligner import align, accept, AlignmentLocus
from .config import RESULTS_DIR
from .results import split_fastq, read_digest, ensure_read, result_dir_for, alignment_status

app = FastAPI()


class AlignRequest(BaseModel):
    fastq_path: str
    context_bam_paths: list[str] = []


class ReadResult(BaseModel):
    digest: str
    status: str
    bam_url: str
    bai_url: str
    expected_bam_url: str | None
    expected_bai_url: str | None
    records: list[AlignmentLocus]


class AlignResponse(BaseModel):
    results: list[ReadResult]
    context_tracks: list[dict]


@app.post("/api/align", response_model=AlignResponse)
async def align_reads(req: AlignRequest):
    fastq_path = Path(req.fastq_path)
    if not fastq_path.exists():
        raise HTTPException(status_code=400, detail=f"FASTQ not found: {fastq_path}")

    for p in req.context_bam_paths:
        if not Path(p).exists():
            raise HTTPException(status_code=400, detail=f"BAM not found: {p}")
        if not Path(p + ".bai").exists():
            raise HTTPException(status_code=400, detail=f"BAI not found: {p}.bai")

    with open(fastq_path, "rb") as f:
        raw = f.read()

    records_raw = split_fastq(raw)
    if not records_raw:
        raise HTTPException(status_code=400, detail="No FASTQ records found in file")

    log.info("Processing %d FASTQ record(s) from %s", len(records_raw), fastq_path)

    results: list[ReadResult] = []
    for record in records_raw:
        log.info("Aligning record: %s", record[:50])
        digest = read_digest(record)
        result_dir = ensure_read(digest, record)
        read_path = result_dir / "read.fastq.gz"

        try:
            _, loci = await align(read_path, result_dir)
        except Exception as e:
            raise HTTPException(status_code=500, detail=str(e))

        status = alignment_status(digest)

        has_expected = (result_dir / "expected.bam").exists()
        results.append(ReadResult(
            digest=digest,
            status=status,
            bam_url=f"/results/{digest}/alignment.bam",
            bai_url=f"/results/{digest}/alignment.bam.bai",
            expected_bam_url=f"/results/{digest}/expected.bam" if has_expected else None,
            expected_bai_url=f"/results/{digest}/expected.bam.bai" if has_expected else None,
            records=loci,
        ))

    context_tracks = [
        {
            "name": Path(p).name,
            "url": f"/api/context-bam?path={p}",
            "indexURL": f"/api/context-bam?path={p}.bai",
            "format": "bam",
        }
        for p in req.context_bam_paths
    ]

    return AlignResponse(results=results, context_tracks=context_tracks)


@app.post("/api/accept/{digest}")
async def accept_alignment(digest: str):
    result_dir = result_dir_for(digest)
    if not (result_dir / "alignment.bam").exists():
        raise HTTPException(status_code=404, detail="No alignment found for this digest")
    try:
        accept(result_dir)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    return {"status": "accepted"}


@app.get("/api/context-bam")
async def context_bam(path: str):
    p = Path(path)
    if not p.exists():
        raise HTTPException(status_code=404, detail="file not found")
    return FileResponse(p)


@app.get("/results/{digest}/{filename}")
async def result_file(digest: str, filename: str):
    path = RESULTS_DIR / digest / filename
    if not path.exists():
        raise HTTPException(status_code=404, detail="not found")
    return FileResponse(path)
