from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from .aligner import align, AlignmentLocus
from .results import new_result_dir

app = FastAPI()


class AlignRequest(BaseModel):
    fastq_path: str
    context_bam_paths: list[str] = []


class AlignResponse(BaseModel):
    result_id: str
    bam_url: str
    bai_url: str
    context_tracks: list[dict]
    records: list[AlignmentLocus]


@app.post("/api/align", response_model=AlignResponse)
async def align_read(req: AlignRequest):
    fastq_path = Path(req.fastq_path)
    if not fastq_path.exists():
        raise HTTPException(status_code=400, detail=f"FASTQ not found: {fastq_path}")

    for p in req.context_bam_paths:
        if not Path(p).exists():
            raise HTTPException(status_code=400, detail=f"BAM not found: {p}")
        if not Path(p + ".bai").exists():
            raise HTTPException(status_code=400, detail=f"BAI not found: {p}.bai")

    result_id, result_dir = new_result_dir()

    try:
        _, loci = await align(fastq_path, result_dir)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

    context_tracks = [
        {
            "name": Path(p).name,
            "url": f"/api/context-bam?path={p}",
            "indexURL": f"/api/context-bam?path={p}.bai",
            "format": "bam",
        }
        for p in req.context_bam_paths
    ]

    return AlignResponse(
        result_id=result_id,
        bam_url=f"/results/{result_id}/alignment.bam",
        bai_url=f"/results/{result_id}/alignment.bam.bai",
        context_tracks=context_tracks,
        records=loci,
    )


@app.get("/api/context-bam")
async def context_bam(path: str):
    p = Path(path)
    if not p.exists():
        raise HTTPException(status_code=404, detail="file not found")
    return FileResponse(p)


@app.get("/results/{result_id}/{filename}")
async def result_file(result_id: str, filename: str):
    from .config import RESULTS_DIR
    path = RESULTS_DIR / result_id / filename
    if not path.exists():
        raise HTTPException(status_code=404, detail="not found")
    return FileResponse(path)
