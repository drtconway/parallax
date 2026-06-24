"""Curation web UI for parallax vs minimap2 alignment review.

Usage:
    curate.py [options]

Options:
    --config PATH    Path to curation YAML config [default: curation.yaml]
    --port=PORT      Port to listen on [default: 8000]
    --deep           Also present 'neither' reads with unchanged fingerprint

Startup:
    1. Loads curation.tsv and identifies uncurated reads (+ 'neither' if --deep)
    2. Filters single-read BAMs for all queued reads from curation BAMs
    3. Serves IGV.js UI for sequential review

    Verdict buttons write immediately to curation.tsv and advance to next read.
"""

import csv
import os
import sys
from pathlib import Path

import docopt
import pysam
import uvicorn
import yaml
from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel


# ─── Curation database ──────────────────────────────────────────────────────

COLUMNS = ['read_name', 'verdict', 'plx_fingerprint', 'mm2_fingerprint', 'notes']


def load_curation(path: Path) -> dict[str, dict]:
    rows = {}
    if not path.exists():
        raise FileNotFoundError(f"Curation database not found: {path}\nRun merge_curation.py first.")
    with open(path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rows[row['read_name']] = dict(row)
    return rows


def save_curation(path: Path, rows: dict[str, dict]) -> None:
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter='\t')
        writer.writeheader()
        for row in sorted(rows.values(), key=lambda r: r['read_name']):
            writer.writerow(row)


# ─── BAM utilities ──────────────────────────────────────────────────────────

def extract_single_read(source_bam: str, read_name: str, output_bam: str) -> None:
    """Extract all records for read_name from source_bam into output_bam,
    coordinate-sort and index so IGV.js can serve range requests."""
    with pysam.AlignmentFile(source_bam, 'rb') as src:
        header = src.header
        tmp_unsorted = output_bam + '.tmp.bam'
        with pysam.AlignmentFile(tmp_unsorted, 'wb', header=header) as out:
            for rec in src:
                if rec.is_secondary:
                    continue
                if rec.query_name == read_name:
                    out.write(rec)
    pysam.sort('-o', output_bam, tmp_unsorted)
    pysam.index(output_bam)
    os.unlink(tmp_unsorted)


def pregenerate_single_read_bams(
    queue: list[str],
    curation_bam: str,
    seeds_bam: str,
    tmpdir: Path,
) -> tuple[dict[str, Path], dict[str, Path]]:
    """Pre-extract per-read BAMs for all queued reads at startup.

    Returns (read_bams, seed_bams) — dicts mapping read_name -> Path.
    Uses a single pass over each source BAM to collect all records, then
    writes and indexes them individually.
    """
    print(f"Pre-collecting records from {curation_bam}...", file=sys.stderr)
    read_records: dict[str, list] = {}
    with pysam.AlignmentFile(curation_bam, 'rb') as bam:
        curation_header = bam.header
        for rec in bam:
            if rec.is_secondary:
                continue
            name = rec.query_name
            if name in set(queue):
                read_records.setdefault(name, []).append(rec)

    print(f"Pre-collecting records from {seeds_bam}...", file=sys.stderr)
    seed_records: dict[str, list] = {}
    with pysam.AlignmentFile(seeds_bam, 'rb') as bam:
        seeds_header = bam.header
        for rec in bam:
            if rec.is_secondary:
                continue
            name = rec.query_name
            if name in set(queue):
                seed_records.setdefault(name, []).append(rec)

    print(f"Writing {len(queue)} per-read BAMs...", file=sys.stderr)
    read_bams: dict[str, Path] = {}
    seed_bams: dict[str, Path] = {}

    for read_name in queue:
        safe = read_name.replace('/', '_').replace(' ', '_')

        # Alignment BAM
        unsorted = tmpdir / f"{safe}.plx.unsorted.bam"
        sorted_bam = tmpdir / f"{safe}.plx.bam"
        with pysam.AlignmentFile(str(unsorted), 'wb', header=curation_header) as out:
            for rec in read_records.get(read_name, []):
                out.write(rec)
        pysam.sort('-o', str(sorted_bam), str(unsorted))
        pysam.index(str(sorted_bam))
        unsorted.unlink()
        read_bams[read_name] = sorted_bam

        # Seeds BAM — only written if there are records (empty BAMs stall IGV.js)
        recs = seed_records.get(read_name, [])
        if recs:
            unsorted = tmpdir / f"{safe}.seeds.unsorted.bam"
            sorted_bam = tmpdir / f"{safe}.seeds.bam"
            with pysam.AlignmentFile(str(unsorted), 'wb', header=seeds_header) as out:
                for rec in recs:
                    out.write(rec)
            pysam.sort('-o', str(sorted_bam), str(unsorted))
            pysam.index(str(sorted_bam))
            unsorted.unlink()
            seed_bams[read_name] = sorted_bam

    return read_bams, seed_bams


# ─── Auto-detect sample ID ──────────────────────────────────────────────────

def detect_sample_id(bam_dir: Path) -> str:
    # Exclude *.curation.plx.sorted.bam — only want the full-sample BAM
    bams = [p for p in bam_dir.glob('*.plx.sorted.bam') if '.curation.' not in p.name]
    if not bams:
        raise FileNotFoundError(f"No *.plx.sorted.bam found in {bam_dir}")
    return bams[0].name.replace('.plx.sorted.bam', '')


# ─── App state ──────────────────────────────────────────────────────────────

class AppState:
    def __init__(
        self,
        results_dir: Path,
        sample_id: str,
        deep: bool,
        reference: str | None,
        tmpdir: Path,
    ):
        self.results_dir = results_dir
        self.bam_dir = results_dir / 'bam'
        self.curation_tsv = results_dir / 'curation.tsv'
        self.sample_id = sample_id
        self.reference = reference
        self.tmpdir = tmpdir
        self.extra_tracks: list[dict] = []

        self.plx_sorted_bam = str(self.bam_dir / f"{sample_id}.plx.sorted.bam")
        self.mm2_sorted_bam = str(self.bam_dir / f"{sample_id}.mm2.sorted.bam")
        self.curation_plx_bam = str(self.bam_dir / f"{sample_id}.curation.plx.sorted.bam")
        self.curation_seeds_bam = str(self.bam_dir / f"{sample_id}.curation.seeds.sorted.bam")

        self.rows = load_curation(self.curation_tsv)

        # Build queue of reads needing curation
        self.queue: list[str] = []
        for name, row in sorted(self.rows.items()):
            v = row['verdict']
            if v == 'uncurated':
                self.queue.append(name)
            elif v == 'neither' and deep:
                self.queue.append(name)
        self.index = 0

        print(f"Queue: {len(self.queue)} reads to curate", file=sys.stderr)

        # Pre-generate all single-read BAMs
        self.read_bams, self.seed_bams = pregenerate_single_read_bams(
            self.queue,
            self.curation_plx_bam,
            self.curation_seeds_bam,
            tmpdir,
        )

    def current_read(self) -> str | None:
        if self.index < len(self.queue):
            return self.queue[self.index]
        return None

    def set_verdict(self, read_name: str, verdict: str, notes: str = '') -> None:
        if read_name not in self.rows:
            raise KeyError(f"Unknown read: {read_name}")
        self.rows[read_name]['verdict'] = verdict
        if notes:
            self.rows[read_name]['notes'] = notes
        save_curation(self.curation_tsv, self.rows)

    def advance(self) -> None:
        self.index += 1

    def stats(self) -> dict:
        counts: dict[str, int] = {}
        for row in self.rows.values():
            v = row['verdict']
            counts[v] = counts.get(v, 0) + 1
        return {
            'total_queue': len(self.queue),
            'position': self.index,
            'remaining': max(0, len(self.queue) - self.index),
            'counts': counts,
        }


state: AppState | None = None
app = FastAPI()


def _path_for_url(url_path: str) -> Path | None:
    """Resolve a /bam/, /tmp_bam/, or /reference/ URL to a filesystem path."""
    assert state is not None
    if url_path.startswith('/bam/'):
        return state.bam_dir / url_path[5:]
    if url_path.startswith('/tmp_bam/'):
        return state.tmpdir / url_path[9:]
    if url_path.startswith('/reference/'):
        return Path(state.reference).parent / url_path[11:]
    return None


@app.middleware('http')
async def log_and_fix_range(request, call_next):
    range_hdr = request.headers.get('range', '')
    print(f"  >> {request.method} {request.url.path} Range={range_hdr!r}", file=sys.stderr)
    response = await call_next(request)
    print(f"  << {response.status_code} {request.url.path}", file=sys.stderr)

    # IGV.js probes file size with Range: bytes=MAX_SAFE_INT-MAX_SAFE_INT+25.
    # Starlette returns 416 without Content-Range, so IGV.js never learns the
    # real file size and stalls. Add Content-Range: bytes */size so it recovers.
    if response.status_code == 416:
        p = _path_for_url(request.url.path)
        if p and p.exists():
            size = p.stat().st_size
            response.headers['Content-Range'] = f'bytes */{size}'
            print(f"     added Content-Range: bytes */{size}", file=sys.stderr)

    return response


# ─── API endpoints ──────────────────────────────────────────────────────────

class VerdictRequest(BaseModel):
    read_name: str
    verdict: str
    notes: str = ''


@app.get('/api/current')
def get_current():
    assert state is not None
    read_name = state.current_read()
    if read_name is None:
        return JSONResponse({'done': True, 'stats': state.stats()})

    row = state.rows[read_name]
    read_bam = state.read_bams.get(read_name)
    seed_bam = state.seed_bams.get(read_name)

    return JSONResponse({
        'done': False,
        'read_name': read_name,
        'verdict': row['verdict'],
        'plx_fingerprint': row['plx_fingerprint'],
        'mm2_fingerprint': row['mm2_fingerprint'],
        'notes': row['notes'],
        'tracks': {
            'mm2_full': f'/bam/{state.sample_id}.mm2.sorted.bam',
            'plx_full': f'/bam/{state.sample_id}.plx.sorted.bam',
            'plx_read': f'/tmp_bam/{read_bam.name}' if read_bam else None,
            'plx_seeds': f'/tmp_bam/{seed_bam.name}' if seed_bam else None,
        },
        'extra_tracks': state.extra_tracks,
        'stats': state.stats(),
    })


@app.post('/api/verdict')
def post_verdict(req: VerdictRequest):
    assert state is not None
    valid = {'mm2', 'plx', 'neither', 'skip'}
    if req.verdict not in valid:
        raise HTTPException(400, f"Invalid verdict '{req.verdict}'. Must be one of {valid}")
    current = state.current_read()
    if current is None:
        raise HTTPException(400, 'No current read')
    if req.read_name != current:
        raise HTTPException(409, f"Expected verdict for '{current}', got '{req.read_name}'")
    if req.verdict != 'skip':
        state.set_verdict(req.read_name, req.verdict, req.notes)
    state.advance()
    return JSONResponse({'ok': True, 'stats': state.stats()})


@app.get('/api/stats')
def get_stats():
    assert state is not None
    return JSONResponse(state.stats())


# File serving is handled via StaticFiles mounts added at startup in main()
# (StaticFiles correctly handles HTTP Range requests that IGV.js requires)


# ─── HTML UI ────────────────────────────────────────────────────────────────

HTML = """<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8"/>
<title>Parallax Curation</title>
<script src="https://cdn.jsdelivr.net/npm/igv@3.0.2/dist/igv.min.js"></script>
<style>
  body { font-family: sans-serif; margin: 0; display: flex; flex-direction: column; height: 100vh; }
  #header { padding: 8px 16px; background: #1a1a2e; color: white; display: flex; align-items: center; gap: 16px; }
  #header h1 { margin: 0; font-size: 1.1em; }
  #progress { font-size: 0.9em; opacity: 0.8; }
  #read-info { font-size: 0.85em; opacity: 0.9; font-family: monospace; }
  #igv-container { flex: 1; overflow: hidden; }
  #controls { padding: 10px 16px; background: #f5f5f5; border-top: 1px solid #ddd;
              display: flex; align-items: center; gap: 12px; flex-wrap: wrap; }
  button { padding: 8px 20px; font-size: 1em; border: none; border-radius: 4px; cursor: pointer; }
  #btn-mm2  { background: #2196F3; color: white; }
  #btn-plx  { background: #4CAF50; color: white; }
  #btn-neither { background: #FF5722; color: white; }
  #btn-skip { background: #9E9E9E; color: white; }
  button:hover { opacity: 0.85; }
  #notes { flex: 1; padding: 7px; border: 1px solid #ccc; border-radius: 4px; font-size: 0.95em; }
  #done-banner { display: none; padding: 40px; text-align: center; font-size: 1.4em; color: #333; }
  #seg-nav { display: none; align-items: center; gap: 6px; font-size: 0.9em; }
  #seg-nav button { padding: 4px 10px; font-size: 0.9em; background: #455a64; color: white; }
  #seg-label { min-width: 120px; text-align: center; font-family: monospace; font-size: 0.8em; }
  #loading-overlay { display: none; position: absolute; inset: 0; background: rgba(255,255,255,0.75);
                     align-items: center; justify-content: center; font-size: 1.3em; color: #555;
                     z-index: 100; pointer-events: none; }
  #igv-wrapper { position: relative; flex: 1; overflow: hidden; }
</style>
</head>
<body>
<div id="header">
  <h1>Parallax Curation</h1>
  <span id="progress"></span>
  <span id="read-info"></span>
</div>
<div id="igv-wrapper">
  <div id="igv-container"></div>
  <div id="loading-overlay">Loading&hellip;</div>
</div>
<div id="controls">
  <button id="btn-mm2"     onclick="submitVerdict('mm2')">minimap2 correct</button>
  <button id="btn-plx"     onclick="submitVerdict('plx')">parallax correct</button>
  <button id="btn-neither" onclick="submitVerdict('neither')">neither</button>
  <button id="btn-skip"    onclick="submitVerdict('skip')">skip</button>
  <div id="seg-nav">
    <button onclick="goToSegment(0)" title="Return to first segment">&#8962;</button>
    <button id="seg-prev" onclick="stepSegment(-1)">&#8592;</button>
    <span id="seg-label"></span>
    <button id="seg-next" onclick="stepSegment(+1)">&#8594;</button>
  </div>
  <input  id="notes" type="text" placeholder="Notes (optional)"/>
</div>
<div id="done-banner" id="done-banner">
  All reads curated. Run build_truth_bam.py to generate the truth BAM.
</div>

<script>
const REFERENCE = %REFERENCE_JSON%;
let browser = null;
let currentRead = null;
let segments = [];   // [{chrom, start, end, label}] — one per distinct locus
let segIndex = 0;

// Parse all distinct loci from a fingerprint string.
// Each segment becomes a padded locus (10% of segment width on each side);
// segments on the same chrom that overlap after padding are merged into one stop.
function parseLoci(fp) {
  if (!fp) return [];
  const loci = [];
  for (const seg of fp.split(';')) {
    const m = seg.match(/^(\\S+):(\\d+)-(\\d+):/);
    if (!m) continue;
    const start = parseInt(m[2]), end = parseInt(m[3]);
    const pad = Math.round((end - start) * 0.1);
    loci.push({ chrom: m[1], start: start - pad, end: end + pad });
  }
  return loci;
}

// Merge loci on the same chrom only when they overlap or are within MAX_MERGE_GAP bp.
// Segments further apart are kept as separate navigation stops.
const MAX_MERGE_GAP = 50000;
function mergeLoci(loci) {
  const merged = [];
  for (const l of loci) {
    const last = merged[merged.length - 1];
    if (last && last.chrom === l.chrom && l.start <= last.end + MAX_MERGE_GAP) {
      last.start = Math.min(last.start, l.start);
      last.end   = Math.max(last.end,   l.end);
    } else {
      merged.push({...l});
    }
  }
  return merged;
}

function locusStr(seg) {
  return `${seg.chrom}:${Math.max(0, seg.start)}-${seg.end}`;
}

async function goToSegment(idx) {
  segIndex = Math.max(0, Math.min(idx, segments.length - 1));
  const seg = segments[segIndex];
  document.getElementById('seg-label').textContent =
    `seg ${segIndex + 1}/${segments.length}  ${seg.chrom}`;
  await browser.search(locusStr(seg));
}

async function stepSegment(delta) {
  await goToSegment(segIndex + delta);
}

async function loadCurrent() {
  const resp = await fetch('/api/current');
  const data = await resp.json();

  if (data.done) {
    document.getElementById('controls').style.display = 'none';
    document.getElementById('igv-container').style.display = 'none';
    document.getElementById('done-banner').style.display = 'block';
    document.getElementById('progress').textContent = 'Done!';
    return;
  }

  currentRead = data.read_name;
  document.getElementById('notes').value = '';
  document.getElementById('progress').textContent =
    `${data.stats.position + 1} / ${data.stats.total_queue}`;
  document.getElementById('read-info').textContent = currentRead;

  // Build segment list from both fingerprints, merge overlaps
  const raw = [
    ...parseLoci(data.plx_fingerprint),
    ...parseLoci(data.mm2_fingerprint),
  ];
  // Deduplicate by chrom:start-end string before merging
  const seen = new Set();
  const deduped = raw.filter(l => {
    const k = locusStr(l);
    if (seen.has(k)) return false;
    seen.add(k);
    return true;
  });
  deduped.sort((a, b) => a.chrom < b.chrom ? -1 : a.chrom > b.chrom ? 1 : a.start - b.start);
  segments = mergeLoci(deduped);
  segIndex = 0;

  // Always show the nav (for the home button); show arrows only with multiple segments
  const nav = document.getElementById('seg-nav');
  nav.style.display = 'flex';
  const multi = segments.length > 1;
  document.getElementById('seg-prev').style.display = multi ? '' : 'none';
  document.getElementById('seg-next').style.display = multi ? '' : 'none';
  document.getElementById('seg-label').style.display = multi ? '' : 'none';

  const tracks = [
    { name: 'minimap2',             url: data.tracks.mm2_full,  indexURL: data.tracks.mm2_full  + '.bai', format: 'bam', color: '#2196F3', height: 200, visibilityWindow: -1 },
    { name: 'parallax',             url: data.tracks.plx_full,  indexURL: data.tracks.plx_full  + '.bai', format: 'bam', color: '#4CAF50', height: 200, visibilityWindow: -1 },
    { name: 'parallax (this read)', url: data.tracks.plx_read,  indexURL: data.tracks.plx_read  + '.bai', format: 'bam', color: '#FF9800', height: 100, visibilityWindow: -1 },
  ];
  if (data.tracks.plx_seeds) {
    tracks.push({ name: 'seeds (this read)', url: data.tracks.plx_seeds, indexURL: data.tracks.plx_seeds + '.bai', format: 'bam', color: '#9C27B0', height: 150, visibilityWindow: -1 });
  }
  for (const t of data.extra_tracks || []) {
    tracks.push({ visibilityWindow: -1, ...t, name: t.name || t.url.split('/').pop() });
  }

  if (!browser) {
    const options = { reference: REFERENCE, locus: locusStr(segments[0]), tracks: tracks, showSVGButton: false, search: false, blat: false };
    browser = await igv.createBrowser(document.getElementById('igv-container'), options);
  } else {
    const existing = browser.trackViews.map(tv => tv.track).filter(t => t.type !== 'sequence' && t.id !== 'ruler');
    for (const t of existing) { await browser.removeTrack(t); }
    for (const t of tracks)   { await browser.loadTrack(t); }
  }
  await goToSegment(0);
  setLoading(false);
}

function setLoading(on) {
  const overlay = document.getElementById('loading-overlay');
  overlay.style.display = on ? 'flex' : 'none';
  for (const id of ['btn-mm2', 'btn-plx', 'btn-neither', 'btn-skip'])
    document.getElementById(id).disabled = on;
}

async function submitVerdict(verdict) {
  if (!currentRead) return;
  setLoading(true);
  const notes = document.getElementById('notes').value.trim();
  await fetch('/api/verdict', {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify({ read_name: currentRead, verdict, notes }),
  });
  loadCurrent();
}

// Keyboard shortcuts
document.addEventListener('keydown', e => {
  if (e.target.id === 'notes') return;
  if (e.key === '1') submitVerdict('mm2');
  if (e.key === '2') submitVerdict('plx');
  if (e.key === '3') submitVerdict('neither');
  if (e.key === '4') submitVerdict('skip');
  if (e.key === 'ArrowLeft')  stepSegment(-1);
  if (e.key === 'ArrowRight') stepSegment(+1);
});

loadCurrent();
</script>
</body>
</html>
"""


@app.get('/', response_class=HTMLResponse)
def index():
    assert state is not None
    ref_name = Path(state.reference).name
    ref_json = f'{{"fastaURL": "/reference/{ref_name}", "indexURL": "/reference/{ref_name}.fai"}}'
    return HTML.replace('%REFERENCE_JSON%', ref_json)


# ─── Entry point ────────────────────────────────────────────────────────────

def main(args):
    global state

    config_path = Path(args['--config'] or 'curation.yaml')
    if not config_path.exists():
        print(f"Error: config file not found: {config_path}", file=sys.stderr)
        sys.exit(1)
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    results_dir = Path(cfg['outdir'])
    if not results_dir.is_dir():
        print(f"Error: results directory not found: {results_dir}", file=sys.stderr)
        sys.exit(1)

    reference = cfg['reference']
    bam_dir = results_dir / 'bam'
    sample_id = detect_sample_id(bam_dir)
    port = int(args['--port'] or 8000)
    deep = bool(args['--deep'])

    if not Path(reference + '.fai').exists():
        print(f"Error: reference index not found — run: samtools faidx {reference}", file=sys.stderr)
        sys.exit(1)

    print(f"Sample: {sample_id}", file=sys.stderr)

    tmpdir = results_dir.resolve() / 'tmp'
    tmpdir.mkdir(parents=True, exist_ok=True)
    print(f"Temp BAMs: {tmpdir}", file=sys.stderr)

    state = AppState(
        results_dir=results_dir,
        sample_id=sample_id,
        deep=deep,
        reference=reference,
        tmpdir=tmpdir,
    )

    # Create a reference dir containing symlinks to the FASTA and .fai so
    # StaticFiles can serve them with proper Range request support
    ref_path = Path(reference)
    refdir = tmpdir / 'reference'
    refdir.mkdir(exist_ok=True)
    for target, src in [
        (refdir / ref_path.name,             ref_path),
        (refdir / (ref_path.name + '.fai'),  Path(reference + '.fai')),
    ]:
        if not target.exists():
            target.symlink_to(src)

    # Built-in gene tracks for known genomes
    GENE_TRACKS = {
        'hg38': {
            'name': 'Genes',
            'type': 'annotation',
            'format': 'refgene',
            'url': 'https://s3.amazonaws.com/igv.org.genomes/hg38/ncbiRefSeq.sorted.txt.gz',
            'indexURL': 'https://s3.amazonaws.com/igv.org.genomes/hg38/ncbiRefSeq.sorted.txt.gz.tbi',
            'visibilityWindow': 300000000,
            'displayMode': 'COLLAPSED',
            'height': 60,
            'order': -1,
        },
    }

    # Symlink extra track files into tmpdir/extra/ so StaticFiles can serve them
    extradir = tmpdir / 'extra'
    extradir.mkdir(exist_ok=True)
    extra_tracks = []
    for t in cfg.get('tracks') or []:
        src = Path(t['path']).resolve()
        if not src.exists():
            print(f"Warning: extra track file not found, skipping: {src}", file=sys.stderr)
            continue
        link = extradir / src.name
        if not link.exists():
            link.symlink_to(src)
        # Symlink companion index files (.bai, .tbi, .csi) if present
        for ext in ('.bai', '.tbi', '.csi'):
            idx_src = Path(str(src) + ext)
            if idx_src.exists():
                idx_link = extradir / (src.name + ext)
                if not idx_link.exists():
                    idx_link.symlink_to(idx_src)
        entry = {'url': f'/extra/{src.name}'}
        if 'name' in t:
            entry['name'] = t['name']
        if 'color' in t:
            entry['color'] = t['color']
        extra_tracks.append(entry)

    genome = cfg.get('genome')
    if genome:
        gene_track = GENE_TRACKS.get(genome.lower())
        if gene_track:
            extra_tracks.insert(0, gene_track)
        else:
            print(f"Warning: no built-in gene track for genome '{genome}'", file=sys.stderr)

    state.extra_tracks = extra_tracks

    # Mount static file directories — StaticFiles handles Range requests
    # correctly, which is required by IGV.js for BAM fetching
    app.mount('/bam',       StaticFiles(directory=str(state.bam_dir)), name='bam')
    app.mount('/tmp_bam',   StaticFiles(directory=str(tmpdir)),        name='tmp_bam')
    app.mount('/reference', StaticFiles(directory=str(refdir), follow_symlink=True), name='reference')
    app.mount('/extra',     StaticFiles(directory=str(extradir), follow_symlink=True), name='extra')

    print(f"\nOpen http://localhost:{port}/ in your browser", file=sys.stderr)
    print("Keyboard shortcuts: 1=mm2  2=plx  3=neither  4=skip", file=sys.stderr)
    uvicorn.run(app, host='0.0.0.0', port=port, log_level='info')


if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    main(args)
