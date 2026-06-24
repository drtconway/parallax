# Alignment Curation

This document describes the process for curating parallax alignments against
minimap2, building a truth BAM, and keeping that truth BAM up to date as
parallax changes.

## Overview

The curation pipeline has four stages:

1. **Compare** — run both aligners and identify reads where they disagree
2. **Merge** — update the curation database with new disagreements
3. **Curate** — visually inspect disagreeing reads and record verdicts
4. **Build truth** — derive a truth BAM from the curated verdicts

Stages 1–2 are re-run each time parallax changes. Stage 3 only presents reads
that are new or whose parallax alignment has changed. Stage 4 is run on demand.

---

## Configuration

All stages are driven by a single YAML file. Copy and edit the example:

```bash
cp curation.yaml my_curation.yaml   # or just edit curation.yaml in place
```

Key fields:

| Field | Description |
|---|---|
| `reference` | Reference FASTA path (must have a `.fai` index: `samtools faidx`) |
| `fastq` | Input FASTQ (glob ok); or use `samplesheet` instead |
| `samplesheet` | CSV with columns `sample,file`; alternative to `fastq` |
| `index` | Pre-built parallax index directory |
| `parallax` | Path to the parallax binary |
| `parallax_config` | Optional TOML config for parallax |
| `outdir` | Directory for all outputs (BAMs, `curation.tsv`, `truth.bam`) |
| `threads` | Threads per process |

---

## Prerequisites

Install Python dependencies (once, into the project venv):

```bash
pip install fastapi uvicorn pydantic docopt pysam pyyaml shapely
```

Build the parallax release binary if not already done:

```bash
cargo build --release
```

---

## Stage 1: Run the comparison pipeline

```bash
nextflow run compare.nf -params-file curation.yaml
```

To override individual parameters on the command line:

```bash
nextflow run compare.nf -params-file curation.yaml --threads 16
```

The pipeline produces under `outdir/bam/`:

| File | Description |
|---|---|
| `<id>.plx.nsorted.bam` | Parallax alignments, name-sorted |
| `<id>.mm2.nsorted.bam` | Minimap2 alignments, name-sorted |
| `<id>.plx.sorted.bam[.bai]` | Parallax alignments, coordinate-sorted + indexed |
| `<id>.mm2.sorted.bam[.bai]` | Minimap2 alignments, coordinate-sorted + indexed |
| `<id>.curation.plx.sorted.bam[.bai]` | Parallax alignments for disagreeing reads only |
| `<id>.curation.seeds.sorted.bam[.bai]` | Parallax seeds for disagreeing reads |

And under `outdir/`:

| File | Description |
|---|---|
| `<id>.compare.txt` | Per-read comparison results (TSV) |

---

## Stage 2: Update the curation database

```bash
python3 scripts/merge_curation.py
# or with a non-default config:
python3 scripts/merge_curation.py --config my_curation.yaml
```

`outdir/curation.tsv` is created on first run. On subsequent runs, the merge script:

- Adds new disagreements as `uncurated`
- Resets any previously curated read to `uncurated` if parallax has changed
  its alignment for that read (detected via alignment fingerprint)
- Keeps `neither` verdicts where parallax is unchanged (these are skipped in
  normal curation mode; use `--deep` to revisit them)
- Marks reads that now agree as `agree` (excluded from the curation queue)

### Curation database format

`outdir/curation.tsv` is a tab-separated file with columns:

| Column | Description |
|---|---|
| `read_name` | Read identifier |
| `verdict` | One of `uncurated`, `mm2`, `plx`, `neither`, `agree` |
| `plx_fingerprint` | Parallax alignment fingerprint at time of curation |
| `mm2_fingerprint` | Minimap2 alignment fingerprint |
| `notes` | Free-text notes |

**Verdicts:**

| Verdict | Meaning |
|---|---|
| `uncurated` | Not yet inspected |
| `mm2` | Minimap2 alignment is correct |
| `plx` | Parallax alignment is correct |
| `neither` | Neither alignment is acceptable |
| `agree` | Both aligners agree (not shown in curation UI) |

**Fingerprint format:** one segment per primary/supplementary alignment,
semicolon-separated, sorted. Each segment is
`chrom:ref_start-ref_end:strand:query_start-query_end`.
Query coordinates are derived from the CIGAR string. If parallax produces a
different fingerprint on re-run, the read is reset to `uncurated`.

---

## Stage 3: Curate

Start the curation web UI:

```bash
python3 scripts/curate.py
# or with a non-default config:
python3 scripts/curate.py --config my_curation.yaml
```

Options:

| Option | Description |
|---|---|
| `--config PATH` | Config file [default: curation.yaml] |
| `--port PORT` | Port to listen on [default: 8000] |
| `--deep` | Also present `neither` reads with unchanged parallax alignment |

Then open `http://localhost:8000/` in a browser.

### IGV tracks

For each read under review, four tracks are shown:

1. **minimap2** — full coordinate-sorted minimap2 BAM (context)
2. **parallax** — full coordinate-sorted parallax BAM (context)
3. **parallax (this read)** — single-read parallax BAM (highlighted in orange)
4. **seeds (this read)** — parallax seeds for this read (purple)

The view is centred on the locus of the read under review with 10% padding.

### Verdict buttons and keyboard shortcuts

| Button | Key | Meaning |
|---|---|---|
| minimap2 correct | `1` | Minimap2 alignment is correct |
| parallax correct | `2` | Parallax alignment is correct |
| neither | `3` | Neither alignment is acceptable |
| skip | `4` | Leave as uncurated, move to next read |

Verdicts are written to `curation.tsv` immediately. The UI survives
interruption — restart it and it resumes from the first uncurated read.

---

## Stage 4: Build the truth BAM

```bash
python3 scripts/build_truth_bam.py
# or with a non-default config:
python3 scripts/build_truth_bam.py --config my_curation.yaml
```

The output is written to `outdir/truth.bam` — a name-sorted BAM containing:

- Reads with verdict `plx` — alignments from the parallax BAM
- Reads with verdict `mm2` — alignments from the minimap2 BAM
- Reads with verdict `agree` — alignments from the parallax BAM

Reads with verdict `neither` or `uncurated` are omitted. The script prints a
warning for any omitted uncurated reads.

---

## Iterating after a parallax change

1. Rebuild the parallax binary (`cargo build --release`)
2. Re-run the comparison pipeline (`nextflow run compare.nf -params-file curation.yaml`)
3. Re-run the merge script (`python3 scripts/merge_curation.py`)
   - Reads where parallax changed are automatically reset to `uncurated`
   - Previously curated reads that are unchanged are kept
4. Run the curation UI — only changed/new reads are presented
5. Rebuild the truth BAM (`python3 scripts/build_truth_bam.py`)

---

## File layout

```
curation.yaml                         shared configuration

outdir/                               (value of outdir in curation.yaml)
  bam/
    <id>.plx.nsorted.bam              parallax, name-sorted
    <id>.mm2.nsorted.bam              minimap2, name-sorted
    <id>.plx.sorted.bam[.bai]         parallax, coordinate-sorted
    <id>.mm2.sorted.bam[.bai]         minimap2, coordinate-sorted
    <id>.curation.plx.sorted.bam[.bai]    parallax subset for curation
    <id>.curation.seeds.sorted.bam[.bai]  parallax seeds for curation
  <id>.compare.txt                    per-read comparison results
  curation.tsv                        curation database (edit carefully)
  truth.bam                           derived truth BAM (regenerate on demand)
  tmp/                                transient per-read BAMs for the UI
```
