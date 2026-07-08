"""Merge new comparison output into the curation database.

Usage:
    merge_curation.py [--config PATH]

Options:
    --config PATH   Path to curation YAML config [default: curation.yaml]

Reads source paths from the config file (reference, outdir, fastq/samplesheet).
Expects the compare pipeline to have already run, with its outputs under outdir/.

The curation TSV has columns:
    read_name  verdict  plx_fingerprint  mm2_fingerprint  notes

Verdicts:
    uncurated   not yet inspected
    mm2         minimap2 alignment accepted as correct
    plx         parallax alignment accepted as correct
    neither     neither alignment acceptable (fingerprinted to detect change)
    agree       alignments now agree -- kept for audit, not shown in UI

Merge rules:
    - Reads in compare_tsv with state != 'agree' and jaccard < 1.0 are
      disagreements and must appear in the curation database.
    - New disagreements are added as 'uncurated'.
    - Existing 'mm2' or 'plx' verdicts where plx_fingerprint has changed
      are reset to 'uncurated' (parallax changed its answer).
    - Existing 'neither' verdicts where plx_fingerprint has changed are
      reset to 'uncurated'.
    - Existing 'neither' verdicts with unchanged plx_fingerprint are kept
      (skipped by default in curate.py unless --deep).
    - Reads that now agree (jaccard == 1.0) have their verdict set to
      'agree' so they are excluded from future curation queues.
"""

import csv
import hashlib
import sys
import tempfile
import zipfile
from pathlib import Path

import docopt
import pysam
import yaml

from fingerprint import fingerprints_from_bam

COLUMNS = ['read_name', 'verdict', 'plx_fingerprint', 'mm2_fingerprint', 'notes']
VERDICTS = {'uncurated', 'mm2', 'plx', 'neither', 'agree', 'clip'}


def load_config(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def resolve_sample_id(bam_dir: Path) -> str:
    bams = list(bam_dir.glob('*.plx.nsorted.bam'))
    if not bams:
        raise FileNotFoundError(f"No *.plx.nsorted.bam found in {bam_dir}")
    return bams[0].name.replace('.plx.nsorted.bam', '')


def load_curation(path: Path) -> dict[str, dict]:
    rows = {}
    if not path.exists():
        return rows
    with open(path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rows[row['read_name']] = row
    return rows


def save_curation(path: Path, rows: dict[str, dict]) -> None:
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter='\t')
        writer.writeheader()
        for row in sorted(rows.values(), key=lambda r: r['read_name']):
            writer.writerow(row)


def load_disagreements(compare_tsv: Path) -> tuple[set[str], set[str], set[str]]:
    """Return (disagreeing_reads, agreeing_reads, clip_reads) from compare output.

    State values from compare_alignments_2.py:
      'both'  — both aligners produced output; jaccard measures overlap
      'lhs'   — only parallax aligned the read (minimap2 missed it)
      'rhs'   — only minimap2 aligned the read (parallax missed it)
      'clip'  — one or both alignments are heavily clipped

    'lhs'-only reads are not disagreements that need curation (parallax found
    something minimap2 didn't — that may well be correct).  'rhs'-only reads
    are the most important case: parallax missed the read entirely and should
    be reviewed.  Both are treated as disagreements here.
    """
    disagreeing = set()
    agreeing = set()
    clipping = set()
    with open(compare_tsv, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            read = row['read_id']
            state = row['state']
            jaccard = float(row['jaccard'])
            if state == 'clip':
                clipping.add(read)
            elif state == 'agree' or (state == 'both' and jaccard >= 1.0):
                agreeing.add(read)
            else:
                # 'both' with jaccard < 1.0, 'lhs'-only, or 'rhs'-only
                disagreeing.add(read)
    return disagreeing, agreeing, clipping


def zip_entry_path(safe: str, suffix: str) -> str:
    """Return the zip-internal path for a file, e.g. 'a3/SRR29147690.1.plx.bam'."""
    bucket = hashlib.md5(safe.encode()).hexdigest()[:2]
    return f"{bucket}/{safe}{suffix}"


def _write_bam_to_zip(
    zf: zipfile.ZipFile,
    read_name: str,
    records: list,
    header: pysam.AlignmentHeader,
    suffix: str,
    tmpdir: Path,
) -> None:
    """Sort and index a BAM into a temp dir, then add both files to the zip."""
    safe = read_name.replace('/', '_').replace(' ', '_')
    unsorted = tmpdir / f"{safe}{suffix}.unsorted.bam"
    sorted_bam = tmpdir / f"{safe}{suffix}.bam"
    try:
        with pysam.AlignmentFile(str(unsorted), 'wb', header=header) as out:
            for rec in records:
                out.write(rec)
        pysam.sort('-o', str(sorted_bam), str(unsorted))
        pysam.index(str(sorted_bam))
        bai = Path(str(sorted_bam) + '.bai')
        zf.write(sorted_bam, zip_entry_path(safe, f"{suffix}.bam"))
        zf.write(bai, zip_entry_path(safe, f"{suffix}.bam.bai"))
    finally:
        for p in (unsorted, sorted_bam, Path(str(sorted_bam) + '.bai')):
            if p.exists():
                p.unlink()


def _collect_records(bam_path: str, needed: set[str]) -> tuple[dict[str, list], pysam.AlignmentHeader | None]:
    """Scan a name-sorted BAM and collect all non-secondary records for reads in `needed`."""
    if not Path(bam_path).exists():
        return {}, None
    records: dict[str, list] = {}
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        header = bam.header
        for rec in bam:
            if rec.is_secondary:
                continue
            if rec.query_name in needed:
                records.setdefault(rec.query_name, []).append(rec)
    return records, header


def generate_per_read_bams(
    queue: list[str],
    curation_bam: str,
    full_plx_bam: str,
    seeds_bam: str,
    zip_path: Path,
) -> None:
    """Extract per-read BAMs for all queued reads into a zip archive.

    Reads are taken from `curation_bam` (the subset re-aligned for curation)
    when available, falling back to `full_plx_bam` (the full-sample BAM) for
    any reads not present in the curation BAM.  Always rebuilds the zip from
    scratch.
    """
    if not queue:
        print("No reads to process.", file=sys.stderr)
        return

    needed_set = set(queue)

    print(f"Collecting records for {len(queue)} reads from {curation_bam}...", file=sys.stderr)
    read_records, curation_header = _collect_records(curation_bam, needed_set)

    # Fall back to the full BAM for reads not found in the curation BAM.
    missing_from_curation = needed_set - set(read_records)
    if missing_from_curation and Path(full_plx_bam).exists():
        print(f"  {len(missing_from_curation)} reads not in curation BAM, "
              f"falling back to {full_plx_bam}...", file=sys.stderr)
        for r in sorted(missing_from_curation):
            print(f"    missing from curation BAM: {r}", file=sys.stderr)
        fallback_records, fallback_header = _collect_records(full_plx_bam, missing_from_curation)
        read_records.update(fallback_records)
        if curation_header is None and fallback_header is not None:
            curation_header = fallback_header

    still_missing = needed_set - set(read_records)
    if still_missing:
        print(f"  Warning: {len(still_missing)} reads not found in any BAM, "
              f"they will be absent from reads.zip", file=sys.stderr)

    print(f"Collecting seed records from {seeds_bam}...", file=sys.stderr)
    seed_records, seeds_header = _collect_records(seeds_bam, needed_set)

    print(f"Writing {len(queue)} per-read BAMs to {zip_path}...", file=sys.stderr)
    tmp_path = zip_path.with_suffix('.zip.tmp')
    with tempfile.TemporaryDirectory() as tmpdir_str:
        tmpdir = Path(tmpdir_str)
        with zipfile.ZipFile(tmp_path, 'w', compression=zipfile.ZIP_STORED) as zf:
            for i, read_name in enumerate(queue):
                if (i + 1) % 500 == 0:
                    print(f"  {i + 1}/{len(queue)}...", file=sys.stderr)
                if read_records.get(read_name):
                    _write_bam_to_zip(zf, read_name, read_records[read_name],
                                      curation_header, '.plx', tmpdir)
                if seeds_header and seed_records.get(read_name):
                    _write_bam_to_zip(zf, read_name, seed_records[read_name],
                                      seeds_header, '.seeds', tmpdir)
    tmp_path.replace(zip_path)
    print(f"Done. Wrote {zip_path}", file=sys.stderr)


def main(args):
    config_path = Path(args['--config'] or 'curation.yaml')
    if not config_path.exists():
        print(f"Error: config file not found: {config_path}", file=sys.stderr)
        sys.exit(1)
    cfg = load_config(config_path)

    outdir = Path(cfg['outdir'])
    bam_dir = outdir / 'bam'
    sample_id = resolve_sample_id(bam_dir)

    compare_tsv = outdir / f"{sample_id}.compare.txt"
    plx_bam = str(bam_dir / f"{sample_id}.plx.nsorted.bam")
    mm2_bam = str(bam_dir / f"{sample_id}.mm2.nsorted.bam")
    curation_tsv = outdir / 'curation.tsv'

    print(f"Loading comparison results from {compare_tsv}", file=sys.stderr)
    disagreeing, agreeing, clipping = load_disagreements(compare_tsv)
    print(f"  {len(disagreeing)} disagreements, {len(agreeing)} agreements, {len(clipping)} clip-only", file=sys.stderr)

    all_reads = disagreeing | agreeing | clipping
    print(f"Fingerprinting {len(all_reads)} reads from parallax BAM...", file=sys.stderr)
    plx_fps = fingerprints_from_bam(plx_bam)
    print(f"Fingerprinting {len(all_reads)} reads from minimap2 BAM...", file=sys.stderr)
    mm2_fps = fingerprints_from_bam(mm2_bam)

    existing = load_curation(curation_tsv)
    print(f"Loaded {len(existing)} existing curation entries", file=sys.stderr)

    added = reset = agreed = kept = 0

    for read in disagreeing:
        plx_fp = plx_fps.get(read, '')
        mm2_fp = mm2_fps.get(read, '')

        if read not in existing:
            existing[read] = {
                'read_name': read,
                'verdict': 'uncurated',
                'plx_fingerprint': plx_fp,
                'mm2_fingerprint': mm2_fp,
                'notes': '',
            }
            added += 1
        else:
            row = existing[read]
            old_plx_fp = row['plx_fingerprint']
            row['plx_fingerprint'] = plx_fp
            row['mm2_fingerprint'] = mm2_fp
            if row['verdict'] == 'agree':
                row['verdict'] = 'uncurated'
                reset += 1
            elif old_plx_fp != plx_fp and row['verdict'] != 'uncurated':
                row['verdict'] = 'uncurated'
                reset += 1
            else:
                kept += 1

    for read in agreeing:
        plx_fp = plx_fps.get(read, '')
        mm2_fp = mm2_fps.get(read, '')
        if read in existing:
            existing[read]['verdict'] = 'agree'
            existing[read]['plx_fingerprint'] = plx_fp
            existing[read]['mm2_fingerprint'] = mm2_fp
        else:
            existing[read] = {
                'read_name': read,
                'verdict': 'agree',
                'plx_fingerprint': plx_fp,
                'mm2_fingerprint': mm2_fp,
                'notes': '',
            }
        agreed += 1

    clipped = 0
    for read in clipping:
        plx_fp = plx_fps.get(read, '')
        mm2_fp = mm2_fps.get(read, '')
        if read in existing and existing[read]['verdict'] not in ('uncurated', 'clip'):
            # Keep a manually-set verdict (mm2/plx/neither) even if now detected as clip
            existing[read]['plx_fingerprint'] = plx_fp
            existing[read]['mm2_fingerprint'] = mm2_fp
            kept += 1
        else:
            if read in existing:
                existing[read]['verdict'] = 'clip'
                existing[read]['plx_fingerprint'] = plx_fp
                existing[read]['mm2_fingerprint'] = mm2_fp
            else:
                existing[read] = {
                    'read_name': read,
                    'verdict': 'clip',
                    'plx_fingerprint': plx_fp,
                    'mm2_fingerprint': mm2_fp,
                    'notes': '',
                }
            clipped += 1

    save_curation(curation_tsv, existing)

    uncurated = sum(1 for r in existing.values() if r['verdict'] == 'uncurated')
    print(f"Merge complete:", file=sys.stderr)
    print(f"  {added} new disagreements added as uncurated", file=sys.stderr)
    print(f"  {reset} existing entries reset to uncurated (parallax changed)", file=sys.stderr)
    print(f"  {kept} existing curated entries retained", file=sys.stderr)
    print(f"  {agreed} reads now agree", file=sys.stderr)
    print(f"  {clipped} reads auto-classified as clip-only", file=sys.stderr)
    print(f"  {uncurated} reads awaiting curation", file=sys.stderr)

    # Pre-generate per-read BAMs for all reads that may need UI review
    uncurated_reads = [r['read_name'] for r in existing.values()
                       if r['verdict'] in ('uncurated', 'clip')]
    curation_plx_bam = str(bam_dir / f"{sample_id}.curation.plx.sorted.bam")
    full_plx_bam = str(bam_dir / f"{sample_id}.plx.nsorted.bam")
    curation_seeds_bam = str(bam_dir / f"{sample_id}.curation.seeds.nsorted.bam")
    reads_zip = outdir / 'reads.zip'
    generate_per_read_bams(uncurated_reads, curation_plx_bam, full_plx_bam, curation_seeds_bam, reads_zip)


if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    main(args)
