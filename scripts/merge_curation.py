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
import sys
from pathlib import Path

import docopt
import pysam
import yaml

from fingerprint import fingerprints_from_bam

COLUMNS = ['read_name', 'verdict', 'plx_fingerprint', 'mm2_fingerprint', 'notes']
VERDICTS = {'uncurated', 'mm2', 'plx', 'neither', 'agree'}


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


def load_disagreements(compare_tsv: Path) -> tuple[set[str], set[str]]:
    """Return (disagreeing_reads, agreeing_reads) from compare output."""
    disagreeing = set()
    agreeing = set()
    with open(compare_tsv, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            read = row['read_id']
            state = row['state']
            jaccard = float(row['jaccard'])
            if state == 'both' and jaccard >= 1.0:
                agreeing.add(read)
            else:
                disagreeing.add(read)
    return disagreeing, agreeing


def generate_per_read_bams(
    queue: list[str],
    curation_bam: str,
    seeds_bam: str,
    reads_dir: Path,
) -> None:
    """Extract per-read BAMs for all queued reads into reads_dir.

    Skips reads whose BAMs already exist and are up to date.
    """
    reads_dir.mkdir(parents=True, exist_ok=True)
    needed = [r for r in queue if not (reads_dir / f"{r.replace('/', '_')}.plx.bam").exists()]
    if not needed:
        print(f"Per-read BAMs already up to date in {reads_dir}", file=sys.stderr)
        return

    needed_set = set(needed)
    print(f"Collecting records for {len(needed)} reads from {curation_bam}...", file=sys.stderr)
    read_records: dict[str, list] = {}
    with pysam.AlignmentFile(curation_bam, 'rb') as bam:
        curation_header = bam.header
        for rec in bam:
            if rec.is_secondary:
                continue
            if rec.query_name in needed_set:
                read_records.setdefault(rec.query_name, []).append(rec)

    seed_records: dict[str, list] = {}
    seeds_header = None
    if Path(seeds_bam).exists():
        print(f"Collecting records for {len(needed)} reads from {seeds_bam}...", file=sys.stderr)
        with pysam.AlignmentFile(seeds_bam, 'rb') as bam:
            seeds_header = bam.header
            for rec in bam:
                if rec.is_secondary:
                    continue
                if rec.query_name in needed_set:
                    seed_records.setdefault(rec.query_name, []).append(rec)

    print(f"Writing {len(needed)} per-read BAMs to {reads_dir}...", file=sys.stderr)
    for read_name in needed:
        safe = read_name.replace('/', '_').replace(' ', '_')

        unsorted = reads_dir / f"{safe}.plx.unsorted.bam"
        sorted_bam = reads_dir / f"{safe}.plx.bam"
        with pysam.AlignmentFile(str(unsorted), 'wb', header=curation_header) as out:
            for rec in read_records.get(read_name, []):
                out.write(rec)
        pysam.sort('-o', str(sorted_bam), str(unsorted))
        pysam.index(str(sorted_bam))
        unsorted.unlink()

        if seeds_header and seed_records.get(read_name):
            unsorted = reads_dir / f"{safe}.seeds.unsorted.bam"
            sorted_bam = reads_dir / f"{safe}.seeds.bam"
            with pysam.AlignmentFile(str(unsorted), 'wb', header=seeds_header) as out:
                for rec in seed_records[read_name]:
                    out.write(rec)
            pysam.sort('-o', str(sorted_bam), str(unsorted))
            pysam.index(str(sorted_bam))
            unsorted.unlink()


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
    disagreeing, agreeing = load_disagreements(compare_tsv)
    print(f"  {len(disagreeing)} disagreements, {len(agreeing)} agreements", file=sys.stderr)

    all_reads = disagreeing | agreeing
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

    save_curation(curation_tsv, existing)

    uncurated = sum(1 for r in existing.values() if r['verdict'] == 'uncurated')
    print(f"Merge complete:", file=sys.stderr)
    print(f"  {added} new disagreements added as uncurated", file=sys.stderr)
    print(f"  {reset} existing entries reset to uncurated (parallax changed)", file=sys.stderr)
    print(f"  {kept} existing curated entries retained", file=sys.stderr)
    print(f"  {agreed} reads now agree", file=sys.stderr)
    print(f"  {uncurated} reads awaiting curation", file=sys.stderr)

    # Pre-generate per-read BAMs for all uncurated reads so curate.py starts instantly
    uncurated_reads = [r['read_name'] for r in existing.values() if r['verdict'] == 'uncurated']
    curation_plx_bam = str(bam_dir / f"{sample_id}.curation.plx.sorted.bam")
    curation_seeds_bam = str(bam_dir / f"{sample_id}.curation.seeds.sorted.bam")
    reads_dir = outdir / 'reads'
    generate_per_read_bams(uncurated_reads, curation_plx_bam, curation_seeds_bam, reads_dir)


if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    main(args)
