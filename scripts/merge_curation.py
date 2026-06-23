"""Merge new comparison output into the curation database.

Usage:
    merge_curation.py <compare_tsv> <plx_bam> <mm2_bam> <curation_tsv>

Arguments:
    <compare_tsv>    Output of compare_alignments_2.py (tsv with header)
    <plx_bam>        Name-sorted parallax BAM (for fingerprinting)
    <mm2_bam>        Name-sorted minimap2 BAM (for fingerprinting)
    <curation_tsv>   Curation database (created if absent, updated in place)

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

from fingerprint import fingerprints_from_bam

COLUMNS = ['read_name', 'verdict', 'plx_fingerprint', 'mm2_fingerprint', 'notes']
VERDICTS = {'uncurated', 'mm2', 'plx', 'neither', 'agree'}


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


def main(args):
    compare_tsv = Path(args['<compare_tsv>'])
    plx_bam = args['<plx_bam>']
    mm2_bam = args['<mm2_bam>']
    curation_tsv = Path(args['<curation_tsv>'])

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
                # Was previously agreeing, now disagrees again
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


if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    main(args)
