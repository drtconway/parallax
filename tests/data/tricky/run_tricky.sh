#!/bin/bash
#
# Diagnostic harness for tricky alignment cases.
#
# Runs parallax on each .fq file in this directory, optionally using a
# per-case TOML config (e.g. del_softclip.toml) for debug output.
#
# Usage:
#   ./run_tricky.sh                  # run all cases
#   ./run_tricky.sh del_softclip     # run one case by name
#
# Environment overrides:
#   PARALLAX_REF   path to hg38 reference FASTA
#   PARALLAX_IDX   path to parallax index directory
#   PARALLAX_BIN   path to parallax binary
#
set -euo pipefail

DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT="$(cd "$DIR/../../.." && pwd)"
REF="${PARALLAX_REF:-/Users/tom.conway/data/hg38/hg38_primary.fasta}"
IDX="${PARALLAX_IDX:-${PROJECT}/hg38_idx}"
PLX="${PARALLAX_BIN:-${PROJECT}/target/release/parallax}"

# Optional filter: run only one case
FILTER="${1:-}"

# ── Helpers ──────────────────────────────────────────────────────────────────

show_alignments() {
    local sam="$1"
    grep -v '^@' "$sam" | while IFS=$'\t' read -r qname flag rname pos mapq cigar rest; do
        printf "    flag=%-5s  chr=%-6s  pos=%-12s  mapq=%-4s  cigar=%s\n" \
            "$flag" "$rname" "$pos" "$mapq" "$cigar"
    done
}

parse_expected() {
    # Extract expected segments from read name (@sim_NNNN:seg1,seg2,...)
    local fq="$1"
    head -1 "$fq" | sed 's/^@[^:]*://' | tr ',' '\n' | while read -r seg; do
        local chr pos_start pos_end strand
        chr=$(echo "$seg" | cut -d_ -f1)
        pos_start=$(echo "$seg" | rev | cut -d_ -f3 | rev)
        pos_end=$(echo "$seg" | rev | cut -d_ -f2 | rev)
        strand=$(echo "$seg" | rev | cut -d_ -f1 | rev)
        printf "    %s:%s-%s (%s)\n" "$chr" "$pos_start" "$pos_end" "$strand"
    done
}

# ── Config generation ────────────────────────────────────────────────────────

ensure_config() {
    local fq="$1"
    local label
    label=$(basename "${fq%.fq}")
    local toml="${fq%.fq}.toml"

    if [[ -f "$toml" ]]; then
        return
    fi

    cat > "$toml" <<EOF
# Config for tricky case: ${label}
# Edit to enable debug outputs or tune parameters.

[seeding]
debug_seeds_sam = "${DIR}/${label}_seeds.sam"
debug_chains_tsv = "${DIR}/${label}_chains.tsv"
debug_gap_fills_tsv = "${DIR}/${label}_gaps.tsv"
debug_split_decisions_tsv = "${DIR}/${label}_splits.tsv"
EOF
    echo "  Created ${label}.toml (edit to customise)"
}

# ── Main loop ────────────────────────────────────────────────────────────────

run_case() {
    local fq="$1"
    local label
    label=$(basename "${fq%.fq}")
    local sam="${fq%.fq}.sam"
    local toml="${fq%.fq}.toml"

    # Create config if missing
    ensure_config "$fq"

    echo "═══════════════════════════════════════════════════════════════"
    echo "  $label"
    echo "═══════════════════════════════════════════════════════════════"
    echo ""
    echo "  Expected segments:"
    parse_expected "$fq"
    echo ""

    # Run parallax
    "$PLX" align "$REF" "$fq" "$sam" -x "$IDX" -t 1 -c "$toml" 2>&1 \
        | grep -v '^\[' || true   # strip log lines

    echo "  Parallax alignments:"
    show_alignments "$sam"

    local n_alns
    n_alns=$(grep -cv '^@' "$sam" || true)
    echo ""
    echo "  Total alignments: $n_alns"
    echo ""
}

count=0
for fq in "$DIR"/*.fq; do
    label=$(basename "${fq%.fq}")
    if [[ -n "$FILTER" && "$label" != "$FILTER" ]]; then
        continue
    fi
    run_case "$fq"
    count=$((count + 1))
done

if [[ $count -eq 0 ]]; then
    echo "No matching cases found${FILTER:+ for filter '$FILTER'}."
    exit 1
fi

echo "═══════════════════════════════════════════════════════════════"
echo "  Ran $count case(s). SAM files are alongside the .fq files."
echo "═══════════════════════════════════════════════════════════════"
