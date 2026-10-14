#!/usr/bin/env bash
set -euo pipefail

# HFD Meta-analysis frozen FASTQ downloader v1
# Run from the HFD_MetaAnalysis project root.
#
# Environment overrides:
#   PREFETCH_MAX_SIZE=200G
#   FASTERQ_THREADS=8
#   CONCURRENT_RUNS=1
#   KEEP_SRA=0

MANIFEST="${1:-00_metadata/contrasts/HFD_primary_FASTQ_download_manifest_frozen_v1.tsv}"

PREFETCH_MAX_SIZE="${PREFETCH_MAX_SIZE:-200G}"
FASTERQ_THREADS="${FASTERQ_THREADS:-8}"
CONCURRENT_RUNS="${CONCURRENT_RUNS:-1}"
KEEP_SRA="${KEEP_SRA:-0}"

SRA_DIR="01_raw/sra"
FASTQ_DIR="01_raw/fastq"
TMP_DIR="01_raw/tmp"
LOG_DIR="logs/fastq_download"
STATUS_DIR="00_metadata/fastq_download_status"

mkdir -p "$SRA_DIR" "$FASTQ_DIR" "$TMP_DIR" "$LOG_DIR" "$STATUS_DIR"

for exe in prefetch fasterq-dump vdb-validate; do
    command -v "$exe" >/dev/null 2>&1 || {
        echo "ERROR: $exe not found in PATH" >&2
        exit 1
    }
done

echo "=== TOOL VERSIONS ==="
prefetch --version || true
fasterq-dump --version || true
vdb-validate --version || true
echo

[[ -s "$MANIFEST" ]] || { echo "ERROR: manifest not found: $MANIFEST" >&2; exit 1; }

SRR_LIST="$STATUS_DIR/SRR_to_download_frozen_v1.txt"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="SRR") c=i; next} {print $c}' "$MANIFEST" \
  | sed '/^$/d' | sort -u > "$SRR_LIST"

EXPECTED=$(wc -l < "$SRR_LIST")
echo "Frozen SRR count: $EXPECTED"
if [[ "$EXPECTED" -ne 1817 ]]; then
    echo "ERROR: expected 1817 frozen SRR, found $EXPECTED" >&2
    exit 1
fi

run_one () {
    srr="$1"
    log="$LOG_DIR/${srr}.log"
    doneflag="$STATUS_DIR/${srr}.done"
    failflag="$STATUS_DIR/${srr}.failed"

    if [[ -f "$doneflag" ]]; then
        echo "[$srr] already complete"
        return 0
    fi

    rm -f "$failflag"
    exec > >(tee -a "$log") 2>&1

    echo "[$srr] START $(date --iso-8601=seconds)"

    # Prefetch is resumable.
    (
      cd "$SRA_DIR"
      prefetch --max-size "$PREFETCH_MAX_SIZE" "$srr"
    )

    # Validate local SRA package.
    vdb-validate "$SRA_DIR/$srr"

    # Convert. --split-files is safe for both paired and single-end runs.
    rm -rf "$TMP_DIR/$srr"
    mkdir -p "$TMP_DIR/$srr" "$FASTQ_DIR/$srr"

    fasterq-dump \
      "$SRA_DIR/$srr" \
      --split-files \
      -e "$FASTERQ_THREADS" \
      -t "$TMP_DIR/$srr" \
      --outdir "$FASTQ_DIR/$srr"

    # Compress generated FASTQ files.
    shopt -s nullglob
    files=("$FASTQ_DIR/$srr"/*.fastq)
    if [[ ${#files[@]} -eq 0 ]]; then
        echo "ERROR: no FASTQ produced for $srr"
        touch "$failflag"
        return 1
    fi

    if command -v pigz >/dev/null 2>&1; then
        pigz -p "$FASTERQ_THREADS" "${files[@]}"
    else
        gzip "${files[@]}"
    fi

    sha256sum "$FASTQ_DIR/$srr"/*.fastq.gz \
      > "$STATUS_DIR/${srr}.fastq.sha256"

    rm -rf "$TMP_DIR/$srr"

    if [[ "$KEEP_SRA" == "0" ]]; then
        rm -rf "$SRA_DIR/$srr"
    fi

    touch "$doneflag"
    echo "[$srr] COMPLETE $(date --iso-8601=seconds)"
}

export -f run_one
export SRA_DIR FASTQ_DIR TMP_DIR LOG_DIR STATUS_DIR PREFETCH_MAX_SIZE FASTERQ_THREADS KEEP_SRA

# Parallelism across SRR is intentionally conservative by default.
xargs -r -n1 -P "$CONCURRENT_RUNS" bash -c 'run_one "$0"' < "$SRR_LIST"

echo
echo "=== FINAL STATUS ==="
DONE=$(find "$STATUS_DIR" -maxdepth 1 -name 'SRR*.done' | wc -l)
FAILED=$(find "$STATUS_DIR" -maxdepth 1 -name 'SRR*.failed' | wc -l)
echo "Expected: $EXPECTED"
echo "Done:     $DONE"
echo "Failed:   $FAILED"

if [[ "$DONE" -ne "$EXPECTED" || "$FAILED" -ne 0 ]]; then
    echo "Download/conversion incomplete. Re-running this same script is safe." >&2
    exit 2
fi

cat "$STATUS_DIR"/SRR*.fastq.sha256 \
  > "$STATUS_DIR/HFD_primary_FASTQ_sha256_frozen_v1.tsv"

echo "All frozen SRR converted and checksummed."
