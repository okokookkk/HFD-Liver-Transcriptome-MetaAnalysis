#!/usr/bin/env bash
set -euo pipefail

ROOT="${HFD_PROJECT_ROOT:-$(pwd)}"
cd "$ROOT"

QUEUE=00_metadata/processing_production/HFD_primary64_production_queue_v1.tsv
DRIVER=scripts/30_process_one_GSE_production_v1_3_2.py
MASTER_STATUS=00_metadata/processing_production/HFD_primary64_production_status_v1.tsv

test -s "$QUEUE"
test -s "$DRIVER"

mkdir -p logs/production
mkdir -p 00_metadata/processing_production

printf "GSE\tExit_code\tProduction_status\tLog\n" > "$MASTER_STATUS"

tail -n +2 "$QUEUE" | cut -f2 | while read -r GSE
do
    test -n "$GSE"

    LOG="logs/production/${GSE}_production_driver_v1_3_2.log"

    echo
    echo "============================================================"
    echo "PRODUCTION: $GSE"
    echo "============================================================"

    set +e
    python3 "$DRIVER" \
      --gse "$GSE" \
      --download-workers 4 \
      --fastqc-threads 16 \
      --star-threads 16 \
      --fc-threads 4 \
      --cleanup \
      2>&1 | tee "$LOG"

    rc=${PIPESTATUS[0]}
    set -e

    case "$rc" in
        0)
            status="PASS_OR_ALREADY_COMPLETE"
            ;;
        20)
            status="PROCESSING_QC_EXCEPTION_REPRESENTATIVE_READLEN_LT50"
            ;;
        21)
            status="PROCESSING_QC_EXCEPTION_R1R2_COUNT_MISMATCH"
            ;;
        22)
            status="REVIEW_REQUIRED_ADAPTER"
            ;;
        23)
            status="REVIEW_REQUIRED_AMBIGUOUS_STRANDEDNESS"
            ;;
        24)
            status="REVIEW_REQUIRED_WITHIN_GSM_STRAND_DISCORDANCE"
            ;;
        25)
            status="SRR_TO_GSM_AGGREGATION_FAIL"
            ;;
        26)
            status="REVIEW_REQUIRED_READ_LENGTH_BORDERLINE"
            ;;
        27)
            status="DEFERRED_ACQUISITION_FAILURE"
            ;;
        *)
            status="UNEXPECTED_FATAL_ERROR"
            printf "%s\t%s\t%s\t%s\n" "$GSE" "$rc" "$status" "$LOG" >> "$MASTER_STATUS"
            echo "Unexpected fatal error for $GSE (exit=$rc). Stopping master production."
            exit "$rc"
            ;;
    esac

    printf "%s\t%s\t%s\t%s\n" "$GSE" "$rc" "$status" "$LOG" >> "$MASTER_STATUS"

    if [[ "$rc" -ne 0 ]]; then
        echo "$GSE recorded as $status; continuing to next frozen GSE."
    fi

done

echo
echo "HFD PRIMARY64 PRODUCTION LOOP: COMPLETE"
echo "Status table: $MASTER_STATUS"
