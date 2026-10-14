#!/usr/bin/env bash
set -euo pipefail

OUT="00_metadata/systematic_search"
mkdir -p "$OUT"

RAW="$OUT/geo_search_strict_all_hits_v2.tsv"
UNIQ="$OUT/geo_candidate_GSE_strict_v2.tsv"
LOG="$OUT/search_queries_strict_v2.tsv"

: > "$RAW"
: > "$LOG"

queries=(

'Mus musculus[Organism] AND liver[Sample Source] AND "high fat diet"[Description] AND 6:10000[Number of Samples] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]'

'Mus musculus[Organism] AND liver[Sample Source] AND "high-fat diet"[Description] AND 6:10000[Number of Samples] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]'

'Mus musculus[Organism] AND liver[Sample Source] AND HFD[Description] AND 6:10000[Number of Samples] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]'

'Mus musculus[Organism] AND liver[Sample Source] AND "diet induced obesity"[Description] AND 6:10000[Number of Samples] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]'

'Mus musculus[Organism] AND liver[Sample Source] AND "diet-induced obesity"[Description] AND 6:10000[Number of Samples] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]'

'Mus musculus[Organism] AND liver[Sample Source] AND DIO[Description] AND 6:10000[Number of Samples] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]'

'Mus musculus[Organism] AND liver[Description] AND "high fat diet"[Description] AND 6:10000[Number of Samples] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]'

'Mus musculus[Organism] AND hepatic[Description] AND "high fat diet"[Description] AND 6:10000[Number of Samples] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]'
)

i=0

for q in "${queries[@]}"; do

    i=$((i+1))

    printf "%02d\t%s\n" "$i" "$q" >> "$LOG"

    echo "========================================"
    echo "STRICT QUERY $i"
    echo "$q"
    echo "========================================"

    esearch \
        -db gds \
        -query "$q" \
    | esummary \
    | xtract \
        -pattern DocumentSummary \
        -element \
            DocumentSummary/Accession \
            DocumentSummary/title \
            DocumentSummary/summary \
    | awk -F'\t' -v Q="$i" '
        BEGIN {
            OFS="\t"
        }

        $1 ~ /^GSE[0-9]+$/ {
            print Q,$1,$2,$3
        }
    ' \
    >> "$RAW"

done


## ============================================================
## Deduplicate at GSE level
##
## RAW columns:
## 1 Query_ID
## 2 GSE
## 3 Series_title
## 4 Series_summary
##
## UNIQ columns:
## 1 GSE
## 2 Series_title
## 3 Series_summary
## ============================================================

awk -F'\t' '
BEGIN {
    OFS="\t"
}

$2 ~ /^GSE[0-9]+$/ {

    accession=$2
    title=$3
    summary=$4

    if (!(accession in seen)) {

        seen[accession]=1

        print \
            accession, \
            title, \
            summary
    }
}
' "$RAW" \
| LC_ALL=C sort -t$'\t' -k1,1 \
> "$UNIQ"


echo
echo "========================================"
echo "STRICT DISCOVERY V2 COMPLETE"
echo "========================================"

echo -n "Raw GSE query hits : "
wc -l < "$RAW"

echo -n "Unique GSE records : "
wc -l < "$UNIQ"

echo
echo "Output:"
echo "  $RAW"
echo "  $UNIQ"
echo "  $LOG"
