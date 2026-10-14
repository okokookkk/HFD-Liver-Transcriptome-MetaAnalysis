#!/usr/bin/env bash
set -euo pipefail

OUT="00_metadata/systematic_search"
mkdir -p "$OUT"

RAW="$OUT/geo_search_all_hits.tsv"
UNIQ="$OUT/geo_candidate_GSE.tsv"
LOG="$OUT/search_queries.tsv"

: > "$RAW"
: > "$LOG"

queries=(
'Mus musculus[Organism] AND liver[All Fields] AND "high fat diet"[All Fields] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]'
'Mus musculus[Organism] AND liver[All Fields] AND "high-fat diet"[All Fields] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]'
'Mus musculus[Organism] AND liver[All Fields] AND HFD[All Fields] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]'
'Mus musculus[Organism] AND liver[All Fields] AND "diet induced obesity"[All Fields] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]'
'Mus musculus[Organism] AND liver[All Fields] AND "diet-induced obesity"[All Fields] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]'
'Mus musculus[Organism] AND liver[All Fields] AND DIO[All Fields] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]'
'Mus musculus[Organism] AND hepatic[All Fields] AND "high fat"[All Fields] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]'
'Mus musculus[Organism] AND liver[All Fields] AND "fatty liver"[All Fields] AND "Expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]'
)

i=0

for q in "${queries[@]}"; do

    i=$((i+1))

    printf "%02d\t%s\n" "$i" "$q" >> "$LOG"

    echo "========================================"
    echo "QUERY $i"
    echo "$q"
    echo "========================================"

    esearch \
        -db gds \
        -query "$q" \
    | esummary \
    | xtract \
        -pattern DocumentSummary \
        -element Accession Title \
    | awk -v Q="$i" 'BEGIN{OFS="\t"} {print Q,$0}' \
    >> "$RAW"

done


## ------------------------------------------------------------
## Retain GSE records and deduplicate
## ------------------------------------------------------------

awk -F'\t' '
BEGIN { OFS="\t" }
$2 ~ /^GSE[0-9]+$/ {
    accession=$2
    title=$3

    if (!(accession in seen)) {
        seen[accession]=1
        print accession,title
    }
}
' "$RAW" \
| sort -V \
> "$UNIQ"


echo
echo "========================================"
echo "DISCOVERY COMPLETE"
echo "========================================"

echo -n "Raw query hits     : "
wc -l < "$RAW"

echo -n "Unique GSE records : "
wc -l < "$UNIQ"

echo
echo "Output:"
echo "  $RAW"
echo "  $UNIQ"
echo "  $LOG"
