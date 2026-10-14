#!/usr/bin/env bash
set -euo pipefail

OUT="00_metadata/systematic_search"
mkdir -p "$OUT"

RAW="$OUT/geo_search_strict_all_hits.tsv"
UNIQ="$OUT/geo_candidate_GSE_strict.tsv"
LOG="$OUT/search_queries_strict.tsv"

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
        -element Accession Title \
    | awk -v Q="$i" 'BEGIN{OFS="\t"} {print Q,$0}' \
    >> "$RAW"

done


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
| LC_ALL=C sort -t$'\t' -k1,1 \
> "$UNIQ"


echo
echo "========================================"
echo "STRICT DISCOVERY COMPLETE"
echo "========================================"

echo -n "Raw strict hits      : "
wc -l < "$RAW"

echo -n "Unique strict GSE    : "
wc -l < "$UNIQ"

echo
echo "Output:"
echo "  $RAW"
echo "  $UNIQ"
echo "  $LOG"
