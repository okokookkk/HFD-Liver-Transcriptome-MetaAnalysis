#!/usr/bin/env python3
from pathlib import Path
import gzip
import re
import pandas as pd

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()

GTF = (
    ROOT
    / "00_metadata/reference/GENCODE_M39_GRCm39/"
      "gencode.vM39.primary_assembly.annotation.gtf"
)

STABLE = (
    ROOT
    / "05_meta/stable_core_v1/"
      "HFD_stable_core_final_v1.tsv"
)

META = (
    ROOT
    / "05_meta/gene_level_v1_2/"
      "HFD_gene_multilevel_meta_primary_v1_2.tsv.gz"
)

ANN_OUT = (
    ROOT
    / "00_metadata/reference/GENCODE_M39_GRCm39/"
      "gencode.vM39_gene_annotation_v1.tsv"
)

STABLE_OUT = (
    ROOT
    / "05_meta/stable_core_v1/"
      "HFD_stable_core_final_annotated_v1.tsv"
)

META_OUT = (
    ROOT
    / "05_meta/gene_level_v1_2/"
      "HFD_gene_multilevel_meta_primary_annotated_v1_2.tsv.gz"
)

SUMMARY = (
    ROOT
    / "00_metadata/analysis_outputs/stable_core_v1/"
      "HFD_stable_core_annotation_summary_v1.txt"
)

for q in [GTF, STABLE, META]:
    if not q.exists():
        raise RuntimeError(f"Missing required input: {q}")

def parse_attrs(s):
    out = {}
    # GTF attributes: key "value"; key2 "value2";
    for m in re.finditer(r'(\S+)\s+"([^"]*)"', s):
        out[m.group(1)] = m.group(2)
    return out

rows = []
with GTF.open() as fh:
    for line in fh:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 9:
            continue
        chrom, source, feature, start, end, score, strand, frame, attrs = parts
        if feature != "gene":
            continue

        a = parse_attrs(attrs)
        gene_id = a.get("gene_id", "")
        if not gene_id:
            continue

        gene_name = a.get("gene_name", "")
        gene_type = a.get("gene_type", a.get("gene_biotype", ""))

        rows.append({
            "Geneid": gene_id,
            "Geneid_no_version": gene_id.split(".")[0],
            "gene_symbol": gene_name,
            "gene_type": gene_type,
            "chromosome": chrom,
            "gene_start": int(start),
            "gene_end": int(end),
            "strand": strand,
            "source": source,
        })

ann = pd.DataFrame(rows)

if ann.empty:
    raise RuntimeError("No gene records parsed from GTF")
if ann["Geneid"].duplicated().any():
    dup = ann.loc[ann["Geneid"].duplicated(keep=False), "Geneid"].unique()
    raise RuntimeError(f"Duplicate versioned Geneid in GTF annotation: {dup[:10]}")

ann = ann.sort_values(["chromosome", "gene_start", "Geneid"]).reset_index(drop=True)
ann.to_csv(ANN_OUT, sep="\t", index=False)

# Secondary fallback by version-stripped ID, only when unique in the GTF.
ann_nv_counts = ann.groupby("Geneid_no_version").size()
unique_nv = set(ann_nv_counts[ann_nv_counts == 1].index)
ann_nv = (
    ann[ann["Geneid_no_version"].isin(unique_nv)]
    .drop(columns="Geneid")
    .rename(columns={"Geneid_no_version": "_join_no_version"})
)

def annotate_table(x, label):
    if "Geneid" not in x.columns:
        raise RuntimeError(f"{label}: Geneid column missing")
    if x["Geneid"].duplicated().any():
        raise RuntimeError(f"{label}: duplicate Geneid")

    y = x.merge(
        ann,
        on="Geneid",
        how="left",
        validate="one_to_one",
        indicator="_versioned_join"
    )

    missing = y["gene_symbol"].isna()

    if missing.any():
        fallback = y.loc[missing, ["Geneid"]].copy()
        fallback["_join_no_version"] = fallback["Geneid"].str.split(".").str[0]
        fallback = fallback.merge(
            ann_nv,
            on="_join_no_version",
            how="left",
            validate="many_to_one"
        )

        for c in [
            "gene_symbol", "gene_type", "chromosome",
            "gene_start", "gene_end", "strand", "source"
        ]:
            y.loc[missing, c] = fallback[c].values

    y["Annotation_status"] = "MATCHED"
    y.loc[y["gene_symbol"].isna(), "Annotation_status"] = "UNMATCHED"

    # Move annotation fields immediately after Geneid.
    annotation_cols = [
        "Geneid",
        "gene_symbol",
        "gene_type",
        "chromosome",
        "gene_start",
        "gene_end",
        "strand",
        "Annotation_status",
    ]
    remaining = [
        c for c in y.columns
        if c not in annotation_cols
        and c not in {
            "Geneid_no_version", "source", "_versioned_join"
        }
    ]
    y = y[annotation_cols + remaining]

    return y

stable = pd.read_csv(STABLE, sep="\t", dtype={"Geneid": str})
meta = pd.read_csv(META, sep="\t", compression="gzip", dtype={"Geneid": str})

stable_a = annotate_table(stable, "stable_core")
meta_a = annotate_table(meta, "primary_meta")

stable_a.to_csv(STABLE_OUT, sep="\t", index=False)
meta_a.to_csv(META_OUT, sep="\t", index=False, compression="gzip")

stable_matched = int((stable_a["Annotation_status"] == "MATCHED").sum())
meta_matched = int((meta_a["Annotation_status"] == "MATCHED").sum())

if len(stable_a) != 59:
    raise RuntimeError(f"Expected 59 stable-core genes, observed {len(stable_a)}")
if stable_matched != len(stable_a):
    bad = stable_a.loc[
        stable_a["Annotation_status"] != "MATCHED",
        ["Geneid", "Annotation_status"]
    ]
    raise RuntimeError(
        "Not all stable-core genes were annotated:\n"
        + bad.to_string(index=False)
    )

# Direction counts and gene-type composition.
stable_a["Direction"] = stable_a["pooled_log2FC"].apply(
    lambda z: "UP" if z > 0 else ("DOWN" if z < 0 else "ZERO")
)

type_counts = (
    stable_a["gene_type"]
    .fillna("UNANNOTATED")
    .replace("", "UNANNOTATED")
    .value_counts()
)

summary = [
    "===== HFD STABLE-CORE GENCODE M39 ANNOTATION SUMMARY v1 =====",
    "",
    "STATUS=PASS",
    f"GENCODE_gene_records={len(ann)}",
    f"Stable_core_genes={len(stable_a)}",
    f"Stable_core_annotated={stable_matched}/{len(stable_a)}",
    f"Primary_meta_genes={len(meta_a)}",
    f"Primary_meta_annotated={meta_matched}/{len(meta_a)}",
    "",
    f"Stable_core_UP={(stable_a['Direction'] == 'UP').sum()}",
    f"Stable_core_DOWN={(stable_a['Direction'] == 'DOWN').sum()}",
    "",
    "STABLE_CORE_GENE_TYPES",
]

for k, v in type_counts.items():
    summary.append(f"{k}={v}")

summary += [
    "",
    f"Gene_annotation={ANN_OUT.relative_to(ROOT)}",
    f"Stable_core_annotated={STABLE_OUT.relative_to(ROOT)}",
    f"Primary_meta_annotated={META_OUT.relative_to(ROOT)}",
    "",
    "GENE_IDENTITIES_CHANGED=NO",
    "META_ANALYSIS_RESULTS_CHANGED=NO",
]

SUMMARY.parent.mkdir(parents=True, exist_ok=True)
SUMMARY.write_text("\n".join(summary) + "\n", encoding="utf-8")
print("\n".join(summary))

print("\n===== FINAL STABLE CORE: SYMBOL + EFFECT =====")
show = stable_a[
    [
        "gene_symbol",
        "gene_type",
        "Geneid",
        "pooled_log2FC",
        "pooled_FDR_BH",
        "N_independent_GSE",
        "I2_total_percent",
        "prediction95_lower_new_GSE",
        "prediction95_upper_new_GSE",
        "GSE_balanced_direction_consistency",
    ]
].sort_values("pooled_log2FC", ascending=False)

print(show.to_string(index=False))

