#!/usr/bin/env python3
from pathlib import Path
import gzip
import hashlib
import pandas as pd
from collections import defaultdict

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()

INDEX = ROOT / "00_metadata/analysis_outputs/contrast_de_v1/HFD_primary105_DESeq2_effect_file_index_v1.tsv"
STATUS = ROOT / "00_metadata/analysis_outputs/contrast_de_v1/HFD_primary105_DESeq2_status_v1.tsv"
CONTRAST_AUDIT = ROOT / "00_metadata/analysis_inputs/frozen_v1/HFD_primary105_contrast_analysis_input_audit_v1_1.tsv"

OUTDIR = ROOT / "00_metadata/analysis_inputs/meta_analysis_v1"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_LONG = OUTDIR / "HFD_primary105_meta_effects_primaryQC_v1.tsv.gz"
OUT_GENE = OUTDIR / "HFD_gene_meta_eligibility_v1.tsv"
OUT_ELIGIBLE = OUTDIR / "HFD_gene_meta_eligible_ge10GSE_v1.tsv"
OUT_COVERAGE = OUTDIR / "HFD_gene_meta_GSE_coverage_distribution_v1.tsv"
OUT_SUMMARY = OUTDIR / "HFD_primary105_meta_input_summary_v1.txt"

for q in [INDEX, STATUS, CONTRAST_AUDIT]:
    if not q.exists():
        raise RuntimeError(f"Missing required input: {q}")

idx = pd.read_csv(INDEX, sep="\t", dtype=str).fillna("")
status = pd.read_csv(STATUS, sep="\t", dtype=str).fillna("")
audit = pd.read_csv(CONTRAST_AUDIT, sep="\t", dtype=str).fillna("")

if len(idx) != 105 or idx["Contrast_id"].nunique() != 105:
    raise RuntimeError("Effect-file index must contain 105 unique contrasts")
if len(status) != 105 or status["Contrast_id"].nunique() != 105:
    raise RuntimeError("DESeq2 status must contain 105 unique contrasts")
if not (status["Status"] == "PASS").all():
    raise RuntimeError("Not all contrast-level DESeq2 runs are PASS")
if len(audit) != 105 or audit["Contrast_id"].nunique() != 105:
    raise RuntimeError("Contrast audit must contain 105 unique contrasts")

# Contrast-level moderator aliases already copied from canonical *_final values.
moderators = [
    "hfd_duration_weeks",
    "sex",
    "strain",
    "hfd_fat_kcal_percent",
    "control_diet_type",
    "collection_ZT",
]
mods_present = [m for m in moderators if m in audit.columns]

contrast_meta_cols = ["Contrast_id", "GSE"] + mods_present
contrast_meta = audit[contrast_meta_cols].drop_duplicates("Contrast_id")
contrast_meta = contrast_meta.set_index("Contrast_id")

# Global audit containers. Sets are kept only for ~20k genes and are tractable.
all_contrasts = defaultdict(set)
all_gses = defaultdict(set)
pass_contrasts = defaultdict(set)
pass_gses = defaultdict(set)
cooks_fail_contrasts = defaultdict(set)
positive_pass = defaultdict(int)
negative_pass = defaultdict(int)
zero_pass = defaultdict(int)

n_all_effects = 0
n_primary_pass = 0
n_cooks_flag = 0
n_nonfinite = 0

selected_effect_cols = [
    "Geneid",
    "log2FoldChange",
    "lfcSE",
    "sampling_variance",
    "baseMean",
    "mean_normalized_CONTROL",
    "mean_normalized_HFD",
    "N_samples_count_ge_10",
    "max_Cooks_distance",
    "DESeq2_Cooks_cutoff",
    "DESeq2_Cooks_flag",
    "Finite_effect_SE",
    "Primary_effect_QC",
]

# Write primary-QC long table incrementally to avoid loading ~1.5M effects at once.
tmp_long = Path(str(OUT_LONG) + ".tmp")
if tmp_long.exists():
    tmp_long.unlink()

header_written = False

with gzip.open(tmp_long, "wt") as gzout:
    for i, r in idx.sort_values(["GSE", "Contrast_id"]).iterrows():
        cid = r["Contrast_id"]
        gse = r["GSE"]

        if cid not in contrast_meta.index:
            raise RuntimeError(f"{cid}: absent from contrast audit")
        if contrast_meta.loc[cid, "GSE"] != gse:
            raise RuntimeError(f"{cid}: GSE mismatch between index and contrast audit")

        f = ROOT / r["Effect_file"]
        done = ROOT / r["DONE_file"]

        if not f.exists():
            raise RuntimeError(f"{cid}: missing effect file {f}")
        if not done.exists():
            raise RuntimeError(f"{cid}: missing DONE file {done}")

        x = pd.read_csv(f, sep="\t", compression="gzip")

        missing = [c for c in selected_effect_cols if c not in x.columns]
        if missing:
            raise RuntimeError(f"{cid}: effect file missing columns {missing}")

        if x["Geneid"].duplicated().any():
            raise RuntimeError(f"{cid}: duplicate Geneid in effect file")

        n_all_effects += len(x)

        finite_mask = (
            pd.to_numeric(x["log2FoldChange"], errors="coerce").notna()
            & pd.to_numeric(x["lfcSE"], errors="coerce").notna()
            & (pd.to_numeric(x["lfcSE"], errors="coerce") > 0)
        )
        n_nonfinite += int((~finite_mask).sum())

        cooks_mask = x["DESeq2_Cooks_flag"].astype(str).eq("YES")
        n_cooks_flag += int(cooks_mask.sum())

        # All modeled effects contribute to coverage diagnostics.
        for gene in x["Geneid"].astype(str):
            all_contrasts[gene].add(cid)
            all_gses[gene].add(gse)

        for gene in x.loc[cooks_mask, "Geneid"].astype(str):
            cooks_fail_contrasts[gene].add(cid)

        p = x[x["Primary_effect_QC"].astype(str).eq("PASS")].copy()
        n_primary_pass += len(p)

        if len(p):
            for gene, yi in zip(
                p["Geneid"].astype(str),
                pd.to_numeric(p["log2FoldChange"], errors="raise")
            ):
                pass_contrasts[gene].add(cid)
                pass_gses[gene].add(gse)
                if yi > 0:
                    positive_pass[gene] += 1
                elif yi < 0:
                    negative_pass[gene] += 1
                else:
                    zero_pass[gene] += 1

            out = p[selected_effect_cols].copy()
            out.insert(0, "GSE", gse)
            out.insert(0, "Contrast_id", cid)

            for m in mods_present:
                out[m] = contrast_meta.loc[cid, m]

            out.to_csv(
                gzout,
                sep="\t",
                index=False,
                header=not header_written,
                na_rep="NA"
            )
            header_written = True

if OUT_LONG.exists():
    OUT_LONG.unlink()
tmp_long.rename(OUT_LONG)

# Build per-gene coverage/eligibility table.
genes = sorted(all_contrasts.keys())

rows = []
for gene in genes:
    n_all_c = len(all_contrasts[gene])
    n_all_g = len(all_gses[gene])
    n_pass_c = len(pass_contrasts[gene])
    n_pass_g = len(pass_gses[gene])
    n_cooks_c = len(cooks_fail_contrasts[gene])

    rows.append({
        "Geneid": gene,
        "N_modeled_contrasts": n_all_c,
        "N_modeled_independent_GSE": n_all_g,
        "N_primaryQC_PASS_contrasts": n_pass_c,
        "N_primaryQC_PASS_independent_GSE": n_pass_g,
        "N_Cooks_flagged_contrasts": n_cooks_c,
        "N_positive_primaryQC_effects": positive_pass[gene],
        "N_negative_primaryQC_effects": negative_pass[gene],
        "N_zero_primaryQC_effects": zero_pass[gene],
        "Primary_meta_eligible_ge10_independent_GSE": "YES" if n_pass_g >= 10 else "NO",
    })

gene = pd.DataFrame(rows)
gene = gene.sort_values(
    ["N_primaryQC_PASS_independent_GSE", "N_primaryQC_PASS_contrasts", "Geneid"],
    ascending=[False, False, True]
)
gene.to_csv(OUT_GENE, sep="\t", index=False)

eligible = gene[gene["Primary_meta_eligible_ge10_independent_GSE"] == "YES"].copy()
eligible.to_csv(OUT_ELIGIBLE, sep="\t", index=False)

coverage = (
    gene["N_primaryQC_PASS_independent_GSE"]
    .value_counts()
    .sort_index()
    .rename_axis("N_primaryQC_PASS_independent_GSE")
    .reset_index(name="N_genes")
)
coverage["Cumulative_genes_ge_this_GSE"] = coverage.iloc[::-1]["N_genes"].cumsum().iloc[::-1].values
coverage.to_csv(OUT_COVERAGE, sep="\t", index=False)

# Hash the actual primary meta-analysis long table for downstream identity checks.
sha = hashlib.sha256()
with open(OUT_LONG, "rb") as fh:
    for block in iter(lambda: fh.read(1024 * 1024), b""):
        sha.update(block)
long_hash = sha.hexdigest()

def n_ge(k):
    return int((gene["N_primaryQC_PASS_independent_GSE"] >= k).sum())

summary = [
    "===== HFD PRIMARY105 META-ANALYSIS INPUT SUMMARY v1 =====",
    "",
    "FINAL_STATUS=PASS",
    "Contrasts=105",
    "Independent_GSE=64",
    f"Genes_modeled_in_at_least_one_contrast={len(gene)}",
    f"All_modeled_contrast_gene_effects={n_all_effects}",
    f"Nonfinite_or_nonpositive_SE_effects={n_nonfinite}",
    f"Cooks_flagged_effects={n_cooks_flag}",
    f"Primary_effect_QC_PASS_effects={n_primary_pass}",
    "",
    "PRIMARY_META_ELIGIBILITY=min_10_independent_GSE_with_primary_QC_PASS_effect",
    f"Genes_ge_5_independent_GSE={n_ge(5)}",
    f"Genes_ge_10_independent_GSE={n_ge(10)}",
    f"Genes_ge_20_independent_GSE={n_ge(20)}",
    f"Genes_ge_30_independent_GSE={n_ge(30)}",
    f"Genes_ge_40_independent_GSE={n_ge(40)}",
    f"Genes_ge_50_independent_GSE={n_ge(50)}",
    f"Genes_ge_60_independent_GSE={n_ge(60)}",
    "",
    f"Primary_meta_eligible_genes={len(eligible)}",
    f"Primary_effect_long_table_SHA256={long_hash}",
    "",
    f"Primary_effect_long_table={OUT_LONG.relative_to(ROOT)}",
    f"Gene_eligibility={OUT_GENE.relative_to(ROOT)}",
    f"Eligible_gene_list={OUT_ELIGIBLE.relative_to(ROOT)}",
    f"GSE_coverage_distribution={OUT_COVERAGE.relative_to(ROOT)}",
    "",
    "WITHIN_CONTRAST_PVALUE_USED_FOR_ELIGIBILITY=NO",
    "COOKS_FLAGGED_EFFECTS_INCLUDED_IN_PRIMARY_META_INPUT=NO",
    "SHARED_GSM_ACROSS_CONTRASTS=0",
    "META_ANALYSIS_PERFORMED=NO",
]

OUT_SUMMARY.write_text("\n".join(summary) + "\n", encoding="utf-8")
print("\n".join(summary))

