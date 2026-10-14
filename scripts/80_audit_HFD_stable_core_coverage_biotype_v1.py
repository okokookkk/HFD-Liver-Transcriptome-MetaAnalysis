#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()

STABLE = ROOT / "05_meta/stable_core_v1/HFD_stable_core_final_annotated_v1.tsv"
ELIG = ROOT / "00_metadata/analysis_inputs/meta_analysis_v1/HFD_gene_meta_eligibility_v1.tsv"
LONG = ROOT / "00_metadata/analysis_inputs/meta_analysis_v1/HFD_primary105_meta_effects_primaryQC_v1.tsv.gz"

OUTDIR = ROOT / "00_metadata/analysis_outputs/stable_core_v1"
OUT_GENE = OUTDIR / "HFD_stable_core_coverage_biotype_expression_audit_v1.tsv"
OUT_BIOTYPE = OUTDIR / "HFD_stable_core_biotype_by_GSE_coverage_v1.tsv"
OUT_SENS = OUTDIR / "HFD_stable_core_coverage_sensitivity_sets_v1.tsv"
OUT_SUMMARY = OUTDIR / "HFD_stable_core_coverage_biotype_summary_v1.txt"

for q in [STABLE, ELIG, LONG]:
    if not q.exists():
        raise RuntimeError(f"Missing required input: {q}")

s = pd.read_csv(STABLE, sep="\t", dtype={"Geneid": str})
e = pd.read_csv(ELIG, sep="\t", dtype={"Geneid": str})

if len(s) != 59 or s["Geneid"].nunique() != 59:
    raise RuntimeError(f"Expected 59 unique stable-core genes, observed rows={len(s)} unique={s['Geneid'].nunique()}")

needed_e = [
    "Geneid",
    "N_modeled_contrasts",
    "N_modeled_independent_GSE",
    "N_primaryQC_PASS_contrasts",
    "N_primaryQC_PASS_independent_GSE",
    "N_Cooks_flagged_contrasts",
]
missing = [c for c in needed_e if c not in e.columns]
if missing:
    raise RuntimeError(f"Eligibility table missing columns: {missing}")

a = s.merge(
    e[needed_e],
    on="Geneid",
    how="left",
    validate="one_to_one"
)

if a[needed_e[1:]].isna().any().any():
    raise RuntimeError("Stable-core gene missing eligibility/QC information")

# Read only stable-core rows from the long effect table.
usecols = [
    "Contrast_id", "GSE", "Geneid",
    "baseMean",
    "mean_normalized_CONTROL",
    "mean_normalized_HFD",
    "N_samples_count_ge_10",
]
chunks = []
stable_ids = set(a["Geneid"])

for ch in pd.read_csv(
    LONG,
    sep="\t",
    compression="gzip",
    usecols=usecols,
    chunksize=250000,
    dtype={"Geneid": str, "Contrast_id": str, "GSE": str},
):
    z = ch[ch["Geneid"].isin(stable_ids)]
    if len(z):
        chunks.append(z)

if not chunks:
    raise RuntimeError("No stable-core effects found in primary meta-analysis long table")

x = pd.concat(chunks, ignore_index=True)

for c in [
    "baseMean",
    "mean_normalized_CONTROL",
    "mean_normalized_HFD",
    "N_samples_count_ge_10",
]:
    x[c] = pd.to_numeric(x[c], errors="coerce")

expr = (
    x.groupby("Geneid", as_index=False)
    .agg(
        Effects_in_primary_long=("Contrast_id", "size"),
        Independent_GSE_in_primary_long=("GSE", "nunique"),
        Median_baseMean_across_contrasts=("baseMean", "median"),
        Min_baseMean_across_contrasts=("baseMean", "min"),
        Median_normalized_CONTROL=("mean_normalized_CONTROL", "median"),
        Median_normalized_HFD=("mean_normalized_HFD", "median"),
        Median_N_samples_count_ge_10=("N_samples_count_ge_10", "median"),
    )
)

a = a.merge(expr, on="Geneid", how="left", validate="one_to_one")

# Sanity check against the final meta model coverage.
bad = a[
    a["Independent_GSE_in_primary_long"].astype(int)
    != a["N_independent_GSE"].astype(int)
]
if len(bad):
    raise RuntimeError(
        "Stable-core primary-long GSE coverage differs from final meta coverage:\n"
        + bad[["Geneid","N_independent_GSE","Independent_GSE_in_primary_long"]].to_string(index=False)
    )

def coverage_bin(n):
    n = int(n)
    if n <= 14:
        return "10-14"
    if n <= 19:
        return "15-19"
    if n <= 29:
        return "20-29"
    if n <= 39:
        return "30-39"
    if n <= 49:
        return "40-49"
    if n <= 63:
        return "50-63"
    if n == 64:
        return "64"
    return "OTHER"

a["GSE_coverage_fraction"] = a["N_independent_GSE"] / 64.0
a["GSE_coverage_bin"] = a["N_independent_GSE"].map(coverage_bin)
a["Protein_coding"] = np.where(a["gene_type"].eq("protein_coding"), "YES", "NO")
a["Any_Cooks_flagged_contrast_before_primaryQC"] = np.where(
    a["N_Cooks_flagged_contrasts"] > 0, "YES", "NO"
)

# Sensitivity flags only. These DO NOT redefine the primary stable core.
for k in [20, 30, 40, 50, 60, 64]:
    a[f"Sensitivity_ge{k}_independent_GSE"] = np.where(
        a["N_independent_GSE"] >= k, "YES", "NO"
    )

a["Sensitivity_protein_coding_and_ge30GSE"] = np.where(
    (a["gene_type"] == "protein_coding") &
    (a["N_independent_GSE"] >= 30),
    "YES", "NO"
)

# Stable sort: broadest coverage, then effect magnitude.
a["_abs_effect"] = a["pooled_log2FC"].abs()
a = a.sort_values(
    ["N_independent_GSE", "_abs_effect", "Geneid"],
    ascending=[False, False, True]
).drop(columns="_abs_effect")

a.to_csv(OUT_GENE, sep="\t", index=False)

# Biotype x coverage bin table.
biotype = (
    a.groupby(["gene_type", "GSE_coverage_bin"], dropna=False)
    .size()
    .reset_index(name="N_genes")
)
biotype.to_csv(OUT_BIOTYPE, sep="\t", index=False)

# Compact sensitivity set table.
sens_rows = []
for k in [10, 20, 30, 40, 50, 60, 64]:
    z = a[a["N_independent_GSE"] >= k]
    sens_rows.append({
        "Minimum_independent_GSE": k,
        "N_stable_core_genes": len(z),
        "N_protein_coding": int((z["gene_type"] == "protein_coding").sum()),
        "N_non_protein_coding": int((z["gene_type"] != "protein_coding").sum()),
        "N_upregulated": int((z["pooled_log2FC"] > 0).sum()),
        "N_downregulated": int((z["pooled_log2FC"] < 0).sum()),
        "Median_abs_pooled_log2FC": float(z["pooled_log2FC"].abs().median()) if len(z) else np.nan,
        "Median_I2_total_percent": float(z["I2_total_percent"].median()) if len(z) else np.nan,
    })

sens = pd.DataFrame(sens_rows)
sens.to_csv(OUT_SENS, sep="\t", index=False)

n_pc = int((a["gene_type"] == "protein_coding").sum())
n_nonpc = len(a) - n_pc
n_cooks_any = int((a["N_Cooks_flagged_contrasts"] > 0).sum())

summary = [
    "===== HFD STABLE-CORE COVERAGE / BIOTYPE AUDIT v1 =====",
    "",
    "STATUS=PASS",
    f"Primary_stable_core_genes={len(a)}",
    f"Protein_coding={n_pc}",
    f"Non_protein_coding={n_nonpc}",
    f"Stable_core_with_any_Cooks_flagged_contrast_before_primaryQC={n_cooks_any}",
    "",
    "COVERAGE_SENSITIVITY_DOES_NOT_REDEFINE_PRIMARY_STABLE_CORE=YES",
]

for r in sens.itertuples(index=False):
    summary.append(
        f"Stable_core_ge_{r.Minimum_independent_GSE}_GSE="
        f"{r.N_stable_core_genes}; protein_coding={r.N_protein_coding}"
    )

summary += [
    "",
    f"Median_GSE_coverage={a['N_independent_GSE'].median():g}/64",
    f"Min_GSE_coverage={int(a['N_independent_GSE'].min())}/64",
    f"Max_GSE_coverage={int(a['N_independent_GSE'].max())}/64",
    f"Median_baseMean_across_stable_core_contrast_profiles={a['Median_baseMean_across_contrasts'].median():.6g}",
    "",
    f"Gene_level_audit={OUT_GENE.relative_to(ROOT)}",
    f"Biotype_by_coverage={OUT_BIOTYPE.relative_to(ROOT)}",
    f"Coverage_sensitivity_sets={OUT_SENS.relative_to(ROOT)}",
    "",
    "PRIMARY_STABLE_CORE_DEFINITION_CHANGED=NO",
    "META_ANALYSIS_RESULTS_CHANGED=NO",
    "PATHWAY_ANALYSIS_PERFORMED=NO",
]

OUT_SUMMARY.write_text("\n".join(summary) + "\n", encoding="utf-8")
print("\n".join(summary))

print("\n===== PROTEIN-CODING STABLE CORE =====")
cols = [
    "gene_symbol","Geneid","pooled_log2FC","pooled_FDR_BH",
    "N_independent_GSE","I2_total_percent",
    "prediction95_lower_new_GSE","prediction95_upper_new_GSE",
    "GSE_balanced_direction_consistency",
    "Median_baseMean_across_contrasts",
    "N_Cooks_flagged_contrasts",
]
print(
    a[a["gene_type"] == "protein_coding"][cols]
    .sort_values(["N_independent_GSE","pooled_log2FC"], ascending=[False,False])
    .to_string(index=False)
)

