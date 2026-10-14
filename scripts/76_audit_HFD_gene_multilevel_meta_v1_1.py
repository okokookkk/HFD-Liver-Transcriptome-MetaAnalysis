#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()

META = ROOT / "05_meta/gene_level_v1/HFD_gene_multilevel_meta_primary_v1_1.tsv.gz"
FAIL = ROOT / "00_metadata/analysis_outputs/gene_meta_v1_1/HFD_gene_multilevel_meta_failures_v1_1.tsv"

OUTDIR = ROOT / "00_metadata/analysis_outputs/gene_meta_v1_1"
OUT_FAILURE_CAT = OUTDIR / "HFD_gene_multilevel_meta_failure_categories_v1_1.tsv"
OUT_BOUNDARY = OUTDIR / "HFD_gene_multilevel_meta_variance_boundary_audit_v1_1.tsv"
OUT_SUMMARY = OUTDIR / "HFD_gene_multilevel_meta_model_audit_summary_v1_1.txt"

for q in [META, FAIL]:
    if not q.exists():
        raise RuntimeError(f"Missing required input: {q}")

m = pd.read_csv(META, sep="\t", compression="gzip")
f = pd.read_csv(FAIL, sep="\t", dtype=str).fillna("")

required = [
    "Geneid",
    "N_contrasts",
    "N_independent_GSE",
    "tau2_GSE",
    "tau2_contrast_within_GSE",
    "tau2_total",
    "I2_total_percent",
    "convergence",
]
missing = [c for c in required if c not in m.columns]
if missing:
    raise RuntimeError(f"Primary meta result missing columns: {missing}")

for c in [
    "N_contrasts", "N_independent_GSE",
    "tau2_GSE", "tau2_contrast_within_GSE",
    "tau2_total", "I2_total_percent"
]:
    m[c] = pd.to_numeric(m[c], errors="coerce")

def categorize_error(s):
    t = str(s).lower()
    if "optimizer" in t or "converg" in t:
        return "OPTIMIZER_OR_CONVERGENCE"
    if "positive definite" in t or "not positive" in t or "hessian" in t:
        return "HESSIAN_OR_POSDEF"
    if "singular" in t:
        return "SINGULAR"
    if "sigma" in t or "variance" in t or "tau" in t:
        return "VARIANCE_COMPONENT"
    if "iteration" in t or "maxit" in t:
        return "ITERATION_LIMIT"
    if "non-finite" in t or "nan" in t or "inf" in t:
        return "NONFINITE_NUMERIC"
    if not t.strip():
        return "EMPTY_ERROR"
    return "OTHER"

if len(f):
    f["Failure_category"] = f["Error"].map(categorize_error)
    cat = (
        f.groupby("Failure_category", dropna=False)
        .size()
        .reset_index(name="N_genes")
        .sort_values(["N_genes", "Failure_category"], ascending=[False, True])
    )
else:
    cat = pd.DataFrame(columns=["Failure_category", "N_genes"])

cat.to_csv(OUT_FAILURE_CAT, sep="\t", index=False)

eps = 1e-10

m["tau2_GSE_at_zero"] = m["tau2_GSE"].abs().le(eps)
m["tau2_contrast_at_zero"] = m["tau2_contrast_within_GSE"].abs().le(eps)
m["both_tau2_zero"] = m["tau2_GSE_at_zero"] & m["tau2_contrast_at_zero"]
m["either_tau2_zero"] = m["tau2_GSE_at_zero"] | m["tau2_contrast_at_zero"]
m["both_tau2_positive"] = (m["tau2_GSE"] > eps) & (m["tau2_contrast_within_GSE"] > eps)

m["GSE_fraction_of_tau2_total"] = np.where(
    m["tau2_total"] > eps,
    m["tau2_GSE"] / m["tau2_total"],
    np.nan
)

boundary_cols = [
    "Geneid",
    "N_contrasts",
    "N_independent_GSE",
    "tau2_GSE",
    "tau2_contrast_within_GSE",
    "tau2_total",
    "GSE_fraction_of_tau2_total",
    "I2_total_percent",
    "tau2_GSE_at_zero",
    "tau2_contrast_at_zero",
    "both_tau2_zero",
    "either_tau2_zero",
    "both_tau2_positive",
    "convergence",
]
m[boundary_cols].to_csv(OUT_BOUNDARY, sep="\t", index=False)

n = len(m)
n_gse_zero = int(m["tau2_GSE_at_zero"].sum())
n_contrast_zero = int(m["tau2_contrast_at_zero"].sum())
n_both_zero = int(m["both_tau2_zero"].sum())
n_either_zero = int(m["either_tau2_zero"].sum())
n_both_pos = int(m["both_tau2_positive"].sum())
n_conv_review = int((m["convergence"].astype(str) != "PASS").sum())

gse_frac = m.loc[
    np.isfinite(m["GSE_fraction_of_tau2_total"]),
    "GSE_fraction_of_tau2_total"
]

summary = [
    "===== HFD GENE MULTILEVEL META MODEL AUDIT v1.1 =====",
    "",
    "STATUS=PASS",
    f"Successful_models={n}",
    f"Failed_models={len(f)}",
    f"Successful_models_with_convergence_review={n_conv_review}",
    "",
    "VARIANCE_COMPONENT_BOUNDARY",
    f"tau2_GSE_at_zero={n_gse_zero}/{n}",
    f"tau2_contrast_within_GSE_at_zero={n_contrast_zero}/{n}",
    f"either_variance_component_at_zero={n_either_zero}/{n}",
    f"both_variance_components_at_zero={n_both_zero}/{n}",
    f"both_variance_components_positive={n_both_pos}/{n}",
    "",
    f"Median_tau2_GSE={m['tau2_GSE'].median():.8g}",
    f"Median_tau2_contrast_within_GSE={m['tau2_contrast_within_GSE'].median():.8g}",
    f"Median_tau2_total={m['tau2_total'].median():.8g}",
    f"Median_I2_total_percent={m['I2_total_percent'].median():.6g}",
    f"Median_GSE_fraction_of_tau2_total={gse_frac.median():.6g}",
    "",
    "FAILURE_CATEGORIES",
]

if len(cat):
    for r in cat.itertuples(index=False):
        summary.append(f"{r.Failure_category}={r.N_genes}")
else:
    summary.append("NONE=0")

summary += [
    "",
    f"Failure_categories={OUT_FAILURE_CAT.relative_to(ROOT)}",
    f"Variance_boundary_audit={OUT_BOUNDARY.relative_to(ROOT)}",
    "",
    "NOTE=Variance-component decomposition is secondary; total heterogeneity is primary.",
    "NO_MODEL_REFIT_PERFORMED=YES",
]

OUT_SUMMARY.write_text("\n".join(summary) + "\n", encoding="utf-8")
print("\n".join(summary))

