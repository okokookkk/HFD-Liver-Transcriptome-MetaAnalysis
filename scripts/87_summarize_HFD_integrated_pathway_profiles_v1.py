#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()

MAP = ROOT / "06_pathway/heterogeneity_v1/HFD_Hallmark_mean_duration_heterogeneity_map_v1.tsv"

OUTDIR = ROOT / "00_metadata/analysis_outputs/integrated_pathway_v1"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_CLASS = OUTDIR / "HFD_Hallmark_integrated_profile_classes_v1.tsv"
OUT_SUM = OUTDIR / "HFD_Hallmark_integrated_profile_summary_v1.txt"

if not MAP.exists():
    raise RuntimeError(f"Missing integrated map: {MAP}")

x = pd.read_csv(MAP, sep="\t")

need = [
    "pathway",
    "mean_response_NES", "mean_response_FDR",
    "duration_modifier_NES", "duration_modifier_FDR",
    "heterogeneity_Q_NES", "heterogeneity_Q_FDR",
    "heterogeneity_I2_NES", "heterogeneity_I2_FDR",
    "heterogeneity_Q_ge30GSE_NES", "heterogeneity_Q_ge30GSE_FDR",
]
missing = [c for c in need if c not in x.columns]
if missing:
    raise RuntimeError(f"Missing columns: {missing}")

for c in need[1:]:
    x[c] = pd.to_numeric(x[c], errors="coerce")

x["Mean_FDR05"] = x["mean_response_FDR"] < 0.05
x["Duration_FDR05"] = x["duration_modifier_FDR"] < 0.05
x["Heterogeneity_Q_FDR05"] = x["heterogeneity_Q_FDR"] < 0.05
x["Heterogeneity_I2_FDR05"] = x["heterogeneity_I2_FDR"] < 0.05
x["Heterogeneity_Q_ge30_FDR05"] = x["heterogeneity_Q_ge30GSE_FDR"] < 0.05

# Descriptive robustness for high heterogeneity.
x["Robust_high_heterogeneity"] = (
    x["Heterogeneity_Q_FDR05"]
    & (x["heterogeneity_Q_NES"] > 0)
    & (x["heterogeneity_I2_NES"] > 0)
    & (x["heterogeneity_Q_ge30GSE_NES"] > 0)
    & (x["Heterogeneity_I2_FDR05"] | x["Heterogeneity_Q_ge30_FDR05"])
)

def sign_label(v):
    if pd.isna(v): return "NA"
    if v > 0: return "POS"
    if v < 0: return "NEG"
    return "ZERO"

x["Mean_direction"] = x["mean_response_NES"].map(sign_label)
x["Duration_direction"] = x["duration_modifier_NES"].map(sign_label)

def classify(r):
    m = r["Mean_FDR05"]
    d = r["Duration_FDR05"]
    h = r["Robust_high_heterogeneity"]
    if m and d and h:
        return "MEAN+HETEROGENEITY+DURATION"
    if h and d and not m:
        return "HETEROGENEITY+DURATION_ONLY"
    if m and h and not d:
        return "MEAN+HETEROGENEITY_ONLY"
    if m and d and not h:
        return "MEAN+DURATION_ONLY"
    if h and not m and not d:
        return "HETEROGENEITY_ONLY"
    if d and not m and not h:
        return "DURATION_ONLY"
    if m and not d and not h:
        return "MEAN_ONLY"
    return "NONE"

x["Integrated_class"] = x.apply(classify, axis=1)

# Additional descriptive relation: do duration-sensitive pathways tend to be heterogeneous?
x["abs_duration_NES"] = x["duration_modifier_NES"].abs()
x["abs_mean_NES"] = x["mean_response_NES"].abs()

pearson_absdur_het = x["abs_duration_NES"].corr(x["heterogeneity_Q_NES"], method="pearson")
spearman_absdur_het = x["abs_duration_NES"].corr(x["heterogeneity_Q_NES"], method="spearman")
pearson_absmean_het = x["abs_mean_NES"].corr(x["heterogeneity_Q_NES"], method="pearson")
spearman_absmean_het = x["abs_mean_NES"].corr(x["heterogeneity_Q_NES"], method="spearman")

# Sort by integrated interest: robust high heterogeneity first, then duration FDR, then heterogeneity FDR.
x["_rank_h"] = (~x["Robust_high_heterogeneity"]).astype(int)
x["_rank_d"] = x["duration_modifier_FDR"].fillna(1)
x["_rank_q"] = x["heterogeneity_Q_FDR"].fillna(1)
x = x.sort_values(["_rank_h", "_rank_d", "_rank_q", "pathway"]).drop(
    columns=["_rank_h","_rank_d","_rank_q"]
)
x.to_csv(OUT_CLASS, sep="\t", index=False)

counts = x["Integrated_class"].value_counts()

def pathway_list(df, n=20):
    if df.empty:
        return ["NONE"]
    out = []
    for r in df.head(n).itertuples(index=False):
        out.append(
            f"{r.pathway}"
            f"|meanNES={r.mean_response_NES:.3f}|meanFDR={r.mean_response_FDR:.3g}"
            f"|durNES={r.duration_modifier_NES:.3f}|durFDR={r.duration_modifier_FDR:.3g}"
            f"|hetNES={r.heterogeneity_Q_NES:.3f}|hetFDR={r.heterogeneity_Q_FDR:.3g}"
        )
    return out

triple = x[x["Integrated_class"] == "MEAN+HETEROGENEITY+DURATION"].copy()
het_dur = x[
    x["Robust_high_heterogeneity"] & x["Duration_FDR05"]
].sort_values(["duration_modifier_FDR","heterogeneity_Q_FDR"])

het_not_dur = x[
    x["Robust_high_heterogeneity"] & ~x["Duration_FDR05"]
].sort_values("heterogeneity_Q_FDR")

mean_not_het = x[
    x["Mean_FDR05"] & ~x["Robust_high_heterogeneity"]
].sort_values("mean_response_FDR")

summary = [
    "===== HFD INTEGRATED HALLMARK PROFILE SUMMARY v1 =====",
    "",
    "STATUS=PASS",
    f"Hallmark_pathways={len(x)}",
    f"Mean_response_FDR05={int(x['Mean_FDR05'].sum())}",
    f"Duration_modifier_FDR05={int(x['Duration_FDR05'].sum())}",
    f"Robust_high_heterogeneity={int(x['Robust_high_heterogeneity'].sum())}",
    "",
    "INTEGRATED_CLASS_COUNTS",
]
for k, v in counts.items():
    summary.append(f"{k}={v}")

summary += [
    "",
    f"Pearson_absDurationNES_vs_heterogeneityNES={pearson_absdur_het:.6g}",
    f"Spearman_absDurationNES_vs_heterogeneityNES={spearman_absdur_het:.6g}",
    f"Pearson_absMeanNES_vs_heterogeneityNES={pearson_absmean_het:.6g}",
    f"Spearman_absMeanNES_vs_heterogeneityNES={spearman_absmean_het:.6g}",
    "",
    "MEAN_PLUS_HETEROGENEITY_PLUS_DURATION",
    *pathway_list(triple, 20),
    "",
    "ROBUST_HETEROGENEITY_PLUS_DURATION",
    *pathway_list(het_dur, 20),
    "",
    "ROBUST_HETEROGENEITY_WITHOUT_DURATION_FDR",
    *pathway_list(het_not_dur, 20),
    "",
    "MEAN_RESPONSE_FDR_WITHOUT_ROBUST_HETEROGENEITY",
    *pathway_list(mean_not_het, 20),
    "",
    f"Integrated_table={OUT_CLASS.relative_to(ROOT)}",
    "",
    "ANALYSIS_STATUS=DESCRIPTIVE_EXPLORATORY",
    "NO_CAUSAL_INTERPRETATION=YES",
]
OUT_SUM.write_text("\n".join(summary) + "\n", encoding="utf-8")
print("\n".join(summary))

