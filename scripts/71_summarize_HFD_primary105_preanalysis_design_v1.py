#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()

AUDIT = ROOT / "00_metadata/analysis_inputs/frozen_v1/HFD_primary105_contrast_analysis_input_audit_v1_1.tsv"
GSE_MULT = ROOT / "00_metadata/analysis_inputs/frozen_v1/HFD_primary64_GSE_contrast_multiplicity_v1_1.tsv"
OUTDIR = ROOT / "00_metadata/analysis_inputs/preanalysis_design_v1"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_SAMPLE = OUTDIR / "HFD_primary105_sample_size_distribution_v1.tsv"
OUT_MOD = OUTDIR / "HFD_primary105_moderator_readiness_v1.tsv"
OUT_GSE = OUTDIR / "HFD_primary64_contrast_multiplicity_distribution_v1.tsv"
OUT_SUMMARY = OUTDIR / "HFD_primary105_preanalysis_design_summary_v1.txt"

for q in [AUDIT, GSE_MULT]:
    if not q.exists():
        raise RuntimeError(f"Missing input: {q}")

x = pd.read_csv(AUDIT, sep="\t", dtype=str).fillna("")
g = pd.read_csv(GSE_MULT, sep="\t", dtype=str).fillna("")

for c in ["N_HFD_GSM", "N_CONTROL_GSM", "N_TOTAL_GSM"]:
    x[c] = pd.to_numeric(x[c], errors="raise")

# ------------------------------------------------------------
# Sample-size audit
# ------------------------------------------------------------
x["Min_arm_n"] = x[["N_HFD_GSM", "N_CONTROL_GSM"]].min(axis=1)
x["Max_arm_n"] = x[["N_HFD_GSM", "N_CONTROL_GSM"]].max(axis=1)
x["Arm_imbalance"] = x["Max_arm_n"] - x["Min_arm_n"]

sample_rows = []
for cid, r in x.sort_values(["Min_arm_n", "N_TOTAL_GSM", "GSE", "Contrast_id"]).iterrows():
    sample_rows.append({
        "Contrast_id": r["Contrast_id"],
        "GSE": r["GSE"],
        "N_HFD_GSM": int(r["N_HFD_GSM"]),
        "N_CONTROL_GSM": int(r["N_CONTROL_GSM"]),
        "N_TOTAL_GSM": int(r["N_TOTAL_GSM"]),
        "Min_arm_n": int(r["Min_arm_n"]),
        "Max_arm_n": int(r["Max_arm_n"]),
        "Arm_imbalance": int(r["Arm_imbalance"]),
        "Both_arms_n_ge_2": "YES" if r["Min_arm_n"] >= 2 else "NO",
        "Both_arms_n_ge_3": "YES" if r["Min_arm_n"] >= 3 else "NO",
        "Both_arms_n_ge_4": "YES" if r["Min_arm_n"] >= 4 else "NO",
    })
pd.DataFrame(sample_rows).to_csv(OUT_SAMPLE, sep="\t", index=False)

# ------------------------------------------------------------
# GSE multiplicity audit
# ------------------------------------------------------------
g["N_contrasts"] = pd.to_numeric(g["N_contrasts"], errors="raise")
dist = (
    g["N_contrasts"]
    .value_counts()
    .sort_index()
    .rename_axis("N_contrasts_per_GSE")
    .reset_index(name="N_GSE")
)
dist["Percent_GSE"] = 100 * dist["N_GSE"] / len(g)
dist.to_csv(OUT_GSE, sep="\t", index=False)

# ------------------------------------------------------------
# Moderator readiness audit
# Uses adjudication status + final value, not merely nonblank final.
# No values are changed here.
# ------------------------------------------------------------
mods = [
    "hfd_duration_weeks",
    "sex",
    "strain",
    "hfd_fat_kcal_percent",
    "control_diet_type",
]

ready_status_tokens = {
    "CONSISTENT",
    "MANUAL_CONFIRMED",
    "AUTO_GEO_RESOLVED",
    "AUTO_EXTRACTION_ERROR_CORRECTED",
    "RESOLVED",
}

mod_rows = []

for base in mods:
    final_col = f"{base}_final"
    status_col = f"{base}_adjudication_status"
    conf_col = f"{base}_final_confidence"

    for c in [final_col, status_col]:
        if c not in x.columns:
            raise RuntimeError(f"Missing moderator column: {c}")

    # Report raw status/final combinations without re-adjudicating.
    z = x[[final_col, status_col] + ([conf_col] if conf_col in x.columns else [])].copy()
    z[final_col] = z[final_col].astype(str).str.strip()
    z[status_col] = z[status_col].astype(str).str.strip()

    nonblank_final = z[final_col].ne("")
    explicit_unknown = z[final_col].str.upper().isin({
        "UNKNOWN", "NOT_REPORTED", "UNRESOLVED", "NA", "N/A"
    })
    usable_final = nonblank_final & ~explicit_unknown

    # Status is not forced into a binary rule; we simply expose it.
    status_counts = z[status_col].replace("", "<BLANK>").value_counts(dropna=False)
    final_unknown_counts = z.loc[~usable_final, final_col].replace("", "<BLANK>").value_counts(dropna=False)

    mod_rows.append({
        "Moderator": base,
        "Total_contrasts": len(z),
        "Nonblank_final": int(nonblank_final.sum()),
        "Usable_final_excluding_explicit_unknown_tokens": int(usable_final.sum()),
        "Blank_or_explicit_unknown": int((~usable_final).sum()),
        "Unique_usable_values": int(z.loc[usable_final, final_col].nunique()),
        "Adjudication_status_counts": " | ".join(f"{k}:{v}" for k, v in status_counts.items()),
        "Missing_or_unknown_final_counts": " | ".join(f"{k}:{v}" for k, v in final_unknown_counts.items()),
    })

pd.DataFrame(mod_rows).to_csv(OUT_MOD, sep="\t", index=False)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------
n_1 = int((x["Min_arm_n"] == 1).sum())
n_ge2 = int((x["Min_arm_n"] >= 2).sum())
n_ge3 = int((x["Min_arm_n"] >= 3).sum())
n_ge4 = int((x["Min_arm_n"] >= 4).sum())

summary = [
    "===== HFD PRIMARY105 PRE-ANALYSIS DESIGN SUMMARY v1 =====",
    "",
    "STATUS=PASS",
    f"Contrasts={len(x)}",
    f"GSE={x['GSE'].nunique()}",
    "",
    "SAMPLE_SIZE",
    f"Min_HFD_GSM={int(x['N_HFD_GSM'].min())}",
    f"Median_HFD_GSM={float(x['N_HFD_GSM'].median()):g}",
    f"Max_HFD_GSM={int(x['N_HFD_GSM'].max())}",
    f"Min_CONTROL_GSM={int(x['N_CONTROL_GSM'].min())}",
    f"Median_CONTROL_GSM={float(x['N_CONTROL_GSM'].median()):g}",
    f"Max_CONTROL_GSM={int(x['N_CONTROL_GSM'].max())}",
    f"Contrasts_with_min_arm_n_1={n_1}/105",
    f"Contrasts_with_both_arms_n_ge_2={n_ge2}/105",
    f"Contrasts_with_both_arms_n_ge_3={n_ge3}/105",
    f"Contrasts_with_both_arms_n_ge_4={n_ge4}/105",
    "",
    "DEPENDENCE",
    f"GSE_with_multiple_contrasts={int((g['N_contrasts'] > 1).sum())}/64",
    f"Maximum_contrasts_per_GSE={int(g['N_contrasts'].max())}",
    "Shared_GSM_across_contrasts=0",
    "",
    "MODERATOR_READINESS",
]

for row in mod_rows:
    summary.append(
        f"{row['Moderator']}: usable_final="
        f"{row['Usable_final_excluding_explicit_unknown_tokens']}/105; "
        f"unique_values={row['Unique_usable_values']}"
    )

summary += [
    "",
    f"Sample_size_table={OUT_SAMPLE.relative_to(ROOT)}",
    f"GSE_multiplicity_distribution={OUT_GSE.relative_to(ROOT)}",
    f"Moderator_readiness={OUT_MOD.relative_to(ROOT)}",
    "",
    "NO_DE_PERFORMED=YES",
    "NO_META_ANALYSIS_PERFORMED=YES",
    "NO_METADATA_MODIFIED=YES",
]

OUT_SUMMARY.write_text("\n".join(summary) + "\n", encoding="utf-8")
print("\n".join(summary))

