#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()
INFILE = ROOT / "00_metadata/analysis_inputs/frozen_v1/HFD_primary105_contrast_analysis_input_audit_v1_1.tsv"

OUTDIR = ROOT / "00_metadata/analysis_inputs/preanalysis_design_v1"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_LEVELS = OUTDIR / "HFD_primary105_moderator_level_GSE_coverage_v1.tsv"
OUT_JOINT = OUTDIR / "HFD_primary105_moderator_joint_coverage_v1.tsv"
OUT_WITHIN = OUTDIR / "HFD_primary64_within_GSE_moderator_variation_v1.tsv"
OUT_SUMMARY = OUTDIR / "HFD_primary105_moderator_modelability_summary_v1.txt"

if not INFILE.exists():
    raise RuntimeError(f"Missing input: {INFILE}")

x = pd.read_csv(INFILE, sep="\t", dtype=str).fillna("")

mods = [
    ("hfd_duration_weeks", "continuous"),
    ("sex", "categorical"),
    ("strain", "categorical"),
    ("hfd_fat_kcal_percent", "continuous"),
    ("control_diet_type", "categorical"),
]

unknown_tokens = {"", "UNKNOWN", "NOT_REPORTED", "UNRESOLVED", "NA", "N/A"}

def usable(s):
    return ~s.astype(str).str.strip().str.upper().isin(unknown_tokens)

level_rows = []
summary_lines = []

for mod, kind in mods:
    if mod not in x.columns:
        raise RuntimeError(f"Missing analysis alias: {mod}")

    mask = usable(x[mod])
    z = x.loc[mask, ["Contrast_id", "GSE", mod]].copy()
    z[mod] = z[mod].astype(str).str.strip()

    summary_lines.append(
        f"{mod}: usable_contrasts={len(z)}/105; independent_GSE={z['GSE'].nunique()}/64"
    )

    if kind == "continuous":
        vals = pd.to_numeric(z[mod], errors="coerce")
        bad = z.loc[vals.isna(), mod].unique().tolist()
        if bad:
            raise RuntimeError(f"{mod}: nonnumeric usable values: {bad}")
        z["_numeric"] = vals

        # Exact values are retained because duration/fat can vary within study.
        for val, q in z.groupby("_numeric", sort=True):
            level_rows.append({
                "Moderator": mod,
                "Type": kind,
                "Level_or_value": f"{val:g}",
                "N_contrasts": len(q),
                "N_independent_GSE": q["GSE"].nunique(),
            })

        summary_lines.append(
            f"  range={vals.min():g}..{vals.max():g}; "
            f"median={vals.median():g}; unique_values={vals.nunique()}"
        )

    else:
        for val, q in z.groupby(mod, sort=True):
            level_rows.append({
                "Moderator": mod,
                "Type": kind,
                "Level_or_value": val,
                "N_contrasts": len(q),
                "N_independent_GSE": q["GSE"].nunique(),
            })

        lv = (
            z.groupby(mod)
            .agg(N_contrasts=("Contrast_id", "size"),
                 N_GSE=("GSE", "nunique"))
            .sort_values(["N_GSE", "N_contrasts"], ascending=False)
        )
        summary_lines.append(
            "  levels=" + "; ".join(
                f"{idx}[C={int(r.N_contrasts)},GSE={int(r.N_GSE)}]"
                for idx, r in lv.iterrows()
            )
        )

levels = pd.DataFrame(level_rows)
levels.to_csv(OUT_LEVELS, sep="\t", index=False)

# Joint complete-case coverage for all subsets of practical interest.
joint_specs = [
    ["hfd_duration_weeks"],
    ["sex"],
    ["strain"],
    ["hfd_fat_kcal_percent"],
    ["control_diet_type"],
    ["hfd_duration_weeks", "hfd_fat_kcal_percent"],
    ["sex", "strain"],
    ["sex", "strain", "control_diet_type"],
    ["hfd_duration_weeks", "hfd_fat_kcal_percent", "sex"],
    ["hfd_duration_weeks", "hfd_fat_kcal_percent", "sex", "strain", "control_diet_type"],
]

joint_rows = []
for spec in joint_specs:
    mask = pd.Series(True, index=x.index)
    for mod in spec:
        mask &= usable(x[mod])
    q = x.loc[mask]
    joint_rows.append({
        "Moderators": "+".join(spec),
        "N_contrasts": len(q),
        "N_independent_GSE": q["GSE"].nunique(),
    })

joint = pd.DataFrame(joint_rows)
joint.to_csv(OUT_JOINT, sep="\t", index=False)

# Within-GSE moderator variation: useful because multi-contrast GSEs may
# contribute within-study moderator contrasts.
within_rows = []
for gse, q in x.groupby("GSE", sort=True):
    row = {"GSE": gse, "N_contrasts": len(q)}
    for mod, _ in mods:
        vals = q.loc[usable(q[mod]), mod].astype(str).str.strip()
        row[f"{mod}_N_unique_usable"] = vals.nunique()
        row[f"{mod}_varies_within_GSE"] = "YES" if vals.nunique() > 1 else "NO"
    within_rows.append(row)

within = pd.DataFrame(within_rows)
within.to_csv(OUT_WITHIN, sep="\t", index=False)

summary = [
    "===== HFD PRIMARY105 MODERATOR MODELABILITY AUDIT v1 =====",
    "",
    "STATUS=PASS",
    *summary_lines,
    "",
    "JOINT_COVERAGE",
]

for r in joint.itertuples(index=False):
    summary.append(
        f"{r.Moderators}: contrasts={r.N_contrasts}/105; "
        f"independent_GSE={r.N_independent_GSE}/64"
    )

summary += [
    "",
    "WITHIN_GSE_VARIATION",
]

for mod, _ in mods:
    n = int((within[f"{mod}_varies_within_GSE"] == "YES").sum())
    summary.append(f"{mod}: varies_within_GSE={n}/64")

summary += [
    "",
    f"Level_coverage={OUT_LEVELS.relative_to(ROOT)}",
    f"Joint_coverage={OUT_JOINT.relative_to(ROOT)}",
    f"Within_GSE_variation={OUT_WITHIN.relative_to(ROOT)}",
    "",
    "NO_METADATA_MODIFIED=YES",
    "NO_MODEL_FIT=YES",
]

OUT_SUMMARY.write_text("\n".join(summary) + "\n", encoding="utf-8")
print("\n".join(summary))

