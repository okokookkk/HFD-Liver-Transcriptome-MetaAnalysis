#!/usr/bin/env python3
from pathlib import Path
import argparse
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument(
    "--root",
    default=__import__("os").environ.get("HFD_PROJECT_ROOT", ".")
)
args = ap.parse_args()

ROOT = Path(args.root)
INFILE = (
    ROOT /
    "00_metadata/contrasts/moderators/manual_adjudication/"
    "HFD_primary105_moderator_manual_adjudication_batch02_v1.tsv"
)
OUTDIR = ROOT / "00_metadata/contrasts/moderators/manual_adjudication"
OUTDIR.mkdir(parents=True, exist_ok=True)

if not INFILE.exists():
    raise FileNotFoundError(INFILE)

x = pd.read_csv(INFILE, sep="\t", dtype=str).fillna("")

assert x["Contrast_id"].nunique() == 105
assert x["GSE"].nunique() == 64

variables = [
    "hfd_duration_weeks",
    "sex",
    "strain",
    "hfd_fat_kcal_percent",
    "control_diet_type",
]

# ------------------------------------------------------------
# Build effective "best available" value with strict precedence:
# manual final > GEO-auto resolved > UNKNOWN.
#
# For analysis readiness:
# - unresolved/missing -> not ready
# - sex=MIXED is known but not usable in binary male/female subgroup
# - suspicious auto HFD fat <20 or >80 is not analysis-ready until reviewed
# - typo-like auto strain strings are not analysis-ready until reviewed
# ------------------------------------------------------------

def auto_ready(var, value, status):
    if status != "RESOLVED":
        return False, "AUTO_NOT_RESOLVED"

    if var == "hfd_fat_kcal_percent":
        try:
            fv = float(value)
        except Exception:
            return False, "AUTO_NONNUMERIC"
        if fv < 20 or fv > 80:
            return False, "AUTO_SUSPICIOUS_FAT_PERCENT"

    if var == "strain" and value in {"C56BL/6J", "C56B6/J"}:
        return False, "AUTO_TYPO_LIKE_STRAIN"

    return True, "AUTO_RESOLVED"

rows = []

for _, r in x.iterrows():
    out = {
        "GSE": r["GSE"],
        "Contrast_id": r["Contrast_id"],
    }

    for var in variables:
        manual_final = r.get(f"{var}_final", "")
        manual_status = r.get(f"{var}_adjudication_status", "")
        auto_value = r.get(f"{var}_geo_auto_value", "")
        auto_status = r.get(f"{var}_geo_auto_status", "")

        if manual_final:
            value = manual_final
            source_level = "MANUAL_FINAL"
            known = True
            ready = True
            reason = manual_status or "MANUAL_FINAL"

            if var == "sex" and value == "MIXED":
                # Known metadata, but not usable for binary sex subgroup.
                ready = False
                reason = "KNOWN_MIXED_SEX"

        else:
            ok, why = auto_ready(var, auto_value, auto_status)
            if auto_status == "RESOLVED" and auto_value:
                value = auto_value
                source_level = "GEO_AUTO"
                known = True
                ready = ok
                reason = why
            else:
                value = "UNKNOWN"
                source_level = "UNRESOLVED"
                known = False
                ready = False
                reason = auto_status or "MISSING"

        out[f"{var}_effective"] = value
        out[f"{var}_source_level"] = source_level
        out[f"{var}_known"] = int(known)
        out[f"{var}_analysis_ready"] = int(ready)
        out[f"{var}_readiness_reason"] = reason

    rows.append(out)

e = pd.DataFrame(rows)

# ------------------------------------------------------------
# Coverage summary
# ------------------------------------------------------------
coverage_rows = []

for var in variables:
    known = int(e[f"{var}_known"].sum())
    ready = int(e[f"{var}_analysis_ready"].sum())

    manual_n = int(
        (e[f"{var}_source_level"] == "MANUAL_FINAL").sum()
    )
    auto_n = int(
        (e[f"{var}_source_level"] == "GEO_AUTO").sum()
    )
    unresolved_n = int(
        (e[f"{var}_source_level"] == "UNRESOLVED").sum()
    )

    # independent GSE with at least one ready contrast
    ready_gse = e.loc[
        e[f"{var}_analysis_ready"].eq(1), "GSE"
    ].nunique()

    coverage_rows.append({
        "Moderator": var,
        "Total_contrasts": len(e),
        "Known_contrasts": known,
        "Known_percent": 100.0 * known / len(e),
        "Analysis_ready_contrasts": ready,
        "Analysis_ready_percent": 100.0 * ready / len(e),
        "Independent_GSE_with_ready_value": ready_gse,
        "Manual_final_contrasts": manual_n,
        "GEO_auto_only_contrasts": auto_n,
        "Unresolved_contrasts": unresolved_n,
    })

coverage = pd.DataFrame(coverage_rows)

# ------------------------------------------------------------
# Distribution by contrast and independent GSE
# ------------------------------------------------------------
distribution_rows = []

for var in variables:
    z = e[e[f"{var}_known"].eq(1)].copy()

    for value, q in z.groupby(f"{var}_effective", sort=True):
        distribution_rows.append({
            "Moderator": var,
            "Value": value,
            "Contrasts": q["Contrast_id"].nunique(),
            "Independent_GSE": q["GSE"].nunique(),
            "Analysis_ready_contrasts": int(
                q[f"{var}_analysis_ready"].sum()
            ),
        })

distribution = pd.DataFrame(distribution_rows)

# ------------------------------------------------------------
# Remaining manual curation queue by GSE
# Use analysis-readiness rather than raw GEO missingness.
# ------------------------------------------------------------
weights = {
    "hfd_duration_weeks": 5,
    "sex": 4,
    "strain": 3,
    "hfd_fat_kcal_percent": 2,
    "control_diet_type": 1,
}

gse_rows = []

for gse, z in e.groupby("GSE", sort=True):
    # Skip fully manually adjudicated GSE only if all priority fields are ready.
    row = {
        "GSE": gse,
        "n_contrasts": z["Contrast_id"].nunique(),
    }
    score = 0

    for var in variables:
        not_ready_n = int(
            (z[f"{var}_analysis_ready"] == 0).sum()
        )
        row[f"{var}_not_ready_n"] = not_ready_n
        score += weights[var] * not_ready_n

    row["manual_priority_score"] = score

    # how much is already manually adjudicated at the GSE level?
    original = x[x["GSE"].eq(gse)]
    row["already_manual_adjudicated"] = int(
        original["overall_adjudication_status"]
        .str.contains("ADJUDICATED", na=False)
        .all()
    )

    gse_rows.append(row)

queue = pd.DataFrame(gse_rows)
queue = queue[queue["manual_priority_score"] > 0].copy()

queue = queue.sort_values(
    [
        "already_manual_adjudicated",
        "manual_priority_score",
        "hfd_duration_weeks_not_ready_n",
        "sex_not_ready_n",
        "strain_not_ready_n",
        "n_contrasts",
        "GSE",
    ],
    ascending=[True, False, False, False, False, False, True]
)

# ------------------------------------------------------------
# Manual progress summary
# ------------------------------------------------------------
manual_mask = x["overall_adjudication_status"].str.contains(
    "ADJUDICATED", na=False
)

progress = pd.DataFrame([{
    "Total_GSE": x["GSE"].nunique(),
    "Total_contrasts": x["Contrast_id"].nunique(),
    "Manually_adjudicated_GSE": x.loc[manual_mask, "GSE"].nunique(),
    "Manually_adjudicated_contrasts": x.loc[
        manual_mask, "Contrast_id"
    ].nunique(),
    "Remaining_GSE_without_complete_batch_adjudication":
        x["GSE"].nunique() - x.loc[manual_mask, "GSE"].nunique(),
    "Remaining_contrasts_without_complete_batch_adjudication":
        x["Contrast_id"].nunique() - x.loc[
            manual_mask, "Contrast_id"
        ].nunique(),
}])

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
effective_out = OUTDIR / "HFD_primary105_moderator_effective_after_batch02_v1.tsv"
coverage_out = OUTDIR / "HFD_primary105_moderator_coverage_after_batch02_v1.tsv"
distribution_out = OUTDIR / "HFD_primary105_moderator_distribution_after_batch02_v1.tsv"
queue_out = OUTDIR / "HFD_primary64_remaining_manual_curation_queue_after_batch02_v1.tsv"
progress_out = OUTDIR / "HFD_moderator_manual_curation_progress_after_batch02_v1.tsv"

e.to_csv(effective_out, sep="\t", index=False)
coverage.to_csv(coverage_out, sep="\t", index=False)
distribution.to_csv(distribution_out, sep="\t", index=False)
queue.to_csv(queue_out, sep="\t", index=False)
progress.to_csv(progress_out, sep="\t", index=False)

print("===== HFD MODERATOR COVERAGE AFTER BATCH02 =====")
print(progress.to_string(index=False))

print("\n===== EFFECTIVE COVERAGE =====")
print(
    coverage[
        [
            "Moderator",
            "Known_contrasts",
            "Known_percent",
            "Analysis_ready_contrasts",
            "Analysis_ready_percent",
            "Independent_GSE_with_ready_value",
            "Manual_final_contrasts",
            "GEO_auto_only_contrasts",
            "Unresolved_contrasts",
        ]
    ].to_string(index=False)
)

print("\n===== KEY DISTRIBUTIONS =====")
for var in variables:
    print(f"\n--- {var} ---")
    z = distribution[distribution["Moderator"].eq(var)]
    print(
        z.sort_values(
            ["Independent_GSE", "Contrasts"],
            ascending=False
        ).to_string(index=False)
    )

print("\n===== TOP 20 REMAINING GSE FOR MANUAL CURATION =====")
cols = [
    "GSE",
    "n_contrasts",
    "manual_priority_score",
    "hfd_duration_weeks_not_ready_n",
    "sex_not_ready_n",
    "strain_not_ready_n",
    "hfd_fat_kcal_percent_not_ready_n",
    "control_diet_type_not_ready_n",
    "already_manual_adjudicated",
]
print(queue[cols].head(20).to_string(index=False))

print("\nSaved:")
for p in [
    effective_out,
    coverage_out,
    distribution_out,
    queue_out,
    progress_out,
]:
    print(p.relative_to(ROOT))

print("\nInterpretation rule:")
print("- Known != necessarily analysis-ready.")
print("- Mixed sex is known metadata but not ready for binary male/female subgroup.")
print("- Suspicious auto HFD fat values and typo-like auto strain strings remain")
print("  non-ready until manually adjudicated.")
