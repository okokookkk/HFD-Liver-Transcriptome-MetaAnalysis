#!/usr/bin/env python3
from pathlib import Path
import argparse
import pandas as pd

# ============================================================
# HFD Meta-analysis
# Initialize manual/publication biological moderator adjudication
#
# INPUT:
#   GEO-auto contrast moderator table v1.1.1
#
# PRINCIPLES:
# - Never overwrite GEO-auto extraction.
# - Preserve raw/auto values verbatim.
# - Manual/publication values are entered in separate columns.
# - UNKNOWN/CONFLICT are acceptable final states.
# - Moderator uncertainty does NOT remove an otherwise valid
#   HFD-vs-control contrast from the primary meta-analysis.
# - If CONTROL/HFD sample identity itself is contradicted, flag
#   SAMPLE_ROLE_IDENTITY_REVIEW_REQUIRED; do not silently relabel.
# ============================================================

ap = argparse.ArgumentParser()
ap.add_argument(
    "--root",
    default=__import__("os").environ.get("HFD_PROJECT_ROOT", ".")
)
args = ap.parse_args()

ROOT = Path(args.root)

AUTO = (
    ROOT /
    "00_metadata/contrasts/moderators/"
    "HFD_primary105_contrast_moderators_GEOauto_v1_1_1.tsv"
)

OUTDIR = ROOT / "00_metadata/contrasts/moderators/manual_adjudication"
OUTDIR.mkdir(parents=True, exist_ok=True)

if not AUTO.exists():
    raise FileNotFoundError(AUTO)

x = pd.read_csv(AUTO, sep="\t", dtype=str).fillna("")

required = {
    "GSE", "Contrast_id",
    "n_CONTROL_GSM", "n_HFD_GSM",
    "sex_value", "sex_status",
    "strain_value", "strain_status",
    "hfd_duration_weeks_value", "hfd_duration_weeks_status",
    "hfd_fat_kcal_percent_value", "hfd_fat_kcal_percent_status",
    "control_diet_type_value", "control_diet_type_status",
}
missing = required - set(x.columns)
if missing:
    raise RuntimeError(
        f"Auto moderator table missing required columns: {sorted(missing)}"
    )

assert x["Contrast_id"].nunique() == 105
assert x["GSE"].nunique() == 64

# ------------------------------------------------------------
# Variables to adjudicate
# ------------------------------------------------------------

variables = [
    ("hfd_duration_weeks", 1),
    ("sex", 2),
    ("strain", 3),
    ("hfd_fat_kcal_percent", 4),
    ("control_diet_type", 5),
]

# Publication/manual source categories.
source_categories = [
    "GEO_GSM",
    "GEO_SERIES",
    "SRA_BIOSAMPLE",
    "PUBLICATION_METHODS",
    "PUBLICATION_RESULTS",
    "PUBLICATION_FIGURE_LEGEND",
    "SUPPLEMENTARY_TABLE",
    "SUPPLEMENTARY_METHODS",
    "DIET_MANUFACTURER_CODE",
    "OTHER",
]

adjudication_status_categories = [
    "NOT_REVIEWED",
    "CONSISTENT",
    "MANUAL_CONFIRMED",
    "GEO_TYPO_CORRECTED",
    "PAPER_GEO_CONFLICT_RESOLVED",
    "CONFLICT_UNRESOLVED",
    "NOT_REPORTED",
    "NOT_APPLICABLE",
    "QC_DISCORDANT",
    "SAMPLE_ROLE_IDENTITY_REVIEW_REQUIRED",
]

confidence_categories = [
    "HIGH",
    "MEDIUM",
    "LOW",
    "UNRESOLVED",
]

# ------------------------------------------------------------
# Build contrast-level adjudication template
# ------------------------------------------------------------

base_cols = [
    "GSE",
    "Contrast_id",
    "n_CONTROL_GSM",
    "n_HFD_GSM",
]

out = x[base_cols].copy()

for var, rank in variables:
    # Preserve auto evidence.
    for suffix in [
        "value",
        "status",
        "raw",
        "source",
        "confidence",
    ]:
        src = f"{var}_{suffix}"
        if src in x.columns:
            out[f"{var}_geo_auto_{suffix}"] = x[src]
        else:
            out[f"{var}_geo_auto_{suffix}"] = ""

    # Manual/publication adjudication fields.
    out[f"{var}_manual_value"] = ""
    out[f"{var}_manual_source_type"] = ""
    out[f"{var}_manual_source_citation"] = ""
    out[f"{var}_manual_evidence"] = ""
    out[f"{var}_manual_confidence"] = ""
    out[f"{var}_adjudication_status"] = "NOT_REVIEWED"

    # Final fields intentionally blank until adjudicated.
    out[f"{var}_final"] = ""
    out[f"{var}_final_confidence"] = ""
    out[f"{var}_final_source"] = ""
    out[f"{var}_final_note"] = ""

# Contrast/sample identity guardrail fields.
out["sample_role_identity_status"] = "NOT_REVIEWED"
out["sample_role_identity_evidence"] = ""
out["sample_role_identity_note"] = ""

# General curation provenance.
out["publication_pmid"] = ""
out["publication_doi"] = ""
out["publication_title"] = ""
out["supplement_reviewed"] = ""
out["curator"] = ""
out["curation_date"] = ""
out["overall_adjudication_status"] = "NOT_REVIEWED"
out["general_note"] = ""

# ------------------------------------------------------------
# Priority scoring
#
# Highest priority:
# 1. duration conflict/missing
# 2. sex missing
# 3. strain missing
# 4. fat %
# 5. control diet
#
# Score at GSE level later to avoid repeatedly opening the same paper.
# ------------------------------------------------------------

weight = {
    "hfd_duration_weeks": 5,
    "sex": 4,
    "strain": 3,
    "hfd_fat_kcal_percent": 2,
    "control_diet_type": 1,
}

def needs_review(status):
    return status != "RESOLVED"

score = []
issues = []

for _, r in x.iterrows():
    s = 0
    iss = []

    for var, _rank in variables:
        st = r.get(f"{var}_status", "")
        if needs_review(st):
            s += weight[var]
            iss.append(f"{var}:{st or 'UNKNOWN'}")

        # Even resolved auto values with suspicious diet-fat values
        # should be manually verified before analysis.
        if var == "hfd_fat_kcal_percent":
            val = r.get(f"{var}_value", "")
            try:
                fv = float(val)
            except Exception:
                fv = None
            if fv is not None and (fv < 20 or fv > 80):
                s += 3
                iss.append(f"hfd_fat_kcal_percent:SUSPICIOUS_{fv:g}")

    # Low-confidence duration conflicts are especially important.
    if r.get("hfd_duration_weeks_status", "") == "CONFLICT":
        s += 5
        iss.append("duration:CONFLICT_PRIORITY")

    score.append(s)
    issues.append(";".join(iss))

out.insert(4, "manual_priority_score", score)
out.insert(5, "manual_review_issues", issues)

# ------------------------------------------------------------
# GSE-level publication queue
# ------------------------------------------------------------

gse_rows = []

for gse, z in out.groupby("GSE", sort=True):
    row = {
        "GSE": gse,
        "n_contrasts": z["Contrast_id"].nunique(),
        "manual_priority_score": int(z["manual_priority_score"].sum()),
    }

    # Count contrast-level needs by variable using original auto table.
    zx = x[x["GSE"] == gse]

    for var, _rank in variables:
        row[f"{var}_needs_review_n"] = int(
            (zx[f"{var}_status"] != "RESOLVED").sum()
        )

    # Additional suspicious fat values.
    suspicious_fat = 0
    for v in zx["hfd_fat_kcal_percent_value"]:
        try:
            fv = float(v)
        except Exception:
            continue
        if fv < 20 or fv > 80:
            suspicious_fat += 1
    row["hfd_fat_suspicious_value_n"] = suspicious_fat

    # Collect contrast IDs for paper-level review.
    row["Contrast_ids"] = ",".join(
        sorted(z["Contrast_id"].unique())
    )

    gse_rows.append(row)

gq = pd.DataFrame(gse_rows)

gq = gq.sort_values(
    [
        "manual_priority_score",
        "hfd_duration_weeks_needs_review_n",
        "sex_needs_review_n",
        "strain_needs_review_n",
        "n_contrasts",
        "GSE",
    ],
    ascending=[False, False, False, False, False, True]
)

# ------------------------------------------------------------
# Conflict-first queue
# ------------------------------------------------------------

conflict_mask = (
    (x["hfd_duration_weeks_status"] == "CONFLICT")
    | (x["sex_status"] == "CONFLICT")
    | (x["strain_status"] == "CONFLICT")
    | (x["hfd_fat_kcal_percent_status"] == "CONFLICT")
    | (x["control_diet_type_status"] == "CONFLICT")
)

conflict = out.loc[conflict_mask].copy()

# ------------------------------------------------------------
# Suspicious auto values queue
# ------------------------------------------------------------

suspicious_rows = []

for _, r in x.iterrows():
    reasons = []

    # HFD fat values biologically suspicious for a conventional HFD label.
    val = r.get("hfd_fat_kcal_percent_value", "")
    try:
        fv = float(val)
    except Exception:
        fv = None
    if fv is not None and fv < 20:
        reasons.append(f"HFD_FAT_PERCENT_LT20:{fv:g}")

    # Strain strings likely typo-like but must NOT be normalized silently.
    strain = r.get("strain_value", "")
    if strain in {"C56BL/6J", "C56B6/J"}:
        reasons.append(f"STRAIN_TYPO_LIKE:{strain}")

    if reasons:
        suspicious_rows.append({
            "GSE": r["GSE"],
            "Contrast_id": r["Contrast_id"],
            "reason": ";".join(reasons),
            "sex_geo_auto": r.get("sex_value", ""),
            "strain_geo_auto": r.get("strain_value", ""),
            "hfd_duration_weeks_geo_auto": r.get(
                "hfd_duration_weeks_value", ""
            ),
            "hfd_fat_kcal_percent_geo_auto": r.get(
                "hfd_fat_kcal_percent_value", ""
            ),
            "control_diet_type_geo_auto": r.get(
                "control_diet_type_value", ""
            ),
        })

suspicious = pd.DataFrame(suspicious_rows)

# ------------------------------------------------------------
# Controlled vocabulary / instructions
# ------------------------------------------------------------

dictionary_rows = []

for var, rank in variables:
    dictionary_rows.append({
        "field_group": var,
        "priority_rank": rank,
        "definition": {
            "hfd_duration_weeks":
                "Time from initiation of HFD feeding to tissue collection "
                "for the specific frozen HFD-vs-control contrast.",
            "sex":
                "Biological sex of animals in the specific contrast.",
            "strain":
                "Mouse strain/substrain as explicitly supported by source evidence.",
            "hfd_fat_kcal_percent":
                "Percent of dietary energy/kcal from fat for the HFD arm; "
                "do not substitute weight percent.",
            "control_diet_type":
                "Diet type used for the frozen control arm "
                "(e.g., chow, matched purified low-fat diet).",
        }[var],
        "curation_rule": {
            "hfd_duration_weeks":
                "Do not confuse age, intervention duration, washout duration, "
                "or total study duration with HFD exposure duration.",
            "sex":
                "If repository and publication disagree, preserve both; "
                "expression-based sex may flag discordance but must not silently relabel.",
            "strain":
                "Preserve typo-like repository strings verbatim; normalize only "
                "after publication/supplement confirms the intended strain.",
            "hfd_fat_kcal_percent":
                "Prefer explicit kcal/energy percent or manufacturer diet code. "
                "Do not infer from generic '% fat' without unit context.",
            "control_diet_type":
                "Distinguish standard chow from matched purified low-fat/control diets "
                "when source evidence permits.",
        }[var],
    })

dictionary = pd.DataFrame(dictionary_rows)

vocab = pd.DataFrame({
    "category": (
        ["manual_source_type"] * len(source_categories)
        + ["adjudication_status"] * len(adjudication_status_categories)
        + ["confidence"] * len(confidence_categories)
    ),
    "allowed_value": (
        source_categories
        + adjudication_status_categories
        + confidence_categories
    )
})

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

template_out = OUTDIR / "HFD_primary105_moderator_manual_adjudication_template_v1.tsv"
gse_queue_out = OUTDIR / "HFD_primary64_publication_curation_queue_v1.tsv"
conflict_out = OUTDIR / "HFD_primary105_moderator_conflict_first_queue_v1.tsv"
suspicious_out = OUTDIR / "HFD_primary105_moderator_suspicious_auto_values_v1.tsv"
dict_out = OUTDIR / "HFD_moderator_manual_adjudication_dictionary_v1.tsv"
vocab_out = OUTDIR / "HFD_moderator_manual_adjudication_vocab_v1.tsv"

out.to_csv(template_out, sep="\t", index=False)
gq.to_csv(gse_queue_out, sep="\t", index=False)
conflict.to_csv(conflict_out, sep="\t", index=False)
suspicious.to_csv(suspicious_out, sep="\t", index=False)
dictionary.to_csv(dict_out, sep="\t", index=False)
vocab.to_csv(vocab_out, sep="\t", index=False)

# ------------------------------------------------------------
# QA / report
# ------------------------------------------------------------

print("===== HFD MODERATOR MANUAL ADJUDICATION INITIALIZATION PASS =====")
print("GSE:", out["GSE"].nunique())
print("Contrasts:", out["Contrast_id"].nunique())
print("Conflict-first contrasts:", len(conflict))
print("Suspicious auto-value contrasts:", len(suspicious))

print("\n===== TOP 20 GSE CURATION PRIORITY =====")
cols = [
    "GSE",
    "n_contrasts",
    "manual_priority_score",
    "hfd_duration_weeks_needs_review_n",
    "sex_needs_review_n",
    "strain_needs_review_n",
    "hfd_fat_kcal_percent_needs_review_n",
    "control_diet_type_needs_review_n",
    "hfd_fat_suspicious_value_n",
]
print(gq[cols].head(20).to_string(index=False))

print("\n===== OUTPUTS =====")
for p in [
    template_out,
    gse_queue_out,
    conflict_out,
    suspicious_out,
    dict_out,
    vocab_out,
]:
    print(p.relative_to(ROOT))

print("\nIMPORTANT:")
print("- Do NOT edit the GEO-auto v1.1.1 source table.")
print("- Fill only the manual_* / final_* / adjudication fields in the template.")
print("- Unresolved moderator metadata does not remove a valid primary contrast.")
print("- If CONTROL/HFD identity itself conflicts across sources, flag it for")
print("  separate sample-role adjudication; do not silently relabel.")
