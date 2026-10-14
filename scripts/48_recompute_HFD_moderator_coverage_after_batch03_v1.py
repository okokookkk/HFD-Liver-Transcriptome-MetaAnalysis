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
MASTER = (
    ROOT /
    "00_metadata/contrasts/moderators/manual_adjudication/"
    "HFD_primary105_moderator_manual_adjudication_batch03_v1.tsv"
)
OUTDIR = ROOT / "00_metadata/contrasts/moderators/manual_adjudication"

if not MASTER.exists():
    raise FileNotFoundError(MASTER)

x = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")

assert x["Contrast_id"].nunique() == 105
assert x["GSE"].nunique() == 64

variables = [
    "hfd_duration_weeks",
    "sex",
    "strain",
    "hfd_fat_kcal_percent",
    "control_diet_type",
]

weights = {
    "hfd_duration_weeks": 5,
    "sex": 4,
    "strain": 3,
    "hfd_fat_kcal_percent": 2,
    "control_diet_type": 1,
}

def canon_number(v):
    v = str(v).strip()
    if v in {"", "UNKNOWN", "CONFLICT", "MIXED"}:
        return v
    try:
        f = float(v)
    except Exception:
        return v
    if f.is_integer():
        return str(int(f))
    return f"{f:g}"

def safe_auto(var, value, status):
    if status != "RESOLVED" or value == "":
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
    row = {
        "GSE": r["GSE"],
        "Contrast_id": r["Contrast_id"],
        "overall_adjudication_status": r.get(
            "overall_adjudication_status", ""
        ),
    }

    for var in variables:
        final = str(r.get(f"{var}_final", "")).strip()
        adj = str(r.get(f"{var}_adjudication_status", "")).strip()
        auto_val = str(r.get(f"{var}_geo_auto_value", "")).strip()
        auto_status = str(r.get(f"{var}_geo_auto_status", "")).strip()

        if var in {"hfd_duration_weeks", "hfd_fat_kcal_percent"}:
            final = canon_number(final)
            auto_val = canon_number(auto_val)

        # 1) Manual final wins.
        if final:
            value = final
            source = "MANUAL_FINAL"
            known = 1
            reviewed = 1
            curation_state = "REVIEWED_RESOLVED"

            # Mixed sex is known but not usable for binary subgroup.
            if var == "sex" and final == "MIXED":
                ready = 0
                readiness_reason = "KNOWN_MIXED_SEX"
            else:
                ready = 1
                readiness_reason = adj or "MANUAL_FINAL"

        # 2) Explicit manual NOT_REPORTED / unresolved review.
        elif adj in {
            "NOT_REPORTED",
            "CONFLICT_UNRESOLVED",
            "QC_DISCORDANT",
        }:
            value = "UNKNOWN"
            source = "MANUAL_REVIEW_UNRESOLVED"
            known = 0
            ready = 0
            reviewed = 1
            curation_state = "REVIEWED_UNRESOLVED"
            readiness_reason = adj

        # 3) Otherwise use safe resolved GEO-auto evidence.
        else:
            ok, why = safe_auto(var, auto_val, auto_status)
            if auto_status == "RESOLVED" and auto_val:
                value = auto_val
                source = "GEO_AUTO"
                known = 1
                ready = int(ok)
                reviewed = int(adj not in {"", "NOT_REVIEWED"})
                curation_state = (
                    "AUTO_READY"
                    if ok else "AUTO_REVIEW_NEEDED"
                )
                readiness_reason = why
            else:
                value = "UNKNOWN"
                source = "UNRESOLVED"
                known = 0
                ready = 0
                reviewed = int(adj not in {"", "NOT_REVIEWED"})
                curation_state = (
                    "REVIEWED_UNRESOLVED"
                    if reviewed else "NOT_YET_REVIEWED"
                )
                readiness_reason = auto_status or "MISSING"

        row[f"{var}_effective"] = value
        row[f"{var}_source_level"] = source
        row[f"{var}_known"] = known
        row[f"{var}_analysis_ready"] = ready
        row[f"{var}_reviewed"] = reviewed
        row[f"{var}_curation_state"] = curation_state
        row[f"{var}_readiness_reason"] = readiness_reason

    rows.append(row)

e = pd.DataFrame(rows)

# ------------------------------------------------------------
# Coverage
# ------------------------------------------------------------
coverage_rows = []
for var in variables:
    known = int(e[f"{var}_known"].sum())
    ready = int(e[f"{var}_analysis_ready"].sum())
    reviewed_unresolved = int(
        (e[f"{var}_curation_state"] == "REVIEWED_UNRESOLVED").sum()
    )
    pending_review = int(
        e[f"{var}_curation_state"].isin(
            ["NOT_YET_REVIEWED", "AUTO_REVIEW_NEEDED"]
        ).sum()
    )

    coverage_rows.append({
        "Moderator": var,
        "Total_contrasts": len(e),
        "Known_contrasts": known,
        "Known_percent": 100.0 * known / len(e),
        "Analysis_ready_contrasts": ready,
        "Analysis_ready_percent": 100.0 * ready / len(e),
        "Independent_GSE_with_ready_value": e.loc[
            e[f"{var}_analysis_ready"].eq(1), "GSE"
        ].nunique(),
        "Reviewed_but_unresolved_contrasts": reviewed_unresolved,
        "Pending_review_contrasts": pending_review,
    })
coverage = pd.DataFrame(coverage_rows)

# ------------------------------------------------------------
# Canonical distributions
# ------------------------------------------------------------
dist_rows = []
for var in variables:
    z = e[e[f"{var}_known"].eq(1)].copy()
    for value, q in z.groupby(f"{var}_effective", sort=True):
        dist_rows.append({
            "Moderator": var,
            "Value": value,
            "Contrasts": q["Contrast_id"].nunique(),
            "Independent_GSE": q["GSE"].nunique(),
            "Analysis_ready_contrasts": int(
                q[f"{var}_analysis_ready"].sum()
            ),
        })
dist = pd.DataFrame(dist_rows)

# ------------------------------------------------------------
# Manual review progress
# Full adjudication vs partial reviewed.
# ------------------------------------------------------------
gse_status_rows = []
for gse, z in x.groupby("GSE", sort=True):
    statuses = set(
        s for s in z["overall_adjudication_status"].tolist()
        if s
    )

    if any("PARTIAL_REVIEWED" in s for s in statuses):
        review_status = "PARTIAL_REVIEWED"
    elif statuses and all("ADJUDICATED" in s for s in statuses):
        review_status = "FULLY_ADJUDICATED"
    elif any("ADJUDICATED" in s for s in statuses):
        review_status = "PARTIALLY_ADJUDICATED"
    else:
        review_status = "NOT_YET_REVIEWED"

    gse_status_rows.append({
        "GSE": gse,
        "n_contrasts": z["Contrast_id"].nunique(),
        "review_status": review_status,
    })

gse_status = pd.DataFrame(gse_status_rows)

# ------------------------------------------------------------
# Batch04 queue:
# Only NOT-YET-REVIEWED GSE are candidates.
# A moderator scores once per GSE if any contrast still needs review.
# Reviewed-but-unresolved fields do not re-enter the queue.
# ------------------------------------------------------------
queue_rows = []

for gse, z in e.groupby("GSE", sort=True):
    review_status = gse_status.loc[
        gse_status["GSE"].eq(gse), "review_status"
    ].iloc[0]

    if review_status != "NOT_YET_REVIEWED":
        continue

    row = {
        "GSE": gse,
        "n_contrasts": z["Contrast_id"].nunique(),
    }

    score = 0
    unresolved_types = []

    for var in variables:
        states = set(z[f"{var}_curation_state"].tolist())

        # Needs manual curation if any contrast in this GSE is either
        # not reviewed or has suspicious auto-only metadata.
        needs = any(
            s in {"NOT_YET_REVIEWED", "AUTO_REVIEW_NEEDED"}
            for s in states
        )

        row[f"{var}_needs_review"] = int(needs)

        if needs:
            score += weights[var]
            unresolved_types.append(var)

    row["gse_marginal_gain_score"] = score
    row["n_unresolved_moderator_types"] = len(unresolved_types)
    row["unresolved_moderators"] = ",".join(unresolved_types)

    queue_rows.append(row)

queue = pd.DataFrame(queue_rows)
if len(queue):
    queue = queue[queue["gse_marginal_gain_score"] > 0].copy()
    queue = queue.sort_values(
        [
            "gse_marginal_gain_score",
            "n_unresolved_moderator_types",
            "n_contrasts",
            "GSE",
        ],
        ascending=[False, False, True, True],
    )

batch04 = queue.head(12).copy()
if len(batch04):
    batch04.insert(
        0, "Batch04_order", range(1, len(batch04) + 1)
    )

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------
summary = pd.DataFrame([{
    "Total_GSE": x["GSE"].nunique(),
    "Total_contrasts": x["Contrast_id"].nunique(),
    "Fully_adjudicated_GSE": int(
        (gse_status["review_status"] == "FULLY_ADJUDICATED").sum()
    ),
    "Partial_reviewed_GSE": int(
        (gse_status["review_status"] == "PARTIAL_REVIEWED").sum()
    ),
    "Not_yet_reviewed_GSE": int(
        (gse_status["review_status"] == "NOT_YET_REVIEWED").sum()
    ),
    "Fully_adjudicated_contrasts": int(
        gse_status.loc[
            gse_status["review_status"].eq("FULLY_ADJUDICATED"),
            "n_contrasts"
        ].sum()
    ),
    "Partial_reviewed_contrasts": int(
        gse_status.loc[
            gse_status["review_status"].eq("PARTIAL_REVIEWED"),
            "n_contrasts"
        ].sum()
    ),
}])

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
eff_out = OUTDIR / "HFD_primary105_moderator_effective_after_batch03_v1.tsv"
cov_out = OUTDIR / "HFD_primary105_moderator_coverage_after_batch03_v1.tsv"
dist_out = OUTDIR / "HFD_primary105_moderator_distribution_after_batch03_v1.tsv"
gse_status_out = OUTDIR / "HFD_primary64_manual_review_status_after_batch03_v1.tsv"
queue_out = OUTDIR / "HFD_primary64_batch04_curation_queue_v1.tsv"
batch_out = OUTDIR / "HFD_batch04_recommended_GSE_queue_v1.tsv"
summary_out = OUTDIR / "HFD_moderator_manual_curation_progress_after_batch03_v1.tsv"

e.to_csv(eff_out, sep="\t", index=False)
coverage.to_csv(cov_out, sep="\t", index=False)
dist.to_csv(dist_out, sep="\t", index=False)
gse_status.to_csv(gse_status_out, sep="\t", index=False)
queue.to_csv(queue_out, sep="\t", index=False)
batch04.to_csv(batch_out, sep="\t", index=False)
summary.to_csv(summary_out, sep="\t", index=False)

print("===== HFD MODERATOR COVERAGE AFTER BATCH03 =====")
print(summary.to_string(index=False))

print("\n===== EFFECTIVE COVERAGE =====")
print(coverage.to_string(index=False))

print("\n===== KEY DISTRIBUTIONS =====")
for var in variables:
    print(f"\n--- {var} ---")
    z = dist[dist["Moderator"].eq(var)]
    print(
        z.sort_values(
            ["Independent_GSE", "Contrasts"],
            ascending=False
        ).to_string(index=False)
    )

print("\n===== REVIEW STATUS =====")
print(
    gse_status["review_status"]
    .value_counts()
    .to_string()
)

print("\n===== RECOMMENDED BATCH04 =====")
if len(batch04):
    print(
        batch04[
            [
                "Batch04_order",
                "GSE",
                "n_contrasts",
                "gse_marginal_gain_score",
                "n_unresolved_moderator_types",
                "unresolved_moderators",
            ]
        ].to_string(index=False)
    )
else:
    print("No unreviewed GSE remain.")

print("\nSaved:")
for p in [
    eff_out,
    cov_out,
    dist_out,
    gse_status_out,
    queue_out,
    batch_out,
    summary_out,
]:
    print(p.relative_to(ROOT))

print("\nRULE:")
print("- Reviewed-but-unresolved fields are NOT requeued automatically.")
print("- Unknown remains unknown; no biological inference is added.")
print("- Numeric normalization is reporting-only.")
print("- No frozen eligibility or contrast membership is changed.")
