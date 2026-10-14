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
OUTDIR = ROOT / "00_metadata/contrasts/moderators/manual_adjudication"

MASTER = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch05_v1.tsv"
if not MASTER.exists():
    raise FileNotFoundError(MASTER)

x = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")

assert x["GSE"].nunique() == 64
assert x["Contrast_id"].nunique() == 105

variables = [
    "hfd_duration_weeks",
    "sex",
    "strain",
    "hfd_fat_kcal_percent",
    "control_diet_type",
]

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
        return False

    if var == "hfd_fat_kcal_percent":
        try:
            fv = float(value)
        except Exception:
            return False
        if fv < 20 or fv > 80:
            return False

    if var == "strain" and value in {"C56BL/6J", "C56B6/J"}:
        return False

    return True

def norm_conf(v):
    v = str(v).strip().upper()
    return v if v in {"HIGH", "MEDIUM", "LOW", "UNRESOLVED"} else ""

rows = []

for _, r in x.iterrows():
    row = {
        "GSE": r["GSE"],
        "Contrast_id": r["Contrast_id"],
        "overall_adjudication_status": r.get("overall_adjudication_status", ""),
    }

    for var in variables:
        final = str(r.get(f"{var}_final", "")).strip()
        final_conf = norm_conf(r.get(f"{var}_final_confidence", ""))
        adj = str(r.get(f"{var}_adjudication_status", "")).strip()

        auto_val = str(r.get(f"{var}_geo_auto_value", "")).strip()
        auto_status = str(r.get(f"{var}_geo_auto_status", "")).strip()
        auto_conf = norm_conf(r.get(f"{var}_geo_auto_confidence", ""))

        if var in {"hfd_duration_weeks", "hfd_fat_kcal_percent"}:
            final = canon_number(final)
            auto_val = canon_number(auto_val)

        if final:
            value = final
            source = "MANUAL_FINAL"
            known = 1
            inclusive_ready = 1
            strict_ready = int(final_conf == "HIGH")
            confidence = final_conf or "UNSPECIFIED"
            curation_state = "REVIEWED_RESOLVED"

            if var == "sex" and final == "MIXED":
                inclusive_ready = 0
                strict_ready = 0

        elif adj in {"NOT_REPORTED", "CONFLICT_UNRESOLVED", "QC_DISCORDANT"}:
            value = "UNKNOWN"
            source = "MANUAL_REVIEW_UNRESOLVED"
            known = 0
            inclusive_ready = 0
            strict_ready = 0
            confidence = "UNRESOLVED"
            curation_state = "REVIEWED_UNRESOLVED"

        else:
            ok = safe_auto(var, auto_val, auto_status)

            if auto_status == "RESOLVED" and auto_val:
                value = auto_val
                source = "GEO_AUTO"
                known = 1
                inclusive_ready = int(ok)
                strict_ready = int(ok and auto_conf not in {"LOW", "UNRESOLVED"})
                confidence = auto_conf or "AUTO_UNSPECIFIED"
                curation_state = "AUTO_READY" if ok else "AUTO_REVIEW_NEEDED"
            else:
                value = "UNKNOWN"
                source = "UNRESOLVED"
                known = 0
                inclusive_ready = 0
                strict_ready = 0
                confidence = "UNRESOLVED"
                curation_state = "NOT_YET_REVIEWED"

        row[f"{var}_effective"] = value
        row[f"{var}_source_level"] = source
        row[f"{var}_confidence"] = confidence
        row[f"{var}_known"] = known
        row[f"{var}_inclusive_ready"] = inclusive_ready
        row[f"{var}_strict_ready"] = strict_ready
        row[f"{var}_curation_state"] = curation_state

    rows.append(row)

e = pd.DataFrame(rows)

# ------------------------------------------------------------
# GSE review status
# ------------------------------------------------------------
gse_status_rows = []

for gse, z in x.groupby("GSE", sort=True):
    statuses = set(s for s in z["overall_adjudication_status"].tolist() if s)

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
# Coverage
# ------------------------------------------------------------
coverage_rows = []

for var in variables:
    coverage_rows.append({
        "Moderator": var,
        "Total_contrasts": len(e),
        "Known_contrasts": int(e[f"{var}_known"].sum()),
        "Inclusive_ready_contrasts": int(e[f"{var}_inclusive_ready"].sum()),
        "Independent_GSE_inclusive_ready": e.loc[
            e[f"{var}_inclusive_ready"].eq(1), "GSE"
        ].nunique(),
        "Strict_ready_contrasts": int(e[f"{var}_strict_ready"].sum()),
        "Independent_GSE_strict_ready": e.loc[
            e[f"{var}_strict_ready"].eq(1), "GSE"
        ].nunique(),
        "Reviewed_but_unresolved_contrasts": int(
            (e[f"{var}_curation_state"] == "REVIEWED_UNRESOLVED").sum()
        ),
        "Pending_review_contrasts": int(
            e[f"{var}_curation_state"].isin(
                ["NOT_YET_REVIEWED", "AUTO_REVIEW_NEEDED"]
            ).sum()
        ),
    })

coverage = pd.DataFrame(coverage_rows)

# ------------------------------------------------------------
# Duration strict attrition after Batch05
# ------------------------------------------------------------
duration_gap = e[
    (e["hfd_duration_weeks_inclusive_ready"] == 1)
    & (e["hfd_duration_weeks_strict_ready"] == 0)
].copy()

duration_gap_gse = (
    duration_gap.groupby("GSE")
    .agg(
        n_gap_contrasts=("Contrast_id", "nunique"),
        duration_values=(
            "hfd_duration_weeks_effective",
            lambda s: ",".join(sorted(set(s)))
        ),
        confidence=(
            "hfd_duration_weeks_confidence",
            lambda s: ",".join(sorted(set(s)))
        ),
        source_level=(
            "hfd_duration_weeks_source_level",
            lambda s: ",".join(sorted(set(s)))
        ),
    )
    .reset_index()
)

duration_gap_gse = duration_gap_gse.merge(
    gse_status[["GSE", "review_status"]],
    on="GSE",
    how="left",
    validate="one_to_one",
)

# Remaining GEO_AUTO_LOW duration strict-upgrade candidates.
strict_upgrade = duration_gap_gse[
    (duration_gap_gse["confidence"] == "LOW")
    & (duration_gap_gse["source_level"] == "GEO_AUTO")
    & (duration_gap_gse["review_status"] == "NOT_YET_REVIEWED")
].copy()

# Add "other unresolved moderator" score only as tiebreaker.
other_vars = ["sex", "strain", "hfd_fat_kcal_percent", "control_diet_type"]
score_rows = []

for _, r in strict_upgrade.iterrows():
    gse = r["GSE"]
    z = e[e["GSE"].eq(gse)]

    other_pending = 0
    unresolved = []

    for var in other_vars:
        needs = (
            (z[f"{var}_curation_state"].isin(
                ["NOT_YET_REVIEWED", "AUTO_REVIEW_NEEDED"]
            ))
        ).any()

        if needs:
            other_pending += 1
            unresolved.append(var)

    score_rows.append({
        "GSE": gse,
        "other_pending_moderator_types": other_pending,
        "other_pending_moderators": ",".join(unresolved),
    })

if score_rows:
    score_df = pd.DataFrame(score_rows)
    strict_upgrade = strict_upgrade.merge(
        score_df, on="GSE", how="left", validate="one_to_one"
    )
else:
    strict_upgrade["other_pending_moderator_types"] = 0
    strict_upgrade["other_pending_moderators"] = ""

strict_upgrade = strict_upgrade.sort_values(
    [
        "other_pending_moderator_types",
        "n_gap_contrasts",
        "GSE",
    ],
    ascending=[False, False, True],
)

strict_upgrade.insert(
    0, "Next_duration_upgrade_order",
    range(1, len(strict_upgrade) + 1)
)

# ------------------------------------------------------------
# Batch05 delta audit: verify the intended four GSE are now strict-ready
# ------------------------------------------------------------
batch05_gses = ["GSE137365", "GSE207281", "GSE247985", "GSE282595"]

delta_rows = []
for gse in batch05_gses:
    z = e[e["GSE"].eq(gse)]
    delta_rows.append({
        "GSE": gse,
        "duration_value": ",".join(sorted(set(z["hfd_duration_weeks_effective"]))),
        "duration_confidence": ",".join(sorted(set(z["hfd_duration_weeks_confidence"]))),
        "strict_ready": int((z["hfd_duration_weeks_strict_ready"] == 1).all()),
        "review_status": gse_status.loc[
            gse_status["GSE"].eq(gse), "review_status"
        ].iloc[0],
    })

delta = pd.DataFrame(delta_rows)

assert (delta["strict_ready"] == 1).all(), delta

# Expected strict duration independent GSE = previous 20 + 4 upgrades = 24.
strict_duration_gse = coverage.loc[
    coverage["Moderator"].eq("hfd_duration_weeks"),
    "Independent_GSE_strict_ready"
].iloc[0]

assert strict_duration_gse == 24, (
    "Expected strict duration GSE=24 after four Batch05 upgrades, "
    f"observed={strict_duration_gse}"
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
eff_out = OUTDIR / "HFD_primary105_moderator_effective_after_batch05_v1.tsv"
cov_out = OUTDIR / "HFD_primary105_moderator_coverage_after_batch05_v1.tsv"
status_out = OUTDIR / "HFD_primary64_manual_review_status_after_batch05_v1.tsv"
gap_out = OUTDIR / "HFD_duration_strict_gap_after_batch05_v1.tsv"
queue_out = OUTDIR / "HFD_duration_next_strict_upgrade_queue_after_batch05_v1.tsv"
delta_out = OUTDIR / "HFD_batch05_duration_strict_upgrade_audit_v1.tsv"

e.to_csv(eff_out, sep="\t", index=False)
coverage.to_csv(cov_out, sep="\t", index=False)
gse_status.to_csv(status_out, sep="\t", index=False)
duration_gap_gse.to_csv(gap_out, sep="\t", index=False)
strict_upgrade.to_csv(queue_out, sep="\t", index=False)
delta.to_csv(delta_out, sep="\t", index=False)

print("===== HFD MODERATOR COVERAGE AFTER BATCH05 =====")
print(coverage.to_string(index=False))

print("\n===== BATCH05 DURATION STRICT-UPGRADE AUDIT =====")
print(delta.to_string(index=False))

print("\n===== DURATION STRICT STATUS =====")
dur = coverage[coverage["Moderator"].eq("hfd_duration_weeks")].iloc[0]
print("Inclusive-ready contrasts:", dur["Inclusive_ready_contrasts"])
print("Inclusive-ready independent GSE:", dur["Independent_GSE_inclusive_ready"])
print("Strict-ready contrasts:", dur["Strict_ready_contrasts"])
print("Strict-ready independent GSE:", dur["Independent_GSE_strict_ready"])
print("Remaining inclusive-not-strict GSE:", duration_gap_gse["GSE"].nunique())

print("\n===== NEXT DURATION STRICT-UPGRADE QUEUE =====")
if len(strict_upgrade):
    print(
        strict_upgrade[
            [
                "Next_duration_upgrade_order",
                "GSE",
                "n_gap_contrasts",
                "duration_values",
                "other_pending_moderator_types",
                "other_pending_moderators",
            ]
        ].to_string(index=False)
    )
else:
    print("No unreviewed GEO_AUTO_LOW duration GSE remain.")

print("\nSaved:")
for p in [
    eff_out, cov_out, status_out, gap_out, queue_out, delta_out
]:
    print(p.relative_to(ROOT))

print("\nINTERPRETATION:")
print("- Batch05 must raise strict duration independent GSE from 20 to 24.")
print("- Remaining strict gap excludes already-reviewed GSE from auto-requeue.")
print("- No moderator value or frozen eligibility is changed by this audit.")
