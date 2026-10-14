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

MASTER = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch04_v1.tsv"
EFFECTIVE = OUTDIR / "HFD_primary105_moderator_effective_after_batch04_v1.tsv"
STATUS = OUTDIR / "HFD_primary64_manual_review_status_after_batch04_v1.tsv"

for p in [MASTER, EFFECTIVE, STATUS]:
    if not p.exists():
        raise FileNotFoundError(p)

m = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")
e = pd.read_csv(EFFECTIVE, sep="\t", dtype=str).fillna("")
s = pd.read_csv(STATUS, sep="\t", dtype=str).fillna("")

assert m["Contrast_id"].nunique() == 105
assert e["Contrast_id"].nunique() == 105
assert s["GSE"].nunique() == 64

# ------------------------------------------------------------
# Contrast-level audit of inclusive-ready but strict-not-ready duration.
# ------------------------------------------------------------
gap = e[
    (e["hfd_duration_weeks_inclusive_ready"] == "1")
    & (e["hfd_duration_weeks_strict_ready"] == "0")
].copy()

gap = gap.merge(
    m[
        [
            "Contrast_id",
            "GSE",
            "hfd_duration_weeks_final",
            "hfd_duration_weeks_final_confidence",
            "hfd_duration_weeks_final_source",
            "hfd_duration_weeks_adjudication_status",
            "hfd_duration_weeks_geo_auto_value",
            "hfd_duration_weeks_geo_auto_status",
            "hfd_duration_weeks_geo_auto_confidence",
            "hfd_duration_weeks_geo_auto_source",
            "hfd_duration_weeks_geo_auto_raw",
        ]
    ],
    on=["Contrast_id", "GSE"],
    how="left",
    validate="one_to_one",
)

def classify_reason(r):
    if r["hfd_duration_weeks_source_level"] == "MANUAL_FINAL":
        conf = r["hfd_duration_weeks_confidence"]
        return f"MANUAL_{conf or 'UNSPECIFIED'}"
    if r["hfd_duration_weeks_source_level"] == "GEO_AUTO":
        conf = r["hfd_duration_weeks_confidence"]
        return f"GEO_AUTO_{conf or 'UNSPECIFIED'}"
    return r["hfd_duration_weeks_source_level"] or "OTHER"

gap["strict_gap_reason"] = gap.apply(classify_reason, axis=1)

# ------------------------------------------------------------
# GSE-level strict-gap summary.
# ------------------------------------------------------------
gse_rows = []

for gse, z in gap.groupby("GSE", sort=True):
    review_status = s.loc[s["GSE"].eq(gse), "review_status"].iloc[0]

    values = sorted(set(z["hfd_duration_weeks_effective"]))
    reasons = sorted(set(z["strict_gap_reason"]))
    auto_conf = sorted(
        set(v for v in z["hfd_duration_weeks_geo_auto_confidence"] if v)
    )
    auto_src = sorted(
        set(v for v in z["hfd_duration_weeks_geo_auto_source"] if v)
    )

    gse_rows.append({
        "GSE": gse,
        "review_status": review_status,
        "n_gap_contrasts": z["Contrast_id"].nunique(),
        "duration_values": ",".join(values),
        "strict_gap_reasons": ",".join(reasons),
        "geo_auto_confidences": ",".join(auto_conf),
        "geo_auto_sources": " || ".join(auto_src),
    })

gse_gap = pd.DataFrame(gse_rows)

# ------------------------------------------------------------
# Strict attrition summary by reason.
# ------------------------------------------------------------
reason_summary = (
    gap.groupby("strict_gap_reason")
    .agg(
        Contrasts=("Contrast_id", "nunique"),
        Independent_GSE=("GSE", "nunique"),
    )
    .reset_index()
    .sort_values(["Independent_GSE", "Contrasts"], ascending=False)
)

# ------------------------------------------------------------
# Duration-focused manual curation queue.
#
# Priority:
# A) Unreviewed GSE with inclusive-ready but strict-not-ready duration:
#    one manual source check can convert a usable but low-confidence duration
#    into HIGH-confidence strict-ready evidence.
#
# B) Unreviewed GSE with duration not inclusive-ready:
#    source review may add a new duration GSE.
#
# Reviewed partial GSE are listed separately, never automatically requeued.
# ------------------------------------------------------------
all_gse = []

for gse, z in e.groupby("GSE", sort=True):
    review_status = s.loc[s["GSE"].eq(gse), "review_status"].iloc[0]

    incl_any = (z["hfd_duration_weeks_inclusive_ready"] == "1").any()
    strict_any = (z["hfd_duration_weeks_strict_ready"] == "1").any()

    pending_other = 0
    for var in ["sex", "strain", "hfd_fat_kcal_percent", "control_diet_type"]:
        pending_other += int(
            z[f"{var}_curation_state"].isin(
                ["NOT_YET_REVIEWED", "AUTO_REVIEW_NEEDED"]
            ).any()
        )

    if incl_any and not strict_any:
        duration_target_class = "A_STRICT_UPGRADE"
        base_score = 100
    elif not incl_any:
        duration_target_class = "B_NEW_DURATION"
        base_score = 50
    else:
        duration_target_class = "C_ALREADY_STRICT"
        base_score = 0

    all_gse.append({
        "GSE": gse,
        "review_status": review_status,
        "n_contrasts": z["Contrast_id"].nunique(),
        "duration_target_class": duration_target_class,
        "duration_inclusive_any": int(incl_any),
        "duration_strict_any": int(strict_any),
        "other_pending_moderator_types": pending_other,
        "duration_curation_priority_score": base_score + pending_other,
    })

duration_queue = pd.DataFrame(all_gse)

unreviewed_queue = duration_queue[
    (duration_queue["review_status"] == "NOT_YET_REVIEWED")
    & (duration_queue["duration_target_class"] != "C_ALREADY_STRICT")
].copy()

unreviewed_queue = unreviewed_queue.sort_values(
    [
        "duration_curation_priority_score",
        "other_pending_moderator_types",
        "n_contrasts",
        "GSE",
    ],
    ascending=[False, False, True, True],
)

# Reviewed-but-unresolved / partial records shown separately only.
reviewed_not_strict = duration_queue[
    (duration_queue["review_status"] != "NOT_YET_REVIEWED")
    & (duration_queue["duration_strict_any"] == 0)
].copy()

# Recommended next batch: strict-upgrade targets first.
batch05_duration = unreviewed_queue.head(12).copy()
if len(batch05_duration):
    batch05_duration.insert(
        0, "DurationBatch05_order",
        range(1, len(batch05_duration) + 1)
    )

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
gap_out = OUTDIR / "HFD_duration_inclusive_vs_strict_gap_contrasts_after_batch04_v1.tsv"
gse_gap_out = OUTDIR / "HFD_duration_inclusive_vs_strict_gap_GSE_after_batch04_v1.tsv"
reason_out = OUTDIR / "HFD_duration_strict_attrition_reason_summary_after_batch04_v1.tsv"
queue_out = OUTDIR / "HFD_duration_focused_curation_queue_after_batch04_v1.tsv"
reviewed_out = OUTDIR / "HFD_duration_reviewed_not_strict_after_batch04_v1.tsv"
batch_out = OUTDIR / "HFD_duration_focused_batch05_queue_v1.tsv"

gap.to_csv(gap_out, sep="\t", index=False)
gse_gap.to_csv(gse_gap_out, sep="\t", index=False)
reason_summary.to_csv(reason_out, sep="\t", index=False)
unreviewed_queue.to_csv(queue_out, sep="\t", index=False)
reviewed_not_strict.to_csv(reviewed_out, sep="\t", index=False)
batch05_duration.to_csv(batch_out, sep="\t", index=False)

print("===== HFD DURATION STRICT-ATTRITION AUDIT AFTER BATCH04 =====")
print("Inclusive-ready duration contrasts:",
      int((e["hfd_duration_weeks_inclusive_ready"] == "1").sum()))
print("Strict-ready duration contrasts:",
      int((e["hfd_duration_weeks_strict_ready"] == "1").sum()))
print("Gap contrasts:", gap["Contrast_id"].nunique())
print("Gap independent GSE:", gap["GSE"].nunique())

print("\n===== STRICT ATTRITION BY REASON =====")
print(reason_summary.to_string(index=False))

print("\n===== GSE LOST FROM STRICT DURATION SET =====")
print(
    gse_gap.sort_values(
        ["review_status", "n_gap_contrasts", "GSE"],
        ascending=[True, False, True]
    ).to_string(index=False)
)

print("\n===== DURATION-FOCUSED BATCH05 =====")
cols = [
    "DurationBatch05_order",
    "GSE",
    "n_contrasts",
    "duration_target_class",
    "duration_inclusive_any",
    "duration_strict_any",
    "other_pending_moderator_types",
    "duration_curation_priority_score",
]
print(batch05_duration[cols].to_string(index=False))

print("\n===== REVIEWED GSE NOT STRICT-READY FOR DURATION =====")
if len(reviewed_not_strict):
    print(
        reviewed_not_strict[
            [
                "GSE", "review_status", "n_contrasts",
                "duration_target_class",
                "duration_inclusive_any", "duration_strict_any",
            ]
        ].to_string(index=False)
    )
else:
    print("None")

print("\nSaved:")
for p in [
    gap_out, gse_gap_out, reason_out,
    queue_out, reviewed_out, batch_out
]:
    print(p.relative_to(ROOT))

print("\nINTERPRETATION:")
print("- The strict-vs-inclusive duration gap is audited before further curation.")
print("- A_STRICT_UPGRADE means duration already exists but is not strict-ready.")
print("- B_NEW_DURATION means duration is not currently analysis-ready.")
print("- Reviewed partial GSE are never automatically requeued.")
print("- No metadata value, eligibility, or frozen contrast membership is changed.")
