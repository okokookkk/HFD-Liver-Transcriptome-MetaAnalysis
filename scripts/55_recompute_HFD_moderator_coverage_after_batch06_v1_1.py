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

MASTER = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch06_v1_1.tsv"
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
        "overall_adjudication_status": r.get(
            "overall_adjudication_status", ""
        ),
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

        # Manual final has precedence.
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

        # Explicit reviewed-but-unresolved.
        elif adj in {
            "NOT_REPORTED",
            "CONFLICT_UNRESOLVED",
            "QC_DISCORDANT",
        }:
            value = "UNKNOWN"
            source = "MANUAL_REVIEW_UNRESOLVED"
            known = 0
            inclusive_ready = 0
            strict_ready = 0
            confidence = "UNRESOLVED"
            curation_state = "REVIEWED_UNRESOLVED"

        # Safe GEO-auto fallback.
        else:
            ok = safe_auto(var, auto_val, auto_status)

            if auto_status == "RESOLVED" and auto_val:
                value = auto_val
                source = "GEO_AUTO"
                known = 1
                inclusive_ready = int(ok)
                strict_ready = int(
                    ok and auto_conf not in {"LOW", "UNRESOLVED"}
                )
                confidence = auto_conf or "AUTO_UNSPECIFIED"
                curation_state = (
                    "AUTO_READY" if ok else "AUTO_REVIEW_NEEDED"
                )
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
# GSE-level manual review status
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
# Coverage
# ------------------------------------------------------------
coverage_rows = []

for var in variables:
    coverage_rows.append({
        "Moderator": var,
        "Total_contrasts": len(e),
        "Known_contrasts": int(e[f"{var}_known"].sum()),
        "Inclusive_ready_contrasts": int(
            e[f"{var}_inclusive_ready"].sum()
        ),
        "Independent_GSE_inclusive_ready": e.loc[
            e[f"{var}_inclusive_ready"].eq(1), "GSE"
        ].nunique(),
        "Strict_ready_contrasts": int(
            e[f"{var}_strict_ready"].sum()
        ),
        "Independent_GSE_strict_ready": e.loc[
            e[f"{var}_strict_ready"].eq(1), "GSE"
        ].nunique(),
        "Reviewed_but_unresolved_contrasts": int(
            (
                e[f"{var}_curation_state"]
                == "REVIEWED_UNRESOLVED"
            ).sum()
        ),
        "Pending_review_contrasts": int(
            e[f"{var}_curation_state"].isin(
                ["NOT_YET_REVIEWED", "AUTO_REVIEW_NEEDED"]
            ).sum()
        ),
    })

coverage = pd.DataFrame(coverage_rows)

# ------------------------------------------------------------
# Batch06 duration audit
# ------------------------------------------------------------
batch06_expected = {
    "GSE293219": ("16", "HIGH", 1),
    "GSE268200": ("27", "HIGH", 1),
    "GSE191004": ("4", "HIGH", 1),
    "GSE205846": ("15", "HIGH", 1),
    "GSE285614": ("13.0446", "MEDIUM", 0),
}

audit_rows = []

for gse, (expected_value, expected_conf, expected_strict) in batch06_expected.items():
    z = e[e["GSE"].eq(gse)].copy()
    if z.empty:
        raise RuntimeError(f"{gse}: missing from effective table")

    vals = sorted(set(z["hfd_duration_weeks_effective"]))
    confs = sorted(set(z["hfd_duration_weeks_confidence"]))
    strict_all = int(
        (z["hfd_duration_weeks_strict_ready"] == 1).all()
    )

    if vals != [expected_value]:
        raise AssertionError(
            f"{gse}: expected duration {expected_value}, observed {vals}"
        )
    if confs != [expected_conf]:
        raise AssertionError(
            f"{gse}: expected confidence {expected_conf}, observed {confs}"
        )
    if strict_all != expected_strict:
        raise AssertionError(
            f"{gse}: expected strict={expected_strict}, observed={strict_all}"
        )

    audit_rows.append({
        "GSE": gse,
        "duration_value": vals[0],
        "duration_confidence": confs[0],
        "strict_ready": strict_all,
        "review_status": gse_status.loc[
            gse_status["GSE"].eq(gse), "review_status"
        ].iloc[0],
    })

batch06_audit = pd.DataFrame(audit_rows)

# ------------------------------------------------------------
# Duration gap after Batch06
# ------------------------------------------------------------
duration_gap = e[
    (e["hfd_duration_weeks_inclusive_ready"] == 1)
    & (e["hfd_duration_weeks_strict_ready"] == 0)
].copy()

gap_rows = []

for gse, z in duration_gap.groupby("GSE", sort=True):
    gap_rows.append({
        "GSE": gse,
        "n_gap_contrasts": z["Contrast_id"].nunique(),
        "duration_values": ",".join(
            sorted(set(z["hfd_duration_weeks_effective"]))
        ),
        "confidence": ",".join(
            sorted(set(z["hfd_duration_weeks_confidence"]))
        ),
        "source_level": ",".join(
            sorted(set(z["hfd_duration_weeks_source_level"]))
        ),
        "review_status": gse_status.loc[
            gse_status["GSE"].eq(gse), "review_status"
        ].iloc[0],
    })

duration_gap_gse = pd.DataFrame(gap_rows)

# ------------------------------------------------------------
# Strict attrition reason summary
# ------------------------------------------------------------
reason_rows = []

if len(duration_gap):
    tmp = duration_gap.copy()
    tmp["gap_reason"] = (
        tmp["hfd_duration_weeks_source_level"]
        + "_"
        + tmp["hfd_duration_weeks_confidence"]
    )

    reason_summary = (
        tmp.groupby("gap_reason")
        .agg(
            Contrasts=("Contrast_id", "nunique"),
            Independent_GSE=("GSE", "nunique"),
        )
        .reset_index()
        .sort_values(
            ["Independent_GSE", "Contrasts"],
            ascending=False
        )
    )
else:
    reason_summary = pd.DataFrame(
        columns=["gap_reason", "Contrasts", "Independent_GSE"]
    )

# ------------------------------------------------------------
# Next unreviewed GEO_AUTO_LOW strict-upgrade queue
# ------------------------------------------------------------
remaining_low = duration_gap_gse[
    (duration_gap_gse["confidence"] == "LOW")
    & (duration_gap_gse["source_level"] == "GEO_AUTO")
    & (duration_gap_gse["review_status"] == "NOT_YET_REVIEWED")
].copy()

other_vars = [
    "sex",
    "strain",
    "hfd_fat_kcal_percent",
    "control_diet_type",
]

score_rows = []

for _, rr in remaining_low.iterrows():
    gse = rr["GSE"]
    z = e[e["GSE"].eq(gse)]

    pending = []
    for var in other_vars:
        needs = z[f"{var}_curation_state"].isin(
            ["NOT_YET_REVIEWED", "AUTO_REVIEW_NEEDED"]
        ).any()
        if needs:
            pending.append(var)

    score_rows.append({
        "GSE": gse,
        "other_pending_moderator_types": len(pending),
        "other_pending_moderators": ",".join(pending),
    })

if score_rows:
    score_df = pd.DataFrame(score_rows)
    remaining_low = remaining_low.merge(
        score_df,
        on="GSE",
        how="left",
        validate="one_to_one",
    )

remaining_low = remaining_low.sort_values(
    [
        "other_pending_moderator_types",
        "n_gap_contrasts",
        "GSE",
    ],
    ascending=[False, False, True],
)

if len(remaining_low):
    remaining_low.insert(
        0,
        "Next_duration_upgrade_order",
        range(1, len(remaining_low) + 1),
    )

# ------------------------------------------------------------
# Reviewed MEDIUM duration set
# ------------------------------------------------------------
manual_medium = duration_gap_gse[
    (duration_gap_gse["confidence"] == "MEDIUM")
    & (duration_gap_gse["source_level"] == "MANUAL_FINAL")
].copy()

# ------------------------------------------------------------
# Hard QA expectations after Batch06
# ------------------------------------------------------------
dur_cov = coverage[
    coverage["Moderator"].eq("hfd_duration_weeks")
].iloc[0]

assert int(dur_cov["Independent_GSE_inclusive_ready"]) == 36, dur_cov
assert int(dur_cov["Independent_GSE_strict_ready"]) == 28, dur_cov
assert int(dur_cov["Strict_ready_contrasts"]) == 65, dur_cov

assert duration_gap_gse["GSE"].nunique() == 8, duration_gap_gse

# Expected composition:
# 6 remaining unreviewed GEO_AUTO_LOW + 2 reviewed MANUAL_MEDIUM.
assert len(remaining_low) == 6, remaining_low
assert manual_medium["GSE"].nunique() == 2, manual_medium

expected_medium = {"GSE285380", "GSE285614"}
observed_medium = set(manual_medium["GSE"])
assert observed_medium == expected_medium, (
    expected_medium, observed_medium
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
eff_out = OUTDIR / "HFD_primary105_moderator_effective_after_batch06_v1_1.tsv"
cov_out = OUTDIR / "HFD_primary105_moderator_coverage_after_batch06_v1_1.tsv"
status_out = OUTDIR / "HFD_primary64_manual_review_status_after_batch06_v1_1.tsv"
audit_out = OUTDIR / "HFD_batch06_duration_strict_upgrade_audit_v1_1.tsv"
gap_out = OUTDIR / "HFD_duration_strict_gap_after_batch06_v1_1.tsv"
reason_out = OUTDIR / "HFD_duration_strict_attrition_reason_summary_after_batch06_v1_1.tsv"
queue_out = OUTDIR / "HFD_duration_next_strict_upgrade_queue_after_batch06_v1_1.tsv"
medium_out = OUTDIR / "HFD_duration_manual_medium_set_after_batch06_v1_1.tsv"

e.to_csv(eff_out, sep="\t", index=False)
coverage.to_csv(cov_out, sep="\t", index=False)
gse_status.to_csv(status_out, sep="\t", index=False)
batch06_audit.to_csv(audit_out, sep="\t", index=False)
duration_gap_gse.to_csv(gap_out, sep="\t", index=False)
reason_summary.to_csv(reason_out, sep="\t", index=False)
remaining_low.to_csv(queue_out, sep="\t", index=False)
manual_medium.to_csv(medium_out, sep="\t", index=False)

print("===== HFD MODERATOR COVERAGE AFTER BATCH06 v1.1 =====")
print(coverage.to_string(index=False))

print("\n===== BATCH06 DURATION AUDIT =====")
print(batch06_audit.to_string(index=False))

print("\n===== DURATION STRICT STATUS =====")
print(
    "Inclusive-ready independent GSE:",
    int(dur_cov["Independent_GSE_inclusive_ready"]),
)
print(
    "Strict-ready independent GSE:",
    int(dur_cov["Independent_GSE_strict_ready"]),
)
print(
    "Strict-ready contrasts:",
    int(dur_cov["Strict_ready_contrasts"]),
)
print(
    "Remaining inclusive-not-strict GSE:",
    duration_gap_gse["GSE"].nunique(),
)

print("\n===== DURATION STRICT ATTRITION BY REASON =====")
print(reason_summary.to_string(index=False))

print("\n===== MANUAL MEDIUM DURATION SET =====")
print(
    manual_medium[
        [
            "GSE",
            "duration_values",
            "confidence",
            "review_status",
        ]
    ].to_string(index=False)
)

print("\n===== NEXT DURATION STRICT-UPGRADE QUEUE =====")
if len(remaining_low):
    print(
        remaining_low[
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
    eff_out,
    cov_out,
    status_out,
    audit_out,
    gap_out,
    reason_out,
    queue_out,
    medium_out,
]:
    print(p.relative_to(ROOT))

print("\nINTERPRETATION:")
print("- Batch06 must raise strict duration independent GSE from 24 to 28.")
print("- GSE268200 final duration must be 27 weeks, not 21 or 28.")
print("- GSE285614 must remain MEDIUM and non-strict for duration.")
print("- Remaining duration gap must be 6 GEO_AUTO_LOW + 2 MANUAL_MEDIUM GSE.")
print("- No moderator value, eligibility, or frozen contrast membership is changed by this audit.")
