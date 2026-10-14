#!/usr/bin/env python3
from pathlib import Path
import argparse
import pandas as pd

# ============================================================
# Recompute moderator coverage after canonical Batch08 v1.
# Also writes:
#   - inclusive-but-not-strict duration queue
#   - reviewed-unresolved duration source-gap table
# ============================================================

ap = argparse.ArgumentParser()
ap.add_argument(
    "--root",
    default=__import__("os").environ.get("HFD_PROJECT_ROOT", "."),
)
args = ap.parse_args()

ROOT = Path(args.root)
OUTDIR = ROOT / "00_metadata/contrasts/moderators/manual_adjudication"

MASTER = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch08_v1.tsv"
COVERAGE_OUT = OUTDIR / "HFD_primary105_moderator_coverage_after_batch08_v1.tsv"
QUEUE_OUT = OUTDIR / "HFD_duration_next_strict_upgrade_queue_after_batch08_v1.tsv"
UNRESOLVED_OUT = OUTDIR / "HFD_duration_reviewed_unresolved_source_gaps_after_batch08_v1.tsv"

if not MASTER.exists():
    raise FileNotFoundError(MASTER)

x = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")
assert x["GSE"].nunique() == 64
assert x["Contrast_id"].nunique() == 105

MODS = [
    "hfd_duration_weeks",
    "sex",
    "strain",
    "hfd_fat_kcal_percent",
    "control_diet_type",
]

VALID_SEX = {"MALE", "FEMALE"}
VALID_STRAINS = {
    "C57BL/6J", "C57BL/6", "C57BL/6NHsd",
    "C57BL/6N", "C57BL/6NTac",
    "129S1/SvImJ", "FVB.129P2",
}
VALID_CONTROL = {
    "CHOW", "LFD", "NORMAL_DIET", "CONTROL_DIET_UNSPECIFIED",
}

def effective(row, var):
    # Manual adjudication supersedes GEO-auto whenever reviewed.
    if row[f"{var}_adjudication_status"] != "NOT_REVIEWED":
        return (
            row[f"{var}_final"].strip(),
            row[f"{var}_final_confidence"].strip(),
            "MANUAL",
        )
    if row[f"{var}_geo_auto_status"] == "RESOLVED":
        return (
            row[f"{var}_geo_auto_value"].strip(),
            row[f"{var}_geo_auto_confidence"].strip(),
            "GEO_AUTO",
        )
    return "", "", "NONE"

def inclusive_ready(var, val):
    if not val or val in {"UNKNOWN", "UNRESOLVED", "NOT_REPORTED", "NA", "N/A"}:
        return False

    if var == "hfd_duration_weeks":
        try:
            return float(val) > 0
        except ValueError:
            return False

    if var == "hfd_fat_kcal_percent":
        try:
            z = float(val)
            return 0 < z <= 100
        except ValueError:
            return False

    if var == "sex":
        # MIXED is biologically known but is not ready for the binary
        # male-vs-female moderator subgroup.
        return val in VALID_SEX

    if var == "strain":
        return val in VALID_STRAINS

    if var == "control_diet_type":
        return val in VALID_CONTROL

    return False

def strict_ready(var, val, confidence, origin):
    if not inclusive_ready(var, val):
        return False
    if origin == "MANUAL":
        return confidence == "HIGH"
    if origin == "GEO_AUTO":
        return confidence in {"MEDIUM", "HIGH"}
    return False

coverage = []
for var in MODS:
    known_contrasts = 0
    inclusive_contrasts = 0
    strict_contrasts = 0
    reviewed_unresolved = 0
    pending_review = 0
    inclusive_gse = set()
    strict_gse = set()

    for _, row in x.iterrows():
        val, conf, origin = effective(row, var)
        inc = inclusive_ready(var, val)
        strict = strict_ready(var, val, conf, origin)

        if val:
            known_contrasts += 1
        if inc:
            inclusive_contrasts += 1
            inclusive_gse.add(row["GSE"])
        if strict:
            strict_contrasts += 1
            strict_gse.add(row["GSE"])

        if row[f"{var}_manual_confidence"] == "UNRESOLVED":
            reviewed_unresolved += 1

        if row[f"{var}_adjudication_status"] == "NOT_REVIEWED" and not inc:
            pending_review += 1

    coverage.append({
        "Moderator": var,
        "Total_contrasts": len(x),
        "Known_contrasts": known_contrasts,
        "Inclusive_ready_contrasts": inclusive_contrasts,
        "Independent_GSE_inclusive_ready": len(inclusive_gse),
        "Strict_ready_contrasts": strict_contrasts,
        "Independent_GSE_strict_ready": len(strict_gse),
        "Reviewed_but_unresolved_contrasts": reviewed_unresolved,
        "Pending_review_contrasts": pending_review,
    })

cov = pd.DataFrame(coverage)
cov.to_csv(COVERAGE_OUT, sep="\t", index=False)

# ------------------------------------------------------------
# Duration inclusive-but-not-strict queue
# ------------------------------------------------------------
queue_rows = []
for gse, z in x.groupby("GSE"):
    gap_rows = []
    for _, row in z.iterrows():
        val, conf, origin = effective(row, "hfd_duration_weeks")
        if (
            inclusive_ready("hfd_duration_weeks", val)
            and not strict_ready("hfd_duration_weeks", val, conf, origin)
        ):
            gap_rows.append((row, val, conf, origin))

    if gap_rows:
        vals = sorted({a[1] for a in gap_rows}, key=float)
        confs = sorted({a[2] for a in gap_rows})
        origins = sorted({a[3] for a in gap_rows})
        queue_rows.append({
            "GSE": gse,
            "n_gap_contrasts": len(gap_rows),
            "duration_values": ",".join(vals),
            "confidence": ",".join(confs),
            "source_level": ",".join(origins),
            "review_status": (
                "REVIEWED_MEDIUM_STRICT_GAP"
                if origins == ["MANUAL"]
                else "GEO_AUTO_LOW_STRICT_GAP"
            ),
        })

queue = pd.DataFrame(queue_rows)
expected_queue = ["GSE285614", "GSE324134"]
observed_queue = sorted(queue["GSE"].tolist()) if len(queue) else []
assert observed_queue == expected_queue, (
    "Unexpected duration strict-gap queue",
    observed_queue,
    expected_queue,
)

order = {g: i + 1 for i, g in enumerate(expected_queue)}
queue.insert(
    0,
    "Next_duration_upgrade_order",
    queue["GSE"].map(order).astype(int),
)
queue = queue.sort_values("Next_duration_upgrade_order")
queue.to_csv(QUEUE_OUT, sep="\t", index=False)

# ------------------------------------------------------------
# Reviewed-unresolved duration source gaps
# ------------------------------------------------------------
unresolved_rows = []
for _, row in x.iterrows():
    if row["hfd_duration_weeks_manual_confidence"] == "UNRESOLVED":
        unresolved_rows.append({
            "GSE": row["GSE"],
            "Contrast_id": row["Contrast_id"],
            "adjudication_status": row["hfd_duration_weeks_adjudication_status"],
            "source_type": row["hfd_duration_weeks_manual_source_type"],
            "source_citation": row["hfd_duration_weeks_manual_source_citation"],
            "evidence": row["hfd_duration_weeks_manual_evidence"],
            "final_note": row["hfd_duration_weeks_final_note"],
        })

unresolved = pd.DataFrame(unresolved_rows)
observed_unresolved = sorted(unresolved["GSE"].unique().tolist())
expected_unresolved = ["GSE267072", "GSE285380"]
assert observed_unresolved == expected_unresolved, (
    "Unexpected reviewed-unresolved duration GSE",
    observed_unresolved,
    expected_unresolved,
)
unresolved.to_csv(UNRESOLVED_OUT, sep="\t", index=False)

# ------------------------------------------------------------
# Canonical Batch08 expected independent-GSE coverage QA
# ------------------------------------------------------------
c = cov.set_index("Moderator")

expected = {
    "hfd_duration_weeks": {"inclusive": 35, "strict": 33},
    "sex": {"strict": 37},
    "strain": {"strict": 40},
    "hfd_fat_kcal_percent": {"strict": 31},
    "control_diet_type": {"strict": 50},
}

for var, exp in expected.items():
    if "inclusive" in exp:
        got = int(c.loc[var, "Independent_GSE_inclusive_ready"])
        assert got == exp["inclusive"], (var, "inclusive", got, exp["inclusive"])
    if "strict" in exp:
        got = int(c.loc[var, "Independent_GSE_strict_ready"])
        assert got == exp["strict"], (var, "strict", got, exp["strict"])

# Explicit semantic QA for the two remaining strict gaps.
r614 = x[x["GSE"].eq("GSE285614")].iloc[0]
assert float(r614["hfd_duration_weeks_final"]) == 13.0446
assert r614["hfd_duration_weeks_final_confidence"] == "MEDIUM"
assert r614["hfd_duration_weeks_adjudication_status"] != "NOT_REVIEWED"

r324 = x[x["GSE"].eq("GSE324134")].iloc[0]
assert r324["hfd_duration_weeks_adjudication_status"] == "NOT_REVIEWED"
assert r324["hfd_duration_weeks_final"] == ""
assert float(r324["hfd_duration_weeks_geo_auto_value"]) == 12.0
assert r324["hfd_duration_weeks_geo_auto_confidence"] == "LOW"

print("===== HFD MODERATOR COVERAGE AFTER BATCH08 v1 =====")
print(cov.to_string(index=False))

print("\n===== DURATION INCLUSIVE-BUT-NOT-STRICT QUEUE =====")
print(queue.to_string(index=False))

print("\n===== DURATION REVIEWED-UNRESOLVED SOURCE GAPS =====")
print(unresolved[["GSE", "Contrast_id", "adjudication_status"]].to_string(index=False))

print("\nQA PASS")
print("- duration inclusive independent GSE = 35")
print("- duration strict independent GSE    = 33")
print("- sex strict independent GSE         = 37")
print("- strain strict independent GSE      = 40")
print("- HFD fat% strict independent GSE    = 31")
print("- control diet strict independent GSE= 50")
print("- strict-gap queue = GSE285614, GSE324134")
print("- reviewed-unresolved = GSE267072, GSE285380")

print("\nSaved:")
print(" ", COVERAGE_OUT.relative_to(ROOT))
print(" ", QUEUE_OUT.relative_to(ROOT))
print(" ", UNRESOLVED_OUT.relative_to(ROOT))
