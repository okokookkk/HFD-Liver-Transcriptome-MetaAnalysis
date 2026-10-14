#!/usr/bin/env python3
from pathlib import Path
import argparse
import pandas as pd

# ============================================================
# Recompute moderator coverage after canonical Batch09 v1.
# Writes:
#   - coverage table
#   - next HFD fat strict-upgrade queue
#   - duration strict-gap queue (for continuity)
# ============================================================

ap = argparse.ArgumentParser()
ap.add_argument(
    "--root",
    default=__import__("os").environ.get("HFD_PROJECT_ROOT", "."),
)
args = ap.parse_args()

ROOT = Path(args.root)
OUTDIR = ROOT / "00_metadata/contrasts/moderators/manual_adjudication"

MASTER = OUTDIR / "HFD_primary105_moderator_manual_adjudication_batch09_v1.tsv"
COVERAGE_OUT = OUTDIR / "HFD_primary105_moderator_coverage_after_batch09_v1.tsv"
FAT_QUEUE_OUT = OUTDIR / "HFD_fat_next_strict_upgrade_queue_after_batch09_v1.tsv"
DURATION_QUEUE_OUT = OUTDIR / "HFD_duration_next_strict_upgrade_queue_after_batch09_v1.tsv"

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
        return val in VALID_SEX
    if var == "strain":
        return val in VALID_STRAINS
    if var == "control_diet_type":
        return val in VALID_CONTROL
    return False

def strict_ready(var, val, conf, origin):
    if not inclusive_ready(var, val):
        return False
    if origin == "MANUAL":
        return conf == "HIGH"
    if origin == "GEO_AUTO":
        return conf in {"MEDIUM", "HIGH"}
    return False

coverage = []
for var in MODS:
    known = inclusive = strict = unresolved = pending = 0
    inclusive_gse, strict_gse = set(), set()

    for _, row in x.iterrows():
        val, conf, origin = effective(row, var)
        inc = inclusive_ready(var, val)
        st = strict_ready(var, val, conf, origin)

        if val:
            known += 1
        if inc:
            inclusive += 1
            inclusive_gse.add(row["GSE"])
        if st:
            strict += 1
            strict_gse.add(row["GSE"])
        if row[f"{var}_manual_confidence"] == "UNRESOLVED":
            unresolved += 1
        if row[f"{var}_adjudication_status"] == "NOT_REVIEWED" and not inc:
            pending += 1

    coverage.append({
        "Moderator": var,
        "Total_contrasts": len(x),
        "Known_contrasts": known,
        "Inclusive_ready_contrasts": inclusive,
        "Independent_GSE_inclusive_ready": len(inclusive_gse),
        "Strict_ready_contrasts": strict,
        "Independent_GSE_strict_ready": len(strict_gse),
        "Reviewed_but_unresolved_contrasts": unresolved,
        "Pending_review_contrasts": pending,
    })

cov = pd.DataFrame(coverage)
cov.to_csv(COVERAGE_OUT, sep="\t", index=False)
c = cov.set_index("Moderator")

# Expected Batch09 independent-GSE coverage.
expected = {
    "hfd_duration_weeks": {"inclusive": 40, "strict": 38},
    "sex": {"strict": 40},
    "strain": {"strict": 41},
    "hfd_fat_kcal_percent": {"inclusive": 37, "strict": 36},
    "control_diet_type": {"strict": 50},
}
for var, exp in expected.items():
    if "inclusive" in exp:
        got = int(c.loc[var, "Independent_GSE_inclusive_ready"])
        assert got == exp["inclusive"], (var, "inclusive", got, exp["inclusive"])
    if "strict" in exp:
        got = int(c.loc[var, "Independent_GSE_strict_ready"])
        assert got == exp["strict"], (var, "strict", got, exp["strict"])

# Target semantic QA.
checks = {
    "GSE124777": {"hfd_fat_kcal_percent": "45", "hfd_duration_weeks": "16", "sex": "MALE"},
    "GSE137022": {"hfd_fat_kcal_percent": "45", "hfd_duration_weeks": "14", "strain": "C57BL/6"},
    "GSE145665": {"hfd_fat_kcal_percent": "60", "hfd_duration_weeks": "24", "sex": "MALE"},
    "GSE146627": {"hfd_fat_kcal_percent": "45", "hfd_duration_weeks": "8"},
    "GSE193785": {"hfd_fat_kcal_percent": "60", "hfd_duration_weeks": "10", "sex": "MALE"},
}
for gse, vals in checks.items():
    z = x[x["GSE"].eq(gse)]
    assert len(z) == 1
    r = z.iloc[0]
    for var, val in vals.items():
        assert r[f"{var}_final"] == val, (gse, var, r[f"{var}_final"], val)
        assert r[f"{var}_final_confidence"] == "HIGH", (gse, var)

# GSE145840 remains the sole numeric manual non-strict fat candidate expected
# from the pre-Batch09 queue.
z145840 = x[x["GSE"].eq("GSE145840")]
assert len(z145840) == 1
r145840 = z145840.iloc[0]
assert float(r145840["hfd_fat_kcal_percent_final"]) == 60.0
assert r145840["hfd_fat_kcal_percent_final_confidence"] == "MEDIUM"

# ------------------------------------------------------------
# Next fat strict-upgrade queue: one row per independent GSE not
# already contributing strict fat coverage.
# ------------------------------------------------------------
fat_strict_gses = set()
for _, row in x.iterrows():
    val, conf, origin = effective(row, "hfd_fat_kcal_percent")
    if strict_ready("hfd_fat_kcal_percent", val, conf, origin):
        fat_strict_gses.add(row["GSE"])

fat_queue_rows = []
for gse, z in x.groupby("GSE", sort=True):
    if gse in fat_strict_gses:
        continue

    eff = [effective(r, "hfd_fat_kcal_percent") for _, r in z.iterrows()]
    vals = sorted({v for v,_,_ in eff if v})
    confs = sorted({q for _,q,_ in eff if q})
    origins = sorted({o for _,_,o in eff if o != "NONE"})

    if any(o == "MANUAL" and q in {"LOW", "MEDIUM"} and v for v,q,o in eff):
        state = "NUMERIC_MANUAL_NONSTRICT"
        priority = 1
    elif any(r["hfd_fat_kcal_percent_manual_confidence"] == "UNRESOLVED"
             for _, r in z.iterrows()):
        state = "REVIEWED_UNRESOLVED"
        priority = 2
    else:
        state = "PENDING_MISSING"
        priority = 3

    fat_queue_rows.append({
        "priority_tier": priority,
        "GSE": gse,
        "n_contrasts": z["Contrast_id"].nunique(),
        "fat_state": state,
        "effective_values": ",".join(vals),
        "effective_confidences": ",".join(confs),
        "effective_sources": ",".join(origins),
    })

fatq = pd.DataFrame(fat_queue_rows).sort_values(
    ["priority_tier", "GSE"], kind="stable"
).reset_index(drop=True)
fatq.insert(0, "upgrade_order", range(1, len(fatq) + 1))

# 64 total - 36 fat strict GSE = 28 remaining.
assert len(fatq) == 28, len(fatq)
assert fatq.iloc[0]["GSE"] == "GSE145840"
assert fatq.iloc[0]["fat_state"] == "NUMERIC_MANUAL_NONSTRICT"
assert float(fatq.iloc[0]["effective_values"]) == 60.0
assert fatq.iloc[0]["effective_confidences"] == "MEDIUM"
fatq.to_csv(FAT_QUEUE_OUT, sep="\t", index=False)

# ------------------------------------------------------------
# Duration inclusive-but-not-strict queue should remain the same
# two GSE as after Batch08.
# ------------------------------------------------------------
dur_rows = []
for gse, z in x.groupby("GSE", sort=True):
    gaps = []
    for _, row in z.iterrows():
        val, conf, origin = effective(row, "hfd_duration_weeks")
        if (
            inclusive_ready("hfd_duration_weeks", val)
            and not strict_ready("hfd_duration_weeks", val, conf, origin)
        ):
            gaps.append((val, conf, origin))
    if gaps:
        dur_rows.append({
            "GSE": gse,
            "n_gap_contrasts": len(gaps),
            "duration_values": ",".join(sorted({v for v,_,_ in gaps}, key=float)),
            "confidence": ",".join(sorted({q for _,q,_ in gaps})),
            "source_level": ",".join(sorted({o for _,_,o in gaps})),
        })

durq = pd.DataFrame(dur_rows)
observed = sorted(durq["GSE"].tolist())
assert observed == ["GSE285614", "GSE324134"], observed
order = {"GSE285614": 1, "GSE324134": 2}
durq.insert(0, "Next_duration_upgrade_order", durq["GSE"].map(order))
durq = durq.sort_values("Next_duration_upgrade_order")
durq.to_csv(DURATION_QUEUE_OUT, sep="\t", index=False)

print("===== HFD MODERATOR COVERAGE AFTER BATCH09 v1 =====")
print(cov.to_string(index=False))

print("\n===== NEXT HFD FAT STRICT-UPGRADE QUEUE (TOP 15) =====")
print(fatq.head(15).to_string(index=False))

print("\n===== DURATION INCLUSIVE-BUT-NOT-STRICT =====")
print(durq.to_string(index=False))

print("\nQA PASS")
print("- HFD duration inclusive/strict GSE = 40 / 38")
print("- sex strict GSE                  = 40")
print("- strain strict GSE               = 41")
print("- HFD fat strict GSE              = 36")
print("- control diet strict GSE         = 50")
print("- remaining fat non-strict GSE    = 28")
print("- GSE145840 remains 60 MEDIUM")
print("- duration strict gaps remain GSE285614 + GSE324134")

print("\nSaved:")
print(" ", COVERAGE_OUT.relative_to(ROOT))
print(" ", FAT_QUEUE_OUT.relative_to(ROOT))
print(" ", DURATION_QUEUE_OUT.relative_to(ROOT))
