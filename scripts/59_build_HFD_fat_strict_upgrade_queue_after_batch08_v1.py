#!/usr/bin/env python3
from pathlib import Path
import argparse
import re
import pandas as pd

# ============================================================
# Build HFD fat kcal% strict-upgrade queue after canonical Batch08 v1.
#
# READ-ONLY with respect to the master.
# No adjudication is performed here.
#
# Goal:
#   Identify GSEs whose upgrade can increase the number of
#   independent GSEs with strict-ready HFD fat kcal%.
#
# Important:
#   - Multiple contrasts within one GSE are dependent effects.
#   - If ANY contrast in a GSE is already strict-ready for fat%,
#     that GSE already contributes to independent-GSE strict coverage
#     and is excluded from this upgrade queue.
#   - GEO auto metadata are scouting clues, not ground truth.
#   - Ingredient wt% is never converted to kcal% fat.
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
GSM_AUTO = (
    ROOT /
    "00_metadata/contrasts/moderators/"
    "HFD_primary837_GSM_moderator_extraction_GEOauto_v1_1_1.tsv"
)
QUEUE_OUT = OUTDIR / "HFD_fat_next_strict_upgrade_queue_after_batch08_v1.tsv"
SCOUT_OUT = OUTDIR / "HFD_fat_source_scout_after_batch08_v1.tsv"

for p in [MASTER, GSM_AUTO]:
    if not p.exists():
        raise FileNotFoundError(p)

x = pd.read_csv(MASTER, sep="\t", dtype=str).fillna("")
g = pd.read_csv(GSM_AUTO, sep="\t", dtype=str).fillna("")

assert x["GSE"].nunique() == 64, x["GSE"].nunique()
assert x["Contrast_id"].nunique() == 105, x["Contrast_id"].nunique()
assert g["GSM"].nunique() == 837, g["GSM"].nunique()

VAR = "hfd_fat_kcal_percent"

required_master = [
    "GSE", "Contrast_id",
    f"{VAR}_geo_auto_value",
    f"{VAR}_geo_auto_status",
    f"{VAR}_geo_auto_raw",
    f"{VAR}_geo_auto_source",
    f"{VAR}_geo_auto_confidence",
    f"{VAR}_manual_value",
    f"{VAR}_manual_source_type",
    f"{VAR}_manual_source_citation",
    f"{VAR}_manual_evidence",
    f"{VAR}_manual_confidence",
    f"{VAR}_adjudication_status",
    f"{VAR}_final",
    f"{VAR}_final_confidence",
    f"{VAR}_final_source",
    f"{VAR}_final_note",
]
missing = [c for c in required_master if c not in x.columns]
if missing:
    raise RuntimeError(f"Master missing expected columns: {missing}")

for c in ["GSE", "Contrast_id", "GSM"]:
    if c not in g.columns:
        raise RuntimeError(f"GSM-auto file missing required column: {c}")

def effective(row):
    if row[f"{VAR}_adjudication_status"] != "NOT_REVIEWED":
        return (
            row[f"{VAR}_final"].strip(),
            row[f"{VAR}_final_confidence"].strip(),
            "MANUAL",
        )
    if row[f"{VAR}_geo_auto_status"] == "RESOLVED":
        return (
            row[f"{VAR}_geo_auto_value"].strip(),
            row[f"{VAR}_geo_auto_confidence"].strip(),
            "GEO_AUTO",
        )
    return "", "", "NONE"

def numeric_ready(val):
    if not val:
        return False
    try:
        z = float(val)
    except ValueError:
        return False
    return 0 < z <= 100

def strict_ready(row):
    val, conf, origin = effective(row)
    if not numeric_ready(val):
        return False
    if origin == "MANUAL":
        return conf == "HIGH"
    if origin == "GEO_AUTO":
        return conf in {"MEDIUM", "HIGH"}
    return False

# Canonical coverage anchor: Batch08 must currently have 31 independent
# GSE strict-ready for HFD fat kcal%.
strict_gses = set()
for _, row in x.iterrows():
    if strict_ready(row):
        strict_gses.add(row["GSE"])

assert len(strict_gses) == 31, (
    "Canonical Batch08 fat strict independent-GSE coverage mismatch",
    len(strict_gses), 31
)

# ------------------------------------------------------------
# Source-scout helpers
# ------------------------------------------------------------

def clean_join(vals, max_chars=5000):
    vals = [str(v).strip() for v in vals if str(v).strip()]
    s = " || ".join(vals)
    return s[:max_chars]

def selected_text(gse):
    z = g[g["GSE"].eq(gse)]
    vals = []
    for c in [
        "sample_title_raw",
        "sample_characteristics_raw",
        "sample_source_name_raw",
        "sample_description_raw",
    ]:
        if c in z.columns:
            vals.extend(z[c].tolist())
    return clean_join(vals)

def detect_diet_codes(text):
    # Conservative scout only. These tokens are review leads, not adjudicated facts.
    pats = [
        r"\bD\d{6,}[A-Z]?\b",        # e.g. D12492, D09100310
        r"\bRD\d{4,}[A-Z]?\b",       # e.g. RD12079B
        r"\bTD[.\-_ ]?\d{4,}\b",     # Teklad-style TD.88137
        r"\bF\d{4,}\b",              # generic formula-like token; scout only
    ]
    found = set()
    for p in pats:
        found.update(m.group(0) for m in re.finditer(p, text, flags=re.I))
    return sorted(found)

def detect_percent_fat_clues(text):
    # Capture explicit textual percentages near "fat".
    clues = set()
    patterns = [
        r"\b\d+(?:\.\d+)?\s*%\s*(?:kcal(?:ories)?\s*)?(?:from\s*)?fat\b",
        r"\bfat\s*[:=]?\s*\d+(?:\.\d+)?\s*%\b",
        r"\b\d+(?:\.\d+)?\s*kcal%\s*fat\b",
        r"\b\d+(?:\.\d+)?\s*%\s*of\s*(?:total\s*)?calories?\s*from\s*fat\b",
    ]
    for p in patterns:
        for m in re.finditer(p, text, flags=re.I):
            clues.add(m.group(0))
    return sorted(clues)

def detect_weight_percent_warning(text):
    # Presence of ingredient wt% wording means manual review must not convert it.
    return bool(
        re.search(
            r"\b(?:w/w|wt\.?\s*%|weight\s*%|by\s+weight|grams?\s+per\s+100\s*g)\b",
            text,
            flags=re.I,
        )
    )

def state_class(rows):
    # rows = all master contrasts for one GSE, none strict-ready by construction.
    eff = [effective(r) for _, r in rows.iterrows()]
    numeric = [(v,c,o) for v,c,o in eff if numeric_ready(v)]

    manual_unresolved = any(
        r[f"{VAR}_manual_confidence"] == "UNRESOLVED"
        for _, r in rows.iterrows()
    )
    pending = any(
        r[f"{VAR}_adjudication_status"] == "NOT_REVIEWED"
        for _, r in rows.iterrows()
    )

    if numeric:
        origins = {o for _,_,o in numeric}
        confs = {c for _,c,_ in numeric}
        if "MANUAL" in origins:
            return "NUMERIC_MANUAL_NONSTRICT"
        if "GEO_AUTO" in origins:
            return "NUMERIC_GEO_AUTO_LOW"
    if manual_unresolved:
        return "REVIEWED_UNRESOLVED"
    if pending:
        return "PENDING_MISSING"
    return "OTHER_NONSTRICT"

# ------------------------------------------------------------
# Build one row per independent GSE that is not already strict-ready.
# ------------------------------------------------------------

rows = []
scout_rows = []

for gse, z in x.groupby("GSE", sort=True):
    if gse in strict_gses:
        continue

    txt = selected_text(gse)
    diet_codes = detect_diet_codes(txt)
    pct_clues = detect_percent_fat_clues(txt)
    wt_warning = detect_weight_percent_warning(txt)

    effs = [effective(r) for _, r in z.iterrows()]
    eff_vals = sorted({v for v,_,_ in effs if v})
    eff_confs = sorted({c for _,c,_ in effs if c})
    eff_origins = sorted({o for _,_,o in effs if o != "NONE"})

    auto_vals = sorted({
        r[f"{VAR}_geo_auto_value"]
        for _, r in z.iterrows()
        if r[f"{VAR}_geo_auto_value"]
    })
    auto_confs = sorted({
        r[f"{VAR}_geo_auto_confidence"]
        for _, r in z.iterrows()
        if r[f"{VAR}_geo_auto_confidence"]
    })
    auto_raw = clean_join(
        [r[f"{VAR}_geo_auto_raw"] for _, r in z.iterrows()],
        max_chars=1800,
    )

    st = state_class(z)

    # Priority:
    # 1 = numeric but only LOW/MEDIUM confidence -> usually easiest strict upgrade
    # 2 = missing/unresolved but sample metadata exposes diet code + explicit kcal-fat clue
    # 3 = missing/unresolved with a diet code
    # 4 = missing/unresolved with explicit percent-fat wording only
    # 5 = no direct scouting clue; publication search likely required
    if st in {"NUMERIC_MANUAL_NONSTRICT", "NUMERIC_GEO_AUTO_LOW"}:
        priority = 1
        reason = "Numeric fat% exists but confidence is below strict threshold."
    elif diet_codes and pct_clues:
        priority = 2
        reason = "Selected GSM metadata contain diet-code and explicit fat-percentage clues."
    elif diet_codes:
        priority = 3
        reason = "Selected GSM metadata contain diet-code clue(s); verify manufacturer/publication composition."
    elif pct_clues:
        priority = 4
        reason = "Selected GSM metadata contain explicit fat-percentage wording; verify kcal% vs wt%."
    else:
        priority = 5
        reason = "No direct selected-GSM diet-code/fat-percentage clue; publication/supplement search likely needed."

    if wt_warning:
        reason += " Weight-percent wording detected: do not convert ingredient wt% to kcal%."

    # Other unresolved/pending moderators are useful for opportunistic co-adjudication.
    other_pending = []
    for other in ["sex", "strain", "hfd_duration_weeks", "control_diet_type"]:
        vals = []
        for _, rr in z.iterrows():
            if rr[f"{other}_adjudication_status"] != "NOT_REVIEWED":
                v = rr[f"{other}_final"].strip()
            elif rr[f"{other}_geo_auto_status"] == "RESOLVED":
                v = rr[f"{other}_geo_auto_value"].strip()
            else:
                v = ""
            vals.append(v)
        if not any(vals):
            other_pending.append(other)

    rows.append({
        "priority_tier": priority,
        "GSE": gse,
        "n_contrasts": z["Contrast_id"].nunique(),
        "fat_state": st,
        "effective_values": ",".join(eff_vals),
        "effective_confidences": ",".join(eff_confs),
        "effective_sources": ",".join(eff_origins),
        "geo_auto_values": ",".join(auto_vals),
        "geo_auto_confidences": ",".join(auto_confs),
        "diet_code_clues": ",".join(diet_codes),
        "explicit_fat_percent_clues": " || ".join(pct_clues),
        "weight_percent_warning": "YES" if wt_warning else "",
        "other_pending_moderators": ",".join(other_pending),
        "priority_reason": reason,
    })

    scout_rows.append({
        "GSE": gse,
        "n_contrasts": z["Contrast_id"].nunique(),
        "fat_state": st,
        "geo_auto_raw": auto_raw,
        "selected_GSM_text_preview": txt[:3500],
        "diet_code_clues": ",".join(diet_codes),
        "explicit_fat_percent_clues": " || ".join(pct_clues),
        "weight_percent_warning": "YES" if wt_warning else "",
    })

queue = pd.DataFrame(rows)

# There are 64 total GSE and 31 already strict-ready.
# Therefore exactly 33 independent GSE should remain in the queue.
assert len(queue) == 33, (
    "Unexpected number of independent GSE outside fat strict coverage",
    len(queue), 33
)

state_order = {
    "NUMERIC_MANUAL_NONSTRICT": 0,
    "NUMERIC_GEO_AUTO_LOW": 1,
    "REVIEWED_UNRESOLVED": 2,
    "PENDING_MISSING": 3,
    "OTHER_NONSTRICT": 4,
}
queue["_state_order"] = queue["fat_state"].map(state_order).fillna(9)
queue = queue.sort_values(
    ["priority_tier", "_state_order", "GSE"],
    kind="stable",
).drop(columns="_state_order").reset_index(drop=True)

queue.insert(0, "upgrade_order", range(1, len(queue) + 1))
queue.to_csv(QUEUE_OUT, sep="\t", index=False)

scout = pd.DataFrame(scout_rows).sort_values("GSE").reset_index(drop=True)
scout.to_csv(SCOUT_OUT, sep="\t", index=False)

print("===== HFD FAT KCAL% STRICT-UPGRADE QUEUE AFTER BATCH08 v1 =====")
print("Current strict independent GSE:", len(strict_gses))
print("Candidate independent GSE:", len(queue))
print()

show_cols = [
    "upgrade_order", "priority_tier", "GSE", "n_contrasts",
    "fat_state", "effective_values", "effective_confidences",
    "geo_auto_values", "geo_auto_confidences",
    "diet_code_clues", "explicit_fat_percent_clues",
    "weight_percent_warning", "other_pending_moderators",
]
print(queue[show_cols].head(15).to_string(index=False))

print("\n===== PRIORITY-TIER COUNTS =====")
print(
    queue.groupby(["priority_tier", "fat_state"])
    .size()
    .rename("n_GSE")
    .reset_index()
    .to_string(index=False)
)

print("\nQA PASS")
print("- master unchanged/read-only")
print("- canonical Batch08 HFD fat strict independent GSE = 31")
print("- independent GSE outside strict coverage = 33")
print("- GSEs already contributing strict fat% coverage are excluded")
print("- multiple contrasts remain grouped as dependent effects")
print("- wt% clues are flagged and never converted to kcal%")

print("\nSaved:")
print(" ", QUEUE_OUT.relative_to(ROOT))
print(" ", SCOUT_OUT.relative_to(ROOT))
