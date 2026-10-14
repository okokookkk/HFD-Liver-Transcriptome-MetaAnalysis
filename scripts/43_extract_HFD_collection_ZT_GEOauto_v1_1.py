#!/usr/bin/env python3
from pathlib import Path
import argparse, re
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument(
    "--root",
    default=__import__("os").environ.get("HFD_PROJECT_ROOT", ".")
)
args = ap.parse_args()

ROOT = Path(args.root)

GSM_AUTO = (
    ROOT /
    "00_metadata/contrasts/moderators/"
    "HFD_primary837_GSM_moderator_extraction_GEOauto_v1_1_1.tsv"
)

BATCH01 = (
    ROOT /
    "00_metadata/contrasts/moderators/manual_adjudication/"
    "HFD_primary105_moderator_manual_adjudication_batch01_v1.tsv"
)

OUTDIR = ROOT / "00_metadata/contrasts/moderators/time_of_day"
OUTDIR.mkdir(parents=True, exist_ok=True)

if not GSM_AUTO.exists():
    raise FileNotFoundError(GSM_AUTO)
if not BATCH01.exists():
    raise FileNotFoundError(BATCH01)

g = pd.read_csv(GSM_AUTO, sep="\t", dtype=str).fillna("")
c = pd.read_csv(BATCH01, sep="\t", dtype=str).fillna("")

assert g["GSM"].nunique() == 837
assert g["Contrast_id"].nunique() == 105
assert c["Contrast_id"].nunique() == 105

def norm_space(x):
    return re.sub(r"\s+", " ", str(x or "")).strip()

# Regression guard: underscore-prefixed ZT tokens must be detectable.
assert re.search(
    r"(?<![A-Za-z0-9])ZT\s*0*(\d{1,2})(?:\.\d+)?\b",
    "RNA-Seq_ZT00 LFD WT REP 1",
    re.I,
).group(1) == "0"

def extract_zt(row):
    """
    Conservative collection-time extraction.
    Only explicit Zeitgeber Time / ZT tokens are accepted automatically.

    We do NOT infer ZT from arbitrary clock time, fasting duration,
    treatment timing, or series-level text.
    """
    texts = [
        row.get("sample_title_raw", ""),
        row.get("sample_characteristics_raw", ""),
    ]
    text = " | ".join(t for t in texts if t)

    hits = []

    # Examples: ZT0, ZT00, ZT 6, ZT12, zeitgeber time: ZT18
    for m in re.finditer(r"(?<![A-Za-z0-9])ZT\s*0*(\d{1,2})(?:\.\d+)?\b", text, re.I):
        val = int(m.group(1))
        if 0 <= val <= 24:
            # Normalize ZT24 to ZT0 only as a clock-cycle equivalent,
            # but preserve the raw string.
            norm = 0 if val == 24 else val
            hits.append((norm, m.group(0), "EXPLICIT_ZT_TOKEN"))

    vals = sorted(set(x[0] for x in hits))

    if len(vals) == 1:
        return {
            "collection_ZT_value": f"ZT{vals[0]:02d}",
            "collection_ZT_hour": str(vals[0]),
            "collection_ZT_status": "RESOLVED",
            "collection_ZT_raw": " | ".join(sorted(set(x[1] for x in hits))),
            "collection_ZT_source": "GEO_GSM_SAMPLE_TEXT",
            "collection_ZT_confidence": "HIGH",
        }

    if len(vals) > 1:
        return {
            "collection_ZT_value": "CONFLICT",
            "collection_ZT_hour": "",
            "collection_ZT_status": "CONFLICT",
            "collection_ZT_raw": " | ".join(sorted(set(x[1] for x in hits))),
            "collection_ZT_source": "GEO_GSM_SAMPLE_TEXT",
            "collection_ZT_confidence": "HIGH",
        }

    return {
        "collection_ZT_value": "UNKNOWN",
        "collection_ZT_hour": "",
        "collection_ZT_status": "MISSING",
        "collection_ZT_raw": "",
        "collection_ZT_source": "",
        "collection_ZT_confidence": "",
    }

rows = []
for _, r in g.iterrows():
    z = {
        "GSE": r["GSE"],
        "Contrast_id": r["Contrast_id"],
        "GSM": r["GSM"],
        "Role": r["Role"],
    }
    z.update(extract_zt(r))
    rows.append(z)

sg = pd.DataFrame(rows)

# GSM-level QA
assert len(sg) == 837
assert sg["GSM"].nunique() == 837
assert sg["Contrast_id"].nunique() == 105

# Contrast-level collapse.
contrast_rows = []
for cid, z in sg.groupby("Contrast_id", sort=True):
    vals = sorted(set(
        v for v in z["collection_ZT_value"]
        if v not in {"", "UNKNOWN"}
    ))

    if len(vals) == 0:
        status = "MISSING"
        value = "UNKNOWN"
        hour = ""
        confidence = ""
    elif "CONFLICT" in vals or len(vals) > 1:
        status = "CONFLICT"
        value = "CONFLICT"
        hour = ""
        confidence = "HIGH"
    else:
        status = "RESOLVED"
        value = vals[0]
        hour = z.loc[
            z["collection_ZT_value"] == value,
            "collection_ZT_hour"
        ].iloc[0]
        confidence = "HIGH"

    # A valid contrast should have the same collection ZT in HFD/control.
    control_vals = sorted(set(
        v for v in z.loc[z["Role"] == "CONTROL", "collection_ZT_value"]
        if v not in {"", "UNKNOWN"}
    ))
    hfd_vals = sorted(set(
        v for v in z.loc[z["Role"] == "HFD", "collection_ZT_value"]
        if v not in {"", "UNKNOWN"}
    ))

    if (
        status == "RESOLVED"
        and control_vals
        and hfd_vals
        and control_vals != hfd_vals
    ):
        status = "ROLE_DISCORDANT"
        value = "CONFLICT"
        hour = ""

    contrast_rows.append({
        "GSE": z["GSE"].iloc[0],
        "Contrast_id": cid,
        "collection_ZT_value": value,
        "collection_ZT_hour": hour,
        "collection_ZT_status": status,
        "collection_ZT_raw": " || ".join(
            sorted(set(v for v in z["collection_ZT_raw"] if v))
        ),
        "collection_ZT_source": (
            "GEO_GSM_SAMPLE_TEXT"
            if status != "MISSING" else ""
        ),
        "collection_ZT_confidence": confidence,
        "n_GSM": z["GSM"].nunique(),
        "n_CONTROL_GSM": z.loc[z["Role"] == "CONTROL", "GSM"].nunique(),
        "n_HFD_GSM": z.loc[z["Role"] == "HFD", "GSM"].nunique(),
    })

ct = pd.DataFrame(contrast_rows).sort_values(["GSE", "Contrast_id"])

assert ct["Contrast_id"].nunique() == 105
assert ct["GSE"].nunique() == 64

# GSE-level coverage / distribution.
gse_summary = (
    ct.groupby("GSE")
      .agg(
          n_contrasts=("Contrast_id", "nunique"),
          resolved_ZT_contrasts=(
              "collection_ZT_status",
              lambda s: int((s == "RESOLVED").sum())
          ),
          missing_ZT_contrasts=(
              "collection_ZT_status",
              lambda s: int((s == "MISSING").sum())
          ),
          conflict_ZT_contrasts=(
              "collection_ZT_status",
              lambda s: int((s.isin(["CONFLICT", "ROLE_DISCORDANT"])).sum())
          ),
      )
      .reset_index()
)

# Distribution by independent GSE.
resolved = ct[ct["collection_ZT_status"] == "RESOLVED"].copy()
if len(resolved):
    dist = (
        resolved.groupby("collection_ZT_value")
        .agg(
            contrasts=("Contrast_id", "nunique"),
            GSE=("GSE", "nunique"),
        )
        .reset_index()
        .sort_values(["GSE", "contrasts"], ascending=False)
    )
else:
    dist = pd.DataFrame(
        columns=["collection_ZT_value", "contrasts", "GSE"]
    )

# Merge ZT auto columns onto current manual adjudication master.
m = c.merge(
    ct[
        [
            "Contrast_id",
            "collection_ZT_value",
            "collection_ZT_hour",
            "collection_ZT_status",
            "collection_ZT_raw",
            "collection_ZT_source",
            "collection_ZT_confidence",
        ]
    ],
    on="Contrast_id",
    how="left",
    validate="one_to_one",
)

# Add manual/final ZT fields for later paper-level adjudication.
m["collection_ZT_manual_value"] = ""
m["collection_ZT_manual_source_type"] = ""
m["collection_ZT_manual_source_citation"] = ""
m["collection_ZT_manual_evidence"] = ""
m["collection_ZT_manual_confidence"] = ""
m["collection_ZT_adjudication_status"] = "NOT_REVIEWED"
m["collection_ZT_final"] = ""
m["collection_ZT_final_confidence"] = ""
m["collection_ZT_final_source"] = ""
m["collection_ZT_final_note"] = ""

# For explicit resolved GSM-level ZT, prefill final as auto-confirmed,
# but leave adjudication status as AUTO_GEO_RESOLVED so later paper
# review can override if necessary.
mask = m["collection_ZT_status"] == "RESOLVED"
m.loc[mask, "collection_ZT_final"] = m.loc[mask, "collection_ZT_value"]
m.loc[mask, "collection_ZT_final_confidence"] = "HIGH"
m.loc[mask, "collection_ZT_final_source"] = "GEO_GSM"
m.loc[mask, "collection_ZT_adjudication_status"] = "AUTO_GEO_RESOLVED"

# Save.
gsm_out = OUTDIR / "HFD_primary837_GSM_collection_ZT_GEOauto_v1_1.tsv"
contrast_out = OUTDIR / "HFD_primary105_contrast_collection_ZT_GEOauto_v1_1.tsv"
gse_out = OUTDIR / "HFD_primary64_collection_ZT_coverage_by_GSE_v1_1.tsv"
dist_out = OUTDIR / "HFD_primary105_collection_ZT_distribution_v1_1.tsv"
master_out = (
    ROOT /
    "00_metadata/contrasts/moderators/manual_adjudication/"
    "HFD_primary105_moderator_manual_adjudication_batch01_plus_ZT_v1_1.tsv"
)

sg.to_csv(gsm_out, sep="\t", index=False)
ct.to_csv(contrast_out, sep="\t", index=False)
gse_summary.to_csv(gse_out, sep="\t", index=False)
dist.to_csv(dist_out, sep="\t", index=False)
m.to_csv(master_out, sep="\t", index=False)

print("===== HFD COLLECTION ZT GEO-AUTO v1.1 PASS =====")
print("GSE:", ct["GSE"].nunique())
print("Contrasts:", ct["Contrast_id"].nunique())
print("GSM:", sg["GSM"].nunique())

print("\n===== CONTRAST ZT STATUS =====")
print(ct["collection_ZT_status"].value_counts().to_string())

print("\n===== RESOLVED ZT DISTRIBUTION =====")
if len(dist):
    print(dist.to_string(index=False))
else:
    print("No resolved ZT contrasts")

print("\n===== GSE WITH MOST RESOLVED ZT CONTRASTS =====")
print(
    gse_summary.sort_values(
        ["resolved_ZT_contrasts", "n_contrasts"],
        ascending=False
    ).head(20).to_string(index=False)
)

print("\nSaved:")
for p in [gsm_out, contrast_out, gse_out, dist_out, master_out]:
    print(p.relative_to(ROOT))

print("\nRULE:")
print("Only explicit GSM-level ZT tokens are auto-resolved.")
print("Clock times, age, intervention timing, and fasting duration are not inferred as ZT.")
