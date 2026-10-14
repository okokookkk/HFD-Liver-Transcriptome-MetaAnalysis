#!/usr/bin/env python3
from pathlib import Path
import argparse, gzip, re
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument(
    "--root",
    default=__import__("os").environ.get("HFD_PROJECT_ROOT", ".")
)
args = ap.parse_args()

ROOT = Path(args.root)
MANIFEST = (
    ROOT /
    "00_metadata/technical_harmonization/frozen_v1/"
    "HFD_primary_ENA_R1R2_download_manifest_frozen_v1.tsv"
)
SOFTDIR = ROOT / "00_metadata/geo_soft"
OUT = ROOT / "00_metadata/contrasts/moderators"
OUT.mkdir(parents=True, exist_ok=True)

def norm_space(x):
    return re.sub(r"\s+", " ", str(x or "")).strip()

def unique_nonempty(values):
    out = []
    seen = set()
    for v in values:
        v = norm_space(v)
        if not v:
            continue
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out

def characteristic_pairs(chars):
    out = []
    for raw in chars:
        raw = norm_space(raw)
        if ":" in raw:
            k, v = raw.split(":", 1)
            out.append((norm_space(k).lower(), norm_space(v), raw))
    return out

def parse_soft(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    samples = {}
    series = {"title": [], "summary": [], "overall_design": []}
    current = None
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                current = norm_space(line.split("=", 1)[1])
                samples.setdefault(current, {
                    "title": [],
                    "source_name_ch1": [],
                    "characteristics_ch1": [],
                    "extract_protocol_ch1": [],
                    "library_strategy": [],
                    "library_source": [],
                    "library_selection": [],
                })
                continue
            if line.startswith("!Series_title = "):
                series["title"].append(norm_space(line.split("=", 1)[1]))
                continue
            if line.startswith("!Series_summary = "):
                series["summary"].append(norm_space(line.split("=", 1)[1]))
                continue
            if line.startswith("!Series_overall_design = "):
                series["overall_design"].append(norm_space(line.split("=", 1)[1]))
                continue
            if current is None:
                continue
            prefix_map = {
                "!Sample_title = ": "title",
                "!Sample_source_name_ch1 = ": "source_name_ch1",
                "!Sample_characteristics_ch1 = ": "characteristics_ch1",
                "!Sample_extract_protocol_ch1 = ": "extract_protocol_ch1",
                "!Sample_library_strategy = ": "library_strategy",
                "!Sample_library_source = ": "library_source",
                "!Sample_library_selection = ": "library_selection",
            }
            for prefix, key in prefix_map.items():
                if line.startswith(prefix):
                    samples[current][key].append(
                        norm_space(line.split("=", 1)[1])
                    )
                    break
    return samples, series

def make_result(value="", raw="", source="", confidence="", note=""):
    return {
        "value": norm_space(value),
        "raw": norm_space(raw),
        "source": norm_space(source),
        "confidence": norm_space(confidence),
        "note": norm_space(note),
    }

def extract_sex(sample):
    pairs = characteristic_pairs(sample["characteristics_ch1"])
    vals, raws = [], []
    for k, v, raw in pairs:
        if k in {"sex", "gender"}:
            s = v.lower()
            if re.search(r"\bmale\b", s):
                vals.append("MALE"); raws.append(raw)
            elif re.search(r"\bfemale\b", s):
                vals.append("FEMALE"); raws.append(raw)
    vals = unique_nonempty(vals)
    if len(vals) == 1:
        return make_result(vals[0], " | ".join(raws),
                           "GEO_SOFT_SAMPLE_CHARACTERISTIC_KEYED", "HIGH")
    if len(vals) > 1:
        return make_result("CONFLICT", " | ".join(raws),
                           "GEO_SOFT_SAMPLE_CHARACTERISTIC_KEYED", "HIGH",
                           "Multiple explicit sex values within GSM")
    text = " | ".join(sample["title"] + sample["characteristics_ch1"])
    hits = []
    if re.search(r"\bmale\b", text, re.I):
        hits.append("MALE")
    if re.search(r"\bfemale\b", text, re.I):
        hits.append("FEMALE")
    if len(hits) == 1:
        return make_result(hits[0], text, "GEO_SOFT_SAMPLE_TEXT", "MEDIUM")
    if len(hits) > 1:
        return make_result("CONFLICT", text, "GEO_SOFT_SAMPLE_TEXT", "LOW")
    return make_result()

def normalize_strain(v):
    s = norm_space(v)
    low = s.lower().replace(" ", "")
    patterns = [
        (r"c57bl/?6j", "C57BL/6J"),
        (r"c57bl/?6n", "C57BL/6N"),
        (r"c57bl/?6", "C57BL/6"),
        (r"balb/?c", "BALB/c"),
        (r"129s\d*", "129"),
        (r"db/?db", "db/db"),
    ]
    for pat, label in patterns:
        if re.search(pat, low, re.I):
            return label
    return s

def extract_strain(sample):
    pairs = characteristic_pairs(sample["characteristics_ch1"])
    vals, raws = [], []
    for k, v, raw in pairs:
        if "strain" in k or "background" in k:
            vals.append(normalize_strain(v))
            raws.append(raw)
    vals = unique_nonempty(vals)
    if len(vals) == 1:
        return make_result(vals[0], " | ".join(raws),
                           "GEO_SOFT_SAMPLE_CHARACTERISTIC_KEYED", "HIGH")
    if len(vals) > 1:
        return make_result("CONFLICT", " | ".join(raws),
                           "GEO_SOFT_SAMPLE_CHARACTERISTIC_KEYED", "HIGH")
    return make_result()

def parse_weeks(text, key=""):
    """
    Parse an explicit duration/age into weeks.

    Handles:
      value='16 weeks'
      value='112 days'
      key='weeks on diet', value='16'
      key='age (weeks)', value='12'

    It deliberately does NOT convert months automatically.
    """
    t = norm_space(text).lower()
    k = norm_space(key).lower()

    m = re.search(
        r"(?<![\d.])(\d+(?:\.\d+)?)\s*(?:week|weeks|wk|wks)\b",
        t
    )
    if m:
        return float(m.group(1)), "WEEK_EXPLICIT_VALUE"

    m = re.search(
        r"(?<![\d.])(\d+(?:\.\d+)?)\s*(?:day|days)\b",
        t
    )
    if m:
        return float(m.group(1)) / 7.0, "DAY_EXPLICIT_VALUE_CONVERTED_TO_WEEK"

    # Numeric-only value when the unit is explicitly encoded in the key.
    m = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*", t)
    if m:
        value = float(m.group(1))
        if re.search(r"\b(?:week|weeks|wk|wks)\b", k):
            return value, "WEEK_FROM_KEY_UNIT"
        if re.search(r"\b(?:day|days)\b", k):
            return value / 7.0, "DAY_FROM_KEY_UNIT_CONVERTED_TO_WEEK"

    return None, ""

def extract_age(sample):
    pairs = characteristic_pairs(sample["characteristics_ch1"])
    candidates = []
    for k, v, raw in pairs:
        if "age" in k:
            w, how = parse_weeks(v, k)
            if w is not None:
                candidates.append((w, raw, how))
    vals = sorted(set(round(x[0], 6) for x in candidates))
    if len(vals) == 1:
        return make_result(str(vals[0]), " | ".join(x[1] for x in candidates),
                           "GEO_SOFT_SAMPLE_CHARACTERISTIC_KEYED", "HIGH")
    if len(vals) > 1:
        return make_result("CONFLICT", " | ".join(x[1] for x in candidates),
                           "GEO_SOFT_SAMPLE_CHARACTERISTIC_KEYED", "HIGH")
    return make_result()

def extract_duration(sample):
    """
    Extract HFD/diet exposure duration from GSM-level metadata.

    Priority:
      1) explicit duration characteristic keys
      2) sample title/characteristics phrases tied to diet/HFD

    Multiple distinct durations within one GSM are CONFLICT, never guessed.
    """
    pairs = characteristic_pairs(sample["characteristics_ch1"])
    keys = (
        "weeks on diet", "week on diet", "wks on diet",
        "days on diet", "day on diet",
        "diet duration", "duration on diet",
        "feeding duration", "hfd duration",
        "high fat diet duration", "high-fat diet duration",
        "treatment duration", "duration",
    )

    candidates = []

    for k, v, raw in pairs:
        if any(sk == k or sk in k for sk in keys):
            w, how = parse_weeks(v, k)
            if w is not None:
                candidates.append((w, raw, how, "KEYED"))

    # Conservative text fallback: require diet/HFD wording near duration.
    text_chunks = sample["title"] + sample["characteristics_ch1"]
    pats = [
        r"(?:hfd|high[- ]fat diet|high fat diet|diet)[^|;]{0,35}?"
        r"(\d+(?:\.\d+)?)\s*(week|weeks|wk|wks)\b",
        r"(\d+(?:\.\d+)?)\s*(week|weeks|wk|wks)\b"
        r"[^|;]{0,35}?(?:hfd|high[- ]fat diet|high fat diet|diet)",
        r"(?:hfd|high[- ]fat diet|high fat diet|diet)[^|;]{0,35}?"
        r"(\d+(?:\.\d+)?)\s*(day|days)\b",
        r"(\d+(?:\.\d+)?)\s*(day|days)\b"
        r"[^|;]{0,35}?(?:hfd|high[- ]fat diet|high fat diet|diet)",
    ]

    for chunk in text_chunks:
        low = norm_space(chunk)
        for pat in pats:
            for m in re.finditer(pat, low, re.I):
                number = float(m.group(1))
                unit = m.group(2).lower()
                w = number if unit.startswith(("week", "wk")) else number / 7.0
                candidates.append((w, m.group(0), "TEXT_DIET_CONTEXT", "TEXT"))

    vals = sorted(set(round(x[0], 6) for x in candidates))

    if len(vals) == 1:
        keyed = any(x[3] == "KEYED" for x in candidates)
        return make_result(
            str(vals[0]),
            " | ".join(unique_nonempty(x[1] for x in candidates)),
            "GEO_SOFT_SAMPLE_CHARACTERISTIC_KEYED"
            if keyed else "GEO_SOFT_SAMPLE_DIET_CONTEXT_TEXT",
            "HIGH" if keyed else "MEDIUM"
        )

    if len(vals) > 1:
        return make_result(
            "CONFLICT",
            " | ".join(unique_nonempty(x[1] for x in candidates)),
            "GEO_SOFT_SAMPLE_TEXT",
            "HIGH",
            "Multiple distinct diet-duration values within GSM"
        )

    return make_result()

def extract_series_duration(series):
    """
    Series-level fallback used only when exactly one explicit diet/HFD duration
    exists across Series title/summary/overall_design.
    """
    candidates = []
    chunks = series["title"] + series["summary"] + series["overall_design"]

    pats = [
        r"(?:hfd|high[- ]fat diet|high fat diet|diet)[^.;|]{0,50}?"
        r"(\d+(?:\.\d+)?)\s*(week|weeks|wk|wks)\b",
        r"(\d+(?:\.\d+)?)\s*(week|weeks|wk|wks)\b"
        r"[^.;|]{0,50}?(?:hfd|high[- ]fat diet|high fat diet|diet)",
        r"(?:hfd|high[- ]fat diet|high fat diet|diet)[^.;|]{0,50}?"
        r"(\d+(?:\.\d+)?)\s*(day|days)\b",
        r"(\d+(?:\.\d+)?)\s*(day|days)\b"
        r"[^.;|]{0,50}?(?:hfd|high[- ]fat diet|high fat diet|diet)",
    ]

    for chunk in chunks:
        for pat in pats:
            for m in re.finditer(pat, norm_space(chunk), re.I):
                number = float(m.group(1))
                unit = m.group(2).lower()
                w = number if unit.startswith(("week", "wk")) else number / 7.0
                candidates.append((round(w, 6), m.group(0)))

    vals = sorted(set(x[0] for x in candidates))
    if len(vals) == 1:
        return make_result(
            str(vals[0]),
            " | ".join(unique_nonempty(x[1] for x in candidates)),
            "GEO_SOFT_SERIES_DIET_CONTEXT_TEXT",
            "LOW",
            "Series-level fallback; verify manually if used in analysis"
        )
    if len(vals) > 1:
        return make_result(
            "CONFLICT",
            " | ".join(unique_nonempty(x[1] for x in candidates)),
            "GEO_SOFT_SERIES_DIET_CONTEXT_TEXT",
            "LOW",
            "Multiple diet durations detected at series level"
        )
    return make_result()

def extract_fasting(sample):
    pairs = characteristic_pairs(sample["characteristics_ch1"])
    keyed = []
    for k, v, raw in pairs:
        if any(word in k for word in ["fast", "feed", "fed", "feeding", "diet", "treatment"]):
            keyed.append(raw)
    text = " | ".join(keyed + sample["title"])
    low = text.lower()
    fasted = bool(re.search(r"\bfasted\b|\bfasting\b|\b\d+\s*h(?:our)?s?\s*fast", low))
    fed = bool(re.search(r"\bad\s*lib(?:itum)?\b|\bfed\b|\bfeeding\b", low))
    if fasted and not fed:
        return make_result("FASTED", text, "GEO_SOFT_SAMPLE_TEXT", "HIGH")
    if fed and not fasted:
        return make_result("FED_AD_LIB_OR_FED", text, "GEO_SOFT_SAMPLE_TEXT", "HIGH")
    if fasted and fed:
        return make_result("CONFLICT", text, "GEO_SOFT_SAMPLE_TEXT", "MEDIUM",
                           "Both fasting and fed terms detected in GSM text")
    return make_result()

def extract_control_diet(sample):
    text = " | ".join(sample["title"] + sample["characteristics_ch1"])
    low = text.lower()
    hits = []
    if re.search(r"\bchow\b", low):
        hits.append("CHOW")
    if re.search(r"\blfd\b|low[- ]fat", low):
        hits.append("LFD")
    if re.search(r"\bnormal diet\b|\bnd\b", low):
        hits.append("NORMAL_DIET")
    if re.search(r"\bcontrol diet\b", low):
        hits.append("CONTROL_DIET_UNSPECIFIED")
    hits = unique_nonempty(hits)
    if len(hits) == 1:
        return make_result(hits[0], text, "GEO_SOFT_SAMPLE_TEXT", "MEDIUM")
    if len(hits) > 1:
        return make_result("|".join(hits), text, "GEO_SOFT_SAMPLE_TEXT", "LOW",
                           "Multiple control-diet descriptors detected")
    return make_result()

def fat_percent_from_text(text):
    """
    Extract dietary fat energy percentage conservatively.
    """
    t = norm_space(text)

    pats_high = [
        r"(\d+(?:\.\d+)?)\s*%\s*(?:of\s*)?(?:kcal|calories?|cal\.?|energy)"
        r"\s*(?:from\s*)?fat",
        r"(\d+(?:\.\d+)?)\s*(?:kcal|cal)\s*%\s*(?:from\s*)?fat",
        r"fat\s*[:=]?\s*(\d+(?:\.\d+)?)\s*%\s*(?:kcal|calories?|energy)",
        r"(?:high[- ]fat diet|high fat diet|hfd)\s*\(?\s*"
        r"(\d+(?:\.\d+)?)\s*%\s*(?:kcal|calories?|energy)?\s*\)?",
        r"(\d+(?:\.\d+)?)\s*%\s*(?:kcal|calories?|energy)\s*"
        r"(?:high[- ]fat diet|high fat diet|hfd)",
    ]
    for pat in pats_high:
        m = re.search(pat, t, re.I)
        if m:
            return float(m.group(1)), "HIGH", m.group(0)

    # Medium confidence only if fat/HFD context is explicit.
    pats_medium = [
        r"(\d+(?:\.\d+)?)\s*%\s*(?:high[- ]?)?fat(?:\s+diet)?",
        r"(?:high[- ]fat diet|high fat diet|hfd)[^|;]{0,20}?"
        r"(\d+(?:\.\d+)?)\s*%",
    ]
    for pat in pats_medium:
        m = re.search(pat, t, re.I)
        if m:
            return float(m.group(1)), "MEDIUM", m.group(0)

    return None, "", ""

def extract_hfd_fat(sample, series):
    sample_text = " | ".join(sample["title"] + sample["characteristics_ch1"])
    v, conf, raw = fat_percent_from_text(sample_text)
    if v is not None:
        return make_result(str(v), raw, "GEO_SOFT_SAMPLE_TEXT", conf)
    candidates = []
    for chunk in series["title"] + series["summary"] + series["overall_design"]:
        vv, cc, rr = fat_percent_from_text(chunk)
        if vv is not None:
            candidates.append((vv, cc, rr))
    vals = sorted(set(x[0] for x in candidates))
    if len(vals) == 1:
        return make_result(str(vals[0]),
                           " | ".join(unique_nonempty(x[2] for x in candidates)),
                           "GEO_SOFT_SERIES_TEXT", "MEDIUM")
    if len(vals) > 1:
        return make_result("CONFLICT",
                           " | ".join(unique_nonempty(x[2] for x in candidates)),
                           "GEO_SOFT_SERIES_TEXT", "LOW",
                           "Multiple fat percentages in series text")
    return make_result()

def extract_hfd_description(sample):
    pairs = characteristic_pairs(sample["characteristics_ch1"])
    parts = []
    for k, v, raw in pairs:
        if any(word in k for word in ["diet", "treatment", "feeding", "fat"]):
            parts.append(raw)
    if not parts:
        parts = sample["title"]
    parts = unique_nonempty(parts)
    return make_result(" | ".join(parts), " | ".join(parts),
                       "GEO_SOFT_SAMPLE_TEXT", "RAW")

def collapse_field(sub, field_prefix, subset_role=None):
    z = sub if subset_role is None else sub[sub["Role"] == subset_role]
    value_col = f"{field_prefix}_value"
    raw_col = f"{field_prefix}_raw"
    source_col = f"{field_prefix}_source"
    conf_col = f"{field_prefix}_confidence"
    vals = unique_nonempty(z[value_col].tolist())
    if len(vals) == 0:
        return {"value": "UNKNOWN", "raw": "", "source": "", "confidence": "", "status": "MISSING"}
    if "CONFLICT" in vals:
        return {
            "value": "CONFLICT",
            "raw": " || ".join(unique_nonempty(z[raw_col].tolist())),
            "source": " || ".join(unique_nonempty(z[source_col].tolist())),
            "confidence": " || ".join(unique_nonempty(z[conf_col].tolist())),
            "status": "CONFLICT"
        }
    if len(vals) == 1:
        return {
            "value": vals[0],
            "raw": " || ".join(unique_nonempty(z[raw_col].tolist())),
            "source": " || ".join(unique_nonempty(z[source_col].tolist())),
            "confidence": " || ".join(unique_nonempty(z[conf_col].tolist())),
            "status": "RESOLVED"
        }
    return {
        "value": "MIXED",
        "raw": " || ".join(unique_nonempty(z[raw_col].tolist())),
        "source": " || ".join(unique_nonempty(z[source_col].tolist())),
        "confidence": " || ".join(unique_nonempty(z[conf_col].tolist())),
        "status": "MIXED"
    }

if not MANIFEST.exists():
    raise FileNotFoundError(MANIFEST)

m = pd.read_csv(MANIFEST, sep="\t", dtype=str).fillna("")
required = {"GSE", "GSM", "SRR", "Contrast_id", "Role"}
missing = required - set(m.columns)
if missing:
    raise RuntimeError(f"Manifest missing required columns: {sorted(missing)}")

runmap = m[["GSE", "GSM", "SRR", "Contrast_id", "Role"]].drop_duplicates().copy()
assert runmap["GSE"].nunique() == 64
assert runmap["Contrast_id"].nunique() == 105
assert runmap["GSM"].nunique() == 837
assert runmap["SRR"].nunique() == 951
assert runmap.groupby("GSM")["Role"].nunique().max() == 1
assert runmap.groupby("GSM")["Contrast_id"].nunique().max() == 1

gsmap = runmap[["GSE", "GSM", "Contrast_id", "Role"]].drop_duplicates().sort_values(
    ["GSE", "Contrast_id", "Role", "GSM"]
)

sample_rows = []
soft_status_rows = []

for gse, gz in gsmap.groupby("GSE"):
    candidates = sorted(SOFTDIR.rglob(f"{gse}_family.soft.gz"))
    if not candidates:
        candidates = sorted(SOFTDIR.rglob(f"{gse}_family.soft"))

    if not candidates:
        soft_status_rows.append({"GSE": gse, "SOFT_status": "MISSING", "SOFT_path": ""})
        for _, r in gz.iterrows():
            sample_rows.append({
                "GSE": gse, "GSM": r["GSM"], "Contrast_id": r["Contrast_id"],
                "Role": r["Role"], "SOFT_status": "MISSING",
            })
        continue

    soft = candidates[0]
    samples, series = parse_soft(soft)
    soft_status_rows.append({
        "GSE": gse, "SOFT_status": "FOUND", "SOFT_path": str(soft.relative_to(ROOT))
    })

    for _, r in gz.iterrows():
        gsm = r["GSM"]
        sample = samples.get(gsm)
        base = {
            "GSE": gse, "GSM": gsm, "Contrast_id": r["Contrast_id"],
            "Role": r["Role"], "SOFT_status": "FOUND" if sample else "GSM_NOT_FOUND_IN_SOFT",
        }
        if sample is None:
            sample_rows.append(base)
            continue

        duration_res = extract_duration(sample)
        if r["Role"] == "HFD" and not duration_res["value"]:
            duration_res = extract_series_duration(series)

        fields = {
            "sex": extract_sex(sample),
            "strain": extract_strain(sample),
            "age_weeks": extract_age(sample),
            "duration_weeks": duration_res if r["Role"] == "HFD" else make_result(),
            "fasting_state": extract_fasting(sample),
            "control_diet_type": extract_control_diet(sample) if r["Role"] == "CONTROL" else make_result(),
            "hfd_fat_kcal_percent": extract_hfd_fat(sample, series) if r["Role"] == "HFD" else make_result(),
            "hfd_diet_description": extract_hfd_description(sample) if r["Role"] == "HFD" else make_result(),
        }

        for name, res in fields.items():
            for k in ["value", "raw", "source", "confidence", "note"]:
                base[f"{name}_{k}"] = res[k]

        base["sample_title_raw"] = " | ".join(unique_nonempty(sample["title"]))
        base["sample_characteristics_raw"] = " | ".join(unique_nonempty(sample["characteristics_ch1"]))
        base["series_title_raw"] = " | ".join(unique_nonempty(series["title"]))
        base["series_overall_design_raw"] = " | ".join(unique_nonempty(series["overall_design"]))
        sample_rows.append(base)

s = pd.DataFrame(sample_rows)
assert len(s) == 837
assert s["GSM"].nunique() == 837
assert s["Contrast_id"].nunique() == 105
assert s["GSE"].nunique() == 64

field_names = [
    "sex", "strain", "age_weeks", "duration_weeks", "fasting_state",
    "control_diet_type", "hfd_fat_kcal_percent", "hfd_diet_description",
]
for f in field_names:
    for suffix in ["value", "raw", "source", "confidence", "note"]:
        col = f"{f}_{suffix}"
        if col not in s.columns:
            s[col] = ""
s = s.fillna("")

srr_per_gsm = runmap.groupby("GSM")["SRR"].nunique().rename("SRR_count")
s = s.merge(srr_per_gsm, left_on="GSM", right_index=True, how="left", validate="one_to_one")

contrast_rows = []

for cid, z in s.groupby("Contrast_id", sort=True):
    gses = unique_nonempty(z["GSE"].tolist())
    assert len(gses) == 1
    gse = gses[0]

    row = {
        "Contrast_id": cid,
        "GSE": gse,
        "n_GSM": z["GSM"].nunique(),
        "n_CONTROL_GSM": z.loc[z["Role"] == "CONTROL", "GSM"].nunique(),
        "n_HFD_GSM": z.loc[z["Role"] == "HFD", "GSM"].nunique(),
        "n_SRR": int(runmap.loc[runmap["Contrast_id"] == cid, "SRR"].nunique()),
    }

    for outname in ["sex", "strain", "age_weeks", "fasting_state"]:
        q = collapse_field(z, outname, subset_role=None)
        for k, v in q.items():
            row[f"{outname}_{k}"] = v

    role_specific = {
        "hfd_duration_weeks": ("duration_weeks", "HFD"),
        "hfd_fat_kcal_percent": ("hfd_fat_kcal_percent", "HFD"),
        "control_diet_type": ("control_diet_type", "CONTROL"),
        "hfd_diet_description": ("hfd_diet_description", "HFD"),
    }

    for outname, (srcname, role) in role_specific.items():
        q = collapse_field(z, srcname, subset_role=role)
        for k, v in q.items():
            row[f"{outname}_{k}"] = v

    review_fields = [
        "sex_status", "strain_status", "age_weeks_status",
        "fasting_state_status", "hfd_duration_weeks_status",
        "hfd_fat_kcal_percent_status", "control_diet_type_status",
    ]
    statuses = [row[x] for x in review_fields]
    if any(x == "CONFLICT" for x in statuses):
        overall = "CONFLICT"
    elif any(x == "MIXED" for x in statuses):
        overall = "MIXED"
    elif all(x == "RESOLVED" for x in statuses):
        overall = "COMPLETE_GEO_AUTO"
    else:
        overall = "INCOMPLETE_GEO_AUTO"

    row["moderator_auto_status"] = overall
    row["missing_or_review_fields"] = ",".join(
        name.replace("_status", "")
        for name in review_fields
        if row[name] != "RESOLVED"
    )
    contrast_rows.append(row)

c = pd.DataFrame(contrast_rows).sort_values(["GSE", "Contrast_id"])
assert len(c) == 105
assert c["Contrast_id"].nunique() == 105
assert c["GSE"].nunique() == 64
assert c["n_GSM"].sum() == 837

coverage_fields = [
    "sex", "strain", "age_weeks", "fasting_state",
    "hfd_duration_weeks", "hfd_fat_kcal_percent", "control_diet_type",
]

coverage_rows = []
for f in coverage_fields:
    status = c[f"{f}_status"]
    resolved = int((status == "RESOLVED").sum())
    mixed = int((status == "MIXED").sum())
    conflict = int((status == "CONFLICT").sum())
    missing = int((status == "MISSING").sum())
    coverage_rows.append({
        "Moderator": f,
        "Total_contrasts": len(c),
        "Resolved": resolved,
        "Resolved_percent": resolved / len(c) * 100.0,
        "Mixed": mixed,
        "Conflict": conflict,
        "Missing": missing,
        "Recommended_next_step": (
            "CANDIDATE_FOR_ANALYSIS"
            if resolved / len(c) >= 0.70
            else "PUBLICATION_OR_MANUAL_CURATION_NEEDED"
        )
    })
coverage = pd.DataFrame(coverage_rows)

review = c[c["moderator_auto_status"] != "COMPLETE_GEO_AUTO"].copy()
review_cols = [
    "GSE", "Contrast_id", "n_CONTROL_GSM", "n_HFD_GSM",
    "moderator_auto_status", "missing_or_review_fields",
]
for f in coverage_fields:
    review_cols += [f"{f}_value", f"{f}_status", f"{f}_source", f"{f}_confidence"]
review = review[review_cols]

dictionary = pd.DataFrame([
    {
        "Moderator": "sex",
        "Definition": "Sex explicitly reported for GSMs in the contrast",
        "Primary_source": "GEO sample characteristics/title",
        "Notes": "MIXED retained if >1 sex occurs in the contrast"
    },
    {
        "Moderator": "strain",
        "Definition": "Mouse strain/substrain explicitly reported",
        "Primary_source": "GEO sample characteristics",
        "Notes": "No substrain guessed when not explicit"
    },
    {
        "Moderator": "age_weeks",
        "Definition": "Explicit sample age normalized to weeks for week/day units",
        "Primary_source": "GEO sample characteristics",
        "Notes": "Months are not automatically converted"
    },
    {
        "Moderator": "hfd_duration_weeks",
        "Definition": "Duration of HFD exposure in the HFD arm",
        "Primary_source": "GEO sample characteristics",
        "Notes": "Only explicit duration fields parsed automatically"
    },
    {
        "Moderator": "hfd_fat_kcal_percent",
        "Definition": "Percent dietary energy/kcal from fat for HFD arm",
        "Primary_source": "GEO sample text, then GEO series text",
        "Notes": "Explicit kcal/energy wording HIGH; generic '% fat' MEDIUM"
    },
    {
        "Moderator": "fasting_state",
        "Definition": "Sampling state: fasted versus fed/ad libitum",
        "Primary_source": "GEO sample characteristics/title",
        "Notes": "MIXED/CONFLICT not collapsed"
    },
    {
        "Moderator": "control_diet_type",
        "Definition": "Control diet descriptor: chow/LFD/normal/control-unspecified",
        "Primary_source": "CONTROL GSM GEO text",
        "Notes": "Multiple descriptors retained rather than guessed"
    },
])

sample_out = OUT / "HFD_primary837_GSM_moderator_extraction_GEOauto_v1_1_1.tsv"
contrast_out = OUT / "HFD_primary105_contrast_moderators_GEOauto_v1_1_1.tsv"
coverage_out = OUT / "HFD_primary105_moderator_coverage_GEOauto_v1_1_1.tsv"
review_out = OUT / "HFD_primary105_moderator_review_queue_v1_1_1.tsv"
soft_out = OUT / "HFD_primary64_GEO_SOFT_availability_v1_1_1.tsv"
dict_out = OUT / "HFD_moderator_dictionary_v1_1_1.tsv"

s.sort_values(["GSE", "Contrast_id", "Role", "GSM"]).to_csv(sample_out, sep="\t", index=False)
c.to_csv(contrast_out, sep="\t", index=False)
coverage.to_csv(coverage_out, sep="\t", index=False)
review.to_csv(review_out, sep="\t", index=False)
pd.DataFrame(soft_status_rows).sort_values("GSE").to_csv(soft_out, sep="\t", index=False)
dictionary.to_csv(dict_out, sep="\t", index=False)

print("===== HFD BIOLOGICAL MODERATOR MASTER: GEO AUTO v1.1.1 PASS =====")
print("GSE:", c["GSE"].nunique())
print("Contrasts:", c["Contrast_id"].nunique())
print("GSM:", c["n_GSM"].sum())
print("CONTROL GSM:", c["n_CONTROL_GSM"].sum())
print("HFD GSM:", c["n_HFD_GSM"].sum())
print("SRR:", runmap["SRR"].nunique())

print("\n===== GEO SOFT AVAILABILITY =====")
print(pd.DataFrame(soft_status_rows)["SOFT_status"].value_counts().to_string())

print("\n===== MODERATOR COVERAGE =====")
print(
    coverage[
        ["Moderator", "Resolved", "Resolved_percent", "Mixed",
         "Conflict", "Missing", "Recommended_next_step"]
    ].to_string(index=False)
)

print("\n===== AUTO STATUS =====")
print(c["moderator_auto_status"].value_counts().to_string())

print("\nSaved:")
for p in [sample_out, contrast_out, coverage_out, review_out, soft_out, dict_out]:
    print(p.relative_to(ROOT))
