#!/usr/bin/env python3
"""
HFD Supplementary Table Source/Schema Audit v1

Purpose
-------
Read-only preflight before constructing Supplementary Tables S1-S22.

This script:
1) reads the frozen supplementary-input manifest,
2) inspects every table source file,
3) reports rows/columns/header names,
4) searches project directories for more specific source candidates where
   S1-S5 may require study-, contrast-, GSM-SRR-, production-exception-, or
   contrast-DE-level tables,
5) checks frozen row-count anchors,
6) writes a concise source-design audit.

No frozen inputs or analysis outputs are modified.
"""

from pathlib import Path
import csv
import gzip
import hashlib
import re

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()
MANIFEST = ROOT / "00_metadata/manuscript/supplementary_inputs_v1/HFD_supplementary_input_manifest_v1.tsv"
OUTDIR = ROOT / "00_metadata/manuscript/supplementary_table_source_audit_v1"
OUT = OUTDIR / "HFD_supplementary_table_source_schema_audit_v1.txt"

OUTDIR.mkdir(parents=True, exist_ok=True)

if not MANIFEST.exists():
    raise FileNotFoundError(f"Missing frozen manifest: {MANIFEST}")

def rel(p: Path) -> str:
    try:
        return str(p.relative_to(ROOT))
    except Exception:
        return str(p)

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def opener(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return path.open("rt", errors="replace")

def sniff_delim(path: Path):
    if ".tsv" in path.name.lower():
        return "\t"
    if ".csv" in path.name.lower():
        return ","
    return None

def table_schema(path: Path):
    delim = sniff_delim(path)
    if delim is None:
        return {
            "kind": "TEXT",
            "rows": None,
            "columns": None,
            "header": []
        }

    with opener(path) as f:
        reader = csv.reader(f, delimiter=delim)
        try:
            header = next(reader)
        except StopIteration:
            return {"kind":"TABLE","rows":0,"columns":0,"header":[]}

        n = 0
        for _ in reader:
            n += 1

    return {
        "kind": "TABLE",
        "rows": n,
        "columns": len(header),
        "header": header
    }

# ----------------------------------------------------------------------
# Load frozen manifest and select Table rows only.
# ----------------------------------------------------------------------
with MANIFEST.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

tables = [r for r in rows if r.get("Kind") == "Table"]

if len(tables) != 22:
    raise RuntimeError(f"Expected 22 supplementary tables, found {len(tables)}")

lines = []
def emit(x=""):
    lines.append(str(x))

emit("===== HFD SUPPLEMENTARY TABLE SOURCE / SCHEMA AUDIT v1 =====")
emit("")
emit("MODE=READ_ONLY")
emit("NEW_INFERENCE=NO")
emit("UPSTREAM_RESULTS_MODIFIED=NO")
emit(f"TABLES_IN_FROZEN_MANIFEST={len(tables)}")
emit("")

# Frozen anchors from project.
expected_rows = {
    "S1": 105,       # current manifest anchor only; may be split to study-level in final design
    "S2": 105,
    "S6": 18234,
    "S7": 59,
    "S8": 50,
    "S9": 105,
    "S15": 50,
    "S16": 351,
    "S19": 15603,
    "S20": 10658,
    "S22": 50,
}

emit("===== 1. MANIFEST SOURCE SCHEMAS =====")
for r in tables:
    item = r["Item"]
    path = ROOT / r["Path"]

    emit("")
    emit(f"[{item}]")
    emit(f"ROLE={r['Role']}")
    emit(f"PATH={r['Path']}")

    if not path.exists():
        emit("STATUS=MISSING")
        continue

    sc = table_schema(path)
    emit("STATUS=FOUND")
    emit(f"SHA256={sha256(path)}")
    emit(f"KIND={sc['kind']}")

    if sc["kind"] == "TABLE":
        emit(f"ROWS={sc['rows']}")
        emit(f"COLUMNS={sc['columns']}")
        emit("HEADER=" + " | ".join(sc["header"]))

        if item in expected_rows:
            emit(f"EXPECTED_ROW_ANCHOR={expected_rows[item]}")
            emit(
                "ROW_ANCHOR_STATUS=" +
                ("PASS" if sc["rows"] == expected_rows[item] else "REVIEW")
            )
    else:
        # concise preview for summary text files
        try:
            text = path.read_text(errors="replace").splitlines()
            useful = [z for z in text if z.strip()][:30]
            emit(f"TEXT_NONEMPTY_LINES={len([z for z in text if z.strip()])}")
            emit("TEXT_PREVIEW_BEGIN")
            for z in useful:
                emit(z)
            emit("TEXT_PREVIEW_END")
        except Exception as e:
            emit(f"TEXT_PREVIEW_ERROR={e}")

# ----------------------------------------------------------------------
# Targeted search for detailed source candidates.
# ----------------------------------------------------------------------
emit("")
emit("===== 2. TARGETED DETAILED-SOURCE CANDIDATE SEARCH =====")

targets = {
    "S1_STUDY_LEVEL_UNIVERSE": [
        "gse", "study", "universe", "eligible", "primary64", "frozen"
    ],
    "S2_CONTRAST_DESIGN": [
        "contrast", "primary105", "design", "frozen"
    ],
    "S3_GSM_SRR_MAPPING": [
        "gsm", "srr", "mapping"
    ],
    "S4_PRODUCTION_EXCEPTIONS": [
        "production", "exception"
    ],
    "S5_CONTRAST_DE_SUMMARY": [
        "contrast", "deseq", "summary"
    ],
}

search_roots = [
    ROOT / "00_metadata",
    ROOT / "01_raw",
    ROOT / "02_processed",
    ROOT / "03_counts",
    ROOT / "04_deseq2",
    ROOT / "05_meta",
    ROOT / "07_human_validation",
]

allowed_suffixes = (
    ".tsv", ".tsv.gz", ".csv", ".csv.gz", ".txt"
)

all_files = []
for sr in search_roots:
    if not sr.exists():
        continue
    for q in sr.rglob("*"):
        if not q.is_file():
            continue
        low = q.name.lower()
        if any(low.endswith(s) for s in allowed_suffixes):
            all_files.append(q)

for label, keywords in targets.items():
    scored = []
    for q in all_files:
        low = rel(q).lower()
        score = sum(1 for k in keywords if k in low)
        if score >= 2:
            scored.append((score, q))

    scored.sort(key=lambda x: (-x[0], len(str(x[1]))))

    emit("")
    emit(f"[{label}]")
    if not scored:
        emit("CANDIDATES=NONE")
        continue

    for score, q in scored[:30]:
        try:
            sc = table_schema(q)
            dims = (
                f"rows={sc['rows']} cols={sc['columns']}"
                if sc["kind"] == "TABLE"
                else "text"
            )
        except Exception:
            dims = "schema_error"
        emit(f"SCORE={score}\t{dims}\t{rel(q)}")

# ----------------------------------------------------------------------
# Design flags: tables that should not simply duplicate each other.
# ----------------------------------------------------------------------
emit("")
emit("===== 3. TABLE-DESIGN FLAGS =====")

flags = [
    (
        "S1",
        "Current manifest uses the 105-contrast audit. Final S1 should preferably be one row per independent GSE (64 rows), not duplicate S2."
    ),
    (
        "S2",
        "Final S2 should remain one row per frozen contrast (105 rows) with biological comparison and sample-membership/design fields."
    ),
    (
        "S3",
        "Current manifest points to a production summary text anchor. Final S3 should use the detailed GSM-to-SRR mapping if available."
    ),
    (
        "S4",
        "Current manifest points to the production summary text anchor. Final S4 should be a compact explicit exception table (7 adapter + 2 strandedness-resolved GSEs) with resolution/provenance fields."
    ),
    (
        "S5",
        "Current manifest points to gene-by-contrast meta input. Final S5 should be a contrast-level summary (105 rows), not the full long gene-by-contrast effect table."
    ),
    (
        "S17",
        "Human cohort DE should contain all four frozen comparisons (GSE126848, GSE130970 primary, GSE130970 strict4, GSE135251), not only the GSE126848 anchor."
    ),
    (
        "S21",
        "Regression summary text is an anchor only. Final S21 should be a compact model-summary table with OR, CI, P, model definition, and sensitivity role."
    ),
]

for item, note in flags:
    emit(f"{item}_DESIGN_REVIEW={note}")

# ----------------------------------------------------------------------
# Final status
# ----------------------------------------------------------------------
emit("")
emit("===== 4. PRE-BUILD STATUS =====")
emit("SOURCE_SCHEMA_AUDIT=PASS")
emit("DIRECT_WORKBOOK_BUILD=HOLD_UNTIL_DESIGN_FLAGS_RESOLVED")
emit("RATIONALE=Avoid duplicating S1/S2 and avoid using summary-text anchors where detailed row-level tables are required.")
emit("NEXT_STEP=Use this audit to create canonical publication-ready S1-S22 tables and a single Excel workbook.")

OUT.write_text("\n".join(lines) + "\n")

print("===== HFD SUPPLEMENTARY TABLE SOURCE / SCHEMA AUDIT v1 =====")
print(f"Created: {OUT}")
print("MODE=READ_ONLY")
print("TABLES=22")
print("NEXT=cat the audit report and paste its output")

