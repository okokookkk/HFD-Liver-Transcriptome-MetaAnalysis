#!/usr/bin/env python3
from pathlib import Path
import re
import hashlib

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()

OUTDIR = ROOT / "00_metadata/manuscript/methods_provenance_v1"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT = OUTDIR / "HFD_methods_provenance_collector_v1.txt"

# High-value roots only; read-only scan.
SCAN_ROOTS = [
    ROOT / "00_metadata",
    ROOT / "scripts",
]

TEXT_EXT = {
    ".txt",".tsv",".csv",".md",".log",".out",".err",".sh",".R",".r",".py",
    ".yaml",".yml",".json"
}
MAX_BYTES = 5_000_000

# Patterns are intentionally broad. The report preserves source lines and paths;
# nothing is inferred if an exact item cannot be found.
PATTERNS = {
    "SYSTEMATIC_SEARCH": [
        r"search strateg", r"search date", r"pubmed", r"geo", r"query",
        r"systematic review", r"eligibility", r"screening", r"PRISMA",
        r"records identified", r"records screened"
    ],
    "REFERENCE_GENOME_ANNOTATION": [
        r"GRCm39", r"mm39", r"GRCm38", r"mm10", r"Mus_musculus",
        r"Ensembl.*release", r"GENCODE", r"\.gtf", r"\.fa(sta)?\b",
        r"reference genome", r"annotation"
    ],
    "ALIGNER_QUANTIFIER": [
        r"\bSTAR\b", r"featureCounts", r"Subread", r"salmon", r"kallisto",
        r"hisat2", r"RSEM", r"HTSeq", r"align", r"quantif"
    ],
    "READ_PROCESSING": [
        r"fastp", r"cutadapt", r"Trim Galore", r"trimmomatic",
        r"adapter", r"quality trim", r"read trimming"
    ],
    "STRANDEDNESS": [
        r"stranded", r"strandedness", r"featureCounts.*-s\s*[012]",
        r"TruSeq Stranded", r"GSE137022", r"GSE278920"
    ],
    "SOFTWARE_VERSION": [
        r"version", r"--version", r"packageVersion", r"sessionInfo",
        r"DESeq2", r"metafor", r"edgeR", r"limma", r"fgsea", r"msigdbr",
        r"featureCounts", r"\bSTAR\b"
    ],
    "DESEQ2_MODEL": [
        r"DESeq2", r"DESeqDataSet", r"DESeq\(", r"results\(",
        r"counts?.*10", r"at least 3", r">=\s*10", r"min.*count"
    ],
    "META_MODEL": [
        r"rma\.mv", r"REML", r"random.?effects", r"multilevel",
        r"prediction interval", r"I2", r"tau2", r"direction consistency",
        r"leave-one-GSE", r"LOO"
    ],
    "MODERATORS": [
        r"hfd_duration", r"duration", r"fat_kcal", r"control_diet",
        r"meta-regression", r"moderator", r"pseudo.?R2"
    ],
    "HUMAN_VALIDATION": [
        r"GSE126848", r"GSE130970", r"GSE135251", r"voom", r"filterByExpr",
        r"calcNormFactors", r"lmFit", r"eBayes", r"org\.Hs\.eg\.db",
        r"Compara", r"ortholog_one2one", r"release.?116"
    ],
}

# Files that are especially important and should be shown in full if small.
PRIORITY_FILES = [
    ROOT / "00_metadata/processing_production/final_audit_v1/HFD_primary64_final_production_audit_summary_v1.txt",
    ROOT / "00_metadata/human_validation/analysis_outputs/final_summary_v1/human_validation_integrated_summary_v1.txt",
]

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()

def safe_read(path):
    try:
        if path.stat().st_size > MAX_BYTES:
            return None
        return path.read_text(errors="replace")
    except Exception:
        return None

# Collect candidate files.
files = []
for base in SCAN_ROOTS:
    if not base.exists():
        continue
    for q in base.rglob("*"):
        if not q.is_file():
            continue
        if q.suffix.lower() not in TEXT_EXT:
            continue
        try:
            if q.stat().st_size <= MAX_BYTES:
                files.append(q)
        except Exception:
            pass

# Newest / most relevant scripts first in source display.
files = sorted(set(files), key=lambda x: (str(x.parent), x.name))

report = []
report.append("===== HFD METHODS PROVENANCE COLLECTOR v1 =====")
report.append("")
report.append(f"ROOT={ROOT}")
report.append(f"FILES_SCANNED={len(files)}")
report.append("MODE=READ_ONLY")
report.append("INFERENCE_POLICY=source_lines_only_no_silent_inference")
report.append("")

# Priority frozen audits.
report.append("===== PRIORITY FROZEN AUDITS =====")
for q in PRIORITY_FILES:
    report.append("")
    report.append(f"--- FILE: {q.relative_to(ROOT) if q.exists() else q} ---")
    if not q.exists():
        report.append("STATUS=MISSING")
        continue
    report.append(f"SHA256={sha256(q)}")
    txt = safe_read(q)
    if txt is None:
        report.append("STATUS=TOO_LARGE_OR_UNREADABLE")
    else:
        report.extend(txt.splitlines())

# Pattern matches.
for section, pats in PATTERNS.items():
    regs = [re.compile(p, re.I) for p in pats]
    hits = []
    for q in files:
        txt = safe_read(q)
        if txt is None:
            continue
        for i, line in enumerate(txt.splitlines(), 1):
            if any(r.search(line) for r in regs):
                # Keep useful context but avoid giant generated tables.
                s = line.strip()
                if len(s) > 800:
                    s = s[:800] + " ...[truncated]"
                hits.append((str(q.relative_to(ROOT)), i, s))
    report.append("")
    report.append(f"===== {section} =====")
    report.append(f"MATCHES={len(hits)}")
    # Prefer scripts and audit/summary files; cap report size.
    def score(hit):
        path = hit[0].lower()
        sc = 0
        if path.startswith("scripts/"): sc -= 5
        if "audit" in path or "summary" in path: sc -= 4
        if "final" in path or "frozen" in path: sc -= 3
        return (sc, path, hit[1])
    for path, ln, s in sorted(hits, key=score)[:250]:
        report.append(f"{path}:L{ln}: {s}")

# Targeted script inventory around production / main analyses.
report.append("")
report.append("===== TARGETED SCRIPT INVENTORY =====")
script_dir = ROOT / "scripts"
if script_dir.exists():
    wanted = []
    for q in sorted(script_dir.iterdir()):
        if not q.is_file():
            continue
        name = q.name
        if re.search(
            r"(production|STAR|featureCount|DESeq2|meta|moderator|Hallmark|human|MASLD|ortholog|audit)",
            name, re.I
        ):
            wanted.append(name)
    for name in wanted:
        report.append(name)

# Explicit unresolved checklist.
report.append("")
report.append("===== MANUSCRIPT METHODS FIELDS TO RESOLVE =====")
fields = [
    "SYSTEMATIC_SEARCH_DATABASES",
    "SYSTEMATIC_SEARCH_DATES",
    "SYSTEMATIC_SEARCH_EXACT_QUERY",
    "REFERENCE_GENOME_BUILD",
    "REFERENCE_FASTA_SOURCE_AND_VERSION",
    "GENE_ANNOTATION_SOURCE_AND_RELEASE",
    "ALIGNER_NAME_AND_VERSION",
    "ALIGNER_CRITICAL_OPTIONS",
    "GENE_COUNTING_SOFTWARE_AND_VERSION",
    "GENE_COUNTING_CRITICAL_OPTIONS",
    "READ_TRIMMER_NAME_AND_VERSION",
    "STANDARD_ADAPTER_TRIMMING_OPTIONS",
    "DESEQ2_VERSION",
    "MOUSE_HALLMARK_SOURCE_EXACT_PACKAGE_FILE",
]
for f in fields:
    report.append(f"{f}=RESOLVE_FROM_MATCHES_ABOVE")

OUT.write_text("\n".join(report) + "\n")
print(f"Created: {OUT}")
print("Read-only provenance scan complete.")
print("Please print the final report with:")
print(f"cat {OUT}")

