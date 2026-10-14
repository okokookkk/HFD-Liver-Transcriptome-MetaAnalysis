#!/usr/bin/env python3
from pathlib import Path
import argparse
import gzip
import json
import re
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime
import pandas as pd

ROOT_DEFAULT = __import__("os").environ.get("HFD_PROJECT_ROOT", ".")

TARGETS = [
    "GSE145665",
    "GSE229697",
    "GSE239490",
    "GSE248297",
    "GSE285614",
    "GSE296500",
    "GSE324134",
]

ADAPTER = "AGATCGGAAGAG"

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=ROOT_DEFAULT)
ap.add_argument(
    "--gse",
    action="append",
    choices=TARGETS,
    help="Run only selected GSE; may be supplied multiple times. Default: all 7.",
)
ap.add_argument("--cutadapt-threads", type=int, default=4)
args = ap.parse_args()

ROOT = Path(args.root)
SELECTED = args.gse if args.gse else TARGETS

MANIFEST = (
    ROOT
    / "00_metadata/technical_harmonization/frozen_v1/"
      "HFD_primary_ENA_R1R2_download_manifest_frozen_v1.tsv"
)
STATUS = (
    ROOT
    / "00_metadata/processing_production/"
      "HFD_primary64_production_status_v1.tsv"
)

OUT_META = ROOT / "00_metadata/processing_production/exception_review"
OUT_META.mkdir(parents=True, exist_ok=True)

CUTADAPT = shutil.which("cutadapt")
FASTQC = shutil.which("fastqc")
MULTIQC = shutil.which("multiqc")

if not CUTADAPT:
    raise RuntimeError("cutadapt not found")
if not FASTQC:
    raise RuntimeError("fastqc not found")
if not MULTIQC:
    raise RuntimeError("multiqc not found")

ver = subprocess.run(
    [CUTADAPT, "--version"], text=True, capture_output=True, check=True
).stdout.strip()
if not ver.startswith("4."):
    raise RuntimeError(f"Unexpected cutadapt version: {ver}")

x = pd.read_csv(MANIFEST, sep="\t", dtype=str).fillna("")
s = pd.read_csv(STATUS, sep="\t", dtype=str).fillna("")

assert x["GSE"].nunique() == 64
assert s["GSE"].nunique() == 64

for gse in SELECTED:
    z = s[s["GSE"].eq(gse)]
    if len(z) != 1:
        raise RuntimeError(f"{gse}: missing/non-unique master production status")
    observed = z.iloc[0]["Production_status"]
    if observed != "REVIEW_REQUIRED_ADAPTER":
        raise RuntimeError(
            f"{gse}: expected REVIEW_REQUIRED_ADAPTER, observed {observed!r}"
        )

def gzip_test(path: Path):
    with gzip.open(path, "rb") as fh:
        while fh.read(16 * 1024 * 1024):
            pass

def fastqc_module(text, name):
    m = re.search(
        rf">>{re.escape(name)}\t([^\n]+)\n(.*?)>>END_MODULE",
        text,
        re.S,
    )
    return (m.group(1), m.group(2)) if m else ("MISSING", "")

def parse_fastqc_zip(zpath: Path):
    with zipfile.ZipFile(zpath) as zz:
        candidates = [n for n in zz.namelist() if n.endswith("/fastqc_data.txt")]
        if len(candidates) != 1:
            raise RuntimeError(
                f"{zpath}: expected one fastqc_data.txt; found {len(candidates)}"
            )
        text = zz.read(candidates[0]).decode("utf-8", errors="replace")

    astat, abody = fastqc_module(text, "Adapter Content")
    ostat, obody = fastqc_module(text, "Overrepresented sequences")

    adapter_max = {}
    header = next(
        (line for line in abody.splitlines() if line.startswith("#Position")),
        None,
    )
    if header:
        cols = header.lstrip("#").split("\t")
        adapters = cols[1:]
        maxima = [0.0] * len(adapters)
        for line in abody.splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            for i, v in enumerate(fields[1:]):
                try:
                    maxima[i] = max(maxima[i], float(v))
                except Exception:
                    pass
        adapter_max = dict(zip(adapters, maxima))

    overrep_rows = [
        line
        for line in obody.splitlines()
        if line and not line.startswith("#")
    ]

    return {
        "Adapter_status": astat,
        "Max_adapter_percent": max(adapter_max.values()) if adapter_max else 0.0,
        "Overrepresented_status": ostat,
        "N_overrepresented_sequences": len(overrep_rows),
    }

def build_pair_table(gse):
    z = x[x["GSE"].eq(gse)].copy()
    if z.empty:
        raise RuntimeError(f"{gse}: absent from frozen manifest")

    rows = []
    for srr, q in z.groupby("SRR", sort=True):
        comp = dict(zip(q["FASTQ_component"], q["FASTQ_filename"]))
        if set(comp) != {"R1", "R2"}:
            raise RuntimeError(
                f"{gse}/{srr}: expected frozen R1/R2 pair, observed {set(comp)}"
            )
        rows.append({"GSE": gse, "SRR": srr, "R1": comp["R1"], "R2": comp["R2"]})
    return pd.DataFrame(rows)

all_audit = []
all_qc = []

for gse in SELECTED:
    print("\n" + "=" * 90)
    print(f"{gse}: adapter-only trimming + QC")

    raw_dir = ROOT / "01_raw/fastq" / gse
    trim_dir = ROOT / "01_processed/adapter_trimmed_fastq" / gse
    fastqc_dir = ROOT / "02_qc/fastqc_adapter_trimmed" / gse
    multiqc_dir = ROOT / "02_qc/multiqc_adapter_trimmed" / gse
    log_dir = ROOT / "logs/production_exception_review" / gse
    meta_dir = OUT_META / gse

    for d in [trim_dir, fastqc_dir, multiqc_dir, log_dir, meta_dir]:
        d.mkdir(parents=True, exist_ok=True)

    pairs = build_pair_table(gse)
    expected_raw = set(pairs["R1"]) | set(pairs["R2"])
    observed_raw = {p.name for p in raw_dir.glob("*.fastq.gz")}
    missing = sorted(expected_raw - observed_raw)
    if missing:
        raise RuntimeError(f"{gse}: missing raw FASTQ: {missing}")

    for row in pairs.itertuples(index=False):
        in1 = raw_dir / row.R1
        in2 = raw_dir / row.R2
        out1 = trim_dir / row.R1
        out2 = trim_dir / row.R2
        part1 = trim_dir / (row.R1 + ".part.gz")
        part2 = trim_dir / (row.R2 + ".part.gz")
        js = meta_dir / f"{row.SRR}_cutadapt_v1_1.json"
        log = log_dir / f"{row.SRR}_cutadapt_v1_1.log"

        # v1 bug cleanup: temporary filenames ended in ".part" rather than ".gz",
        # so cutadapt wrote plain FASTQ. Remove only those known temporary files.
        legacy_part1 = trim_dir / (row.R1 + ".part")
        legacy_part2 = trim_dir / (row.R2 + ".part")
        legacy_part1.unlink(missing_ok=True)
        legacy_part2.unlink(missing_ok=True)

        complete = out1.exists() and out2.exists() and js.exists()
        if complete:
            try:
                gzip_test(out1)
                gzip_test(out2)
                print(f"{row.SRR}: existing trimmed pair validated; skip")
            except Exception:
                complete = False

        if not complete:
            part1.unlink(missing_ok=True)
            part2.unlink(missing_ok=True)
            cmd = [
                CUTADAPT,
                "-j", str(args.cutadapt_threads),
                "-a", ADAPTER,
                "-A", ADAPTER,
                "-O", "5",
                "-m", "20",
                "--json", str(js),
                "-o", str(part1),
                "-p", str(part2),
                str(in1),
                str(in2),
            ]
            print(f"{row.SRR}: cutadapt")
            with log.open("w") as fh:
                proc = subprocess.run(
                    cmd, stdout=fh, stderr=subprocess.STDOUT, text=True
                )
            if proc.returncode != 0:
                raise RuntimeError(
                    f"{gse}/{row.SRR}: cutadapt failed rc={proc.returncode}; see {log}"
                )
            gzip_test(part1)
            gzip_test(part2)
            part1.replace(out1)
            part2.replace(out2)

        dat = json.loads(js.read_text())
        rc = dat.get("read_counts", {})
        all_audit.append(
            {
                "GSE": gse,
                "SRR": row.SRR,
                "Raw_R1": str(in1.relative_to(ROOT)),
                "Raw_R2": str(in2.relative_to(ROOT)),
                "Trimmed_R1": str(out1.relative_to(ROOT)),
                "Trimmed_R2": str(out2.relative_to(ROOT)),
                "Input_pairs": rc.get("input"),
                "Output_pairs": rc.get("output"),
                "Pairs_too_short": rc.get("filtered", {}).get("too_short"),
                "Adapter_sequence": ADAPTER,
                "Minimum_overlap": 5,
                "Minimum_length": 20,
                "Quality_trimming": "NO",
                "Cutadapt_version": ver,
            }
        )

    for pth in sorted(trim_dir.glob("*.fastq.gz")):
        stem = pth.name[:-9]
        z = fastqc_dir / f"{stem}_fastqc.zip"
        h = fastqc_dir / f"{stem}_fastqc.html"
        if z.exists() and h.exists():
            continue
        log = log_dir / f"{pth.name}_FastQC_adapter_trimmed_v1_1.log"
        with log.open("w") as fh:
            proc = subprocess.run(
                [FASTQC, "-t", "1", "-o", str(fastqc_dir), str(pth)],
                stdout=fh,
                stderr=subprocess.STDOUT,
                text=True,
            )
        if proc.returncode != 0:
            raise RuntimeError(f"{gse}: FastQC failed for {pth.name}; see {log}")

    trimmed = sorted(trim_dir.glob("*.fastq.gz"))
    zips = sorted(fastqc_dir.glob("*_fastqc.zip"))
    if len(zips) != len(trimmed):
        raise RuntimeError(
            f"{gse}: trimmed FASTQ/FastQC ZIP mismatch: {len(trimmed)} vs {len(zips)}"
        )

    mq_name = f"{gse}_adapter_trimmed_MultiQC_v1_1.html"
    mq_log = log_dir / f"{gse}_MultiQC_adapter_trimmed_v1_1.log"
    with mq_log.open("w") as fh:
        proc = subprocess.run(
            [
                MULTIQC,
                str(fastqc_dir),
                "-o", str(multiqc_dir),
                "-n", mq_name,
                "--force",
            ],
            stdout=fh,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if proc.returncode != 0:
        raise RuntimeError(f"{gse}: MultiQC failed; see {mq_log}")

    raw_fastqc_dir = ROOT / "02_qc/fastqc" / gse
    raw_map = {p.name: p for p in raw_fastqc_dir.glob("*_fastqc.zip")}
    trim_map = {p.name: p for p in fastqc_dir.glob("*_fastqc.zip")}

    for name, tpath in sorted(trim_map.items()):
        if name not in raw_map:
            raise RuntimeError(f"{gse}: missing raw FastQC counterpart for {name}")
        before = parse_fastqc_zip(raw_map[name])
        after = parse_fastqc_zip(tpath)
        all_qc.append(
            {
                "GSE": gse,
                "FastQC_zip": name,
                "Raw_adapter_status": before["Adapter_status"],
                "Raw_max_adapter_percent": before["Max_adapter_percent"],
                "Trimmed_adapter_status": after["Adapter_status"],
                "Trimmed_max_adapter_percent": after["Max_adapter_percent"],
                "Raw_overrepresented_status": before["Overrepresented_status"],
                "Raw_n_overrepresented": before["N_overrepresented_sequences"],
                "Trimmed_overrepresented_status": after["Overrepresented_status"],
                "Trimmed_n_overrepresented": after["N_overrepresented_sequences"],
            }
        )

    marker = meta_dir / f"{gse}_ADAPTER_TRIM_QC_DONE_v1_1.txt"
    marker.write_text(
        "\n".join(
            [
                f"GSE={gse}",
                "STATUS=ADAPTER_TRIM_QC_DONE",
                "SCRIPT_VERSION=v1.1",
                f"TIMESTAMP={datetime.now().astimezone().isoformat()}",
                f"SRR={len(pairs)}",
                f"FASTQ={2 * len(pairs)}",
                f"ADAPTER={ADAPTER}",
                "MINIMUM_OVERLAP=5",
                "MINIMUM_LENGTH=20",
                "QUALITY_TRIMMING=NO",
                "RAW_FASTQ_MODIFIED=NO",
                "FROZEN_TECHNICAL_ELIGIBILITY_CHANGED=NO",
                "DOWNSTREAM_STAR_FEATURECOUNTS_RUN=NO",
            ]
        ) + "\n"
    )

audit = pd.DataFrame(all_audit)
qc = pd.DataFrame(all_qc)

audit_out = OUT_META / "HFD_adapter_exception_cutadapt_audit_v1_1.tsv"
qc_out = OUT_META / "HFD_adapter_exception_FastQC_before_after_v1_1.tsv"
audit.to_csv(audit_out, sep="\t", index=False)
qc.to_csv(qc_out, sep="\t", index=False)

rank = {"pass": 0, "warn": 1, "fail": 2, "MISSING": 3}
def worst_status(vals):
    vals = list(vals)
    return max(vals, key=lambda v: rank.get(v, 99)) if vals else "MISSING"

gse_summary = []
for gse, z in qc.groupby("GSE", sort=True):
    gse_summary.append(
        {
            "GSE": gse,
            "N_FASTQ": len(z),
            "Raw_worst_adapter_status": worst_status(z["Raw_adapter_status"]),
            "Raw_max_adapter_percent": z["Raw_max_adapter_percent"].max(),
            "Trimmed_worst_adapter_status": worst_status(z["Trimmed_adapter_status"]),
            "Trimmed_max_adapter_percent": z["Trimmed_max_adapter_percent"].max(),
            "Trimmed_adapter_PASS_n": int((z["Trimmed_adapter_status"] == "pass").sum()),
            "Trimmed_adapter_WARN_n": int((z["Trimmed_adapter_status"] == "warn").sum()),
            "Trimmed_adapter_FAIL_n": int((z["Trimmed_adapter_status"] == "fail").sum()),
        }
    )

gse_summary = pd.DataFrame(gse_summary)
gse_summary_out = OUT_META / "HFD_adapter_exception_FastQC_GSE_summary_v1_1.tsv"
gse_summary.to_csv(gse_summary_out, sep="\t", index=False)

print("\n===== ADAPTER TRIM + QC COMPLETE =====")
print(gse_summary.to_string(index=False))
print("\nSaved:")
print(" ", audit_out.relative_to(ROOT))
print(" ", qc_out.relative_to(ROOT))
print(" ", gse_summary_out.relative_to(ROOT))
print("\nNo STAR/featureCounts were run. Review trimmed FastQC before downstream.")

