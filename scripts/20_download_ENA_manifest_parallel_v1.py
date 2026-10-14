#!/usr/bin/env python3

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import pandas as pd
import subprocess
import hashlib
import gzip
import threading

ap = argparse.ArgumentParser()
ap.add_argument("--manifest", required=True)
ap.add_argument("--gse", required=True)
ap.add_argument("--workers", type=int, default=4)
args = ap.parse_args()

MANIFEST = Path(args.manifest)
GSE = args.gse
N_WORKERS = args.workers

RAW = Path("01_raw/fastq") / GSE
QCOUT = Path("00_metadata/processing_pilot") / GSE

RAW.mkdir(parents=True, exist_ok=True)
QCOUT.mkdir(parents=True, exist_ok=True)

x = pd.read_csv(
    MANIFEST,
    sep="\t",
    dtype=str
).fillna("")

x = x[x["GSE"] == GSE].copy()

if x.empty:
    raise SystemExit(f"No rows for {GSE}")

x["FASTQ_bytes"] = pd.to_numeric(
    x["FASTQ_bytes"],
    errors="raise"
)

assert not x["FASTQ_filename"].duplicated().any()

# Strict PE frozen manifest check
pair = (
    x.groupby("SRR")["FASTQ_component"]
    .apply(lambda z: set(z))
)

assert all(z == {"R1", "R2"} for z in pair)

lock = threading.Lock()

def say(*a):
    with lock:
        print(*a, flush=True)

def urlfix(url):
    if url.startswith(("http://", "https://")):
        return url
    return "https://" + url.lstrip("/")

def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(
            lambda: f.read(16 * 1024 * 1024),
            b""
        ):
            h.update(b)
    return h.hexdigest()

def gzip_test(path):
    with gzip.open(path, "rb") as f:
        while f.read(16 * 1024 * 1024):
            pass

def process(row):

    fn = row["FASTQ_filename"]
    dst = RAW / fn

    expected_size = int(row["FASTQ_bytes"])
    expected_md5 = row["ENA_MD5"].lower()
    url = urlfix(row["ENA_URL"])

    say(f"\n===== {fn} =====")

    # Existing complete file
    if dst.exists() and dst.stat().st_size == expected_size:

        say(f"{fn}: full size exists; checking MD5")

        actual_md5 = md5sum(dst)

        if actual_md5 == expected_md5:

            gzip_test(dst)

            say(f"{fn}: EXISTING VALID -> SKIP")

            status = "EXISTING_VALID"

            return {
                "GSE": GSE,
                "GSM": row["GSM"],
                "SRR": row["SRR"],
                "Contrast_id": row["Contrast_id"],
                "Role": row["Role"],
                "FASTQ_component": row["FASTQ_component"],
                "FASTQ_filename": fn,
                "Expected_bytes": expected_size,
                "Observed_bytes": dst.stat().st_size,
                "Expected_MD5": expected_md5,
                "Observed_MD5": actual_md5,
                "Download_status": status,
                "Size_QC": "PASS",
                "MD5_QC": "PASS",
                "Gzip_QC": "PASS"
            }

        say(f"{fn}: full size but MD5 mismatch -> restart")
        dst.unlink()

    if dst.exists():
        say(
            f"{fn}: resume from "
            f"{dst.stat().st_size / 1024**2:.1f} MiB"
        )
    else:
        say(f"{fn}: new download")

    cmd = [
        "curl",
        "-L",
        "--fail",
        "--retry", "10",
        "--retry-delay", "5",
        "--connect-timeout", "30",
        "--continue-at", "-",
        "--output", str(dst),
        url
    ]

    rc = subprocess.run(cmd).returncode

    if rc != 0:
        raise RuntimeError(
            f"{fn}: curl exit={rc}"
        )

    actual_size = dst.stat().st_size

    if actual_size != expected_size:
        raise RuntimeError(
            f"{fn}: SIZE FAIL "
            f"{actual_size} != {expected_size}"
        )

    actual_md5 = md5sum(dst)

    if actual_md5 != expected_md5:
        raise RuntimeError(
            f"{fn}: MD5 FAIL"
        )

    gzip_test(dst)

    say(f"{fn}: SIZE / MD5 / GZIP PASS")

    return {
        "GSE": GSE,
        "GSM": row["GSM"],
        "SRR": row["SRR"],
        "Contrast_id": row["Contrast_id"],
        "Role": row["Role"],
        "FASTQ_component": row["FASTQ_component"],
        "FASTQ_filename": fn,
        "Expected_bytes": expected_size,
        "Observed_bytes": actual_size,
        "Expected_MD5": expected_md5,
        "Observed_MD5": actual_md5,
        "Download_status": "DOWNLOADED",
        "Size_QC": "PASS",
        "MD5_QC": "PASS",
        "Gzip_QC": "PASS"
    }

results = []
errors = []

with ThreadPoolExecutor(
    max_workers=N_WORKERS
) as pool:

    futures = {
        pool.submit(process, row):
            row["FASTQ_filename"]
        for _, row in x.iterrows()
    }

    for fut in as_completed(futures):

        fn = futures[fut]

        try:
            results.append(fut.result())

        except Exception as e:
            errors.append((fn, str(e)))
            say(f"ERROR: {fn}: {e}")

if errors:
    print("\n===== ERRORS =====")
    for fn, e in errors:
        print(fn, e)
    raise SystemExit(
        f"{len(errors)} FASTQ file(s) failed"
    )

q = pd.DataFrame(results).sort_values(
    ["GSM","SRR","FASTQ_component"]
)

out = (
    QCOUT /
    f"{GSE}_FASTQ_acquisition_QC_v1.tsv"
)

q.to_csv(
    out,
    sep="\t",
    index=False
)

print(f"\n===== {GSE} ACQUISITION SUMMARY =====")
print("Workers:", N_WORKERS)
print("GSM:", q["GSM"].nunique())
print("SRR:", q["SRR"].nunique())
print("FASTQ:", len(q))
print(
    "GiB:",
    q["Observed_bytes"].sum() / 1024**3
)

print("\nALL DOWNLOAD / SIZE / MD5 / GZIP QC: PASS")
print("Saved:", out)
