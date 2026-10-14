#!/usr/bin/env python3
import csv, gzip, io, re, subprocess, time, urllib.request
from pathlib import Path

ROOT = Path(".")
LIST = ROOT / "00_metadata" / "frozen_search" / "HFD_primary_analysis_GSE_frozen_v1.txt"
SOFT_DIR = ROOT / "00_metadata" / "geo_soft"
OUT_DIR = ROOT / "00_metadata" / "run_mapping"
LOG_DIR = ROOT / "logs"

for p in [SOFT_DIR, OUT_DIR, LOG_DIR]:
    p.mkdir(parents=True, exist_ok=True)

GSE_FAIL = OUT_DIR / "primary107_geo_failed_v1.txt"
SRX_FAIL = OUT_DIR / "primary107_srx_failed_v1.txt"
SAMPLE_OUT = OUT_DIR / "primary107_geo_samples_v1.tsv"
RUN_OUT = OUT_DIR / "primary107_sra_runs_v1.tsv"
SUMMARY = OUT_DIR / "primary107_mapping_summary_v1.txt"

def geo_bucket(gse):
    return re.sub(r"\d{3}$", "nnn", gse)

def get_soft(gse):
    dest = SOFT_DIR / f"{gse}_family.soft.gz"
    if dest.exists() and dest.stat().st_size > 100:
        return dest
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_bucket(gse)}/{gse}/soft/{gse}_family.soft.gz"
    last = None
    for attempt in range(1, 4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent":"HFD-meta-analysis/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as w:
                while True:
                    b = r.read(1024*1024)
                    if not b:
                        break
                    w.write(b)
            if dest.stat().st_size <= 100:
                raise RuntimeError("downloaded file too small")
            return dest
        except Exception as e:
            last = str(e)
            if dest.exists():
                dest.unlink()
            time.sleep(2 * attempt)
    raise RuntimeError(last or "unknown GEO download failure")

def first(vals, key):
    x = vals.get(key, [])
    return x[0] if x else ""

def parse_soft(gse, path):
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        current = None
        vals = {}
        def flush():
            nonlocal current, vals
            if not current:
                return
            characteristics = vals.get("!Sample_characteristics_ch1", [])
            char_pairs = {}
            for c in characteristics:
                if ":" in c:
                    k, v = c.split(":", 1)
                    k = re.sub(r"\s+", " ", k.strip().lower())
                    char_pairs.setdefault(k, []).append(v.strip())
            def cfind(*keys):
                for key in keys:
                    for ck, cv in char_pairs.items():
                        if ck == key or key in ck:
                            return " | ".join(cv)
                return ""
            relations = vals.get("!Sample_relation", [])
            srx = []
            for rel in relations:
                srx += re.findall(r"\bSRX\d+\b", rel)
            rows.append({
                "GSE": gse,
                "GSM": current,
                "Sample_title": first(vals, "!Sample_title"),
                "Source_name": first(vals, "!Sample_source_name_ch1"),
                "Organism": first(vals, "!Sample_organism_ch1"),
                "Library_strategy": first(vals, "!Sample_library_strategy"),
                "Library_source": first(vals, "!Sample_library_source"),
                "Library_selection": first(vals, "!Sample_library_selection"),
                "Molecule": first(vals, "!Sample_molecule_ch1"),
                "Tissue": cfind("tissue", "source"),
                "Diet": cfind("diet"),
                "Genotype": cfind("genotype"),
                "Sex": cfind("sex", "gender"),
                "Age": cfind("age"),
                "Strain": cfind("strain"),
                "Treatment": cfind("treatment", "drug treatment", "compound"),
                "SRX": ";".join(sorted(set(srx))),
                "All_characteristics": " || ".join(characteristics),
                "All_relations": " || ".join(relations),
            })
            current, vals = None, {}
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                flush()
                current = line.split(" = ", 1)[1].strip()
                vals = {}
            elif current and line.startswith("!Sample_") and " = " in line:
                k, v = line.split(" = ", 1)
                vals.setdefault(k, []).append(v)
        flush()
    return rows

def fetch_runinfo(srx):
    query = f"{srx}[ACCN]"
    cmd = f"esearch -db sra -query '{query}' | efetch -format runinfo"
    last = None
    for attempt in range(1, 4):
        try:
            cp = subprocess.run(
                ["bash", "-lc", cmd],
                text=True, capture_output=True, timeout=180
            )
            if cp.returncode != 0:
                raise RuntimeError(cp.stderr.strip() or f"return code {cp.returncode}")
            txt = cp.stdout.strip()
            if not txt or "\n" not in txt:
                raise RuntimeError("empty RunInfo")
            reader = csv.DictReader(io.StringIO(txt))
            rows = [r for r in reader if r.get("Run")]
            if not rows:
                raise RuntimeError("no SRR rows in RunInfo")
            return rows
        except Exception as e:
            last = str(e)
            time.sleep(2 * attempt)
    raise RuntimeError(last or "unknown SRA failure")

def main():
    gses = [x.strip() for x in LIST.read_text().splitlines() if x.strip()]
    if len(gses) != 107 or len(set(gses)) != 107:
        raise SystemExit(f"Expected 107 unique GSE, found {len(gses)} rows / {len(set(gses))} unique")

    print("========================================")
    print("PRIMARY 107 SAMPLE/RUN MAPPING V1")
    print("========================================")
    print("Target GSE:", len(gses))

    all_samples = []
    geo_failed = []

    for i, gse in enumerate(gses, 1):
        print(f"[GEO {i:03d}/107] {gse}")
        try:
            path = get_soft(gse)
            rows = parse_soft(gse, path)
            all_samples.extend(rows)
            n_rna = sum(r["Library_strategy"].lower() == "rna-seq" for r in rows)
            n_srx = len({s for r in rows for s in r["SRX"].split(";") if s})
            print(f"    samples={len(rows)} RNA-Seq={n_rna} unique_SRX={n_srx}")
        except Exception as e:
            geo_failed.append((gse, str(e)))
            print("    FAILED:", e)

    sample_fields = [
        "GSE","GSM","Sample_title","Source_name","Organism",
        "Library_strategy","Library_source","Library_selection","Molecule",
        "Tissue","Diet","Genotype","Sex","Age","Strain","Treatment",
        "SRX","All_characteristics","All_relations"
    ]
    with SAMPLE_OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=sample_fields, delimiter="\t")
        w.writeheader()
        w.writerows(all_samples)

    with GSE_FAIL.open("w", encoding="utf-8") as f:
        for gse, err in geo_failed:
            f.write(f"{gse}\t{err}\n")

    # Map unique SRX only from RNA-seq samples.
    srx_to_samples = {}
    rnaseq_no_srx = []
    for r in all_samples:
        if r["Library_strategy"].lower() != "rna-seq":
            continue
        srxs = [x for x in r["SRX"].split(";") if x]
        if not srxs:
            rnaseq_no_srx.append((r["GSE"], r["GSM"], r["Sample_title"]))
        for srx in srxs:
            srx_to_samples.setdefault(srx, []).append(r)

    run_rows = []
    srx_failed = []

    unique_srx = sorted(srx_to_samples, key=lambda x:int(x[3:]))
    for i, srx in enumerate(unique_srx, 1):
        print(f"[SRA {i:04d}/{len(unique_srx)}] {srx}")
        try:
            infos = fetch_runinfo(srx)
            for info in infos:
                for s in srx_to_samples[srx]:
                    run_rows.append({
                        "GSE": s["GSE"],
                        "GSM": s["GSM"],
                        "SRX": srx,
                        "SRR": info.get("Run",""),
                        "BioSample": info.get("BioSample",""),
                        "BioProject": info.get("BioProject",""),
                        "LibraryLayout": info.get("LibraryLayout",""),
                        "LibraryStrategy": info.get("LibraryStrategy",""),
                        "LibrarySource": info.get("LibrarySource",""),
                        "Platform": info.get("Platform",""),
                        "Model": info.get("Model",""),
                        "spots": info.get("spots",""),
                        "bases": info.get("bases",""),
                        "avgLength": info.get("avgLength",""),
                        "size_MB": info.get("size_MB",""),
                        "ReleaseDate": info.get("ReleaseDate",""),
                        "Sample_title": s["Sample_title"],
                        "Source_name": s["Source_name"],
                        "Tissue": s["Tissue"],
                        "Diet": s["Diet"],
                        "Genotype": s["Genotype"],
                        "Sex": s["Sex"],
                        "Age": s["Age"],
                        "Strain": s["Strain"],
                        "Treatment": s["Treatment"],
                    })
        except Exception as e:
            srx_failed.append((srx, str(e)))
            print("    FAILED:", e)

    run_fields = [
        "GSE","GSM","SRX","SRR","BioSample","BioProject",
        "LibraryLayout","LibraryStrategy","LibrarySource","Platform","Model",
        "spots","bases","avgLength","size_MB","ReleaseDate",
        "Sample_title","Source_name","Tissue","Diet","Genotype","Sex","Age","Strain","Treatment"
    ]
    with RUN_OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=run_fields, delimiter="\t")
        w.writeheader()
        w.writerows(run_rows)

    with SRX_FAIL.open("w", encoding="utf-8") as f:
        for srx, err in srx_failed:
            f.write(f"{srx}\t{err}\n")

    unique_srr = {r["SRR"] for r in run_rows if r["SRR"]}
    srr_counts = {}
    for r in run_rows:
        if r["SRR"]:
            srr_counts[r["SRR"]] = srr_counts.get(r["SRR"], 0) + 1
    duplicate_srr = {k:v for k,v in srr_counts.items() if v > 1}

    summary = [
        "PRIMARY 107 SAMPLE/RUN MAPPING V1",
        "=================================",
        f"Target GSE: {len(gses)}",
        f"GEO completed: {len(gses)-len(geo_failed)}",
        f"GEO failed: {len(geo_failed)}",
        f"Total GSM parsed: {len(all_samples)}",
        f"RNA-Seq GSM: {sum(r['Library_strategy'].lower() == 'rna-seq' for r in all_samples)}",
        f"RNA-Seq GSM without SRX relation: {len(rnaseq_no_srx)}",
        f"Unique SRX queried: {len(unique_srx)}",
        f"SRX failed: {len(srx_failed)}",
        f"Run rows: {len(run_rows)}",
        f"Unique SRR: {len(unique_srr)}",
        f"SRR appearing in >1 GSE/GSM mapping row: {len(duplicate_srr)}",
        "",
        "Outputs:",
        f"  {SAMPLE_OUT}",
        f"  {RUN_OUT}",
        f"  {GSE_FAIL}",
        f"  {SRX_FAIL}",
    ]
    SUMMARY.write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n" + "\n".join(summary))

if __name__ == "__main__":
    main()
