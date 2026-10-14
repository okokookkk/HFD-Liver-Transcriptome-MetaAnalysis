#!/usr/bin/env python3
"""
Build canonical Supplementary Tables S1-S22 v1_1.

Supersedes failed v1.
Key correction:
- S3 is built from the 64 final per-GSE sample maps:
    03_counts/gsm/GSE*/GSE*_sample_map_v1.tsv
  This family exactly reproduces:
    64 GSE / 837 GSM / 951 SRR / 951 unique GSE-GSM-SRR

Important:
- Frozen upstream inputs are unchanged.
- No global version replacement is used.
- New outputs are written under tables_v1_1 / supplementary_tables_v1_1.
- v1 partial outputs are preserved.
"""

from pathlib import Path
import csv
import gzip
import hashlib
import shutil
import statistics
import tarfile
import math
import re

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()

OUTDIR = ROOT / "09_supplementary" / "tables_v1_1"
METADIR = ROOT / "00_metadata" / "manuscript" / "supplementary_tables_v1_1"

if OUTDIR.exists() and any(OUTDIR.iterdir()):
    raise RuntimeError(
        f"Refusing to overwrite non-empty corrective output directory: {OUTDIR}. "
        "Create a superseding version instead."
    )
if METADIR.exists() and any(METADIR.iterdir()):
    raise RuntimeError(
        f"Refusing to overwrite non-empty corrective metadata directory: {METADIR}. "
        "Create a superseding version instead."
    )

OUTDIR.mkdir(parents=True, exist_ok=True)
METADIR.mkdir(parents=True, exist_ok=True)

MANIFEST = METADIR / "HFD_supplementary_tables_manifest_v1_1.tsv"
SUMMARY = METADIR / "HFD_supplementary_tables_freeze_summary_v1_1.txt"
README = OUTDIR / "README_Supplementary_Tables_v1_1.txt"
BUNDLE = METADIR / "HFD_supplementary_tables_bundle_v1_1.tar.gz"

# ---------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------
def require(rel):
    p = ROOT / rel
    if not p.exists():
        raise FileNotFoundError(f"Missing required frozen source: {rel}")
    return p

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def text_open(path, mode="rt"):
    if path.name.endswith(".gz"):
        return gzip.open(path, mode, newline="")
    return path.open(mode, newline="")

def read_tsv_dicts(path):
    with text_open(path, "rt") as f:
        r = csv.DictReader(f, delimiter="\t")
        rows = list(r)
        return r.fieldnames or [], rows

def write_tsv(path, header, rows, gzip_output=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    if gzip_output:
        with gzip.open(path, "wt", newline="") as f:
            w = csv.DictWriter(
                f, delimiter="\t", fieldnames=header,
                extrasaction="ignore", lineterminator="\n"
            )
            w.writeheader()
            for row in rows:
                w.writerow({k: row.get(k, "") for k in header})
    else:
        with path.open("w", newline="") as f:
            w = csv.DictWriter(
                f, delimiter="\t", fieldnames=header,
                extrasaction="ignore", lineterminator="\n"
            )
            w.writeheader()
            for row in rows:
                w.writerow({k: row.get(k, "") for k in header})

def count_rows(path):
    with text_open(path, "rt") as f:
        r = csv.reader(f, delimiter="\t")
        try:
            next(r)
        except StopIteration:
            return 0
        return sum(1 for _ in r)

def normalize(s):
    return re.sub(r"[^a-z0-9]+", "", str(s).lower())

def find_col(headers, exact_candidates=(), contains_all=(), required=True):
    norm_to_orig = {normalize(h): h for h in headers}
    for c in exact_candidates:
        if normalize(c) in norm_to_orig:
            return norm_to_orig[normalize(c)]
    if contains_all:
        for h in headers:
            low = h.lower()
            if all(tok.lower() in low for tok in contains_all):
                return h
    if required:
        raise KeyError(
            f"Could not find column. exact={exact_candidates}, "
            f"contains_all={contains_all}\nAvailable={headers}"
        )
    return None

def safe_float(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except Exception:
        return None

def truthy(x):
    return str(x).strip().upper() in {"TRUE","YES","1","PASS"}

manifest_rows = []

def register(item, title, role, outpath, source_paths, expected_rows=None, note=""):
    n = count_rows(outpath)
    if expected_rows is not None and n != expected_rows:
        raise RuntimeError(f"{item}: expected {expected_rows} rows, got {n}")
    manifest_rows.append({
        "Table": item,
        "Title": title,
        "Role": role,
        "Canonical_file": str(outpath.relative_to(ROOT)),
        "Rows": n,
        "SHA256": sha256(outpath),
        "Source_files": " | ".join(str(p.relative_to(ROOT)) for p in source_paths),
        "Note": note,
    })

# ---------------------------------------------------------------------
# S1. 64-GSE frozen study universe
# ---------------------------------------------------------------------
s1_src = require(
    "00_metadata/analysis_inputs/frozen_v1/HFD_primary64_GSE_contrast_multiplicity_v1_1.tsv"
)
s1_out = OUTDIR / "TableS1_frozen_64GSE_study_universe_v1_1.tsv"
shutil.copy2(s1_src, s1_out)
register(
    "S1",
    "Frozen 64-GSE mouse liver HFD study universe",
    "One row per independent GSE; distinct from contrast-level Table S2.",
    s1_out, [s1_src], expected_rows=64
)

# ---------------------------------------------------------------------
# S2. 105 frozen contrast design
# ---------------------------------------------------------------------
s2_src = require(
    "00_metadata/analysis_inputs/frozen_v1/HFD_primary105_contrast_analysis_input_audit_v1_1.tsv"
)
hdr, rows = read_tsv_dicts(s2_src)

s2_cols = [
    "Contrast_id","GSE",
    "N_HFD_GSM","N_CONTROL_GSM","N_TOTAL_GSM",
    "N_HFD_SRR","N_CONTROL_SRR","N_TOTAL_SRR",
    "Any_shared_GSM","HFD_GSMs","CONTROL_GSMs",
    "Constructable",
    "hfd_duration_weeks","sex","strain",
    "hfd_fat_kcal_percent","control_diet_type","collection_ZT",
    "publication_pmid","publication_doi","publication_title",
    "overall_adjudication_status"
]
missing = [c for c in s2_cols if c not in hdr]
if missing:
    raise KeyError(f"S2 missing expected columns: {missing}")

s2_out = OUTDIR / "TableS2_frozen_105_contrast_design_v1_1.tsv"
write_tsv(s2_out, s2_cols, rows)
register(
    "S2",
    "Frozen design of 105 mouse HFD-control contrasts",
    "One row per contrast with biological/technical sample counts, final moderators, and publication identifiers.",
    s2_out, [s2_src], expected_rows=105,
    note="Detailed moderator evidence/provenance is provided separately in Table S9."
)

# ---------------------------------------------------------------------
# S3. Correct final GSM-SRR mapping from 64 per-GSE sample maps.
# ---------------------------------------------------------------------
sample_map_files = sorted(
    (ROOT / "03_counts/gsm").glob("GSE*/GSE*_sample_map_v1.tsv")
)

if len(sample_map_files) != 64:
    raise RuntimeError(
        f"S3: expected 64 final per-GSE sample maps, found {len(sample_map_files)}"
    )

s3_rows = []
triples = set()
pairs = set()
srrs = set()
gses = set()

expected_header = ["GSE","GSM","SRR","Contrast_id","Role"]

for src in sample_map_files:
    h, rr = read_tsv_dicts(src)
    if h != expected_header:
        raise RuntimeError(
            f"S3 schema mismatch in {src.relative_to(ROOT)}\n"
            f"Expected={expected_header}\nObserved={h}"
        )

    for row in rr:
        gse = row["GSE"].strip()
        gsm = row["GSM"].strip()
        srr = row["SRR"].strip()

        if not (gse and gsm and srr):
            raise RuntimeError(f"S3 blank identifier in {src.relative_to(ROOT)}")

        key = (gse,gsm,srr)
        if key in triples:
            raise RuntimeError(f"S3 duplicate GSE/GSM/SRR triplet: {key}")

        triples.add(key)
        pairs.add((gse,gsm))
        srrs.add(srr)
        gses.add(gse)

        x = dict(row)
        x["Source_sample_map"] = str(src.relative_to(ROOT))
        s3_rows.append(x)

if len(gses) != 64:
    raise RuntimeError(f"S3: expected 64 GSE, got {len(gses)}")
if len(pairs) != 837:
    raise RuntimeError(f"S3: expected 837 unique GSE/GSM pairs, got {len(pairs)}")
if len(srrs) != 951:
    raise RuntimeError(f"S3: expected 951 unique SRR, got {len(srrs)}")
if len(triples) != 951:
    raise RuntimeError(
        f"S3: expected 951 unique GSE/GSM/SRR triplets, got {len(triples)}"
    )
if len(s3_rows) != 951:
    raise RuntimeError(f"S3: expected 951 rows, got {len(s3_rows)}")

s3_rows.sort(key=lambda r: (r["GSE"],r["GSM"],r["SRR"]))
s3_header = expected_header + ["Source_sample_map"]

s3_out = OUTDIR / "TableS3_frozen_GSM_SRR_mapping_v1_1.tsv"
write_tsv(s3_out, s3_header, s3_rows)
register(
    "S3",
    "Frozen GSM-to-SRR technical-run mapping",
    "One row per final production sequencing run; derived from the 64 final per-GSE sample maps.",
    s3_out, sample_map_files, expected_rows=951,
    note="Exactly reproduces 64 GSE, 837 GSM, and 951 SRR from the final production audit."
)

# ---------------------------------------------------------------------
# S4. Production exceptions: 7 adapter + 2 strandedness resolutions.
# ---------------------------------------------------------------------
adapter_src = require(
    "00_metadata/processing_production/exception_review/HFD_adapter_exception_downstream_status_v1.tsv"
)
strand_src = require(
    "00_metadata/processing_production/exception_review/HFD_strandedness_exception_resolution_status_v1.tsv"
)
ah, ar = read_tsv_dicts(adapter_src)
sh, sr = read_tsv_dicts(strand_src)

union = ["Exception_type","Source_file"]
for c in ah + sh:
    if c not in union:
        union.append(c)

s4_rows = []
for r in ar:
    x = {
        "Exception_type":"Adapter-trim resolution",
        "Source_file":str(adapter_src.relative_to(ROOT))
    }
    x.update(r)
    s4_rows.append(x)

for r in sr:
    x = {
        "Exception_type":"Strandedness adjudication",
        "Source_file":str(strand_src.relative_to(ROOT))
    }
    x.update(r)
    s4_rows.append(x)

if len(ar) != 7 or len(sr) != 2:
    raise RuntimeError(
        f"S4: expected 7 adapter + 2 strandedness rows; got {len(ar)} + {len(sr)}"
    )

s4_out = OUTDIR / "TableS4_production_exceptions_v1_1.tsv"
write_tsv(s4_out, union, s4_rows)
register(
    "S4",
    "Resolved production exceptions",
    "Explicit technical exceptions and their frozen downstream/resolution status.",
    s4_out, [adapter_src, strand_src], expected_rows=9
)

# ---------------------------------------------------------------------
# S5. 105-row contrast-level DE summary.
# ---------------------------------------------------------------------
status_src = require(
    "00_metadata/analysis_outputs/contrast_de_v1/HFD_primary105_DESeq2_status_v1.tsv"
)
meta_src = require(
    "00_metadata/analysis_inputs/meta_analysis_v1/HFD_primary105_meta_effects_primaryQC_v1.tsv.gz"
)

status_h, status_rows = read_tsv_dicts(status_src)
contrast_status_col = find_col(status_h, ("Contrast_id","contrast_id"))

with gzip.open(meta_src, "rt", newline="") as f:
    rr = csv.DictReader(f, delimiter="\t")

    needed = {
        "Contrast_id","GSE","log2FoldChange","lfcSE",
        "Primary_effect_QC","DESeq2_Cooks_flag"
    }
    missing = needed - set(rr.fieldnames or [])
    if missing:
        raise KeyError(f"S5 meta input missing: {sorted(missing)}")

    agg = {}
    for row in rr:
        cid = row["Contrast_id"]
        a = agg.setdefault(cid, {
            "GSE": row["GSE"],
            "N_gene_effect_rows": 0,
            "N_primary_effect_QC_PASS": 0,
            "N_Cooks_flag_TRUE": 0,
            "abs_effects": [],
            "ses": [],
        })

        a["N_gene_effect_rows"] += 1

        if truthy(row["Primary_effect_QC"]):
            a["N_primary_effect_QC_PASS"] += 1
        if truthy(row["DESeq2_Cooks_flag"]):
            a["N_Cooks_flag_TRUE"] += 1

        e = safe_float(row["log2FoldChange"])
        se = safe_float(row["lfcSE"])

        if e is not None:
            a["abs_effects"].append(abs(e))
        if se is not None:
            a["ses"].append(se)

derived = {}
for cid, a in agg.items():
    derived[cid] = {
        "GSE_effect_input": a["GSE"],
        "N_gene_effect_rows": a["N_gene_effect_rows"],
        "N_primary_effect_QC_PASS": a["N_primary_effect_QC_PASS"],
        "N_Cooks_flag_TRUE": a["N_Cooks_flag_TRUE"],
        "median_abs_log2FC":
            statistics.median(a["abs_effects"]) if a["abs_effects"] else "",
        "median_lfcSE":
            statistics.median(a["ses"]) if a["ses"] else "",
    }

s5_header = status_h + [
    "GSE_effect_input","N_gene_effect_rows","N_primary_effect_QC_PASS",
    "N_Cooks_flag_TRUE","median_abs_log2FC","median_lfcSE"
]

s5_rows = []
for r in status_rows:
    cid = r[contrast_status_col]
    if cid not in derived:
        raise RuntimeError(f"S5: no frozen effect rows found for contrast {cid}")
    x = dict(r)
    x.update(derived[cid])
    s5_rows.append(x)

if len(s5_rows) != 105:
    raise RuntimeError(f"S5 expected 105 contrast rows, got {len(s5_rows)}")

s5_out = OUTDIR / "TableS5_contrast_level_DE_summary_v1_1.tsv"
write_tsv(s5_out, s5_header, s5_rows)
register(
    "S5",
    "Contrast-level differential-expression and effect-input summary",
    "One row per frozen contrast with DESeq2 status and compact aggregate effect/QC metrics.",
    s5_out, [status_src, meta_src], expected_rows=105
)

# ---------------------------------------------------------------------
# S6-S16. Straight copies from frozen canonical outputs.
# ---------------------------------------------------------------------
copy_specs = [
    ("S6","Complete mouse gene-level multilevel meta-analysis",
     "Complete mouse gene meta-analysis.",
     "05_meta/gene_level_v1_2/HFD_gene_multilevel_meta_primary_annotated_v1_2.tsv.gz",
     "TableS6_complete_mouse_gene_meta_v1_1.tsv.gz",18234),

    ("S7","Final 59-gene statistical stable core",
     "Final stable core with exact leave-one-GSE-out robustness fields.",
     "05_meta/stable_core_v1/HFD_stable_core_final_v1.tsv",
     "TableS7_final_59_gene_stable_core_v1_1.tsv",59),

    ("S8","Hallmark enrichment of pooled mouse HFD response",
     "All 50 Hallmark mean-response GSEA results.",
     "06_pathway/hallmark_mean_response_v1_2/HFD_Hallmark_GSEA_pooled_meta_zstat_all_v1_2.tsv",
     "TableS8_mouse_Hallmark_mean_response_v1_1.tsv",50),

    ("S9","Moderator curation and provenance",
     "Full 105-contrast manual adjudication/provenance table.",
     "00_metadata/contrasts/moderators/manual_adjudication/HFD_primary105_moderator_manual_adjudication_batch11_v1.tsv",
     "TableS9_moderator_curation_provenance_v1_1.tsv",105),

    ("S10","HFD-duration meta-regression",
     "Gene-wise duration meta-regression results.",
     "05_meta/moderators_v1_2/HFD_duration_meta_regression_v1_2.tsv.gz",
     "TableS10_duration_meta_regression_v1_1.tsv.gz",14421),

    ("S11","Dietary-fat-percentage meta-regression",
     "Gene-wise dietary-fat-percentage meta-regression results.",
     "05_meta/moderators_v1_2/HFD_fat_percent_meta_regression_v1_2.tsv.gz",
     "TableS11_fat_percent_meta_regression_v1_1.tsv.gz",14588),

    ("S12","Control-diet-type meta-regression",
     "Gene-wise control-diet omnibus and category contrasts.",
     "05_meta/moderators_v1_2/HFD_control_diet_meta_regression_v1_2.tsv.gz",
     "TableS12_control_diet_meta_regression_v1_1.tsv.gz",13300),

    ("S13","Duration Hallmark enrichment and fat-adjustment sensitivity",
     "Primary and fat-adjusted duration Hallmark results.",
     "06_pathway/duration_modifier_v1/HFD_duration_Hallmark_primary_vs_adjusted_comparison_v1.tsv",
     "TableS13_duration_Hallmark_results_v1_1.tsv",50),

    ("S14","Hallmark enrichment of transcriptional heterogeneity",
     "Primary standardized excess-Q heterogeneity Hallmark results.",
     "06_pathway/heterogeneity_v1/HFD_Hallmark_heterogeneity_standardized_Q_GSEA_v1.tsv",
     "TableS14_heterogeneity_Hallmark_results_v1_1.tsv",50),

    ("S15","Integrated Hallmark response classes",
     "Integrated mean-response, duration-sensitivity, and heterogeneity profile for all 50 Hallmarks.",
     "00_metadata/analysis_outputs/integrated_pathway_v1/HFD_Hallmark_integrated_profile_classes_v1.tsv",
     "TableS15_integrated_Hallmark_classes_v1_1.tsv",50),

    ("S16","Frozen human MASLD sample groups",
     "Frozen primary human validation sample grouping and phenotype fields.",
     "00_metadata/human_validation/frozen_groups_v1_2/human_MASLD_primary_sample_groups_frozen_v1_2.tsv",
     "TableS16_human_frozen_sample_groups_v1_1.tsv",351),
]

for item, title, role, relsrc, outname, nrows in copy_specs:
    src = require(relsrc)
    out = OUTDIR / outname
    shutil.copy2(src, out)
    register(item,title,role,out,[src],expected_rows=nrows)

# ---------------------------------------------------------------------
# S17. All four human cohort DE comparisons.
# ---------------------------------------------------------------------
s17_sources = [
    ("GSE126848_primary",
     require("07_human_validation/cohort_DE_v1_2/GSE126848_MASLD_vs_control_voom_DE_v1_2.tsv.gz")),

    ("GSE130970_primary_control6",
     require("07_human_validation/cohort_DE_v1_2/GSE130970_MASLD_vs_control6_voom_DE_v1_2.tsv.gz")),

    ("GSE130970_strict_control4",
     require("07_human_validation/cohort_DE_v1_2/GSE130970_MASLD_vs_strict_control4_voom_DE_v1_2.tsv.gz")),

    ("GSE135251_primary",
     require("07_human_validation/cohort_DE_v1_2/GSE135251_MASLD_vs_control_voom_DE_v1_2.tsv.gz")),
]

s17_out = OUTDIR / "TableS17_human_cohort_DE_all_comparisons_v1_1.tsv.gz"

with gzip.open(s17_out, "wt", newline="") as out_f:
    writer = None
    total = 0
    common_header = None

    for label, src in s17_sources:
        with gzip.open(src, "rt", newline="") as in_f:
            r = csv.DictReader(in_f, delimiter="\t")

            if common_header is None:
                common_header = r.fieldnames or []
                writer = csv.DictWriter(
                    out_f,
                    delimiter="\t",
                    fieldnames=["Cohort_comparison"] + common_header,
                    lineterminator="\n"
                )
                writer.writeheader()

            elif (r.fieldnames or []) != common_header:
                raise RuntimeError(f"S17 schema mismatch for {src}")

            for row in r:
                writer.writerow({"Cohort_comparison": label, **row})
                total += 1

if total != 62012:
    raise RuntimeError(f"S17 expected 62,012 total DE rows, got {total}")

register(
    "S17",
    "Human cohort-level differential expression",
    "All four frozen human MASLD case-control comparisons in one long table.",
    s17_out, [p for _,p in s17_sources], expected_rows=62012
)

# ---------------------------------------------------------------------
# S18-S20. Frozen copies.
# ---------------------------------------------------------------------
for item, title, role, relsrc, outname, nrows in [
    ("S18","Human three-cohort gene meta-analysis",
     "Random-effects and fixed-effect human meta-analysis across the three external cohorts.",
     "07_human_validation/gene_meta_v1/human_MASLD_three_cohort_random_effects_meta_v1.tsv.gz",
     "TableS18_human_three_cohort_gene_meta_v1_1.tsv.gz",12477),

    ("S19","Frozen Ensembl 116 high-confidence one-to-one ortholog map",
     "Mouse-human high-confidence one-to-one ortholog universe.",
     "00_metadata/human_validation/analysis_inputs/ortholog_map_v1_4/Ensembl116_mouse_human_highconfidence_one2one_orthologs_v1_4.tsv.gz",
     "TableS19_Ensembl116_mouse_human_ortholog_map_v1_1.tsv.gz",15603),

    ("S20","Complete mouse-human gene-level concordance table",
     "All 10,658 translational genes with mouse reproducibility and human concordance endpoints/sensitivities.",
     "07_human_validation/mouse_human_concordance_v1_6/mouse_HFD_human_MASLD_one2one_ortholog_concordance_v1_6.tsv.gz",
     "TableS20_mouse_human_gene_concordance_v1_1.tsv.gz",10658),
]:
    src = require(relsrc)
    out = OUTDIR / outname
    shutil.copy2(src,out)
    register(item,title,role,out,[src],expected_rows=nrows)

# ---------------------------------------------------------------------
# S21. Compact model-summary table from frozen integrated summaries.
# ---------------------------------------------------------------------
final_sum = require(
    "00_metadata/human_validation/analysis_outputs/final_summary_v1/human_validation_integrated_summary_v1.txt"
)
shape_sum = require(
    "00_metadata/human_validation/analysis_outputs/mouse_human_concordance_followup_v1_1/mouse_human_concordance_shape_summary_v1_1.txt"
)

def parse_kv(path):
    d = {}
    for line in path.read_text(errors="replace").splitlines():
        if "=" in line:
            k,v = line.split("=",1)
            d[k.strip()] = v.strip()
    return d

kv = {}
kv.update(parse_kv(final_sum))
kv.update(parse_kv(shape_sum))

def get_first(keys, default=""):
    for k in keys:
        if k in kv and kv[k] != "":
            return kv[k]
    return default

def split_ci(s):
    if not s:
        return "",""
    s = (
        s.replace("_to_"," ")
         .replace(" to "," ")
         .replace("–"," ")
    )
    nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
    return (nums[0],nums[1]) if len(nums)>=2 else ("","")

primary_lo, primary_hi = split_ci(
    get_first(["Primary_95CI","Primary_CI"])
)
q4_lo, q4_hi = split_ci(
    get_first(["Q4_vs_Q1_95CI","Q4_vs_Q1_unadjusted_95CI","Q4_vs_Q1_CI"])
)

s21_rows = [
    {
        "Model":"Primary continuous",
        "Analysis_role":"Primary prespecified",
        "Endpoint":"Human same-direction pooled RE effect",
        "Predictor":"Mouse GSE-balanced direction consistency, per +10 percentage points",
        "OR":get_first(["Primary_OR_per_10pp"]),
        "CI95_lower":primary_lo,
        "CI95_upper":primary_hi,
        "P_value":get_first(["Primary_P"]),
        "Notes":"Primary translational inference."
    },
    {
        "Model":"Adjusted",
        "Analysis_role":"Supportive sensitivity",
        "Endpoint":"Human same-direction pooled RE effect",
        "Predictor":"Mouse consistency + |mouse log2FC| + log10(mouse nGSE)",
        "OR":get_first(["Adjusted_OR_per_10pp","Adjusted_OR"]),
        "CI95_lower":"",
        "CI95_upper":"",
        "P_value":get_first(["Adjusted_P"]),
        "Notes":"Adjusted for effect magnitude and mouse study coverage."
    },
    {
        "Model":"Mouse coverage ≥20 GSE",
        "Analysis_role":"Coverage sensitivity",
        "Endpoint":"Human same-direction pooled RE effect",
        "Predictor":"Mouse consistency, per +10 percentage points",
        "OR":get_first(["Coverage_ge20_OR_per_10pp","Coverage_ge20_OR"]),
        "CI95_lower":"",
        "CI95_upper":"",
        "P_value":get_first(["Coverage_ge20_P"]),
        "Notes":"Restricted to mouse genes observed in at least 20 independent GSEs."
    },
    {
        "Model":"Human fixed-effect sensitivity",
        "Analysis_role":"Human model sensitivity",
        "Endpoint":"Human same-direction fixed-effect pooled effect",
        "Predictor":"Mouse consistency, per +10 percentage points",
        "OR":get_first(["Human_FE_sensitivity_OR_per_10pp","Human_FE_sensitivity_OR"]),
        "CI95_lower":"",
        "CI95_upper":"",
        "P_value":get_first(["Human_FE_sensitivity_P","Human_FE_P"]),
        "Notes":"Human fixed-effect sensitivity analysis."
    },
    {
        "Model":"GSE130970 strict4 sensitivity",
        "Analysis_role":"Human control-definition sensitivity",
        "Endpoint":"Human same-direction pooled RE effect",
        "Predictor":"Mouse consistency, per +10 percentage points",
        "OR":get_first(["Strict4_sensitivity_OR_per_10pp","Strict4_sensitivity_OR"]),
        "CI95_lower":"",
        "CI95_upper":"",
        "P_value":get_first(["Strict4_sensitivity_P","Strict4_P"]),
        "Notes":"Uses GSE130970 strict NAS=0 four-control sensitivity."
    },
    {
        "Model":"Q4 vs Q1",
        "Analysis_role":"Secondary post hoc shape",
        "Endpoint":"Human same-direction pooled RE effect",
        "Predictor":"Highest vs lowest mouse reproducibility quartile",
        "OR":get_first(["Q4_vs_Q1_unadjusted_OR"]),
        "CI95_lower":q4_lo,
        "CI95_upper":q4_hi,
        "P_value":get_first(["Q4_vs_Q1_unadjusted_P","Q4_vs_Q1_P"]),
        "Notes":"Secondary characterization; does not replace continuous primary model."
    },
    {
        "Model":"Q4 vs lower 75%",
        "Analysis_role":"Secondary post hoc shape",
        "Endpoint":"Human same-direction pooled RE effect",
        "Predictor":"Highest mouse reproducibility quartile vs Q1-Q3",
        "OR":get_first(["Q4_vs_lower75_adjusted_OR"]),
        "CI95_lower":"",
        "CI95_upper":"",
        "P_value":get_first(["Q4_vs_lower75_adjusted_P"]),
        "Notes":"Adjusted secondary characterization."
    },
    {
        "Model":"PI excludes zero: unadjusted",
        "Analysis_role":"Secondary sparse subgroup",
        "Endpoint":"Human same-direction pooled RE effect",
        "Predictor":"Mouse 95% prediction interval excludes zero",
        "OR":get_first(["PI_excludes_zero_unadjusted_OR","PI_exclusion_unadjusted_OR"]),
        "CI95_lower":"",
        "CI95_upper":"",
        "P_value":get_first(["PI_excludes_zero_unadjusted_P","PI_exclusion_unadjusted_P"]),
        "Notes":"Only 31 translational genes met PI-exclusion criterion; underpowered."
    },
    {
        "Model":"PI excludes zero: adjusted",
        "Analysis_role":"Secondary sparse subgroup",
        "Endpoint":"Human same-direction pooled RE effect",
        "Predictor":"Mouse 95% prediction interval excludes zero",
        "OR":get_first(["PI_excludes_zero_adjusted_OR","PI_exclusion_adjusted_OR"]),
        "CI95_lower":"",
        "CI95_upper":"",
        "P_value":get_first(["PI_excludes_zero_adjusted_P","PI_exclusion_adjusted_P"]),
        "Notes":"Only 31 translational genes met PI-exclusion criterion; underpowered."
    },
]

s21_header = [
    "Model","Analysis_role","Endpoint","Predictor",
    "OR","CI95_lower","CI95_upper","P_value","Notes"
]

s21_out = OUTDIR / "TableS21_mouse_human_concordance_models_v1_1.tsv"
write_tsv(s21_out, s21_header, s21_rows)
register(
    "S21",
    "Mouse-human concordance regression and sensitivity models",
    "Compact primary, supportive, and secondary model summary.",
    s21_out, [final_sum,shape_sum], expected_rows=9,
    note="Blank CI/P fields indicate values not present in the frozen integrated summaries; no values were reconstructed."
)

# ---------------------------------------------------------------------
# S22. Frozen Hallmark joint classes.
# ---------------------------------------------------------------------
s22_src = require(
    "07_human_validation/final_summary_v1/mouse_human_Hallmark_joint_significance_classes_v1.tsv"
)
s22_out = OUTDIR / "TableS22_mouse_human_Hallmark_joint_classes_v1_1.tsv"
shutil.copy2(s22_src,s22_out)
register(
    "S22",
    "Mouse-human Hallmark concordance and joint significance classes",
    "All 50 Hallmarks with mouse, human RE/FE/strict4 results and joint class.",
    s22_out, [s22_src], expected_rows=50
)

# ---------------------------------------------------------------------
# Final manifest / README / bundle
# ---------------------------------------------------------------------
with MANIFEST.open("w",newline="") as f:
    fieldnames = [
        "Table","Title","Role","Canonical_file","Rows","SHA256",
        "Source_files","Note"
    ]
    w = csv.DictWriter(
        f,delimiter="\t",fieldnames=fieldnames,lineterminator="\n"
    )
    w.writeheader()
    w.writerows(manifest_rows)

titles = [
    ("S1","Frozen 64-GSE mouse liver HFD study universe"),
    ("S2","Frozen design of 105 mouse HFD-control contrasts"),
    ("S3","Frozen GSM-to-SRR technical-run mapping"),
    ("S4","Resolved production exceptions"),
    ("S5","Contrast-level differential-expression and effect-input summary"),
    ("S6","Complete mouse gene-level multilevel meta-analysis"),
    ("S7","Final 59-gene statistical stable core"),
    ("S8","Hallmark enrichment of pooled mouse HFD response"),
    ("S9","Moderator curation and provenance"),
    ("S10","HFD-duration meta-regression"),
    ("S11","Dietary-fat-percentage meta-regression"),
    ("S12","Control-diet-type meta-regression"),
    ("S13","Duration Hallmark enrichment and fat-adjustment sensitivity"),
    ("S14","Hallmark enrichment of transcriptional heterogeneity"),
    ("S15","Integrated Hallmark response classes"),
    ("S16","Frozen human MASLD sample groups"),
    ("S17","Human cohort-level differential expression"),
    ("S18","Human three-cohort gene meta-analysis"),
    ("S19","Frozen Ensembl 116 high-confidence one-to-one ortholog map"),
    ("S20","Complete mouse-human gene-level concordance table"),
    ("S21","Mouse-human concordance regression and sensitivity models"),
    ("S22","Mouse-human Hallmark concordance and joint significance classes"),
]

README.write_text(
    "HFD Meta-analysis Paper — Canonical Supplementary Tables v1.1\n\n"
    "NEW_INFERENCE=NO\n"
    "UPSTREAM_RESULTS_MODIFIED=NO\n\n" +
    "\n".join(f"Table {k}. {v}" for k,v in titles) +
    "\n\nSee HFD_supplementary_tables_manifest_v1_1.tsv "
    "for frozen source paths and SHA256 hashes.\n"
)

if len(manifest_rows) != 22:
    raise RuntimeError(f"Expected 22 manifest rows, got {len(manifest_rows)}")

with tarfile.open(BUNDLE,"w:gz") as tar:
    for q in sorted(OUTDIR.iterdir()):
        tar.add(q,arcname=f"tables_v1_1/{q.name}")
    tar.add(MANIFEST,arcname=MANIFEST.name)
    tar.add(README,arcname=f"tables_v1_1/{README.name}")

summary = [
    "===== HFD CANONICAL SUPPLEMENTARY TABLES FREEZE v1.1 =====",
    "",
    "STATUS=PASS",
    "ROLE=CANONICAL_SUPPLEMENTARY_TABLE_BUILD",
    "NEW_INFERENCE=NO",
    "UPSTREAM_RESULTS_MODIFIED=NO",
    "",
    "CORRECTION_FROM_v1=S3_final_GSM_SRR_source_family",
    "S3_SOURCE_FAMILY=03_counts/gsm/GSE*/GSE*_sample_map_v1.tsv",
    "S3_ASSERT_GSE=64",
    "S3_ASSERT_GSM=837",
    "S3_ASSERT_SRR=951",
    "",
    "TABLES=22",
    "S1_ROWS=64",
    "S2_ROWS=105",
    "S3_ROWS=951",
    "S4_ROWS=9",
    "S5_ROWS=105",
    "S6_ROWS=18234",
    "S7_ROWS=59",
    "S8_ROWS=50",
    "S9_ROWS=105",
    "S10_ROWS=14421",
    "S11_ROWS=14588",
    "S12_ROWS=13300",
    "S13_ROWS=50",
    "S14_ROWS=50",
    "S15_ROWS=50",
    "S16_ROWS=351",
    "S17_ROWS=62012",
    "S18_ROWS=12477",
    "S19_ROWS=15603",
    "S20_ROWS=10658",
    "S21_ROWS=9",
    "S22_ROWS=50",
    "",
    f"MANIFEST={MANIFEST.relative_to(ROOT)}",
    f"TABLE_DIR={OUTDIR.relative_to(ROOT)}",
    f"BUNDLE={BUNDLE.relative_to(ROOT)}",
    "",
    "NEXT_STEP=Upload_bundle_for_publication_ready_Excel_workbook_build"
]

SUMMARY.write_text("\n".join(summary) + "\n")

print("\n".join(summary))

