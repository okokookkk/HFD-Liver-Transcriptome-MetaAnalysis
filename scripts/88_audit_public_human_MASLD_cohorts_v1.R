#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("GEOquery", quietly = TRUE))
    stop("Required R package not installed: GEOquery")
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Required R package not installed: data.table")
})

library(GEOquery)
library(data.table)

ROOT <- normalizePath(Sys.getenv("HFD_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)

OUTDIR <- file.path(ROOT, "00_metadata", "human_validation", "cohort_audit_v1")
SUPPDIR <- file.path(ROOT, "00_metadata", "human_validation", "supplementary_v1")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPPDIR, recursive = TRUE, showWarnings = FALSE)

OUT_SAMPLE <- file.path(OUTDIR, "human_MASLD_public_cohort_sample_metadata_v1.tsv")
OUT_FILES <- file.path(OUTDIR, "human_MASLD_public_cohort_supplementary_files_v1.tsv")
OUT_PLAN <- file.path(OUTDIR, "human_MASLD_validation_design_v1.tsv")
OUT_SUMMARY <- file.path(OUTDIR, "human_MASLD_public_cohort_audit_summary_v1.txt")

GSES <- c("GSE126848", "GSE130970", "GSE135251")

message("===== PUBLIC HUMAN MASLD COHORT AUDIT v1 =====")

sample_rows <- list()
file_rows <- list()

for (gse in GSES) {
  message("Reading GEO metadata: ", gse)

  gset <- getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE)
  if (length(gset) < 1L) stop(gse, ": GEOquery returned no ExpressionSet")

  # Usually one platform here.
  eset <- gset[[1]]
  ph <- as.data.table(pData(eset), keep.rownames = "GSM")
  ph[, GSE := gse]

  # Preserve all GEO phenotype columns; this is an audit artifact.
  setcolorder(ph, c("GSE", "GSM", setdiff(names(ph), c("GSE", "GSM"))))
  sample_rows[[gse]] <- ph

  # Download supplementary files only; no SRA/raw FASTQ.
  gdir <- file.path(SUPPDIR, gse)
  dir.create(gdir, recursive = TRUE, showWarnings = FALSE)

  supp <- getGEOSuppFiles(
    gse,
    makeDirectory = FALSE,
    baseDir = gdir,
    fetch_files = TRUE
  )

  sf <- as.data.table(supp, keep.rownames = "Downloaded_path")
  sf[, GSE := gse]
  sf[, Filename := basename(Downloaded_path)]
  if ("size" %in% names(sf)) {
    sf[, Size := size]
  } else {
    sf[, Size := NA_real_]
  }

  # Capture URL if GEOquery exposes it.
  url_col <- intersect(c("url", "URL"), names(sf))
  if (length(url_col)) {
    sf[, Source_URL := get(url_col[1])]
  } else {
    sf[, Source_URL := NA_character_]
  }

  file_rows[[gse]] <- sf[, .(
    GSE, Filename, Downloaded_path, Size, Source_URL
  )]
}

samples <- rbindlist(sample_rows, fill = TRUE, use.names = TRUE)
files <- rbindlist(file_rows, fill = TRUE, use.names = TRUE)

fwrite(samples, OUT_SAMPLE, sep = "\t", quote = FALSE, na = "NA")
fwrite(files, OUT_FILES, sep = "\t", quote = FALSE, na = "NA")

# ------------------------------------------------------------------
# Prespecified external-validation design.
# This does not infer phenotype from arbitrary fields; it records the
# intended contrast and requires manual confirmation from downloaded
# GEO sample metadata before DE is run.
# ------------------------------------------------------------------
plan <- data.table(
  GSE = GSES,
  Human_role = c(
    "PRIMARY_EXTERNAL_VALIDATION",
    "PRIMARY_EXTERNAL_VALIDATION",
    "PRIMARY_EXTERNAL_VALIDATION"
  ),
  Intended_primary_contrast = c(
    "NAFL_or_NASH_vs_histologically_normal",
    "NAFLD_spectrum_vs_histologically_normal",
    "NAFLD_vs_control"
  ),
  Intended_control_definition = c(
    "Normal-weight + obese subjects with normal liver histology",
    "Histologically normal liver",
    "Control liver samples"
  ),
  Intended_case_definition = c(
    "NAFL + NASH",
    "NAFLD spectrum",
    "NAFLD cases"
  ),
  Expected_n_control_from_GEO = c(26L, 6L, 10L),
  Expected_n_case_from_GEO = c(31L, 72L, 206L),
  Processed_count_strategy = c(
    "single public raw gene-count matrix",
    "single public tximport gene-count matrix",
    "public per-GSM gene-count files / RAW tar"
  ),
  Raw_FASTQ_required = "NO",
  Controlled_access_required = "NO",
  Manual_sample_group_confirmation_required = "YES"
)

fwrite(plan, OUT_PLAN, sep = "\t", quote = FALSE, na = "NA")

# File audit by cohort.
file_summary <- files[, .(
  N_downloaded_files = .N,
  Total_downloaded_bytes = sum(as.numeric(Size), na.rm = TRUE),
  Example_files = paste(head(Filename, 6L), collapse = ";")
), by = GSE]

# Sample counts from GEO Matrix metadata, not yet biological classification.
sample_summary <- samples[, .(
  GEO_matrix_samples = uniqueN(GSM)
), by = GSE]

audit <- merge(plan, sample_summary, by = "GSE", all.x = TRUE)
audit <- merge(audit, file_summary, by = "GSE", all.x = TRUE)

expected_files <- list(
  GSE126848 = "GSE126848_Gene_counts_raw.txt.gz",
  GSE130970 = "GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz"
)

g126_ok <- any(files[GSE == "GSE126848", Filename] ==
               expected_files$GSE126848)
g130_ok <- any(files[GSE == "GSE130970", Filename] ==
               expected_files$GSE130970)
g135_any <- nrow(files[GSE == "GSE135251"]) > 0L

status <- if (all(c(g126_ok, g130_ok, g135_any))) "PASS" else "REVIEW"

summary <- c(
  "===== PUBLIC HUMAN MASLD COHORT AUDIT SUMMARY v1 =====",
  "",
  paste0("STATUS=", status),
  paste0("GEOquery_version=", as.character(packageVersion("GEOquery"))),
  "",
  "COHORTS=GSE126848;GSE130970;GSE135251",
  "CONTROLLED_ACCESS_REQUIRED=NO",
  "RAW_FASTQ_REPROCESSING_PLANNED=NO",
  "PROCESSED_COUNTS_USED_FOR_EXTERNAL_VALIDATION=YES",
  "",
  "PRESPECIFIED_PRIMARY_HUMAN_VALIDATION:",
  "GSE126848=NAFL_or_NASH_vs_histologically_normal",
  "GSE130970=NAFLD_spectrum_vs_histologically_normal",
  "GSE135251=NAFLD_vs_control",
  "",
  "EXPECTED_SAMPLE_COUNTS_FROM_GEO_DESIGN:",
  "GSE126848=31_case;26_histologically_normal_control",
  "GSE130970=72_case;6_histologically_normal_control",
  "GSE135251=206_case;10_control",
  "",
  paste0("GSE126848_expected_count_matrix_found=", ifelse(g126_ok, "YES", "NO")),
  paste0("GSE130970_expected_count_matrix_found=", ifelse(g130_ok, "YES", "NO")),
  paste0("GSE135251_supplementary_files_downloaded=", ifelse(g135_any, "YES", "NO")),
  "",
  "HUMAN_DE_ANALYSIS_PERFORMED=NO",
  "MOUSE_HUMAN_ORTHOLOG_MAPPING_PERFORMED=NO",
  "HUMAN_META_ANALYSIS_PERFORMED=NO",
  "CONCORDANCE_TEST_PERFORMED=NO",
  "",
  paste0("Sample_metadata=", sub(paste0("^", ROOT, "/"), "", OUT_SAMPLE)),
  paste0("Supplementary_files=", sub(paste0("^", ROOT, "/"), "", OUT_FILES)),
  paste0("Validation_design=", sub(paste0("^", ROOT, "/"), "", OUT_PLAN))
)

writeLines(summary, OUT_SUMMARY)
cat(paste(summary, collapse = "\n"), "\n")

if (status != "PASS") quit(status = 2L)

