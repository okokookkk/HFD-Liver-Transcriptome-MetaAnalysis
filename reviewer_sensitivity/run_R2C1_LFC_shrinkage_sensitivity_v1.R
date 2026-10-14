#!/usr/bin/env Rscript
# Reviewer 2 Comment 1: DESeq2 apeglm LFC-shrinkage sensitivity.
# Reproduces the frozen DESeq2 analysis, does not alter any primary files.
# Shrunk estimates are used only for descriptive sign-stability checks,
# NEVER as sampling-error inputs for a new inferential meta-analysis.

suppressPackageStartupMessages({
  for (p in c("DESeq2", "apeglm", "data.table")) {
    if (!requireNamespace(p, quietly = TRUE)) stop("Missing required R package: ", p)
  }
  library(DESeq2)
  library(data.table)
})

ROOT <- normalizePath(Sys.getenv("HFD_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
LIMIT <- if (length(args)) as.integer(args[1L]) else 2L
if (length(LIMIT) != 1L || is.na(LIMIT) || LIMIT < 1L || LIMIT > 105L) {
  stop("Usage: Rscript run_R2C1_LFC_shrinkage_sensitivity_v1.R [2|105]")
}

AUDIT <- file.path(ROOT, "00_metadata/analysis_inputs/frozen_v1",
                   "HFD_primary105_contrast_analysis_input_audit_v1_1.tsv")
CORE <- file.path(ROOT, "05_meta/stable_core_v1",
                  "HFD_stable_core_final_v1.tsv")
OUTDIR <- file.path(ROOT, "00_metadata/analysis_outputs/reviewer2_comment1_shrinkage_v1")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
for (f in c(AUDIT, CORE)) if (!file.exists(f)) stop("Missing required file: ", f)

audit <- fread(AUDIT)
core <- fread(CORE)
if (nrow(audit) != 105L || uniqueN(audit$Contrast_id) != 105L ||
    nrow(core) != 59L || uniqueN(core$Geneid) != 59L) {
  stop("Frozen input sizes differ from 105 contrasts / 59 core genes")
}
if (!all(c("Geneid", "pooled_log2FC") %in% names(core))) {
  stop("Core file missing Geneid or pooled_log2FC")
}
if (!all(audit$Constructable == "PASS")) stop("Some frozen contrasts are not constructable")

split_ids <- function(s) {
  v <- trimws(strsplit(as.character(s), ",", fixed = TRUE)[[1L]])
  v[nzchar(v)]
}

cat("===== R2C1 LFC SHRINKAGE SENSITIVITY =====\n")
cat("Contrasts requested:", LIMIT, "/105\n")
cat("R:", R.version.string, "\n")
cat("DESeq2:", as.character(packageVersion("DESeq2")), "\n")
cat("apeglm:", as.character(packageVersion("apeglm")), "\n")
cat("Output:", OUTDIR, "\n\n")

for (i in seq_len(LIMIT)) {
  r <- audit[i]
  cid <- r$Contrast_id
  gse <- r$GSE
  outfile <- file.path(OUTDIR, paste0(cid, "_core59_apeglm_v1.tsv"))
  cat(sprintf("[%03d/%03d] %s (%s)\n", i, LIMIT, cid, gse))

  if (file.exists(outfile)) {
    old <- tryCatch(fread(outfile), error = function(e) NULL)
    if (!is.null(old) && nrow(old) == 59L && "LFC_shrunk_apeglm" %in% names(old)) {
      cat("  SKIP: existing 59-gene output\n")
      next
    }
    stop("Existing output is incomplete: ", outfile)
  }

  ctl <- split_ids(r$CONTROL_GSMs)
  hfd <- split_ids(r$HFD_GSMs)
  samples <- c(ctl, hfd)
  if (length(ctl) < 3L || length(hfd) < 3L || anyDuplicated(samples)) {
    stop("Frozen sample assignment invalid: ", cid)
  }
  matrix_file <- file.path(ROOT, r$Final_matrix)
  if (!file.exists(matrix_file)) stop("Missing GSM count matrix: ", matrix_file)
  x <- fread(matrix_file, check.names = FALSE)
  if (!"Geneid" %in% names(x) || anyDuplicated(x$Geneid)) {
    stop("Geneid missing or duplicated: ", matrix_file)
  }
  missing_samples <- setdiff(samples, names(x))
  if (length(missing_samples)) stop("Missing sample columns: ", paste(missing_samples, collapse=","))
  mat <- as.matrix(x[, ..samples])
  storage.mode(mat) <- "numeric"
  if (any(!is.finite(mat)) || any(mat < 0) || any(abs(mat - round(mat)) > 1e-8)) {
    stop("Invalid raw counts: ", cid)
  }
  rownames(mat) <- x$Geneid
  keep <- rowSums(mat >= 10) >= 3L  # identical filter to original script 73
  mat <- mat[keep, , drop = FALSE]
  storage.mode(mat) <- "integer"

  group <- factor(c(rep("CONTROL", length(ctl)), rep("HFD", length(hfd))),
                  levels = c("CONTROL", "HFD"))
  coldata <- data.frame(row.names = samples, group = group)
  dds <- DESeqDataSetFromMatrix(countData = mat, colData = coldata, design = ~ group)
  dds <- DESeq(dds, test = "Wald", betaPrior = FALSE,
               minReplicatesForReplace = Inf, quiet = TRUE)
  res <- results(dds, contrast = c("group", "HFD", "CONTROL"),
                 independentFiltering = FALSE, cooksCutoff = FALSE, alpha = 0.05)
  coef_name <- "group_HFD_vs_CONTROL"
  if (!coef_name %in% resultsNames(dds)) {
    stop("Unexpected coefficient name: ", paste(resultsNames(dds), collapse=", "))
  }
  shr <- lfcShrink(dds, coef = coef_name, res = res, type = "apeglm")

  # Compare to original, frozen gene-by-contrast DESeq2 output and QC flags.
  frozen <- file.path(ROOT, "04_de/contrast_level_v1", gse,
                      paste0(cid, "_DESeq2_effects_v1.tsv.gz"))
  if (!file.exists(frozen)) stop("Missing original DESeq2 effects: ", frozen)
  original <- fread(frozen, select = c("Geneid", "log2FoldChange", "lfcSE",
                                       "sampling_variance", "Primary_effect_QC"))
  k_new <- match(core$Geneid, rownames(dds))
  k_old <- match(core$Geneid, original$Geneid)
  k_has <- which(!is.na(k_new))
  mle_new <- se_new <- shrink <- shrink_se <- rep(NA_real_, nrow(core))
  if (length(k_has)) {
    mle_new[k_has] <- res$log2FoldChange[k_new[k_has]]
    se_new[k_has] <- res$lfcSE[k_new[k_has]]
    shrink[k_has] <- shr$log2FoldChange[k_new[k_has]]
    shrink_se[k_has] <- shr$lfcSE[k_new[k_has]]
  }
  mle_orig <- original$log2FoldChange[k_old]
  se_orig <- original$lfcSE[k_old]
  var_orig <- original$sampling_variance[k_old]
  qc_orig <- original$Primary_effect_QC[k_old]
  delta <- abs(mle_new - mle_orig)
  delta <- delta[is.finite(delta)]
  maxdiff <- if (length(delta)) max(delta) else NA_real_
  cat("  Core modeled:", length(k_has), "; original-MLE max difference:",
      format(maxdiff, digits=5), "\n")
  if (is.finite(maxdiff) && maxdiff > 0.01) {
    stop("New MLE differs materially from frozen input (max >0.01); check DESeq2 version")
  }
  out <- data.table(
    Geneid = core$Geneid,
    GSE = gse,
    Contrast_id = cid,
    Original_QC = qc_orig,
    Original_MLE_log2FC = mle_orig,
    Original_Wald_SE = se_orig,
    Original_sampling_variance = var_orig,
    Refit_MLE_log2FC = mle_new,
    Refit_Wald_SE = se_new,
    LFC_shrunk_apeglm = shrink,
    LFC_shrunk_posterior_SD = shrink_se,
    Original_pooled_log2FC = core$pooled_log2FC,
    Original_MLE_reproduction_abs_diff = abs(mle_new - mle_orig)
  )
  tmp <- paste0(outfile, ".tmp")
  fwrite(out, tmp, sep="\t", quote=FALSE, na="NA")
  if (!file.rename(tmp, outfile)) stop("Could not write output: ", outfile)
  cat("  WROTE:", outfile, "\n")
}

# Summarize only when all 105 contrasts have been processed.
if (LIMIT == 105L) {
  paths <- file.path(OUTDIR, paste0(audit$Contrast_id, "_core59_apeglm_v1.tsv"))
  if (!all(file.exists(paths))) stop("105-contrast run incomplete")
  dt <- rbindlist(lapply(paths, fread), use.names=TRUE)
  pass <- dt[Original_QC == "PASS" &
             is.finite(Original_MLE_log2FC) &
             is.finite(LFC_shrunk_apeglm) &
             is.finite(Original_sampling_variance) & Original_sampling_variance > 0]
  pass[, wt := 1 / Original_sampling_variance]
  by_gse <- pass[, .(
    N_contrasts = .N,
    Original_GSE_effect = sum(wt * Original_MLE_log2FC) / sum(wt),
    Shrunk_GSE_effect = sum(wt * LFC_shrunk_apeglm) / sum(wt),
    Primary_pooled_effect = first(Original_pooled_log2FC)
  ), by = .(Geneid, GSE)]
  by_gse[, Original_direction_agrees := sign(Original_GSE_effect) == sign(Primary_pooled_effect)]
  by_gse[, Shrunk_direction_agrees := sign(Shrunk_GSE_effect) == sign(Primary_pooled_effect)]
  by_gse[, GSE_direction_flipped := sign(Original_GSE_effect) != sign(Shrunk_GSE_effect)]
  summary <- by_gse[, .(
    N_GSE_paired = .N,
    Original_GSE_direction_consistency = mean(Original_direction_agrees),
    Shrunk_GSE_direction_consistency = mean(Shrunk_direction_agrees),
    N_GSE_direction_flips = sum(GSE_direction_flipped)
  ), by = Geneid]
  summary <- merge(core[, .(Geneid, pooled_log2FC, GSE_balanced_direction_consistency)],
                   summary, by="Geneid", all.x=TRUE, sort=FALSE)
  summary[, Shrunk_consistency_ge_0_80 := is.finite(Shrunk_GSE_direction_consistency) &
            Shrunk_GSE_direction_consistency >= 0.80]
  summary_file <- file.path(OUTDIR, "R2C1_59core_apeglm_directional_sensitivity_summary_v1.tsv")
  fwrite(summary, summary_file, sep="\t", quote=FALSE, na="NA")

  report <- c(
    "===== R2C1 APEGLM CORE59 DIRECTIONAL SENSITIVITY =====",
    paste("Contrasts:", uniqueN(dt$Contrast_id)),
    paste("Independent GSE:", uniqueN(dt$GSE)),
    paste("Core gene–contrast pairs passing primary QC + shrinkage:", nrow(pass)),
    paste("Core genes with >=10 paired GSEs:", sum(summary$N_GSE_paired >= 10, na.rm=TRUE)),
    paste("Core genes with shrunk GSE consistency >=0.80:",
          sum(summary$Shrunk_consistency_ge_0_80)),
    paste("Core genes with no GSE-level sign flips:",
          sum(summary$N_GSE_direction_flips == 0, na.rm=TRUE)),
    paste("Total GSE-level sign flips:", sum(summary$N_GSE_direction_flips, na.rm=TRUE)),
    paste("Detailed summary:", summary_file),
    "NOTE: This is a descriptive direction-stability analysis, not a refit of the primary multilevel meta-analysis."
  )
  report_file <- file.path(OUTDIR, "R2C1_apeglm_summary_v1.txt")
  writeLines(report, report_file)
  cat("\n", paste(report, collapse="\n"), "\n", sep="")
}
cat("\nDONE\n")
