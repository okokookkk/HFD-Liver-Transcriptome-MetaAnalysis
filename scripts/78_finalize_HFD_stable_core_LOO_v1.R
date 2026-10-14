#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("metafor", quietly = TRUE)) stop("Required package: metafor")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Required package: data.table")
})

library(metafor)
library(data.table)

ROOT <- normalizePath(Sys.getenv("HFD_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)

META <- file.path(
  ROOT, "05_meta/gene_level_v1_2",
  "HFD_gene_multilevel_meta_primary_v1_2.tsv.gz"
)

LONG <- file.path(
  ROOT, "00_metadata/analysis_inputs/meta_analysis_v1",
  "HFD_primary105_meta_effects_primaryQC_v1.tsv.gz"
)

OUTDIR <- file.path(ROOT, "05_meta", "stable_core_v1")
METADIR <- file.path(
  ROOT, "00_metadata/analysis_outputs/stable_core_v1"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)

OUT_DETAIL <- file.path(
  OUTDIR,
  "HFD_stable_core_preLOO_leave_one_GSE_out_detail_v1.tsv.gz"
)
OUT_FINAL <- file.path(
  OUTDIR,
  "HFD_stable_core_final_v1.tsv"
)
OUT_STATUS <- file.path(
  METADIR,
  "HFD_stable_core_LOO_status_v1.tsv"
)
OUT_SUMMARY <- file.path(
  METADIR,
  "HFD_stable_core_LOO_summary_v1.txt"
)

for (q in c(META, LONG)) {
  if (!file.exists(q)) stop("Missing required input: ", q)
}

meta <- fread(META)
cand <- meta[stable_core_preLOO == "YES"]

if (nrow(cand) != 59L) {
  warning("Expected 59 preLOO candidates; observed ", nrow(cand))
}

candidate_genes <- cand$Geneid

need_cols <- c(
  "Contrast_id", "GSE", "Geneid",
  "log2FoldChange", "sampling_variance"
)

dt <- fread(
  LONG,
  select = need_cols,
  na.strings = c("NA", "")
)
dt <- dt[Geneid %in% candidate_genes]

if (any(!is.finite(dt$log2FoldChange))) stop("Non-finite effect in candidate input")
if (any(!is.finite(dt$sampling_variance))) stop("Non-finite variance in candidate input")
if (any(dt$sampling_variance <= 0)) stop("Non-positive variance in candidate input")

fit_meta <- function(d, optimizer = NULL) {
  args <- list(
    yi = d$log2FoldChange,
    V = d$sampling_variance,
    random = list(~ 1 | GSE, ~ 1 | Contrast_id),
    data = d,
    method = "REML",
    test = "z",
    sparse = TRUE
  )
  if (!is.null(optimizer)) {
    args$control <- list(optimizer = optimizer)
  }
  do.call(rma.mv, args)
}

fit_with_rescue <- function(d) {
  fit <- tryCatch(fit_meta(d), error = function(e) e)
  opt <- "DEFAULT_nlminb"

  if (inherits(fit, "error")) {
    fit <- tryCatch(fit_meta(d, "BFGS"), error = function(e) e)
    opt <- "BFGS"
  }

  if (inherits(fit, "error")) {
    fit <- tryCatch(fit_meta(d, "Nelder-Mead"), error = function(e) e)
    opt <- "Nelder-Mead"
  }

  list(fit = fit, optimizer = opt)
}

loo_rows <- list()
status_rows <- list()

cat("===== HFD STABLE-CORE LEAVE-ONE-GSE-OUT v1 =====\n")
cat("Candidates:", nrow(cand), "\n\n")

for (ii in seq_len(nrow(cand))) {

  gene <- cand$Geneid[ii]
  pooled_full <- cand$pooled_log2FC[ii]
  full_sign <- sign(pooled_full)

  d <- dt[Geneid == gene]
  gses <- sort(unique(d$GSE))

  if (length(gses) < 10L) {
    stop(gene, ": candidate unexpectedly has <10 independent GSE")
  }

  cat(sprintf("[%02d/%02d] %s : %d GSE\n", ii, nrow(cand), gene, length(gses)))

  gene_fail <- 0L
  gene_sign_flip <- 0L
  gene_zero <- 0L

  gene_effects <- numeric()
  gene_ses <- numeric()
  gene_pvals <- numeric()

  for (gse_drop in gses) {

    z <- d[GSE != gse_drop]

    rr <- fit_with_rescue(z)
    fit <- rr$fit

    if (inherits(fit, "error")) {

      gene_fail <- gene_fail + 1L

      loo_rows[[length(loo_rows) + 1L]] <- data.frame(
        Geneid = gene,
        Full_pooled_log2FC = pooled_full,
        Dropped_GSE = gse_drop,
        Remaining_GSE = uniqueN(z$GSE),
        Remaining_contrasts = nrow(z),
        LOO_status = "FAIL",
        Optimizer_used = rr$optimizer,
        LOO_pooled_log2FC = NA_real_,
        LOO_SE = NA_real_,
        LOO_CI95_lower = NA_real_,
        LOO_CI95_upper = NA_real_,
        LOO_p_value = NA_real_,
        LOO_tau2_GSE = NA_real_,
        LOO_tau2_contrast_within_GSE = NA_real_,
        LOO_tau2_total = NA_real_,
        Direction_same_as_full = NA_character_,
        Error = conditionMessage(fit),
        stringsAsFactors = FALSE
      )

      next
    }

    mu <- as.numeric(coef(fit)[1])
    se <- as.numeric(fit$se[1])
    pval <- as.numeric(fit$pval[1])

    tau_gse <- as.numeric(fit$sigma2[1])
    tau_contrast <- as.numeric(fit$sigma2[2])
    tau_total <- tau_gse + tau_contrast

    same <- sign(mu) == full_sign

    if (!same) gene_sign_flip <- gene_sign_flip + 1L
    if (sign(mu) == 0) gene_zero <- gene_zero + 1L

    gene_effects <- c(gene_effects, mu)
    gene_ses <- c(gene_ses, se)
    gene_pvals <- c(gene_pvals, pval)

    loo_rows[[length(loo_rows) + 1L]] <- data.frame(
      Geneid = gene,
      Full_pooled_log2FC = pooled_full,
      Dropped_GSE = gse_drop,
      Remaining_GSE = uniqueN(z$GSE),
      Remaining_contrasts = nrow(z),
      LOO_status = "PASS",
      Optimizer_used = rr$optimizer,
      LOO_pooled_log2FC = mu,
      LOO_SE = se,
      LOO_CI95_lower = as.numeric(fit$ci.lb[1]),
      LOO_CI95_upper = as.numeric(fit$ci.ub[1]),
      LOO_p_value = pval,
      LOO_tau2_GSE = tau_gse,
      LOO_tau2_contrast_within_GSE = tau_contrast,
      LOO_tau2_total = tau_total,
      Direction_same_as_full = ifelse(same, "YES", "NO"),
      Error = "",
      stringsAsFactors = FALSE
    )
  }

  if (length(gene_effects)) {
    min_effect <- min(gene_effects)
    max_effect <- max(gene_effects)
    max_abs_change <- max(abs(gene_effects - pooled_full))
    min_abs_effect <- min(abs(gene_effects))
    max_p <- max(gene_pvals)
  } else {
    min_effect <- max_effect <- max_abs_change <- min_abs_effect <- max_p <- NA_real_
  }

  final_stable <- (
    gene_fail == 0L &&
    gene_sign_flip == 0L &&
    gene_zero == 0L
  )

  status_rows[[length(status_rows) + 1L]] <- data.frame(
    Geneid = gene,
    Full_pooled_log2FC = pooled_full,
    Full_pooled_FDR_BH = cand$pooled_FDR_BH[ii],
    Full_prediction95_lower = cand$prediction95_lower_new_GSE[ii],
    Full_prediction95_upper = cand$prediction95_upper_new_GSE[ii],
    Full_GSE_direction_consistency = cand$GSE_balanced_direction_consistency[ii],
    N_full_independent_GSE = cand$N_independent_GSE[ii],
    N_LOO_fits_expected = length(gses),
    N_LOO_fits_PASS = length(gses) - gene_fail,
    N_LOO_fits_FAIL = gene_fail,
    N_LOO_direction_flip = gene_sign_flip,
    N_LOO_zero_effect = gene_zero,
    LOO_min_pooled_log2FC = min_effect,
    LOO_max_pooled_log2FC = max_effect,
    LOO_min_abs_pooled_log2FC = min_abs_effect,
    LOO_max_abs_change_from_full = max_abs_change,
    LOO_max_nominal_p_value = max_p,
    LOO_direction_stable = ifelse(
      gene_fail == 0L && gene_sign_flip == 0L && gene_zero == 0L,
      "YES", "NO"
    ),
    Final_stable_core = ifelse(final_stable, "YES", "NO"),
    stringsAsFactors = FALSE
  )
}

detail <- rbindlist(loo_rows, fill = TRUE)
status <- rbindlist(status_rows, fill = TRUE)

fwrite(
  detail,
  OUT_DETAIL,
  sep = "\t",
  quote = FALSE,
  na = "NA",
  compress = "gzip"
)

fwrite(
  status,
  OUT_STATUS,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# Final stable-core table retains the full primary meta-analysis statistics
# and adds the LOO robustness statistics.
final <- merge(
  cand,
  status,
  by = "Geneid",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("", "_LOO")
)

final <- final[Final_stable_core == "YES"]

# Stable output ordering: direction first, then pooled FDR/effect.
final[, abs_effect_sort := abs(pooled_log2FC)]
setorder(final, -pooled_log2FC, pooled_FDR_BH, -abs_effect_sort, Geneid)
final[, abs_effect_sort := NULL]

fwrite(
  final,
  OUT_FINAL,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

n_final <- nrow(final)
n_fail_candidates <- sum(status$N_LOO_fits_FAIL > 0)
n_flip_candidates <- sum(status$N_LOO_direction_flip > 0)
n_direction_stable <- sum(status$LOO_direction_stable == "YES")

summary <- c(
  "===== HFD FINAL STABLE-CORE LOO SUMMARY v1 =====",
  "",
  "STATUS=PASS",
  paste0("PreLOO_candidates=", nrow(cand)),
  paste0("Candidates_with_any_LOO_fit_failure=", n_fail_candidates),
  paste0("Candidates_with_any_direction_flip=", n_flip_candidates),
  paste0("Candidates_LOO_direction_stable=", n_direction_stable),
  paste0("Final_stable_core_genes=", n_final),
  "",
  "FINAL_STABLE_CORE_CRITERIA:",
  "1_pooled_BH_FDR_lt_0.05=YES",
  "2_prediction_interval_excludes_zero=YES",
  "3_GSE_balanced_direction_consistency_ge_0.80=YES",
  "4_leave_one_GSE_out_direction_stable=YES",
  "",
  paste0(
    "Final_stable_core_upregulated=",
    sum(final$pooled_log2FC > 0, na.rm = TRUE)
  ),
  paste0(
    "Final_stable_core_downregulated=",
    sum(final$pooled_log2FC < 0, na.rm = TRUE)
  ),
  "",
  paste0("LOO_detail=", sub(paste0("^", ROOT, "/"), "", OUT_DETAIL)),
  paste0("LOO_status=", sub(paste0("^", ROOT, "/"), "", OUT_STATUS)),
  paste0("Final_stable_core=", sub(paste0("^", ROOT, "/"), "", OUT_FINAL)),
  "",
  "MODERATOR_ANALYSIS_PERFORMED=NO",
  "PATHWAY_ANALYSIS_PERFORMED=NO"
)

writeLines(summary, OUT_SUMMARY)
cat("\n", paste(summary, collapse = "\n"), "\n", sep = "")

