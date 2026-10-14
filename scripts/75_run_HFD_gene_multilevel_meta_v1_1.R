#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Required R package not installed: metafor")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Required R package not installed: data.table")
  }
})

library(metafor)
library(data.table)

ROOT <- normalizePath(Sys.getenv("HFD_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)

LONG <- file.path(
  ROOT,
  "00_metadata/analysis_inputs/meta_analysis_v1",
  "HFD_primary105_meta_effects_primaryQC_v1.tsv.gz"
)

ELIG <- file.path(
  ROOT,
  "00_metadata/analysis_inputs/meta_analysis_v1",
  "HFD_gene_meta_eligible_ge10GSE_v1.tsv"
)

OUTDIR <- file.path(ROOT, "05_meta", "gene_level_v1")
METADIR <- file.path(
  ROOT,
  "00_metadata/analysis_outputs/gene_meta_v1_1"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)

OUT <- file.path(OUTDIR, "HFD_gene_multilevel_meta_primary_v1_1.tsv.gz")
FAIL <- file.path(METADIR, "HFD_gene_multilevel_meta_failures_v1_1.tsv")
SUMMARY <- file.path(METADIR, "HFD_gene_multilevel_meta_summary_v1_1.txt")

for (q in c(LONG, ELIG)) {
  if (!file.exists(q)) stop("Missing required input: ", q)
}

cat("===== HFD GENE MULTILEVEL META-ANALYSIS v1.1 =====\n")
cat("R:", R.version.string, "\n")
cat("metafor:", as.character(packageVersion("metafor")), "\n")
cat("data.table:", as.character(packageVersion("data.table")), "\n\n")

# Read only columns needed for the primary intercept-only model.
need_cols <- c(
  "Contrast_id", "GSE", "Geneid",
  "log2FoldChange", "lfcSE", "sampling_variance"
)

dt <- fread(
  LONG,
  sep = "\t",
  select = need_cols,
  na.strings = c("NA", "")
)

elig <- fread(ELIG, sep = "\t")
eligible_genes <- elig[
  Primary_meta_eligible_ge10_independent_GSE == "YES",
  Geneid
]

if (length(eligible_genes) != 18234L) {
  warning(
    "Expected 18,234 eligible genes from audit, observed ",
    length(eligible_genes)
  )
}

dt <- dt[Geneid %in% eligible_genes]

if (any(!is.finite(dt$log2FoldChange))) stop("Non-finite log2FC in primary input")
if (any(!is.finite(dt$sampling_variance))) stop("Non-finite variance in primary input")
if (any(dt$sampling_variance <= 0)) stop("Non-positive variance in primary input")

setkey(dt, Geneid)
gene_split <- split(dt, by = "Geneid", keep.by = TRUE, drop = TRUE)

typical_sampling_variance <- function(vi, p = 1L) {
  wi <- 1 / vi
  denom <- sum(wi)^2 - sum(wi^2)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  (length(vi) - p) * sum(wi) / denom
}

direction_metrics <- function(d, pooled) {
  s <- sign(pooled)

  # Contrast-level directional consistency.
  contrast_sign <- sign(d$log2FoldChange)
  contrast_cons <- mean(contrast_sign == s)

  # GSE-balanced directional consistency:
  # one inverse-variance weighted summary per independent GSE.
  gse_eff <- d[, .(
    GSE_effect_IVW = sum(log2FoldChange / sampling_variance) /
                     sum(1 / sampling_variance)
  ), by = GSE]

  gse_cons <- mean(sign(gse_eff$GSE_effect_IVW) == s)

  list(
    contrast_consistency = contrast_cons,
    gse_consistency = gse_cons,
    n_positive_contrasts = sum(d$log2FoldChange > 0),
    n_negative_contrasts = sum(d$log2FoldChange < 0),
    n_zero_contrasts = sum(d$log2FoldChange == 0),
    n_positive_GSE = sum(gse_eff$GSE_effect_IVW > 0),
    n_negative_GSE = sum(gse_eff$GSE_effect_IVW < 0),
    n_zero_GSE = sum(gse_eff$GSE_effect_IVW == 0)
  )
}

fit_one <- function(d) {
  gene <- d$Geneid[1L]

  k_contrast <- nrow(d)
  k_gse <- uniqueN(d$GSE)

  if (k_gse < 10L) {
    return(list(
      ok = FALSE, gene = gene,
      error = paste0("Independent GSE <10: ", k_gse)
    ))
  }

  # Unique contrast identifiers are retained as the lower random-effect level.
  # Because Contrast_id is globally unique, this estimates residual
  # between-contrast heterogeneity beyond known sampling variance.
  fit <- tryCatch(
    rma.mv(
      yi = log2FoldChange,
      V = sampling_variance,
      random = list(~ 1 | GSE, ~ 1 | Contrast_id),
      data = d,
      method = "REML",
      test = "z",
      sparse = TRUE
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(ok = FALSE, gene = gene, error = conditionMessage(fit)))
  }

  mu <- as.numeric(coef(fit)[1])
  se_mu <- as.numeric(fit$se[1])
  pval <- as.numeric(fit$pval[1])
  ci_lb <- as.numeric(fit$ci.lb[1])
  ci_ub <- as.numeric(fit$ci.ub[1])

  if (length(fit$sigma2) != 2L) {
    return(list(
      ok = FALSE, gene = gene,
      error = paste0("Expected 2 sigma2 components; got ", length(fit$sigma2))
    ))
  }

  tau_gse <- as.numeric(fit$sigma2[1])
  tau_contrast <- as.numeric(fit$sigma2[2])
  tau_total <- tau_gse + tau_contrast

  vt <- typical_sampling_variance(d$sampling_variance, p = 1L)

  denom_i2 <- tau_total + vt
  if (is.finite(denom_i2) && denom_i2 > 0) {
    I2_gse <- 100 * tau_gse / denom_i2
    I2_contrast <- 100 * tau_contrast / denom_i2
    I2_total <- 100 * tau_total / denom_i2
  } else {
    I2_gse <- NA_real_
    I2_contrast <- NA_real_
    I2_total <- NA_real_
  }

  # Prediction interval for a new contrast from a new independent GSE:
  # pooled-mean uncertainty + both estimated heterogeneity components.
  pred_se <- sqrt(se_mu^2 + tau_total)
  pred_lb <- mu - qnorm(0.975) * pred_se
  pred_ub <- mu + qnorm(0.975) * pred_se

  dm <- direction_metrics(d, mu)

  out <- data.frame(
    Geneid = gene,
    N_contrasts = k_contrast,
    N_independent_GSE = k_gse,
    pooled_log2FC = mu,
    pooled_SE = se_mu,
    CI95_lower = ci_lb,
    CI95_upper = ci_ub,
    z_value = as.numeric(fit$zval[1]),
    p_value = pval,
    tau2_GSE = tau_gse,
    tau2_contrast_within_GSE = tau_contrast,
    tau2_total = tau_total,
    typical_sampling_variance = vt,
    I2_GSE_percent = I2_gse,
    I2_contrast_within_GSE_percent = I2_contrast,
    I2_total_percent = I2_total,
    QE = as.numeric(fit$QE),
    QE_df = as.integer(fit$k - fit$p),
    QE_p_value = as.numeric(fit$QEp),
    prediction95_lower_new_GSE = pred_lb,
    prediction95_upper_new_GSE = pred_ub,
    prediction_interval_excludes_zero = ifelse(
      pred_lb > 0 | pred_ub < 0, "YES", "NO"
    ),
    contrast_direction_consistency = dm$contrast_consistency,
    GSE_balanced_direction_consistency = dm$gse_consistency,
    N_positive_contrasts = dm$n_positive_contrasts,
    N_negative_contrasts = dm$n_negative_contrasts,
    N_zero_contrasts = dm$n_zero_contrasts,
    N_positive_GSE = dm$n_positive_GSE,
    N_negative_GSE = dm$n_negative_GSE,
    N_zero_GSE = dm$n_zero_GSE,
    convergence = ifelse(is.null(fit$converged) || isTRUE(fit$converged),
                         "PASS", "REVIEW"),
    stringsAsFactors = FALSE
  )

  list(ok = TRUE, gene = gene, result = out, error = "")
}

# Parallel execution across genes. On Linux, mclapply uses forked workers.
# Default cores are deliberately conservative; override with HFD_META_CORES.
cores_env <- Sys.getenv("HFD_META_CORES", unset = "")
if (nzchar(cores_env)) {
  NCORES <- as.integer(cores_env)
} else {
  NCORES <- max(1L, min(8L, parallel::detectCores(logical = FALSE) - 1L))
}
if (!is.finite(NCORES) || NCORES < 1L) NCORES <- 1L

cat("Eligible genes:", length(gene_split), "\n")
cat("Parallel workers:", NCORES, "\n\n")

fits <- parallel::mclapply(
  gene_split,
  fit_one,
  mc.cores = NCORES,
  mc.preschedule = TRUE
)

ok <- vapply(fits, function(z) isTRUE(z$ok), logical(1))

fail_df <- if (any(!ok)) {
  rbindlist(lapply(fits[!ok], function(z) {
    data.frame(Geneid = z$gene, Error = z$error, stringsAsFactors = FALSE)
  }), fill = TRUE)
} else {
  data.table(Geneid = character(), Error = character())
}

fwrite(fail_df, FAIL, sep = "\t", quote = FALSE, na = "NA")

if (!any(ok)) {
  stop("All gene-level meta-analysis models failed")
}

res <- rbindlist(lapply(fits[ok], function(z) z$result), fill = TRUE)

# Gene-wise multiplicity correction for the primary pooled effect.
res[, pooled_FDR_BH := p.adjust(p_value, method = "BH")]
res[, pooled_significant_FDR05 := ifelse(pooled_FDR_BH < 0.05, "YES", "NO")]

# Stable-core components available at this stage.
# Leave-one-GSE-out stability is intentionally added in the next stage.
res[, stable_core_preLOO := ifelse(
  pooled_FDR_BH < 0.05 &
    prediction_interval_excludes_zero == "YES" &
    GSE_balanced_direction_consistency >= 0.80,
  "YES", "NO"
)]

res[, abs_pooled_log2FC_sort := abs(pooled_log2FC)]
setorder(res, pooled_FDR_BH, -abs_pooled_log2FC_sort, Geneid)
res[, abs_pooled_log2FC_sort := NULL]

# Atomic gzip write.
tmp <- paste0(OUT, ".tmp")
fwrite(res, tmp, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
if (file.exists(OUT)) file.remove(OUT)
if (!file.rename(tmp, OUT)) stop("Failed to rename primary meta result")

n_fail <- nrow(fail_df)
n_fit <- nrow(res)

summary <- c(
  "===== HFD GENE MULTILEVEL META-ANALYSIS SUMMARY v1.1 =====",
  "",
  paste0("FINAL_STATUS=", ifelse(n_fail == 0L, "PASS", "PASS_WITH_MODEL_FAILURES")),
  paste0("Eligible_genes_input=", length(eligible_genes)),
  paste0("Genes_successfully_fit=", n_fit),
  paste0("Genes_failed=", n_fail),
  paste0("Independent_GSE_total=64"),
  paste0("Contrasts_total=105"),
  "",
  paste0("R_version=", R.version.string),
  paste0("metafor_version=", as.character(packageVersion("metafor"))),
  paste0("Parallel_workers=", NCORES),
  "",
  "MODEL=rma.mv_REML_random_GSE_plus_Contrast",
  "Effect_size=DESeq2_unshrunken_log2FC",
  "Sampling_variance=lfcSE_squared",
  "GSE_dependence_modeled=YES",
  "Shared_GSM_covariance_required=NO",
  "",
  paste0("Pooled_FDR_lt_0.05=", sum(res$pooled_FDR_BH < 0.05, na.rm = TRUE)),
  paste0("Prediction_interval_excludes_zero=",
         sum(res$prediction_interval_excludes_zero == "YES", na.rm = TRUE)),
  paste0("GSE_direction_consistency_ge_0.80=",
         sum(res$GSE_balanced_direction_consistency >= 0.80, na.rm = TRUE)),
  paste0("Stable_core_preLOO=",
         sum(res$stable_core_preLOO == "YES", na.rm = TRUE)),
  "",
  paste0("Median_I2_total_percent=",
         format(median(res$I2_total_percent, na.rm = TRUE), digits = 6)),
  paste0("Median_tau2_total=",
         format(median(res$tau2_total, na.rm = TRUE), digits = 6)),
  "",
  paste0("Primary_result=", sub(paste0("^", ROOT, "/"), "", OUT)),
  paste0("Failure_table=", sub(paste0("^", ROOT, "/"), "", FAIL)),
  "",
  "LEAVE_ONE_GSE_OUT_PERFORMED=NO",
  "MODERATOR_ANALYSIS_PERFORMED=NO",
  "PATHWAY_ANALYSIS_PERFORMED=NO"
)

writeLines(summary, SUMMARY)
cat("\n", paste(summary, collapse = "\n"), "\n", sep = "")

if (n_fail > 0L) {
  quit(status = 2L)
}

