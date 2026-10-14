#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("metafor", quietly = TRUE)) stop("Required package: metafor")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Required package: data.table")
})

library(metafor)
library(data.table)

ROOT <- normalizePath(Sys.getenv("HFD_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)

LONG <- file.path(ROOT, "00_metadata/analysis_inputs/meta_analysis_v1",
                  "HFD_primary105_meta_effects_primaryQC_v1.tsv.gz")
PREV <- file.path(ROOT, "05_meta/gene_level_v1",
                  "HFD_gene_multilevel_meta_primary_v1_1.tsv.gz")
FAIL <- file.path(ROOT, "00_metadata/analysis_outputs/gene_meta_v1_1",
                  "HFD_gene_multilevel_meta_failures_v1_1.tsv")

OUTDIR <- file.path(ROOT, "05_meta/gene_level_v1_2")
METADIR <- file.path(ROOT, "00_metadata/analysis_outputs/gene_meta_v1_2")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(METADIR, recursive=TRUE, showWarnings=FALSE)

OUT <- file.path(OUTDIR, "HFD_gene_multilevel_meta_primary_v1_2.tsv.gz")
RETRY <- file.path(METADIR, "HFD_gene_multilevel_meta_retry_status_v1_2.tsv")
SUMMARY <- file.path(METADIR, "HFD_gene_multilevel_meta_summary_v1_2.txt")

for (q in c(LONG, PREV, FAIL)) if (!file.exists(q)) stop("Missing: ", q)

prev <- fread(PREV)
fail <- fread(FAIL)

if (nrow(fail) != 23L) warning("Expected 23 failed genes; observed ", nrow(fail))
failed_genes <- fail$Geneid

need_cols <- c("Contrast_id","GSE","Geneid","log2FoldChange","lfcSE","sampling_variance")
dt <- fread(LONG, select=need_cols)
dt <- dt[Geneid %in% failed_genes]

typical_sampling_variance <- function(vi, p=1L) {
  wi <- 1/vi
  denom <- sum(wi)^2 - sum(wi^2)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  (length(vi)-p) * sum(wi) / denom
}

direction_metrics <- function(d, pooled) {
  s <- sign(pooled)
  gse_eff <- d[, .(
    GSE_effect_IVW = sum(log2FoldChange/sampling_variance) /
                     sum(1/sampling_variance)
  ), by=GSE]
  list(
    contrast_consistency = mean(sign(d$log2FoldChange) == s),
    gse_consistency = mean(sign(gse_eff$GSE_effect_IVW) == s),
    n_positive_contrasts = sum(d$log2FoldChange > 0),
    n_negative_contrasts = sum(d$log2FoldChange < 0),
    n_zero_contrasts = sum(d$log2FoldChange == 0),
    n_positive_GSE = sum(gse_eff$GSE_effect_IVW > 0),
    n_negative_GSE = sum(gse_eff$GSE_effect_IVW < 0),
    n_zero_GSE = sum(gse_eff$GSE_effect_IVW == 0)
  )
}

fit_with_optimizer <- function(d, optimizer) {
  ctrl <- list(optimizer=optimizer)
  rma.mv(
    yi=log2FoldChange,
    V=sampling_variance,
    random=list(~1|GSE, ~1|Contrast_id),
    data=d,
    method="REML",
    test="z",
    sparse=TRUE,
    control=ctrl
  )
}

make_result <- function(d, fit, optimizer_used) {
  mu <- as.numeric(coef(fit)[1])
  se_mu <- as.numeric(fit$se[1])
  tau_gse <- as.numeric(fit$sigma2[1])
  tau_contrast <- as.numeric(fit$sigma2[2])
  tau_total <- tau_gse + tau_contrast
  vt <- typical_sampling_variance(d$sampling_variance)

  denom <- tau_total + vt
  if (is.finite(denom) && denom > 0) {
    I2_gse <- 100*tau_gse/denom
    I2_contrast <- 100*tau_contrast/denom
    I2_total <- 100*tau_total/denom
  } else {
    I2_gse <- I2_contrast <- I2_total <- NA_real_
  }

  pred_se <- sqrt(se_mu^2 + tau_total)
  pred_lb <- mu - qnorm(.975)*pred_se
  pred_ub <- mu + qnorm(.975)*pred_se
  dm <- direction_metrics(d, mu)

  data.frame(
    Geneid=d$Geneid[1],
    N_contrasts=nrow(d),
    N_independent_GSE=uniqueN(d$GSE),
    pooled_log2FC=mu,
    pooled_SE=se_mu,
    CI95_lower=as.numeric(fit$ci.lb[1]),
    CI95_upper=as.numeric(fit$ci.ub[1]),
    z_value=as.numeric(fit$zval[1]),
    p_value=as.numeric(fit$pval[1]),
    tau2_GSE=tau_gse,
    tau2_contrast_within_GSE=tau_contrast,
    tau2_total=tau_total,
    typical_sampling_variance=vt,
    I2_GSE_percent=I2_gse,
    I2_contrast_within_GSE_percent=I2_contrast,
    I2_total_percent=I2_total,
    QE=as.numeric(fit$QE),
    QE_df=as.integer(fit$k-fit$p),
    QE_p_value=as.numeric(fit$QEp),
    prediction95_lower_new_GSE=pred_lb,
    prediction95_upper_new_GSE=pred_ub,
    prediction_interval_excludes_zero=ifelse(pred_lb>0 | pred_ub<0, "YES","NO"),
    contrast_direction_consistency=dm$contrast_consistency,
    GSE_balanced_direction_consistency=dm$gse_consistency,
    N_positive_contrasts=dm$n_positive_contrasts,
    N_negative_contrasts=dm$n_negative_contrasts,
    N_zero_contrasts=dm$n_zero_contrasts,
    N_positive_GSE=dm$n_positive_GSE,
    N_negative_GSE=dm$n_negative_GSE,
    N_zero_GSE=dm$n_zero_GSE,
    convergence="PASS",
    rescue_optimizer=optimizer_used,
    stringsAsFactors=FALSE
  )
}

retry_rows <- list()
rescued_rows <- list()

for (gene in failed_genes) {
  d <- dt[Geneid == gene]
  if (nrow(d) == 0L) {
    retry_rows[[length(retry_rows)+1L]] <- data.frame(
      Geneid=gene, Retry_status="FAIL", Optimizer_used="NONE",
      Error="No primary-QC effects found", stringsAsFactors=FALSE)
    next
  }

  fit <- tryCatch(fit_with_optimizer(d, "BFGS"), error=function(e)e)
  optimizer <- "BFGS"

  if (inherits(fit, "error")) {
    fit2 <- tryCatch(fit_with_optimizer(d, "Nelder-Mead"), error=function(e)e)
    if (!inherits(fit2, "error")) {
      fit <- fit2
      optimizer <- "Nelder-Mead"
    }
  }

  if (inherits(fit, "error")) {
    retry_rows[[length(retry_rows)+1L]] <- data.frame(
      Geneid=gene, Retry_status="FAIL", Optimizer_used=optimizer,
      Error=conditionMessage(fit), stringsAsFactors=FALSE)
  } else {
    rescued_rows[[length(rescued_rows)+1L]] <- make_result(d, fit, optimizer)
    retry_rows[[length(retry_rows)+1L]] <- data.frame(
      Geneid=gene, Retry_status="PASS", Optimizer_used=optimizer,
      Error="", stringsAsFactors=FALSE)
  }
}

retry <- rbindlist(retry_rows, fill=TRUE)
fwrite(retry, RETRY, sep="\t", quote=FALSE, na="NA")

rescued <- if (length(rescued_rows)) rbindlist(rescued_rows, fill=TRUE) else data.table()

# Preserve previous successful fits. Add a provenance column.
if (!"rescue_optimizer" %in% names(prev)) prev[, rescue_optimizer := "DEFAULT_nlminb"]

# Align columns before append.
all_cols <- union(names(prev), names(rescued))
for (cc in setdiff(all_cols, names(prev))) prev[, (cc) := NA]
for (cc in setdiff(all_cols, names(rescued))) rescued[, (cc) := NA]
prev <- prev[, ..all_cols]
rescued <- rescued[, ..all_cols]

final <- rbind(prev, rescued, fill=TRUE)

# Recompute multiplicity correction after rescuing any failed models.
final[, pooled_FDR_BH := p.adjust(p_value, method="BH")]
final[, pooled_significant_FDR05 := ifelse(pooled_FDR_BH < .05, "YES","NO")]
final[, stable_core_preLOO := ifelse(
  pooled_FDR_BH < .05 &
  prediction_interval_excludes_zero == "YES" &
  GSE_balanced_direction_consistency >= .80,
  "YES","NO"
)]

final[, abs_pooled_log2FC_sort := abs(pooled_log2FC)]
setorder(final, pooled_FDR_BH, -abs_pooled_log2FC_sort, Geneid)
final[, abs_pooled_log2FC_sort := NULL]

fwrite(final, OUT, sep="\t", quote=FALSE, na="NA", compress="gzip")

n_rescued <- sum(retry$Retry_status == "PASS")
n_still_fail <- sum(retry$Retry_status == "FAIL")

summary <- c(
  "===== HFD GENE MULTILEVEL META-ANALYSIS FINALIZATION v1.2 =====",
  "",
  paste0("FINAL_STATUS=", ifelse(n_still_fail==0L, "PASS", "PASS_WITH_REMAINING_FAILURES")),
  paste0("Previous_successful_models=", nrow(prev)),
  paste0("Previous_failed_models=", nrow(fail)),
  paste0("Retry_rescued=", n_rescued),
  paste0("Retry_still_failed=", n_still_fail),
  paste0("Final_successful_models=", nrow(final)),
  "",
  paste0("Pooled_FDR_lt_0.05=", sum(final$pooled_FDR_BH < .05, na.rm=TRUE)),
  paste0("Prediction_interval_excludes_zero=",
         sum(final$prediction_interval_excludes_zero=="YES", na.rm=TRUE)),
  paste0("GSE_direction_consistency_ge_0.80=",
         sum(final$GSE_balanced_direction_consistency >= .80, na.rm=TRUE)),
  paste0("Stable_core_preLOO=",
         sum(final$stable_core_preLOO=="YES", na.rm=TRUE)),
  "",
  paste0("Median_I2_total_percent=",
         format(median(final$I2_total_percent, na.rm=TRUE), digits=6)),
  paste0("Median_tau2_total=",
         format(median(final$tau2_total, na.rm=TRUE), digits=6)),
  "",
  paste0("Final_primary_result=", sub(paste0("^",ROOT,"/"),"",OUT)),
  paste0("Retry_status=", sub(paste0("^",ROOT,"/"),"",RETRY)),
  "",
  "PRIMARY_MODEL_STRUCTURE_CHANGED=NO",
  "ONLY_PREVIOUS_OPTIMIZER_FAILURES_RETRIED=YES",
  "BH_FDR_RECOMPUTED_AFTER_RESCUE=YES",
  "LEAVE_ONE_GSE_OUT_PERFORMED=NO"
)

writeLines(summary, SUMMARY)
cat(paste(summary, collapse="\n"), "\n")

