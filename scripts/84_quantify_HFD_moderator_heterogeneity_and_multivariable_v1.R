#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("metafor", quietly = TRUE)) stop("Required package: metafor")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Required package: data.table")
})
library(metafor)
library(data.table)

ROOT <- normalizePath(Sys.getenv("HFD_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)

LONG <- file.path(
  ROOT, "00_metadata/analysis_inputs/meta_analysis_v1",
  "HFD_primary105_meta_effects_primaryQC_v1.tsv.gz"
)

DUR <- file.path(
  ROOT, "05_meta/moderators_v1_2",
  "HFD_duration_meta_regression_v1_2.tsv.gz"
)
FAT <- file.path(
  ROOT, "05_meta/moderators_v1_2",
  "HFD_fat_percent_meta_regression_v1_2.tsv.gz"
)
CTRL <- file.path(
  ROOT, "05_meta/moderators_v1_2",
  "HFD_control_diet_meta_regression_v1_2.tsv.gz"
)

OUTDIR <- file.path(ROOT, "05_meta", "moderator_followup_v1")
METADIR <- file.path(
  ROOT, "00_metadata/analysis_outputs/moderator_followup_v1"
)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)

OUT_DUR_R2 <- file.path(
  OUTDIR, "HFD_duration_matched_subset_heterogeneity_reduction_v1.tsv.gz"
)
OUT_FAT_R2 <- file.path(
  OUTDIR, "HFD_fat_percent_matched_subset_heterogeneity_reduction_v1.tsv.gz"
)
OUT_CTRL_R2 <- file.path(
  OUTDIR, "HFD_control_diet_matched_subset_heterogeneity_reduction_v1.tsv.gz"
)
OUT_MULTI <- file.path(
  OUTDIR, "HFD_duration_plus_fat_multivariable_meta_regression_v1.tsv.gz"
)
OUT_MULTI_FAIL <- file.path(
  METADIR, "HFD_duration_plus_fat_multivariable_failures_v1.tsv"
)
OUT_NULL_FAIL <- file.path(
  METADIR, "HFD_matched_null_model_failures_v1.tsv"
)
OUT_SUMMARY <- file.path(
  METADIR, "HFD_moderator_heterogeneity_and_multivariable_summary_v1.txt"
)

for (q in c(LONG, DUR, FAT, CTRL)) {
  if (!file.exists(q)) stop("Missing required input: ", q)
}

need <- c(
  "Contrast_id", "GSE", "Geneid",
  "log2FoldChange", "sampling_variance",
  "hfd_duration_weeks", "hfd_fat_kcal_percent",
  "control_diet_type"
)
dt <- fread(LONG, select = need, na.strings = c("NA", ""))

dt[, log2FoldChange := as.numeric(log2FoldChange)]
dt[, sampling_variance := as.numeric(sampling_variance)]
dt[, hfd_duration_weeks := suppressWarnings(as.numeric(hfd_duration_weeks))]
dt[, hfd_fat_kcal_percent := suppressWarnings(as.numeric(hfd_fat_kcal_percent))]
dt[, control_diet_type := trimws(as.character(control_diet_type))]
dt[control_diet_type == "", control_diet_type := NA_character_]

if (any(!is.finite(dt$log2FoldChange))) stop("Non-finite effect")
if (any(!is.finite(dt$sampling_variance)) || any(dt$sampling_variance <= 0))
  stop("Invalid sampling variance")

dt[, duration_per4 := (hfd_duration_weeks - 12) / 4]
dt[, fat_per10 := (hfd_fat_kcal_percent - 60) / 10]

dur <- fread(DUR)
fat <- fread(FAT)
ctrl <- fread(CTRL)

if (!nrow(dur) || !nrow(fat) || !nrow(ctrl))
  stop("One or more primary moderator result tables are empty")

cores_env <- Sys.getenv("HFD_META_CORES", unset = "")
if (nzchar(cores_env)) {
  NCORES <- as.integer(cores_env)
} else {
  NCORES <- max(1L, min(8L, parallel::detectCores(logical = FALSE) - 1L))
}
if (!is.finite(NCORES) || NCORES < 1L) NCORES <- 1L

cat("===== HFD MODERATOR HETEROGENEITY / MULTIVARIABLE FOLLOW-UP v1 =====\n")
cat("metafor:", as.character(packageVersion("metafor")), "\n")
cat("Workers:", NCORES, "\n\n")

fit_null_rescue <- function(d) {
  dofit <- function(opt = NULL) {
    if (is.null(opt)) {
      rma.mv(
        yi = d$log2FoldChange,
        V = d$sampling_variance,
        random = list(~1 | GSE, ~1 | Contrast_id),
        data = d,
        method = "REML",
        test = "z",
        sparse = TRUE
      )
    } else {
      rma.mv(
        yi = d$log2FoldChange,
        V = d$sampling_variance,
        random = list(~1 | GSE, ~1 | Contrast_id),
        data = d,
        method = "REML",
        test = "z",
        sparse = TRUE,
        control = list(optimizer = opt)
      )
    }
  }

  fit <- tryCatch(dofit(), error = function(e) e)
  opt <- "DEFAULT_nlminb"
  if (inherits(fit, "error")) {
    fit <- tryCatch(dofit("BFGS"), error = function(e) e)
    opt <- "BFGS"
  }
  if (inherits(fit, "error")) {
    fit <- tryCatch(dofit("Nelder-Mead"), error = function(e) e)
    opt <- "Nelder-Mead"
  }
  list(fit = fit, optimizer = opt)
}

fit_multi_rescue <- function(d) {
  dofit <- function(opt = NULL) {
    if (is.null(opt)) {
      rma.mv(
        yi = d$log2FoldChange,
        V = d$sampling_variance,
        mods = ~ duration_per4 + fat_per10,
        random = list(~1 | GSE, ~1 | Contrast_id),
        data = d,
        method = "REML",
        test = "z",
        sparse = TRUE
      )
    } else {
      rma.mv(
        yi = d$log2FoldChange,
        V = d$sampling_variance,
        mods = ~ duration_per4 + fat_per10,
        random = list(~1 | GSE, ~1 | Contrast_id),
        data = d,
        method = "REML",
        test = "z",
        sparse = TRUE,
        control = list(optimizer = opt)
      )
    }
  }

  fit <- tryCatch(dofit(), error = function(e) e)
  opt <- "DEFAULT_nlminb"
  if (inherits(fit, "error")) {
    fit <- tryCatch(dofit("BFGS"), error = function(e) e)
    opt <- "BFGS"
  }
  if (inherits(fit, "error")) {
    fit <- tryCatch(dofit("Nelder-Mead"), error = function(e) e)
    opt <- "Nelder-Mead"
  }
  list(fit = fit, optimizer = opt)
}

gene_split <- split(dt, by = "Geneid", keep.by = TRUE, drop = TRUE)

subset_for <- function(d, moderator) {
  if (moderator == "duration") {
    return(d[is.finite(hfd_duration_weeks)])
  }
  if (moderator == "fat") {
    return(d[is.finite(hfd_fat_kcal_percent)])
  }
  if (moderator == "control") {
    return(d[
      !is.na(control_diet_type) &
      control_diet_type %in% c("CHOW", "LFD", "NORMAL_DIET")
    ])
  }
  stop("Unknown moderator")
}

matched_null_one <- function(gene, moderator, mod_row) {
  d <- gene_split[[gene]]
  if (is.null(d)) {
    return(list(
      status = "FAIL", Geneid = gene,
      reason = "Gene absent from primary long table"
    ))
  }

  x <- subset_for(d, moderator)
  rr <- fit_null_rescue(x)
  fit <- rr$fit

  if (inherits(fit, "error")) {
    return(list(
      status = "FAIL", Geneid = gene,
      reason = conditionMessage(fit)
    ))
  }

  tau_null_gse <- as.numeric(fit$sigma2[1])
  tau_null_contrast <- as.numeric(fit$sigma2[2])
  tau_null_total <- tau_null_gse + tau_null_contrast

  tau_mod_total <- as.numeric(mod_row$tau2_total)

  r2 <- if (is.finite(tau_null_total) && tau_null_total > 0) {
    100 * (tau_null_total - tau_mod_total) / tau_null_total
  } else {
    NA_real_
  }

  out <- data.frame(
    Geneid = gene,
    N_contrasts_matched = nrow(x),
    N_independent_GSE_matched = uniqueN(x$GSE),
    tau2_null_GSE = tau_null_gse,
    tau2_null_contrast_within_GSE = tau_null_contrast,
    tau2_null_total = tau_null_total,
    tau2_moderator_total = tau_mod_total,
    tau2_absolute_reduction = tau_null_total - tau_mod_total,
    pseudo_R2_total_percent = r2,
    null_optimizer = rr$optimizer,
    stringsAsFactors = FALSE
  )
  list(status = "PASS", result = out, Geneid = gene, reason = "")
}

run_matched_null <- function(modres, moderator, label) {
  cat("Matched-null heterogeneity audit:", label, "\n")

  genes <- modres$Geneid

  # One-gene preflight.
  pf_gene <- genes[1]
  pf <- matched_null_one(
    pf_gene, moderator,
    modres[Geneid == pf_gene][1]
  )
  if (pf$status != "PASS")
    stop("PRECHECK_FAIL matched null ", label, ": ", pf$reason)

  ans <- parallel::mclapply(
    genes,
    function(g) matched_null_one(
      g, moderator,
      modres[Geneid == g][1]
    ),
    mc.cores = NCORES,
    mc.preschedule = TRUE
  )

  ok <- vapply(ans, function(z) z$status == "PASS", logical(1))
  fail <- !ok

  r <- rbindlist(lapply(ans[ok], `[[`, "result"), fill = TRUE)

  # Merge primary moderator statistics.
  keep <- intersect(
    c(
      "Geneid", "beta", "SE", "z_value", "p_value", "FDR_BH",
      "QM_p_value", "omnibus_FDR_BH",
      "beta_LFD_vs_CHOW", "p_LFD_vs_CHOW", "LFD_vs_CHOW_FDR_BH",
      "beta_NORMAL_DIET_vs_CHOW", "p_NORMAL_DIET_vs_CHOW",
      "NORMAL_DIET_vs_CHOW_FDR_BH",
      "tau2_GSE", "tau2_contrast_within_GSE", "tau2_total"
    ),
    names(modres)
  )

  r <- merge(
    r,
    modres[, ..keep],
    by = "Geneid",
    all.x = TRUE,
    sort = FALSE,
    suffixes = c("", "_moderator")
  )

  f <- if (any(fail)) {
    rbindlist(lapply(ans[fail], function(z) {
      data.frame(
        Geneid = z$Geneid,
        Moderator = label,
        Error = z$reason,
        stringsAsFactors = FALSE
      )
    }), fill = TRUE)
  } else {
    data.table(Geneid = character(), Moderator = character(), Error = character())
  }

  list(result = r, failures = f)
}

dur_r2 <- run_matched_null(dur, "duration", "duration")
fat_r2 <- run_matched_null(fat, "fat", "fat_percent")
ctrl_r2 <- run_matched_null(ctrl, "control", "control_diet")

fwrite(
  dur_r2$result, OUT_DUR_R2,
  sep = "\t", quote = FALSE, na = "NA", compress = "gzip"
)
fwrite(
  fat_r2$result, OUT_FAT_R2,
  sep = "\t", quote = FALSE, na = "NA", compress = "gzip"
)
fwrite(
  ctrl_r2$result, OUT_CTRL_R2,
  sep = "\t", quote = FALSE, na = "NA", compress = "gzip"
)

null_fail <- rbindlist(
  list(dur_r2$failures, fat_r2$failures, ctrl_r2$failures),
  fill = TRUE
)
fwrite(null_fail, OUT_NULL_FAIL, sep = "\t", quote = FALSE, na = "NA")

# ------------------------------------------------------------------
# Prespecified multivariable sensitivity: duration + fat%
# ------------------------------------------------------------------
fit_multi_gene <- function(d) {
  x <- d[is.finite(hfd_duration_weeks) & is.finite(hfd_fat_kcal_percent)]

  if (uniqueN(x$GSE) < 20L) {
    return(list(
      status = "INELIGIBLE", Geneid = d$Geneid[1],
      reason = paste0("known_joint_GSE_lt20:", uniqueN(x$GSE))
    ))
  }
  if (uniqueN(x$hfd_duration_weeks) < 2L ||
      uniqueN(x$hfd_fat_kcal_percent) < 2L) {
    return(list(
      status = "INELIGIBLE", Geneid = d$Geneid[1],
      reason = "insufficient_joint_moderator_variation"
    ))
  }

  rr <- fit_multi_rescue(x)
  fit <- rr$fit

  if (inherits(fit, "error")) {
    return(list(
      status = "FAIL", Geneid = d$Geneid[1],
      reason = conditionMessage(fit)
    ))
  }

  cn <- names(coef(fit))
  id <- which(cn == "duration_per4")
  iff <- which(cn == "fat_per10")
  if (length(id) != 1L || length(iff) != 1L) {
    return(list(
      status = "FAIL", Geneid = d$Geneid[1],
      reason = paste0("unexpected_coefficients:", paste(cn, collapse = ","))
    ))
  }

  out <- data.frame(
    Geneid = d$Geneid[1],
    N_contrasts = nrow(x),
    N_independent_GSE = uniqueN(x$GSE),

    beta_duration_per4 = as.numeric(coef(fit)[id]),
    SE_duration_per4 = as.numeric(fit$se[id]),
    z_duration_per4 = as.numeric(fit$zval[id]),
    p_duration_per4 = as.numeric(fit$pval[id]),

    beta_fat_per10 = as.numeric(coef(fit)[iff]),
    SE_fat_per10 = as.numeric(fit$se[iff]),
    z_fat_per10 = as.numeric(fit$zval[iff]),
    p_fat_per10 = as.numeric(fit$pval[iff]),

    QM = as.numeric(fit$QM),
    QM_p_value = as.numeric(fit$QMp),

    tau2_GSE = as.numeric(fit$sigma2[1]),
    tau2_contrast_within_GSE = as.numeric(fit$sigma2[2]),
    tau2_total = as.numeric(sum(fit$sigma2)),
    optimizer = rr$optimizer,
    stringsAsFactors = FALSE
  )

  list(status = "PASS", result = out, Geneid = d$Geneid[1], reason = "")
}

# Preflight multivariable model.
pf_multi <- NULL
for (d in gene_split) {
  z <- fit_multi_gene(d)
  if (z$status == "PASS") {
    pf_multi <- z
    break
  }
}
if (is.null(pf_multi))
  stop("PRECHECK_FAIL: no duration+fat multivariable model could be fit")
cat("PRECHECK duration+fat PASS:", pf_multi$Geneid, "\n")

multi_ans <- parallel::mclapply(
  gene_split,
  fit_multi_gene,
  mc.cores = NCORES,
  mc.preschedule = TRUE
)

multi_pass <- vapply(multi_ans, function(z) z$status == "PASS", logical(1))
multi_fail <- vapply(multi_ans, function(z) z$status == "FAIL", logical(1))
multi_inel <- vapply(multi_ans, function(z) z$status == "INELIGIBLE", logical(1))

multi <- rbindlist(
  lapply(multi_ans[multi_pass], `[[`, "result"),
  fill = TRUE
)

multi[, duration_FDR_BH := p.adjust(p_duration_per4, method = "BH")]
multi[, fat_FDR_BH := p.adjust(p_fat_per10, method = "BH")]
multi[, omnibus_FDR_BH := p.adjust(QM_p_value, method = "BH")]

multi[, abs_duration_sort := abs(beta_duration_per4)]
setorder(multi, duration_FDR_BH, -abs_duration_sort, Geneid)
multi[, abs_duration_sort := NULL]

fwrite(
  multi, OUT_MULTI,
  sep = "\t", quote = FALSE, na = "NA", compress = "gzip"
)

multi_fail_df <- if (any(multi_fail)) {
  rbindlist(lapply(multi_ans[multi_fail], function(z) {
    data.frame(
      Geneid = z$Geneid,
      Status = "FAIL",
      Reason = z$reason,
      stringsAsFactors = FALSE
    )
  }), fill = TRUE)
} else {
  data.table(Geneid = character(), Status = character(), Reason = character())
}
fwrite(
  multi_fail_df, OUT_MULTI_FAIL,
  sep = "\t", quote = FALSE, na = "NA"
)

# ------------------------------------------------------------------
# Moderator correlation in jointly observed design.
# ------------------------------------------------------------------
design <- unique(dt[, .(
  Contrast_id, GSE, hfd_duration_weeks, hfd_fat_kcal_percent
)])

joint <- design[
  is.finite(hfd_duration_weeks) & is.finite(hfd_fat_kcal_percent)
]

gse_joint <- joint[, .(
  duration_weeks = median(hfd_duration_weeks),
  fat_kcal_percent = median(hfd_fat_kcal_percent)
), by = GSE]

cor_contrast_pearson <- cor(
  joint$hfd_duration_weeks,
  joint$hfd_fat_kcal_percent,
  method = "pearson"
)
cor_contrast_spearman <- cor(
  joint$hfd_duration_weeks,
  joint$hfd_fat_kcal_percent,
  method = "spearman"
)
cor_gse_pearson <- cor(
  gse_joint$duration_weeks,
  gse_joint$fat_kcal_percent,
  method = "pearson"
)
cor_gse_spearman <- cor(
  gse_joint$duration_weeks,
  gse_joint$fat_kcal_percent,
  method = "spearman"
)

summ_r2 <- function(r, sig_col = NULL) {
  vals <- r$pseudo_R2_total_percent
  vals <- vals[is.finite(vals)]

  if (is.null(sig_col)) {
    sigvals <- numeric()
  } else {
    sigvals <- r[get(sig_col) < 0.05, pseudo_R2_total_percent]
    sigvals <- sigvals[is.finite(sigvals)]
  }

  list(
    median_all = if (length(vals)) median(vals) else NA_real_,
    positive_all = if (length(vals)) mean(vals > 0) else NA_real_,
    median_sig = if (length(sigvals)) median(sigvals) else NA_real_,
    positive_sig = if (length(sigvals)) mean(sigvals > 0) else NA_real_
  )
}

sdur <- summ_r2(dur_r2$result, "FDR_BH")
sfat <- summ_r2(fat_r2$result, "FDR_BH")
sctrl <- summ_r2(ctrl_r2$result, "omnibus_FDR_BH")

status <- if (
  nrow(null_fail) == 0L &&
  sum(multi_fail) == 0L &&
  nrow(multi) > 0L
) "PASS" else "PASS_WITH_MODEL_FAILURES"

summary <- c(
  "===== HFD MODERATOR HETEROGENEITY / MULTIVARIABLE SUMMARY v1 =====",
  "",
  paste0("STATUS=", status),
  paste0("R_version=", R.version.string),
  paste0("metafor_version=", as.character(packageVersion("metafor"))),
  paste0("Parallel_workers=", NCORES),
  "",
  "MATCHED_SUBSET_HETEROGENEITY_REDUCTION=descriptive_pseudo_R2_from_total_sigma2",
  "PSEUDO_R2_FORMULA=100*(tau2_null_total-tau2_moderator_total)/tau2_null_total",
  "NEGATIVE_PSEUDO_R2_VALUES_RETAINED=YES",
  "RMA_MV_BUILTIN_R2_USED=NO",
  "",
  paste0("Duration_matched_null_models=", nrow(dur_r2$result)),
  paste0("Duration_matched_null_failures=", nrow(dur_r2$failures)),
  paste0("Duration_median_pseudo_R2_all_percent=", format(sdur$median_all, digits = 6)),
  paste0("Duration_fraction_positive_R2_all=", format(sdur$positive_all, digits = 6)),
  paste0("Duration_median_pseudo_R2_FDR05_percent=", format(sdur$median_sig, digits = 6)),
  paste0("Duration_fraction_positive_R2_FDR05=", format(sdur$positive_sig, digits = 6)),
  "",
  paste0("Fat_matched_null_models=", nrow(fat_r2$result)),
  paste0("Fat_matched_null_failures=", nrow(fat_r2$failures)),
  paste0("Fat_median_pseudo_R2_all_percent=", format(sfat$median_all, digits = 6)),
  paste0("Fat_fraction_positive_R2_all=", format(sfat$positive_all, digits = 6)),
  paste0("Fat_median_pseudo_R2_FDR05_percent=", format(sfat$median_sig, digits = 6)),
  paste0("Fat_fraction_positive_R2_FDR05=", format(sfat$positive_sig, digits = 6)),
  "",
  paste0("Control_matched_null_models=", nrow(ctrl_r2$result)),
  paste0("Control_matched_null_failures=", nrow(ctrl_r2$failures)),
  paste0("Control_median_pseudo_R2_all_percent=", format(sctrl$median_all, digits = 6)),
  paste0("Control_fraction_positive_R2_all=", format(sctrl$positive_all, digits = 6)),
  paste0("Control_median_pseudo_R2_FDR05_percent=", format(sctrl$median_sig, digits = 6)),
  paste0("Control_fraction_positive_R2_FDR05=", format(sctrl$positive_sig, digits = 6)),
  "",
  "MULTIVARIABLE_SENSITIVITY=duration_plus_fat_percent",
  "MULTIVARIABLE_MIN_JOINT_KNOWN_GSE=20",
  paste0("Multivariable_genes_fit=", nrow(multi)),
  paste0("Multivariable_genes_ineligible=", sum(multi_inel)),
  paste0("Multivariable_model_failures=", sum(multi_fail)),
  paste0("Multivariable_duration_FDR_lt_0.05=",
         sum(multi$duration_FDR_BH < 0.05, na.rm = TRUE)),
  paste0("Multivariable_fat_FDR_lt_0.05=",
         sum(multi$fat_FDR_BH < 0.05, na.rm = TRUE)),
  paste0("Multivariable_omnibus_FDR_lt_0.05=",
         sum(multi$omnibus_FDR_BH < 0.05, na.rm = TRUE)),
  "",
  paste0("Joint_design_contrasts=", nrow(joint)),
  paste0("Joint_design_GSE=", nrow(gse_joint)),
  paste0("Duration_fat_contrast_Pearson_r=", format(cor_contrast_pearson, digits = 6)),
  paste0("Duration_fat_contrast_Spearman_rho=", format(cor_contrast_spearman, digits = 6)),
  paste0("Duration_fat_GSEbalanced_Pearson_r=", format(cor_gse_pearson, digits = 6)),
  paste0("Duration_fat_GSEbalanced_Spearman_rho=", format(cor_gse_spearman, digits = 6)),
  "",
  paste0("Duration_R2_result=", sub(paste0("^", ROOT, "/"), "", OUT_DUR_R2)),
  paste0("Fat_R2_result=", sub(paste0("^", ROOT, "/"), "", OUT_FAT_R2)),
  paste0("Control_R2_result=", sub(paste0("^", ROOT, "/"), "", OUT_CTRL_R2)),
  paste0("Duration_plus_fat_result=", sub(paste0("^", ROOT, "/"), "", OUT_MULTI)),
  paste0("Matched_null_failures=", sub(paste0("^", ROOT, "/"), "", OUT_NULL_FAIL)),
  paste0("Multivariable_failures=", sub(paste0("^", ROOT, "/"), "", OUT_MULTI_FAIL)),
  "",
  "INTERPRETATION_CAUSAL=NO",
  "MODERATORS_ARE_PRIMARILY_BETWEEN_GSE=YES",
  "MODERATOR_PATHWAY_ANALYSIS_PERFORMED=NO"
)

writeLines(summary, OUT_SUMMARY)
cat("\n", paste(summary, collapse = "\n"), "\n", sep = "")

if (status != "PASS") quit(status = 2L)

