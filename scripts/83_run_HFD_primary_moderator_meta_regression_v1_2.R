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

OUTDIR <- file.path(ROOT, "05_meta", "moderators_v1_2")
METADIR <- file.path(ROOT, "00_metadata/analysis_outputs/moderators_v1")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(METADIR, recursive=TRUE, showWarnings=FALSE)

OUT_DUR <- file.path(OUTDIR, "HFD_duration_meta_regression_v1_2.tsv.gz")
OUT_FAT <- file.path(OUTDIR, "HFD_fat_percent_meta_regression_v1_2.tsv.gz")
OUT_CTRL <- file.path(OUTDIR, "HFD_control_diet_meta_regression_v1_2.tsv.gz")
OUT_FAIL <- file.path(METADIR, "HFD_primary_moderator_model_failures_v1_2.tsv")
OUT_SUM <- file.path(METADIR, "HFD_primary_moderator_meta_regression_summary_v1_2.txt")

if (!file.exists(LONG)) stop("Missing input: ", LONG)

need <- c(
  "Contrast_id","GSE","Geneid","log2FoldChange","sampling_variance",
  "hfd_duration_weeks","hfd_fat_kcal_percent","control_diet_type"
)

dt <- fread(LONG, select=need, na.strings=c("NA",""))
missing <- setdiff(need, names(dt))
if (length(missing)) stop("Long table missing columns: ", paste(missing, collapse=", "))

dt[, log2FoldChange := as.numeric(log2FoldChange)]
dt[, sampling_variance := as.numeric(sampling_variance)]
dt[, hfd_duration_weeks := suppressWarnings(as.numeric(hfd_duration_weeks))]
dt[, hfd_fat_kcal_percent := suppressWarnings(as.numeric(hfd_fat_kcal_percent))]
dt[, control_diet_type := trimws(as.character(control_diet_type))]
dt[control_diet_type == "", control_diet_type := NA_character_]

if (any(!is.finite(dt$log2FoldChange))) stop("Non-finite effect")
if (any(!is.finite(dt$sampling_variance)) || any(dt$sampling_variance <= 0))
  stop("Invalid sampling variance")

# Prespecified scaling/centering.
# Coefficients are per +4 weeks and +10 kcal%-points.
dt[, duration_per4 := (hfd_duration_weeks - 12) / 4]
dt[, fat_per10 := (hfd_fat_kcal_percent - 60) / 10]

# Primary control-diet levels and reference.
allowed_ctrl <- c("CHOW","LFD","NORMAL_DIET")
bad_ctrl <- setdiff(unique(na.omit(dt$control_diet_type)), allowed_ctrl)
if (length(bad_ctrl))
  stop("Unexpected control_diet_type values: ", paste(bad_ctrl, collapse=", "))

cores_env <- Sys.getenv("HFD_META_CORES", unset="")
if (nzchar(cores_env)) {
  NCORES <- as.integer(cores_env)
} else {
  NCORES <- max(1L, min(8L, parallel::detectCores(logical=FALSE)-1L))
}
if (!is.finite(NCORES) || NCORES < 1L) NCORES <- 1L

cat("===== HFD PRIMARY MODERATOR META-REGRESSION v1.2 =====\n")
cat("metafor:", as.character(packageVersion("metafor")), "\n")
cat("Workers:", NCORES, "\n\n")

fit_mv_rescue <- function(d, mods_formula) {
  dofit <- function(opt = NULL) {
    if (is.null(opt)) {
      rma.mv(
        yi = d$log2FoldChange,
        V = d$sampling_variance,
        mods = mods_formula,
        random = list(~ 1 | GSE, ~ 1 | Contrast_id),
        data = d,
        method = "REML",
        test = "z",
        sparse = TRUE
      )
    } else {
      rma.mv(
        yi = d$log2FoldChange,
        V = d$sampling_variance,
        mods = mods_formula,
        random = list(~ 1 | GSE, ~ 1 | Contrast_id),
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

# ------------------------------
# Continuous moderator fitter
# ------------------------------
fit_cont_gene <- function(d, moderator, scaled_name, label) {
  x <- d[is.finite(get(moderator))]
  k_gse <- uniqueN(x$GSE)

  if (k_gse < 20L) {
    return(list(
      status="INELIGIBLE", Geneid=d$Geneid[1],
      reason=paste0("known_GSE_lt20:", k_gse)
    ))
  }
  if (uniqueN(x[[moderator]]) < 2L) {
    return(list(
      status="INELIGIBLE", Geneid=d$Geneid[1],
      reason="no_moderator_variation"
    ))
  }

  rr <- fit_mv_rescue(x, as.formula(paste0("~ ", scaled_name)))
  fit <- rr$fit
  if (inherits(fit,"error")) {
    return(list(
      status="FAIL", Geneid=d$Geneid[1],
      reason=conditionMessage(fit), optimizer=rr$optimizer
    ))
  }

  # coefficient 2 is moderator slope.
  b <- as.numeric(coef(fit)[2])
  se <- as.numeric(fit$se[2])
  z <- as.numeric(fit$zval[2])
  p <- as.numeric(fit$pval[2])

  out <- data.frame(
    Geneid=d$Geneid[1],
    Moderator=label,
    N_contrasts=nrow(x),
    N_independent_GSE=k_gse,
    Moderator_min=min(x[[moderator]], na.rm=TRUE),
    Moderator_median=median(x[[moderator]], na.rm=TRUE),
    Moderator_max=max(x[[moderator]], na.rm=TRUE),
    beta=b,
    SE=se,
    CI95_lower=b-qnorm(.975)*se,
    CI95_upper=b+qnorm(.975)*se,
    z_value=z,
    p_value=p,
    QM=as.numeric(fit$QM),
    QM_df=as.integer(fit$m),
    QM_p_value=as.numeric(fit$QMp),
    tau2_GSE=as.numeric(fit$sigma2[1]),
    tau2_contrast_within_GSE=as.numeric(fit$sigma2[2]),
    tau2_total=as.numeric(sum(fit$sigma2)),
    optimizer=rr$optimizer,
    stringsAsFactors=FALSE
  )
  list(status="PASS", result=out, Geneid=d$Geneid[1], reason="")
}

# ------------------------------
# Control-diet categorical fitter
# Primary eligibility:
# >=20 GSE overall, and each of CHOW/LFD/NORMAL_DIET >=5 independent GSE.
# This freezes the same 3-level model across genes.
# ------------------------------
fit_ctrl_gene <- function(d) {
  x <- d[!is.na(control_diet_type) & control_diet_type %in% allowed_ctrl]
  k_gse <- uniqueN(x$GSE)

  bycat <- x[, .(N_GSE=uniqueN(GSE), N_contrasts=.N), by=control_diet_type]
  gse_counts <- setNames(rep(0L, length(allowed_ctrl)), allowed_ctrl)
  con_counts <- setNames(rep(0L, length(allowed_ctrl)), allowed_ctrl)
  if (nrow(bycat)) {
    gse_counts[bycat$control_diet_type] <- bycat$N_GSE
    con_counts[bycat$control_diet_type] <- bycat$N_contrasts
  }

  if (k_gse < 20L) {
    return(list(
      status="INELIGIBLE", Geneid=d$Geneid[1],
      reason=paste0("known_GSE_lt20:", k_gse)
    ))
  }
  if (any(gse_counts < 5L)) {
    return(list(
      status="INELIGIBLE", Geneid=d$Geneid[1],
      reason=paste0(
        "category_GSE_lt5:",
        paste(names(gse_counts), gse_counts, sep="=", collapse=",")
      )
    ))
  }

  x[, control_diet_type := factor(
    control_diet_type,
    levels=c("CHOW","LFD","NORMAL_DIET")
  )]

  rr <- fit_mv_rescue(x, ~ control_diet_type)
  fit <- rr$fit
  if (inherits(fit,"error")) {
    return(list(
      status="FAIL", Geneid=d$Geneid[1],
      reason=conditionMessage(fit), optimizer=rr$optimizer
    ))
  }

  cn <- names(coef(fit))
  i_lfd <- which(cn == "control_diet_typeLFD")
  i_norm <- which(cn == "control_diet_typeNORMAL_DIET")
  if (length(i_lfd)!=1L || length(i_norm)!=1L)
    stop("Unexpected control-diet coefficient names for ", d$Geneid[1],
         ": ", paste(cn, collapse=", "))

  out <- data.frame(
    Geneid=d$Geneid[1],
    Moderator="control_diet_type",
    N_contrasts=nrow(x),
    N_independent_GSE=k_gse,
    N_GSE_CHOW=gse_counts["CHOW"],
    N_GSE_LFD=gse_counts["LFD"],
    N_GSE_NORMAL_DIET=gse_counts["NORMAL_DIET"],
    N_contrasts_CHOW=con_counts["CHOW"],
    N_contrasts_LFD=con_counts["LFD"],
    N_contrasts_NORMAL_DIET=con_counts["NORMAL_DIET"],
    beta_LFD_vs_CHOW=as.numeric(coef(fit)[i_lfd]),
    SE_LFD_vs_CHOW=as.numeric(fit$se[i_lfd]),
    z_LFD_vs_CHOW=as.numeric(fit$zval[i_lfd]),
    p_LFD_vs_CHOW=as.numeric(fit$pval[i_lfd]),
    beta_NORMAL_DIET_vs_CHOW=as.numeric(coef(fit)[i_norm]),
    SE_NORMAL_DIET_vs_CHOW=as.numeric(fit$se[i_norm]),
    z_NORMAL_DIET_vs_CHOW=as.numeric(fit$zval[i_norm]),
    p_NORMAL_DIET_vs_CHOW=as.numeric(fit$pval[i_norm]),
    QM=as.numeric(fit$QM),
    QM_df=as.integer(fit$m),
    QM_p_value=as.numeric(fit$QMp),
    tau2_GSE=as.numeric(fit$sigma2[1]),
    tau2_contrast_within_GSE=as.numeric(fit$sigma2[2]),
    tau2_total=as.numeric(sum(fit$sigma2)),
    optimizer=rr$optimizer,
    stringsAsFactors=FALSE
  )
  list(status="PASS", result=out, Geneid=d$Geneid[1], reason="")
}


# ------------------------------------------------------------------
# PRE-FLIGHT: fit one genuinely eligible gene for each moderator
# sequentially before launching thousands of parallel models.
# ------------------------------------------------------------------
gene_split <- split(dt, by="Geneid", keep.by=TRUE, drop=TRUE)

find_cont_preflight <- function(moderator) {
  for (d in gene_split) {
    x <- d[is.finite(get(moderator))]
    if (uniqueN(x$GSE) >= 20L && uniqueN(x[[moderator]]) >= 2L)
      return(d)
  }
  NULL
}

find_ctrl_preflight <- function() {
  for (d in gene_split) {
    x <- d[!is.na(control_diet_type) & control_diet_type %in% allowed_ctrl]
    if (uniqueN(x$GSE) < 20L) next
    bycat <- x[, .(N_GSE=uniqueN(GSE)), by=control_diet_type]
    counts <- setNames(rep(0L, length(allowed_ctrl)), allowed_ctrl)
    if (nrow(bycat)) counts[bycat$control_diet_type] <- bycat$N_GSE
    if (all(counts >= 5L)) return(d)
  }
  NULL
}

pf_dur_d <- find_cont_preflight("hfd_duration_weeks")
if (is.null(pf_dur_d)) stop("PRECHECK_FAIL: no duration-eligible gene found")
pf_dur <- fit_cont_gene(
  pf_dur_d, "hfd_duration_weeks", "duration_per4", "duration_per_4_weeks"
)
if (pf_dur$status != "PASS")
  stop("PRECHECK_FAIL duration: ", pf_dur$reason)

pf_fat_d <- find_cont_preflight("hfd_fat_kcal_percent")
if (is.null(pf_fat_d)) stop("PRECHECK_FAIL: no fat-eligible gene found")
pf_fat <- fit_cont_gene(
  pf_fat_d, "hfd_fat_kcal_percent", "fat_per10", "fat_percent_per_10_points"
)
if (pf_fat$status != "PASS")
  stop("PRECHECK_FAIL fat: ", pf_fat$reason)

pf_ctrl_d <- find_ctrl_preflight()
if (is.null(pf_ctrl_d)) stop("PRECHECK_FAIL: no control-diet-eligible gene found")
pf_ctrl <- fit_ctrl_gene(pf_ctrl_d)
if (pf_ctrl$status != "PASS")
  stop("PRECHECK_FAIL control diet: ", pf_ctrl$reason)

cat("PRECHECK duration PASS:", pf_dur_d$Geneid[1], "\n")
cat("PRECHECK fat PASS:", pf_fat_d$Geneid[1], "\n")
cat("PRECHECK control diet PASS:", pf_ctrl_d$Geneid[1], "\n\n")

run_cont <- function(moderator, scaled_name, label) {
  cat("Running:", label, "\n")
  ans <- parallel::mclapply(
    gene_split,
    fit_cont_gene,
    moderator=moderator,
    scaled_name=scaled_name,
    label=label,
    mc.cores=NCORES,
    mc.preschedule=TRUE
  )

  pass <- vapply(ans, function(z) z$status=="PASS", logical(1))
  fail <- vapply(ans, function(z) z$status=="FAIL", logical(1))
  inel <- vapply(ans, function(z) z$status=="INELIGIBLE", logical(1))

  res <- if (any(pass)) rbindlist(lapply(ans[pass], `[[`, "result"), fill=TRUE) else data.table()
  if (nrow(res)) {
    res[, FDR_BH := p.adjust(p_value, method="BH")]
    res[, significant_FDR05 := ifelse(FDR_BH < .05, "YES","NO")]
  }

  failures <- rbindlist(lapply(ans[fail], function(z)
    data.frame(Geneid=z$Geneid, Moderator=label, Status="FAIL",
               Reason=z$reason, stringsAsFactors=FALSE)), fill=TRUE)
  ineligible <- rbindlist(lapply(ans[inel], function(z)
    data.frame(Geneid=z$Geneid, Moderator=label, Status="INELIGIBLE",
               Reason=z$reason, stringsAsFactors=FALSE)), fill=TRUE)

  list(res=res, failures=failures, ineligible=ineligible)
}

# data.table cannot sort by expression inside setorder; helper fixes that.
safe_sort_cont <- function(res) {
  if (!nrow(res)) return(res)
  res[, abs_beta_sort := abs(beta)]
  setorder(res, FDR_BH, -abs_beta_sort, Geneid)
  res[, abs_beta_sort := NULL]
  res
}

dur <- run_cont("hfd_duration_weeks", "duration_per4", "duration_per_4_weeks")
dur$res <- safe_sort_cont(dur$res)
fwrite(dur$res, OUT_DUR, sep="\t", quote=FALSE, na="NA", compress="gzip")

fat <- run_cont("hfd_fat_kcal_percent", "fat_per10", "fat_percent_per_10_points")
fat$res <- safe_sort_cont(fat$res)
fwrite(fat$res, OUT_FAT, sep="\t", quote=FALSE, na="NA", compress="gzip")

cat("Running: control_diet_type\n")
ctrl_ans <- parallel::mclapply(
  gene_split,
  fit_ctrl_gene,
  mc.cores=NCORES,
  mc.preschedule=TRUE
)
ctrl_pass <- vapply(ctrl_ans, function(z) z$status=="PASS", logical(1))
ctrl_fail <- vapply(ctrl_ans, function(z) z$status=="FAIL", logical(1))
ctrl_inel <- vapply(ctrl_ans, function(z) z$status=="INELIGIBLE", logical(1))

ctrl <- if (any(ctrl_pass)) rbindlist(lapply(ctrl_ans[ctrl_pass], `[[`, "result"), fill=TRUE) else data.table()
if (nrow(ctrl)) {
  ctrl[, omnibus_FDR_BH := p.adjust(QM_p_value, method="BH")]
  ctrl[, omnibus_significant_FDR05 := ifelse(omnibus_FDR_BH < .05, "YES","NO")]
  ctrl[, LFD_vs_CHOW_FDR_BH := p.adjust(p_LFD_vs_CHOW, method="BH")]
  ctrl[, NORMAL_DIET_vs_CHOW_FDR_BH := p.adjust(p_NORMAL_DIET_vs_CHOW, method="BH")]
  setorder(ctrl, omnibus_FDR_BH, Geneid)
}
fwrite(ctrl, OUT_CTRL, sep="\t", quote=FALSE, na="NA", compress="gzip")

ctrl_fail_df <- rbindlist(lapply(ctrl_ans[ctrl_fail], function(z)
  data.frame(Geneid=z$Geneid, Moderator="control_diet_type", Status="FAIL",
             Reason=z$reason, stringsAsFactors=FALSE)), fill=TRUE)
ctrl_inel_df <- rbindlist(lapply(ctrl_ans[ctrl_inel], function(z)
  data.frame(Geneid=z$Geneid, Moderator="control_diet_type", Status="INELIGIBLE",
             Reason=z$reason, stringsAsFactors=FALSE)), fill=TRUE)

all_problem <- rbindlist(
  list(
    dur$failures, dur$ineligible,
    fat$failures, fat$ineligible,
    ctrl_fail_df, ctrl_inel_df
  ),
  fill=TRUE
)
fwrite(all_problem, OUT_FAIL, sep="\t", quote=FALSE, na="NA")

count_sig <- function(x, col) {
  if (!nrow(x)) return(0L)
  sum(x[[col]] < .05, na.rm=TRUE)
}
count_pos <- function(x) {
  if (!nrow(x)) return(0L)
  sum(x$FDR_BH < .05 & x$beta > 0, na.rm=TRUE)
}
count_neg <- function(x) {
  if (!nrow(x)) return(0L)
  sum(x$FDR_BH < .05 & x$beta < 0, na.rm=TRUE)
}

summary <- c(
  "===== HFD PRIMARY MODERATOR META-REGRESSION SUMMARY v1.2 =====",
  "",
  paste0(
    "STATUS=",
    ifelse(
      nrow(dur$res) == 0L || nrow(fat$res) == 0L || nrow(ctrl) == 0L,
      "FAIL",
      ifelse(
        nrow(dur$failures) + nrow(fat$failures) + sum(ctrl_fail) > 0L,
        "PASS_WITH_MODEL_FAILURES",
        "PASS"
      )
    )
  ),
  paste0("R_version=", R.version.string),
  paste0("metafor_version=", as.character(packageVersion("metafor"))),
  paste0("Parallel_workers=", NCORES),
  "",
  "FROZEN_META_INPUT=00_metadata/analysis_inputs/meta_analysis_v1/HFD_primary105_meta_effects_primaryQC_v1.tsv.gz",
  "MODEL=rma.mv_REML_random_GSE_plus_Contrast",
  "PRIMARY_MODERATORS=duration_continuous;fat_kcal_percent_continuous;control_diet_type_categorical",
  "DURATION_EFFECT_UNIT=per_4_weeks",
  "DURATION_CENTER=12_weeks",
  "FAT_EFFECT_UNIT=per_10_kcal_percent_points",
  "FAT_CENTER=60_percent",
  "CONTROL_DIET_REFERENCE=CHOW",
  "CONTROL_DIET_LEVELS=CHOW;LFD;NORMAL_DIET",
  "",
  "GENE_LEVEL_ELIGIBILITY_CONTINUOUS=min_20_independent_GSE_with_known_moderator",
  "GENE_LEVEL_ELIGIBILITY_CONTROL_DIET=min_20_GSE_total_and_min_5_GSE_each_CHOW_LFD_NORMAL_DIET",
  "",
  paste0("Duration_genes_fit=", nrow(dur$res)),
  paste0("Duration_genes_ineligible=", nrow(dur$ineligible)),
  paste0("Duration_model_failures=", nrow(dur$failures)),
  paste0("Duration_FDR_lt_0.05=", count_sig(dur$res, "FDR_BH")),
  paste0("Duration_FDR05_positive_beta=", count_pos(dur$res)),
  paste0("Duration_FDR05_negative_beta=", count_neg(dur$res)),
  "",
  paste0("Fat_genes_fit=", nrow(fat$res)),
  paste0("Fat_genes_ineligible=", nrow(fat$ineligible)),
  paste0("Fat_model_failures=", nrow(fat$failures)),
  paste0("Fat_FDR_lt_0.05=", count_sig(fat$res, "FDR_BH")),
  paste0("Fat_FDR05_positive_beta=", count_pos(fat$res)),
  paste0("Fat_FDR05_negative_beta=", count_neg(fat$res)),
  "",
  paste0("Control_diet_genes_fit=", nrow(ctrl)),
  paste0("Control_diet_genes_ineligible=", sum(ctrl_inel)),
  paste0("Control_diet_model_failures=", sum(ctrl_fail)),
  paste0("Control_diet_omnibus_FDR_lt_0.05=",
         ifelse(nrow(ctrl), sum(ctrl$omnibus_FDR_BH < .05, na.rm=TRUE), 0L)),
  paste0("Control_diet_LFD_vs_CHOW_FDR_lt_0.05=",
         ifelse(nrow(ctrl), sum(ctrl$LFD_vs_CHOW_FDR_BH < .05, na.rm=TRUE), 0L)),
  paste0("Control_diet_NORMAL_DIET_vs_CHOW_FDR_lt_0.05=",
         ifelse(nrow(ctrl), sum(ctrl$NORMAL_DIET_vs_CHOW_FDR_BH < .05, na.rm=TRUE), 0L)),
  "",
  paste0("Duration_result=", sub(paste0("^",ROOT,"/"),"",OUT_DUR)),
  paste0("Fat_result=", sub(paste0("^",ROOT,"/"),"",OUT_FAT)),
  paste0("Control_diet_result=", sub(paste0("^",ROOT,"/"),"",OUT_CTRL)),
  paste0("Eligibility_and_failure_audit=", sub(paste0("^",ROOT,"/"),"",OUT_FAIL)),
  "",
  "SEX_STRAIN_ANALYSIS_PERFORMED=NO",
  "MULTIVARIABLE_DURATION_PLUS_FAT_PERFORMED=NO",
  "MODERATOR_PATHWAY_ANALYSIS_PERFORMED=NO"
)

writeLines(summary, OUT_SUM)
cat("\n", paste(summary, collapse="\n"), "\n", sep="")

if (nrow(dur$res) == 0L || nrow(fat$res) == 0L || nrow(ctrl) == 0L) {
  quit(status = 1L)
}
if (nrow(dur$failures) + nrow(fat$failures) + sum(ctrl_fail) > 0L) {
  quit(status = 2L)
}

