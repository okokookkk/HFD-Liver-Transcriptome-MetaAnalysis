#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  req <- c("data.table", "metafor")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) stop("Missing required R packages: ", paste(miss, collapse=", "))
})
library(data.table)
library(metafor)

ROOT <- normalizePath(Sys.getenv("HFD_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)

D126 <- file.path(
  ROOT, "07_human_validation/cohort_DE_v1_2",
  "GSE126848_MASLD_vs_control_voom_DE_v1_2.tsv.gz"
)
D130 <- file.path(
  ROOT, "07_human_validation/cohort_DE_v1_2",
  "GSE130970_MASLD_vs_control6_voom_DE_v1_2.tsv.gz"
)
D130S <- file.path(
  ROOT, "07_human_validation/cohort_DE_v1_2",
  "GSE130970_MASLD_vs_strict_control4_voom_DE_v1_2.tsv.gz"
)
D135 <- file.path(
  ROOT, "07_human_validation/cohort_DE_v1_2",
  "GSE135251_MASLD_vs_control_voom_DE_v1_2.tsv.gz"
)

OUTDIR <- file.path(ROOT, "07_human_validation", "gene_meta_v1")
METADIR <- file.path(
  ROOT, "00_metadata/human_validation/analysis_outputs/gene_meta_v1"
)
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(METADIR, recursive=TRUE, showWarnings=FALSE)

OUT_PRIMARY <- file.path(
  OUTDIR, "human_MASLD_three_cohort_random_effects_meta_v1.tsv.gz"
)
OUT_STRICT <- file.path(
  OUTDIR, "human_MASLD_three_cohort_strict4_sensitivity_meta_v1.tsv.gz"
)
OUT_PAIR <- file.path(
  OUTDIR, "human_MASLD_cohort_effect_concordance_v1.tsv"
)
OUT_COMPARE <- file.path(
  OUTDIR, "human_MASLD_primary_vs_strict4_meta_comparison_v1.tsv.gz"
)
OUT_FAIL <- file.path(
  METADIR, "human_MASLD_gene_meta_model_failures_v1.tsv"
)
OUT_SUM <- file.path(
  METADIR, "human_MASLD_gene_meta_summary_v1.txt"
)

for (q in c(D126, D130, D130S, D135)) {
  if (!file.exists(q)) stop("Missing required input: ", q)
}

read_de <- function(path, cohort) {
  x <- fread(path)
  need <- c("ENTREZID","human_log2FC","human_SE","human_sampling_variance")
  miss <- setdiff(need, names(x))
  if (length(miss)) stop(cohort, " missing columns: ", paste(miss, collapse=","))
  x <- x[, .(
    ENTREZID=as.character(ENTREZID),
    yi=as.numeric(human_log2FC),
    sei=as.numeric(human_SE),
    vi=as.numeric(human_sampling_variance)
  )]
  if (anyDuplicated(x$ENTREZID)) stop(cohort, ": duplicate ENTREZID")
  if (any(!is.finite(x$yi)) || any(!is.finite(x$sei)) ||
      any(!is.finite(x$vi)) || any(x$sei <= 0) || any(x$vi <= 0))
    stop(cohort, ": invalid effect/SE/variance")
  x[, Cohort := cohort]
  x
}

a <- read_de(D126, "GSE126848")
b <- read_de(D130, "GSE130970")
bs <- read_de(D130S, "GSE130970_STRICT4")
c <- read_de(D135, "GSE135251")

cat("===== HUMAN MASLD THREE-COHORT GENE META-ANALYSIS v1 =====\n")
cat("metafor:", as.character(packageVersion("metafor")), "\n\n")

# ------------------------------------------------------------------
# Pairwise cohort effect concordance before pooling.
# Uses all genes jointly tested in each pair.
# ------------------------------------------------------------------
pair_concordance <- function(x, y, name_x, name_y) {
  z <- merge(
    x[,.(ENTREZID, x=yi)],
    y[,.(ENTREZID, y=yi)],
    by="ENTREZID"
  )
  data.table(
    Cohort_A=name_x,
    Cohort_B=name_y,
    N_common_genes=nrow(z),
    Pearson_r=cor(z$x,z$y,method="pearson"),
    Spearman_rho=cor(z$x,z$y,method="spearman"),
    Same_direction_fraction=mean(sign(z$x)==sign(z$y)),
    Same_direction_nonzero_fraction=mean(
      sign(z$x[z$x != 0 & z$y != 0]) ==
      sign(z$y[z$x != 0 & z$y != 0])
    )
  )
}

pairs <- rbindlist(list(
  pair_concordance(a,b,"GSE126848","GSE130970"),
  pair_concordance(a,c,"GSE126848","GSE135251"),
  pair_concordance(b,c,"GSE130970","GSE135251"),
  pair_concordance(b,bs,"GSE130970_PRIMARY6","GSE130970_STRICT4")
))
fwrite(pairs, OUT_PAIR, sep="\t", quote=FALSE, na="NA")

# ------------------------------------------------------------------
# Primary eligibility: gene must be estimable in ALL THREE cohorts.
# Human meta-analysis is external validation, not discovery by missingness.
# ------------------------------------------------------------------
common_primary <- Reduce(
  intersect,
  list(a$ENTREZID, b$ENTREZID, c$ENTREZID)
)
common_strict <- Reduce(
  intersect,
  list(a$ENTREZID, bs$ENTREZID, c$ENTREZID)
)

fit_one <- function(gene, x130) {
  da <- a[ENTREZID==gene]
  db <- x130[ENTREZID==gene]
  dc <- c[ENTREZID==gene]

  d <- rbindlist(list(da,db,dc), use.names=TRUE, fill=TRUE)
  if (nrow(d) != 3L) {
    return(list(status="FAIL", ENTREZID=gene, reason="not_three_cohorts"))
  }

  fit_re <- tryCatch(
    rma.uni(yi=yi, vi=vi, data=d, method="REML", test="z"),
    error=function(e)e
  )
  if (inherits(fit_re,"error")) {
    return(list(
      status="FAIL", ENTREZID=gene,
      reason=paste0("REML:",conditionMessage(fit_re))
    ))
  }

  fit_fe <- tryCatch(
    rma.uni(yi=yi, vi=vi, data=d, method="FE", test="z"),
    error=function(e)e
  )
  if (inherits(fit_fe,"error")) {
    return(list(
      status="FAIL", ENTREZID=gene,
      reason=paste0("FE:",conditionMessage(fit_fe))
    ))
  }

  # Leave-one-cohort-out sign stability using FE pooled effects on 2 cohorts.
  loo_signs <- numeric(3)
  for (j in 1:3) {
    dd <- d[-j]
    ff <- tryCatch(
      rma.uni(yi=yi, vi=vi, data=dd, method="FE", test="z"),
      error=function(e)e
    )
    if (inherits(ff,"error")) {
      loo_signs[j] <- NA_real_
    } else {
      loo_signs[j] <- sign(as.numeric(coef(ff)[1]))
    }
  }

  re_beta <- as.numeric(coef(fit_re)[1])
  cohort_signs <- sign(d$yi)

  out <- data.table(
    ENTREZID=gene,

    human_RE_log2FC=re_beta,
    human_RE_SE=as.numeric(fit_re$se),
    human_RE_z=as.numeric(fit_re$zval),
    human_RE_p=as.numeric(fit_re$pval),
    human_tau2=as.numeric(fit_re$tau2),
    human_I2_percent=as.numeric(fit_re$I2),
    human_QE=as.numeric(fit_re$QE),
    human_QE_p=as.numeric(fit_re$QEp),

    human_FE_log2FC=as.numeric(coef(fit_fe)[1]),
    human_FE_SE=as.numeric(fit_fe$se),
    human_FE_z=as.numeric(fit_fe$zval),
    human_FE_p=as.numeric(fit_fe$pval),

    effect_GSE126848=da$yi,
    effect_GSE130970=db$yi,
    effect_GSE135251=dc$yi,

    direction_consistency_3cohort=max(
      mean(cohort_signs > 0),
      mean(cohort_signs < 0)
    ),
    all_three_same_direction=(
      all(cohort_signs > 0) || all(cohort_signs < 0)
    ),
    primary_RE_sign=sign(re_beta),
    all_LOO_FE_same_sign=(
      all(is.finite(loo_signs)) &&
      all(loo_signs == sign(re_beta))
    )
  )

  list(status="PASS", result=out, ENTREZID=gene, reason="")
}

run_meta <- function(genes, x130, label) {
  cat("Running:", label, "genes=", length(genes), "\n")
  ans <- lapply(genes, fit_one, x130=x130)
  ok <- vapply(ans, function(z) z$status=="PASS", logical(1))
  fail <- !ok

  res <- rbindlist(lapply(ans[ok], `[[`, "result"), fill=TRUE)
  res[, human_RE_FDR_BH := p.adjust(human_RE_p, method="BH")]
  res[, human_FE_FDR_BH := p.adjust(human_FE_p, method="BH")]

  res[, abs_RE_sort := abs(human_RE_log2FC)]
  setorder(res, human_RE_FDR_BH, -abs_RE_sort, ENTREZID)
  res[, abs_RE_sort := NULL]

  failures <- if (any(fail)) {
    rbindlist(lapply(ans[fail], function(z) data.table(
      Analysis=label, ENTREZID=z$ENTREZID, Reason=z$reason
    )))
  } else {
    data.table(Analysis=character(), ENTREZID=character(), Reason=character())
  }

  list(res=res, failures=failures)
}

primary <- run_meta(common_primary, b, "PRIMARY_CONTROL6")
strict <- run_meta(common_strict, bs, "STRICT_CONTROL4")

fwrite(primary$res, OUT_PRIMARY, sep="\t", quote=FALSE, na="NA", compress="gzip")
fwrite(strict$res, OUT_STRICT, sep="\t", quote=FALSE, na="NA", compress="gzip")

fails <- rbindlist(list(primary$failures, strict$failures), fill=TRUE)
fwrite(fails, OUT_FAIL, sep="\t", quote=FALSE, na="NA")

# ------------------------------------------------------------------
# Primary vs strict4 sensitivity comparison.
# ------------------------------------------------------------------
cmp <- merge(
  primary$res[, .(
    ENTREZID,
    primary_RE_log2FC=human_RE_log2FC,
    primary_RE_FDR=human_RE_FDR_BH,
    primary_all3_same=all_three_same_direction
  )],
  strict$res[, .(
    ENTREZID,
    strict4_RE_log2FC=human_RE_log2FC,
    strict4_RE_FDR=human_RE_FDR_BH,
    strict4_all3_same=all_three_same_direction
  )],
  by="ENTREZID"
)
cmp[, Same_RE_direction :=
      sign(primary_RE_log2FC)==sign(strict4_RE_log2FC)]
fwrite(cmp, OUT_COMPARE, sep="\t", quote=FALSE, na="NA", compress="gzip")

sens_pearson <- cor(cmp$primary_RE_log2FC, cmp$strict4_RE_log2FC, method="pearson")
sens_spearman <- cor(cmp$primary_RE_log2FC, cmp$strict4_RE_log2FC, method="spearman")
sens_same <- mean(cmp$Same_RE_direction)

summ <- function(r) {
  list(
    n=nrow(r),
    fdr05=sum(r$human_RE_FDR_BH < .05, na.rm=TRUE),
    up=sum(r$human_RE_FDR_BH < .05 & r$human_RE_log2FC > 0, na.rm=TRUE),
    down=sum(r$human_RE_FDR_BH < .05 & r$human_RE_log2FC < 0, na.rm=TRUE),
    all3=sum(r$all_three_same_direction, na.rm=TRUE),
    loo=sum(r$all_LOO_FE_same_sign, na.rm=TRUE),
    medianI2=median(r$human_I2_percent, na.rm=TRUE),
    medianTau=median(r$human_tau2, na.rm=TRUE)
  )
}
sp <- summ(primary$res)
ss <- summ(strict$res)

summary <- c(
  "===== HUMAN MASLD GENE META-ANALYSIS SUMMARY v1 =====",
  "",
  "STATUS=PASS",
  paste0("R_version=", R.version.string),
  paste0("metafor_version=", as.character(packageVersion("metafor"))),
  "",
  "PRIMARY_ELIGIBILITY=gene_estimable_in_all_3_human_cohorts",
  "PRIMARY_MODEL=random_effects_REML",
  "PRIMARY_TEST=z",
  "FIXED_EFFECT_SENSITIVITY=YES",
  "HUMAN_PREDICTION_INTERVAL_USED=NO_k_equals_3",
  "LOO_SIGN_STABILITY=two_cohort_fixed_effect_sign",
  "",
  paste0("Primary_common_genes=", sp$n),
  paste0("Primary_RE_FDR_lt_0.05=", sp$fdr05),
  paste0("Primary_RE_FDR05_up=", sp$up),
  paste0("Primary_RE_FDR05_down=", sp$down),
  paste0("Primary_all3_same_direction=", sp$all3),
  paste0("Primary_all_LOO_FE_same_sign=", sp$loo),
  paste0("Primary_median_I2_percent=", format(sp$medianI2,digits=6)),
  paste0("Primary_median_tau2=", format(sp$medianTau,digits=6)),
  "",
  paste0("Strict4_common_genes=", ss$n),
  paste0("Strict4_RE_FDR_lt_0.05=", ss$fdr05),
  paste0("Strict4_RE_FDR05_up=", ss$up),
  paste0("Strict4_RE_FDR05_down=", ss$down),
  paste0("Strict4_all3_same_direction=", ss$all3),
  paste0("Strict4_all_LOO_FE_same_sign=", ss$loo),
  "",
  paste0("Primary_vs_strict4_common_genes=", nrow(cmp)),
  paste0("Primary_vs_strict4_RE_effect_Pearson_r=", format(sens_pearson,digits=6)),
  paste0("Primary_vs_strict4_RE_effect_Spearman_rho=", format(sens_spearman,digits=6)),
  paste0("Primary_vs_strict4_same_direction_fraction=", format(sens_same,digits=6)),
  "",
  "PAIRWISE_COHORT_EFFECT_CONCORDANCE:"
)
for (i in seq_len(nrow(pairs))) {
  rr <- pairs[i]
  summary <- c(summary, paste0(
    rr$Cohort_A,"_vs_",rr$Cohort_B,
    "|N=",rr$N_common_genes,
    "|Pearson=",format(rr$Pearson_r,digits=5),
    "|Spearman=",format(rr$Spearman_rho,digits=5),
    "|SameDirection=",format(rr$Same_direction_fraction,digits=5)
  ))
}

summary <- c(
  summary, "",
  paste0("Primary_result=", sub(paste0("^",ROOT,"/"),"",OUT_PRIMARY)),
  paste0("Strict4_result=", sub(paste0("^",ROOT,"/"),"",OUT_STRICT)),
  paste0("Pairwise_concordance=", sub(paste0("^",ROOT,"/"),"",OUT_PAIR)),
  paste0("Primary_vs_strict4=", sub(paste0("^",ROOT,"/"),"",OUT_COMPARE)),
  paste0("Model_failures=", sub(paste0("^",ROOT,"/"),"",OUT_FAIL)),
  "",
  "MOUSE_HUMAN_ORTHOLOG_MAPPING_PERFORMED=NO",
  "MOUSE_HUMAN_CONCORDANCE_TEST_PERFORMED=NO"
)

writeLines(summary, OUT_SUM)
cat("\n", paste(summary, collapse="\n"), "\n", sep="")

if (nrow(primary$failures) > 0L || nrow(strict$failures) > 0L)
  quit(status=2L)

